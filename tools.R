library(here)
library(googledrive)
library(googlesheets)
library(readr)
library(tools)
library(stringr)
library(stats)
library(dplyr)

source(file.path(here(), "configuration.R"))


read_fasta <- function(file_path) {
  # Read raw string
  raw_data <- read_file(file_path)
  
  # Return an empty vector an a warning if no sequences are found
  if (raw_data == "") {
    warning(paste0("No sequences found in the file: ", file_path))
    return(character(0))
  }
  
  # Find location of every header start 
  split_data <- str_split(raw_data, pattern = "\n>", simplify = TRUE)
  
  # Split the data for each sequence into lines
  split_data <- str_split(split_data, pattern = "\n")
  
  # The first lines are headers, so remvove those
  headers <- vapply(split_data, FUN = `[`, FUN.VALUE = character(1), 1)
  split_data <- lapply(split_data, FUN = `[`, -1)
  
  # Remove the > from the first sequence. The others were removed by the split
  headers[1] <- sub(headers[1], pattern = "^>", replacement = "")
  
  # Combine multiple lines into single sequences
  seqs <- vapply(split_data, FUN = paste0, FUN.VALUE = character(1), collapse = "")
  
  # Combine and return results 
  return(setNames(seqs, headers))
}


get_release_data <- function() {
  data <- read_csv(local_release_spreadsheet_path, col_type = cols())
  data$release_name <- paste0(release_name_prefix, data$release_number)
  data$release_path <- file.path(local_release_dir, paste0(data$release_name, ".fa"))
  data$blast_path <- file.path(blast_database_dir, data$release_name)
  count_data <- system2("grep", c("-c", "'^>'", file.path(local_release_dir, "*.fa")), stdout = TRUE) %>%
    paste0(collapse = "\n") %>%
    read_delim(":", col_names = c("release_name", "seq_count")) %>%
    mutate(release_name = file_path_sans_ext(basename(release_name)))
  data <- left_join(data, count_data, by = "release_name")
  return(data)
}


get_public_release_names <- function() {
  release_data <- get_release_data()
  release_data <- release_data[release_data$public, ]
  return(release_data$release_name)
}


make_blast_database <- function(fasta_path, out_dir_path) {
  # Make temporary version of the database with indexes instead of full headers
  #   This is needed because BLAST will not return header info after the first space..
  ref_seqs <- read_fasta(fasta_path)
  names(ref_seqs) <- as.character(seq_along(ref_seqs))
  temp_database_path <- tempfile()
  write_lines(paste0(">", names(ref_seqs), "\n", ref_seqs), temp_database_path)
  
  # Make blast database with temporary file
  blast_db_name <- file_path_sans_ext(basename(fasta_path))
  blast_db_path <- file.path(out_dir_path, blast_db_name)
  makeblastdb_command <- paste("makeblastdb",
                               "-in", temp_database_path,
                               "-parse_seqids",
                               "-dbtype nucl",
                               "-out", blast_db_path)
  system(makeblastdb_command)
}

update_releases <- function() {
  # turn off authentication for google drive
  drive_auth_config(active = FALSE)
  
  # Update local spreadsheet
  drive_files <- drive_ls(path = as_id(googledrive_root_id))
  if (! release_spreadsheet_name %in% drive_files$name) {
    stop('Cant find release spreadsheet named "', release_spreadsheet_name,
         '". Check that the file exists in Google Drive or change the name of the file to look for in "configuration.R".')
  }
  releases_spreadsheet_id <- drive_files$id[drive_files$name == release_spreadsheet_name]
  release_spreadsheet_obj <- gs_key(releases_spreadsheet_id)
  gs_download(release_spreadsheet_obj, to = local_release_spreadsheet_path, overwrite = TRUE)
  release_data <- read_csv(local_release_spreadsheet_path)
 
  # Check for new releases
  local_release_names <- list.files(local_release_dir, pattern = paste0(release_name_prefix, "[0-9]+\\.fa"))
  local_release_nums <- str_match(local_release_names, pattern = paste0(release_name_prefix, "([0-9]+)\\.fa"))[, 2]
  new_release_indexes <- which(! release_data$release_number %in% local_release_nums)
  message('Adding ', length(new_release_indexes), ' releases.')
  
  # Process new releases
  db_file_data <- drive_ls(path = as_id(googledrive_database_id))
  make_one_release <- function(index) {
    # Store a local copy of the database
    new_release_remote <- release_data$source_file[index]
    if (! new_release_remote %in% db_file_data$name) {
      stop('Cannot find remote release file "', new_release_remote, '".')
    }
    new_release_id <- db_file_data$id[db_file_data$name == new_release_remote]
    new_release_file_name <- paste0(release_name_prefix, release_data$release_number[index], ".fa")
    new_release_file_path <- file.path(local_release_dir, new_release_file_name)
    drive_download(as_id(new_release_id), path = new_release_file_path)
    
    # Create a blast database for it
    make_blast_database(fasta_path = new_release_file_path, out_dir_path = blast_database_dir)
  }
  
  not_used <- lapply(new_release_indexes, make_one_release)
}



get_blast_databases <- function(db_dir) {
  # Get possible blast database names
  db_file_names <- file_path_sans_ext(list.files(db_dir))
  output <- unique(db_file_names)
  
  # Remove any that do not appear the correct number of times
  is_database <- vapply(output, function(x) as.list(table(db_file_names))[[x]] == 6, logical(1))
  output <- output[is_database]
  
  return(rev(output))
}


get_latest_release_fa <- function() {
  local_release_names <- list.files(local_release_dir, pattern = paste0(release_name_prefix, "[0-9]+\\.fa"))
  local_release_nums <- str_match(local_release_names, pattern = paste0(release_name_prefix, "([0-9]+)\\.fa"))[, 2]
  local_release_nums <- as.numeric(local_release_nums)
  file.path(local_release_dir, local_release_names[which.max(local_release_nums)])
}

convert_oomydb_headers_to_taxmap <- function(headers, ...) {
  extract_tax_data(headers,
                   class_sep = ";",
                   regex = "name=(.+)\\|strain=(.+)\\|ncbi_acc=(.+)\\|ncbi_taxid=(.+)\\|oodb_id=(.+)\\|taxonomy=(.+)$",
                   key = c(name = "info", strain = "info", ncbi_acc = "info", ncbi_taxid = "info", oodb_id = "info", taxonomy = "class"),
                   ...)
}
