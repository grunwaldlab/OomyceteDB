library(here)
library(googledrive)
library(googlesheets)

source(file.path(here(), "configuration.R"))


make_blast_database <- function(fasta_path, out_dir_path) {
  # Make temporary version of the database with indexes instead of full headers
  #   This is needed because BLAST will not return header info after the first space..
  ref_seqs <- metacoder::read_fasta(fasta_path)
  names(ref_seqs) <- as.character(seq_along(ref_seqs))
  temp_database_path <- tempfile()
  readr::write_lines(paste0(">", names(ref_seqs), "\n", ref_seqs), temp_database_path)
  
  # Make blast database with temporary file
  blast_db_name <- tools::file_path_sans_ext(basename(fasta_path))
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
  googledrive::drive_auth_config(active = FALSE)
  
  # Update local spreadsheet
  drive_files <- drive_ls(path = as_id(googledrive_root_id))
  if (! release_spreadsheet_name %in% drive_files$name) {
    stop('Cant find release spreadsheet named "', release_spreadsheet_name,
         '". Check that the file exists in Google Drive or change the name of the file to look for in "configuration.R".')
  }
  releases_spreadsheet_id <- drive_files$id[drive_files$name == release_spreadsheet_name]
  release_spreadsheet_obj <- googlesheets::gs_key(releases_spreadsheet_id)
  googlesheets::gs_download(release_spreadsheet_obj, to = local_release_spreadsheet_path, overwrite = TRUE)
  release_data <- readr::read_csv(local_release_spreadsheet_path)
 
  # Check for new releases
  local_release_names <- list.files(local_release_dir, pattern = paste0(release_name_prefix, "[0-9]+\\.fa"))
  local_release_nums <- stringr::str_match(local_release_names, pattern = paste0(release_name_prefix, "([0-9]+)\\.fa"))[, 2]
  remote_release_names <- release_data$release_number
  new_release_indexes <- which(! remote_release_names %in% local_release_nums)
  message('Adding ', length(new_release_indexes), ' releases.')
  
  # Process new releases
  make_one_release <- function(index) {
    # Store a local copy of the database
    new_release_remote <- release_data$source_file[index]
    if (! new_release_remote %in% drive_files$name) {
      stop('Cannot find remote release file "', new_release_remote, '".')
    }
    new_release_id <- drive_files$id[drive_files$name == new_release_remote]
    new_release_file_name <- paste0(release_name_prefix, release_data$release_number[index], ".fa")
    new_release_file_path <- file.path(local_release_dir, new_release_file_name)
    googledrive::drive_download(as_id(new_release_id), path = new_release_file_path)
    
    # Create a blast database for it
    make_blast_database(fasta_path = new_release_file_path, out_dir_path = blast_database_dir)
  }
  
  not_used <- lapply(new_release_indexes, make_one_release)
}



get_blast_databases <- function(db_dir) {
  # Get possible blast database names
  db_file_names <- tools::file_path_sans_ext(list.files(db_dir))
  output <- unique(db_file_names)
  
  # Remove any that do not appear the correct number of times
  is_database <- vapply(output, function(x) as.list(table(db_file_names))[[x]] == 6, logical(1))
  output <- output[is_database]
  
  return(rev(output))
}


get_latest_release_fa <- function() {
  local_release_names <- list.files(local_release_dir, pattern = paste0(release_name_prefix, "[0-9]+\\.fa"))
  local_release_nums <- stringr::str_match(local_release_names, pattern = paste0(release_name_prefix, "([0-9]+)\\.fa"))[, 2]
  local_release_nums <- as.numeric(local_release_nums)
  file.path(local_release_dir, local_release_names[which.max(local_release_nums)])
}
