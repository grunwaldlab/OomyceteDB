#!/usr/bin/env Rscript

# Import dependencies
library(argparser, quietly = TRUE)
library(readxl, quietly = TRUE)

# Create command line argument parser and parse input
parser <- arg_parser("Convert oomycetedb excel file to FASTA for release. This assumes the following columns are present: oodb_id, sequence, genus, species, isolate/strain, genbank_id, taxon_id, source, Order, Family")
parser <- add_argument(parser, "xlsx_path", help="Path to excel file with data to convert to FASTA.", type="character")
parser <- add_argument(parser, "out_path", help="Path to write the FASTA file to.", type="character")
arguments <- parse_args(parser)
# arguments <- list(xlsx_path = "~/downloads/DB_release_1.2.xlsx", out_path = 'deleteme.fasta')

# Parse excel file
input_data <- readxl::read_xlsx(arguments$xlsx_path)
colnames(input_data) <- tolower(colnames(input_data))
colnames(input_data)[colnames(input_data) == 'isolate/strain'] <- 'strain'

# Convert to FASTA
input_data$binomial <- paste0(input_data$genus, '_', input_data$species)
input_data$taxonomy <- paste('cellular_organisms;Eukaryota;Stramenopiles;Oomycetes', input_data$order, input_data$family, input_data$binomial, sep = ';')
col_key <- c(
  binomial = 'name',
  strain = 'strain',
  genbank_id = 'ncbi_acc',
  taxon_id = 'ncbi_taxid',
  oodb_id = 'oodb_id',
  taxonomy = 'taxonomy'
)
headers <- vapply(seq_len(nrow(input_data)), function(i) {
  paste0(col_key, '=', input_data[i, names(col_key)], collapse = '|')  
}, FUN.VALUE = character(1))
output <- paste0('>', headers, '\n', input_data$sequence)
writeLines(output, arguments$out_path)
