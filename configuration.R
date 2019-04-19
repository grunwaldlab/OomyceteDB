library(here)

googledrive_root_id <- "https://drive.google.com/drive/u/1/folders/1XJLFW_S9EzVOGgHCqg42H88dNtD91gRq"
release_spreadsheet_name <- "releases"
local_release_spreadsheet_path <- file.path(here(), "data/releases.csv")
release_name_prefix <- "release_"
local_release_dir <- file.path(here(), "data/releases")
blast_database_dir <- file.path(here(), "data/blast_databases")