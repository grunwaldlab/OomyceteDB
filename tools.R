make_blast_database <- function() {
  # Make temporary version of the database with indexes instead of full headers
  #   This is needed because BLAST will not return header info after the first space..
  ref_seqs <- metacoder::read_fasta("data/source/rps10_database.fa")
  names(ref_seqs) <- as.character(seq_along(ref_seqs))
  temp_database_path <- tempfile()
  readr::write_lines(paste0(">", names(ref_seqs), "\n", ref_seqs), temp_database_path)
  
  # Make blast database with temporary file
  makeblastdb_command <- paste("makeblastdb",
                               "-in", temp_database_path,
                               "-parse_seqids",
                               "-dbtype nucl",
                               "-out data/blast_databases/current")
  system(makeblastdb_command)
}