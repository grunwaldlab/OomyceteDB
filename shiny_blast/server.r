require(XML)
library(plyr)
library(dplyr)
library(DT)

server <- function(input, output, session) {
  
  blast_results <- eventReactive(
    eventExpr = input$blast,
    ignoreNULL = TRUE,
    {
      
      # Choose database (NOT ACTUALLY CHOOSING YET)
      db <- file.path("../data/blast_databases", input$db)
      remote <- c("")
      
      # Format query
      query_temp_path <- tempfile(fileext = ".fa")
      if (startsWith(input$query, ">")){
        writeLines(input$query, query_temp_path)
      } else {
        writeLines(paste0(">Query\n", input$query), query_temp_path)
      }
      
      #calls the blast
      blast_command <- paste(input$program,
                             "-query", query_temp_path,
                             "-db", db ,
                             "-dust no",
                             "-evalue", input$eval,
                             "-outfmt 5",
                             "-max_hsps 1",
                             "-max_target_seqs 10",
                             remote)
      data <- system(blast_command, intern = TRUE)
      xmlParse(data)
    })

  # Now to parse the results...
  parsed_results <- reactive({
    if (is.null(blast_results()))
    {
      
    } else {
      xmltop = xmlRoot(blast_results())
      
      # the first chunk is for multi-fastas
      results <- xpathApply(blast_results(), '//Iteration', function(row) {
        
        get_hit_info <- function(info) {
          getNodeSet(row, paste0('Iteration_hits//Hit//Hit_hsps//Hsp//', info)) %>% sapply(., xmlValue)
        }
        
        query_id <-   getNodeSet(row, 'Iteration_query-def') %>% sapply(., xmlValue)
        query_len <-  getNodeSet(row, 'Iteration_query-len') %>% sapply(., xmlValue)
        hit_ids <-    getNodeSet(row, 'Iteration_hits//Hit//Hit_id') %>% sapply(., xmlValue)
        hit_length <- getNodeSet(row, 'Iteration_hits//Hit//Hit_len') %>% sapply(., xmlValue)
        bit_score <-  get_hit_info('Hsp_bit-score')
        score <-      get_hit_info('Hsp_score') 
        evalue <-     get_hit_info('Hsp_evalue')
        query_from <- get_hit_info('Hsp_query-from')
        query_to <-   get_hit_info('Hsp_query-to')
        hit_from <-   get_hit_info('Hsp_hit-from')
        hit_to <-     get_hit_info('Hsp_hit-to')
        identity <-   get_hit_info('Hsp_identity')
        positive <-   get_hit_info('Hsp_positive')
        gaps <-       get_hit_info('Hsp_gaps')
        align_len <-  get_hit_info('Hsp_align-len')
        query_seq <-  get_hit_info('Hsp_qseq')
        hit_seq <-    get_hit_info('Hsp_hseq')
        midline <-    get_hit_info('Hsp_midline')
        
        cbind(query_id, query_len, hit_ids, hit_length, bit_score, score, evalue,
              query_from, query_to, hit_from, hit_to, identity, positive, gaps,
              align_len, query_seq, hit_seq, midline)
      })
      
      # this ensures that NAs get added for no hits
      results <- rbind.fill(lapply(results, function(y) {as.data.frame(y, stringsAsFactors = FALSE)}))
      
      # Convert numeric cols to numeric
      numeric_cols <- c("query_len", "hit_length", "bit_score", "score", "evalue", "query_from", 
                        "query_to", "hit_from", "hit_to", "identity", "positive", "gaps", 
                        "align_len")
      results[numeric_cols] <- lapply(results[numeric_cols], as.numeric)
      
      # Replace database index with header
      database_seqs <- metacoder::read_fasta("../data/source/rps10_database.fa") # SHOULD BE SELECTED BY INPUT OPTION!!!
      results$hit_ids <- names(database_seqs)[as.numeric(results$hit_ids)]
      
      # Calculate derived columns
      results$prop_identity <- results$identity / results$align_len
      results$prop_match_len <- (results$query_to - results$query_from + 1) / results$query_len
      
      return(results)
    }
  })
  
  # makes the datatable
  output$blast_results <- renderDataTable({
    if (is.null(blast_results())) {
    } else {
      results <- parsed_results()
      tm_obj <- taxa::extract_tax_data(results$hit_ids,
                                       class_sep = ";",
                                       regex = "id=(.+)\\|name=(.+)\\|source=(.+)\\|tax_id=(.+)\\|taxonomy=(.+)$",
                                       key = c(id = "info", name = "info", source = "info", tax_id = "info", taxonomy = "class"))
      results$classification <- tm_obj %>% 
        filter_taxa(taxon_names == "Oomycetes", subtaxa = TRUE) %>%
        filter_taxa(!is_root) %>%
        classifications() %>% 
        `[`(tm_obj$data$tax_data$taxon_id)
        
      
      results %>% 
        transmute("Query ID" = query_id,
                  "Hit taxonomic classification" = classification,
                  "Identity" = round(prop_identity, digits = 3),
                  "Query Matched" = round(prop_match_len, digits = 3))
    }
  }, selection = "single")
  
  # this chunk gets the alignemnt information from a clicked row
  output$clicked <- renderTable({
    if(is.null(input$blast_results_rows_selected)){}
    else{
      xmltop = xmlRoot(blast_results())
      clicked = input$blast_results_rows_selected
      tableout<- data.frame(parsed_results()[clicked,])
      
      tableout <- t(tableout)
      names(tableout) <- c("")
      rownames(tableout) <- c("Query ID","Hit ID", "Length", "Bit Score", "e-value")
      colnames(tableout) <- NULL
      data.frame(tableout)
    }
  },rownames =T,colnames =F)
  
  #this chunk makes the alignments for clicked rows
  output$alignment <- renderText({
    if(is.null(input$blast_results_rows_selected)){}
    else{
      xmltop = xmlRoot(blast_results())
      
      clicked = input$blast_results_rows_selected
        
      #loop over the xml to get the alignments
      align <- xpathApply(blast_results(), '//Iteration',function(row){
        top <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq') %>% sapply(., xmlValue)
        mid <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_midline') %>% sapply(., xmlValue)
        bottom <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq') %>% sapply(., xmlValue)
        rbind(top,mid,bottom)
      })
        
      #split the alignments every 40 carachters to get a "wrapped look"
      alignx <- do.call("cbind", align)
      splits <- strsplit(gsub("(.{40})", "\\1,", alignx[1:3,clicked]),",")
      
      #paste them together with returns '\n' on the breaks
      split_out <- lapply(1:length(splits[[1]]),function(i){
        rbind(paste0("Q-",splits[[1]][i],"\n"),paste0("M-",splits[[2]][i],"\n"),paste0("H-",splits[[3]][i],"\n"))
      })
      unlist(split_out)
    }
  })
}
#phew...
