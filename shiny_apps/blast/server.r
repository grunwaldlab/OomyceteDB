require(XML)
library(plyr)
library(dplyr)
library(DT)
library(metacoder)
library(shinyjs)

outfmt_options <- c(" 0: pairwise",
                    " 1: query-anchored showing identities", 
                    " 2: query-anchored no identities",
                    " 3: flat query-anchored show identities", 
                    " 4: flat query-anchored no identities",
                    " 5: XML Blast output", 
                    " 6: tabular",
                    " 7: tabular with comment lines",
                    " 8: Text ASN.1", 
                    " 9: Binary ASN.1",
                    "10: Comma-separated values",
                    "11: BLAST archive format (ASN.1)", 
                    "12: JSON Seqalign output",
                    "13: JSON Blast output",
                    "14: XML2 Blast output")

server <- function(input, output, session) {
  
  raw_blast_results <- eventReactive(
    eventExpr = input$blast,
    ignoreNULL = TRUE,
    {
      
      # Choose database (NOT ACTUALLY CHOOSING YET)
      db <- file.path("../../data/blast_databases", input$db)
      remote <- c("")
      
      # Format query
      query_temp_path <- tempfile(fileext = ".fa")
      if (startsWith(input$query, ">")){
        writeLines(input$query, query_temp_path)
      } else {
        writeLines(paste0(">Query\n", input$query), query_temp_path)
      }
      
      # Calls BLAST
      blast_command <- paste(input$program,
                             "-query", query_temp_path,
                             "-db", db ,
                             "-dust no",
                             "-evalue", input$eval,
                             "-outfmt 11",
                             "-max_hsps 1",
                             "-max_target_seqs", input$max_target_seqs,
                             remote)
      system(blast_command, intern = TRUE)
    })
  
  
  blast_results <- eventReactive(
    eventExpr = input$blast,
    ignoreNULL = TRUE,
    {
      # Convert to XML format
      raw_results <- raw_blast_results()
      tmp_file <- tempfile()
      readr::write_lines(raw_results, path = tmp_file)
      blast_formatter_cmd <- paste("blast_formatter",
                                   "-archive", tmp_file,
                                   "-outfmt 5")
      xmlParse(system(blast_formatter_cmd, intern = TRUE))
    })
  

  blast_results_for_download <- eventReactive(
    eventExpr = input$outfmt,
    ignoreNULL = TRUE,
    {
      # Convert to XML format
      raw_results <- raw_blast_results()
      tmp_file <- tempfile()
      readr::write_lines(raw_results, path = tmp_file)
      blast_formatter_cmd <- paste("blast_formatter",
                                   "-archive", tmp_file,
                                   "-outfmt", stringr::str_match(input$outfmt, " *([0-9]+):.+")[, 2])
      system(blast_formatter_cmd, intern = TRUE)
    })
  
  # Now to parse the results...
  parsed_results <- reactive({
    if (is.null(blast_results()))
    {
      
    } else {

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
      
      # Convert to taxmap
      tm_obj <- taxa::extract_tax_data(results$hit_ids,
                                       class_sep = ";",
                                       regex = "id=(.+)\\|name=(.+)\\|source=(.+)\\|tax_id=(.+)\\|taxonomy=(.+)$",
                                       key = c(id = "info", name = "info", source = "info", tax_id = "info", taxonomy = "class"))
      
      # Add on results table to taxmap
      tm_obj$data$tax_data <- bind_cols(tm_obj$data$tax_data, results)
      
      return(tm_obj)
    }
  })
  
  # makes the datatable
  output$blast_results <- renderDataTable({
    if (is.null(blast_results())) {
    } else {
      results <- parsed_results()
      classification <- results %>% 
        filter_taxa(taxon_names == "Oomycetes", subtaxa = TRUE) %>%
        filter_taxa(!is_root) %>%
        classifications() %>% 
        `[`(results$data$tax_data$taxon_id)
        
      
      results$data$tax_data %>% 
        transmute("Query ID" = query_id,
                  "Hit taxonomic classification" = classification,
                  "Identity (%)" = round(prop_identity * 100, digits = 3),
                  "Query Coverage (%)" = round(prop_match_len * 100, digits = 3))
    }
  }, selection = "single")
  
  # this chunk gets the alignemnt information from a clicked row
  output$clicked <- renderTable({
    if(is.null(input$blast_results_rows_selected)){}
    else{
      clicked <- input$blast_results_rows_selected
      results <- parsed_results()
      tableout <- data.frame(results$data$tax_data[clicked,])
      tableout <- transmute(tableout,
                            "Hit taxonomic classification" = taxonomy,
                            "Hit NCBI taxon ID" = tax_id,
                            # accesion number when we have it
                            "Identity (%)" = round(prop_identity * 100, digits = 3),
                            "Alignment length" = align_len,
                            "Mismatches" = align_len - identity, #check 
                            "Gaps" =  gaps,
                            "Query aligned range" = paste(query_from, "-", query_to),
                            "Hit aligned range" = paste(hit_from, "-", hit_to),
                            "E value" = evalue,
                            "Bit score"  = bit_score,
                            "Query ID" = query_id,
                            # "Hit ID" = hit_ids,
                            "Query Coverage (%)" = round(prop_match_len * 100, digits = 3),
                            "Matching base pairs" = identity,
                            "Hit length" = hit_length
                            # "Hit sequence" = hit_seq
                            )
      
      tableout <- t(tableout)
      names(tableout) <- c("")
      colnames(tableout) <- NULL
      data.frame(tableout)
    }
  },rownames = TRUE,colnames = FALSE)
  
  #this chunk makes the alignments for clicked rows
  output$alignment <- renderText({
    if(is.null(input$blast_results_rows_selected)){}
    else{
      clicked = input$blast_results_rows_selected
        
      #loop over the xml to get the alignments
      align <- xpathApply(blast_results(), '//Iteration',function(row){
        top <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq') %>% sapply(., xmlValue)
        mid <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_midline') %>% sapply(., xmlValue)
        bottom <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq') %>% sapply(., xmlValue)
        rbind(top,mid,bottom)
      })
        
      #split the alignments every 80 carachters to get a "wrapped look"
      alignx <- do.call("cbind", align)
      splits <- strsplit(gsub("(.{80})", "\\1,", alignx[1:3,clicked]),",")
      
      #paste them together with returns '\n' on the breaks
      split_out <- lapply(seq_along(splits[[1]]), function(i){
        rbind(paste0("Q-", splits[[1]][i], "\n"),
              paste0("  ", splits[[2]][i], "\n"),
              paste0("H-", splits[[3]][i], "\n"),
              "\n")
      })
      output <- unlist(split_out)
      output[1] <- paste0(" ", output[1])
      return(output)
    }
  })
  
  # Downloadable csv of selected dataset ----
  output$download_data <- downloadHandler(
    contentType = "text/csv",
    filename = function() {
      paste0("blast_results_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".txt")
    },
    content = function(file) {
      print("heres")
      not_used <- readr::write_lines(blast_results_for_download(), path = file)
    }
  )
  
  output$download_data_form <- renderUI({
    if (is.null(blast_results())) {
    } else {
      list(
        h4("Download BLAST output"),
        selectInput("outfmt", "Output format:", choices = outfmt_options, width = "300px",
                    selected = "10: Comma-separated values"),
        downloadButton(outputId = "download_data", label = "Download Results")
      )
    }
  })
  
  output$hit_list <- renderUI({
    if (is.null(blast_results())) {
    } else {
      list(
        h4("Hits"),
        DT::dataTableOutput("blast_results")
        )
    }
  })
  
  output$hit_details <- renderUI({
    if (is.null(blast_results())) {
    } else {
      if (is.null(input$blast_results_rows_selected)) {
        return(list(
          h4("Hit details"),
          p("Click on a row for details...")
        ))
      }
      return(list(
        h4("Hit details"),
        p(tableOutput("clicked")),
        verbatimTextOutput("alignment")
      ))
    }
  })
  
}

