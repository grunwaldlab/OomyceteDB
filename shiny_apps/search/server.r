require(XML)
library(plyr)
library(dplyr)
library(DT)
library(metacoder)
library(shinyjs)

server <- function(input, output, session) {
  
  selected_subset <- eventReactive(
    eventExpr = input$search,
    ignoreNULL = TRUE,
    {
      # Parse whole database
      database_seqs <- metacoder::read_fasta("../../data/source/rps10_database.fa") # SHOULD BE SELECTED BY INPUT OPTION!!!
      tm_obj <- taxa::extract_tax_data(names(database_seqs),
                                       include_match = FALSE,
                                       class_sep = ";",
                                       regex = "id=(.+)\\|name=(.+)\\|source=(.+)\\|tax_id=(.+)\\|taxonomy=(.+)$",
                                       key = c(id = "info", name = "info", source = "info", tax_id = "info", taxonomy = "class"))
      tm_obj$data$tax_data$sequence <- database_seqs
      
      # Subset to taxon
      if (input$taxon_subset != "") {
        taxa_to_use <- trimws(strsplit(input$taxon_subset, split = ",+")[[1]])
        tm_obj <- filter_taxa(tm_obj, tolower(taxon_names) %in% tolower(taxa_to_use),
                              subtaxa = TRUE, supertaxa = TRUE, drop_obs = TRUE,
                              reassign_obs = FALSE)
      }
      
      
      return(tm_obj)
    })
  
  
  # makes the datatable
  output$database_list <- renderDataTable({
    if (is.null(selected_subset())) {
    } else {
      results <- selected_subset()
      
      classification <- results %>% 
        filter_taxa(taxon_names == "Oomycetes", subtaxa = TRUE) %>%
        filter_taxa(!is_root) %>%
        classifications() %>% 
        `[`(results$data$tax_data$taxon_id)
      
      
      results$data$tax_data %>% 
        transmute("Sequence ID" = id,
                  "Organism name" = name,
                  # "NCBI taxon ID" = tax_id,
                  "Taxonomy" = classification)
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
  }, rownames = TRUE, colnames = FALSE)
  
  
  database_for_download <- eventReactive(
    eventExpr = input$search,
    ignoreNULL = TRUE,
    {
      results <- selected_subset()
      paste0(">", results$data$tax_data$input, "\n", results$data$tax_data$sequence)
    })
  
  
  output$download_data <- downloadHandler(
    contentType = "text/csv",
    filename = function() {
      paste0("oomycetedb_subset_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".fa")
    },
    content = function(file) {
      print("heres")
      not_used <- readr::write_lines(database_for_download(), path = file)
    }
  )
  
  
  output$download_data_form <- renderUI({
    if (is.null(selected_subset())) {
    } else {
      list(
        h4("Download database subset"),
        downloadButton(outputId = "download_data", label = "Download FASTA")
      )
    }
  })
  
  
}

