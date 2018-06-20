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
  }, selection = "multiple")
  
  
  output$selected_seq_printout <- renderText({
    if(is.null(input$database_list_rows_selected)){}
    else{
      results <- selected_subset()
      clicked <- input$database_list_rows_selected
      output <- paste0(paste0(">", results$data$tax_data$input[clicked], "\n",
                       results$data$tax_data$sequence[clicked]), collapse = "\n")
      return(output)
    }
  })
  
  
  database_for_download <- eventReactive(
    eventExpr = input$search,
    ignoreNULL = TRUE,
    {
      results <- selected_subset()
      paste0(">", results$data$tax_data$input, "\n",
             results$data$tax_data$sequence)
    })
  
  
  database_for_download_selected <- eventReactive(
    eventExpr = input$database_list_rows_selected,
    ignoreNULL = TRUE,
    {
      results <- selected_subset()
      paste0(">", results$data$tax_data$input[input$database_list_rows_selected], "\n",
             results$data$tax_data$sequence[input$database_list_rows_selected])
    })
  
  
  output$download_data <- downloadHandler(
    contentType = "text/csv",
    filename = function() {
      paste0("oomycetedb_subset_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".fa")
    },
    content = function(file) {
      not_used <- readr::write_lines(database_for_download(), path = file)
    }
  )
  
  output$download_data_selected <- downloadHandler(
    contentType = "text/csv",
    filename = function() {
      paste0("oomycetedb_subset_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".fa")
    },
    content = function(file) {
      not_used <- readr::write_lines(database_for_download_selected(), path = file)
    }
  )
  
  
  output$download_data_form <- renderUI({
    if (is.null(selected_subset())) {
    } else {
      if (is.null(input$database_list_rows_selected)) {
        return(list(
          # h4("Download database subset"),
          downloadButton(outputId = "download_data", label = "Download database subset")
        ))
      } else {
        return(list(
          # h4("Download database subset"),
          downloadButton(outputId = "download_data", label = "Download database subset"),
          downloadButton(outputId = "download_data_selected", label = "Download selected sequences")
        ))
        
      }
    }
  })
  
  output$seq_list_table <- renderUI({
    if (is.null(selected_subset())) {
    } else {
      list(
        h4("Sequences in subset"),
        DT::dataTableOutput("database_list")
      )
    }
  })
  
  output$seq_details <- renderUI({
    x = 1
    if (is.null(selected_subset())) {
    } else {
      if (is.null(input$database_list_rows_selected)) {
        return(list(
          h4("Selected sequences"),
          p("Click on one or more rows to get the FASTA entry. Click again to deselect.")
        ))
      } else {
        return(list(
          h4("Selected sequences"),
          verbatimTextOutput("selected_seq_printout")
        ))
      }
    }
  })
}

