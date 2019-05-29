library(XML)
library(dplyr)
library(DT)
library(taxa)
library(shinyjs)
library(here)

source(file.path(here(), "tools.R"))

options(shiny.sanitize.errors = FALSE)

server <- function(input, output, session) {
  
  
  selected_database_name <- reactive({
    data <- public_release_data()
    data$release_name[input$database_table_rows_selected]
  })
  
  selected_database <- reactive(
    {
      database_path <- file.path(local_release_dir, paste0(selected_database_name(), ".fa"))
      database_seqs <- read_fasta(database_path)
      tm_obj <- taxa::extract_tax_data(trimws(names(database_seqs)),
                                       include_match = FALSE,
                                       class_sep = ";",
                                       regex = "id=(.*)\\|name=(.*)\\|source=(.*)\\|tax_id=(.*)\\|taxonomy=(.*)$",
                                       key = c(id = "info", name = "info", source = "info", tax_id = "info", taxonomy = "class"))
      tm_obj$data$tax_data$sequence <- database_seqs
      return(tm_obj)
    })

  
  selected_subset <- eventReactive(
    eventExpr = input$search,
    ignoreNULL = TRUE,
    {
      # Parse whole database
      tm_obj <- selected_database()
      
      # Subset to taxon
      if (input$taxon_subset != "") {
        taxa_to_use <- trimws(strsplit(input$taxon_subset, split = ",+")[[1]])
        tm_obj <- filter_taxa(tm_obj, tolower(taxon_names) %in% tolower(taxa_to_use),
                              subtaxa = TRUE, supertaxa = TRUE, drop_obs = TRUE,
                              reassign_obs = FALSE)
      }
      
      
      return(tm_obj)
    })
  
  
  # makes the table of sequences displayed as search results
  output$sequence_table <- renderDataTable({
    if (is.null(selected_subset())) {
    } else {
      results <- selected_subset()
      
      classification <- results %>% 
        filter_taxa(taxon_names == "Oomycetes", subtaxa = TRUE) %>%
        filter_taxa(!is_root) %>%
        classifications() %>% 
        `[`(results$data$tax_data$taxon_id)
      
      family <- unlist(results$supertaxa(results$data$tax_data$taxon_id, recursive = FALSE, value = "taxon_names"))
      species <- taxon_names(results)[results$data$tax_data$taxon_id]
      
      
      results$data$tax_data %>% 
        transmute("Sequence ID" = id,
                  # "Organism name" = name,
                  # "NCBI taxon ID" = tax_id,
                  # "Taxonomy" = classification,
                  "Family" = family,
                  "Species" =  name)
    }
  }, selection = "multiple")
  
  
  public_release_data <- function() {
    get_release_data() %>%
      filter(public == TRUE) %>%
      arrange(-row_number())
  }
  
  # makes the table of sequences displayed as search results
  output$database_table <- renderDataTable({
    public_release_data() %>%
      transmute("Release" = release_number,
                "Date released" = release_date,
                "Notes" =  release_notes)
  }, selection = list(mode = 'single', selected = 1, target = 'row'))
  
  
  output$selected_seq_printout <- renderText({
    if(is.null(input$sequence_table_rows_selected)){}
    else{
      results <- selected_subset()
      clicked <- input$sequence_table_rows_selected
      output <- paste0(paste0(">", results$data$tax_data$input[clicked], "\n",
                       results$data$tax_data$sequence[clicked]), collapse = "\n")
      return(output)
    }
  })
  
  database_for_download_whole <- reactive({
      results <- selected_database()
      paste0(">", results$data$tax_data$input, "\n",
             results$data$tax_data$sequence)
  })
  
  database_for_download_subset <- eventReactive(
    eventExpr = input$search,
    ignoreNULL = TRUE,
    {
      results <- selected_subset()
      paste0(">", results$data$tax_data$input, "\n",
             results$data$tax_data$sequence)
    })
  
  
  database_for_download_selected <- eventReactive(
    eventExpr = input$sequence_table_rows_selected,
    ignoreNULL = TRUE,
    {
      results <- selected_subset()
      paste0(">", results$data$tax_data$input[input$sequence_table_rows_selected], "\n",
             results$data$tax_data$sequence[input$sequence_table_rows_selected])
    })
  
  output$download_data_whole <- downloadHandler(
    contentType = "text/csv",
    filename = function() {
      paste0("oomycetedb_whole_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".fa")
    },
    content = function(file) {
      not_used <- readr::write_lines(database_for_download_whole(), path = file)
    }
  )
  
  output$download_data_subset <- downloadHandler(
    contentType = "text/csv",
    filename = function() {
      paste0("oomycetedb_subset_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".fa")
    },
    content = function(file) {
      not_used <- readr::write_lines(database_for_download_subset(), path = file)
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
      if (is.null(input$sequence_table_rows_selected)) {
        return(list(
          # h4("Download database subset"),
          downloadButton(outputId = "download_data_whole", label = "Download database"),
          downloadButton(outputId = "download_data_subset", label = "Download database subset")
        ))
      } else {
        return(list(
          # h4("Download database subset"),
          downloadButton(outputId = "download_data_whole", label = "Download database"),
          downloadButton(outputId = "download_data_subset", label = "Download database subset"),
          downloadButton(outputId = "download_data_selected", label = "Download selected sequences")
        ))
        
      }
    }
  })
  
  output$database_table_ui <- renderUI({
    list(
      h3("Releases"),
      p('Click on one or more rows to select the release of the database to download or search'),
      DT::dataTableOutput("database_table")
    )
  })
  
  output$sequence_table_ui <- renderUI({
    if (is.null(selected_subset())) {
    } else {
      list(
        h3("Search Results"),
        p('Click on one or more rows to get the FASTA entry. Click again to deselect. You can also download the search results as a FASTA file by pressing the "Download database subset" button below.'),
        DT::dataTableOutput("sequence_table")
      )
    }
  })
  
  output$seq_details <- renderUI({
    x = 1
    if (is.null(selected_subset())) {
    } else {
      if (is.null(input$sequence_table_rows_selected)) {
        return(list(
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

