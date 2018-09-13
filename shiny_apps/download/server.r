require(XML)
library(plyr)
library(dplyr)
library(DT)
library(taxa)
library(shinyjs)

local_release_dir = "data/releases"
blast_database_dir = "data/blast_databases"

server <- function(input, output, session) {
  

  # makes the datatable
  output$database_list <- renderDataTable({
    readr::read_csv("../../data/releases.csv") %>%
      filter(public) %>%
      arrange(desc(release_number)) %>%
      transmute("Release" = paste("Release", release_number),
                "Date released" = release_date,
                "Notes" = release_notes)
  }, selection = "single")
  

  database_for_download_selected <- eventReactive(
    eventExpr = input$database_list_rows_selected,
    ignoreNULL = TRUE,
    {
      release_data <- readr::read_csv("../../data/releases.csv")
      selected_release_num <- release_data$release_number[input$database_list_rows_selected]
      selected_release_path <- file.path("..", "..", local_release_dir, paste0("release_", selected_release_num, ".fa"))
      results <- metacoder::read_fasta(selected_release_path)
      paste0(">", names(results), "\n", results)
    })
  
  
  output$download_data <- downloadHandler(
    contentType = "text/csv",
    filename = function() {
      release_data <- readr::read_csv("../../data/releases.csv")
      selected_release_num <- release_data$release_number[input$database_list_rows_selected]
      paste0("oomycetedb_release_", selected_release_num, "_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".fa")
    },
    content = function(file) {
      not_used <- readr::write_lines(database_for_download_selected(), path = file)
    }
  )
  

  output$download_data_form <- renderUI({
    if (is.null(input$database_list_rows_selected)) {
      return(list(
        p("Click on a release to download.")
      ))
    } else {
      return(list(
        # h4("Download database subset"),
        downloadButton(outputId = "download_data", label = "Download")
      ))
      
    }
  })
  
}

