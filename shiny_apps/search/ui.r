library(shinythemes)
library(DT)
library(shinyjs)

local_release_dir = "data/releases"
blast_database_dir = "data/blast_databases"
option_width <- "130px"

get_blast_databases <- function(db_dir) {
  # Get possible blast database names
  db_file_names <- tools::file_path_sans_ext(list.files(db_dir))
  output <- unique(db_file_names)
  
  # Remove any that do not appear the correct number of times
  is_database <- vapply(output, function(x) as.list(table(db_file_names))[[x]] == 6, logical(1))
  output <- output[is_database]
  
  return(rev(output))
}



ui <- fluidPage(
  
  theme = shinytheme("cerulean"),
  
  tagList(
    tags$head(
      tags$link(rel="stylesheet", type="text/css", href="style.css"),
      tags$script(type="text/javascript", src = "busy.js")
    )
  ),
  
  includeCSS("www/style.css"),
  
  mainPanel(
    h3("Subsetting options"),
    selectInput("db", "Database:", choices = get_blast_databases(file.path("..", "..", blast_database_dir)), width = option_width),
    textInput("taxon_subset", "Taxa to subset the database to:", value = "",
              placeholder = "pythium, phytophthora, Albuginaceae, etc ...", width = "500px"),
    br(),
    actionButton("search", "Search database")
  ),
  
  
  mainPanel(
    h3("Results"),
    # DT::dataTableOutput("database_list"),
    uiOutput("seq_list_table"),
    uiOutput("download_data_form"),
    uiOutput("seq_details")
  )
  
  
)
