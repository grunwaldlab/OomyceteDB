# library(shinythemes)
library(DT)
library(shinyjs)
library(here)

source(file.path(here(), "tools.R"))

option_width <- "130px"
options(shiny.sanitize.errors = FALSE)

ui <- fluidPage(
  
  # theme = shinytheme("cerulean"),
  
  tagList(
    tags$head(
      tags$link(rel="stylesheet", type="text/css", href="style.css"),
      tags$script(type="text/javascript", src = "busy.js")
    )
  ),
  
  includeCSS("www/style.css"),
  
  mainPanel(
    h3("Releases"),
    p('Click on one or more rows to select the release of the database to download or search'),
    uiOutput("database_table_ui"),
    uiOutput("database_whole_download_ui"),
    uiOutput("search_ui"),
    uiOutput("sequence_table_ui"),
    uiOutput("download_data_form"),
    uiOutput("seq_details")
  )
  
  
)
