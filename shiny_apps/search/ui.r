library(shinythemes)
library(DT)
library(shinyjs)


option_width <- "130px"

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
    div(style="display:inline-block",
        selectInput("db", "Database:", choices = c("current"), width = option_width)),
    div(style="display:inline-block",
        textInput("taxon_subset", "Taxa:", value = "", placeholder = "pythium, phytophthora", width = option_width)),
    br(),
    actionButton("search", "Search database")
  ),
  
  
  mainPanel(
    h3("Results"),
    DT::dataTableOutput("database_list"),
    uiOutput("download_data_form")
    # uiOutput("seq_list"),
    # uiOutput("seq_details")
  )
  
  
)
