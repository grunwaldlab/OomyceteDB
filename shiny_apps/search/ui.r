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
    selectInput("db", "Database:", choices = c("current"), width = option_width),
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
