library(shiny)
library(shinyjs)
library(shinythemes)

ui <- fluidPage(
  mainPanel(
    plotOutput(outputId = "heat_tree"),
    tableOutput(outputId = "count_table")
  )
  
)
