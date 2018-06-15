library(shinythemes)
library(DT)

ui <- fluidPage(theme = shinytheme("cerulean"),
                tagList(
                  tags$head(
                    tags$link(rel="stylesheet", type="text/css", href="style.css"),
                    tags$script(type="text/javascript", src = "busy.js")
                  )
                ),
                
                includeCSS("www/style.css"),
                
                #This block gives us all the inputs:
                mainPanel(
                  textAreaInput('query', 'Input sequence:', value = "", placeholder = "", width = "800px", height="200px"),
                  selectInput("db", "Database:", choices=c("current"), width="120px"),
                  div(style="display:inline-block",
                      selectInput("program", "Program:", choices=c("blastn","tblastn"), width="100px")),
                  div(style="display:inline-block",
                      selectInput("eval", "e-value:", choices=c(1,0.001,1e-4,1e-5,1e-10), width="120px")),
                  actionButton("blast", "BLAST!")
                ),
                
                #this snippet generates a progress indicator for long BLASTs
                div(class = "busy",  
                    p("Calculation in progress.."), 
                    img(src="../images/blast_progress.gif", height = 100, width = 100, align = "center")
                ),
                
                #Basic results output
                mainPanel(
                  h4("Results"),
                  DT::dataTableOutput("blast_results"),
                  p("Alignment:", tableOutput("clicked") ),
                  verbatimTextOutput("alignment")
                )
)
