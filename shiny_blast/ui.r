library(shinythemes)
library(DT)
library(shinyjs)


option_width <- "130px"
outfmt_options <- c(" 0: pairwise",
                    " 1: query-anchored showing identities", 
                    " 2: query-anchored no identities",
                    " 3: flat query-anchored show identities", 
                    " 4: flat query-anchored no identities",
                    " 5: XML Blast output", 
                    " 6: tabular",
                    " 7: tabular with comment lines",
                    " 8: Text ASN.1", 
                    " 9: Binary ASN.1",
                    "10: Comma-separated values",
                    "11: BLAST archive format (ASN.1)", 
                    "12: JSON Seqalign output",
                    "13: JSON Blast output",
                    "14: XML2 Blast output")

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
                  h3("Inputs"),
                  textAreaInput('query', 'Query sequences in FASTA format:', value = "", placeholder = "", width = "800px", height = "200px"),
                  div(style="display:inline-block",
                      selectInput("program", "Program:", choices = c("blastn", "tblastn"), width = option_width)),
                  div(style="display:inline-block",
                      selectInput("db", "Database:", choices = c("current"), width = option_width)),
                  div(style="display:inline-block",
                      selectInput("max_target_seqs", "Hits per query:", choices = c(1, 5, 10, 20, 50, 100), width = option_width, selected = 10)),
                  div(style="display:inline-block",
                      selectInput("eval", "Max e-value:", choices = c(1,0.001,1e-4,1e-5,1e-10), width = option_width)),
                  br(),
                  actionButton("blast", "Run BLAST")
                ),
                
                # this snippet generates a progress indicator for long BLASTs
                div(class = "busy",  
                    p("Calculation in progress.."), 
                    img(src="../images/blast_progress.gif", height = 100, width = 100, align = "center")
                ),
                
                # Basic results output
                mainPanel(
                  h3("Results"),
                  # h4("Download BLAST output"),
                  div(style="display:inline-block", uiOutput("download_data_form")),
                  div(style="display:inline-block", uiOutput("hit_list")),
                  div(style="display:inline-block", uiOutput("hit_details")),
                  uiOutput("before_blast_msg")
                )
)
