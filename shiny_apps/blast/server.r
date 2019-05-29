library(XML)
library(shiny)
library(plyr)
library(dplyr)
library(DT)
library(taxa)
library(shinyjs)
library(stringr)
library(ggplot2)
library(readr)
library(ape)
library(shinythemes)
library(here)

options(shiny.sanitize.errors = FALSE)

source(file.path(here(), "tools.R"))

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


server <- function(input, output, session) {
  
  check_blast_input <- eventReactive(
    eventExpr = input$blast,
    ignoreNULL = TRUE,
    {
      if (input$query == "") {
        return("Please enter a query sequence in FASTA format to run BLAST.")
      }
      if (nchar(input$query) > 100000) {
        return("Query is too large for our servers. Please download the database and run BLAST on your computer instead.")
      }
      return("OK")
    }
  )
  
  observeEvent(
    eventExpr = input$blast,
    {
      check_results <- check_blast_input()
      if (check_results != "OK") {
        showNotification(check_results)
      }
      # validate(need(check_results == "OK", message = check_results))
    }
  )
  
  raw_blast_results <- eventReactive(
    eventExpr = input$blast,
    ignoreNULL = TRUE,
    {
      # Check input
      if (check_blast_input() != "OK") {
        return(NULL)
      }
      
      # Choose database 
      db <- file.path(blast_database_dir, input$db)
      remote <- c("")
      
      # Format query
      query_temp_path <- tempfile(fileext = ".fa")
      if (startsWith(input$query, ">")){
        writeLines(input$query, query_temp_path)
      } else {
        writeLines(paste0(">Query\n", input$query), query_temp_path)
      }
      
      # Calls BLAST
      blast_command <- paste(file.path(blast_path, input$program),
                             "-query", query_temp_path,
                             "-db", db ,
                             ifelse(input$program == "blastn", "-dust no", ""),
                             "-evalue", input$eval,
                             "-outfmt 11",
                             "-max_hsps 1",
                             "-max_target_seqs", input$max_target_seqs,
                             remote)
      system(blast_command, intern = TRUE)
    })
  
  
  blast_results <- reactive(
    {
      # Require raw blast results to be available
      req(raw_blast_results())
      
      # Convert to XML format
      raw_results <- raw_blast_results()
      tmp_file <- tempfile()
      write_lines(raw_results, path = tmp_file)
      blast_formatter_cmd <- paste(file.path(blast_path, "blast_formatter"),
                                   "-archive", tmp_file,
                                   "-outfmt 5")
      xmlParse(system(blast_formatter_cmd, intern = TRUE))
    })
  
  
  blast_results_for_download <- eventReactive(
    eventExpr = input$outfmt,
    ignoreNULL = TRUE,
    {
      # Convert to XML format
      raw_results <- raw_blast_results()
      tmp_file <- tempfile()
      write_lines(raw_results, path = tmp_file)
      blast_formatter_cmd <- paste(file.path(blast_path, "blast_formatter"),
                                   "-archive", tmp_file,
                                   "-outfmt", str_match(input$outfmt, " *([0-9]+):.+")[, 2])
      system(blast_formatter_cmd, intern = TRUE)
    })
  
  # Now to parse the results...
  parsed_results <- reactive({
    # Check that blast results are available
    my_blast_results <- req(blast_results())
    
    # the first chunk is for multi-fastas
    results <- xpathApply(my_blast_results, '//Iteration', function(row) {
      
      get_hit_info <- function(info) {
        getNodeSet(row, paste0('Iteration_hits//Hit//Hit_hsps//Hsp//', info)) %>% sapply(., xmlValue)
      }
      
      query_id <-   getNodeSet(row, 'Iteration_query-def') %>% sapply(., xmlValue)
      query_len <-  getNodeSet(row, 'Iteration_query-len') %>% sapply(., xmlValue)
      hit_ids <-    getNodeSet(row, 'Iteration_hits//Hit//Hit_id') %>% sapply(., xmlValue)
      hit_length <- getNodeSet(row, 'Iteration_hits//Hit//Hit_len') %>% sapply(., xmlValue)
      bit_score <-  get_hit_info('Hsp_bit-score')
      score <-      get_hit_info('Hsp_score') 
      evalue <-     get_hit_info('Hsp_evalue')
      query_from <- get_hit_info('Hsp_query-from')
      query_to <-   get_hit_info('Hsp_query-to')
      hit_from <-   get_hit_info('Hsp_hit-from')
      hit_to <-     get_hit_info('Hsp_hit-to')
      identity <-   get_hit_info('Hsp_identity')
      positive <-   get_hit_info('Hsp_positive')
      gaps <-       get_hit_info('Hsp_gaps')
      align_len <-  get_hit_info('Hsp_align-len')
      query_seq <-  get_hit_info('Hsp_qseq')
      hit_seq <-    get_hit_info('Hsp_hseq')
      midline <-    get_hit_info('Hsp_midline')
      
      cbind(query_id, query_len, hit_ids, hit_length, bit_score, score, evalue,
            query_from, query_to, hit_from, hit_to, identity, positive, gaps,
            align_len, query_seq, hit_seq, midline)
    })
    
    # this ensures that NAs get added for no hits
    results <- rbind.fill(lapply(results, function(y) {as.data.frame(y, stringsAsFactors = FALSE)}))
    
    # If there are no hits, make an empty table
    if (ncol(results) == 2) {
      results <- data.frame(query_id = character(),
                            query_len = character(),
                            hit_ids = character(),
                            hit_length = character(),
                            bit_score = character(),
                            score = character(),
                            evalue = character(),
                            query_from = character(),
                            query_to = character(),
                            hit_from = character(),
                            hit_to = character(),
                            identity = character(),
                            positive = character(),
                            gaps = character(),
                            align_len = character(),
                            query_seq = character(),
                            hit_seq = character(),
                            midline = character())
    }
    
    # Convert numeric cols to numeric
    numeric_cols <- c("query_len", "hit_length", "bit_score", "score", "evalue", "query_from", 
                      "query_to", "hit_from", "hit_to", "identity", "positive", "gaps", 
                      "align_len")
    results[numeric_cols] <- lapply(results[numeric_cols], as.numeric)
    
    # Replace database index with header
    database_path <- file.path(local_release_dir, paste0(isolate(input$db), ".fa"))
    database_seqs <- read.FASTA(database_path)
    results$hit_ids <- names(database_seqs)[as.numeric(results$hit_ids)]
    
    # Calculate derived columns
    results$hit_index <- seq_len(nrow(results))
    results$prop_identity <- results$identity / results$align_len
    results$prop_match_len <- (results$query_to - results$query_from + 1) / results$query_len
    
    if (nrow(results) > 0) {
      # Convert to taxmap
      tm_obj <- extract_tax_data(results$hit_ids,
                                 class_sep = ";",
                                 regex = "id=(.+)\\|name=(.+)\\|source=(.+)\\|tax_id=(.+)\\|taxonomy=(.+)$",
                                 key = c(id = "info", name = "info", source = "info", tax_id = "info", taxonomy = "class"))
      
      # Add on results table to taxmap
      tm_obj$data$tax_data <- bind_cols(tm_obj$data$tax_data, results)
    } else {
      # Convert to taxmap
      tm_obj <- taxmap()
      
      # Add on results table to taxmap
      tm_obj$data$tax_data <- as_tibble(bind_cols(data.frame(taxon_id = character(),
                                                             id = character(),
                                                             name = character(),
                                                             source = character(),
                                                             tax_id = character(),
                                                             taxonomy = character(),
                                                             input = character()),
                                                  results))
    }
    
    return(tm_obj)
    
  })
  
  selected_query_hit_table <- reactive({
    # Check that blast results are available
    my_parsed_results <- req(parsed_results())
    hit_data <- my_parsed_results$data$tax_data
    
    # Get the hit id that is clicked
    clicked_index <- input$blast_results_rows_selected
    hit_selected <- hit_data$hit_ids[clicked_index]
    
    # Subset data for query 
    query_selected <- hit_data$query_id[clicked_index]
    hit_data <- filter(hit_data, query_id == query_selected)
    
    return(hit_data)
  })
  
  # Make hit plot
  output$hit_plot <- renderPlot({
    # Check that blast results are available
    my_parsed_results <- req(parsed_results())
    hit_data <- my_parsed_results$data$tax_data
    
    # Get the hit id that is clicked
    clicked_index <- input$blast_results_rows_selected
    hit_selected <- hit_data$hit_ids[clicked_index]
    
    # Subset data for query 
    query_selected <- hit_data$query_id[clicked_index]
    hit_data <- filter(hit_data, query_id == query_selected)
    
    # Show a maximum number
    max_shown <- 20
    if (nrow(hit_data) > max_shown) {
      hit_data <- hit_data[1:max_shown, ]
      plot_title <- paste0("Distribution of the top ", max_shown, " Blast Hits")
    } else {
      plot_title <- paste0("Distribution of ", nrow(hit_data), " Blast Hits")
    }
    
    # Prepare data for plotting
    hit_data$y <- - seq_len(nrow(hit_data))
    hit_data$color_interval <- cut(hit_data$score,
                                   breaks = c(-99999, 40, 50, 80, 200, 99999),
                                   labels = c("<40", "40-50", "50-80", "80-200", ">=200"))
    hit_data <- mutate(hit_data,
                       start = query_from - hit_from + 1,
                       end = query_to + (hit_length - hit_to))
    
    # Get data for mismatches
    mismatches <- lapply(str_locate_all(hit_data$midline, " +"), as.data.frame)
    align_data <- do.call(bind_rows, mismatches)
    align_data$y <- rep(hit_data$y, sapply(mismatches, nrow))
    gap_correction <- rep((hit_data$align_len - hit_data$gaps) / hit_data$align_len, sapply(mismatches, nrow))
    align_data$start <- align_data$start * gap_correction + rep(hit_data$query_from, sapply(mismatches, nrow)) - 1
    align_data$end <- align_data$end * gap_correction + rep(hit_data$query_from, sapply(mismatches, nrow)) - 1
    
    # Get data for highlight
    x_min <- min(hit_data$query_from)
    x_max <- max(hit_data$query_to)
    highlight_data <- data.frame(x = x_min, xend = x_max, y = - which(hit_data$hit_ids == hit_selected))
    
    # Get data for query
    query_data <- data.frame(x = 1, xend = hit_data$query_len[1], y = 0)
    
    # plot
    ggplot() +
      geom_segment(data = highlight_data,
                   aes(x = x, xend = xend, y = y, yend = y),
                   colour = "yellow", size = 5) +
      geom_segment(data = hit_data,
                   aes(x = start, xend = end, y = y, yend = y, colour = color_interval),
                   size = 2) +
      scale_colour_manual(values = c("#000000", "#0000ff", "#00ff00", "#ff00ff", "#ff0000"),
                          labels = c("<40", "40-50", "50-80", "80-200", ">=200"), 
                          drop = FALSE) +
      geom_segment(data = query_data,
                   aes(x = x, xend = xend, y = y, yend = y),
                   colour = "#58c7c7", size = 5) +
      geom_segment(data = align_data,
                   aes(x = start - 0.5, xend = end + 0.5,
                       y = y, yend = y), color = "black", size = 4) + 
      guides(colour = guide_legend(title = "Alignment scores:")) +
      scale_x_continuous(position = "top") + 
      ggtitle(plot_title) +
      theme(legend.position = "top",
            axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_blank(),
            panel.background = element_blank(), 
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 12))
  })
  
  # makes the datatable
  output$blast_results <- renderDataTable({
    # Check that blast results are available
    my_parsed_results <- req(parsed_results())
    if (nrow(my_parsed_results$data$tax_data) == 0) {
      return(NULL)
    }
    
    # Filter taxonomic classifications for better printing 
    results <- parsed_results()
    classification <- results %>% 
      filter_taxa(taxon_names == "Oomycetes", subtaxa = TRUE) %>%
      filter_taxa(!is_root) %>%
      classifications() %>% 
      `[`(results$data$tax_data$taxon_id)
    
    genus <- unlist(results$supertaxa(results$data$tax_data$taxon_id, recursive = FALSE, value = "taxon_names"))
    species <- taxon_names(results)[results$data$tax_data$taxon_id]
    
    results$data$tax_data %>% 
      arrange(desc(prop_identity)) %>%
      transmute("Query ID" = stringr::str_trunc(query_id, 25),
                # "Hit taxonomic classification" = classification,
                "Genus" = genus,
                "Species" =  name,
                "Identity (%)" = round(prop_identity * 100, digits = 3),
                "Query Coverage (%)" = round(prop_match_len * 100, digits = 3),
                "E value" = evalue)
    
  }, selection = "single")
  
  # this chunk gets the alignemnt information from a clicked row
  output$clicked <- renderTable({
    # Check that a row has been selected
    req(input$blast_results_rows_selected)
    
    # Create table of output information to display
    clicked <- input$blast_results_rows_selected
    results <- parsed_results()
    tableout <- data.frame(results$data$tax_data[clicked,])
    tableout <- transmute(tableout,
                          "Hit taxonomic classification" = taxonomy,
                          "Hit NCBI taxon ID" = tax_id,
                          # accesion number when we have it
                          "Identity (%)" = round(prop_identity * 100, digits = 3),
                          "Alignment length" = align_len,
                          "Mismatches" = align_len - identity, #check 
                          "Gaps" =  gaps,
                          "Query aligned range" = paste(query_from, "-", query_to),
                          "Hit aligned range" = paste(hit_from, "-", hit_to),
                          # "E value" = evalue,
                          "Bit score"  = bit_score,
                          "Query ID" = query_id,
                          # "Hit ID" = hit_ids,
                          "Query Coverage (%)" = round(prop_match_len * 100, digits = 3),
                          "Matching base pairs" = identity,
                          "Hit length" = hit_length
                          # "Hit sequence" = hit_seq
    )
    
    tableout <- t(tableout)
    names(tableout) <- c("")
    colnames(tableout) <- NULL
    data.frame(tableout)
    
  },rownames = TRUE,colnames = FALSE)
  
  # make alignments for clicked rows
  output$alignment <- renderText({
    # Check that a row has been selected
    req(input$blast_results_rows_selected)
    
    clicked = input$blast_results_rows_selected
    
    # loop over the xml to get the alignments
    align <- xpathApply(blast_results(), '//Iteration',function(row){
      top <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq') %>% sapply(., xmlValue)
      mid <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_midline') %>% sapply(., xmlValue)
      bottom <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq') %>% sapply(., xmlValue)
      rbind(top,mid,bottom)
    })
    
    # split the alignments every 80 carachters to get a "wrapped look"
    alignx <- do.call("cbind", align)
    splits <- strsplit(gsub("(.{80})", "\\1,", alignx[1:3,clicked]),",")
    
    # paste them together with returns '\n' on the breaks
    split_out <- lapply(seq_along(splits[[1]]), function(i){
      rbind(paste0("Q-", splits[[1]][i], "\n"),
            paste0("  ", splits[[2]][i], "\n"),
            paste0("H-", splits[[3]][i], "\n"),
            "\n")
    })
    output <- unlist(split_out)
    output[1] <- paste0(" ", output[1])
    return(output)
  })
  
  # Downloadable csv of selected dataset ----
  output$download_data <- downloadHandler(
    contentType = "text/csv",
    filename = function() {
      paste0("blast_results_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".txt")
    },
    content = function(file) {
      print("heres")
      not_used <- write_lines(blast_results_for_download(), path = file)
    }
  )
  
  output$download_data_form <- renderUI({
    # Check that blast results are available
    req(blast_results())
    if (nrow(parsed_results()$data$tax_data) == 0) {
      return(NULL)
    }
    
    list(
      selectInput("outfmt", "Output format:", choices = outfmt_options, width = "300px",
                  selected = "10: Comma-separated values"),
      downloadButton(outputId = "download_data", label = "Download Results")
    )
    
  })
  
  output$hit_list <- renderUI({
    # Check that blast results are available
    req(blast_results())
    if (nrow(parsed_results()$data$tax_data) == 0) {
      return(NULL)
    }    
    list(
      h4("Hits"),
      dataTableOutput("blast_results")
    )
    
  })
  
  output$no_hits <- renderText({
    # Check that blast results are available
    req(blast_results())
    
    if (nrow(parsed_results()$data$tax_data) == 0) {
      return('No hits found.')
    } else {
      return(NULL)
    }
  })
  
  output$hit_details <- renderUI({
    # Check that blast results are available
    req(blast_results())
    if (nrow(parsed_results()$data$tax_data) == 0) {
      return(NULL)
    }
    
    if (is.null(input$blast_results_rows_selected)) {
      return(list(
        h4("Hit details"),
        p("Click on a row for details...")
      ))
    }
    return(list(
      h4("Hit details"),
      plotOutput("hit_plot", click = "clicked_position"),
      p(tableOutput("clicked")),
      verbatimTextOutput("alignment")
    ))
    
  })
  
  observe({
    req(input$clicked_position)
    
    # Get index of hit that was clicked on
    pos_data <- input$clicked_position
    row_index <- abs(round(pos_data$y))
    hit_index <- selected_query_hit_table()$hit_index[row_index]
    
    # Select the row in the hit table
    if (! is.null(hit_index) && hit_index %in% seq_len(nrow(parsed_results()$data$tax_data))) {
      dt_proxy <- dataTableProxy("blast_results")
      selectRows(dt_proxy, hit_index)
    }
  })
  
}

