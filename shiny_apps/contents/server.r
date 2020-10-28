library(shiny)
library(shinyjs)
library(shinythemes)
library(metacoder)
library(ape)
library(here)
library(dplyr)
library(knitr)
library(kableExtra)



source(file.path(here(), "style.R"))
source(file.path(here(), "tools.R"))
source(file.path(here(), "configuration.R"))


server <- function(input, output, session) {
  output$heat_tree <- renderPlot({
    seqs <- read.FASTA(get_latest_release_fa())
    obj <- convert_oomydb_headers_to_taxmap(names(seqs))
    set.seed(2)
    obj %>%
      filter_taxa(taxon_names == "Oomycetes", subtaxa = TRUE) %>%
      remove_redundant_names() %>%
      heat_tree(node_label = taxon_names,
                node_size = n_obs, 
                node_color = n_obs,
                node_size_range = c(0.01, 0.04),
                node_label_size_range = c(0.015, 0.03),
                edge_size_range =  c(0.01, 0.04) * .5,
                layout = "da", initial_layout = "re",
                node_color_axis_label = "Number of sequences")
    
  })
  
  output$count_table <- renderTable({
    taxa_in_table <- c("Phytophthora", "Pythium", "Phytopythium", "Achlya", "Albugo", "Peronospora", "Plasmopara", "Saprolegnia")
    genus_table <- obj %>%
      filter_taxa(taxon_names %in% taxa_in_table, subtaxa = TRUE) %>%
      get_data_frame(c("n_obs", "taxon_names", "n_subtaxa")) %>%
      filter(taxon_names %in% taxa_in_table)
    
    
    genus_table <- rbind(genus_table, list(sum(genus_table$n_obs),
                                           "Total",
                                           sum(genus_table$n_subtaxa)))
    
    table_out <- genus_table %>%
      arrange(desc(as.numeric(n_subtaxa))) %>%
      transmute("Genus" = taxon_names, "Number of species" = as.integer(n_subtaxa), "Number of sequences" = as.integer(n_obs)) # %>%
      # kable() %>%
      # kable_styling(bootstrap_options = "striped", full_width = F, position = "center") %>% 
      # column_spec(1, italic = T) 
    
    table_out
  })
}

