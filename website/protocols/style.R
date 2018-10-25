library(knitr)
library(dplyr)
library(readr)
library(kableExtra)

knitr::opts_chunk$set(echo = FALSE, message = FALSE)


custom_kable <- function(data, ...) {
  kable(data, booktabs = TRUE, linesep = "") %>%
    kable_styling(bootstrap_options = c("striped"),
                  latex_options = c("striped"),
                  full_width = FALSE,
                  position = "left") %>%
    row_spec(0, bold = TRUE)
}