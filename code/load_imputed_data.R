library(tidyverse)

imputed_df <- read.csv('~/projects/metabolomics/data/metabolomics_converted/imputed-data-7-11-2016-samples.txt', sep = '\t')
keys_df <- data.frame(DIAG_code = imputed_df$DIAG_code, CLIENT_IDENTIFIER = imputed_df$CLIENT_IDENTIFIER)

rownames(imputed_df) <- imputed_df$DIAG_code
imputed_df$DIAG_code <- NULL
imputed_df$CLIENT_IDENTIFIER <- NULL

imputed_df <- drop_na(imputed_df)
imputed_df <- imputed_df %>% select_if(is.numeric)

zero_var <- function(v)any(sd(v) != 0)
imputed_df <- imputed_df %>% select_if(zero_var)
