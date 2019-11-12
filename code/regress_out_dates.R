suppressWarnings(library(tidyverse))
suppressWarnings(library(argparse))
# source('~/projects/metabolomics/workflowr/code/load_imputed_data.R')
parser <- ArgumentParser()

parser$add_argument('-result_dir', default = 'data')
parser$add_argument('-result_prefix', default = 'imputed_data_date_controlled')
parser$add_argument('-input_samples')
parser$add_argument('-input_covariates')
h <- "Choose which covariates to include in the regression. Currently not implemented. Currently only PARAM_RUN_DAY is considered"
parser$add_argument('-covariate_types', action = 'append')
parser$add_argument('-inverse_normalize', action = 'store_true', default = FALSE)

args <- parser$parse_args()



main<-function(args){
  dir.create(args$result_dir, 
             recursive = TRUE,
             showWarnings = FALSE)
  samples_data <- load_samples(args$input_samples, args$inverse_normalize)
  cov_data <-load_covariates(args$input_covariates)
  resid_df <- data.frame(row.names = rownames(samples_data))
  coefficient_df <- NULL
  pval_df <- NULL
  for (metab_name in colnames(samples_data)){
    reg_obj <- run_single_regression(samples_data[[metab_name]], cov_data, metab_name)
    reg_output <- summary(reg_obj)$coefficients

    resid_df[[metab_name]] <- reg_obj$residuals
    if (!is.null(pval_df)){
      pval_df[[metab_name]] <- reg_output[,"Pr(>|t|)"]
      coefficient_df[[metab_name]] <-reg_output[,'Estimate']
    }else{
      pval_df <- data.frame(reg_output[,"Pr(>|t|)"])
      colnames(pval_df) <- c(metab_name)

      coefficient_df <- data.frame(reg_output[,'Estimate'])
      colnames(coefficient_df) <- c(metab_name)
    }
  }
  suffixes <- c('_coefficients.txt',
                '_pvals.txt',
                '_resids.txt')
  df_lst <- list(coefficient_df,
                 pval_df,
                 resid_df)
  for(i in 1:3){
    fp = file.path(args$result_dir, paste0(args$result_prefix, suffixes[i]))
    write.csv(df_lst[i], fp,
                quote = TRUE)
  }
  
  
}



load_samples <- function(fp, normalize){
  df <- read.csv(fp, sep = '\t')
  
  rownames(df) <- df$DIAG_code
  df$DIAG_code <- NULL
  df$CLIENT_IDENTIFIER <- NULL
  
  df <- drop_na(df)
  df <- df %>% select_if(is.numeric)
  
  zero_var <- function(v)any(sd(v) != 0)
  df <- df %>% select_if(zero_var)
  if(normalize){
    df <- data.frame(invnorm(df))
  }
  return(df)
}

invnorm <- function(x) {
  if(is.null(dim(x))) res = invnorm.vector(x) else
    res=apply(x,2,invnorm.vector)
  res
}

load_covariates<-function(fp){
  covariates <- read.csv('~/projects/metabolomics/data/metabolomics_converted/imputed-data-7-11-2016-sample-metadata.txt',
                         sep = '\t',
                         colClasses = c('factor'))
  rownames(covariates) <- covariates$DIAG_code
  return(data.frame(PARAM_RUN_DAY = covariates$PARAM_RUN_DAY))
}

run_single_regression<-function(response, cov_data, name){
  cov_data[[name]] <- response
  ff <- paste(name , '~ 0 + .')
  reg_obj  <- lm(ff , data = cov_data)
  return(reg_obj)
}


invnorm.vector = function(x) {yy = rank(x)/(length(x)+1); qnorm(yy)}


main(args)