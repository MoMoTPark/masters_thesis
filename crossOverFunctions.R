#' The following functions investigate the cross-over points between the two models


library(tidyverse)

#' @param dataframe Simulation dataframe that contains the results from both models
#' @param rd_model Chosen RD model from simulated data
#' @param vaf_model Chosen SV-VAF model from simulated data
#' Identifies crossing point between RD and SV-VAF model based on inferred CI width
#' interval and outputs a dataframe only containing datapoints with <= 0.05 CI
#' width difference between two models.   
crossOverDF = function (dataframe,
                       rd_model = c("rd start_in", "rd overlap"),
                       vaf_model = c("vaf split_read", "vaf read_pair")){
  
  rd_model = match.arg(rd_model)
  vaf_model = match.arg(vaf_model)
  
  vaf_df = dataframe %>%
    group_by(model) %>%
    filter(model == vaf_model) %>%
    rename(vaf_ci_width = ci_width) %>%
    ungroup()
  
  rd_df = dataframe %>%
    group_by(model) %>%
    filter(model == rd_model) %>%
    rename(rd_ci_width = ci_width) %>%
    ungroup()
  
  crossOver_df = bind_cols(rd_df, vaf_df$vaf_ci_width) %>%
    rename(vaf_ci_width = ...11) %>%
    #filter(abs(rd_ci_width - vaf_ci_width) <= 0.05)
    # To check for datapoints where performance of two models are significantly different
    filter(abs(rd_ci_width - vaf_ci_width) > 0.05)
  
  return(crossOver_df)
}


combineCI = function (dataframe,
                        rd_model = c("rd start_in", "rd overlap"),
                        vaf_model = c("vaf split_read", "vaf read_pair")){
  
  rd_model = match.arg(rd_model)
  vaf_model = match.arg(vaf_model)
  
  vaf_df = dataframe %>%
    filter(model == vaf_model) %>%
    rename(vaf_ci_width = ci_width) %>%
    ungroup()
  
  rd_df = dataframe %>%
    filter(model == rd_model) %>%
    rename(rd_ci_width = ci_width) %>%
    bind_cols(vaf_ci_width = vaf_df$vaf_ci_width) %>%
    ungroup() %>%
    select(-c(model))
    
    
    return(rd_df)
}