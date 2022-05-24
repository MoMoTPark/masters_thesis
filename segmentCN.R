library(tidyverse)
library(testthat)
source("Simulation_functions.R")

#' Output: n segments with inferred CI width interval
segmentCN = function(n = 10000, 
                     read_length = 100, 
                     coverage = 30,
                     # Only 1st initial value
                     flank_seg_length,
                     # list of n values for n segments
                     cur_seg_length,
                     # Only 1st initial value
                     actual_flank_seg_cn,
                     # list of n values for n segments
                     actual_cur_seg_cn,
                     fragment_size = 600,
                     rd_model,
                     vaf_model) {
  
  inputDF = data.frame(n = n, 
                      read_length = read_length, 
                      coverage = coverage,
                      flank_seg_length = flank_seg_length,
                      cur_seg_length = cur_seg_length,
                      actual_flank_seg_cn = actual_flank_seg_cn,
                      actual_cur_seg_cn = actual_cur_seg_cn,
                      fragment_size = fragment_size,
                      rd_model = rd_model,
                      vaf_model = vaf_model)
  
  # Add the value of current segment from previous row to flanking segment of current row
  n = dim(inputDF)[1]
  for (i in 2:n) {
    inputDF$flank_seg_length[i] = inputDF$cur_seg_length[(i - 1)]
    inputDF$actual_flank_seg_cn[i] = inputDF$actual_cur_seg_cn[(i - 1)]
  }
  
  tempDF = inputDF %>%
    rowwise() %>%
    do({obs_cn_width(n, 
                                .$read_length, 
                                .$coverage,
                                .$flank_seg_length,
                                .$cur_seg_length,
                                .$actual_flank_seg_cn,
                                .$actual_cur_seg_cn,
                                .$fragment_size,
                                .$rd_model,
                                .$vaf_model)}) %>%
     ungroup() %>%
    fixDF(rd_model = paste("rd", rd_model), vaf_model = paste("vaf", vaf_model))
  
  data.frame(fragment_size = tempDF$fragment_size,
             read_length = tempDF$read_length,
             coverage = tempDF$coverage,
             flank_seg_length = tempDF$flank_seg_length,
             cur_seg_length = tempDF$cur_seg_length,
             actual_flank_seg_cn = tempDF$actual_flank_seg_cn,
             actual_cur_seg_cn = tempDF$actual_cur_seg_cn,
             rd_ci_width = tempDF$rd_ci_width,
             vaf_ci_width = tempDF$vaf_ci_width) #%>%
    #select(-c(fragment_size, read_length, coverage))
}


fixDF = function(df, rd_model, vaf_model) {
  vafDF = df %>%
    filter(model == vaf_model)
  
  df %>%
    filter(model == rd_model) %>%
    rename(rd_ci_width = ci_width) %>%
    bind_cols(vafDF$ci_width) %>%
    rename(vaf_ci_width = ...10) %>%
    select(-c(model))
}