library(tidyverse)
library(testthat)

#' The ultimate goal of the following functions is to calculate the copy number 
#' of current segment using either SV VAF (variant allele frequency at SV breakpoint) 
#' or read depth (RD). Two models, RD-based `rd_model` and SV-VAF based `vaf_model`,
#' are used to investigate the effect of current segment's length on copy number
#' inference. Each model takes positional arguments which specifies the method
#' that is used for calculating the relevant signal in order to infer copy number.
#' Sequencing reads are generated using a Poisson distribution for the purposes of
#' simulations.  


#' Calculates copy number using SV-VAF
#' In a given bin with the presence of SV breakpoints, SV-VAF value of each 
#' SV breakpoint can be used for copy number inference.
#' In the context of this function, flanking segment refers to the larger segment 
#' in length where a bin has two segments separated by at SV breakpoints.    
#' @param obs_sup_reads The number of sequencing reads which support SV at an 
#' SV breakpoint
#' @param obs_ref_reads The number of sequencing reads that support the reference 
#' genome for a given segment
#' @param read_length The length of sequencing reads
#' @param flank_seg_reads The number of sequence reads which are in the flanking 
#' segment of a given window (read count in the flanking segment)
#' @param flank_seg_length The length of flanking segment
#' @param coverage Coverage of sequencing reads (NOT sequencing coverage)
#' @param orientation Orientation of SV in relation to its direction within the
#' genome. This parameter specifies whether there is an increase or decrease in 
#' copy number when moving in the direction of flanking to current segment.
#' Takes "increasing" or "decreasing" with reference to whether the copy number 
#' (equivalent to read count if copy number is inferred by read depth) from a 
#' flanking segment is increased when compared to the current segment and "decreasing" 
#' for when the copy number of flanking segment is higher than the current segment (SV).
#' @param rd_model Model that utilises RD signal for copy number inference.
#' "start_in" option specifies counting only reads that start in current segment. 
#' "overlap" option specifies counting reads which either start-in or overlap the 
#' current segment
cn_from_vaf = function(obs_sup_reads, 
                       obs_ref_reads, 
                       read_length, 
                       flank_seg_reads, 
                       flank_seg_length,
                       coverage,
                       orientation = c("increasing", "decreasing"),
                       rd_model = c("start_in", "overlap")){
  orientation = match.arg(orientation)
  rd_model = match.arg(rd_model)
  flank_seg_cn = cn_from_rd(obs_reads = flank_seg_reads, 
                            seg_length = flank_seg_length, 
                            read_length = read_length, 
                            coverage = coverage,
                            rd_model = rd_model)
  cn_from_vaf_with_flank_seg_cn(obs_sup_reads, 
                                obs_ref_reads, 
                                flank_seg_cn, 
                                orientation)
}


#' Calculate copy number of current segment with a given flanking segment copy number
#' @param obs_sup_reads The number of sequencing reads which support SV at an 
#' SV breakpoint
#' @param obs_ref_reads The number of sequencing reads that support the reference 
#' genome for a given segment
#' @param flank_seg_cn The copy number of flanking segment calculated from read depth (RD)
#' by `rd_model`
#' @param orientation Orientation of SV in relation to its direction within the
#' genome. This parameter specifies whether there is an increase or decrease in 
#' copy number when moving in the direction of flanking to current segment.
#' Takes "increasing" or "decreasing" with reference to whether the copy number 
#' (equivalent to read count if copy number is inferred by read depth) from a 
#' flanking segment is increased when compared to the current segment and "decreasing" 
#' for when the copy number of flanking segment is higher than the current segment (SV).
cn_from_vaf_with_flank_seg_cn = function(obs_sup_reads, 
                                         obs_ref_reads, 
                                         flank_seg_cn, 
                                         orientation = c("increasing", "decreasing")) {
  input_arg = match.arg(orientation)
  obs_total_reads = obs_sup_reads + obs_ref_reads
  sv_vaf = obs_sup_reads / obs_total_reads
  if (input_arg == "increasing") {
    cur_seg_cn = flank_seg_cn / (1 - sv_vaf) 
  } else {
    cur_seg_cn = flank_seg_cn * (1 - sv_vaf)
  }
  cur_seg_cn
}

#' `rd_model`: calculates copy number of a segment using read depth (RD) 
#' @param obs_reads Number of reads at a given base position
#' @param seg_length The length of a segment for copy number inference
#' @param read_length The length of sequencing reads 
#' @param coverage Coverage of sequencing reads (NOT sequencing coverage)
#' @param rd_model Model that utilises RD signal for copy number inference.
#' "start_in" option specifies counting only reads that start in current segment. 
#' "overlap" option specifies counting reads which either start-in or overlap the 
#' current segment
cn_from_rd = function(obs_reads, 
                      seg_length, 
                      read_length, 
                      coverage,
                      rd_model = c("start_in", "overlap")){
  rd_arg = match.arg(rd_model)
  eff_length = ifelse(rd_arg == "overlap", read_length + seg_length, seg_length)
  obs_bases = obs_reads * read_length
  obs_bases / (coverage * eff_length)
}


#' Calculates the expected number of read counts in a segment for `rd_model`
#' @param cn Copy number of the current segment
#' @param coverage Coverage of sequencing reads (NOT sequencing coverage)
#' @param read_length The length of sequencing reads 
#' @param seg_length The length of current segment
expected_segment_reads = function(cn, 
                          coverage, 
                          read_length, 
                          seg_length,
                          rd_model = c("start_in", "overlap")){
  model = match.arg(rd_model)
  ifelse(model == "overlap",
    expected_reads(cn, coverage, read_length, read_length + seg_length),
    expected_reads(cn, coverage, read_length, seg_length))
}
#' Calculates the expected number of read counts in a segment for `sv_model`
#' @param cn Copy number of the current segment
#' @param coverage Coverage of sequencing reads (NOT sequencing coverage)
#' @param read_length The length of sequencing reads 
#' @param fragment_size Fragment size of sequencing reads (paired-end sequencing
#' fragment size: 2 x `read_length` + inner-distance between reads). Fragment 
#' size is set explicitely to prevent any ambiguities in fragment size calculation.
#' @param vaf_model copy number inference model which uses variant allele frequency
#' of SV breakpoint to infer copy number of a segment. "read_pair" argument specifies
#' using fragment size as the effective length of an SV, instead of just considering
#' split reads which support an SV. "split_read" argument specified that the length
#' of an SV is only the read length, namely SV is only supported by split reads.
expected_sv_reads = function(cn, 
                          coverage, 
                          read_length, 
                          fragment_size,
                          vaf_model = c("read_pair", "split_read")){
  vaf_model = match.arg(vaf_model)
  ifelse(vaf_model == "read_pair",
    expected_reads(cn, coverage, read_length, fragment_size),
    expected_reads(cn, coverage, read_length, read_length))
}
#' Calculates the expected number of read counts by considering the SV length as
#' `effective_interval_length`
#' @param cn Copy number of a segment
#' @param coverage Coverage of sequencing reads (NOT sequencing coverage)
#' @param read_length The length of sequencing reads 
#' @param effective_interval_length The effective length of the interval over which 
#' the read can start in (in most cases simply SV length)
expected_reads = function(cn, 
                          coverage,
                          read_length, 
                          effective_interval_length) {
  cn * coverage * effective_interval_length / read_length
}


#' Produces inferred copy number values for specified input data and outputs a dataframe 
#' containing copy numbers calculated from both SV-VAF and read depth
#' @param n The number of simulations to perform
#' @param coverage  Coverage of sequencing reads (NOT sequencing coverage)
#' @param read_length The length of sequencing reads 
#' @param flank_seg_length The length of flanking segment
#' @param cur_seg_length The length of current segment (to infer copy number for this segment)
#' @param actual_flank_seg_cn Copy number of flanking segment (a known value for simulation)
#' @param actual_cur_seg_cn Copy number of current segment (a known value for simulation)
#' @param rd_model Model that utilises RD signal for copy number inference.
#' "start_in" option specifies counting only reads that start in current segment. 
#' "overlap" option specifies counting reads which either start-in or overlap the 
#' current segment
#' @param vaf_model copy number inference model which uses variant allele frequency
#' of SV breakpoint to infer copy number of a segment. "read_pair" argument specifies
#' using fragment size as the effective length of an SV, instead of just considering
#' split reads which support an SV. "split_read" argument specified that the length
#' of an SV is only the read length, namely SV is only supported by split reads.    
obs_cn = function(n, 
                  read_length, 
                  coverage,
                  flank_seg_length,
                  cur_seg_length,
                  actual_flank_seg_cn,
                  actual_cur_seg_cn,
                  fragment_size,
                  rd_model = c("start_in", "overlap"),
                  vaf_model = c("read_pair", "split_read")) {
  rd_model = match.arg(rd_model)
  vaf_model = match.arg(vaf_model)
  orientation = ifelse(actual_flank_seg_cn < actual_cur_seg_cn, "increasing", "decreasing")
  sv_cn = abs(actual_flank_seg_cn - actual_cur_seg_cn) # Copy number of SV 
  ref_cn = min(actual_flank_seg_cn, actual_cur_seg_cn) # Copy number of Reference 
  # The number of reads supporting the SV at the SV breakpoint 
  obs_sup_read_count = expected_sv_reads(sv_cn, coverage, read_length, fragment_size, vaf_model)
  # The number of reads supporting the reference genome at the SV breakpoint
  obs_ref_read_count = expected_sv_reads(ref_cn, coverage, read_length, fragment_size, vaf_model)
  # The number of reads observed along the flanking segment
  flank_seg_read_count = expected_segment_reads(actual_flank_seg_cn, coverage, read_length, flank_seg_length, rd_model)
  cur_seg_read_count = expected_segment_reads(actual_cur_seg_cn, coverage, read_length, cur_seg_length, rd_model)
  obs_sup_reads = rpois(n, obs_sup_read_count)
  obs_ref_reads = rpois(n, obs_ref_read_count)
  flank_seg_reads = rpois(n, flank_seg_read_count)
  cur_seg_reads = rpois(n, cur_seg_read_count)
  vaf_cn_estimates = cn_from_vaf(obs_sup_reads, 
              obs_ref_reads, 
              read_length, 
              flank_seg_reads, 
              flank_seg_length,
              coverage,
              orientation,
              rd_model)
  rd_cn_estimates = cn_from_rd(cur_seg_reads, cur_seg_length, read_length, coverage, rd_model)
  bind_rows(
    data.frame(
      model=paste("vaf", vaf_model),
      cn_estimate=vaf_cn_estimates),
    data.frame(
      model=paste("rd", rd_model),
      cn_estimate=rd_cn_estimates)) %>%
    mutate(
      read_length=read_length,
      fragment_size=fragment_size,
      coverage=coverage,
      flank_seg_length=flank_seg_length,
      cur_seg_length=cur_seg_length,
      actual_flank_seg_cn=actual_flank_seg_cn,
      actual_cur_seg_cn=actual_cur_seg_cn)
}


#' Calculates confidence interval width of inferred copy numbers for simulated
#' reads for given segment lengths  
#' @param n The number of simulations to perform
#' @param read_length The length of sequencing reads 
#' @param coverage  Coverage of sequencing reads (NOT sequencing coverage)
#' @param flank_seg_length The length of flanking segment
#' @param cur_seg_length The length of current segment (to infer copy number for this segment)
#' @param actual_flank_seg_cn Copy number of flanking segment (a known value for simulation)
#' @param actual_cur_seg_cn Copy number of current segment (a known value for simulation)
#' @param fragment_size Fragment size of sequencing reads (paired-end sequencing
#' fragment size: 2 x `read_length` + inner-distance between reads). Fragment 
#' size is set explicitly to prevent any ambiguities in fragment size calculation.
#' @param rd_model Model that utilises RD signal for copy number inference.
#' "start_in" option specifies counting only reads that start in current segment. 
#' "overlap" option specifies counting reads which either start-in or overlap the 
#' current segment
#' @param vaf_model copy number inference model which uses variant allele frequency
#' of SV breakpoint to infer copy number of a segment. "read_pair" argument specifies
#' using fragment size as the effective length of an SV, instead of just considering
#' split reads which support an SV. "split_read" argument specified that the length
#' of an SV is only the read length, namely SV is only supported by split reads. 
obs_cn_width = function(n, 
                        read_length, 
                        coverage,
                        flank_seg_length,
                        cur_seg_length,
                        actual_flank_seg_cn,
                        actual_cur_seg_cn,
                        fragment_size,
                        rd_model = c("start_in", "overlap"),
                        vaf_model = c("read_pair", "split_read")){
  rd_model = match.arg(rd_model)
  vaf_model = match.arg(vaf_model)
  obs_cn(n, 
         read_length, 
         coverage,
         flank_seg_length,
         cur_seg_length,
         actual_flank_seg_cn,
         actual_cur_seg_cn,
         fragment_size,
         rd_model = rd_model,
         vaf_model = vaf_model) %>%
    group_by(model, 
             fragment_size, 
             read_length, 
             coverage, 
             cur_seg_length, 
             flank_seg_length, 
             actual_flank_seg_cn, 
             actual_cur_seg_cn) %>%
    summarise(ci_width=estimate_ci_width(cn_estimate)) %>%
    ungroup()
}


#' Calculates the confidence interval of inferred copy number
#' @param observed A range of inferred copy numbers 
estimate_ci_width = function(observed, interval=0.95) {
  outlier_size = (1 - interval) / 2
  x = sort(observed)
  x[round(length(observed) * (1 - outlier_size))] - x[pmax(1, round(length(observed) * outlier_size))]
}


#' Calculate SV-VAF from given actual flanking and current segment copy number values
#' @param actual_cur_seg_cn Copy number of the current segment 
#' @param actual_flank_seg_cn Copy number of the flanking segment
sv_vaf = function(actual_flank_seg_cn, actual_cur_seg_cn) {
  ifelse(actual_flank_seg_cn < actual_cur_seg_cn,
         1 - (actual_flank_seg_cn / actual_cur_seg_cn),
         1 - (actual_cur_seg_cn / actual_flank_seg_cn))
}


#' Create data frames with specified variable values to perform simulation over a range of 
#' inputs. Calculated confidence interval of inferred copy numbers by `rd_model`
#' and `sv_vaf` are used as a metric for performance comparison of the two models.
#' The focus of permutated inputs in this function would be to primarily investigate
#' the effect of current segment length on inferred copy numbers using various 
#' possible actual (known) copy numbers of flanking and current segments. 
sim_permutations = function(
    read_length = 100,
    coverage = 30,
    cur_seg_length,
    flank_seg_length = 1000000,
    actual_flank_seg_cn,
    actual_cur_seg_cn,
    fragment_size,
    n=10000,
    rd_model,
    vaf_model 
    ){
  expand.grid(
    read_length = read_length,
    coverage = coverage,
    cur_seg_length = cur_seg_length,
    flank_seg_length = flank_seg_length,
    actual_flank_seg_cn = actual_flank_seg_cn,
    actual_cur_seg_cn = actual_cur_seg_cn,
    fragment_size = fragment_size,
    rd_model = rd_model,
    vaf_model = vaf_model, 
    stringsAsFactors = FALSE) %>%
  filter(actual_cur_seg_cn != actual_flank_seg_cn) %>%
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
  mutate(vaf=sv_vaf(actual_flank_seg_cn, actual_cur_seg_cn))
}





