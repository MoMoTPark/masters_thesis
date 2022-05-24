source("Simulation_functions.R")


test_that("cn_from_vaf_with_flank_seg_cn", {
  expect_equal(2, cn_from_vaf_with_flank_seg_cn(obs_sup_reads=0, 
                                                obs_ref_reads=100, 
                                                flank_seg_cn=2, 
                                                orientation="increasing"))
  expect_equal(2, cn_from_vaf_with_flank_seg_cn(obs_sup_reads=0, 
                                                obs_ref_reads=100, 
                                                flank_seg_cn=2, 
                                                orientation="decreasing"))
  expect_equal(0, cn_from_vaf_with_flank_seg_cn(obs_sup_reads=100, 
                                                obs_ref_reads=0, 
                                                flank_seg_cn=2, 
                                                orientation="decreasing"))
  test_symmetry = function(seg_a, sv_count, ref_count) {
    seg_b = cn_from_vaf_with_flank_seg_cn(sv_count, ref_count, seg_a, "increasing")
    expect_equal(seg_a, cn_from_vaf_with_flank_seg_cn(sv_count, ref_count, seg_b, "decreasing"))
    seg_b = cn_from_vaf_with_flank_seg_cn(sv_count, ref_count, seg_a, "decreasing")
    expect_equal(seg_a, cn_from_vaf_with_flank_seg_cn(sv_count, ref_count, seg_b, "increasing"))
  }
})

test_that("sv_vaf", {
  expect_equal(sv_vaf(1, 1), 0)
  expect_equal(sv_vaf(1, 2), 0.5)
  expect_equal(sv_vaf(2, 1), 0.5)
  expect_equal(sv_vaf(0, 1), 1)
  expect_equal(sv_vaf(1, 0), 1)
  expect_equal(sv_vaf(4, 1), 0.75)
  expect_equal(sv_vaf(1, 4), 0.75)
  expect_equal(sv_vaf(0, 0), NaN)
  expect_equal(sv_vaf(c(1, 2), c(2, 4)), c(0.5, 0.5))
})


# Test that rd_model input argument behave as expected. 
test_that("cn_from_rd", {
  expect_equal(cn_from_rd(obs_reads = 1000, 
                          seg_length = 250, 
                          read_length = 150, 
                          coverage = 30,
                          rd_model = "start_in"), 20)
  expect_equal(cn_from_rd(obs_reads = 1000, 
                          seg_length = 250, 
                          read_length = 150, 
                          coverage = 30,
                          rd_model = "overlap"), 12.5)
})


# Test that rd_model incorporated into cn_from_vaf input argument behave as expected. 
test_that("cn_from_vaf", {
  expect_equal(cn_from_vaf(obs_sup_reads = 100, 
                           obs_ref_reads = 200, 
                           read_length = 150, 
                           flank_seg_reads = 1015, 
                           flank_seg_length = 10000,
                           coverage = 30,
                           orientation = "increasing",
                           rd_model = "overlap"), 0.75)
  expect_equal(cn_from_vaf(obs_sup_reads = 100, 
                           obs_ref_reads = 200, 
                           read_length = 150, 
                           flank_seg_reads = 1000, 
                           flank_seg_length = 10000,
                           coverage = 30,
                           orientation = "increasing",
                           rd_model = "start_in"), 0.75)
  expect_equal(cn_from_vaf(obs_sup_reads = 100, 
                           obs_ref_reads = 200, 
                           read_length = 150, 
                           flank_seg_reads = 1000, 
                           flank_seg_length = 10000,
                           coverage = 30,
                           orientation = "decreasing",
                           rd_model = "overlap"), 0.329, 
               tolerance = 0.001)
  expect_equal(cn_from_vaf(obs_sup_reads = 100, 
                           obs_ref_reads = 200, 
                           read_length = 150, 
                           flank_seg_reads = 1000, 
                           flank_seg_length = 10000,
                           coverage = 30,
                           orientation = "decreasing",
                           rd_model = "start_in"), 0.333, 
               tolerance = 0.001)
})


# Adding the rd_model argument to expected_reads is causing issue for calculating sv_vaf cn
# This is maybe because my implementation of vaf_model in the code is wrong (lines 132 to 139 in 
# Simulation_functions.R)
test_that("obs_cn", {
  df = obs_cn(n = 1, 
       read_length = 150, 
       coverage = 1000000,
       flank_seg_length = 10000,
       cur_seg_length = 250,
       actual_flank_seg_cn = 2,
       actual_cur_seg_cn = 4,
       fragment_size = 300,
       rd_model = "start_in",
       vaf_model = "read_pair")
  expect_equal(as.data.frame(df), data.frame(
                              model = c("vaf read_pair", "rd start_in"),
                              cn_estimate = 4,
                              read_length = 150,
                              fragment_size = 300,
                              coverage = 1000000,
                              flank_seg_length = 10000,
                              cur_seg_length = 250,
                              actual_flank_seg_cn = 2,
                              actual_cur_seg_cn = 4), tolerance = 0.01)
  })

obs_cn(n = 1, 
       read_length = 150, 
       coverage = 30,
       flank_seg_length = 10000,
       cur_seg_length = 250,
       actual_flank_seg_cn = 4,
       actual_cur_seg_cn = 4,
       rd_model = "overlap",
       vaf_model = "read_pair")


obs_cn(n = 1, 
       read_length = 150, 
       coverage = 30,
       flank_seg_length = 10000,
       cur_seg_length = 250,
       actual_flank_seg_cn = 2,
       actual_cur_seg_cn = 4,
       rd_model = "overlap",
       vaf_model = "split_read")

test_that("sim_permutation cardinality", {
  df = sim_permutations(
    read_length = c(10, 20),
    coverage = 30,
    cur_seg_length = 5000,
    flank_seg_length = 100000,
    actual_flank_seg_cn = 3,
    actual_cur_seg_cn = c(2, 4),
    fragment_size = 300,
    n=50,
    rd_model = "overlap",
    vaf_model = "read_pair")
  expect_equal(2 * 4, nrow(df))
})

test_that("obs_cn_width", {
  df = obs_cn_width(n = 10, 
              read_length = 100, 
              coverage = 1000000,
              flank_seg_length = 10000,
              cur_seg_length = 10000,
              actual_flank_seg_cn = 2,
              actual_cur_seg_cn = 1,
              fragment_size = 300,
              rd_model = "overlap",
              vaf_model = "split_read")
  expect_equal(as.data.frame(df), data.frame(model = c("rd overlap", "vaf split_read"),
                              read_length = 100, 
                              coverage = 1000000,
                              cur_seg_length = 10000,
                              flank_seg_length = 10000,
                              actual_flank_seg_cn = 2,
                              actual_cur_seg_cn = 1,
                              ci_width = c(0,0)), tolerance = 0.01)
  })

test_that("expected_reads", {
  expect_equal(30, expected_reads(cn = 1, 
                                  coverage = 30, 
                                  read_length = 100, 
                                  seg_length = 100, 
                                  rd_model = "start_in"))
  expect_equal(60, expected_reads(cn = 1, 
                                  coverage = 30, 
                                  read_length = 100, 
                                  seg_length = 100, 
                                  rd_model = "overlap"))
               })
# expect_error()

test_that("estimate_ci_width", {
  expect_equal(95, estimate_ci_width(1:100))
  expect_equal(1, estimate_ci_width(c(0, 1)))
})
test_that("expected_segment_reads", {
  expect_equal(2, expected_segment_reads(2, 1, 100, 100, "start_in"))
  expect_equal(6, expected_segment_reads(2, 3, 100, 100, "start_in"))
  expect_equal(20, expected_segment_reads(2, 1, 100, 1000, "start_in"))
  expect_equal(22, expected_segment_reads(2, 1, 100, 1000, "overlap"))
})
test_that("expected_sv_reads", {
  expect_equal(2, expected_sv_reads(2, 1, 100, 300, "split_read"))
  expect_equal(6, expected_sv_reads(2, 3, 100, 300, "split_read"))
  expect_equal(18, expected_sv_reads(2, 3, 100, 300, "read_pair"))
})
test_that("expected_reads", {
  expect_equal(1, expected_reads(1, 1, 10, 10))
  expect_equal(10, expected_reads(1, 1, 10, 100))
  expect_equal(6, expected_reads(2, 3, 10, 10))
})

# test_that("segment_CN_function_check", {
#   expect_equal()
# })