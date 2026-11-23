library(tinytest)

# Statistical properties of the CLUMPY allocators (src/alloc_clumpy.cpp) that Mazy (2022),
# https://theses.hal.science/tel-04382012v1, Ch. 3, claims for uSAM and uPAM. Each check is one
# seeded Monte Carlo run compared with its expectation, within `z_max` standard errors: a
# correct allocator fails a check about once in 16 000 seeds, and the seeds are fixed. The
# deterministic properties are in test_alloc_clumpy.R.

if (!at_home()) {
  exit_file("Monte Carlo tests skipped (not at_home)")
}

z_max <- 4
# TRUE where an observed frequency lies within z_max binomial standard errors of p
within_binomial <- function(frequency, p, n) {
  abs(frequency - p) <= z_max * sqrt(p * (1 - p) / n)
}

# 150 000 cells of class 1 in three strata of potential, in random order over the raster, and a
# transition 1 -> 2 demanded at their mean
strata <- c(0.02, 0.1, 0.5)
n_stratum <- 50000L
set.seed(1L)
potential <- sample(rep(strata, each = n_stratum))
allocate_strata <- function(seed, method, batch_size = 0L, area_mean = 1, rarefy = TRUE) {
  set.seed(seed)
  evoland.plus:::allocate_clumpy_cpp(
    landscape = rep(1L, length(potential)),
    nrow = 300L,
    ncol = 500L,
    trans_from = 1L,
    trans_to = 2L,
    prob_cell = list(seq_along(potential)),
    prob_value = list(potential),
    area_mean = area_mean,
    area_var = 0,
    elongation = 0,
    target_rate = mean(strata),
    method = method,
    batch_size = batch_size,
    rarefy = rarefy,
    shuffle = TRUE,
    avoid_aggregation = FALSE,
    area_dist = 1L
  )
}
# share of changed cells per stratum
changed_by_stratum <- function(posterior) {
  as.numeric(tapply(posterior == 2L, potential, mean))
}

# uSAM is unbiased ####
# Mazy eqs. 3.6 and 3.7, App. 3.B.2: every cell changes with its own potential, so the changed
# cells are distributed over the potentials as the potentials prescribe, and their number is
# the sum of the potentials on average
usam <- allocate_strata(seed = 2L, method = 0L)
expect_true(all(within_binomial(changed_by_stratum(usam), strata, n_stratum)))
expect_true(
  abs(sum(usam == 2L) - sum(potential)) <= z_max * sqrt(sum(potential * (1 - potential)))
)

# with two transitions from one class, each cell goes to either with its potential
set.seed(3L)
two_way <- evoland.plus:::allocate_clumpy_cpp(
  landscape = rep(1L, 40000L),
  nrow = 200L,
  ncol = 200L,
  trans_from = c(1L, 1L),
  trans_to = c(2L, 3L),
  prob_cell = list(1:40000, 1:40000),
  prob_value = list(rep(0.1, 40000L), rep(0.3, 40000L)),
  area_mean = c(1, 1),
  area_var = c(0, 0),
  elongation = c(0, 0),
  target_rate = c(0.1, 0.3),
  method = 0L,
  batch_size = 0L,
  rarefy = TRUE,
  shuffle = TRUE,
  avoid_aggregation = FALSE,
  area_dist = 1L
)
expect_true(all(within_binomial(
  c(mean(two_way == 2L), mean(two_way == 3L)),
  c(0.1, 0.3),
  40000L
)))

# uPAM pivots ####
# Mazy eq. 3.7, App. 3.D: with mono-pixel patches, the pivots are the allocated cells. uPAM
# draws them in repeated MuST passes from a shrinking pool; with the probabilities updated as an
# urn per bin of potential (Mazy & Longaretti 2023, App. D.3, E.2.3), they follow the
# potentials as in uSAM, in one pass over all candidates (batch_size < 0) as in the default
# batches of about 1 % of the pool. Without the update, the cells of high potential, drawn
# early and then gone from the pool, were under-represented (0.4755 instead of 0.5).
upam_one_pass <- allocate_strata(seed = 2L, method = 1L, batch_size = -1L)
expect_true(all(within_binomial(changed_by_stratum(upam_one_pass), strata, n_stratum)))
upam_passes <- allocate_strata(seed = 2L, method = 1L)
expect_true(all(within_binomial(changed_by_stratum(upam_passes), strata, n_stratum)))
# the quantity is exact all the same: 0.62 / 3 of 150 000 cells
expect_equal(sum(upam_passes == 2L), 31000L)

# rarefaction ####
# Mazy eqs. 3.11 and 3.12, Fig. 3.2: a pivot grows into a patch of E(sigma) cells on average,
# so the pivot probability is the potential over E(sigma). Seen on the single MuST pass of
# uSAM, where pivots are the changed cells: with a mean area of 4, a quarter of the potentials
rarefied <- allocate_strata(seed = 4L, method = 0L, area_mean = 4)
expect_true(all(within_binomial(changed_by_stratum(rarefied), strata / 4, n_stratum)))
not_rarefied <- allocate_strata(seed = 4L, method = 0L, area_mean = 4, rarefy = FALSE)
expect_true(all(within_binomial(changed_by_stratum(not_rarefied), strata, n_stratum)))

# patch areas ####
# Mazy eq. 3.10: the areas of the allocated patches are distributed as prescribed. First the
# samplers: mean and variance of 20 000 areas, which are rounded to whole cells (adding 1/12
# to the variance)
set.seed(5L)
areas <- list(
  lognormal = vapply(1:20000, \(i) evoland.plus:::sample_lognorm_area_cpp(10, 9), integer(1)),
  normal = vapply(1:20000, \(i) evoland.plus:::sample_normal_area_cpp(10, 4), integer(1))
)
expect_equal(mean(areas[["lognormal"]]), 10, tolerance = 0.01)
expect_equal(stats::var(areas[["lognormal"]]), 9 + 1 / 12, tolerance = 0.05)
expect_equal(mean(areas[["normal"]]), 10, tolerance = 0.01)
expect_equal(stats::var(areas[["normal"]]), 4 + 1 / 12, tolerance = 0.05)

# then the patches of an allocation, where they do not merge: 300 x 300 cells, 2 % of them
# demanded in patches of 8 cells on average with a variance of 4, elongation 0.5
allocate_patches <- function(seed, elongation) {
  set.seed(seed)
  posterior <- evoland.plus:::allocate_clumpy_cpp(
    landscape = rep(1L, 90000L),
    nrow = 300L,
    ncol = 300L,
    trans_from = 1L,
    trans_to = 2L,
    prob_cell = list(1:90000),
    prob_value = list(rep(0.02, 90000L)),
    area_mean = 8,
    area_var = 4,
    elongation = elongation,
    target_rate = 0.02,
    method = 1L,
    batch_size = 0L,
    rarefy = TRUE,
    shuffle = TRUE,
    avoid_aggregation = TRUE,
    area_dist = 1L
  )
  ifelse(posterior == 2L, 1L, NA_integer_) |>
    matrix(nrow = 300L)
}
patches <- allocate_patches(seed = 6L, elongation = 0.5)
by_rook <- evoland.plus:::calculate_class_stats_cpp(patches, cellsize = 1, directions = 4L)
# about 1800 cells / 8 patches; the standard error of their mean area is sqrt(4 / 225)
expect_true(abs(by_rook[["patch_area_mean"]] - 8) <= z_max * sqrt(4 / by_rook[["patch_count"]]))
expect_equal(by_rook[["patch_area_variance"]], 4 + 1 / 12, tolerance = 0.25)

# estimated patch parameters ####
# what calculate_class_stats_cpp() measures on an allocation is what the allocation was given.
# The elongation the grower steers towards (Mazy App. 3.I.3.3) is reached on average:
expect_equal(by_rook[["patch_elongation_mean"]], 0.5, tolerance = 0.1)
elongated <-
  allocate_patches(seed = 6L, elongation = 0.8) |>
  evoland.plus:::calculate_class_stats_cpp(cellsize = 1, directions = 4L)
expect_equal(elongated[["patch_elongation_mean"]], 0.8, tolerance = 0.1)
# This needs rook adjacency, which create_alloc_params_t() uses: by queen adjacency, patches
# that touch by a corner count as one, which the grower does not prevent, so fewer and larger
# patches are measured than were allocated.
by_queen <- evoland.plus:::calculate_class_stats_cpp(patches, cellsize = 1)
expect_true(by_queen[["patch_count"]] < by_rook[["patch_count"]])
expect_true(by_queen[["patch_area_mean"]] > by_rook[["patch_area_mean"]])
