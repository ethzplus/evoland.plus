library(tinytest)

# The CLUMPY allocation backend (src/alloc_clumpy.cpp) against Mazy (2022),
# https://theses.hal.science/tel-04382012v1, Ch. 3. Everything here is deterministic or
# reproduces the random draws in R: the C++ code draws from R's generator, so after the same
# set.seed() the test knows every uniform the allocator sees. Properties that only hold on
# average are in test_alloc_clumpy_stat.R.

# must_cpp: Multinomial Sampling Test ####
# Mazy App. 3.B.1, eq. 3.B.1: one uniform u per cell picks the first state w whose cumulative
# probability exceeds u; the last column is the stay state, taken when no other is.
row_cumsum <- function(x) {
  for (j in seq_len(ncol(x))[-1L]) {
    x[, j] <- x[, j] + x[, j - 1L]
  }
  x
}
must_oracle <- function(prob, states, u) {
  states[pmin(rowSums(u >= row_cumsum(prob)) + 1L, ncol(prob))]
}
set.seed(1L)
must_prob <- matrix(stats::runif(600L), ncol = 3L) / 3 # rows sum below 1
must_states <- c(10L, 20L, 30L)
set.seed(2L)
must_drawn <- evoland.plus:::must_cpp(must_prob, must_states)
set.seed(2L)
must_u <- stats::runif(200L)
expect_equal(must_drawn, must_oracle(must_prob, must_states, must_u))
# all three states occur, so the comparison exercises every branch
expect_equal(sort(unique(must_drawn)), must_states)
# negative and NaN probabilities count as 0 (as in the reference clumpy), whatever is drawn
expect_equal(
  evoland.plus:::must_cpp(matrix(c(-1, NaN, 1), nrow = 1L), 1:3),
  3L
)

# area samplers ####
# both are parameterised by the mean E and variance V of the patch area in cells, rounded to
# whole cells and at least 1; the draws below are reproduced with R's own samplers
set.seed(1L)
lognorm_drawn <- vapply(1:50, \(i) evoland.plus:::sample_lognorm_area_cpp(4, 2), integer(1))
set.seed(1L)
expect_equal(
  lognorm_drawn,
  # a log-normal with mean E and variance V has sdlog^2 = log(1 + V / E^2) and
  # meanlog = log(E) - sdlog^2 / 2
  stats::rlnorm(50L, meanlog = log(4) - log(1 + 2 / 16) / 2, sdlog = sqrt(log(1 + 2 / 16))) |>
    round() |>
    pmax(1)
)
set.seed(1L)
normal_drawn <- vapply(1:50, \(i) evoland.plus:::sample_normal_area_cpp(4, 2), integer(1))
set.seed(1L)
expect_equal(normal_drawn, pmax(round(pmax(stats::rnorm(50L, 4, sqrt(2)), 1)), 1))
# without variance (0, or NA for a transition observed as a single patch) the area is the mean
expect_equal(evoland.plus:::sample_lognorm_area_cpp(5, 0), 5L)
expect_equal(evoland.plus:::sample_lognorm_area_cpp(1, 0), 1L)
expect_equal(evoland.plus:::sample_lognorm_area_cpp(5, NA_real_), 5L)
expect_equal(evoland.plus:::sample_normal_area_cpp(5, 0), 5L)
# no usable mean: a single cell
expect_equal(evoland.plus:::sample_lognorm_area_cpp(0, 1), 1L)
expect_equal(evoland.plus:::sample_lognorm_area_cpp(NA_real_, 1), 1L)
expect_equal(evoland.plus:::sample_normal_area_cpp(0, 1), 1L)

# patch elongation ####
# Mazy eq. 3.I.12: e = 1 - sqrt(lambda_min / lambda_max) of the second central moments of the
# cell coordinates. Measured through calculate_class_stats_cpp(), which shares the definition
# (src/clumpy_geometry.h) with the patch grower.
patch_stats <- function(cells, nrow = 4L, ncol = 4L, directions = 8L) {
  patch <- matrix(NA_integer_, nrow, ncol)
  patch[cells] <- 1L
  evoland.plus:::calculate_class_stats_cpp(patch, cellsize = 1, directions = directions)
}
elongation <- function(...) patch_stats(...)[["patch_elongation_mean"]]
expect_equal(elongation(cbind(1, 1)), 0) # a single cell is isometric by convention
expect_equal(elongation(cbind(1, 1:2)), 1) # domino: all variance along one axis
expect_equal(elongation(cbind(1, 1:4)), 1)
expect_equal(elongation(cbind(c(1, 1, 2, 2), c(1, 2, 1, 2))), 0) # 2 x 2 square
# L-tromino (1,1), (1,2), (2,1): the moment matrix [[2, -1], [-1, 2]] / 9 has eigenvalues 1/9
# and 3/9
expect_equal(elongation(cbind(c(1, 1, 2), c(1, 2, 1))), 1 - 1 / sqrt(3))
# the eigenvalues do not depend on position or on which axis is which
expect_equal(elongation(cbind(c(3, 3, 4), c(2, 3, 2))), 1 - 1 / sqrt(3))
expect_equal(elongation(cbind(c(1, 2, 1), c(1, 1, 2))), 1 - 1 / sqrt(3))
# two cells touching by a corner are one patch by queen adjacency, two by rook adjacency
expect_equal(patch_stats(cbind(1:2, 1:2))[["patch_count"]], 1L)
expect_equal(patch_stats(cbind(1:2, 1:2), directions = 4L)[["patch_count"]], 2L)
expect_error(patch_stats(cbind(1, 1), directions = 6L), "directions must be 4 or 8")

# raster_neighbors_cpp ####

expect_equal(
  evoland.plus:::raster_neighbors_cpp(3L, 4L), # 3-row, 4-col raster (12 cells)
  list(
    above = c(0L, 0L, 0L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L),
    below = c(5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 0L, 0L, 0L, 0L),
    left = c(0L, 1L, 2L, 3L, 0L, 5L, 6L, 7L, 0L, 9L, 10L, 11L),
    right = c(2L, 3L, 4L, 0L, 6L, 7L, 8L, 0L, 10L, 11L, 12L, 0L)
  )
)

# grow_patch_cpp ####

# 4x4 raster, all class 1, no obstacles
n <- 16L
landscape <- as.integer(rep(1L, n))
ant_land <- as.integer(rep(1L, n))
probs <- rep(0.8, n)
nbrs <- evoland.plus:::raster_neighbors_cpp(4L, 4L)

patch <- evoland.plus:::grow_patch_cpp(
  landscape = landscape, # changes landscape by reference
  ant_landscape = ant_land,
  probs = probs,
  nbr_above = nbrs$above,
  nbr_below = nbrs$below,
  nbr_left = nbrs$left,
  nbr_right = nbrs$right,
  pivot = 1L,
  target_area = 4L,
  from_class = 1L,
  to_class = 2L,
  elongation = 0.5,
  ncol = 4L
)
expect_equal(patch, c(1L, 2L, 5L, 3L))
expect_equal(
  landscape,
  c(2L, 2L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L)
)

# Pivot with wrong class -> empty result
landscape_wrong <- as.integer(rep(2L, n))
patch_empty <- evoland.plus:::grow_patch_cpp(
  landscape = landscape_wrong,
  ant_landscape = ant_land,
  probs = probs,
  nbr_above = nbrs$above,
  nbr_below = nbrs$below,
  nbr_left = nbrs$left,
  nbr_right = nbrs$right,
  pivot = 1L,
  target_area = 4L,
  from_class = 1L,
  to_class = 2L,
  elongation = 0.5,
  ncol = 4L
)
expect_length(patch_empty, 0L)

# avoid_aggregation: a patch that would touch an existing patch fails entirely.
# 1x4 row [1,1,2,1] (cell 3 is an existing class-2 patch, originally class 1).
nbr14 <- evoland.plus:::raster_neighbors_cpp(1L, 4L)
land_adj <- as.integer(c(1L, 1L, 2L, 1L))
ant_adj <- as.integer(c(1L, 1L, 1L, 1L))
probs_adj <- rep(1.0, 4L)

# With avoid_aggregation = TRUE, growing from cell 1 toward area 3 hits cell 3
# (a foreign patch) and fails -> nothing allocated.
patch_agg <- evoland.plus:::grow_patch_cpp(
  landscape = land_adj,
  ant_landscape = ant_adj,
  probs = probs_adj,
  nbr_above = nbr14$above,
  nbr_below = nbr14$below,
  nbr_left = nbr14$left,
  nbr_right = nbr14$right,
  pivot = 1L,
  target_area = 3L,
  from_class = 1L,
  to_class = 2L,
  elongation = 0.0,
  ncol = 4L,
  avoid_aggregation = TRUE
)
expect_length(patch_agg, 0L)

# With avoid_aggregation = FALSE, it grows as far as it can (cells 1 and 2).
land_adj2 <- as.integer(c(1L, 1L, 2L, 1L))
patch_noagg <- evoland.plus:::grow_patch_cpp(
  landscape = land_adj2,
  ant_landscape = ant_adj,
  probs = probs_adj,
  nbr_above = nbr14$above,
  nbr_below = nbr14$below,
  nbr_left = nbr14$left,
  nbr_right = nbr14$right,
  pivot = 1L,
  target_area = 3L,
  from_class = 1L,
  to_class = 2L,
  elongation = 0.0,
  ncol = 4L,
  avoid_aggregation = FALSE
)
expect_equal(sort(patch_noagg), c(1L, 2L))

# allocate_clumpy_cpp ####
# one transition 1 -> 2 on a square of class 1 unless said otherwise; potentials are recycled
# over the cells
allocate_square <- function(
  n_side,
  potential,
  rate,
  area_mean = 1,
  area_var = 0,
  method = 1L,
  trans_to = 2L,
  batch_size = 0L,
  rarefy = TRUE,
  avoid_aggregation = TRUE,
  elongation = 0,
  area_dist = 1L, # normal: with a variance of 0, every patch has exactly area_mean cells
  extra_trans_to = integer(0)
) {
  n_cells <- n_side^2
  n_extra <- length(extra_trans_to)
  evoland.plus:::allocate_clumpy_cpp(
    landscape = rep(1L, n_cells),
    nrow = n_side,
    ncol = n_side,
    trans_from = rep(1L, 1L + n_extra),
    trans_to = c(trans_to, extra_trans_to),
    prob_cell = c(list(seq_len(n_cells)), rep(list(integer(0)), n_extra)),
    prob_value = c(list(rep_len(potential, n_cells)), rep(list(numeric(0)), n_extra)),
    area_mean = rep(area_mean, 1L + n_extra),
    area_var = rep(area_var, 1L + n_extra),
    elongation = rep(elongation, 1L + n_extra),
    target_rate = c(rate, rep(0, n_extra)),
    method = method,
    batch_size = batch_size,
    rarefy = rarefy,
    shuffle = TRUE,
    avoid_aggregation = avoid_aggregation,
    area_dist = area_dist
  )
}
n_changed <- function(seed, ...) {
  set.seed(seed)
  sum(allocate_square(...) == 2L)
}
# sizes of the rook-connected patches of class 2
rook_patches <- function(posterior, n_side) {
  ifelse(posterior == 2L, 1L, NA_integer_) |>
    matrix(nrow = n_side) |>
    evoland.plus:::calculate_class_stats_cpp(cellsize = 1, directions = 4L)
}

## uSAM ####
# Mazy sec. 3.4.1: every cell of an anterior class is tested once with MuST, anterior classes
# in ascending order, cells in raster order. The posterior map follows from the uniforms.
# `potential` is a cells x transitions matrix.
usam_oracle <- function(landscape, trans_from, trans_to, potential) {
  posterior <- landscape
  for (id_lulc in sort(unique(trans_from))) {
    pool <- which(landscape == id_lulc)
    trans <- which(trans_from == id_lulc)
    u <- stats::runif(length(pool))
    selected <- rowSums(u >= row_cumsum(potential[pool, trans, drop = FALSE])) + 1L
    changes <- selected <= length(trans)
    posterior[pool[changes]] <- trans_to[trans[selected[changes]]]
  }
  posterior
}
# three classes with a block without data; class 1 can become 2 or 3, class 3 can become 1,
# class 2 does not change
set.seed(7L)
usam_landscape <- sample(1:3, 400L, replace = TRUE)
usam_landscape[101:140] <- NA_integer_
usam_from <- c(1L, 1L, 3L)
usam_to <- c(2L, 3L, 1L)
usam_potential <- matrix(stats::runif(1200L, 0, 0.4), ncol = 3L)
for (trans in 1:3) {
  usam_potential[!usam_landscape %in% usam_from[trans], trans] <- 0
}
usam_cells <- lapply(1:3, \(trans) which(usam_potential[, trans] > 0))
set.seed(99L)
usam_posterior <- evoland.plus:::allocate_clumpy_cpp(
  landscape = usam_landscape,
  nrow = 20L,
  ncol = 20L,
  trans_from = usam_from,
  trans_to = usam_to,
  prob_cell = usam_cells,
  prob_value = lapply(1:3, \(trans) usam_potential[usam_cells[[trans]], trans]),
  area_mean = c(1, 1, 1),
  area_var = c(0, 0, 0),
  elongation = c(0, 0, 0),
  target_rate = c(0.2, 0.2, 0.2),
  method = 0L,
  batch_size = 0L,
  rarefy = TRUE,
  shuffle = TRUE,
  avoid_aggregation = FALSE,
  area_dist = 0L
)
set.seed(99L)
expect_equal(usam_posterior, usam_oracle(usam_landscape, usam_from, usam_to, usam_potential))
# a cell changes at most once and only along a declared transition: class 2 and the cells
# without data are untouched (implied by the oracle, stated for the reader)
expect_equal(usam_posterior[usam_landscape %in% 2L], rep(2L, sum(usam_landscape %in% 2L)))
expect_true(all(is.na(usam_posterior[101:140])))

# a potential of 1 changes every cell, whatever is drawn; no potential changes none
expect_equal(n_changed(1L, 5L, potential = 1, rate = 1, method = 0L), 25L)
expect_equal(n_changed(1L, 5L, potential = 0, rate = 0.3), 0L)
# only cells that carry a potential can change
set.seed(9L)
expect_equal(
  which(allocate_square(5L, rep(c(1, 0, 0, 0, 0, 0), length.out = 25L), 1, method = 0L) == 2L),
  c(1L, 7L, 13L, 19L, 25L)
)

## uPAM: updating the probabilities ####
# Mazy & Longaretti (2023), App. D.1, E.2.2: an urn per bin, here a bin of potential, with N
# cells left and N_v of them still due to change; a cell's probability is its own potential
# times N_v over the potentials of the N cells. 10 cells at 0.2 and 10 at 0.6, quota 8: two
# bins, as equal potentials are never split, expecting 2 and 6 cells
urn <- function(potential, quota, n_bins, pivot = integer(0), size = integer(0), removed = pivot) {
  evoland.plus:::urn_bins_cpp(potential, quota, n_bins, pivot, size, removed)
}
two_strata <- rep(c(0.2, 0.6), each = 10L)
expect_equal(
  urn(two_strata, 8, n_bins = 10L),
  list(bin = rep(1:2, each = 10L), probability = two_strata, expected = c(2, 6))
)
# three cells of the second bin allocated, two of the first failed: 3 balls of state v among 7
# in the second, 2 among 8 in the first
after <- urn(
  two_strata,
  8,
  n_bins = 10L,
  pivot = 11:13,
  size = rep(1L, 3L),
  removed = c(1:2, 11:13)
)
expect_equal(after[["expected"]], c(2, 3))
expect_equal(after[["probability"]], c(0, 0, rep(2 / 8, 8L), 0, 0, 0, rep(3 / 7, 7L)))
# a patch of 3 cells from a pivot in the second bin, spilling into the first, counts as 3 cells
# of the second bin; its other cells leave the first bin as failed cells do
spilled <- urn(two_strata, 8, n_bins = 10L, pivot = 11L, size = 3L, removed = c(1:2, 11L))
expect_equal(spilled[["expected"]], c(2, 3))
expect_equal(spilled[["probability"]][c(3L, 12L)], c(0.2 * 2 / 1.6, 0.6 * 3 / 5.4))
# the expected counts are scaled to the quota, here half the sum of the potentials
expect_equal(urn(two_strata, 4, n_bins = 10L)[["probability"]], two_strata / 2)
# within a bin, the update is shared in proportion to the potentials: one bin of 0.2 and 0.6,
# expecting 1.6 cells; once one cell is allocated, 0.6 cells are due from a mass of 1
shared <- urn(c(0.2, 0.6, 0.2, 0.6), 1.6, n_bins = 1L, pivot = 2L, size = 1L)
expect_equal(shared[["probability"]], c(0.12, 0, 0.12, 0.36))
# bins of about equal potential mass: 0.01, ..., 0.10 in two bins, the first closed once it
# holds half the mass (0.01 + ... + 0.07 = 0.28 >= 0.275); cells without potential have none
expect_equal(urn(c(1:10 / 100, 0), 0.55, n_bins = 2L)[["bin"]], c(rep(1L, 7L), rep(2L, 3L), 0L))
# once a bin's expected count is used up, its cells are not drawn any more
expect_equal(
  urn(two_strata, 8, n_bins = 10L, pivot = 1L, size = 2L)[["probability"]],
  c(0, rep(0, 9L), rep(0.6, 10L))
)

## uPAM: quantity of change ####
# Mazy sec. 3.4.2: patches are allocated until the quantity of change of the transition is
# used up, so the allocated cells reach the quota and exceed it by less than one patch. The
# quota is round(rate * cells of the anterior class). Checked over seeds, with patches of
# exactly 6 cells on 1600 cells: quota 80, i.e. 14 patches.
seeds <- 1:20
expect_equal(
  vapply(seeds, n_changed, numeric(1), n_side = 40L, potential = 0.05, rate = 0.05, area_mean = 6),
  rep(84, 20L)
)
# the same when patches may merge
expect_equal(
  vapply(
    seeds,
    n_changed,
    numeric(1),
    n_side = 40L,
    potential = 0.05,
    rate = 0.05,
    area_mean = 6,
    avoid_aggregation = FALSE
  ),
  rep(84, 20L)
)
# a quota of a single patch: with 400 cells at potential 0.01 and patches of 4 cells, a MuST
# pass selects no pivot with probability exp(-1); the pass is repeated until one is
expect_equal(
  vapply(seeds, n_changed, numeric(1), n_side = 20L, potential = 0.01, rate = 0.01, area_mean = 4),
  rep(4, 20L)
)
# the quota is a whole number of cells: 0.026 * 400 = 10.4 gives 10, 0.001 * 400 = 0.4 gives 0
expect_equal(
  vapply(seeds, n_changed, numeric(1), n_side = 20L, potential = 0.026, rate = 0.026),
  rep(10, 20L)
)
expect_equal(
  vapply(seeds, n_changed, numeric(1), n_side = 20L, potential = 0.001, rate = 0.001),
  rep(0, 20L)
)
# log-normal areas: the bound holds with the largest patch of the realisation
expect_true(all(vapply(
  seeds,
  function(seed) {
    set.seed(seed)
    posterior <-
      allocate_square(40L, 0.05, 0.05, area_mean = 6, area_var = 4, area_dist = 0L)
    largest <- max(table(terra::values(terra::patches(
      terra::rast(matrix(ifelse(posterior == 2L, 1, NA), 40L, byrow = TRUE)),
      directions = 4L
    ))))
    sum(posterior == 2L) >= 80L && sum(posterior == 2L) < 80L + largest
  },
  logical(1)
)))

## uPAM: patches ####
# Mazy sec. 3.2.4 and App. 3.D.1.2: patches of one period must not merge, or the patch sizes
# shift towards larger patches. With avoid_aggregation and patches of exactly 6 cells at a
# density where patches crowd (a tenth of the map), every rook-connected patch has 6 cells.
set.seed(3L)
crowded <- allocate_square(40L, 0.1, 0.1, area_mean = 6)
expect_equal(
  rook_patches(crowded, 40L)[c("patch_count", "patch_area_mean", "patch_area_variance")],
  data.frame(patch_count = 27L, patch_area_mean = 6, patch_area_variance = 0)
)
expect_equal(sum(crowded == 2L), 27L * 6L)
# without it patches merge, and the same count of cells sits in fewer, larger patches
set.seed(3L)
merged <- allocate_square(40L, 0.1, 0.1, area_mean = 6, avoid_aggregation = FALSE)
expect_true(rook_patches(merged, 40L)[["patch_area_mean"]] > 6)

# a 1 x 5 row with potential 1 (every cell is a pivot, taken left to right), patches of 2
# cells and all cells demanded: without avoidance the patches tile the row; with it, the patch
# from cell 3 would touch the first one, fails and leaves cell 3 unchanged
allocate_row <- function(avoid_aggregation) {
  evoland.plus:::allocate_clumpy_cpp(
    landscape = rep(1L, 5L),
    nrow = 1L,
    ncol = 5L,
    trans_from = 1L,
    trans_to = 2L,
    prob_cell = list(1:5),
    prob_value = list(rep(1, 5L)),
    area_mean = 2,
    area_var = 0,
    elongation = 0,
    target_rate = 1,
    method = 1L,
    batch_size = 1L,
    rarefy = FALSE,
    shuffle = FALSE,
    avoid_aggregation = avoid_aggregation,
    area_dist = 1L
  )
}
expect_equal(allocate_row(avoid_aggregation = FALSE), c(2L, 2L, 2L, 2L, 2L))
expect_equal(allocate_row(avoid_aggregation = TRUE), c(2L, 2L, 1L, 2L, 2L))

## relabeling ####
# the realisation depends on the classes only through which cells belong to which transition:
# another id for the posterior class, or a further transition without potential or demand,
# leaves the allocated cells the same under one seed
set.seed(5L)
reference <- allocate_square(30L, 0.05, 0.05, area_mean = 4, area_var = 2, area_dist = 0L)
set.seed(5L)
relabeled <-
  allocate_square(30L, 0.05, 0.05, area_mean = 4, area_var = 2, area_dist = 0L, trans_to = 7L)
expect_equal(which(relabeled == 7L), which(reference == 2L))
set.seed(5L)
expect_equal(
  allocate_square(
    30L,
    0.05,
    0.05,
    area_mean = 4,
    area_var = 2,
    area_dist = 0L,
    extra_trans_to = 3L
  ),
  reference
)

# alloc_clumpy_spec on a hand-built allocation context ####
# a 10 x 10 grid of class 1 with one viable transition 1 -> 2 at rate 0.2, mono-pixel patches,
# so uSAM is selected
grid_coords <- data.table::CJ(lat = seq(5, 95, 10), lon = seq(5, 95, 10))[,
  .(id_coord = .I, lon, lat)
]
data.table::setattr(grid_coords, "epsg", 2056L)
data.table::setattr(grid_coords, "resolution", 10)
clumpy_ctx <- list(
  anterior = grid_coords[, .(id_coord, id_lulc = 1L)],
  viable_trans = data.table::data.table(
    id_trans = 1L,
    id_lulc_anterior = 1L,
    id_lulc_posterior = 2L
  ),
  quota = data.table::data.table(id_trans = 1L, rate = 0.2, n_anterior = 100L, quota = 20L),
  potentials = grid_coords[, .(id_trans = 1L, id_coord, value = 0.2)],
  alloc_params = data.table::data.table(
    id_trans = 1L,
    mean_patch_size = 1,
    patch_size_variance = 0,
    patch_elongation = NA_real_
  ),
  grid = evoland.plus:::coords_grid(grid_coords)
)
# the uSAM oracle on a context: id_coord counts from the bottom row of the grid, the allocator
# draws in raster order, from the top row
ctx_oracle <- function(ctx) {
  cells <- ctx$grid[["cells"]][ctx$anterior, on = "id_coord"]
  landscape <- cells[order(cell), id_lulc]
  potential <-
    ctx$potentials[cells, on = "id_coord"] |>
    data.table::dcast(cell ~ id_trans, value.var = "value", fill = 0)
  posterior <- usam_oracle(
    landscape,
    ctx$viable_trans[["id_lulc_anterior"]],
    ctx$viable_trans[["id_lulc_posterior"]],
    as.matrix(potential[, -"cell"])
  )
  cells[, .(id_coord, id_lulc = posterior[cell])]
}
set.seed(3L)
expect_message(
  clumpy_posterior <- alloc_clumpy_spec()$one_period(clumpy_ctx),
  "CLUMPY method: uSAM"
)
set.seed(3L)
# one row per row of the anterior map, in its order
expect_identical(clumpy_posterior, ctx_oracle(clumpy_ctx))
# uSAM meets the quota of 20 cells on average only (Mazy sec. 3.4.1); this realisation has 19
expect_equal(sum(clumpy_posterior[["id_lulc"]] == 2L), 19L)

# a viable transition at a rate of 0 has no positive potential anywhere; it still gets its
# (empty) element, so the other transition is allocated and it changes no cell
zero_rate_ctx <- clumpy_ctx
zero_rate_ctx$anterior <- grid_coords[, .(id_coord, id_lulc = rep(1:2, 50L))]
zero_rate_ctx$viable_trans <- data.table::data.table(
  id_trans = 1:2,
  id_lulc_anterior = 1:2,
  id_lulc_posterior = 2:1
)
zero_rate_ctx$quota <- data.table::data.table(
  id_trans = 1:2,
  rate = c(0, 0.2),
  n_anterior = 50L,
  quota = c(0L, 10L)
)
zero_rate_ctx$potentials <- rbind(
  grid_coords[, .(id_trans = 1L, id_coord, value = 0)],
  zero_rate_ctx$anterior[id_lulc == 2L, .(id_trans = 2L, id_coord, value = 0.2)]
)
zero_rate_ctx$alloc_params <- data.table::data.table(
  id_trans = 1:2,
  mean_patch_size = 1,
  patch_size_variance = 0,
  patch_elongation = NA_real_
)
set.seed(3L)
expect_message(
  zero_rate_posterior <- alloc_clumpy_spec()$one_period(zero_rate_ctx),
  "CLUMPY method: uSAM"
)
set.seed(3L)
expect_identical(zero_rate_posterior, ctx_oracle(zero_rate_ctx))
# class 1 keeps its 50 cells and gains those that left class 2
expect_true(sum(zero_rate_posterior[["id_lulc"]] == 1L) > 50L)

# patches larger than a cell select uPAM, which takes its quota from the rate: 20 cells
patch_ctx <- clumpy_ctx
patch_ctx$alloc_params <- data.table::data.table(
  id_trans = 1L,
  mean_patch_size = 2,
  patch_size_variance = 0,
  patch_elongation = NA_real_
)
set.seed(3L)
expect_message(
  patch_posterior <- alloc_clumpy_spec(area_dist = "normal")$one_period(patch_ctx),
  "CLUMPY method: uPAM"
)
expect_equal(sum(patch_posterior[["id_lulc"]] == 2L), 20L)

expect_error(alloc_clumpy_spec(default_elongation = 2), "default_elongation")
expect_error(alloc_clumpy_spec(n_bins = 0), "n_bins")
