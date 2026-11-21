library(tinytest)

overrides <- data.table::data.table(id_run = 1L, id_trans = 1:2, is_viable = c(FALSE, TRUE))
expect_stdout(print(as_trans_viability_t(overrides)), "allowed: 1, disallowed: 1")

# an override says either way; NA would fall through to the default unnoticed
expect_error(
  as_trans_viability_t(list(id_run = 1L, id_trans = 1L, is_viable = NA)),
  "is_viable must be TRUE or FALSE"
)

# create_trans_viability_t() derives the defaults from thresholds and exclusions
trans_meta <- as_trans_meta_t(data.table::data.table(
  id_trans = 1:4,
  id_lulc_anterior = c(1L, 1L, 2L, 3L),
  id_lulc_posterior = c(2L, 3L, 1L, 1L),
  cardinality = c(10L, 2L, 6L, 2L),
  frequency_rel = c(0.5, 0.1, 0.3, 0.1),
  frequency_abs = c(0.1, 0.02, 0.06, 0.02)
))
viable_by <- function(...) create_trans_viability_t(trans_meta, ...)[["is_viable"]]
expect_equal(viable_by(), rep(TRUE, 4L))
expect_equal(viable_by(min_cardinality_abs = 5L), c(TRUE, FALSE, TRUE, FALSE))
expect_equal(viable_by(min_frequency_rel = 0.4), c(TRUE, FALSE, FALSE, FALSE))
# with both thresholds, meeting either keeps a transition viable
expect_equal(
  viable_by(min_cardinality_abs = 5L, min_frequency_rel = 0.4),
  c(TRUE, FALSE, TRUE, FALSE)
)
expect_equal(viable_by(exclude_anterior = 3L, exclude_posterior = 2L), c(FALSE, TRUE, TRUE, FALSE))
expect_equal(unique(create_trans_viability_t(trans_meta, id_run = 2L)[["id_run"]]), 2L)
