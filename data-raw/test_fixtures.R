# Stores the test fixtures in R/sysdata.rda, for use in tests, examples, and vignettes.
#
# The fixtures are built by inst/tinytest/helper_fixtures.R, from one synthetic landscape.
# This script is not run during checks; the constructor tests rebuild each fixture with the
# same helpers and compare it against the stored one. Changing a fixture therefore requires
# rerunning this script, so that the stored data cannot go stale relative to the
# constructors it is built from.
library("evoland.plus")
library("data.table")

source("inst/tinytest/helper_fixtures.R")

landscape <- fixture_landscape()
test_coords_t <- fixture_coords_t(landscape)
test_periods_t <- fixture_periods_t()
test_lulc_meta_t <- fixture_lulc_meta_t()
test_lulc_data_t <- fixture_lulc_data_t(landscape, test_coords_t)
test_pred_spec <- fixture_pred_spec()
test_pred_meta_t <- fixture_pred_meta_t()
test_pred_data_t <- fixture_pred_data_t(landscape, test_coords_t)

usethis::use_data(
  test_coords_t,
  test_periods_t,
  test_lulc_meta_t,
  test_lulc_data_t,
  test_pred_spec,
  test_pred_meta_t,
  test_pred_data_t,
  overwrite = TRUE,
  compress = "xz",
  internal = TRUE
)
