# Creates canonical test fixture data for use in tests, examples, and vignettes.
#
# This script is not run during checks. Each fixture below is rebuilt in the
# corresponding inst/tinytest/test_*.R file (the synthetic ones by helper_fixtures.R),
# which compares the freshly built object against the one stored in R/sysdata.rda. Changing a fixture here therefore requires
# rerunning this script, so that the stored data cannot go stale relative to the
# constructors it is built from.
library("evoland.plus")
library("data.table")

# Coordinates: 30×30 grid in Swiss LV95
test_coords_t <- create_coords_t_square(
  epsg = 2056,
  extent = terra::ext(c(xmin = 2697000, xmax = 2700000, ymin = 1252000, ymax = 1255000)),
  resolution = 100
)

# Time periods
test_periods_t <- create_periods_t(
  period_length_str = "P10Y",
  start_observed = "1985-01-01",
  end_observed = "2005-01-01",
  end_extrapolated = "2015-01-01"
)

# Land use with a known change process, and predictors derived from its drivers, on the same
# grid; built by inst/tinytest/helper_fixtures.R, which the constructor tests use to check the
# stored data
source("inst/tinytest/helper_fixtures.R")
landscape <- fixture_landscape()
test_lulc_meta_t <- fixture_lulc_meta_t()
test_lulc_data_t <- fixture_lulc_data_t(landscape, test_coords_t)

# Predictor metadata (one for each data type)
test_pred_spec <- list(
  elevation = list(
    unit = "m",
    pretty_name = "Elevation",
    description = "Elevation above sea level",
    data_type = "float",
    sources = list(list(url = "https://example.com/elevation.tif", md5sum = "abc123"))
  ),
  population = list(
    unit = "count",
    pretty_name = "Population",
    description = "Population count",
    data_type = "int",
    sources = list(list(url = "https://example.com/pop.tif", md5sum = "def456"))
  ),
  is_protected = list(
    unit = "boolean",
    pretty_name = "Protected Area",
    description = "Is the area protected",
    data_type = "bool",
    sources = list(list(url = "https://example.com/protected.gpkg", md5sum = "ghi789"))
  ),
  soil_type = list(
    unit = "category",
    pretty_name = "Soil Type",
    description = "Type of soil",
    data_type = "factor",
    factor_levels = c("Sandy", "Loamy", "Clay", "Peaty", "Chalky"),
    sources = list(list(url = "https://example.com/soil.tif", md5sum = "jkl012"))
  )
)
test_pred_meta_t <- create_pred_meta_t(test_pred_spec)

test_pred_data_t <- fixture_pred_data_t(landscape, test_coords_t)

# Save all fixtures
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
