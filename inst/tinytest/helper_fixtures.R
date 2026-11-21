# The test fixtures, built from one synthetic_landscape(): the grid of the landscape gives the
# coordinates, its land use the land use data, and the drivers of its change the predictors.
# data-raw/test_fixtures.R stores them in R/sysdata.rda; the constructor tests rebuild them
# here and compare, so that the stored data cannot go stale.

# the landscape: three observed periods of four classes, and the drivers of their change, on
# a 30 x 30 grid of 100 m cells in Swiss LV95
fixture_landscape <- function() {
  set.seed(42)
  synthetic_landscape(
    n_periods = 3L,
    n_cells = 30L,
    resolution = 100,
    xmin = 2697000,
    ymin = 1252000,
    crs = "EPSG:2056"
  )
}

# the coordinates: one point per cell of the landscape's grid
fixture_coords_t <- function(landscape) {
  grid <- landscape$lulc
  create_coords_t_square(
    epsg = as.integer(terra::crs(grid, describe = TRUE)$code),
    extent = terra::ext(grid),
    resolution = terra::res(grid)[1]
  )
}

# the landscape's periods 1 to 3 are observed, period 4 is extrapolated
fixture_periods_t <- function() {
  create_periods_t(
    period_length_str = "P10Y",
    start_observed = "1985-01-01",
    end_observed = "2005-01-01",
    end_extrapolated = "2015-01-01"
  )
}

# the synthetic classes, so that id_lulc equals the source class
fixture_lulc_meta_t <- function() {
  create_lulc_meta_t(list(
    forest = list(pretty_name = "Forest", src_classes = 1L),
    arable = list(pretty_name = "Arable land", src_classes = 2L),
    urban = list(pretty_name = "Urban areas", src_classes = 3L),
    lake = list(pretty_name = "Lake", src_classes = 4L)
  ))
}

fixture_lulc_data_t <- function(landscape, coords) {
  extract_using_coords_t(landscape$lulc, coords)[, .(
    id_run = 0L,
    id_coord,
    id_lulc = as.integer(value),
    id_period = as.integer(sub(".*=", "", layer))
  )] |>
    as_lulc_data_t()
}

# one predictor per data type; elevation and soil type are static, population and protection
# change from period to period
fixture_pred_spec <- function() {
  list(
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
      description = "Population count, per period",
      data_type = "int",
      sources = list(list(url = "https://example.com/pop.tif", md5sum = "def456"))
    ),
    is_protected = list(
      unit = "boolean",
      pretty_name = "Protected Area",
      description = "Is the area protected, per period",
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
}

fixture_pred_meta_t <- function() {
  create_pred_meta_t(fixture_pred_spec())
}

# the predictors, derived from the landscape so that they carry its change process:
# - elevation (static, period 0) follows site quality;
# - population (periods 1 to 3) follows accessibility, grows by a fifth each period, and adds
#   500 in cells that are urban in that period;
# - protection (periods 1 to 3) covers the least accessible cells, and its threshold rises
#   each period, so that the protected area grows;
# - soil type (static, period 0) is the nuisance field, which the process ignores.
fixture_pred_data_t <- function(landscape, coords) {
  drivers <- data.table::dcast(
    extract_using_coords_t(landscape$drivers, coords),
    id_coord ~ layer,
    value.var = "value"
  )
  lulc <- fixture_lulc_data_t(landscape, coords)
  dynamic <- lulc[drivers, on = "id_coord"]

  as_pred_data_t(rbind(
    drivers[, .(id_pred = 1L, id_period = 0L, id_coord, value = 400 + 400 * site_quality)],
    dynamic[, .(
      id_pred = 2L,
      id_period,
      id_coord,
      value = round(1000 * accessibility * (1 + 0.2 * (id_period - 1L))) +
        500 * (id_lulc == 3L)
    )],
    dynamic[, .(
      id_pred = 3L,
      id_period,
      id_coord,
      value = as.numeric(accessibility < 0.2 + 0.05 * id_period)
    )],
    drivers[, .(
      id_pred = 4L,
      id_period = 0L,
      id_coord,
      value = ceiling(5 * pmax(random_nuisance, 1e-9))
    )]
  )[, id_run := 0L])
}
