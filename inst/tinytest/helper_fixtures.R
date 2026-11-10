# The synthetic test fixtures, built from synthetic_landscape() on the grid of test_coords_t:
# land use with a known change process, and predictors derived from the drivers of that
# process. data-raw/test_fixtures.R stores them in R/sysdata.rda; the constructor tests rebuild
# them here and compare, so that the stored data cannot go stale.

# the landscape: three observed periods of four classes, and the drivers of their change
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

# static predictors (period 0), one per data type, derived from the drivers so that they carry
# the change process: elevation follows site quality, population and protection follow
# accessibility; the soil type is the nuisance field, which the process ignores
fixture_pred_data_t <- function(landscape, coords) {
  drivers <- data.table::dcast(
    extract_using_coords_t(landscape$drivers, coords),
    id_coord ~ layer,
    value.var = "value"
  )
  as_pred_data_t(rbind(
    drivers[, .(id_pred = 1L, id_coord, value = 400 + 400 * site_quality)],
    drivers[, .(id_pred = 2L, id_coord, value = round(1000 * accessibility))],
    drivers[, .(id_pred = 3L, id_coord, value = as.numeric(accessibility < 0.3))],
    drivers[, .(id_pred = 4L, id_coord, value = ceiling(5 * pmax(random_nuisance, 1e-9)))]
  )[, `:=`(id_run = 0L, id_period = 0L)])
}
