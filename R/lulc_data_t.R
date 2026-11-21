LULC_DATA_T <- table_spec(
  "lulc_data_t",
  columns = c(id_run = "int", id_coord = "int", id_lulc = "int", id_period = "int"),
  descriptions = c(
    id_run = "Foreign key to [runs_t]",
    id_coord = "Foreign key to [coords_t]",
    id_lulc = "Foreign key to [lulc_meta_t]",
    id_period = "Foreign key to [periods_t]"
  ),
  key_cols = c("id_run", "id_coord", "id_period"),
  partition_cols = "id_run",
  # run-scoped: read through the lineage, one slice per period
  inherit_by = "id_period",
  allow_extra_cols = FALSE
)

#' Create LULC Data Table
#'
#' Creates an `lulc_data_t` table and validates that it matches the schema.
#'
#' @name lulc_data_t
#' @include ducklake_db_utils.R
#'
#' @param x An object that can be passed to [data.table::setDT()]
#'
#' @eval table_spec_rd(LULC_DATA_T)
#'
#' @export
as_lulc_data_t <- function(x) {
  as_table_t(x, LULC_DATA_T)
}

#' @export
#' @describeIn lulc_data_t Print an lulc_data_t object, passing params to data.table print
#' @param nrow see [data.table::print.data.table]
#' @param ... passed to [data.table::print.data.table]
print.lulc_data_t <- function(x, nrow = 10, ...) {
  if (nrow(x) > 0) {
    n_coords <- data.table::uniqueN(x[["id_coord"]])
    n_lulc <- data.table::uniqueN(x[["id_lulc"]])
    n_periods <- data.table::uniqueN(x[["id_period"]])

    cat(glue::glue(
      "LULC Data Table\n",
      "Observations: {nrow(x)}\n",
      "Coordinates: {n_coords}, LULC classes: {n_lulc}, Periods: {n_periods}\n\n"
    ))
  } else {
    cat("LULC Data Table (empty)\n")
  }
  NextMethod(nrow = nrow, ...)
  invisible(x)
}

#' @describeIn lulc_data_t Retrieve LULC data as SpatRaster
#' @param id_period Integer or vector of integers specifying which periods to
#' retrieve; returned as layers.
#' @keywords internal
#' @param id_run Integer; the run whose land use is read, through its lineage
lulc_data_as_rast <- function(self, private, id_period = NULL, id_run = 0L) {
  # Build query to join lulc_data_t with coords_t
  where_clause <- NULL
  if (!is.null(id_period)) {
    where_clause <- private$sql("id_period in ({id_period*})")
  }

  data <- self$fetch("lulc_data_t", where = where_clause, id_run = id_run)

  if (nrow(data) == 0L) {
    stop("No LULC data found for the specified periods")
  }

  tabular_to_raster(
    data = data,
    coords = self$coords_minimal,
    value_col = "id_lulc"
  )
}
