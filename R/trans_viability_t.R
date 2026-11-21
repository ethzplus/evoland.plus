#' Create Transition Viability Table
#'
#' Creates a trans_viability_t table, which overrides the viability of transitions for a run.
#' [trans_meta_t] identifies the transitions and gives each a default viability, `is_viable`.
#' A run that allows or disallows transitions differently, e.g. to compare model performance
#' with and without a rare transition, stores a row per transition it changes here. Like every
#' run-dependent table, it is read through the run's lineage: a run without rows of its own
#' takes those of its closest ancestor that has them, and a transition without an override in
#' the lineage keeps the default of [trans_meta_t]. A run's resolved set of transitions is
#' returned by `db$trans_meta_v(id_run)`, and it is the set that the run's predictor
#' selection, model fitting, prediction, allocation parameters and allocation work on.
#'
#' A run that branches off its parent (see [allocate()]) and overrides viability predicts its
#' own transition potentials instead of taking its parent's.
#'
#' To derive a run's overrides from other thresholds, take `is_viable` of
#' [create_trans_meta_t()] for the same transitions, e.g.
#' `create_trans_meta_t(db$trans_v(), min_cardinality_abs = 20L)[, .(id_run = 1L, id_trans,
#' is_viable)]`.
#'
#' @name trans_viability_t
#'
#' @param x A list or data.frame coercible to a data.table
#'
#' @return A data.table of class "trans_viability_t" with columns:
#'   - `id_run`: Foreign key to runs_t
#'   - `id_trans`: Foreign key to trans_meta_t
#'   - `is_viable`: Whether the transition is viable for modelling in this run (and the runs
#'     inheriting from it)
#' @export
as_trans_viability_t <- function(x) {
  if (missing(x)) {
    x <- data.table::data.table(
      id_run = integer(0),
      id_trans = integer(0),
      is_viable = logical(0)
    )
  }

  data.table::setDT(x) |>
    cast_dt_col("id_run", "int") |>
    cast_dt_col("id_trans", "int") |>
    cast_dt_col("is_viable", "bool")

  as_ducklake_db_t(
    x,
    class_name = "trans_viability_t",
    key_cols = c("id_run", "id_trans")
  )
}

#' @export
validate.trans_viability_t <- function(x, ...) {
  NextMethod()

  data.table::setcolorder(x, c("id_run", "id_trans", "is_viable"))

  stopifnot(
    "id_run is not integer" = is.integer(x[["id_run"]]),
    "id_trans is not integer" = is.integer(x[["id_trans"]]),
    "is_viable must be TRUE or FALSE" = is.logical(x[["is_viable"]]) && !anyNA(x[["is_viable"]])
  )

  return(x)
}

#' @describeIn trans_viability_t Print a trans_viability_t object, passing params to
#' data.table print
#' @param nrow see [data.table::print.data.table]
#' @param ... passed to [data.table::print.data.table]
#' @export
print.trans_viability_t <- function(x, nrow = 10, ...) {
  if (nrow(x) > 0) {
    cat(glue::glue(
      "Transition Viability Overrides\n",
      "Runs: {data.table::uniqueN(x[['id_run']])}, ",
      "allowed: {sum(x[['is_viable']])}, disallowed: {sum(!x[['is_viable']])}\n\n"
    ))
  } else {
    cat("Transition Viability Overrides (empty)\n")
  }
  NextMethod(nrow = nrow, ...)
  invisible(x)
}
