TRANS_VIABILITY_T <- table_spec(
  "trans_viability_t",
  columns = c(
    id_run = "int",
    id_trans = "int",
    is_viable = "bool"
  ),
  descriptions = c(
    id_run = "Foreign key to [runs_t]",
    id_trans = "Foreign key to [trans_meta_t]",
    is_viable = paste(
      "Whether the transition is viable for modelling in this run (and the runs inheriting",
      "from it)"
    )
  ),
  key_cols = c("id_run", "id_trans"),
  inherit_by = "id_trans"
)

#' Create Transition Viability Table
#'
#' Creates a trans_viability_t table, which says whether a transition of [trans_meta_t] is
#' viable for modelling in a run. A transition can be non-viable for conceptual or statistical
#' reasons: transitions out of an "other" class say nothing about what comes next
#' (`exclude_anterior`), and a transition may be observed too rarely for statistical inference
#' (`min_cardinality_abs`, `min_frequency_rel`). Which transitions are modelled is also a
#' modelling choice that runs may vary, e.g. to compare model performance with and without a
#' rare transition.
#'
#' The default viability is stored as rows of the root run, 0, usually derived with
#' `create_trans_viability_t()`. A run that allows or disallows transitions differently
#' stores a row per transition it changes. Like every run-dependent table, the table is read
#' through the run's lineage, per transition: a transition takes the row of the closest run in
#' the lineage that has one, and a transition without any row is viable. A run's resolved set
#' of transitions is returned by `db$trans_meta_v(id_run)`, and it is the set that the run's
#' predictor selection, model fitting, prediction, allocation parameters and allocation work
#' on. A run that branches off its parent (see [allocate()]) and holds rows of its own
#' predicts its own transition potentials instead of taking its parent's.
#'
#' ```r
#' # the defaults, as rows of run 0
#' db$trans_viability_t <- create_trans_viability_t(db$trans_meta_t, min_cardinality_abs = 5L)
#' # run 1 also disallows transition 3
#' db$trans_viability_t <- as_trans_viability_t(list(id_run = 1L, id_trans = 3L, is_viable = FALSE))
#' ```
#'
#' @name trans_viability_t
#' @include ducklake_db_utils.R
#'
#' @param x A list or data.frame coercible to a data.table
#' @param trans_meta A [trans_meta_t] table
#' @param min_cardinality_abs Minimum absolute number of transitions for viability (optional)
#' @param min_frequency_rel Minimum relative frequency of transitions for viability (optional)
#' @param exclude_anterior Vector of id_lulc values to exclude as anterior (source) classes
#' @param exclude_posterior Vector of id_lulc values to exclude as posterior (target) classes
#' @param id_run Integer, the run the rows are for; 0, the root run, for the defaults
#'
#' @eval table_spec_rd(TRANS_VIABILITY_T)
#' @export
as_trans_viability_t <- function(x) {
  as_table_t(x, TRANS_VIABILITY_T)
}

#' @describeIn trans_viability_t Mark each transition of `trans_meta` as viable unless it is
#' excluded or falls below the thresholds; with both thresholds, a transition is non-viable
#' only if it falls below both.
#' @export
create_trans_viability_t <- function(
  trans_meta,
  min_cardinality_abs = NULL,
  min_frequency_rel = NULL,
  exclude_anterior = NULL,
  exclude_posterior = NULL,
  id_run = 0L
) {
  stopifnot(
    "trans_meta must be a trans_meta_t" = inherits(trans_meta, "trans_meta_t"),
    "id_run must be a single integer" = is_scalar_integerish(id_run)
  )
  below_cardinality <- if (is.null(min_cardinality_abs)) {
    NULL
  } else {
    trans_meta[["cardinality"]] < min_cardinality_abs
  }
  below_frequency <- if (is.null(min_frequency_rel)) {
    NULL
  } else {
    trans_meta[["frequency_rel"]] < min_frequency_rel
  }
  # with both thresholds given, meeting either one suffices
  thresholds <- Filter(Negate(is.null), list(below_cardinality, below_frequency))
  below_threshold <- if (length(thresholds) > 0L) Reduce(`&`, thresholds) else FALSE

  data.table::data.table(
    id_run = id_run,
    id_trans = trans_meta[["id_trans"]],
    is_viable = !(trans_meta[["id_lulc_anterior"]] %in%
      exclude_anterior |
      trans_meta[["id_lulc_posterior"]] %in% exclude_posterior |
      below_threshold)
  ) |>
    as_trans_viability_t()
}

#' @export
validate.trans_viability_t <- function(x, ...) {
  NextMethod()

  stopifnot("is_viable must be TRUE or FALSE" = !anyNA(x[["is_viable"]]))

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

