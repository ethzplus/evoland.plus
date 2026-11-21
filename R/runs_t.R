RUNS_T <- table_spec(
  "runs_t",
  columns = c(
    id_run = "int",
    parent_id_run = "int",
    description = "char",
    rng_seed = "int",
    id_period_eval_start = "int"
  ),
  descriptions = c(
    id_run = "Unique ID for each run; 0 is the base run",
    parent_id_run = "The run this run inherits data from; `NA` for the base run",
    description = "Human-readable description of the run",
    rng_seed = "Seed for the run's allocation, see [allocate()]; `NA` for none",
    id_period_eval_start = paste(
      "First period the run allocates itself, i.e. where it branches off its parent, see",
      "[allocate()]; `NA` for none"
    )
  ),
  key_cols = "id_run"
)

#' Create Runs Table
#'
#' Creates a runs_t table that stores the identities, hierarchy, and description of
#' runs. Defaults to a table with one row with the base `id_run := 0`. Runs are a
#' feature that allows one to label data that is either overriding base data (e.g.
#' predictors being used for forecasting) or produced by a specific model experiment
#' (e.g. a given allocation of land use change given a set of transition potential
#' models and allocation parameters).
#'
#' Like any evoland.plus table, `runs_t` takes columns beyond the ones listed below, e.g. to
#' document how a run was made (its allocator, model selection, ...): add them to the table
#' you assign or commit.
#'
#' @name runs_t
#' @include table_spec.R
#'
#' @param x A list or data.frame coercible to a data.table
#'
#' @eval table_spec_rd(RUNS_T)
#' @export
as_runs_t <- function(x) {
  if (missing(x)) {
    x <- data.table::data.table(
      id_run = 0L,
      parent_id_run = NA_integer_,
      description = "Base"
    )
  }

  as_table_t(x, RUNS_T)
}

#' @export
validate.runs_t <- function(x, ...) {
  NextMethod()
  validate_table_spec(x, RUNS_T)

  stopifnot(
    "all parent_id_run must be in id_run or NA" = all(
      is.na(x[["parent_id_run"]]) | x[["parent_id_run"]] %in% x[["id_run"]]
    ),
    "no base (0) id_run" = 0L %in% x[["id_run"]]
  )

  return(x)
}

#' @describeIn runs_t Print a runs_t object, passing params to data.table print
#' @param nrow see [data.table::print.data.table]
#' @param ... passed to [data.table::print.data.table]
#' @export
print.runs_t <- function(x, nrow = 10, ...) {
  if (nrow(x) > 0) {
    # could add a small recursive function to look up max hierarchy depth
    cat(glue::glue(
      "Run metadata table\n",
      "Number of runs: {nrow(x)}\n\n"
    ))
  } else {
    cat("Runs Table (empty)\n")
  }
  NextMethod(nrow = nrow, ...)
  invisible(x)
}

# lineage is ordered from most recent to oldest (i.e. id_run, parent_id_run,
# grandparent_id_run, etc.)
get_lineage <- function(runs_t, id_run) {
  lineage <- integer(0)
  current_id <- id_run

  if (nrow(runs_t[id_run == current_id]) == 0L) {
    stop(
      glue::glue(
        "Requested run (id_run = {id_run}) not found in runs_t"
      ),
      call. = FALSE
    )
  }

  repeat {
    parent_id <- runs_t[id_run == current_id, parent_id_run]

    if (length(parent_id) == 0L) {
      stop("parent not found, should never happen with valid runs_t")
    }

    # append to lineage
    lineage <- c(lineage, current_id)
    if (is.na(parent_id)) {
      # encountered root
      break
    }
    current_id <- parent_id
  }

  lineage
}
