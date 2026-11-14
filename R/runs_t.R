#' Create Runs Table
#'
#' Creates a runs_t table that stores the identities, hierarchy, and description of
#' runs. Defaults to a table with one row with the base `id_run := 0`. Runs are a
#' feature that allows one to label data that is either overriding base data (e.g.
#' predictors being used for forecasting) or produced by a specific model experiment
#' (e.g. a given allocation of land use change given a set of transition potential
#' models and allocation parameters).
#'
#' A run also records how its results were made, so that they can be traced and replayed
#' from the database:
#' - `rng_seed` seeds R's random number generator for [allocate()]. It belongs to its run
#'   alone: a run without one draws a seed when it is first allocated, and records it.
#' - `alloc_args` holds the allocator and its arguments, as [allocate()] used them;
#'   `allocate()` without an `allocator` uses them.
#' - `model_selection_criterion` holds `select_score` and `select_maximize`, as
#'   [predict_trans_pot()] and [fit_full_models()] used them; without `select_score`, they
#'   use it.
#' - `id_period_eval_start` is the first period the run allocates itself, i.e. where its
#'   maps branch off its parent's.
#'
#' `alloc_args` and `model_selection_criterion` are inherited through the lineage: a run
#' without them uses those of its closest ancestor that has them. A value passed as an
#' argument is recorded on the run when the run's own row holds none.
#'
#' @name runs_t
#'
#' @param x A list or data.frame coercible to a data.table
#'
#' @return A data.table of class "runs_t" with columns:
#'   - `id_run`: Unique ID for each run; 0 is the base run
#'   - `parent_id_run`: The run this run inherits data from; `NA` for the base run
#'   - `description`: Human-readable description of the run
#'   - `rng_seed`: Seed of R's random number generator for allocation; `NA` for none yet
#'   - `alloc_args`: Allocator and its arguments, as a named list; `NULL` for none
#'   - `model_selection_criterion`: `select_score` and `select_maximize`, as a named list;
#'     `NULL` for none
#'   - `id_period_eval_start`: First period the run allocates itself; `NA` for none
#' @export
as_runs_t <- function(x) {
  if (missing(x)) {
    x <- data.table::data.table(
      id_run = 0L,
      parent_id_run = NA_integer_,
      description = "Base"
    )
  }

  data.table::setDT(x)
  # runs registered before these columns existed have none of them
  for (col in c("rng_seed", "id_period_eval_start")) {
    if (!col %in% names(x)) data.table::set(x, j = col, value = NA_integer_)
  }
  for (col in RUN_PARAM_MAPS) {
    # wrapped in list(): a bare list of NULLs would read as removing the column
    if (!col %in% names(x)) data.table::set(x, j = col, value = list(vector("list", nrow(x))))
  }

  x |>
    cast_dt_col("id_run", "int") |>
    cast_dt_col("parent_id_run", "int") |>
    cast_dt_col("rng_seed", "int") |>
    cast_dt_col("id_period_eval_start", "int")

  as_ducklake_db_t(
    x,
    class_name = "runs_t",
    key_cols = "id_run",
    map_cols = RUN_PARAM_MAPS
  )
}

#' @export
validate.runs_t <- function(x, ...) {
  NextMethod()

  data.table::setcolorder(
    x,
    c(
      "id_run",
      "parent_id_run",
      "description",
      "rng_seed",
      "alloc_args",
      "model_selection_criterion",
      "id_period_eval_start"
    )
  )

  stopifnot(
    "id_run is not integer" = is.integer(x[["id_run"]]),
    "parent_id_run is not integer" = is.integer(x[["parent_id_run"]]),
    "rng_seed is not integer" = is.integer(x[["rng_seed"]]),
    "id_period_eval_start is not integer" = is.integer(x[["id_period_eval_start"]]),
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

# the run parameters stored as MAP columns, see as_runs_t()
RUN_PARAM_MAPS <- c("alloc_args", "model_selection_criterion")

# Create runs_t with the base run, or bring an existing one up to date: add the columns
# an older database lacks, and the base run if it is missing. Never upserts the base run
# itself, which would reset what it records each time the database is opened.
ensure_runs_t <- function(self) {
  if (!"runs_t" %in% self$list_tables()) {
    return(self$commit(as_runs_t(), "runs_t", method = "overwrite"))
  }

  stored <- self$get_table_metadata("runs_t")
  if (!all(RUN_PARAM_MAPS %in% stored[["map_cols"]])) {
    # small table: rewriting it also rewrites the map columns into its metadata
    runs <- self$fetch("runs_t")
    self$commit(as_runs_t(data.table::copy(runs)), "runs_t", method = "overwrite")
  }

  if (self$row_count("runs_t", where = "id_run = 0") == 0L) {
    self$commit(as_runs_t(), "runs_t", method = "upsert")
  }

  invisible(NULL)
}

# A run's parameter: the run's own value, or for an inherited parameter, that of the closest
# ancestor that has one. NULL if there is none.
get_run_param <- function(self, id_run, param, inherit = param %in% RUN_PARAM_MAPS) {
  runs <- self$runs_t
  lineage <- if (inherit) get_lineage(runs, as.integer(id_run)) else as.integer(id_run)
  for (run in lineage) {
    value <- runs[[param]][runs[["id_run"]] == run]
    if (is.list(value)) {
      value <- value[[1L]]
    }
    if (length(value) > 0L && !all(is.na(value))) {
      return(value)
    }
  }
  NULL
}

# Record parameters on a run's own row, without touching its other columns. Map parameters
# are given as named lists.
set_run_params <- function(self, id_run, ...) {
  params <- list(...)
  row <- data.table::data.table(id_run = as.integer(id_run))
  for (param in names(params)) {
    value <- params[[param]]
    data.table::set(row, j = param, value = if (is.list(value)) list(value) else value)
  }
  data.table::setattr(row, "key_cols", "id_run")
  data.table::setattr(row, "map_cols", intersect(names(params), RUN_PARAM_MAPS))
  self$commit(row, "runs_t", method = "upsert")
  invisible(NULL)
}

# Run parameters as stored in a MAP column: one string per argument, vectors joined by ","
encode_run_args <- function(args) {
  args <- Filter(Negate(is.null), args)
  lapply(args, function(value) paste(as.character(value), collapse = ","))
}

decode_run_args <- function(args) {
  lapply(args, function(value) {
    utils::type.convert(strsplit(as.character(value), ",", fixed = TRUE)[[1L]], as.is = TRUE)
  })
}

# The model selection criterion to use for `id_run`: the arguments when given, otherwise
# what the run records (inherited). Recorded on the run when its own row has none.
resolve_model_selection <- function(self, id_run, select_score, select_maximize) {
  if (is.null(select_score)) {
    recorded <- get_run_param(self, id_run, "model_selection_criterion")
    if (is.null(recorded)) {
      stop(
        glue::glue(
          "Run {id_run} records no model selection criterion; pass select_score and ",
          "select_maximize"
        ),
        call. = FALSE
      )
    }
    return(list(
      select_score = as.character(recorded[["select_score"]]),
      select_maximize = as.logical(recorded[["select_maximize"]])
    ))
  }

  stopifnot(
    "select_score must be a character string" = is_scalar_character(select_score),
    "select_maximize must be TRUE or FALSE" = is_scalar_logical(select_maximize)
  )
  criterion <- list(select_score = select_score, select_maximize = select_maximize)
  if (is.null(get_run_param(self, id_run, "model_selection_criterion", inherit = FALSE))) {
    set_run_params(self, id_run, model_selection_criterion = encode_run_args(criterion))
  }
  criterion
}
