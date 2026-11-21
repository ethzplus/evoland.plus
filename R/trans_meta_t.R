TRANS_META_T <- table_spec(
  "trans_meta_t",
  columns = c(
    id_trans = "int",
    id_lulc_anterior = "int",
    id_lulc_posterior = "int",
    cardinality = "int",
    frequency_rel = "float",
    frequency_abs = "float"
  ),
  descriptions = c(
    id_trans = "Unique ID for each transition",
    id_lulc_anterior = "Foreign key to [lulc_meta_t], the class before the transition",
    id_lulc_posterior = "Foreign key to [lulc_meta_t], the class after the transition",
    cardinality = "How many times this transition was observed",
    frequency_rel = "Frequency relative to all observed transitions",
    frequency_abs = "Frequency relative to all coordinate points"
  ),
  key_cols = c("id_lulc_anterior", "id_lulc_posterior"),
  alternate_key_cols = "id_trans"
)

#' Create Transition Metadata Table
#'
#' Creates a trans_meta_t table based on observed land use transitions in the
#' LULC data: one row per observed transition type, with how often it was observed.
#' The table is the same for all runs: `id_trans` identifies a pair of anterior and
#' posterior classes in every run. Whether a transition is viable for modelling is
#' stored per run in [trans_viability_t]; `db$trans_meta_v(id_run)` returns this table
#' with the viability resolved for a run, in a column `is_viable`.
#'
#' The usual workflow derives both tables from the observed transitions:
#'
#' ```r
#' db$trans_meta_t <- create_trans_meta_t(db$trans_v())
#' db$trans_viability_t <- create_trans_viability_t(db$trans_meta_t, min_cardinality_abs = 5L)
#' ```
#'
#' @name trans_meta_t
#' @include table_spec.R
#'
#' @param x A list or data.frame coercible to a trans_meta_t object.
#' @param transitions A trans_v table, with columns id_coord, id_lulc_anterior,
#' id_lulc_posterior, id_period; contains transitions and non-transitions (where
#' anterior == posterior)
#' @param ... Not used. The thresholds and exclusions that `create_trans_meta_t()` used to
#' take are arguments of [create_trans_viability_t()].
#'
#' @eval table_spec_rd(TRANS_META_T)
#' @export
as_trans_meta_t <- function(x) {
  as_table_t(x, TRANS_META_T)
}

#' @describeIn trans_meta_t Calculate the transition metadata from observed transitions.
#' `id_trans` numbers the transitions in the order of their classes.
#' @export
create_trans_meta_t <- function(transitions, ...) {
  if (...length() > 0L) {
    stop(
      "create_trans_meta_t() no longer derives viability; pass ",
      toString(names(list(...))),
      " to create_trans_viability_t() instead",
      call. = FALSE
    )
  }
  # Return empty table if no data provided
  if (missing(transitions) || nrow(transitions) == 0L) {
    return(as_trans_meta_t())
  }

  # Count total coordinate-period pairs that could have transitions
  n_total_pairs <- nrow(transitions)

  # Aggregate transitions by type and calculate frequencies
  transitions[
    # excludes non-transitions
    id_lulc_anterior != id_lulc_posterior
  ][,
    .(cardinality = .N),
    by = .(id_lulc_anterior, id_lulc_posterior)
  ][,
    `:=`(
      frequency_rel = cardinality / sum(cardinality),
      frequency_abs = cardinality / n_total_pairs
    )
  ][
    order(id_lulc_anterior, id_lulc_posterior) # deterministic id_trans
  ][,
    id_trans := seq_len(.N)
  ] |>
    as_trans_meta_t()
}

#' @export
validate.trans_meta_t <- function(x, ...) {
  NextMethod()
  validate_table_spec(x, TRANS_META_T)

  stopifnot(
    all(x[["cardinality"]] >= 0),
    all(x[["frequency_rel"]] >= 0 & x[["frequency_rel"]] <= 1),
    all(x[["frequency_abs"]] >= 0 & x[["frequency_abs"]] <= 1)
  )

  return(x)
}

#' @export
#' @describeIn trans_meta_t Print a trans_meta_t object, passing params to data.table print
#' @param nrow see [data.table::print.data.table]
#' @param ... passed to [data.table::print.data.table]
print.trans_meta_t <- function(x, nrow = 10, ...) {
  if (nrow(x) > 0) {
    cat(glue::glue(
      "Transition Metadata Table\n",
      "Total transition types: {nrow(x)}\n",
      # only as returned by db$trans_meta_v(), resolved for a run
      if ("is_viable" %in% names(x)) "Viable for modelling: {sum(x[['is_viable']])}\n" else "",
      "Total observed transitions: {sum(x[['cardinality']])}\n\n"
    ))
  } else {
    cat("Transition Metadata Table (empty)\n")
  }
  NextMethod(nrow = nrow, ...)
  invisible(x)
}

# trans_meta_t with is_viable resolved for `id_run`: the row in trans_viability_t of the
# closest run in the lineage that has one; a transition without any is viable. A database
# written before trans_viability_t existed, and not yet migrated (i.e. opened read-only), still
# carries is_viable in trans_meta_t, which then serves as the root run's rows.
trans_meta_read_expr <- function(self, id_run) {
  meta_read_expr <- self$get_read_expr("trans_meta_t")
  legacy <- has_legacy_is_viable(self)
  viability_read_expr <- if ("trans_viability_t" %in% self$list_tables()) {
    self$get_read_expr("trans_viability_t", id_run = id_run)
  } else {
    DBI::SQL("(select null::integer as id_trans, null::boolean as is_viable where false)")
  }
  glue::glue_sql(
    "(select m.* {exclude}, coalesce(v.is_viable, {fallback}true) as is_viable
      from {meta_read_expr} m
      left join {viability_read_expr} v using (id_trans))",
    exclude = DBI::SQL(if (legacy) "exclude (is_viable)" else ""),
    fallback = DBI::SQL(if (legacy) "m.is_viable, " else ""),
    .con = self$connection
  )
}

# Whether trans_meta_t still has the is_viable column it had before trans_viability_t
has_legacy_is_viable <- function(self) {
  "trans_meta_t" %in%
    self$list_tables() &&
    "is_viable" %in%
      self$get_query("select column_name from (describe dl_db.trans_meta_t)", as_atomic = TRUE)
}

#' @describeIn trans_meta_t The transition metadata with `is_viable` resolved for a run,
#' see [trans_viability_t]
#' @param self An [evoland_db] instance
#' @param id_run Integer, the run whose viability overrides are read, through its lineage
trans_meta_v <- function(self, id_run = 0L) {
  stopifnot("id_run must be a single integer" = is_scalar_integerish(id_run))
  self$get_query("from {trans_meta_read_expr(self, id_run)} order by id_trans") |>
    as_trans_meta_t()
}

# The viable transitions of a run, optionally narrowed to `id_trans`
viable_trans_meta <- function(self, id_trans = NULL, id_run = 0L) {
  stopifnot("id_trans must be NULL or integers" = is.null(id_trans) || is_integerish(id_trans))
  viable <- trans_meta_v(self, id_run)[is_viable == TRUE]
  if (is.null(id_trans)) {
    return(viable)
  }
  # env renames the argument, which inside `[` would be masked by the column of the same name
  viable[id_trans %in% ids, env = list(ids = I(id_trans))]
}
