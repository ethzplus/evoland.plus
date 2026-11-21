TRANS_META_T <- table_spec(
  "trans_meta_t",
  columns = c(
    id_trans = "int",
    id_lulc_anterior = "int",
    id_lulc_posterior = "int",
    cardinality = "int",
    frequency_rel = "float",
    frequency_abs = "float",
    is_viable = "bool"
  ),
  descriptions = c(
    id_trans = "Unique ID for each transition",
    id_lulc_anterior = "Foreign key to [lulc_meta_t], the class before the transition",
    id_lulc_posterior = "Foreign key to [lulc_meta_t], the class after the transition",
    cardinality = "How many times this transition was observed",
    frequency_rel = "Frequency relative to all observed transitions",
    frequency_abs = "Frequency relative to all coordinate points",
    is_viable = "Whether a transition potential model is fitted for this transition"
  ),
  key_cols = c("id_lulc_anterior", "id_lulc_posterior"),
  alternate_key_cols = "id_trans"
)

#' Create Transition Metadata Table
#'
#' Creates a trans_meta_t table based on observed land use transitions in the
#' LULC data. This function analyzes transition patterns and creates metadata
#' entries for each viable transition type. The metadata is kept constant across
#' runs: `id_trans` identifies a pair of anterior and posterior classes in every run.
#' `is_viable` is the default viability; a run may allow or disallow transitions
#' differently through [trans_viability_t], and `db$trans_meta_v(id_run)` returns the
#' table with the viability resolved for a run.
#'
#' @name trans_meta_t
#' @include table_spec.R
#'
#' @param x A list or data.frame coercible to a trans_meta_t object.
#' @param transitions A trans_v table, with columns id_coord, id_lulc_anterior,
#' id_lulc_posterior, id_period; contains transitions and non-transitions (where
#' anterior == posterior)
#' @param min_cardinality_abs Minimum absolute number of transitions for viability (optional)
#' @param min_frequency_rel Minimum relative frequency of transitions for viability (optional)
#' @param exclude_anterior Vector of id_lulc values to exclude as anterior (source) classes
#' @param exclude_posterior Vector of id_lulc values to exclude as posterior (target) classes
#'
#' @eval table_spec_rd(TRANS_META_T)
#' @export
as_trans_meta_t <- function(x) {
  as_table_t(x, TRANS_META_T)
}

#' @describeIn trans_meta_t Calculate the transition metadata and mark for modelling
#' feasibility. Does not attribute `id_trans`; this only makes sense as part of a DB.
#' @export
create_trans_meta_t <- function(
  transitions,
  min_cardinality_abs = NULL,
  min_frequency_rel = NULL,
  exclude_anterior = NULL,
  exclude_posterior = NULL
) {
  # Return empty table if no data provided
  if (missing(transitions) || nrow(transitions) == 0L) {
    return(as_trans_meta_t())
  }

  # Count total coordinate-period pairs that could have transitions
  n_total_pairs <- nrow(transitions)

  # Aggregate transitions by type and calculate frequencies
  trans_summary <-
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
    ]

  # default: all viable until determined otherwise
  trans_summary[, is_viable := TRUE]
  trans_summary[, id_trans := seq_len(.N)]

  # Apply exclusions
  if (!is.null(exclude_anterior) && length(exclude_anterior) > 0) {
    trans_summary[id_lulc_anterior %in% exclude_anterior, is_viable := FALSE]
  }

  if (!is.null(exclude_posterior) && length(exclude_posterior) > 0) {
    trans_summary[id_lulc_posterior %in% exclude_posterior, is_viable := FALSE]
  }

  # Apply minimum thresholds
  if (!is.null(min_cardinality_abs) || !is.null(min_frequency_rel)) {
    if (!is.null(min_cardinality_abs) && !is.null(min_frequency_rel)) {
      # Both specified: viable if both thresholds are met
      trans_summary[
        cardinality < min_cardinality_abs &
          frequency_rel < min_frequency_rel,
        is_viable := FALSE
      ]
    } else if (!is.null(min_cardinality_abs)) {
      # Only absolute threshold specified
      trans_summary[
        cardinality < min_cardinality_abs,
        is_viable := FALSE
      ]
    } else {
      # Only relative threshold specified
      trans_summary[
        frequency_rel < min_frequency_rel,
        is_viable := FALSE
      ]
    }
  }

  as_trans_meta_t(trans_summary)
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
    n_viable <- sum(x[["is_viable"]])
    n_total <- nrow(x)
    total_transitions <- sum(x[["cardinality"]])

    cat(glue::glue(
      "Transition Metadata Table\n",
      "Total transition types: {n_total}\n",
      "Viable for modelling: {n_viable}\n",
      "Total observed transitions: {total_transitions}\n\n"
    ))
  } else {
    cat("Transition Metadata Table (empty)\n")
  }
  NextMethod(nrow = nrow, ...)
  invisible(x)
}

# trans_meta_t with is_viable resolved for `id_run`: the override in trans_viability_t of the
# closest run in the lineage that has one, else the default. A database without
# trans_viability_t (e.g. one written before it existed) reads the defaults as stored.
trans_meta_read_expr <- function(self, id_run) {
  meta_read_expr <- self$get_read_expr("trans_meta_t")
  if (!"trans_viability_t" %in% self$list_tables()) {
    return(meta_read_expr)
  }
  glue::glue_sql(
    "(select m.* replace (coalesce(v.is_viable, m.is_viable) as is_viable)
      from {meta_read_expr} m
      left join {viability_read_expr} v using (id_trans))",
    viability_read_expr = self$get_read_expr("trans_viability_t", id_run = id_run),
    .con = self$connection
  )
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
