#' Registry of Table Specifications
#'
#' `EVOLAND_SCHEMA` is a named list of the [table_spec] of every table of [evoland_db]. The
#' ER diagram and the column tables in `devdocs/design/database.md` are generated from it by
#' `devdocs/er_diagram.R`, and [get_evoland_db_read_expr()] reads which tables are
#' run-scoped from it.
#'
#' @name EVOLAND_SCHEMA
#' @keywords internal
#' @include alloc_params_t.R coords_t.R lulc_data_t.R lulc_meta_t.R neighbors_t.R periods_t.R
#' @include pred_data_t.R pred_meta_t.R reporting_t.R runs_t.R trans_meta_t.R trans_models_t.R
#' @include trans_pot_t.R trans_preds_t.R trans_rates_t.R trans_viability_t.R
EVOLAND_SCHEMA <- list(
  ALLOC_PARAMS_T,
  COORDS_T,
  LULC_DATA_T,
  LULC_META_T,
  NEIGHBORS_T,
  PERIODS_T,
  PRED_DATA_T,
  PRED_META_T,
  REPORTING_T,
  RUNS_T,
  TRANS_META_T,
  TRANS_MODELS_T,
  TRANS_POT_T,
  TRANS_PREDS_T,
  TRANS_RATES_T,
  TRANS_VIABILITY_T
)
names(EVOLAND_SCHEMA) <- vapply(EVOLAND_SCHEMA, \(spec) spec[["name"]], character(1))
