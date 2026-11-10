# `evoland_db` - R6 Class for lulc data management

An R6 class that provides an interface to a folder-based data storage
system for the evoland.plus package. This class uses DuckDB for
in-memory SQL operations while persisting data to a DuckLake catalog on
disk.

Inherits from
[ducklake_db](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.md)
for generic database operations.

## See also

Additional methods and active bindings are added to this class in
separate files:

- [evoland_db_views](https://ethzplus.github.io/evoland.plus/reference/evoland_db_views.md) -
  View active bindings (lulc_meta_long_v, etc.) and methods

## Super class

[`ducklake_db`](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.md)
-\> `evoland_db`

## Active bindings

- `coords_t`:

  Get or upsert
  [coords_t](https://ethzplus.github.io/evoland.plus/reference/coords_t.md)

- `periods_t`:

  Get or upsert
  [periods_t](https://ethzplus.github.io/evoland.plus/reference/periods_t.md)

- `lulc_meta_t`:

  Get or upsert
  [lulc_meta_t](https://ethzplus.github.io/evoland.plus/reference/lulc_meta_t.md)

- `lulc_data_t`:

  Get or upsert
  [lulc_data_t](https://ethzplus.github.io/evoland.plus/reference/lulc_data_t.md)

- `pred_data_t`:

  Get or upsert
  [pred_data_t](https://ethzplus.github.io/evoland.plus/reference/pred_data_t.md)

- `pred_meta_t`:

  Get or upsert
  [pred_meta_t](https://ethzplus.github.io/evoland.plus/reference/pred_meta_t.md)

- `trans_meta_t`:

  Get or upsert
  [trans_meta_t](https://ethzplus.github.io/evoland.plus/reference/trans_meta_t.md)

- `trans_preds_t`:

  Get or upsert
  [trans_preds_t](https://ethzplus.github.io/evoland.plus/reference/trans_preds_t.md);
  deselect a predictor by setting its `is_selected` to `FALSE`

- `trans_rates_t`:

  Get or upsert
  [trans_rates_t](https://ethzplus.github.io/evoland.plus/reference/trans_rates_t.md)

- `trans_models_t`:

  Get or upsert
  [trans_models_t](https://ethzplus.github.io/evoland.plus/reference/trans_models_t.md)

- `alloc_params_t`:

  Get or upsert
  [alloc_params_t](https://ethzplus.github.io/evoland.plus/reference/alloc_params_t.md)

- `trans_pot_t`:

  Get or upsert raw transition potentials
  [trans_pot_t](https://ethzplus.github.io/evoland.plus/reference/trans_pot_t.md).
  These are per-transition model probabilities stored by
  [`predict_trans_pot()`](https://ethzplus.github.io/evoland.plus/reference/trans_pot_t.md).
  Use
  [`adjusted_trans_pot_v()`](https://ethzplus.github.io/evoland.plus/reference/evoland_db_views.md)
  for allocation-ready values.

- `neighbors_t`:

  Get or upsert
  [neighbors_t](https://ethzplus.github.io/evoland.plus/reference/neighbors_t.md)

- `reporting_t`:

  Get or upsert
  [reporting_t](https://ethzplus.github.io/evoland.plus/reference/reporting_t.md)

- `runs_t`:

  Get or upsert
  [runs_t](https://ethzplus.github.io/evoland.plus/reference/runs_t.md)

- `id_run`:

  Get or set active id_run, see
  [runs_t](https://ethzplus.github.io/evoland.plus/reference/runs_t.md)

- `run_lineage`:

  Get id_run, see
  [runs_t](https://ethzplus.github.io/evoland.plus/reference/runs_t.md)

## Methods

### Public methods

- [`evoland_db$adjusted_trans_pot_v()`](#method-evoland_db-adjusted_trans_pot_v)

- [`evoland_db$new()`](#method-evoland_db-initialize)

- [`evoland_db$get_read_expr()`](#method-evoland_db-get_read_expr)

- [`evoland_db$print()`](#method-evoland_db-print)

- [`evoland_db$set_report()`](#method-evoland_db-set_report)

- [`evoland_db$set_neighbors()`](#method-evoland_db-set_neighbors)

- [`evoland_db$generate_neighbor_predictors()`](#method-evoland_db-generate_neighbor_predictors)

- [`evoland_db$upsert_new_neighbors()`](#method-evoland_db-upsert_new_neighbors)

- [`evoland_db$add_predictor()`](#method-evoland_db-add_predictor)

- [`evoland_db$pred_wide_train_v()`](#method-evoland_db-pred_wide_train_v)

- [`evoland_db$pred_wide_predict_v()`](#method-evoland_db-pred_wide_predict_v)

- [`evoland_db$figure_of_merit_v()`](#method-evoland_db-figure_of_merit_v)

- [`evoland_db$lulc_crosstab_v()`](#method-evoland_db-lulc_crosstab_v)

- [`evoland_db$create_alloc_params_t()`](#method-evoland_db-create_alloc_params_t)

- [`evoland_db$lulc_data_as_rast()`](#method-evoland_db-lulc_data_as_rast)

- [`evoland_db$fit_full_models()`](#method-evoland_db-fit_full_models)

- [`evoland_db$fit_partial_models()`](#method-evoland_db-fit_partial_models)

- [`evoland_db$get_crossval_plots()`](#method-evoland_db-get_crossval_plots)

- [`evoland_db$set_full_trans_preds()`](#method-evoland_db-set_full_trans_preds)

- [`evoland_db$get_pred_filter_score()`](#method-evoland_db-get_pred_filter_score)

- [`evoland_db$predict_trans_pot()`](#method-evoland_db-predict_trans_pot)

- [`evoland_db$get_obs_trans_rates()`](#method-evoland_db-get_obs_trans_rates)

- [`evoland_db$clone()`](#method-evoland_db-clone)

Inherited methods

- [`ducklake_db$checkpoint()`](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.html#method-checkpoint)
- [`ducklake_db$commit()`](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.html#method-commit)
- [`ducklake_db$delete_from()`](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.html#method-delete_from)
- [`ducklake_db$execute()`](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.html#method-execute)
- [`ducklake_db$fetch()`](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.html#method-fetch)
- [`ducklake_db$get_query()`](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.html#method-get_query)
- [`ducklake_db$get_table_metadata()`](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.html#method-get_table_metadata)
- [`ducklake_db$list_tables()`](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.html#method-list_tables)
- [`ducklake_db$next_id()`](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.html#method-next_id)
- [`ducklake_db$row_count()`](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.html#method-row_count)
- [`ducklake_db$transaction()`](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.html#method-transaction)

------------------------------------------------------------------------

### `evoland_db$adjusted_trans_pot_v()`

#### Usage

    evoland_db$adjusted_trans_pot_v(
      id_period_post,
      id_run = self$id_run,
      id_run_trans_pot = id_run
    )

------------------------------------------------------------------------

### `evoland_db$new()`

Initialize a new evoland_db object

#### Usage

    evoland_db$new(
      path = NULL,
      id_run = 0L,
      read_only = FALSE,
      catalog = NULL,
      data_path = NULL,
      sqlite_journal_mode = "wal",
      expire_older_than = NULL,
      delete_older_than = NULL,
      ...
    )

#### Arguments

- `path`:

  Character. Path to the data folder; may be omitted when both `catalog`
  and `data_path` are given, see
  [ducklake_db](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.md).

- `id_run`:

  Atomic integer run ID, defaults to 0. Can be set to NULL

- `read_only`:

  Logical. Whether to update the reporting table upon initialization; if
  TRUE, the catalog is attached read-only.

- `catalog`:

  Character. DuckLake catalog connection, see
  [ducklake_db](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.md).

- `data_path`:

  Character. DuckLake data files location, see
  [ducklake_db](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.md).

- `sqlite_journal_mode`:

  Character. Journal mode for a SQLite catalog, see
  [ducklake_db](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.md).

- `expire_older_than, delete_older_than`:

  Retention for `$checkpoint()`, see
  [ducklake_db](https://ethzplus.github.io/evoland.plus/reference/ducklake_db.md).

- `...`:

  passed on to `set_report`

#### Returns

A new `evoland_db` object

------------------------------------------------------------------------

### `evoland_db$get_read_expr()`

Get SQL expression to read a table, respecting active run hierarchy. For
each data slice, returns data from the closest ancestor run that has it.

#### Usage

    evoland_db$get_read_expr(table_name, id_run = self$id_run)

#### Arguments

- `table_name`:

  Character string table name

- `id_run`:

  Integer, the run whose lineage is read; defaults to the active run.
  Another run can be read without changing the active one.

#### Returns

Character string SQL expression (a composable subquery)

------------------------------------------------------------------------

### `evoland_db$print()`

Print method for evoland_db

#### Usage

    evoland_db$print(...)

#### Arguments

- `...`:

  passed to super\$print

------------------------------------------------------------------------

### `evoland_db$set_report()`

Set reporting metadata, see
[`db_set_report()`](https://ethzplus.github.io/evoland.plus/reference/reporting_t.md)

#### Usage

    evoland_db$set_report(...)

#### Arguments

- `...`:

  each named argument is entered into the table with the argument name
  as its key

------------------------------------------------------------------------

### `evoland_db$set_neighbors()`

Set the neighbors table, see
[`set_neighbors()`](https://ethzplus.github.io/evoland.plus/reference/neighbors_t.md)

#### Usage

    evoland_db$set_neighbors(
      max_distance = 1000,
      distance_breaks = c(0, 100, 500, 1000),
      overwrite = FALSE,
      quiet = FALSE,
      chunksize = 1e+07
    )

#### Arguments

- `max_distance`:

  Maximum distance

- `distance_breaks`:

  Numeric vector of distance breaks

- `overwrite`:

  Logical, whether to overwrite existing neighbors_t table

- `quiet`:

  Logical, whether to suppress messages about progress

- `chunksize`:

  Integer, chunksize for writing; keeps down memory pressure

------------------------------------------------------------------------

### `evoland_db$generate_neighbor_predictors()`

Generate neighbor prediction table, i.e. "how many neighbors within
distance X are of type Y", see
[`generate_neighbor_predictors()`](https://ethzplus.github.io/evoland.plus/reference/neighbors_t.md)

#### Usage

    evoland_db$generate_neighbor_predictors()

------------------------------------------------------------------------

### `evoland_db$upsert_new_neighbors()`

Append new neighbors to predictor for a given period; depends on
generate_neighbor_predictors() having been run.

#### Usage

    evoland_db$upsert_new_neighbors(id_period, id_run = self$id_run)

#### Arguments

- `id_period`:

  Period to calculate predictors for

- `id_run`:

  Run whose map is counted and under which the predictors are written

------------------------------------------------------------------------

### `evoland_db$add_predictor()`

Add a predictor using the currently active id_run, see
[`add_predictor()`](https://ethzplus.github.io/evoland.plus/reference/pred_data_t.md)

#### Usage

    evoland_db$add_predictor(
      pred_data_raw,
      name,
      fill_value,
      pretty_name = name,
      description = NA_character_,
      orig_format = NA_character_,
      sources = list(),
      unit = NA_character_
    )

#### Arguments

- `pred_data_raw`:

  Data frame with columns `id_coord`, `id_period`, and `value`

- `name`:

  Unique short name

- `fill_value`:

  Value to substitute if a coordinate point in
  [coords_t](https://ethzplus.github.io/evoland.plus/reference/coords_t.md)
  does not have an explicit associated value

- `pretty_name`:

  char, Friendly name for use in reporting

- `description`:

  char, For use in reporting.

- `orig_format`:

  char, Format description of the underlying data (raster, vector…)

- `sources`:

  data.frame with url/md5sum columns, used for keeping track of
  underlying raw data, see
  [`download_and_verify()`](https://ethzplus.github.io/evoland.plus/reference/util_download.md)

- `unit`:

  char, SI unit for physical predictors, or descriptions like "bed
  nights/year" as a proxy for touristic activity

------------------------------------------------------------------------

### `evoland_db$pred_wide_train_v()`

Get transitions along with their predictor data in a wide data.table,
see
[`pred_wide_train_v()`](https://ethzplus.github.io/evoland.plus/reference/pred_data_t.md)

#### Usage

    evoland_db$pred_wide_train_v(id_trans, id_pred, id_period_anterior = NULL)

#### Arguments

- `id_trans`:

  Integer transition ID, see
  [trans_meta_t](https://ethzplus.github.io/evoland.plus/reference/trans_meta_t.md)

- `id_pred`:

  Optional integer vector of predictor IDs to include

- `id_period_anterior`:

  Optional integer vector of anterior periods to train on; `NULL`
  (default) uses every observed anterior period

------------------------------------------------------------------------

### `evoland_db$pred_wide_predict_v()`

Retrieve a wide view of the predictor data, see
[`pred_wide_predict_v()`](https://ethzplus.github.io/evoland.plus/reference/pred_data_t.md)

#### Usage

    evoland_db$pred_wide_predict_v(
      id_trans,
      id_period_anterior,
      id_run = self$id_run
    )

#### Arguments

- `id_trans`:

  Integer transition ID, see
  [trans_meta_t](https://ethzplus.github.io/evoland.plus/reference/trans_meta_t.md)

- `id_period_anterior`:

  Integer ID of period to retrieve data for

- `id_run`:

  Integer, the run whose predictors and anterior map are read

------------------------------------------------------------------------

### `evoland_db$figure_of_merit_v()`

Pontius' figure of merit of simulated against observed change, see
[`figure_of_merit_v()`](https://ethzplus.github.io/evoland.plus/reference/figure_of_merit_v.md)

#### Usage

    evoland_db$figure_of_merit_v(
      id_period_anterior,
      id_period_post,
      id_run_reference,
      id_run_simulated = self$id_run,
      by_transition = FALSE,
      exclude_id_lulc = integer(0)
    )

#### Arguments

- `id_period_anterior`:

  Integer, period of the initial map

- `id_period_post`:

  Integer, period of the observed and simulated maps

- `id_run_reference`:

  Integer, run holding the initial and observed maps

- `id_run_simulated`:

  Integer vector, runs holding the simulated maps; defaults to the
  active run

- `by_transition`:

  Logical; if `TRUE`, report per transition instead of overall.

- `exclude_id_lulc`:

  Integer vector of classes; cells whose class in the reference run, at
  `id_period_anterior` or `id_period_post`, is one of them are ignored.
  Default none.

------------------------------------------------------------------------

### `evoland_db$lulc_crosstab_v()`

Cross-tabulation of land use between two periods, see
[`lulc_crosstab_v()`](https://ethzplus.github.io/evoland.plus/reference/lulc_crosstab_v.md)

#### Usage

    evoland_db$lulc_crosstab_v(
      id_period_anterior,
      id_period_post,
      id_run_anterior = self$id_run,
      id_run_post = id_run_anterior
    )

#### Arguments

- `id_period_anterior, id_period_post`:

  Integer, the two periods to compare

- `id_run_anterior`:

  Integer, run to read `id_period_anterior` from; defaults to the active
  run

- `id_run_post`:

  Integer, run to read `id_period_post` from; defaults to
  `id_run_anterior`

------------------------------------------------------------------------

### `evoland_db$create_alloc_params_t()`

Estimates allocation parameters for all viable transitions under the
active `id_run`, see
[`create_alloc_params_t()`](https://ethzplus.github.io/evoland.plus/reference/alloc_params_t.md)

#### Usage

    evoland_db$create_alloc_params_t(id_period_post = NULL)

#### Arguments

- `id_period_post`:

  Optional integer vector; estimate only from the period pairs ending in
  these periods. `NULL` (default) uses every observed pair.

------------------------------------------------------------------------

### `evoland_db$lulc_data_as_rast()`

Retrieve LULC data as a SpatRaster object for a given period. See
[`lulc_data_as_rast()`](https://ethzplus.github.io/evoland.plus/reference/lulc_data_t.md)

#### Usage

    evoland_db$lulc_data_as_rast(id_period = NULL)

#### Arguments

- `id_period`:

  Optional integer vector of period IDs to include. If NULL (default),
  all periods are included.

------------------------------------------------------------------------

### `evoland_db$fit_full_models()`

Fit full models (trained on the complete dataset) for each viable
transition, see
[`fit_full_models()`](https://ethzplus.github.io/evoland.plus/reference/trans_models_t.md).
Two mutually exclusive modes: pass `learner` to train directly, or pass
`select_score` to pick the best partial model by score.

#### Usage

    evoland_db$fit_full_models(
      learner = NULL,
      select_score = NULL,
      select_maximize = TRUE,
      cluster = NULL,
      id_trans = NULL,
      id_period_anterior = NULL
    )

#### Arguments

- `learner`:

  An mlr3 `Learner` or `AutoTuner` for direct-learner mode (`NULL` when
  `select_score` is used).

- `select_score`:

  Measure ID string for score-select mode, e.g. `"classif.auc"` (`NULL`
  when `learner` is used).

- `select_maximize`:

  Logical; maximize (`TRUE`) or minimize (`FALSE`) the score.

- `cluster`:

  Optional cluster object for parallel processing

- `id_trans`:

  Optional integer vector; fit only these transitions.

- `id_period_anterior`:

  Optional integer vector; train only on transitions observed from these
  anterior periods. `NULL` (default) uses all observed periods.

------------------------------------------------------------------------

### `evoland_db$fit_partial_models()`

Fit partial models for each viable transition using stratified sampling.
Models are trained on a subsample and evaluated on held-out data, see
[`fit_partial_models()`](https://ethzplus.github.io/evoland.plus/reference/trans_models_t.md)
for details.

#### Usage

    evoland_db$fit_partial_models(
      learner,
      measures,
      sample_frac = 0.7,
      seed = NULL,
      cluster = NULL,
      id_trans = NULL,
      id_period_anterior = NULL
    )

#### Arguments

- `learner`:

  An mlr3 `Learner` or `AutoTuner` R6 object.

- `measures`:

  A vector of `Measure` names passed to
  [mlr3::msr](https://mlr3.mlr-org.com/reference/mlr_sugar.html) or a
  list of `Measure` objects for scoring the held-out split.

- `sample_frac`:

  Fraction in \\0, 1\\ for stratified sampling.

- `seed`:

  Random seed for reproducible sampling

- `cluster`:

  Optional cluster object for parallel processing

- `id_trans`:

  Optional integer vector; fit only these transitions.

- `id_period_anterior`:

  Optional integer vector; train and score only on transitions observed
  from these anterior periods. `NULL` (default) uses all observed
  periods.

------------------------------------------------------------------------

### `evoland_db$get_crossval_plots()`

Get cross-validation plots for stored predictions, see
[`get_crossval_plots()`](https://ethzplus.github.io/evoland.plus/reference/trans_models_t.md)

#### Usage

    evoland_db$get_crossval_plots(id_run = NULL, id_trans = NULL, ...)

#### Arguments

- `id_run`:

  Optional integer; filter by run ID.

- `id_trans`:

  Optional integer; filter by transition ID.

- `...`:

  Passed to
  [`mlr3viz::autoplot.PredictionClassif()`](https://mlr3viz.mlr-org.com/reference/autoplot.PredictionClassif.html)

------------------------------------------------------------------------

### `evoland_db$set_full_trans_preds()`

Set an initial full set of transition / predictor relations for the
active run, see
[`set_full_trans_preds()`](https://ethzplus.github.io/evoland.plus/reference/trans_preds_t.md)

#### Usage

    evoland_db$set_full_trans_preds(reset = FALSE)

#### Arguments

- `reset`:

  Logical, start the active run's selection over, removing its rows and
  their scores? Default FALSE.

------------------------------------------------------------------------

### `evoland_db$get_pred_filter_score()`

Add filter scores to predictors for each `id_run, id_trans`. See
[`get_pred_filter_score()`](https://ethzplus.github.io/evoland.plus/reference/trans_preds_t.md).

#### Usage

    evoland_db$get_pred_filter_score(
      filter = "correlation",
      id_trans = NULL,
      cluster = NULL,
      ...
    )

#### Arguments

- `filter`:

  Character passed to
  [mlr3filters::flt](https://mlr3filters.mlr-org.com/reference/flt.html)
  or
  [mlr3filters::Filter](https://mlr3filters.mlr-org.com/reference/Filter.html)
  object specifying the filter method to use for feature selection.

- `id_trans`:

  Optional integer vector; score only these transitions.

- `cluster`:

  Optional cluster object for parallel processing

- `...`:

  Additional arguments passed to `flt`.

------------------------------------------------------------------------

### `evoland_db$predict_trans_pot()`

Predict the raw transition potential for a given period and store in
`trans_pot_t`, see
[`predict_trans_pot()`](https://ethzplus.github.io/evoland.plus/reference/trans_pot_t.md).
Raw potentials are per-transition model probabilities (not yet
allocation-ready); use
[`adjusted_trans_pot_v()`](https://ethzplus.github.io/evoland.plus/reference/evoland_db_views.md)
to obtain column-scaled, row-closed values.

#### Usage

    evoland_db$predict_trans_pot(
      id_period_post,
      select_score,
      select_maximize,
      force = FALSE,
      id_run = self$id_run
    )

#### Arguments

- `id_period_post`:

  Integerish, posterior period of the transition potential interval

- `select_score`:

  Character string; mlr3 measure ID (e.g. `"classif.auc"`) used to
  select model for extrapolation

- `select_maximize`:

  Logical; maximize (`TRUE`) or minimize (`FALSE`) the score.

- `force`:

  Logical; predict even if `trans_pot_t` already holds potentials for
  `id_run` and the period.

- `id_run`:

  Integer; the run to predict for, and to write the potentials under.

------------------------------------------------------------------------

### `evoland_db$get_obs_trans_rates()`

Get the transition rates that were observed, see
[trans_rates_t](https://ethzplus.github.io/evoland.plus/reference/trans_rates_t.md)

#### Usage

    evoland_db$get_obs_trans_rates()

------------------------------------------------------------------------

### `evoland_db$clone()`

The objects of this class are cloneable with this method.

#### Usage

    evoland_db$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
