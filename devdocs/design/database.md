# Data model

## Overview

The analytic workflow is centrally managed around a relational database.
The database is a [DuckLake](https://ducklake.select/) catalog, accessed through an in-memory [DuckDB](https://duckdb.org/) instance, holding the tables described below.
By default a database is one folder: the catalog in `<path>/catalog.sqlite`, the data as Parquet files under `<path>/data/`.
The catalog can also live in PostgreSQL and the data in object storage (S3 and similar), so several machines can work on one database.

DuckLake was chosen over plain Parquet files because it gives atomic writes and snapshot isolation: several processes can write to one database concurrently, and a reader never sees a half-written table.
It supports no constraints, keys, or indexes so `ducklake_db$commit()` enforces key uniqueness itself.
Each table is declared once, as a `table_spec()` next to its `as_*_t()` constructor: its columns and their types, its keys, partitioning and map columns, and whether and how it is inherited along the run lineage.
The constructor, the reads through the lineage, the column lists in the help pages and the tables in this document are derived from that declaration; the registry of all of them is `EVOLAND_SCHEMA`.
The keys, partitioning, map columns and inheritance of each table are also stored in the catalog as the table's comment, one `key: "value1", "value2"` line each, so that other tools can read them.
An append or upsert that brings columns a table does not have adds them (`alter table ... add column`, which only touches the catalog; existing rows read `NULL`).
The validators of the large per-coordinate tables (`coords_t`, `lulc_data_t`, `pred_data_t`, `neighbors_t`, `trans_pot_t`) reject columns they do not declare, to keep those tables lean.
DuckLake also does not support ENUMs, so R factors are stored as strings and cast back by the `as_*_t()` constructors.

Because other tools may read and write the catalog, the schema should be considered to be part of the package's API; see [../style/r.md](../style/r.md#schema-names).

## Tables

Grouped by the [workflow](workflow.md) phase that writes them:

- Setup
  - `coords_t`: the spatial domain, one row per coordinate point
  - `periods_t`: the temporal domain, observed and extrapolated periods
  - `lulc_meta_t`: the LULC classes and the source classes they aggregate
- Data
  - `lulc_data_t`: observed (and later, allocated) LULC per coordinate and period
  - `pred_meta_t` + `pred_data_t`: predictor metadata and values
  - `neighbors_t`: neighbor relationships between coordinates, for neighborhood predictors
- Calibration
  - `trans_meta_t`: the transitions between classes and how often each was observed
  - `trans_viability_t`: whether each transition is viable for modelling, per run
  - `trans_preds_t`: which predictors each transition's model may use
  - `trans_models_t`: fitted transition potential models and their scores
- Allocation
  - `trans_rates_t`: transition rates per period, observed or demanded
  - `trans_pot_t`: predicted transition potentials per coordinate
  - `alloc_params_t`: patch size and shape parameters per transition
- Bookkeeping
  - `runs_t`: the run hierarchy every run-dependent table is keyed on
  - `reporting_t`: key-value report metadata

<!-- BEGIN GENERATED ER DIAGRAM: run `Rscript devdocs/er_diagram.R`, do not edit -->
```mermaid
erDiagram
    alloc_params_t {
        int id_run PK, FK
        int id_trans PK, FK
        float mean_patch_size
        float patch_size_variance
        float patch_elongation
        float patch_isometry
        float frac_expander
        float frac_patcher
    }
    coords_t {
        int id_coord PK
        float lon
        float lat
    }
    lulc_data_t {
        int id_run PK, FK
        int id_coord PK, FK
        int id_lulc FK
        int id_period PK, FK
    }
    lulc_meta_t {
        int id_lulc UK
        string name PK
        string pretty_name
        string description
        list src_classes
    }
    neighbors_t {
        int id_coord_origin PK, FK
        int id_coord_neighbor PK, FK
        float distance
    }
    periods_t {
        int id_period UK
        date start_date PK
        date end_date PK
        date mean_date
        int period_length_d
        bool is_extrapolated
    }
    pred_data_t {
        int id_run PK, FK
        int id_period PK, FK
        int id_pred PK, FK
        int id_coord PK, FK
        float value
    }
    pred_meta_t {
        int id_pred UK
        string name PK
        string pretty_name
        string description
        string orig_format
        list sources
        string unit
        list factor_levels
        factor data_type
        string fill_value
    }
    reporting_t {
        string key PK
        string value
    }
    runs_t {
        int id_run PK
        int parent_id_run FK
        string description
        int rng_seed
        int id_period_eval_start
    }
    trans_meta_t {
        int id_trans UK
        int id_lulc_anterior PK, FK
        int id_lulc_posterior PK, FK
        int cardinality
        float frequency_rel
        float frequency_abs
    }
    trans_models_t {
        int id_run PK, FK
        int id_trans PK, FK
        string learner_id PK
        map learner_params
        list learner_spec
        map crossval_score
        list crossval_predictions
        list learner_full
    }
    trans_pot_t {
        int id_run PK, FK
        int id_trans PK, FK
        int id_period_post PK, FK
        int id_coord PK, FK
        float value
    }
    trans_preds_t {
        int id_run PK, FK
        int id_pred PK, FK
        int id_trans PK, FK
        bool is_selected
    }
    trans_rates_t {
        int id_run PK, FK
        int id_period PK, FK
        int id_trans PK, FK
        int count
        float rate
    }
    trans_viability_t {
        int id_run PK, FK
        int id_trans PK, FK
        bool is_viable
    }

    alloc_params_t }o--|| runs_t : "id_run"
    alloc_params_t }o--|| trans_meta_t : "id_trans"
    lulc_data_t }o--|| runs_t : "id_run"
    lulc_data_t }o--|| coords_t : "id_coord"
    lulc_data_t }o--|| lulc_meta_t : "id_lulc"
    lulc_data_t }o--|| periods_t : "id_period"
    neighbors_t }o--|| coords_t : "id_coord_origin"
    neighbors_t }o--|| coords_t : "id_coord_neighbor"
    pred_data_t }o--|| runs_t : "id_run"
    pred_data_t }o--|| periods_t : "id_period"
    pred_data_t }o--|| pred_meta_t : "id_pred"
    pred_data_t }o--|| coords_t : "id_coord"
    runs_t }o--o| runs_t : "parent_id_run"
    trans_meta_t }o--|| lulc_meta_t : "id_lulc_anterior"
    trans_meta_t }o--|| lulc_meta_t : "id_lulc_posterior"
    trans_models_t }o--|| runs_t : "id_run"
    trans_models_t }o--|| trans_meta_t : "id_trans"
    trans_pot_t }o--|| runs_t : "id_run"
    trans_pot_t }o--|| trans_meta_t : "id_trans"
    trans_pot_t }o--|| periods_t : "id_period_post"
    trans_pot_t }o--|| coords_t : "id_coord"
    trans_preds_t }o--|| runs_t : "id_run"
    trans_preds_t }o--|| pred_meta_t : "id_pred"
    trans_preds_t }o--|| trans_meta_t : "id_trans"
    trans_rates_t }o--|| runs_t : "id_run"
    trans_rates_t }o--|| periods_t : "id_period"
    trans_rates_t }o--|| trans_meta_t : "id_trans"
    trans_viability_t }o--|| runs_t : "id_run"
    trans_viability_t }o--|| trans_meta_t : "id_trans"
```
<!-- END GENERATED ER DIAGRAM -->

Throughout this document, a _transition_ takes a coordinate point from its _anterior_ class (the class it had at the start of a period) to its _posterior_ class (the class it has at the end).
Likewise, `id_period_anterior` and `id_period_post` name the periods before and after a transition.

### `reporting_t`

A key/value (string) store for report metadata.

<!-- BEGIN GENERATED COLUMNS reporting_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname | type   | description          |
| ---------- | ------- | ------ | -------------------- |
| \*         | key     | string | Unique character key |
|            | value   | string | Value as character   |
<!-- END GENERATED COLUMNS reporting_t -->

On initialization, `evoland_db` sets `report_name`, `report_name_pretty`, `report_include_date` and `creator_username` if the table is empty, and updates `last_opened` and `last_opened_username` every time the database is opened read-write.

Few functions read this table yet; the plotting of model diagnostics should.
Which standard outputs (tables, graphs, maps) the package generates, and therefore whether this table earns its place, is open in ethzplus/evoland.plus#76.

### `coords_t`

<!-- BEGIN GENERATED COLUMNS coords_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname  | type  | description                                    |
| ---------- | -------- | ----- | ---------------------------------------------- |
| \*         | id_coord | int   | Unique ID for each coordinate point            |
|            | lon      | float | Longitude/x coordinate, in the unit of the CRS |
|            | lat      | float | Latitude/y coordinate, in the unit of the CRS  |
<!-- END GENERATED COLUMNS coords_t -->

`coords_t` specifies the coordinate points at which the model is set up and run.
It is either derived from a reference grid (e.g. existing LULC data) or generated from parameters; `create_coords_t_square()` builds a regular square grid from an extent and a resolution.
Triangular or hexagonal tilings, and a polygon column describing the surface each point stands for, are possible future extensions.

The CRS and grid geometry are stored as table metadata: `epsg`, the bounding box (`xmin`, `xmax`, `ymin`, `ymax`) and `resolution` for regular grids.

Every other spatial table refers to `id_coord` instead of carrying coordinates, so data are stored in long format, one row per coordinate point.
To turn a slice of `lulc_data_t` or `pred_data_t` back into a map, join it to `coords_t` and rasterize the points with `terra::rasterize()` against a template built from the grid metadata (see `lulc_data_as_rast()`).
This works for irregular coordinate points as well as for dense grids.

### `periods_t`

<!-- BEGIN GENERATED COLUMNS periods_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname         | type | description                                                                                |
| ---------- | --------------- | ---- | ------------------------------------------------------------------------------------------ |
| alternate  | id_period       | int  | Unique ID for each period                                                                  |
| \*         | start_date      | date | First day of the period                                                                    |
| \*         | end_date        | date | Last day of the period                                                                     |
|            | mean_date       | date | Midpoint of the period, derived from start and end date                                    |
|            | period_length_d | int  | Days between this and the preceding period's midpoints, derived from start and end date    |
|            | is_extrapolated | bool | Whether the period is used for extrapolation (`TRUE`) or matched to observations (`FALSE`) |
<!-- END GENERATED COLUMNS periods_t -->

This table identifies periods in the past and future.
`create_periods_t()` builds it from a start and end date for the observations, an end date for the extrapolation, and a period length.
Periods do not overlap: each period ends the day before the next one starts.
A date on the boundary between two periods therefore belongs to the later one.

Static predictors (variables that do not change over time) are stored under the special period 0, which starts and ends at the end of the observed range.
They can then be stored and joined like time-varying predictors.

### `runs_t`

<!-- BEGIN GENERATED COLUMNS runs_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname              | type   | description                                                                                                   |
| ---------- | -------------------- | ------ | ------------------------------------------------------------------------------------------------------------- |
| \*         | id_run               | int    | Unique ID for each run; 0 is the base run                                                                     |
|            | parent_id_run        | int    | The run this run inherits data from; `NA` for the base run                                                    |
|            | description          | string | Human-readable description of the run                                                                         |
|            | rng_seed             | int    | Seed for the run's allocation, see `allocate()`; `NA` for none                                                |
|            | id_period_eval_start | int    | First period the run allocates itself, i.e. where it branches off its parent, see `allocate()`; `NA` for none |
<!-- END GENERATED COLUMNS runs_t -->

Runs label data that either overrides base data (e.g. predictors for a scenario) or is produced by a specific experiment (e.g. one allocation in an ensemble).
Run 0 is the base run.
A child run inherits every slice of data it does not override: reading a run-dependent table under a run returns, for each slice, the data of the closest run in its lineage that has it.
For instance, a run that only overrides one predictor reads all other predictors from its parent, and so on up to the base run.
Which columns form a slice is declared per table, as the `inherit_by` of its specification: `id_period` for `lulc_data_t`, the pair of `id_period` and `id_pred` for `pred_data_t`, and so on.
A slice that holds `id_period` falls back to the static period 0 where the period itself is missing.

Users add columns to `runs_t` to document how a run was made (allocator, model selection, anything else); the package does not record these itself.
`rng_seed` and `id_period_eval_start`, `NA` where unset and not inherited, are read by `allocate()`.
Given several runs, `allocate()` goes period by period, each for every run that has reached its `id_period_eval_start`, so that a branch finds the maps of the run it branches off.
At its branching period, a run uses its parent's potentials unless it holds inputs of its own (predictors, the anterior map, potentials or `trans_viability_t` overrides), and shares its parent's allocation context if it also holds no rates or allocation parameters of its own.

### `lulc_meta_t` + `lulc_data_t`

`lulc_meta_t` describes each land use / land cover class.

<!-- BEGIN GENERATED COLUMNS lulc_meta_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname     | type   | description                                                      |
| ---------- | ----------- | ------ | ---------------------------------------------------------------- |
| alternate  | id_lulc     | int    | Unique ID for each land use class                                |
| \*         | name        | string | Concise name for use in code and queries, e.g. `forest_dense`    |
|            | pretty_name | string | Long name for plots and output, e.g. _Dense/Old Growth Forest_   |
|            | description | string | Long description / operationalisation                            |
|            | src_classes | list   | List of the source class IDs that are aggregated into this class |
<!-- END GENERATED COLUMNS lulc_meta_t -->

The `src_classes` column allows aggregating from source data class IDs to the harmonized LULC classes used in the model.

`lulc_data_t` records that at a given place and time, a certain class was observed or allocated.
Every coordinate point of `coords_t` must have a class in every observed period.

<!-- BEGIN GENERATED COLUMNS lulc_data_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname   | type | description                  |
| ---------- | --------- | ---- | ---------------------------- |
| \*         | id_run    | int  | Foreign key to `runs_t`      |
| \*         | id_coord  | int  | Foreign key to `coords_t`    |
|            | id_lulc   | int  | Foreign key to `lulc_meta_t` |
| \*         | id_period | int  | Foreign key to `periods_t`   |
<!-- END GENERATED COLUMNS lulc_data_t -->

The table is partitioned by `id_run`.

Some observational datasets, e.g. the Swiss [_Surface Statistics_](https://www.bfs.admin.ch/bfs/de/home/statistiken/raum-umwelt/nomenklaturen/arealstatistik.html), have exact dates for individual observations.
If they matter, add them as predictors, e.g. "number of years before or after the period's mean date".

### `pred_meta_t` + `pred_data_t`

`pred_meta_t` describes each predictor.

<!-- BEGIN GENERATED COLUMNS pred_meta_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname       | type   | description                                                                                                                                         |
| ---------- | ------------- | ------ | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| alternate  | id_pred       | int    | Unique ID for each predictor                                                                                                                        |
| \*         | name          | string | Concise name for use in code and queries, e.g. `distance_road`                                                                                      |
|            | pretty_name   | string | Name for plots and output, e.g. _Distance To Closest Road_                                                                                          |
|            | description   | string | Long description / operationalisation                                                                                                               |
|            | orig_format   | string | Original format, e.g. "100m raster, daily"                                                                                                          |
|            | sources       | list   | Where the data were fetched: list column of data.frames with columns `url` and `md5sum`                                                             |
|            | unit          | string | SI units for physical properties, or more complex descriptors like "bed nights/year" as a proxy for touristic activity; `NA` for categorical values |
|            | factor_levels | list   | List of character vectors, the ordered level labels when `data_type` is "factor" or "ordered"                                                       |
|            | data_type     | factor | Factor with levels "int", "float", "bool", "factor", "ordered": the type the stored values are cast back to                                         |
|            | fill_value    | string | Value of the `coords_t` coordinate points that are not explicitly stored in `pred_data_t`                                                           |
<!-- END GENERATED COLUMNS pred_meta_t -->

`pred_data_t` stores the predictor values.
All predictor types share one `value` column stored as float; `data_type` and `factor_levels` cast them back when read.

<!-- BEGIN GENERATED COLUMNS pred_data_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname   | type  | description                                                           |
| ---------- | --------- | ----- | --------------------------------------------------------------------- |
| \*         | id_run    | int   | Foreign key to `runs_t`                                               |
| \*         | id_period | int   | Foreign key to `periods_t`                                            |
| \*         | id_pred   | int   | Foreign key to `pred_meta_t`                                          |
| \*         | id_coord  | int   | Foreign key to `coords_t`                                             |
|            | value     | float | Predictor value, stored as double and cast back through `pred_meta_t` |
<!-- END GENERATED COLUMNS pred_data_t -->

Unlike `lulc_data_t`, `pred_data_t` need not hold a row for every coordinate point.
Values may be stored only where they are informative; reads join against `coords_t` and fill the remaining points with the predictor's `fill_value`.
A population predictor, for instance, only needs rows where the population is greater than 0, with a `fill_value` of 0.

The table is partitioned by `(id_run, id_period)`.

### `trans_meta_t`

Relates _anterior_ and _posterior_ classes to each other as _transitions_, with statistics about their frequency in the observed data.
`create_trans_meta_t()` derives it from the observed transitions (`db$trans_v`).
The table is the same for every run: `id_trans` identifies a pair of classes everywhere it is referenced.

<!-- BEGIN GENERATED COLUMNS trans_meta_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname           | type  | description                                                   |
| ---------- | ----------------- | ----- | ------------------------------------------------------------- |
| alternate  | id_trans          | int   | Unique ID for each transition                                 |
| \*         | id_lulc_anterior  | int   | Foreign key to `lulc_meta_t`, the class before the transition |
| \*         | id_lulc_posterior | int   | Foreign key to `lulc_meta_t`, the class after the transition  |
|            | cardinality       | int   | How many times this transition was observed                   |
|            | frequency_rel     | float | Frequency relative to all observed transitions                |
|            | frequency_abs     | float | Frequency relative to all coordinate points                   |
<!-- END GENERATED COLUMNS trans_meta_t -->

> [!NOTE]
> The previous generation of this model, [evoland-plus-legacy](https://github.com/ethzplus/evoland-plus-legacy), determined the viability of a transition per period.
> This meant that a single period for training our transition models had to be picked (usually the most recent), which made us lose out on the data provided by previous periods. We could still restrict our training periods in training now, but this only really makes sense if we suppose that different causal mechanisms drove land use change in different periods.

### `trans_viability_t`

<!-- BEGIN GENERATED COLUMNS trans_viability_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname   | type | description                                                                                  |
| ---------- | --------- | ---- | -------------------------------------------------------------------------------------------- |
| \*         | id_run    | int  | Foreign key to `runs_t`                                                                      |
| \*         | id_trans  | int  | Foreign key to `trans_meta_t`                                                                |
|            | is_viable | bool | Whether the transition is viable for modelling in this run (and the runs inheriting from it) |
<!-- END GENERATED COLUMNS trans_viability_t -->

Whether a transition is viable for modelling, per run.
A transition can be non-viable for conceptual or statistical reasons.
A model may contain an "other" class that simplifies what a class can turn into; since "other" says nothing about what comes next, transitions out of it cannot be modelled (`exclude_anterior`).
And a transition may be observed too rarely for statistical inference (`min_cardinality_abs`, `min_frequency_rel`).
Which transitions are modelled is also a modelling choice: a model may do better or worse when a rare transition is allowed or left out (ethzplus/evoland.plus#142).

The defaults are rows of the root run, derived from the thresholds and exclusions by `create_trans_viability_t()`:

```r
db$trans_meta_t <- create_trans_meta_t(db$trans_v())
db$trans_viability_t <- create_trans_viability_t(db$trans_meta_t, min_cardinality_abs = 5L)
```

A run that allows or disallows transitions differently adds rows of its own.
The table is read through the lineage per `id_trans`, like any run-dependent table: a transition takes the row of the closest run in the lineage that has one, and a transition without any row is viable.
`db$trans_meta_v(id_run)` returns `trans_meta_t` with the resolved `is_viable` of a run; predictor selection, model fitting, prediction, allocation parameters and allocation all work on that set.
A run that branches off its parent and holds rows of its own predicts its own potentials instead of taking its parent's, see `runs_t`.

Viability was a column of `trans_meta_t` until ethzplus/evoland.plus#142; databases written before then are not migrated.
The alternative of making `trans_meta_t` itself run-dependent was rejected: it duplicates the transition statistics per run and lets `id_trans` drift between runs although every other table refers to it.

### `trans_preds_t`

An m to n relation listing the predictors each transition's model may use.
`$set_full_trans_preds()` initializes it with the cartesian product of predictors and viable transitions, all selected.
Feature selection then prunes it: `$get_pred_filter_score()` scores the selected transition-predictor pairs with an [mlr3filters](https://mlr3filters.mlr-org.com) filter, the caller sets `is_selected` from the score, and the scored table is upserted.
Repeating this with different filters prunes the selection step by step, e.g. a cheap correlation filter first and an importance-based one on what remains.

Pruning flips `is_selected` rather than deleting rows, for two reasons.
A run inherits each pair from its closest ancestor that holds it, so a child run can only drop an inherited predictor by holding the pair with `is_selected = FALSE`; a missing row would let the ancestor's show through.
And the scores each filter writes, in a column named after the filter, stay next to the decision they informed.
These score columns are added to the table as they are first written and are not part of the schema readers rely on; readers use only the selected pairs.

<!-- BEGIN GENERATED COLUMNS trans_preds_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname     | type | description                                            |
| ---------- | ----------- | ---- | ------------------------------------------------------ |
| \*         | id_run      | int  | Foreign key to `runs_t`                                |
| \*         | id_pred     | int  | Foreign key to `pred_meta_t`                           |
| \*         | id_trans    | int  | Foreign key to `trans_meta_t`                          |
|            | is_selected | bool | Whether the models of the transition use the predictor |
<!-- END GENERATED COLUMNS trans_preds_t -->

### `trans_rates_t`

Transition rates per transition and period.
Rates are not probabilities: for observed periods they describe the change that happened, and for extrapolated periods they prescribe how much change allocation has to produce.

<!-- BEGIN GENERATED COLUMNS trans_rates_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname   | type  | description                                                                                                         |
| ---------- | --------- | ----- | ------------------------------------------------------------------------------------------------------------------- |
| \*         | id_run    | int   | Foreign key to `runs_t`                                                                                             |
| \*         | id_period | int   | Foreign key to `periods_t`                                                                                          |
| \*         | id_trans  | int   | Foreign key to `trans_meta_t`                                                                                       |
|            | count     | int   | Number of transitioning cells for (id_trans, id_period), for reporting; see Description                             |
|            | rate      | float | Transition rate: count of transitions in (id_trans, id_period) over count of cells of id_lulc_anterior in id_period |
<!-- END GENERATED COLUMNS trans_rates_t -->

Observed rates come from `$get_obs_trans_rates()`.
Rates for extrapolated periods can be extrapolated from the observed ones with a linear regression (`extrapolate_trans_rates()`), or derived from target class shares or areas with a linear program that yields a constrained optimal set of rates (`solve_trans_rates()`, see `trans_rates_lp`).

### `trans_models_t`

Transition potential models: for each viable `id_trans`, one or more fitted models with their scores.
The table is populated in one of two ways:

- Select by score: each call to `$fit_partial_models(learner, ...)` fits one candidate learner per transition on a training split and scores it on the held-out split; `$fit_full_models(select_score = ...)` refits the best-scoring candidate per transition on the full data.
- Direct: `$fit_full_models(learner = ...)` fits the given learner on the full data straight away, with no candidates to compare. An mlr3 `AutoTuner` works in either mode and does its own tuning internally.

<!-- BEGIN GENERATED COLUMNS trans_models_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname              | type   | description                                                                                                                                                                                                                                                         |
| ---------- | -------------------- | ------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| \*         | id_run               | int    | Foreign key to `runs_t`                                                                                                                                                                                                                                             |
| \*         | id_trans             | int    | Foreign key to `trans_meta_t`                                                                                                                                                                                                                                       |
| \*         | learner_id           | string | mlr3 twoclass `LearnerClassif`(https://mlr3.mlr-org.com/reference/LearnerClassif.html) key, e.g. `"classif.ranger"`                                                                                                                                                 |
|            | learner_params       | map    | MAP of atomic scalar learner hyperparameters for querying; complete hyperparameters are captured by `learner_spec`. A fit that failed is recorded as a sentinel row whose `learner_params` holds a single `error_message` entry, and whose `learner_full` is `NULL` |
|            | learner_spec         | list   | BLOB of serialized untrained mlr3 `Learner`; for AutoTuners, this is the optimal inner learner after tuning                                                                                                                                                         |
|            | crossval_score       | map    | MAP of cross-validation performance scores (from `prediction$score(measures)`)                                                                                                                                                                                      |
|            | crossval_predictions | list   | BLOB of serialized mlr3 `PredictionClassif` on the held-out test split                                                                                                                                                                                              |
|            | learner_full         | list   | BLOB of serialized trained mlr3 `Learner` fitted on the full dataset, used for extrapolation                                                                                                                                                                        |
<!-- END GENERATED COLUMNS trans_models_t -->

The table is partitioned by `id_run`.

The blob columns are the only part of the database that depends on R: they hold mlr3 objects serialized with `qs2`, which only R can read.
Every other table is plain data that any DuckLake client can read and write.

Each model is currently a two-class model of a single transition, so the potentials of a cell's transitions are not constrained to sum to 1.
Multiclass learners (ethzplus/evoland.plus#58) would fit one model per anterior class; they could be stored here under `id_trans = 0` or `NULL`.

### `alloc_params_t`

<!-- BEGIN GENERATED COLUMNS alloc_params_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname             | type  | description                                                                                                                                                                                       |
| ---------- | ------------------- | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| \*         | id_run              | int   | Foreign key to `runs_t`                                                                                                                                                                           |
| \*         | id_trans            | int   | Foreign key to `trans_meta_t`                                                                                                                                                                     |
|            | mean_patch_size     | float | Mean area of new patches (in cell units); a patch is a set of cells connected by their edges (rook adjacency)                                                                                     |
|            | patch_size_variance | float | Variance of patch area (in cell units)                                                                                                                                                            |
|            | patch_elongation    | float | Mean patch elongation, `1 - sqrt(lambda_2 / lambda_1)` from the eigenvalues of the patch's covariance, from 0 (isometric) to 1 (linear); the raw shape summary from `calculate_class_stats_cpp()` |
|            | patch_isometry      | float | Dinamica-specific isometry parameter derived from `patch_elongation` via `isometry_from_elongation()`                                                                                             |
|            | frac_expander       | float | Fraction of transition cells adjacent to existing patches, in [0, 1]                                                                                                                              |
|            | frac_patcher        | float | Fraction of transition cells forming new patches, in [0, 1]                                                                                                                                       |
<!-- END GENERATED COLUMNS alloc_params_t -->

Patch statistics measured on the observed transitions, used to parameterize allocation.
`create_alloc_params_t()` estimates one set per viable transition on the active run; perturbed variants for sensitivity analysis go on child runs.

The CLUMPY allocator reads `mean_patch_size`, `patch_size_variance` and `patch_elongation`.
`patch_isometry`, `frac_expander` and `frac_patcher` are only used by the Dinamica allocator; whether CLUMPY should split allocation into expansion and new patches too is open in ethzplus/evoland.plus#78.

### `neighbors_t`

<!-- BEGIN GENERATED COLUMNS neighbors_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname           | type             | description                                                  |
| ---------- | ----------------- | ---------------- | ------------------------------------------------------------ |
| \*         | id_coord_origin   | int              | Foreign key to `coords_t`, the origin coordinate point       |
| \*         | id_coord_neighbor | int              | Foreign key to `coords_t`, the neighbor coordinate point     |
|            | distance          | float            | Distance between origin and neighbor, in the unit of the CRS |
|            | distance_class    | factor, optional | Distance band, if `distance_breaks` were given               |
<!-- END GENERATED COLUMNS neighbors_t -->

Lists, for each coordinate point, the points within `max_distance`.
`$set_neighbors()` builds it with a spatial hash (`distance_neighbors_cpp()`); `$generate_neighbor_predictors()` then counts, per point, the neighbors of each class within each distance band.

On a dense raster this table is unnecessary: neighborhood counts follow from a focal (moving window) operation on the raster itself.
It is also expensive, with one row per pair of nearby points.
What it buys is generality: it works for irregular coordinate points, which have no raster to run a window over.

A `bearing` or `bearing_class` column recording the direction in which a neighbor lies would allow directional neighborhood predictors; see ethzplus/evoland.plus#77.

### `trans_pot_t`

The raw transition potential per coordinate point, transition and posterior period: the probability the transition's model predicts, in [0, 1].
Raw potentials are not calibrated to the demanded rates, and the potentials of a cell's transitions may sum to more than 1.
The `adjusted_trans_pot_v` view makes them allocation-ready (see [Views](#views)).

The table is written by `$predict_trans_pot()`.
It carries `id_run` so that each run can use its own potentials; with `allocate(..., potentials_ancestor = 1)`, sibling runs (e.g. a Monte Carlo ensemble) share the potentials predicted under their parent run.
Potentials can also be computed outside the package and written to the table directly.
It is partitioned by `id_run`.

<!-- BEGIN GENERATED COLUMNS trans_pot_t: run `Rscript devdocs/er_diagram.R`, do not edit -->
| uniqueness | colname        | type  | description                                                                |
| ---------- | -------------- | ----- | -------------------------------------------------------------------------- |
| \*         | id_run         | int   | Foreign key to `runs_t`                                                    |
| \*         | id_trans       | int   | Foreign key to `trans_meta_t`                                              |
| \*         | id_period_post | int   | Foreign key to `periods_t`, the posterior period                           |
| \*         | id_coord       | int   | Foreign key to `coords_t`                                                  |
|            | value          | float | Raw transition potential in [0, 1], as predicted by the transition's model |
<!-- END GENERATED COLUMNS trans_pot_t -->

## Views

Views are computed on demand and not stored.
They are exposed on `evoland_db` as active bindings, or as methods if they take parameters, and read through the active run's lineage.

| View                                                   | Returns                                                                                                                                                                         |
| ------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `trans_v`                                              | Observed transitions: one row per coordinate point and period, with the class in the preceding period as anterior. Input to `create_trans_meta_t()` and the observed rates. |
| `lulc_meta_long_v`                                     | `lulc_meta_t` unrolled to one row per source class, for reclassifying source data.                                                                                              |
| `pred_sources_v`                                       | Distinct predictor source URLs and their MD5 checksums.                                                                                                                         |
| `pred_wide_train_v(id_trans, id_pred)`                 | The design matrix for fitting one transition's model: `id_coord`, `id_period_anterior`, `did_transition`, and one `id_pred_<n>` column per predictor.                            |
| `pred_wide_predict_v(id_trans, id_period_anterior)`    | The same predictors in wide format for one period, as input to prediction.                                                                                                      |
| `adjusted_trans_pot_v(id_period_post)`                 | Allocation-ready potentials: each transition's raw potentials scaled so their mean matches the demanded rate, then divided by the cell's total wherever that exceeds 1.          |
| `figure_of_merit_v(...)`                               | Pontius' figure of merit of simulated against observed change, per run and optionally per transition.                                                                            |
| `coords_minimal`, `extent`                             | `coords_t` reduced to `id_coord`, `lon`, `lat` with grid metadata, and its extent as a `terra::SpatExtent`.                                                                     |

Views longer than a few lines live as SQL in `inst/*.sql`; see [../style/sql.md](../style/sql.md).
