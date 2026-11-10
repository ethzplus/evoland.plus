# Create Transition-Predictor Relationship Table

Creates a trans_preds_t table, which lists the transition / predictor
pairs a run considers. `is_selected` says which of them the transition
models use: feature selection sets it to `FALSE` for the predictors it
prunes. Further columns, such as the scores of
`get_pred_filter_score()`, are added to the table as they are written;
they are for reporting only. When writing a selection back, assign the
scored table with its flags rather than a subset of it: rows left out of
the assignment keep the flags they had.

## Usage

``` r
as_trans_preds_t(x)

# S3 method for class 'trans_preds_t'
print(x, nrow = 10, ...)

set_full_trans_preds(self, reset = FALSE)

get_pred_filter_score(
  self,
  filter = "correlation",
  id_trans = NULL,
  cluster = NULL,
  ...
)
```

## Arguments

- x:

  A list or data.frame coercible to a trans_preds_t object. If missing,
  an empty table will be created. If it has no `is_selected` column,
  every pair is selected.

- nrow:

  see
  [data.table::print.data.table](https://rdrr.io/pkg/data.table/man/print.data.table.html)

- ...:

  Additional arguments passed to `flt` if `filter` is a character string

- self:

  An
  [evoland_db](https://ethzplus.github.io/evoland.plus/reference/evoland_db.md)
  instance with populated trans_meta_t and pred_meta_t tables

- reset:

  Bool, start the active run's selection over? This removes the active
  run's rows, including their score columns' values, before writing the
  full set.

- filter:

  An
  [mlr3filters::Filter](https://mlr3filters.mlr-org.com/reference/Filter.html)
  object or a character string specifying the filter method, retrieved
  via
  [mlr3filters::flt](https://mlr3filters.mlr-org.com/reference/flt.html).
  Note that your filter must be compatible with the feature data types;
  compare your `pred_meta_t` table to <https://mlr3filters.mlr-org.com>
  for filter compatibility.

- id_trans:

  Optional integer vector of transition IDs; only these transitions are
  scored. `NULL` (default) scores every transition with selected
  predictors.

- cluster:

  An optional cluster object, see
  [run_parallel_evoland](https://ethzplus.github.io/evoland.plus/reference/evoland_db_util.md)

## Value

A data.table of class "trans_preds_t" with columns:

- `id_run`: Foreign key to runs_t

- `id_pred`: Foreign key to pred_meta_t

- `id_trans`: Foreign key to trans_meta_t

- `is_selected`: Logical, whether the models of the transition use the
  predictor

- any further columns, e.g. filter scores

## Methods (by generic)

- `print(trans_preds_t)`: Print a trans_preds_t object, passing params
  to data.table print

## Functions

- `set_full_trans_preds()`: Set an initial full set of transition /
  predictor relations for the active run: every predictor for every
  viable transition, all selected

- `get_pred_filter_score()`: Get a filter score for the selected
  transition-predictor relationships based on mlr3filters. Returns
  trans_preds_t with an additional column named after the filter\$id.
  The filter score can be used for feature selection: set `is_selected`
  according to the score and write the table back, e.g.
  `scores[, is_selected := importance > 100]; db$trans_preds_t <- scores`.
  Writing it back also stores the scores. Only the selected pairs are
  scored, so a second filter refines the selection of the first.

## Examples

``` r
if (FALSE) { # \dontrun{
# a typical selection step
scores <- db$get_pred_filter_score("importance")
scores[, is_selected := importance > 100]
db$trans_preds_t <- scores
} # }
```
