# Allocate land use change

*Allocation* places the projected quantity of change on the map: for
each period, it turns the land use at the start of the period (the
*anterior* map) into the land use at its end (the *posterior* map),
given the transition potentials and the demanded transition rates.

`allocate()` runs a contiguous sequence of periods. For each period it

1.  predicts the transition potentials if `trans_pot_t` does not hold
    them yet, see
    [`predict_trans_pot()`](https://ethzplus.github.io/evoland.plus/reference/trans_pot_t.md);

2.  builds the allocation context: the anterior map, the
    allocation-ready potentials
    ([`adjusted_trans_pot_v()`](https://ethzplus.github.io/evoland.plus/reference/evoland_db_views.md)),
    the viable transitions, the quantity of change per transition and
    the allocation parameters;

3.  hands the context to the allocator, which returns the posterior map;

4.  commits the map to `lulc_data_t`, together with the neighbor
    predictors the next period is predicted from, in one transaction.

`allocate_one_period()` runs steps 1-3 for one period and returns the
map without committing it. Use it to write your own loop: edit
`trans_pot_t` between prediction and allocation (e.g. theory-driven
shifts in potential due to policy interventions), edit the allocated map
afterwards, or draw an ensemble of single-period realisations. Commit
the map with `db$commit(x, "lulc_data_t", method = "upsert")` and, if a
later period follows, run `db$upsert_new_neighbors(id_period_post)`.

The allocator is one of
[`alloc_clumpy_spec()`](https://ethzplus.github.io/evoland.plus/reference/alloc_clumpy_spec.md)
(stochastic, unbiased; the default),
[`alloc_greedy_spec()`](https://ethzplus.github.io/evoland.plus/reference/alloc_greedy_spec.md)
(deterministic rank-and-fill) or
[`alloc_dinamica_spec()`](https://ethzplus.github.io/evoland.plus/reference/alloc_dinamica_spec.md)
(Dinamica EGO's expander and patcher).

The quantity of change is the same for every allocator: each viable
transition changes `round(rate * n_anterior)` cells, where `rate` is the
transition's rate in `trans_rates_t` for the posterior period and
`n_anterior` the number of cells of its anterior class in the anterior
map. Transitions without a rate change no cells.

## Usage

``` r
allocate(
  db,
  id_periods,
  allocator = alloc_clumpy_spec(),
  select_score,
  select_maximize,
  id_run = db$id_run,
  potentials_ancestor = 0L,
  force_predict = FALSE,
  update_neighbors = TRUE
)

allocate_one_period(
  db,
  id_period_post,
  allocator = alloc_clumpy_spec(),
  select_score,
  select_maximize,
  id_run = db$id_run,
  potentials_ancestor = 0L,
  force_predict = FALSE
)
```

## Arguments

- db:

  An
  [evoland_db](https://ethzplus.github.io/evoland.plus/reference/evoland_db.md)
  instance.

- id_periods:

  Integer vector of contiguous posterior period IDs to simulate.

- allocator:

  An allocator specification, see Description.

- select_score:

  Character; mlr3 measure ID (e.g. `"classif.auc"`) used to select the
  model for prediction, see
  [`predict_trans_pot()`](https://ethzplus.github.io/evoland.plus/reference/trans_pot_t.md).

- select_maximize:

  Logical; maximize (`TRUE`) or minimize (`FALSE`) `select_score`.

- id_run:

  Integer; the run to allocate. It reads its anterior map, demand and
  allocation parameters through its lineage, and the allocated maps are
  written under it.

- potentials_ancestor:

  Integer; how many generations up the lineage of `id_run` to take the
  transition potentials from: `0` (default) the run itself, `1` its
  parent, `2` its grandparent. They are predicted under that run if
  missing. Use `1` to let sibling runs (e.g. a Monte Carlo ensemble of
  allocations) share their parent's potentials. `2` only makes sense if
  the parent's siblings share the same transition potentials too. Above
  `0` for a single period only: from the second period on, the anterior
  map is the run's own allocation, which its ancestors cannot see.

- force_predict:

  Logical; if `TRUE`, predict the transition potentials even if
  `trans_pot_t` already holds them for that run and the period.

- update_neighbors:

  Logical; whether to recompute the neighbor predictors after the last
  period (default `TRUE`). Intermediate periods are always updated,
  since the next period's prediction reads them. Set `FALSE` when
  nothing is allocated onwards from the last period, e.g. for
  single-period ensembles.

- id_period_post:

  Integer posterior period ID.

## Value

`allocate()`: called for its side effects on `lulc_data_t` and
`pred_data_t`; returns `NULL` invisibly. `allocate_one_period()`: an
[lulc_data_t](https://ethzplus.github.io/evoland.plus/reference/lulc_data_t.md)
with the posterior map.
