/*
Which run allocates which period, on whose potentials and allocation context: one row per run
and period. A run allocates from where it branches off (id_period_eval_start) on, and from the
first of the requested periods otherwise. In the period it branches off, it takes its parent's
potentials unless it holds inputs of its own (the anterior map, predictors before the period,
potentials or overrides of transition viability), and its parent's whole allocation context
unless it also holds its own rates or allocation parameters.

Interpolated by `ducklake_db$get_query()`, requiring
Filters:
- {id_run*}: the runs to allocate
- {id_period_first}, {id_period_last}: the first and last of the requested periods
Data sources, each read as stored, i.e. every run's own rows:
- {runs_read_expr}: runs_t
- {lulc_read_expr}: lulc_data_t
- {pred_read_expr}: pred_data_t
- {pot_read_expr}: trans_pot_t
- {rates_read_expr}: trans_rates_t
- {params_read_expr}: alloc_params_t
- {viability_read_expr}: trans_viability_t
*/
with
  runs as (
    select
      id_run,
      parent_id_run,
      rng_seed,
      id_period_eval_start,
      greatest(
        {id_period_first},
        coalesce(id_period_eval_start, {id_period_first})
      ) as id_period_start
    from
      {runs_read_expr}
    where
      id_run in ({id_run*})
  ),
  schedule as (
    select
      r.*,
      unnest(
        generate_series(r.id_period_start, {id_period_last})
      ) as id_period
    from
      runs r
  ),
  -- where a run branches off, which of its inputs and context it holds itself
  branching as (
    select
      s.id_run,
      s.id_period,
      s.parent_id_run,
      (
        exists (
          select
            1
          from
            {lulc_read_expr} l
          where
            l.id_run = s.id_run
            and l.id_period = s.id_period - 1
        )
        or exists (
          select
            1
          from
            {pred_read_expr} p
          where
            p.id_run = s.id_run
            and p.id_period < s.id_period
        )
        or exists (
          select
            1
          from
            {pot_read_expr} t
          where
            t.id_run = s.id_run
            and t.id_period_post = s.id_period
        )
        -- a run allowing other transitions than its parent needs other potentials
        or exists (
          select
            1
          from
            {viability_read_expr} v
          where
            v.id_run = s.id_run
        )
      ) as has_own_inputs,
      (
        exists (
          select
            1
          from
            {rates_read_expr} r
          where
            r.id_run = s.id_run
            and r.id_period = s.id_period
        )
        or exists (
          select
            1
          from
            {params_read_expr} a
          where
            a.id_run = s.id_run
        )
      ) as has_own_rates_or_params
    from
      schedule s
    where
      s.id_period = s.id_period_eval_start
      and s.parent_id_run is not null
  )
select
  s.id_period,
  s.id_run,
  s.rng_seed,
  case
    when not b.has_own_inputs then b.parent_id_run
    else s.id_run
  end as potentials_run,
  case
    when not (
      b.has_own_inputs
      or b.has_own_rates_or_params
    ) then b.parent_id_run
    else s.id_run
  end as context_run
from
  schedule s
  left join branching b on b.id_run = s.id_run
  and b.id_period = s.id_period
order by
  s.id_period,
  s.id_run
