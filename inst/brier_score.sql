/*
Multi-category Brier score of ensembles of simulated maps, scored as probabilistic forecasts of
the observed map, per ensemble and anterior class; R aggregates over the classes.

Per cell, an ensemble of m runs forecasts each class with the share of its runs that put the
cell in that class. The Brier score of a cell is the squared distance between that forecast
and the observed class as a one-hot vector, summed over the classes; with the forecast
summing to 1, it is sum(p_k^2) - 2 p_observed + 1. The fair score (Ferro 2014) subtracts the
expected penalty of a finite ensemble, sum(p_k (1 - p_k)) / (m - 1) = (1 - sum(p_k^2)) / (m - 1).
The reference forecast is the uniform one: each cell of an anterior class forecasts the
observed shares of that class's transitions, i.e. the observed quantity of change placed at
random within the class.

Interpolated by `ducklake_db$get_query()`, requiring
Filters:
- {id_period_anterior}
- {id_period_post}
- {cell_filter}: a condition on a.id_lulc (anterior class) and o.id_lulc (observed class)
Data sources:
- {reference_read_expr}: lulc_data_t of the reference run, holding the anterior and observed map
- {simulated_read_expr}: rows (ensemble, id_run, id_coord, id_lulc) of each simulated run's map
  at id_period_post
*/
with
  cells as (
    select
      a.id_coord,
      a.id_lulc as anterior,
      o.id_lulc as observed
    from
      {reference_read_expr} a
      inner join {reference_read_expr} o using (id_coord)
    where
      a.id_period = {id_period_anterior}
      and o.id_period = {id_period_post}
      and {cell_filter}
  ),
  simulated as ({simulated_read_expr}),
  members as (
    select
      ensemble,
      count(distinct id_run) as m
    from
      simulated
    group by
      ensemble
  ),
  forecast as (
    -- the share of an ensemble's runs that put each scored cell in each class
    select
      s.ensemble,
      s.id_coord,
      s.id_lulc,
      count(*) / any_value(m.m) as p
    from
      simulated s
      inner join members m using (ensemble)
      semi join cells using (id_coord)
    group by
      s.ensemble,
      s.id_coord,
      s.id_lulc
  ),
  forecast_cell as (
    select
      f.ensemble,
      f.id_coord,
      c.anterior,
      sum(f.p * f.p) as sum_p2,
      sum(f.p) filter(
        where
          f.id_lulc = c.observed
      ) as p_observed
    from
      forecast f
      inner join cells c using (id_coord)
    group by
      f.ensemble,
      f.id_coord,
      c.anterior
  ),
  rates as (
    -- observed shares of each anterior class's transitions, persistence included
    select
      anterior,
      observed as id_lulc,
      count(*) / sum(count(*)) over (
        partition by
          anterior
      ) as p
    from
      cells
    group by
      anterior,
      observed
  ),
  reference as (
    select
      c.anterior,
      avg(r.sum_p2 - 2 * coalesce(r_obs.p, 0) + 1) as brier_reference
    from
      cells c
      inner join (
        select
          anterior,
          sum(p * p) as sum_p2
        from
          rates
        group by
          anterior
      ) r using (anterior)
      left join rates r_obs on r_obs.anterior = c.anterior
      and r_obs.id_lulc = c.observed
    group by
      c.anterior
  )
select
  f.ensemble,
  any_value(m.m)::integer as n_members,
  f.anterior as id_lulc_anterior,
  count(*)::integer as n_cells,
  avg(f.sum_p2 - 2 * coalesce(f.p_observed, 0) + 1) as brier,
  avg(
    f.sum_p2 - 2 * coalesce(f.p_observed, 0) + 1 - case
      when m.m > 1 then (1 - f.sum_p2) / (m.m - 1)
      else 0
    end
  ) as brier_fair,
  any_value(r.brier_reference) as brier_reference
from
  forecast_cell f
  inner join members m using (ensemble)
  inner join reference r using (anterior)
group by
  f.ensemble,
  f.anterior
order by
  f.ensemble,
  f.anterior
