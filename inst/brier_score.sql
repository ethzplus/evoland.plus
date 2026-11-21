/*
Multi-category Brier score of ensembles of simulated maps, scored as probabilistic forecasts of
the observed map, per ensemble and, optionally, per anterior class; and the coverage of the
scored cells by each member run, which the scores rely on.

Per cell, an ensemble of m runs forecasts each class with the share of its runs that put the
cell in that class. The Brier score of a cell is the squared distance between that forecast
and the observed class as a one-hot vector, summed over the classes; with the forecast
summing to 1, it is sum(p_k^2) - 2 p_observed + 1. The fair score (Ferro 2014) subtracts the
expected penalty of a finite ensemble, sum(p_k (1 - p_k)) / (m - 1) = (1 - sum(p_k^2)) / (m - 1).
The reference forecast is the uniform one: each cell of an anterior class forecasts the
observed shares of that class's transitions, i.e. the observed quantity of change placed at
random within the class.

Members are read from their own rows of lulc_data_t, not through their lineage: a member
without a map of its own for id_period_post would otherwise read the observed map. The scores
are only valid if every member covers every scored cell, which "coverage" reports.

Interpolated by `ducklake_db$get_query()`, requiring
Filters:
- {id_period_anterior}
- {id_period_post}
- {changeable_only}: if true, only cells whose anterior class has a viable transition
- {exclude_id_lulc*}: classes whose cells, as anterior or observed class, are left out
Ensembles, one element per member run:
- {ensemble*}: name of the ensemble
- {ensemble_position*}: position of the ensemble in the order given, for the output
- {id_run*}: the member run
Data sources:
- {reference_read_expr}: lulc_data_t of the reference run, holding the anterior and observed map
- {lulc_data_read_expr}: lulc_data_t as stored, every run's own rows
- {trans_meta_read_expr}: trans_meta_t
Output:
- {`result`}: "coverage", "overall" or "per_anterior"
*/
with
  members as (
    select
      unnest([{ensemble*}]) as ensemble,
      unnest([{ensemble_position*}]::integer[]) as ensemble_position,
      unnest([{id_run*}]::integer[]) as id_run
  ),
  ensembles as (
    select
      ensemble,
      any_value(ensemble_position) as ensemble_position,
      count(*) as m
    from
      members
    group by
      ensemble
  ),
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
      and (
        not {changeable_only}
        or a.id_lulc in (
          select
            id_lulc_anterior
          from
            {trans_meta_read_expr}
          where
            is_viable
        )
      )
      -- an empty exclusion arrives as [NULL], which contains no class
      and not list_contains([{exclude_id_lulc*}]::integer[], a.id_lulc)
      and not list_contains([{exclude_id_lulc*}]::integer[], o.id_lulc)
  ),
  simulated as (
    -- each member's own map at id_period_post
    select
      m.ensemble,
      s.id_coord,
      s.id_lulc
    from
      members m
      inner join {lulc_data_read_expr} s on s.id_run = m.id_run
      and s.id_period = {id_period_post}
  ),
  coverage as (
    -- n_rows = 0: no map of its own; n_cells_covered < n_cells: a map with gaps
    select
      r.id_run,
      count(s.id_coord) as n_rows,
      count(c.id_coord) as n_cells_covered,
      (
        select
          count(*)
        from
          cells
      ) as n_cells
    from
      (
        select distinct
          id_run
        from
          members
      ) r
      left join {lulc_data_read_expr} s on s.id_run = r.id_run
      and s.id_period = {id_period_post}
      left join cells c on c.id_coord = s.id_coord
    group by
      r.id_run
    order by
      r.id_run
  ),
  votes as (
    -- how many of an ensemble's runs put each scored cell in each class; the ensemble size is
    -- joined only after this aggregation, which, over all members' rows, is the costly step
    select
      s.ensemble,
      s.id_coord,
      s.id_lulc,
      count(*) as n_runs
    from
      simulated s
      semi join cells using (id_coord)
    group by
      s.ensemble,
      s.id_coord,
      s.id_lulc
  ),
  forecast as (
    -- the share of an ensemble's runs that put each scored cell in each class
    select
      v.ensemble,
      v.id_coord,
      v.id_lulc,
      v.n_runs / e.m as p
    from
      votes v
      inner join ensembles e using (ensemble)
  ),
  forecast_cell as (
    select
      f.ensemble,
      f.id_coord,
      c.anterior,
      c.observed,
      sum(f.p * f.p) as sum_p2,
      coalesce(
        sum(f.p) filter (
          where
            f.id_lulc = c.observed
        ),
        0
      ) as p_observed
    from
      forecast f
      inner join cells c using (id_coord)
    group by
      f.ensemble,
      f.id_coord,
      c.anterior,
      c.observed
  ),
  rates as (
    -- observed shares of each anterior class's transitions, persistence included
    select
      anterior,
      observed,
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
    -- the reference forecast's score on a cell, by its anterior and observed class
    select
      anterior,
      observed,
      sum(p * p) over (
        partition by
          anterior
      ) - 2 * p + 1 as brier_reference
    from
      rates
  ),
  cell_scores as (
    select
      f.ensemble,
      e.ensemble_position,
      e.m,
      f.anterior,
      f.sum_p2 - 2 * f.p_observed + 1 as brier,
      f.sum_p2 - 2 * f.p_observed + 1 - case
        when e.m > 1 then (1 - f.sum_p2) / (e.m - 1)
        else 0
      end as brier_fair,
      r.brier_reference
    from
      forecast_cell f
      inner join ensembles e using (ensemble)
      inner join reference r using (anterior, observed)
  ),
  scores as (
    -- per ensemble, and per ensemble and anterior class; a skill score is undefined where
    -- the reference is perfect, i.e. nothing changed
    select
      ensemble,
      any_value(ensemble_position) as ensemble_position,
      any_value(m)::integer as n_members,
      anterior as id_lulc_anterior,
      grouping(anterior) = 1 as is_overall,
      count(*)::integer as n_cells,
      avg(brier) as brier,
      avg(brier_fair) as brier_fair,
      avg(brier_reference) as brier_reference,
      1 - avg(brier) / nullif(avg(brier_reference), 0) as skill,
      1 - avg(brier_fair) / nullif(avg(brier_reference), 0) as skill_fair
    from
      cell_scores
    group by
      grouping sets ((ensemble), (ensemble, anterior))
  ),
  overall as (
    select
      ensemble,
      n_members,
      n_cells,
      brier,
      brier_fair,
      brier_reference,
      skill,
      skill_fair
    from
      scores
    where
      is_overall
    order by
      ensemble_position
  ),
  per_anterior as (
    select
      ensemble,
      n_members,
      id_lulc_anterior,
      n_cells,
      brier,
      brier_fair,
      brier_reference,
      skill,
      skill_fair
    from
      scores
    where
      not is_overall
    order by
      ensemble_position,
      id_lulc_anterior
  )
select
  *
from
  {`result`}
