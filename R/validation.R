#' Validation of simulated land use change
#'
#' @description
#' Compare a simulated map against the observed one, both starting from the same initial map:
#'
#' - [lulc_crosstab_v()] (`evoland_db$lulc_crosstab_v()`): cross-tabulation of land use between
#'   two periods, of one run or of two;
#' - [figure_of_merit_v()] (`evoland_db$figure_of_merit_v()`): Pontius' figure of merit and its
#'   components, cell by cell, with the value expected from random allocation for comparison;
#' - [brier_score_v()] (`evoland_db$brier_score_v()`): ensembles of simulated maps scored as
#'   probabilistic forecasts of the observed map, with the Brier score, its fair variant for
#'   finite ensembles, and the skill against random allocation of the observed quantities;
#' - [calc_transition_similarity()]: fuzzy similarity of differences for one transition, which
#'   tolerates change placed near, rather than exactly on, the observed change.
#'
#' @name validation
NULL

#' Where land use changed between two maps
#'
#' @param initial_map,final_map SpatRasters on the same grid
#' @param from_class,to_class Integer, optional: only count this transition. Give both or
#'   neither.
#'
#' @return Logical SpatRaster: `TRUE` where the land use changed (from `from_class` to
#'   `to_class`, if given), `FALSE` where it did not, `NA` where either map is `NA`.
#' @keywords internal
create_change_map <- function(initial_map, final_map, from_class = NULL, to_class = NULL) {
  stopifnot(
    "Both from_class and to_class must be specified, or both must be NULL" = {
      is.null(from_class) == is.null(to_class)
    }
  )
  if (is.null(from_class)) {
    return(initial_map != final_map)
  }
  initial_map == from_class & final_map == to_class
}

#' Cross-tabulation of land use between two periods
#'
#' @description
#' Counts the cells by their land use class in one period and in another: the transition
#' matrix, in long form. Each period is read from its own run, through that run's lineage, so
#' the same call gives observed change (both periods from one run), simulated change (a
#' simulated run from its initial period), or agreement between an observed and a simulated
#' map of the same period (two runs, one period).
#'
#' @param self An [evoland_db] instance
#' @param id_period_anterior,id_period_post Integer, the two periods to compare
#' @param id_run_anterior Integer, run to read `id_period_anterior` from; defaults to the active
#'   run
#' @param id_run_post Integer, run to read `id_period_post` from; defaults to `id_run_anterior`
#'
#' @return A [data.table::data.table()] with `id_lulc_anterior`, `id_lulc_posterior` and
#'   `n_cells`, one row per pair of classes that occurs. Cells missing from either map are
#'   ignored. Use [data.table::dcast()] for the matrix form.
#'
#' @keywords internal
lulc_crosstab_v <- function(
  self,
  id_period_anterior,
  id_period_post,
  id_run_anterior = self$id_run,
  id_run_post = id_run_anterior
) {
  stopifnot(
    "id_period_anterior must be a single integer" = is_scalar_integerish(id_period_anterior),
    "id_period_post must be a single integer" = is_scalar_integerish(id_period_post),
    "id_run_anterior must be a single integer; set an active run or pass it" = is_scalar_integerish(
      id_run_anterior
    ),
    "id_run_post must be a single integer" = is_scalar_integerish(id_run_post)
  )

  self$get_query(
    r"{
    select
      a.id_lulc as id_lulc_anterior,
      p.id_lulc as id_lulc_posterior,
      count(*)::integer as n_cells
    from
      {anterior_read_expr} a
      inner join {post_read_expr} p using (id_coord)
    where
      a.id_period = {id_period_anterior}
      and p.id_period = {id_period_post}
    group by
      a.id_lulc,
      p.id_lulc
    order by
      a.id_lulc,
      p.id_lulc
    }",
    anterior_read_expr = self$get_read_expr("lulc_data_t", id_run = id_run_anterior),
    post_read_expr = self$get_read_expr("lulc_data_t", id_run = id_run_post)
  )
}

#' Figure of merit of simulated land use change
#'
#' @description
#' Compares simulated against observed change with Pontius' figure of merit (Pontius et al.
#' 2008). The initial and observed maps are read from a reference run, the simulated maps from
#' one or more other runs, each through its own run lineage. Every cell that changed in either
#' the observed or a simulated map is one of:
#'
#' - **hit**: observed change, simulated as change to the observed class;
#' - **wrong hit**: observed change, simulated as change to another class;
#' - **miss**: observed change, simulated as persistence;
#' - **false alarm**: observed persistence, simulated as change.
#'
#' The figure of merit is `hits / (misses + hits + wrong hits + false alarms)`, the share of the
#' union of observed and simulated change that the simulation got right. Cells that persisted in
#' both maps do not enter it, so, unlike overall agreement, it is not inflated by a mostly
#' unchanging landscape. Its value scales with how much change there is to get right, so compare
#' it against `figure_of_merit_null`, the value expected when the simulated quantity of each
#' transition is placed at random among the cells of its initial class.
#'
#' Cells whose change is not the model's to get right, such as a class whose transitions are
#' imposed deterministically, can be left out with `exclude_id_lulc`. Which cells are left out
#' is decided on the reference run alone: a cell is ignored, in every simulated run, if its class
#' in the reference run's initial or observed map is one of them. Its simulated class does not
#' matter.
#'
#' For an ensemble of runs, the expected counts are the means of the per-run counts; recompute
#' the ratios from those means rather than averaging the per-run ratios.
#'
#' @param self An [evoland_db] instance
#' @param id_period_anterior Integer, period of the initial map
#' @param id_period_post Integer, period of the observed and simulated maps
#' @param id_run_reference Integer, run holding the initial and observed maps
#' @param id_run_simulated Integer vector, runs holding the simulated maps; defaults to the
#'   active run
#' @param by_transition Logical; if `TRUE`, report per transition instead of overall.
#' @param exclude_id_lulc Integer vector of classes; cells whose class in the reference run, at
#'   `id_period_anterior` or `id_period_post`, is one of them are ignored. Default none.
#'
#' @return A [data.table::data.table()] with one row per simulated `id_run` (and transition).
#'   Overall (`by_transition = FALSE`): `hits`, `wrong_hits`, `misses`, `false_alarms`,
#'   `figure_of_merit`, `producers_accuracy` (`hits` over observed change), `users_accuracy`
#'   (`hits` over simulated change) and `figure_of_merit_null`. Per transition:
#'   `id_lulc_anterior`, `id_lulc_posterior`, the `observed` and `simulated` cell counts,
#'   `hits`, `figure_of_merit` (`hits` over the union of observed and simulated change) and
#'   `figure_of_merit_null`. Cells missing from any of the three maps are ignored.
#'
#' @references
#' Pontius, R. G., Boersma, W., Castella, J.-C., et al. (2008). Comparing the input, output, and
#' validation maps for several models of land change. The Annals of Regional Science, 42(1),
#' 11-37. https://doi.org/10.1007/s00168-007-0138-2
#'
#' @keywords internal
figure_of_merit_v <- function(
  self,
  id_period_anterior,
  id_period_post,
  id_run_reference,
  id_run_simulated = self$id_run,
  by_transition = FALSE,
  exclude_id_lulc = integer(0)
) {
  stopifnot(
    "id_period_anterior must be a single integer" = is_scalar_integerish(id_period_anterior),
    "id_period_post must be a single integer" = is_scalar_integerish(id_period_post),
    "id_run_reference must be a single integer" = is_scalar_integerish(id_run_reference),
    "id_run_simulated must be a non-empty vector of integers" = {
      length(id_run_simulated) >= 1L && is_integerish(id_run_simulated)
    },
    "exclude_id_lulc must be a vector of integers" = is_integerish(exclude_id_lulc)
  )

  # each run is read through its own lineage
  lulc_read_expr <- function(id_run) self$get_read_expr("lulc_data_t", id_run = id_run)
  reference_read_expr <- lulc_read_expr(id_run_reference)
  simulated_read_expr <- glue::glue_sql_collapse(
    lapply(id_run_simulated, function(id_run) {
      glue::glue_sql(
        "select {id_run}::integer as id_run, id_coord, id_lulc
        from {lulc_read_expr(id_run)} where id_period = {id_period_post}",
        .con = self$connection
      )
    }),
    sep = " union all "
  )
  exclude_filter <- if (length(exclude_id_lulc) > 0L) {
    glue::glue_sql(
      "and a.id_lulc not in ({exclude_id_lulc*}) and o.id_lulc not in ({exclude_id_lulc*})",
      .con = self$connection
    )
  } else {
    DBI::SQL("")
  }

  self$get_query(
    read_sql("figure_of_merit.sql"),
    result = if (by_transition) "per_transition" else "overall"
  )
}

#' Brier score of ensembles of simulated maps
#'
#' @description
#' Scores ensembles of simulated maps as probabilistic forecasts of the observed map. The
#' figure of merit ([figure_of_merit_v()]) scores one map at a time and rewards allocators that
#' concentrate change on the most probable cells; an ensemble of realisations of a stochastic
#' allocator is better judged by how well the share of its realisations that put a cell in a
#' class forecasts where change happens.
#'
#' For each ensemble of `m` runs and each scored cell, the forecast probability of class `k` is
#' the share of the runs that put the cell in `k`. The multi-category **Brier score** (Brier
#' 1950) is the squared distance between that forecast and the observed class as a one-hot
#' vector, summed over the classes and averaged over the cells. It ranges from 0 to 2; a single
#' map scores 0 on a correct cell and 2 on any miss, false alarm or wrong hit.
#'
#' With few members, the forecast probabilities are coarse, which costs an ensemble for its
#' size rather than its skill. The **fair Brier score** (Ferro 2014) subtracts that expected
#' penalty, `sum_k p_k (1 - p_k) / (m - 1)` per cell, and estimates the score of the ensemble
#' with infinitely many members. It equals the Brier score for `m = 1`, and assumes
#' exchangeable members: realisations of one allocator from one set of potentials, not a mix of
#' learners or parameters.
#'
#' The **skill scores** compare against a reference forecast, `1 - brier / brier_reference`:
#' 0 for no better than the reference, 1 for perfect, negative for worse. The reference is the
#' uniform forecast (Murphy 1988): each cell of an anterior class forecasts the observed shares
#' of that class's transitions, the cell-wise expectation of placing the observed quantity of
#' each transition at random within its class, as `figure_of_merit_null` does.
#'
#' @param self An [evoland_db] instance
#' @param id_period_anterior Integer, period of the initial map
#' @param id_period_post Integer, period of the observed and simulated maps
#' @param id_run_reference Integer, run holding the initial and observed maps
#' @param ensembles Named list of integer vectors, the runs of each ensemble, read through
#'   their own lineages; or an integer vector for a single ensemble. Unnamed ensembles are
#'   numbered.
#' @param by_anterior Logical; if `TRUE`, report per anterior class instead of overall.
#' @param changeable_only Logical; if `TRUE` (default), score only the cells whose anterior
#'   class has a viable transition in [trans_meta_t]. The cells of other classes cannot
#'   change, and every forecast, the reference included, gets them right.
#' @param exclude_id_lulc Integer vector of classes; cells whose class in the reference run,
#'   at `id_period_anterior` or `id_period_post`, is one of them are ignored. Default none.
#'
#' @return A [data.table::data.table()] with one row per ensemble (and anterior class):
#'   `ensemble`, `n_members`, (`id_lulc_anterior`), `n_cells`, `brier`, `brier_fair`,
#'   `brier_reference`, `skill` and `skill_fair`. Cells that a run's map does not cover are
#'   scored on the runs that do.
#'
#' @references
#' Brier, G. W. (1950). Verification of forecasts expressed in terms of probability. Monthly
#' Weather Review, 78(1), 1-3. https://doi.org/10.1175/1520-0493(1950)078<0001:VOFEIT>2.0.CO;2
#'
#' Ferro, C. A. T. (2014). Fair scores for ensemble forecasts. Quarterly Journal of the Royal
#' Meteorological Society, 140(683), 1917-1923. https://doi.org/10.1002/qj.2270
#'
#' Murphy, A. H. (1988). Skill scores based on the mean square error and their relationships to
#' the correlation coefficient. Monthly Weather Review, 116(12), 2417-2424.
#' https://doi.org/10.1175/1520-0493(1988)116<2417:SSBOTM>2.0.CO;2
#'
#' @keywords internal
brier_score_v <- function(
  self,
  id_period_anterior,
  id_period_post,
  id_run_reference,
  ensembles,
  by_anterior = FALSE,
  changeable_only = TRUE,
  exclude_id_lulc = integer(0)
) {
  if (!is.list(ensembles)) {
    ensembles <- list(ensembles)
  }
  if (is.null(names(ensembles))) {
    names(ensembles) <- as.character(seq_along(ensembles))
  }
  stopifnot(
    "id_period_anterior must be a single integer" = is_scalar_integerish(id_period_anterior),
    "id_period_post must be a single integer" = is_scalar_integerish(id_period_post),
    "id_run_reference must be a single integer" = is_scalar_integerish(id_run_reference),
    "ensembles must be non-empty vectors of integers" = {
      length(ensembles) >= 1L &&
        all(vapply(ensembles, function(x) length(x) >= 1L && is_integerish(x), logical(1)))
    },
    "ensembles must have unique names" = !anyDuplicated(names(ensembles)),
    "by_anterior must be TRUE or FALSE" = is_scalar_logical(by_anterior),
    "changeable_only must be TRUE or FALSE" = is_scalar_logical(changeable_only),
    "exclude_id_lulc must be a vector of integers" = is_integerish(exclude_id_lulc)
  )

  sql <- function(...) glue::glue_sql(..., .con = self$connection, .envir = parent.frame())
  lulc_read_expr <- function(id_run) self$get_read_expr("lulc_data_t", id_run = id_run)
  simulated_read_expr <- glue::glue_sql_collapse(
    unlist(lapply(names(ensembles), function(ensemble) {
      lapply(ensembles[[ensemble]], function(id_run) {
        sql(
          "select {ensemble}::varchar as ensemble, {id_run}::integer as id_run, id_coord, id_lulc
          from {lulc_read_expr(id_run)} where id_period = {id_period_post}"
        )
      })
    })),
    sep = " union all "
  )

  conditions <- list(DBI::SQL("true"))
  if (changeable_only) {
    changeable <- unique(viable_trans_meta(self)[["id_lulc_anterior"]])
    conditions <- c(conditions, sql("a.id_lulc in ({changeable*})"))
  }
  if (length(exclude_id_lulc) > 0L) {
    conditions <- c(
      conditions,
      sql("a.id_lulc not in ({exclude_id_lulc*}) and o.id_lulc not in ({exclude_id_lulc*})")
    )
  }

  per_class <- self$get_query(
    read_sql("brier_score.sql"),
    reference_read_expr = lulc_read_expr(id_run_reference),
    simulated_read_expr = simulated_read_expr,
    cell_filter = glue::glue_sql_collapse(conditions, sep = " and ")
  )

  scores <- if (by_anterior) {
    per_class
  } else {
    # the overall mean is the mean of the per-class means, weighted by their cells
    per_class[,
      .(
        n_cells = sum(n_cells),
        brier = stats::weighted.mean(brier, n_cells),
        brier_fair = stats::weighted.mean(brier_fair, n_cells),
        brier_reference = stats::weighted.mean(brier_reference, n_cells)
      ),
      by = .(ensemble, n_members)
    ]
  }
  scores[, `:=`(
    skill = 1 - brier / brier_reference,
    skill_fair = 1 - brier_fair / brier_reference
  )]
  # in the order the ensembles were given
  scores[order(match(ensemble, names(ensembles)))]
}

#' Compute fuzzy similarity of differences for transition validation
#'
#' @description
#' Implements the "similarity of differences" approach from Dinamica EGO (Soares-Filho et al.
#' 2009), after the fuzzy set comparison of Hagen (2003). Compares the spatial pattern of a
#' specific transition in observed vs simulated maps, allowing for spatial tolerance.
#'
#' @details
#' Both maps are reduced to where the transition happened. For each cell that changed in one
#' map, its membership in the other map's change is the distance decay to the nearest cell that
#' changed there, within the window: 1 for a change in the same cell, `exp(-d / decay_divisor)`
#' for one `d` cells away (or 1 anywhere in the window without decay), 0 if there is none. The
#' directional similarity is the mean of these memberships over the changed cells of the first
#' map only; cells that changed in neither map do not enter it. The overall similarity is the
#' minimum of the two directions, so a simulation cannot score well by placing much more (or
#' much less) change than was observed.
#'
#' @param initial_map SpatRaster, initial LULC state
#' @param observed_map SpatRaster, observed final LULC state
#' @param simulated_map SpatRaster, simulated final LULC state
#' @param from_class Integer, initial class of transition
#' @param to_class Integer, final class of transition
#' @param window_size Integer, size of moving window (must be odd). Default 11.
#' @param use_exp_decay Logical, use exponential decay? Default TRUE.
#' @param decay_divisor Numeric, attenuation factor. Default 2.
#'
#' @return List with:
#'   - observed_change: Logical SpatRaster, `TRUE` where the transition was observed
#'   - simulated_change: Logical SpatRaster, `TRUE` where the transition was simulated
#'   - similarity: Minimum of the two directional similarities; `NA` if the transition
#'     happened in neither map
#'   - similarity_observed_to_simulated: Mean membership of observed changes in the simulated change
#'   - similarity_simulated_to_observed: Mean membership of simulated changes in the observed change
#'   - similarity_map: SpatRaster of the observed-to-simulated membership at each observed
#'     change (`NA` elsewhere)
#'   - n_observed: Number of cells with observed transition
#'   - n_simulated: Number of cells with simulated transition
#'
#' @references
#' Hagen, A. (2003). Fuzzy set approach to assessing similarity of categorical maps.
#' International Journal of Geographical Information Science, 17(3), 235-249.
#' https://doi.org/10.1080/13658810210157822
#'
#' Soares-Filho, B. S., Rodrigues, H. O., & Costa, W. L. (2009). Modeling Environmental
#' Dynamics with Dinamica EGO. Centro de Sensoriamento Remoto, UFMG.
#'
#' @export
calc_transition_similarity <- function(
  initial_map,
  observed_map,
  simulated_map,
  from_class,
  to_class,
  window_size = 11L, # TODO insert warning if domain is too small for window
  use_exp_decay = TRUE,
  decay_divisor = 2.0
) {
  stopifnot(
    "window_size must be odd" = (window_size %% 2L) == 1L,
    "window_size must be positive" = window_size > 0L,
    "decay_divisor must be positive" = decay_divisor > 0
  )

  obs_change <- create_change_map(initial_map, observed_map, from_class, to_class)
  sim_change <- create_change_map(initial_map, simulated_map, from_class, to_class)

  n_obs <- terra::global(obs_change, "sum", na.rm = TRUE)[1, 1]
  n_sim <- terra::global(sim_change, "sum", na.rm = TRUE)[1, 1]

  if (n_obs == 0 && n_sim == 0) {
    return(list(
      observed_change = obs_change,
      simulated_change = sim_change,
      similarity = NA_real_,
      similarity_observed_to_simulated = NA_real_,
      similarity_simulated_to_observed = NA_real_,
      similarity_map = NULL,
      n_observed = 0,
      n_simulated = 0
    ))
  }

  weights <- .build_weight_matrix(window_size, use_exp_decay, decay_divisor)
  obs_membership <- .change_membership(obs_change, weights)
  sim_membership <- .change_membership(sim_change, weights)

  obs_to_sim_map <- terra::mask(sim_membership, obs_change, maskvalues = c(NA, FALSE))
  sim_to_obs_map <- terra::mask(obs_membership, sim_change, maskvalues = c(NA, FALSE))

  mean_over_changes <- function(membership_map, n_changed) {
    if (n_changed == 0) {
      return(NA_real_)
    }
    terra::global(membership_map, "mean", na.rm = TRUE)[1, 1]
  }
  similarity_observed_to_simulated <- mean_over_changes(obs_to_sim_map, n_obs)
  similarity_simulated_to_observed <- mean_over_changes(sim_to_obs_map, n_sim)

  list(
    observed_change = obs_change,
    simulated_change = sim_change,
    similarity = min(
      similarity_observed_to_simulated,
      similarity_simulated_to_observed,
      na.rm = TRUE
    ),
    similarity_observed_to_simulated = similarity_observed_to_simulated,
    similarity_simulated_to_observed = similarity_simulated_to_observed,
    similarity_map = obs_to_sim_map,
    n_observed = n_obs,
    n_simulated = n_sim
  )
}

#' Distance-decayed membership in a change map
#'
#' @param change_map Logical SpatRaster, `TRUE` where the change happened, see
#'   [create_change_map()]
#' @param weight_matrix Matrix of distance weights, see [.build_weight_matrix()]
#'
#' @return SpatRaster with, for every cell, the weight of the nearest change within the window
#'   (0 if there is none)
#' @keywords internal
.change_membership <- function(change_map, weight_matrix) {
  changed <- terra::ifel(change_map, 1, 0)
  changed <- terra::classify(changed, cbind(NA, 0))
  terra::focal(changed, w = weight_matrix, fun = "max", na.rm = TRUE, fillvalue = 0)
}

#' Build distance weight matrix for the fuzzy similarity window
#'
#' @param window_size Integer, size of window (must be odd)
#' @param use_exp_decay Logical, use exponential decay?
#' @param decay_divisor Numeric, attenuation factor
#'
#' @return Matrix of weights: `exp(-d / decay_divisor)` for Euclidean distance `d` from the
#'   centre cell, or 1 throughout without decay
#' @keywords internal
.build_weight_matrix <- function(window_size, use_exp_decay, decay_divisor) {
  offsets <- seq_len(window_size) - (window_size + 1) / 2
  distance <- sqrt(outer(offsets^2, offsets^2, "+"))
  if (use_exp_decay) exp(-distance / decay_divisor) else distance * 0 + 1
}
