#' Allocate land use change
#'
#' @description
#' *Allocation* places the projected quantity of change on the map: for each period, it turns
#' the land use at the start of the period (the *anterior* map) into the land use at its end
#' (the *posterior* map), given the transition potentials and the demanded transition rates.
#'
#' `allocate()` runs a contiguous sequence of periods for one or more runs. For each period it
#'
#' 1. predicts the transition potentials if `trans_pot_t` does not hold them yet, see
#'    [predict_trans_pot()];
#' 2. builds the allocation context: the anterior map, the allocation-ready potentials
#'    ([adjusted_trans_pot_v()]), the viable transitions, the quantity of change per
#'    transition and the allocation parameters;
#' 3. hands the context to the allocator, once per run, which returns the posterior map;
#' 4. commits the maps of all runs to `lulc_data_t`, together with the neighbor predictors the
#'    next period is predicted from, in one transaction.
#'
#' `allocate_one_period()` runs steps 1-3 for one period and returns the maps without
#' committing them. Use it to write your own loop: edit `trans_pot_t` between prediction and
#' allocation (e.g. theory-driven shifts in potential due to policy interventions), edit the
#' allocated map afterwards, or draw an ensemble of single-period realisations. Commit the map
#' with `db$commit(x, "lulc_data_t", method = "upsert")` and, if a later period follows, run
#' `db$upsert_new_neighbors(id_period_post)`.
#'
#' The allocator is one of [alloc_clumpy_spec()] (stochastic, unbiased; the default),
#' [alloc_greedy_spec()] (deterministic rank-and-fill) or [alloc_dinamica_spec()] (Dinamica
#' EGO's expander and patcher). It is not recorded: to document it, add a column to
#' [runs_t].
#'
#' # Runs, branches and seeds
#'
#' `allocate()` takes a vector of runs, e.g. an ensemble of stochastic realisations, or runs
#' that branch off one another. It reads two columns of [runs_t]:
#'
#' - `id_period_eval_start`: the period a run *branches off* its parent. The run allocates the
#'   periods of `id_periods` from this one on, and reads the earlier ones from its lineage.
#'   Without it, the run allocates every period of `id_periods`.
#' - `rng_seed`: before allocating a period, R's default random number generator is seeded
#'   with `rng_seed + id_period`, so that a run's maps replay whichever runs it is allocated
#'   with, and on whichever process.
#'   Runs without a seed draw from the caller's generator. Seeding does not change the state
#'   of the caller's generator.
#'
#' Periods are allocated one after the other, and each for all runs that have reached it, so
#' that a branch finds the maps of the run it branches off. In the period it branches off, a
#' run allocates on its parent's transition potentials, unless it holds its own: predictors
#' before that period, the anterior map or potentials. Branches that also hold no
#' rates or allocation parameters of their own share their parent's allocation context: it is
#' built once, and each of them allocated on it in turn. The runs of a period are spread over
#' the mirai daemons of `compute`, see [map_evoland()], and their maps committed once all are
#' done; transition potentials are predicted beforehand, on the calling process.
#'
#' The quantity of change is the same for every allocator: each viable transition changes
#' `round(rate * n_anterior)` cells, where `rate` is the transition's rate in `trans_rates_t`
#' for the posterior period and `n_anterior` the number of cells of its anterior class in the
#' anterior map. Transitions without a rate change no cells.
#'
#' @param db An [evoland_db] instance.
#' @param id_periods Integer vector of contiguous posterior period IDs to simulate.
#' @param id_period_post Integer posterior period ID.
#' @param allocator An allocator specification, see Description.
#' @param select_score Character; mlr3 measure ID (e.g. `"classif.auc"`) used to select the
#'   model for prediction, see [predict_trans_pot()].
#' @param select_maximize Logical; maximize (`TRUE`) or minimize (`FALSE`) `select_score`.
#' @param id_run Integer vector of distinct runs to allocate. Each reads its anterior map,
#'   demand and allocation parameters through its lineage, and its allocated maps are written
#'   under it. For `allocate_one_period()`, several runs must share their ancestor
#'   `potentials_ancestor` generations up, and are all allocated on the context read under that
#'   ancestor: their own rates and allocation parameters are not read.
#' @param potentials_ancestor Integer; how many generations up the lineage of `id_run` to
#'   take the transition potentials from: `0` (default) the run itself, `1` its parent, `2` its
#'   grandparent. They are predicted under that run if missing. Use `1` to let sibling runs
#'   (e.g. a Monte Carlo ensemble of allocations) share their parent's potentials. `2` only
#'   makes sense if the parent's siblings share the same transition potentials too. `allocate()`
#'   takes the parent's potentials where a run branches off, see Description.
#' @param force_predict Logical; if `TRUE`, predict the transition potentials even if
#'   `trans_pot_t` already holds them for that run and the period.
#' @param update_neighbors Logical; whether to recompute the neighbor predictors after the last
#'   period (default `TRUE`). Intermediate periods are always updated, since the next period's
#'   prediction reads them. Set `FALSE` when nothing is allocated onwards from the last period,
#'   e.g. for single-period ensembles.
#' @param checkpoint Logical; if `TRUE`, run the database's `$checkpoint()` once, after the
#'   last period (default `FALSE`). Every allocated period adds a snapshot to the catalog,
#'   and DuckLake neither expires snapshots nor merges small files on its own, so a study
#'   that runs many realisations grows until it is checkpointed. Set it on the last call of
#'   a batch, or call `db$checkpoint()` once the batch is done. Leave it off when other
#'   processes allocate into the same database at the same time: a checkpoint rewrites files
#'   they may be reading. How many snapshots it expires depends on the retention set on the
#'   database, see [ducklake_db].
#' @param compute Name of the mirai compute profile to allocate the runs of a period on, see
#'   [map_evoland()].
#'
#' @return `allocate()`: called for its side effects on `lulc_data_t` and `pred_data_t`;
#'   returns `NULL` invisibly. `allocate_one_period()`: an [lulc_data_t] with the posterior
#'   map of each run.
#'
#' @name allocate
#' @include evoland_db.R
NULL

#' @rdname allocate
#' @export
allocate <- function(
  db,
  id_periods,
  allocator = alloc_clumpy_spec(),
  select_score,
  select_maximize,
  id_run = 0L,
  force_predict = FALSE,
  update_neighbors = TRUE,
  checkpoint = FALSE,
  compute = "default"
) {
  check_allocator(allocator)
  stopifnot(
    "checkpoint must be TRUE or FALSE" = is_scalar_logical(checkpoint),
    "id_periods must be integers" = is_integerish(id_periods),
    "id_periods must be contiguous" = all(diff(id_periods) == 1L),
    "id_periods must be in periods_t" = all(id_periods %in% db$periods_t[["id_period"]]),
    "id_run must be distinct integers" = {
      length(id_run) > 0L && is_integerish(id_run) && !anyDuplicated(id_run)
    }
  )

  schedule <- allocation_schedule(db, as.integer(id_periods), as.integer(id_run))
  for (run in schedule[, .(periods = paste(id_period, collapse = " -> ")), by = id_run][,
    paste0("run ", id_run, ", periods ", periods)
  ]) {
    message(glue::glue("Allocating with {allocator$name}: {run}"))
  }
  n_workers <- n_evoland_workers(compute)

  for (id_period_post in id_periods) {
    todo <- schedule[id_period == id_period_post]

    # on the calling process, so that runs sharing potentials do not predict them twice
    for (potentials_run in unique(todo[["potentials_run"]])) {
      db$predict_trans_pot(
        id_period_post = id_period_post,
        select_score = select_score,
        select_maximize = select_maximize,
        force = force_predict,
        id_run = potentials_run
      )
    }

    # the runs sharing a context, in as many items as there are workers
    items <- todo |>
      split(by = c("context_run", "potentials_run")) |>
      lapply(function(group) {
        chunks <- min(n_workers, nrow(group))
        split(group, rep_len(seq_len(chunks), nrow(group)))
      }) |>
      unlist(recursive = FALSE)

    posterior <-
      map_evoland(
        items = unname(items),
        worker_fun = allocate_on_context,
        parent_db = db,
        compute = compute,
        id_period_post = id_period_post,
        allocator = allocator
      ) |>
      data.table::rbindlist() |>
      as_lulc_data_t()

    # the maps and the neighbor counts derived from them land together or not at all
    db$transaction({
      db$commit(posterior, "lulc_data_t", method = "upsert")
      if (update_neighbors || id_period_post != max(id_periods)) {
        for (run in todo[["id_run"]]) {
          db$upsert_new_neighbors(id_period_post, id_run = run)
        }
      }
    })
  }

  if (checkpoint) {
    db$checkpoint()
  }

  invisible(NULL)
}

#' @rdname allocate
#' @export
allocate_one_period <- function(
  db,
  id_period_post,
  allocator = alloc_clumpy_spec(),
  select_score,
  select_maximize,
  id_run = 0L,
  potentials_ancestor = 0L,
  force_predict = FALSE
) {
  check_allocator(allocator)
  stopifnot(
    "id_period_post must be a single integer" = is_scalar_integerish(id_period_post),
    "id_run must be distinct integers" = {
      length(id_run) > 0L && is_integerish(id_run) && !anyDuplicated(id_run)
    },
    "potentials_ancestor must be a single integer, 0 or more" = {
      is_scalar_integerish(potentials_ancestor) && potentials_ancestor >= 0L
    }
  )
  id_run <- as.integer(id_run)
  runs <- db$runs_t

  # the run potentials_ancestor generations up the lineage, the run itself for 0
  potentials_run <- vapply(
    id_run,
    function(run) get_lineage(runs, run)[potentials_ancestor + 1L],
    integer(1L)
  )
  if (anyNA(potentials_run)) {
    stop(
      glue::glue(
        "Run {id_run[is.na(potentials_run)][1L]} has no ancestor ",
        "{potentials_ancestor} generations up"
      ),
      call. = FALSE
    )
  }
  potentials_run <- unique(potentials_run)
  stopifnot(
    "several runs must share their ancestor potentials_ancestor generations up" = {
      length(potentials_run) == 1L
    }
  )

  db$predict_trans_pot(
    id_period_post = id_period_post,
    select_score = select_score,
    select_maximize = select_maximize,
    force = force_predict,
    id_run = potentials_run
  )

  allocate_on_context(
    item = data.table::data.table(
      id_run = id_run,
      rng_seed = runs[["rng_seed"]][match(id_run, runs[["id_run"]])],
      # a single run reads its own anterior map, rates and allocation parameters
      context_run = if (length(id_run) == 1L) id_run else potentials_run,
      potentials_run = potentials_run
    ),
    db = db,
    id_period_post = id_period_post,
    allocator = allocator
  )
}

check_allocator <- function(allocator) {
  stopifnot(
    "allocator must be an allocator specification, e.g. alloc_clumpy_spec()" = {
      inherits(allocator, "allocator_spec")
    }
  )
  # e.g. Dinamica's setup, before anything is predicted or written
  allocator$check()
  invisible(NULL)
}

# Which run allocates which period, on whose potentials and allocation context: one row per
# run and period, ordered by period. A run allocates from where it branches off
# (id_period_eval_start) on. There, it takes its parent's potentials unless it holds inputs of
# its own: an anterior map, predictors before the period or potentials. If it holds no rates
# or allocation parameters of its own either, it shares its parent's whole context.
allocation_schedule <- function(db, id_periods, id_run) {
  # data.table expressions below would read id_run as the column
  run_ids <- id_run
  runs <- db$runs_t
  unknown <- setdiff(run_ids, runs[["id_run"]])
  if (length(unknown) > 0L) {
    stop(glue::glue("Runs not in runs_t: {toString(unknown)}"), call. = FALSE)
  }

  runs <- runs[(match(run_ids, runs[["id_run"]]))]
  runs[, first := pmax(min(id_periods), id_period_eval_start, na.rm = TRUE)]
  late <- runs[["first"]] > max(id_periods)
  if (any(late)) {
    stop(
      glue::glue(
        "Runs {toString(run_ids[late])} branch off after the last of id_periods ",
        "(id_period_eval_start {toString(runs[['id_period_eval_start']][late])})"
      ),
      call. = FALSE
    )
  }

  schedule <- runs[,
    .(id_period = seq.int(first, max(id_periods))),
    by = .(id_run, parent_id_run, rng_seed, id_period_eval_start, first)
  ][, `:=`(potentials_run = id_run, context_run = id_run)]

  branching <- schedule[id_period == id_period_eval_start & !is.na(parent_id_run)]
  for (period in unique(branching[["id_period"]])) {
    branch_runs <- branching[id_period == period, id_run]
    own_inputs <- c(
      runs_with_own_rows(db, "lulc_data_t", branch_runs, paste("id_period =", period - 1L)),
      runs_with_own_rows(db, "pred_data_t", branch_runs, paste("id_period <", period)),
      runs_with_own_rows(db, "trans_pot_t", branch_runs, paste("id_period_post =", period))
    )
    own_context <- c(
      own_inputs,
      runs_with_own_rows(db, "trans_rates_t", branch_runs, paste("id_period =", period)),
      runs_with_own_rows(db, "alloc_params_t", branch_runs)
    )
    schedule[
      id_period == period & id_run %in% setdiff(branch_runs, own_inputs),
      potentials_run := parent_id_run
    ]
    schedule[
      id_period == period & id_run %in% setdiff(branch_runs, own_context),
      context_run := parent_id_run
    ]
  }

  data.table::setorder(schedule, id_period)
  schedule <- schedule[, .(id_period, id_run, rng_seed, potentials_run, context_run)]
  data.table::setattr(schedule, "class", c("data.table", "data.frame"))
}

# Those of `id_run` that hold rows of their own, matching `where`, in a table
runs_with_own_rows <- function(db, table_name, id_run, where = "true") {
  if (!table_name %in% db$list_tables()) {
    return(integer(0L))
  }
  db$get_query(
    "select distinct id_run from dl_db.{`table_name`}
     where id_run in ({id_run*}) and ({DBI::SQL(where)})",
    as_atomic = TRUE
  )
}

# Allocate each run of `item` (a data.table with id_run, rng_seed, context_run and
# potentials_run, the last two the same for all rows) on one allocation context; the worker
# of allocate(). Returns the posterior maps of all runs as one lulc_data_t.
allocate_on_context <- function(item, db, id_period_post, allocator) {
  ctx <- allocation_context(
    db,
    id_period_post,
    item[["context_run"]][1L],
    item[["potentials_run"]][1L]
  )

  # seeded runs leave the caller's generator as it was; unseeded ones draw from it in turn
  caller_rng <- get_rng_state()
  on.exit(set_rng_state(caller_rng), add = TRUE)

  posteriors <- lapply(seq_len(nrow(item)), function(i) {
    run <- item[["id_run"]][i]
    seed <- item[["rng_seed"]][i]
    if (is.na(seed)) {
      set_rng_state(caller_rng)
    } else {
      # R's default generator, whichever the process uses: mirai daemons use L'Ecuyer-CMRG
      set.seed(
        (seed + id_period_post) %% .Machine$integer.max,
        kind = "default",
        normal.kind = "default",
        sample.kind = "default"
      )
    }

    message(glue::glue(
      "  {allocator$name}: run {run}, period {id_period_post - 1L} -> {id_period_post}"
    ))
    ctx$id_run <- run
    posterior <- allocator$one_period(ctx)
    if (is.na(seed)) {
      caller_rng <<- get_rng_state()
    }

    n_changed <- ctx$anterior[posterior, on = "id_coord", sum(id_lulc != i.id_lulc)]
    message(glue::glue("  Changed {n_changed} of {nrow(posterior)} cells"))

    posterior[, .(id_run = run, id_coord, id_lulc, id_period = id_period_post)]
  })

  as_lulc_data_t(data.table::rbindlist(posteriors))
}

get_rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
}

set_rng_state <- function(state) {
  if (!is.null(state)) {
    assign(".Random.seed", state, envir = globalenv())
  } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    # R seeds a fresh generator on the next draw
    rm(".Random.seed", envir = globalenv())
  }
}

# Everything an allocator reads, gathered in one place. Allocators get this list and never
# touch the database themselves. `grid` serves the raster-based allocators (CLUMPY, Dinamica),
# see coords_grid().
allocation_context <- function(db, id_period_post, id_run, potentials_run) {
  id_period_anterior <- id_period_post - 1L

  anterior <- db$get_query(
    "select id_coord, id_lulc
     from {db$get_read_expr('lulc_data_t', id_run = id_run)}
     where id_period = {id_period_anterior}
     order by id_coord"
  )
  if (nrow(anterior) == 0L) {
    stop(
      glue::glue("No LULC data for run {id_run} in the anterior period {id_period_anterior}"),
      call. = FALSE
    )
  }

  viable_trans <-
    viable_trans_meta(db)[, .(id_trans, id_lulc_anterior, id_lulc_posterior)] |>
    data.table::setorder(id_trans)
  stopifnot("No viable transitions found" = nrow(viable_trans) > 0L)

  rates <- db$get_query(
    "select id_trans, rate
     from {db$get_read_expr('trans_rates_t', id_run = id_run)}
     where id_period = {id_period_post}"
  )

  # one quantity of change for every allocator: the rate's share of the anterior class
  quota <-
    rates[viable_trans, on = "id_trans"][
      anterior[, .(n_anterior = .N), by = .(id_lulc_anterior = id_lulc)],
      on = "id_lulc_anterior"
    ][
      viable_trans[, .(id_trans)],
      on = "id_trans"
    ][, `:=`(
      rate = data.table::fcoalesce(rate, 0),
      n_anterior = data.table::fcoalesce(n_anterior, 0L)
    )][, .(id_trans, rate, n_anterior, quota = as.integer(round(rate * n_anterior)))]

  potentials <- db$adjusted_trans_pot_v(
    id_period_post,
    id_run = id_run,
    id_run_trans_pot = potentials_run
  )

  alloc_params <- db$get_query(
    "from {db$get_read_expr('alloc_params_t', id_run = id_run)}"
  )

  list(
    id_run = id_run,
    id_period_post = id_period_post,
    anterior = anterior,
    viable_trans = viable_trans,
    quota = quota,
    potentials = potentials,
    alloc_params = alloc_params,
    grid = coords_grid(db$coords_minimal)
  )
}

#' Allocator specifications
#'
#' @description
#' An allocator specification names an allocation algorithm and carries its arguments; pass it
#' to [allocate()]. The arguments are checked when the specification is created, so that a
#' mistake surfaces before anything is predicted or written.
#'
#' @details
#' A specification is a list of class `allocator_spec` with
#' - `name`: a label for messages;
#' - `args`: the arguments it was created with;
#' - `check()`: called by [allocate()] before anything else, e.g. to check that an external
#'   program is installed;
#' - `one_period(ctx)`: takes the allocation context [allocate()] builds and returns the
#'   posterior map as a data.table with `id_coord` and `id_lulc`, one row per row of
#'   `ctx$anterior`.
#'
#' @param name Character, label for messages.
#' @param args Named list, the specification's arguments.
#' @param one_period Function of the allocation context, see Details.
#' @param check Function without arguments, see Details.
#' @param x An allocator specification.
#' @param ... Ignored.
#'
#' @name allocator_spec
#' @keywords internal
new_allocator_spec <- function(name, args, one_period, check = function() invisible(NULL)) {
  stopifnot(
    "name must be a single string" = is_scalar_character(name),
    "args must be a list" = is.list(args),
    "one_period must be a function of one argument, the allocation context" = {
      is.function(one_period) && length(formals(one_period)) == 1L
    },
    "check must be a function without arguments" = {
      is.function(check) && length(formals(check)) == 0L
    }
  )
  structure(
    list(name = name, args = args, one_period = one_period, check = check),
    class = "allocator_spec"
  )
}

#' @rdname allocator_spec
#' @export
print.allocator_spec <- function(x, ...) {
  cat("Allocator:", x[["name"]], "\n")
  for (arg in names(x[["args"]])) {
    cat("  ", arg, ": ", deparse1(x[["args"]][[arg]]), "\n", sep = "")
  }
  invisible(x)
}
