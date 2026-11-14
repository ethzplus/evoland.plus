#' Allocate land use change
#'
#' @description
#' *Allocation* places the projected quantity of change on the map: for each period, it turns
#' the land use at the start of the period (the *anterior* map) into the land use at its end
#' (the *posterior* map), given the transition potentials and the demanded transition rates.
#'
#' `allocate()` runs a contiguous sequence of periods. For each period it
#'
#' 1. predicts the transition potentials if `trans_pot_t` does not hold them yet, see
#'    [predict_trans_pot()];
#' 2. builds the allocation context: the anterior map, the allocation-ready potentials
#'    ([adjusted_trans_pot_v()]), the viable transitions, the quantity of change per
#'    transition and the allocation parameters;
#' 3. hands the context to the allocator, which returns the posterior map;
#' 4. commits the map to `lulc_data_t`, together with the neighbor predictors the next period
#'    is predicted from, in one transaction.
#'
#' `allocate_one_period()` runs steps 1-3 for one period and returns the map without
#' committing it. Use it to write your own loop: edit `trans_pot_t` between prediction and
#' allocation (e.g. theory-driven shifts in potential due to policy interventions), edit the
#' allocated map afterwards, or draw an ensemble of single-period realisations. Commit the map
#' with `db$commit(x, "lulc_data_t", method = "upsert")` and, if a later period follows, run
#' `db$upsert_new_neighbors(id_period_post)`.
#'
#' The allocator is one of [alloc_clumpy_spec()] (stochastic, unbiased; the default),
#' [alloc_greedy_spec()] (deterministic rank-and-fill) or [alloc_dinamica_spec()] (Dinamica
#' EGO's expander and patcher).
#'
#' `allocate()` seeds R's random number generator with the run's `rng_seed` in [runs_t] for
#' the whole sequence of periods, and restores the generator's previous state afterwards, so
#' that a run replays from the database. A run without a seed draws one from the generator's
#' current state and records it. Only allocators that draw from R's generator are seeded
#' (CLUMPY); Dinamica draws its own seed, see ethzplus/evoland.plus#104. The first period a
#' run allocates is recorded as its `id_period_eval_start`. `allocate_one_period()` neither
#' seeds nor records these: seed and record them in your own loop.
#'
#' The quantity of change is the same for every allocator: each viable transition changes
#' `round(rate * n_anterior)` cells, where `rate` is the transition's rate in `trans_rates_t`
#' for the posterior period and `n_anterior` the number of cells of its anterior class in the
#' anterior map. Transitions without a rate change no cells.
#'
#' @param db An [evoland_db] instance.
#' @param id_periods Integer vector of contiguous posterior period IDs to simulate.
#' @param id_period_post Integer posterior period ID.
#' @param allocator An allocator specification, see Description. `NULL` (default) uses the
#'   allocator the run records in [runs_t] (inherited from its ancestors), and CLUMPY if none
#'   does. The allocator used is recorded on the run if it records none itself.
#' @param select_score Character; mlr3 measure ID (e.g. `"classif.auc"`) used to select the
#'   model for prediction, see [predict_trans_pot()]. `NULL` (default) uses the criterion the
#'   run of the potentials records in [runs_t].
#' @param select_maximize Logical; maximize (`TRUE`) or minimize (`FALSE`) `select_score`.
#' @param id_run Integer; the run to allocate. It reads its anterior map, demand and
#'   allocation parameters through its lineage, and the allocated maps are written under it.
#' @param potentials_ancestor Integer; how many generations up the lineage of `id_run` to
#'   take the transition potentials from: `0` (default) the run itself, `1` its parent, `2` its
#'   grandparent. They are predicted under that run if missing. Use `1` to let sibling runs
#'   (e.g. a Monte Carlo ensemble of allocations) share their parent's potentials. `2` only
#'   makes sense if the parent's siblings share the same transition potentials too. Above `0`
#'   for a single period only: from the second period on, the anterior map is the run's own
#'   allocation, which its ancestors cannot see.
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
#'
#' @return `allocate()`: called for its side effects on `lulc_data_t` and `pred_data_t`;
#'   returns `NULL` invisibly. `allocate_one_period()`: an [lulc_data_t] with the posterior
#'   map.
#'
#' @name allocate
#' @include evoland_db.R
NULL

#' @rdname allocate
#' @export
allocate <- function(
  db,
  id_periods,
  allocator = NULL,
  select_score = NULL,
  select_maximize = NULL,
  id_run = db$id_run,
  potentials_ancestor = 0L,
  force_predict = FALSE,
  update_neighbors = TRUE,
  checkpoint = FALSE
) {
  allocator <- resolve_allocator(db, id_run, allocator)
  check_allocation_args(allocator, id_run, potentials_ancestor)
  stopifnot(
    "checkpoint must be TRUE or FALSE" = is_scalar_logical(checkpoint),
    "id_periods must be integers" = is_integerish(id_periods),
    "id_periods must be contiguous" = all(diff(id_periods) == 1L),
    "id_periods must be in periods_t" = all(id_periods %in% db$periods_t[["id_period"]]),
    "potentials_ancestor above 0 works for a single period only" = {
      length(id_periods) == 1L || potentials_ancestor == 0L
    }
  )

  message(glue::glue(
    "Allocating with {allocator$name}: run {id_run}, ",
    "periods {paste(id_periods, collapse = ' -> ')}"
  ))

  # the seed and the period the run branches off at belong to the run itself
  own <- list(
    rng_seed = get_run_param(db, id_run, "rng_seed", inherit = FALSE),
    id_period_eval_start = get_run_param(db, id_run, "id_period_eval_start", inherit = FALSE)
  )
  to_record <- list()
  if (is.null(own$id_period_eval_start)) {
    to_record$id_period_eval_start <- as.integer(min(id_periods))
  }
  seed <- own$rng_seed
  if (allocator$seeded && is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
    to_record$rng_seed <- seed
  }
  if (length(to_record) > 0L) {
    do.call(set_run_params, c(list(self = db, id_run = id_run), to_record))
  }

  with_seed(
    if (allocator$seeded) seed,
    allocate_periods(
      db,
      id_periods,
      allocator,
      select_score,
      select_maximize,
      id_run,
      potentials_ancestor,
      force_predict,
      update_neighbors
    )
  )

  if (checkpoint) {
    db$checkpoint()
  }

  invisible(NULL)
}

# the loop of allocate(), one period after the other
allocate_periods <- function(
  db,
  id_periods,
  allocator,
  select_score,
  select_maximize,
  id_run,
  potentials_ancestor,
  force_predict,
  update_neighbors
) {
  for (id_period_post in id_periods) {
    posterior <- allocate_one_period(
      db = db,
      id_period_post = id_period_post,
      allocator = allocator,
      select_score = select_score,
      select_maximize = select_maximize,
      id_run = id_run,
      potentials_ancestor = potentials_ancestor,
      force_predict = force_predict
    )

    # the map and the neighbor counts derived from it land together or not at all
    db$transaction({
      db$commit(posterior, "lulc_data_t", method = "upsert")
      if (update_neighbors || id_period_post != max(id_periods)) {
        db$upsert_new_neighbors(id_period_post, id_run = id_run)
      }
    })
  }

  invisible(NULL)
}

#' @rdname allocate
#' @export
allocate_one_period <- function(
  db,
  id_period_post,
  allocator = NULL,
  select_score = NULL,
  select_maximize = NULL,
  id_run = db$id_run,
  potentials_ancestor = 0L,
  force_predict = FALSE
) {
  allocator <- resolve_allocator(db, id_run, allocator)
  check_allocation_args(allocator, id_run, potentials_ancestor)
  stopifnot(
    "id_period_post must be a single integer" = is_scalar_integerish(id_period_post)
  )

  # the run potentials_ancestor generations up the lineage, the run itself for 0
  potentials_run <- get_lineage(db$runs_t, as.integer(id_run))[potentials_ancestor + 1L]
  if (is.na(potentials_run)) {
    stop(
      glue::glue("Run {id_run} has no ancestor {potentials_ancestor} generations up"),
      call. = FALSE
    )
  }

  db$predict_trans_pot(
    id_period_post = id_period_post,
    select_score = select_score,
    select_maximize = select_maximize,
    force = force_predict,
    id_run = potentials_run
  )

  ctx <- allocation_context(db, id_period_post, id_run, potentials_run)

  message(glue::glue(
    "  {allocator$name}: period {id_period_post - 1L} -> {id_period_post}"
  ))
  posterior <- allocator$one_period(ctx)

  n_changed <- ctx$anterior[posterior, on = "id_coord", sum(id_lulc != i.id_lulc)]
  message(glue::glue("  Changed {n_changed} of {nrow(posterior)} cells"))

  posterior[, .(id_run = id_run, id_coord, id_lulc, id_period = id_period_post)] |>
    as_lulc_data_t()
}

check_allocation_args <- function(allocator, id_run, potentials_ancestor) {
  stopifnot(
    "allocator must be an allocator specification, e.g. alloc_clumpy_spec()" = {
      inherits(allocator, "allocator_spec")
    },
    "id_run must be a single integer" = is_scalar_integerish(id_run),
    "potentials_ancestor must be a single integer, 0 or more" = {
      is_scalar_integerish(potentials_ancestor) && potentials_ancestor >= 0L
    }
  )
  # e.g. Dinamica's setup, before anything is predicted or written
  allocator$check()
  invisible(NULL)
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
#' - `id`: the allocator, e.g. `"clumpy"` for [alloc_clumpy_spec()]; recorded in [runs_t]
#'   together with `args`, from which [allocate()] can recreate the specification;
#' - `args`: the arguments it was created with;
#' - `seeded`: whether it draws from R's random number generator, which [allocate()] seeds
#'   from the run's `rng_seed`;
#' - `check()`: called by [allocate()] before anything else, e.g. to check that an external
#'   program is installed;
#' - `one_period(ctx)`: takes the allocation context [allocate()] builds and returns the
#'   posterior map as a data.table with `id_coord` and `id_lulc`, one row per row of
#'   `ctx$anterior`.
#'
#' @param name Character, label for messages.
#' @param id Character, the allocator: `alloc_<id>_spec()` recreates the specification.
#' @param args Named list, the specification's arguments.
#' @param seeded Logical, whether the allocator draws from R's random number generator.
#' @param one_period Function of the allocation context, see Details.
#' @param check Function without arguments, see Details.
#' @param x An allocator specification.
#' @param ... Ignored.
#'
#' @name allocator_spec
#' @keywords internal
new_allocator_spec <- function(
  name,
  id,
  args,
  one_period,
  check = function() invisible(NULL),
  seeded = FALSE
) {
  stopifnot(
    "name must be a single string" = is_scalar_character(name),
    "id must be a single string" = is_scalar_character(id),
    "seeded must be TRUE or FALSE" = is_scalar_logical(seeded),
    "args must be a list" = is.list(args),
    "one_period must be a function of one argument, the allocation context" = {
      is.function(one_period) && length(formals(one_period)) == 1L
    },
    "check must be a function without arguments" = {
      is.function(check) && length(formals(check)) == 0L
    }
  )
  structure(
    list(
      name = name,
      id = id,
      args = args,
      one_period = one_period,
      check = check,
      seeded = seeded
    ),
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

# The allocator to use for `id_run`: the given one, or the one the run records (inherited),
# or CLUMPY. Recorded on the run if its own row holds none.
resolve_allocator <- function(db, id_run, allocator) {
  if (is.null(allocator)) {
    recorded <- get_run_param(db, id_run, "alloc_args")
    allocator <- if (is.null(recorded)) {
      alloc_clumpy_spec()
    } else {
      spec_fun <- get(paste0("alloc_", recorded[["allocator"]], "_spec"), mode = "function")
      do.call(spec_fun, decode_run_args(recorded[setdiff(names(recorded), "allocator")]))
    }
  }
  stopifnot(
    "allocator must be an allocator specification, e.g. alloc_clumpy_spec()" = {
      inherits(allocator, "allocator_spec")
    }
  )
  if (is.null(get_run_param(db, id_run, "alloc_args", inherit = FALSE))) {
    set_run_params(
      db,
      id_run,
      alloc_args = c(list(allocator = allocator$id), encode_run_args(allocator$args))
    )
  }
  allocator
}

# Evaluate `expr` with R's random number generator seeded with `seed`, and restore the
# generator's previous state afterwards, so that seeding a run does not fix the draws of
# whatever the caller does next. A NULL seed evaluates `expr` as it is.
with_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(expr)
  }
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  previous_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
  on.exit(
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = globalenv())
    } else {
      rm(".Random.seed", envir = globalenv())
    },
    add = TRUE
  )
  set.seed(seed)
  expr
}
