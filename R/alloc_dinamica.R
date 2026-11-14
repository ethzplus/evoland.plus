#' Dinamica EGO allocator
#'
#' @description
#' Allocation with [Dinamica EGO](https://dinamicaego.com/)'s expander and patcher
#' (`AllocateTransitions`), run as an external program, see the `install-dinamica` vignette.
#' Each period's inputs are written to a working directory: the anterior map, one probability
#' map per viable transition (the adjusted potentials), the transition rates, and the
#' expansion and patcher tables from `alloc_params_t` (`frac_expander`, `mean_patch_size`,
#' `patch_size_variance`, `patch_isometry`).
#'
#' [check_dinamica()] runs before anything is predicted or written.
#'
#' @section Replaying a realisation:
#' DinamicaConsole cannot be given a seed value. By default it seeds its random number
#' generators from the current time, so no two realisations are alike and none can be replayed.
#' Its only seed option, `-predefined-seed`, takes no value: every run then uses the same
#' built-in seed. `replayable = TRUE` passes that option and runs Dinamica single-threaded
#' (`-processors=1 -disable-parallel-functors -disable-parallel-map-load`), because parallel
#' execution has made other Dinamica functors non-deterministic. A replayable run reproduces one
#' realisation; it cannot produce an ensemble of distinct ones, since all runs share the seed.
#'
#' @param work_dir Character or `NULL`; directory under which each period gets its own working
#'   directory, `run_<id_run>_period_<id_period_post>`. `NULL` (default) uses a temporary
#'   directory.
#' @param keep_intermediate Logical; keep the working directories after each period (default
#'   `FALSE`, remove them).
#' @param replayable Logical; run Dinamica with its predefined seed and single-threaded, so that
#'   the same inputs give the same allocation (default `FALSE`, a time-derived seed). See
#'   "Replaying a realisation".
#'
#' @return An allocator specification for [allocate()].
#'
#' @include allocate.R util_dinamica.R
#' @export
alloc_dinamica_spec <- function(work_dir = NULL, keep_intermediate = FALSE, replayable = FALSE) {
  stopifnot(
    "work_dir must be NULL or a single string" = is.null(work_dir) || is_scalar_character(work_dir),
    "keep_intermediate must be TRUE or FALSE" = is_scalar_logical(keep_intermediate),
    "replayable must be TRUE or FALSE" = is_scalar_logical(replayable)
  )
  args <- list(work_dir = work_dir, keep_intermediate = keep_intermediate, replayable = replayable)

  new_allocator_spec(
    name = "Dinamica EGO",
    args = args,
    one_period = function(ctx) alloc_dinamica_one_period(ctx, args),
    check = function() require_dinamica()
  )
}

alloc_dinamica_one_period <- function(ctx, args) {
  work_dir <- ensure_dir(
    if (is.null(args[["work_dir"]])) {
      tempfile("dinamica_")
    } else {
      file.path(args[["work_dir"]], glue::glue("run_{ctx$id_run}_period_{ctx$id_period_post}"))
    }
  )
  if (args[["keep_intermediate"]]) {
    message(glue::glue("  Dinamica working directory kept: {work_dir}"))
  } else {
    on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  }

  write_dinamica_inputs(ctx, work_dir)

  run_alloc_dinamica(
    work_dir = work_dir,
    additional_args = dinamica_replay_args(args[["replayable"]]),
    echo = FALSE,
    write_logfile = TRUE
  )

  posterior <- rast_to_grid(ctx$grid, terra::rast(file.path(work_dir, "posterior.tif")))
  posterior[ctx$anterior[, .(id_coord)], on = "id_coord"][
    !is.na(value),
    .(id_coord, id_lulc = as.integer(value))
  ]
}

# Write the input files for one Dinamica allocation into work_dir: rates, expansion and
# patcher tables, anterior map and one probability map per viable transition. The tables
# and the probability maps follow the same order of transitions, which is how Dinamica
# matches them up.
write_dinamica_inputs <- function(ctx, work_dir) {
  # one row per viable transition, with its rate and patch parameters
  trans <-
    ctx$alloc_params[ctx$viable_trans[ctx$quota, on = "id_trans"], on = "id_trans"] |>
    data.table::setorder(id_lulc_anterior, id_lulc_posterior)

  trans[, .(`From*` = id_lulc_anterior, `To*` = id_lulc_posterior, Rate = rate)] |>
    data.table::fwrite(file.path(work_dir, "trans_rates.csv"))

  trans[, .(
    `From*` = id_lulc_anterior,
    `To*` = id_lulc_posterior,
    # Dinamica cannot do full expansion or full patching, hence constrained to (0, 1)
    Frac_expander = pmax(1e-6, pmin(1 - 1e-6, frac_expander))
  )] |>
    data.table::fwrite(file.path(work_dir, "expansion_table.csv"))

  trans[, .(
    `From*` = id_lulc_anterior,
    `To*` = id_lulc_posterior,
    Mean_Patch_Size = data.table::fcoalesce(mean_patch_size, 1),
    Patch_Size_Variance = data.table::fcoalesce(patch_size_variance, 0),
    Patch_Isometry = data.table::fcoalesce(patch_isometry, 1)
  )] |>
    data.table::fwrite(file.path(work_dir, "patcher_table.csv"))

  grid_to_rast(ctx$grid, ctx$anterior, "id_lulc") |>
    terra::writeRaster(
      file.path(work_dir, "anterior.tif"),
      overwrite = TRUE,
      datatype = "INT1U",
      NAflag = 255 # Dinamica cannot handle NaN
    )

  prob_map_dir <- ensure_dir(file.path(work_dir, "probability_map_dir"))
  for (i in seq_len(nrow(trans))) {
    grid_to_rast(ctx$grid, ctx$potentials[id_trans == trans[["id_trans"]][i]], "value") |>
      terra::writeRaster(
        # prefixed with 001, 002, ... so that the files sort like the tables on every
        # filesystem
        file.path(
          prob_map_dir,
          glue::glue("{sprintf('%03d', i)}_id_trans_{trans[['id_trans']][i]}.tif")
        ),
        overwrite = TRUE,
        NAflag = -999 # Dinamica cannot handle NaN
      )
  }

  invisible(work_dir)
}

# DinamicaConsole flags for a replayable run: its built-in seed, and no multithreading
dinamica_replay_args <- function(replayable) {
  if (!replayable) {
    return(NULL)
  }
  c(
    "-predefined-seed",
    "-processors=1",
    "-disable-parallel-functors",
    "-disable-parallel-map-load"
  )
}
