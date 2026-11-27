#' CLUMPY allocator
#'
#' @description
#' Stochastic, unbiased allocation after Mazy (2022, Ch. 3), reimplemented from
#' [CLUMPY](https://github.com/mmyrte/clumpy). The whole pivot-selection and patch-growth
#' routine runs in C++ ([allocate_clumpy_cpp()]). The method is chosen from the patch
#' parameters in `alloc_params_t`:
#'
#' * **uSAM** (Unbiased Simple Allocation Method, Mazy sec. 3.4.1) when every transition is
#'   mono-pixel (`mean_patch_size == 1` and `patch_size_variance == 0`): one MuST (Multinomial
#'   Sampling Test, Mazy App. 3.B; the same test the reference `clumpy` calls "GART") pass per
#'   anterior class, each selected pivot allocated as a single cell. The quantity of change is
#'   met on average: every cell changes with its own probability, so the number of changed
#'   cells varies between realisations.
#' * **uPAM** (Unbiased Patch Allocation Method, Mazy sec. 3.4.2, Fig. 3.2) otherwise: iterative
#'   MuST with a per-transition pixel quota and sampling without replacement, updating the
#'   probabilities as the pool shrinks (see below). Each pivot grows into a patch (Mazy App.
#'   3.I.3). Patches are allocated until the quota,
#'   `round(rate * n_anterior)` cells, is used up, so the changed cells reach it and exceed it by
#'   less than the last patch. They fall short only where the cells of the anterior class run
#'   out, or patches fail (see `avoid_aggregation`).
#'
#' Multi-pixel patches require uPAM; "uSAM with patches larger than one pixel" is not a valid
#' method, hence the automatic selection rather than a user switch.
#'
#' The per-cell pivot probability is divided by the mean patch area (the 1/E(sigma) factor,
#' Mazy Fig. 3.2) so that the allocated quantity of change matches the target transition rate;
#' without it allocation over-shoots by roughly the mean patch size.
#'
#' **Updating the probabilities.** uPAM draws pivots in repeated MuST passes from a pool that
#' shrinks. With potentials that stay fixed, a cell's potential would act as a hazard rate: by
#' the end, cells of low potential are allocated more often than their potential, and cells of
#' high potential less often. Mazy & Longaretti (2023, \doi{10.1016/j.envsoft.2023.105706},
#' sec. 4.2, App. D.3, E.2.3) update the probabilities after each patch, as an urn per bin
#' that holds the cells still in the pool and the cells still expected to change. Here the bins
#' are bins of the potential, joint over the transitions from one anterior class: cells with
#' the same potential for every transition are exchangeable for all of them, so the update needs
#' neither the predictors nor the model. The cells with a potential are cut in two, again and
#' again, at about the median of the potential mass of the transition whose potentials spread most,
#' never between equal potentials, into at most `n_bins` bins. Each transition expects at least
#' ten patches in each of its bins: where a cut would leave it less on one side, it keeps the
#' uncut group as one bin. With one transition, this gives bins of about equal potential mass,
#' at most one per ten patches of the quota. Per transition, a bin starts with the sum of its
#' potentials as its expected count, scaled so that the bins sum to the quota. A cell's
#' probability is its own potential times the bin's remaining expected count over the remaining
#' sum of potentials of its cells, so that cells of different potential within a bin keep their
#' ratio, and the first pass draws from the potentials as they are. A successful patch counts as
#' cells of its pivot's bin, as in the article's urn; all its cells, and those of failed
#' patches, leave the pool and their bin of every transition. Because the bins are joint, the
#' probability that one transition loses where another takes cells goes back to cells with the
#' same potentials for both, not to cells that only share its own potential.
#' The probabilities are recomputed for every MuST pass, so `batch_size` sets how many pivots
#' share one update.
#'
#' **What uPAM is unbiased about.** The pivot cells are drawn in proportion to the potentials,
#' and so are the patch sizes and elongations if same-step patches do not merge (Mazy App. 3.B.2,
#' 3.D.1.2). The other cells of a patch are not drawn: the grower adds the neighbour with the
#' best score, potential over distance to the target elongation. The allocated cells as a whole
#' follow the potentials only where these vary slowly over the extent of a patch (Mazy App.
#' 3.A.2: patch scale much smaller than the scale of the explanatory variables). Where
#' potentials vary from cell to cell, uPAM places change where the potentials are high in a
#' cell's neighbourhood, not cell by cell; with mono-pixel patches, use uSAM. The update does not
#' change this: counting each cell of a patch against its own bin would steer the pivots so
#' that the patch cells follow the potentials, which the grower's deterministic choice of
#' neighbours defeats, and the pivots would no longer follow the potentials.
#'
#' **Differences from the reference CLUMPY.** Checked against clumpy's patch grower on 2000
#' scenarios (evoland-experiments, `2026-09-paper-figures/010-verify-clumpy.r`):
#' * no connectedness step: clumpy first adds, at random, a cell almost enclosed by the patch
#'   (Mazy App. 3.I.3.2), to avoid small holes;
#' * exact elongation: a candidate is scored by the elongation of the patch including it, from
#'   moments about that patch's own centroid (Mazy eqs. 3.I.6-3.I.8); the thesis (App. 3.I.3.3)
#'   and clumpy keep the centroid of the patch without the candidate, for speed;
#' * patch areas are rounded to the nearest cell, where clumpy grows to the next larger one
#'   (half a cell above the prescribed mean on average);
#' * without `avoid_aggregation`, a patch that runs out of neighbours of its class is kept as
#'   grown, where clumpy discards it;
#' * the probabilities are updated on bins of the potential (Mazy & Longaretti 2023, App. D.3,
#'   E.2.3), where clumpy re-estimates the density of the remaining pool in the space of the
#'   explanatory variables (App. D.4).
#'
#' With the connectedness step and the elongation aligned, both grow the same patch, except where
#' equally scored candidates are separated by floating-point round-off. uSAM matches clumpy bit
#' for bit when clumpy's random numbers are supplied through R's user-supplied generator.
#'
#' The patch parameters are read from `alloc_params_t`: `mean_patch_size` and
#' `patch_size_variance` give the patch-area distribution, `patch_elongation` its shape. Areas
#' are rounded to whole cells and are at least one cell. A transition without a variance (0, or
#' missing because it was observed as a single patch) gets patches of exactly its mean size.
#'
#' The allocator works on a raster built from `coords_t`, so it needs coordinate points on a
#' regular grid. R's random number generator drives it, so [set.seed()] makes a realisation
#' replay.
#'
#' @param area_dist Character; patch-area distribution, `"lognormal"` (default) or `"normal"`
#'   (Gaussian with sd = `sqrt(patch_size_variance)`, clamped to >= 1).
#' @param avoid_aggregation Logical; if `TRUE` (default), a uPAM patch that would merge with
#'   another patch allocated in the same period, or cannot reach its drawn area, fails and
#'   allocates nothing (clumpy `GaussianPatcher` semantics, Mazy App. 3.D.1.2: merging lowers the
#'   quantity of change and shifts the patch sizes towards larger patches). Land that had the
#'   target class before the period does not count as a patch: patches may grow next to it, so
#'   existing areas expand with avoidance on. Failed patches remove their cells from the pool, so
#'   where patches crowd, the quota can be left slightly unmet. `FALSE` lets same-period patches
#'   merge and keeps patches that run out of neighbours. Ignored for the mono-pixel uSAM path.
#' @param batch_size Integer; uPAM pivots attempted per MuST pass, i.e. between two updates of
#'   the probabilities. `0` (default) auto-scales to ~1% of each class's source pool, bounding
#'   the number of MuST passes so large rasters stay tractable; `> 0` is an explicit cap (1 =
#'   strict uPAM, an update after every patch); `< 0` processes all candidates of a pass, and
#'   updates only where one pass falls short of the quota. All three are unbiased; the fewer
#'   pivots share an update, the closer the count allocated in each range of potential comes to
#'   the sum of its potentials, at the cost of more passes. Ignored for uSAM.
#' @param default_elongation Numeric in \[0, 1\]; patch elongation for transitions whose
#'   `patch_elongation` is missing (0 isometric, 1 linear).
#' @param n_bins Integer >= 1; uPAM maximum number of joint bins of potential per anterior class
#'   for the update of the probabilities (default 100), and so of bins per transition. Fewer,
#'   wider bins share one update among more different potentials; bins expecting few patches
#'   each get a whole patch or none. Ignored for uSAM.
#'
#' @return An allocator specification for [allocate()].
#'
#' @references Mazy, 2022 (\url{https://theses.hal.science/tel-04382012v1}), Ch. 3.
#'
#' Mazy, F.-R., Longaretti, P.-Y., 2023. Towards a generic theoretical framework for
#' pattern-based LUCC modeling. Environmental Modelling & Software 166, 105706.
#' \doi{10.1016/j.envsoft.2023.105706}
#'
#' @include allocate.R
#' @export
alloc_clumpy_spec <- function(
  area_dist = c("lognormal", "normal"),
  avoid_aggregation = TRUE,
  batch_size = 0L,
  default_elongation = 0.5,
  n_bins = 100L
) {
  area_dist <- match.arg(area_dist)
  stopifnot(
    "avoid_aggregation must be TRUE or FALSE" = is_scalar_logical(avoid_aggregation),
    "batch_size must be a single integer" = is_scalar_integerish(batch_size),
    "default_elongation must be a number in [0, 1]" = {
      is_scalar_numeric(default_elongation) && default_elongation >= 0 && default_elongation <= 1
    },
    "n_bins must be a single integer, 1 or more" = is_scalar_integerish(n_bins) && n_bins >= 1
  )
  args <- list(
    area_dist = area_dist,
    avoid_aggregation = avoid_aggregation,
    batch_size = as.integer(batch_size),
    default_elongation = default_elongation,
    n_bins = as.integer(n_bins)
  )

  new_allocator_spec(
    name = "CLUMPY",
    args = args,
    one_period = function(ctx) alloc_clumpy_one_period(ctx, args)
  )
}

alloc_clumpy_one_period <- function(ctx, args) {
  grid <- ctx$grid
  cells <- grid[["cells"]]
  viable_trans <- ctx$viable_trans

  # the anterior map as a raster in row-major cell order, NA where there is no coordinate point
  landscape <- grid_values(grid, ctx$anterior, "id_lulc")

  # Sparse per-transition potential columns: the cells carrying a nonzero adjusted potential
  # and their values. The adjusted potentials are already sparse, so they are passed as
  # per-transition lists rather than a dense n_cells x n_trans matrix (prohibitive at >1e7
  # cells).
  # one table per viable transition, in their order, empty for a transition without potentials
  potentials_by_trans <-
    ctx$potentials[value > 0][cells, on = "id_coord", cell := i.cell][,
      # as a factor, a transition without any potential (e.g. at a rate of 0) keeps its element
      id_trans := factor(id_trans, levels = viable_trans[["id_trans"]])
    ] |>
    split(by = "id_trans", drop = FALSE, sorted = TRUE)
  prob_cell <- lapply(potentials_by_trans, `[[`, "cell")
  prob_value <- lapply(potentials_by_trans, `[[`, "value")

  # patch parameters and rates in the order of viable_trans
  params <- ctx$alloc_params[viable_trans[, .(id_trans)], on = "id_trans"]
  rates <- ctx$quota[viable_trans[, .(id_trans)], on = "id_trans"]
  area_mean <- as.numeric(params[["mean_patch_size"]])
  area_var <- as.numeric(params[["patch_size_variance"]])
  elongation <- data.table::fcoalesce(
    as.numeric(params[["patch_elongation"]]),
    args[["default_elongation"]]
  )

  # every transition mono-pixel -> uSAM, otherwise uPAM
  is_mono <- all(!is.na(area_mean) & area_mean == 1 & (is.na(area_var) | area_var == 0))
  message(glue::glue("  CLUMPY method: {if (is_mono) 'uSAM' else 'uPAM'}"))

  posterior <- allocate_clumpy_cpp(
    landscape = landscape,
    nrow = terra::nrow(grid[["template"]]),
    ncol = terra::ncol(grid[["template"]]),
    trans_from = viable_trans[["id_lulc_anterior"]],
    trans_to = viable_trans[["id_lulc_posterior"]],
    prob_cell = prob_cell,
    prob_value = prob_value,
    area_mean = area_mean,
    area_var = area_var,
    elongation = elongation,
    target_rate = rates[["rate"]],
    method = if (is_mono) 0L else 1L,
    batch_size = args[["batch_size"]],
    rarefy = TRUE,
    shuffle = TRUE,
    avoid_aggregation = args[["avoid_aggregation"]],
    # integer code used by allocate_clumpy_cpp: 0 = log-normal, 1 = normal
    area_dist = if (args[["area_dist"]] == "lognormal") 0L else 1L,
    n_bins = args[["n_bins"]]
  )

  cells[ctx$anterior[, .(id_coord)], on = "id_coord"][,
    .(id_coord, id_lulc = as.integer(posterior[cell]))
  ]
}
