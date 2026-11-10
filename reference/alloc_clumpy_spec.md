# CLUMPY allocator

Stochastic, unbiased allocation after Mazy (2022, Ch. 3), reimplemented
from [CLUMPY](https://github.com/mmyrte/clumpy). The whole
pivot-selection and patch-growth routine runs in C++
([`allocate_clumpy_cpp()`](https://ethzplus.github.io/evoland.plus/reference/allocate_clumpy_cpp.md)).
The method is chosen from the patch parameters in `alloc_params_t`:

- **uSAM** (Unbiased Simple Allocation Method, Mazy sec. 3.4.1) when
  every transition is mono-pixel (`mean_patch_size == 1` and
  `patch_size_variance == 0`): one MuST (Multinomial Sampling Test, Mazy
  App. 3.B; the same test the reference `clumpy` calls "GART") pass per
  anterior class, each selected pivot allocated as a single cell.
  Quantity of change is enforced in expectation.

- **uPAM** (Unbiased Patch Allocation Method, Mazy sec. 3.4.2, Fig. 3.2)
  otherwise: iterative MuST with a per-transition pixel quota and
  sampling without replacement. Affordable here because evoland.plus's
  potentials come from a fixed fitted model, so the marginal density
  does not need to be re-estimated between patches.

Multi-pixel patches require uPAM; "uSAM with patches larger than one
pixel" is not a valid method, hence the automatic selection rather than
a user switch.

The per-cell pivot probability is divided by the mean patch area (the
1/E(sigma) factor, Mazy Fig. 3.2) so that the allocated quantity of
change matches the target transition rate; without it allocation
over-shoots by roughly the mean patch size.

The patch parameters are read from `alloc_params_t`: `mean_patch_size`
and `patch_size_variance` give the patch-area distribution,
`patch_elongation` its shape.

The allocator works on a raster built from `coords_t`, so it needs
coordinate points on a regular grid. R's random number generator drives
it, so [`set.seed()`](https://rdrr.io/r/base/Random.html) makes a
realisation replay.

## Usage

``` r
alloc_clumpy_spec(
  area_dist = c("lognormal", "normal"),
  avoid_aggregation = TRUE,
  batch_size = 0L,
  default_elongation = 0.5
)
```

## Arguments

- area_dist:

  Character; patch-area distribution, `"lognormal"` (default) or
  `"normal"` (Gaussian with sd = `sqrt(patch_size_variance)`, clamped to
  \>= 1).

- avoid_aggregation:

  Logical; if `TRUE` (default) uPAM patches that would merge with an
  existing patch fail and allocate nothing (clumpy `GaussianPatcher`
  semantics). Ignored for the mono-pixel uSAM path.

- batch_size:

  Integer; uPAM pivots attempted per MuST re-draw. `0` (default)
  auto-scales to ~1% of each class's source pool, bounding the number of
  MuST passes so large rasters stay tractable; `> 0` is an explicit cap
  (1 = strict uPAM); `< 0` processes all candidates in a single pass.
  Ignored for uSAM.

- default_elongation:

  Numeric in \[0, 1\]; patch elongation for transitions whose
  `patch_elongation` is missing (0 isometric, 1 linear).

## Value

An allocator specification for
[`allocate()`](https://ethzplus.github.io/evoland.plus/reference/allocate.md).

## References

Mazy, 2022 (<https://theses.hal.science/tel-04382012v1>), Ch. 3.
