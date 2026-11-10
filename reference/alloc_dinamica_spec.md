# Dinamica EGO allocator

Allocation with [Dinamica EGO](https://dinamicaego.com/)'s expander and
patcher (`AllocateTransitions`), run as an external program, see the
`install-dinamica` vignette. Each period's inputs are written to a
working directory: the anterior map, one probability map per viable
transition (the adjusted potentials), the transition rates, and the
expansion and patcher tables from `alloc_params_t` (`frac_expander`,
`mean_patch_size`, `patch_size_variance`, `patch_isometry`).

[`check_dinamica()`](https://ethzplus.github.io/evoland.plus/reference/util_dinamica.md)
runs before anything is predicted or written. Dinamica draws its own
random seed, so its realisations cannot be replayed (see
ethzplus/evoland.plus#104).

## Usage

``` r
alloc_dinamica_spec(work_dir = NULL, keep_intermediate = FALSE)
```

## Arguments

- work_dir:

  Character or `NULL`; directory under which each period gets its own
  working directory, `run_<id_run>_period_<id_period_post>`. `NULL`
  (default) uses a temporary directory.

- keep_intermediate:

  Logical; keep the working directories after each period (default
  `FALSE`, remove them).

## Value

An allocator specification for
[`allocate()`](https://ethzplus.github.io/evoland.plus/reference/allocate.md).
