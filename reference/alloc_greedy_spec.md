# Greedy (rank-and-fill) allocator

Deterministic allocation that fills each transition's demanded quantity
with the cells of highest transition potential: rank the candidate
(cell, transition) pairs, then walk the ranking and accept a pair when
its cell has not changed yet and its transition still has quantity left.
Every cell changes at most once, and no transition exceeds its demand.
There is no randomness and no patch geometry, so the result is one map,
the same on every run.

This is the allocator behind lulcc's Ordered model (Fuchs et al. 2013)
and SEALS' rank-and-fill, and a useful comparator for the stochastic
allocators: it concentrates change on the most probable cells, which a
single-map score such as the figure of merit rewards, whereas an
ensemble of unbiased realisations
([`alloc_clumpy_spec()`](https://ethzplus.github.io/evoland.plus/reference/alloc_clumpy_spec.md))
reproduces where change occurs.

Two ways to settle the competition for a cell between transitions
(`arbitration`):

- `"joint"` (default): one ranking of all pairs by adjusted potential,
  across transitions; a cell goes to whichever of its transitions comes
  first.

- `"ordered"`: transitions take their cells in turn, in the order given
  by `order`, each taking its best remaining cells; earlier transitions
  win contested cells. This is the Ordered procedure of Fuchs et al.
  (2013), with transitions in place of classes.

**Tie-breaking is arbitrary.** To be deterministic, the ranking needs a
total order, so pairs with equal potential are ranked by secondary keys
that carry no meaning:

- `"joint"`: by `id_trans` (lower first), then `id_coord`. Where two
  transitions have equal potential at a cell, the one with the lower
  `id_trans` gets it; where equally ranked pairs compete for the last
  cells of a quota, those with the lower `id_coord` are taken. The
  `"joint"` ranking is therefore not symmetric in the transitions:
  renumbering them can change the result wherever potentials tie.

- `"ordered"`: within a transition, by `id_coord`.

Ties are rare for continuous potentials, but common for coarse ones
(e.g. classification trees, or potentials rounded on input), and
`id_coord` usually follows the raster's row order, so tied cells are
taken from the top of the map first.

The allocation works on `id_coord` rows only, not on a raster, so it
applies to any regular tessellation `coords_t` can describe. The
quantity of change is a number of cells, though, so cells of unequal
area or productivity are not weighted (see ethzplus/evoland.plus#91).
The ranking uses the allocation-ready potentials of
[`adjusted_trans_pot_v()`](https://ethzplus.github.io/evoland.plus/reference/evoland_db_views.md),
which are scaled to each transition's rate, so that they are comparable
across transitions.

## Usage

``` r
alloc_greedy_spec(arbitration = c("joint", "ordered"), order = NULL)
```

## Arguments

- arbitration:

  Character, `"joint"` (default) or `"ordered"`; how transitions compete
  for a cell. Both break ties in potential by ids, which is arbitrary;
  see Description.

- order:

  Integer vector of `id_trans`, the priority for
  `arbitration = "ordered"`; viable transitions not listed follow in
  `id_trans` order. Ignored for `"joint"`.

## Value

An allocator specification for
[`allocate()`](https://ethzplus.github.io/evoland.plus/reference/allocate.md).

## References

Fuchs, R., Herold, M., Verburg, P. H., & Clevers, J. G. P. W. (2013). A
high-resolution and harmonized model approach for reconstructing and
analysing historic land changes in Europe. Biogeosciences, 10(3),
1543-1559. https://doi.org/10.5194/bg-10-1543-2013
