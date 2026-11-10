# Synthetic land use with a known change process

Simulates a square landscape and its land use over several periods,
driven by a stated change process. It is meant for tutorials, examples
and tests: unlike real data or random noise, the change has a signal a
model can learn, and the process that generated it is known.

There are four classes: 1 forest, 2 arable land, 3 urban areas, 4 a lake
that never changes. Two latent drivers shape both the initial map and
its change:

- `accessibility`, high towards the west and south, plus a smooth random
  field;

- `site_quality`, high towards the north and east, plus a smooth random
  field.

In the initial map, urban land clusters where accessibility is high, the
lake sits in a random patch, arable land takes the better sites and
forest the rest. Each later period applies four transitions, whose
per-period probabilities depend on the drivers and on the share of each
class in the 5 x 5 cell neighbourhood:

|  |  |
|----|----|
| transition | logit of the per-period probability |
| arable -\> urban | -5.5 + 7 \* share_urban + 3 \* accessibility |
| forest -\> urban | -7 + 6 \* share_urban + 3 \* accessibility |
| forest -\> arable | -5 + 4 \* site_quality + 3 \* share_arable |
| arable -\> forest | -6 + 4 \* (1 - site_quality) \* (1 - accessibility) + 3 \* share_forest |

A cell undergoes at most one transition per period; the transitions are
drawn in the order of the table. A third, `random_nuisance` field has no
influence on the change; it is returned as a predictor a model should
learn to ignore.

Arable land, with its higher economic input and output, is cleared on
good sites and abandoned to forest on poor, remote ones. All random
draws go through R's random number generator, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) makes a landscape
reproducible.

The grid lies at the origin of Web Mercator (EPSG:3857) by default: a
coordinate system that belongs to no country, with coordinates that
start near zero and, at the equator, units of true metres.

## Usage

``` r
synthetic_landscape(
  n_periods = 3L,
  n_cells = 30L,
  resolution = 100,
  xmin = 0,
  ymin = 0,
  crs = "EPSG:3857"
)
```

## Arguments

- n_periods:

  Integer, the number of maps: the initial map and `n_periods - 1` steps
  of the process.

- n_cells:

  Integer, cells per side of the square grid.

- resolution:

  Numeric, cell size in the units of `crs`.

- xmin, ymin:

  Numeric, the grid's lower left corner in the units of `crs`.

- crs:

  Character, the grid's coordinate reference system, passed to
  [`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html)

## Value

A list of two
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)s
on the same grid:

- `lulc`: one layer per period, named `id_period=1`, `id_period=2`, ...,
  holding the class of each cell;

- `drivers`: the layers `accessibility`, `site_quality` and
  `random_nuisance`, in \[0, 1\].

## Examples

``` r
set.seed(1)
landscape <- synthetic_landscape(n_periods = 3, n_cells = 30)
terra::plot(landscape$lulc, nc = 3)

terra::plot(landscape$drivers, nc = 3)

```
