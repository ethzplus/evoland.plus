# Raster grid of a coordinate table

Internal helpers for the code that works on rasters (the CLUMPY and
Dinamica allocators,
[`tabular_to_raster()`](https://ethzplus.github.io/evoland.plus/reference/util_terra.md)):
a raster covering the coordinate points, and the raster cell of each
point, to move values between rows keyed by `id_coord` and raster cells.

- `coords_grid()` builds the grid.

- `coords_cells()` finds the cell of every coordinate point in any
  raster `x`; points outside it get `NA`.
  [`extract_using_coords_t()`](https://ethzplus.github.io/evoland.plus/reference/util_terra.md)
  reads rasters through it.

- `grid_values()` lays values given per coordinate point out over every
  cell of the grid, in terra's cell order; cells without a coordinate
  point, or without a value, are `NA`.

- `grid_to_rast()` does the same, but returns a single-layer raster.

- `rast_to_grid()` reads a raster on the grid back into one value per
  coordinate point.

The cells are only meaningful for a regular grid of points: on any other
set of points, several points can fall into one cell.

## Usage

``` r
coords_grid(coords, resolution = NULL)

coords_cells(coords, x)

grid_values(grid, x, value_col)

grid_to_rast(grid, x, value_col)

rast_to_grid(grid, x)
```

## Arguments

- coords:

  A
  [coords_t](https://ethzplus.github.io/evoland.plus/reference/coords_t.md)
  (or `db$coords_minimal`): `id_coord`, `lon`, `lat` and an `epsg`
  attribute.

- resolution:

  Numeric scalar or `NULL`; the cell size in CRS units. `NULL` uses the
  `resolution` attribute of `coords`, or else the smallest distance
  among its first 1000 points.

- x:

  For `grid_values()` and `grid_to_rast()`, a data.table with `id_coord`
  and `value_col`; for `coords_cells()`, any
  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html);
  for `rast_to_grid()`, a single-layer
  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  on the grid.

- grid:

  A list as returned by `coords_grid()`: `template`, an empty
  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  and `cells`, a data.table with `id_coord` and `cell`.

- value_col:

  Character scalar, the column of `x` holding the values.

## Value

`coords_grid()`: the grid, see `grid`. `coords_cells()`: a data.table
with `id_coord` and `cell`. `grid_values()`: a vector of the type of
`x[[value_col]]`, one element per cell. `grid_to_rast()`: a
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
with one layer named `value_col`. `rast_to_grid()`: a data.table with
`id_coord` and `value`.
