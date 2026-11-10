# Allocator specifications

An allocator specification names an allocation algorithm and carries its
arguments; pass it to
[`allocate()`](https://ethzplus.github.io/evoland.plus/reference/allocate.md).
The arguments are checked when the specification is created, so that a
mistake surfaces before anything is predicted or written.

## Usage

``` r
new_allocator_spec(name, args, one_period, check = function() invisible(NULL))

# S3 method for class 'allocator_spec'
print(x, ...)
```

## Arguments

- name:

  Character, label for messages.

- args:

  Named list, the specification's arguments.

- one_period:

  Function of the allocation context, see Details.

- check:

  Function without arguments, see Details.

- x:

  An allocator specification.

- ...:

  Ignored.

## Details

A specification is a list of class `allocator_spec` with

- `name`: a label for messages;

- `args`: the arguments it was created with;

- `check()`: called by
  [`allocate()`](https://ethzplus.github.io/evoland.plus/reference/allocate.md)
  before anything else, e.g. to check that an external program is
  installed;

- `one_period(ctx)`: takes the allocation context
  [`allocate()`](https://ethzplus.github.io/evoland.plus/reference/allocate.md)
  builds and returns the posterior map as a data.table with `id_coord`
  and `id_lulc`, one row per row of `ctx$anterior`.
