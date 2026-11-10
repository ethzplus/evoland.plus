# Dinamica Utility Functions

Interact with Dinamica from R, see **Functions** section below.

## Usage

``` r
check_dinamica(quiet = FALSE)

exec_dinamica(
  model_path,
  disable_parallel = TRUE,
  log_level = NULL,
  additional_args = NULL,
  write_logfile = TRUE,
  echo = FALSE
)

run_alloc_dinamica(
  work_dir = format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss"),
  ...
)

process_dinamica_script(infile, outfile, mode = "encode", check = TRUE)
```

## Arguments

- quiet:

  bool, suppress the message summarising the result?

- model_path:

  Path to the .ego model file to run. Any submodels must be included in
  a directory of the exact form `basename(modelpath)_ego_Submodels`,
  [see
  wiki](https://csr.ufmg.br/dinamica/dokuwiki/doku.php?id=submodels)

- disable_parallel:

  Whether to disable parallel steps (default TRUE)

- log_level:

  Logging level (1-7, default NULL)

- additional_args:

  Additional arguments to pass to DinamicaConsole, see
  `DinamicaConsole -help`

- write_logfile:

  bool, write stdout&stderr to a file?

- echo:

  bool, direct echo to console?

- work_dir:

  Working dir, where to place ego files and control table

- ...:

  passed to `exec_dinamica()`

- infile:

  Input file path. Treated as input if passed AsIs using
  [`base::I()`](https://rdrr.io/r/base/AsIs.html)

- outfile:

  Output file path (optional)

- mode:

  Character, either "encode" or "decode"

- check:

  Default TRUE, simple check to ensure that you're handling what you're
  expecting

## Value

`check_dinamica()`: a list with character vectors `required` and
`recommended` naming the problems found, invisibly; both are empty for a
complete setup.

## Functions

- `check_dinamica()`: Check that Dinamica EGO is installed and
  configured as the `install-dinamica` vignette describes. Problems that
  stop evoland.plus from running Dinamica at all (platform,
  `DinamicaConsole` on the `PATH`, `DINAMICA_EGO_8_HOME` and its
  libraries) are *required*; those that make a run fail inside Dinamica
  (configuration file, Dinamica's R package) are *recommended*.
  `DINAMICA_EGO_8_TEMP_DIR` is not checked: Dinamica sets it for each R
  process it starts, and a value set by the user breaks that coupling.

- `exec_dinamica()`: Execute a Dinamica .ego file using
  `DinamicaConsole`

- `run_alloc_dinamica()`: Set up evoland-specific Dinamica EGO files;
  execute using `exec_dinamica()`

- `process_dinamica_script()`: Encode or decode raw R and Python code
  chunks in .ego files and their submodels to/from base64
