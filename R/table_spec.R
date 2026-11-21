#' Table Specifications
#'
#' @description
#' Every table of [evoland_db] is declared once, as a `table_spec` object next to its
#' `as_<table>_t()` constructor. The constructor, the inheritance through the run lineage
#' (see [get_evoland_db_read_expr()]), the ER diagram and the column lists in the
#' documentation are all derived from it. The registry of all specs is `EVOLAND_SCHEMA`.
#'
#' A spec lists the columns a table guarantees. Further columns pass through the
#' constructor unless `allow_extra_cols` is `FALSE`.
#'
#' @param name Character scalar, the table name and the class of its objects, e.g.
#'   `"lulc_data_t"`.
#' @param columns Named character vector, column name to type, in the canonical column
#'   order. Types are `"int"`, `"float"`, `"bool"`, `"char"`, `"date"`, `"factor"`,
#'   `"list"` (an R list column, e.g. a serialized object) and `"map"` (a list column of
#'   named lists, stored as a DuckDB `MAP`).
#' @param descriptions Named character vector, column name to its meaning, for the
#'   documentation; one entry for every column in `columns` and `optional_cols`.
#' @param key_cols,alternate_key_cols,partition_cols Character vectors, see
#'   [as_ducklake_db_t()].
#' @param inherit_by Character vector or `NULL`. `NULL` for a table that is not run-scoped.
#'   For a run-scoped table, the columns besides `id_run` that form one inheritance slice:
#'   reading the table for a run returns each slice from the closest run in its lineage
#'   that holds it. `character(0)` makes the whole table one slice.
#' @param optional_cols Named character vector, like `columns`, for columns a table may
#'   have: they are cast when present, but not added when missing.
#' @param levels Named list, factor column name to its levels.
#' @param allow_extra_cols Logical; if `FALSE`, the constructor rejects columns that are
#'   neither in `columns` nor in `optional_cols`.
#' @return A `table_spec` object
#' @name table_spec
#' @keywords internal
#' @include util.R
table_spec <- function(
  name,
  columns,
  descriptions,
  key_cols = NULL,
  alternate_key_cols = NULL,
  partition_cols = NULL,
  inherit_by = NULL,
  optional_cols = NULL,
  levels = list(),
  allow_extra_cols = TRUE
) {
  all_cols <- c(columns, optional_cols)
  stopifnot(
    "name must be a character scalar" = is_scalar_character(name),
    "columns must be a named character vector" = is.character(columns) &&
      !is.null(names(columns)),
    "column types must be one of int, float, bool, char, date, factor, list, map" = all(
      all_cols %in% names(TABLE_SPEC_TYPES)
    ),
    "column names must be unique" = !anyDuplicated(names(all_cols)),
    "descriptions must name every column" = setequal(names(descriptions), names(all_cols)),
    "key, partition and inheritance columns must be in columns" = all(
      c(key_cols, alternate_key_cols, partition_cols, inherit_by) %in% names(columns)
    ),
    "a run-scoped table needs an id_run column" = is.null(inherit_by) ||
      "id_run" %in% names(columns),
    "levels must name factor columns" = all(
      names(levels) %in% names(all_cols)[all_cols == "factor"]
    ),
    "allow_extra_cols must be TRUE or FALSE" = is_scalar_logical(allow_extra_cols)
  )

  map_cols <- names(columns)[columns == "map"]

  structure(
    list(
      name = name,
      columns = columns,
      optional_cols = optional_cols,
      descriptions = descriptions[names(all_cols)],
      key_cols = key_cols,
      alternate_key_cols = alternate_key_cols,
      partition_cols = partition_cols,
      map_cols = if (length(map_cols)) map_cols,
      inherit_by = inherit_by,
      levels = levels,
      allow_extra_cols = allow_extra_cols
    ),
    class = "table_spec"
  )
}

# spec type -> R predicate and prototype, and its name in the documentation; columns are
# cast by cast_dt_col(), which knows the same types
TABLE_SPEC_TYPES <- list(
  int = list(is = is.integer, proto = integer(0), doc = "int"),
  float = list(is = is.double, proto = double(0), doc = "float"),
  bool = list(is = is.logical, proto = logical(0), doc = "bool"),
  char = list(is = is.character, proto = character(0), doc = "string"),
  date = list(is = \(x) inherits(x, "Date"), proto = as.Date(character(0)), doc = "date"),
  factor = list(is = is.factor, proto = factor(), doc = "factor"),
  list = list(is = is.list, proto = list(), doc = "list"),
  map = list(is = is.list, proto = list(), doc = "map")
)

#' @describeIn table_spec Coerce `x` to the table that `spec` declares: build an empty table
#' if `x` is missing, cast every column, add the missing ones (filled with `NA`, or `NULL`
#' for list columns), set the canonical column order and finish with [as_ducklake_db_t()],
#' which validates the result.
#' @param x An object accepted by [data.table::setDT()], or missing for an empty table.
#' @param spec A `table_spec` object.
#' @param derive Optional function of the cast table, which computes derived columns by
#'   reference before the table is validated.
as_table_t <- function(x, spec, derive = NULL) {
  if (missing(x)) {
    x <- lapply(spec[["columns"]], \(type) TABLE_SPEC_TYPES[[type]][["proto"]])
    for (col in intersect(names(spec[["levels"]]), names(spec[["columns"]]))) {
      x[[col]] <- factor(character(0), levels = spec[["levels"]][[col]])
    }
  }
  data.table::setDT(x)

  for (col in names(spec[["columns"]])) {
    cast_spec_col(x, col, spec[["columns"]][[col]], spec[["levels"]][[col]])
  }
  for (col in intersect(names(spec[["optional_cols"]]), names(x))) {
    cast_spec_col(x, col, spec[["optional_cols"]][[col]], spec[["levels"]][[col]])
  }

  if (!is.null(derive)) {
    derive(x)
  }

  data.table::setcolorder(
    x,
    intersect(names(c(spec[["columns"]], spec[["optional_cols"]])), names(x))
  )
  data.table::setattr(x, "inherit_by", spec[["inherit_by"]])

  as_ducklake_db_t(
    x,
    class_name = spec[["name"]],
    key_cols = spec[["key_cols"]],
    alternate_key_cols = spec[["alternate_key_cols"]],
    map_cols = spec[["map_cols"]],
    partition_cols = spec[["partition_cols"]]
  )
}

# Cast one column by reference; a missing column is added, filled with NA (NULL for lists)
cast_spec_col <- function(x, col, type, levels = NULL) {
  if (type %in% c("list", "map")) {
    if (!col %in% names(x)) {
      # wrapped once more, or a single NULL would remove the column
      data.table::set(x, j = col, value = list(rep(list(NULL), nrow(x))))
    }
    return(invisible(x))
  }
  if (!col %in% names(x)) {
    value <- if (type == "factor") {
      factor(rep(NA_character_, nrow(x)), levels = levels)
    } else {
      rep(TABLE_SPEC_TYPES[[type]][["proto"]][NA], nrow(x))
    }
    data.table::set(x, j = col, value = value)
    return(invisible(x))
  }
  cast_dt_col(x, col, type, levels = levels)
}

#' @describeIn table_spec Check that `x` has the columns `spec` declares, with their types;
#' for use at the start of the `validate()` method of a table.
validate_table_spec <- function(x, spec) {
  for (col in names(spec[["columns"]])) {
    if (!col %in% names(x)) {
      stop(glue::glue("{spec[['name']]} needs a column {col}"), call. = FALSE)
    }
  }
  all_cols <- c(spec[["columns"]], spec[["optional_cols"]])
  for (col in intersect(names(all_cols), names(x))) {
    type <- all_cols[[col]]
    if (!TABLE_SPEC_TYPES[[type]][["is"]](x[[col]])) {
      stop(
        glue::glue(
          "{col} is not {TABLE_SPEC_R_TYPES[[type]]}, but {class(x[[col]])[1]}"
        ),
        call. = FALSE
      )
    }
  }
  if (!spec[["allow_extra_cols"]]) {
    stop_on_extra_cols(x, names(all_cols))
  }
  invisible(x)
}

# how a type is named in the messages of validate_table_spec()
TABLE_SPEC_R_TYPES <- c(
  int = "integer",
  float = "double",
  bool = "logical",
  char = "character",
  date = "a Date",
  factor = "a factor",
  list = "a list",
  map = "a list"
)

#' @describeIn table_spec The `@return` section of a table's help page, listing its columns;
#' for use as `#' @eval table_spec_rd(SPEC)`.
table_spec_rd <- function(spec) {
  optional <- names(spec[["descriptions"]]) %in% names(spec[["optional_cols"]])
  bullets <- paste0(
    "  - `",
    names(spec[["descriptions"]]),
    "`: ",
    ifelse(optional, "Optional. ", ""),
    gsub("\n", "\n    ", spec[["descriptions"]], fixed = TRUE)
  )
  paste(
    c(glue::glue('@return A data.table of class "{spec[["name"]]}" with columns:'), bullets),
    collapse = "\n"
  )
}

#' @describeIn table_spec The columns of a table as a Markdown table, for
#' `devdocs/design/database.md`: uniqueness, name, type and description.
table_spec_md <- function(spec) {
  all_cols <- c(spec[["columns"]], spec[["optional_cols"]])
  uniqueness <- ifelse(
    names(all_cols) %in% spec[["key_cols"]],
    "\\*",
    ifelse(names(all_cols) %in% spec[["alternate_key_cols"]], "alternate", "")
  )
  type <- vapply(all_cols, \(t) TABLE_SPEC_TYPES[[t]][["doc"]], character(1))
  type <- ifelse(
    names(all_cols) %in% names(spec[["optional_cols"]]),
    paste0(type, ", optional"),
    type
  )
  description <-
    spec[["descriptions"]] |>
    gsub("\n", " ", x = _, fixed = TRUE) |>
    # roxygen links to backticked code
    gsub("\\[`?([A-Za-z0-9_.:]+(\\(\\))?)`?\\]", "`\\1`", x = _) |>
    # roxygen escapes of square brackets
    gsub("\\\\([][])", "\\1", x = _) |>
    gsub("|", "\\|", x = _, fixed = TRUE)

  rows <- list(
    c("uniqueness", "colname", "type", "description"),
    cbind(uniqueness, names(all_cols), type, description)
  )
  cells <- rbind(rows[[1]], rows[[2]])
  widths <- apply(cells, 2, \(col) max(nchar(col, type = "width")))
  format_row <- function(row) {
    padded <- paste0(row, strrep(" ", widths - nchar(row, type = "width")))
    paste0("| ", paste(padded, collapse = " | "), " |")
  }
  c(
    format_row(cells[1, ]),
    paste0("| ", paste(strrep("-", widths), collapse = " | "), " |"),
    apply(cells[-1, , drop = FALSE], 1, format_row)
  ) |>
    unname()
}
