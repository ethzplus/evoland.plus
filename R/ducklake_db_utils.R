#' DuckLake database utility functions
#'
#' @description
#' A set of utility functions for working with [ducklake_db].
#'
#' @name ducklake_db_utils
NULL

#' @describeIn ducklake_db_utils Coerce to ducklake_db_t subclass. It coerces an object to a
#' data.table and adds attributes that [ducklake_db] relies on for database-like
#' operations. See the paramaters for details.
#' @param class_name Optional class name to prepend before "ducklake_db_t". Used
#'        to create more specific subclasses.
#' @param key_cols Optional character vector. Used for data.table's [data.table::setkey]
#'        and used by [ducklake_db] to determine which columns to use for upsert
#'        operations. Together with `alternate_key_cols`, these are used to define
#'        uniqueness constraints.
#' @param alternate_key_cols Optional character vector. These are used in
#'        conjunction with `key_cols` to define uniqueness constraints, but are
#'        not used for upsert operations. Example: numeric primary key whose
#'        distinctness is identical to a character primary key column, but
#'        upsert operations should be based on the character key, leaving the
#'        character<->numeric mapping intact.
#' @param map_cols Optional character vector. Used to coerce columns to DuckDB's
#'        MAP type, used to store named unnested lists from R.
#' @param partition_cols Optional character vector. Used to specify columns for
#'        hive style partitioning of the data files.
#' @export
as_ducklake_db_t <- function(
  x,
  class_name = character(),
  key_cols = NULL,
  alternate_key_cols = NULL,
  map_cols = NULL,
  partition_cols = NULL
) {
  data.table::setDT(x)

  # subclassses for more specific validation methods
  distinct_classes <- unique(c(class_name, "ducklake_db_t", class(x)))
  data.table::setattr(x, "class", distinct_classes)

  # uniqueness, matching, distinctness can be derived from key cols
  data.table::setkeyv(x, key_cols)
  data.table::setattr(x, "key_cols", key_cols)
  data.table::setattr(x, "alternate_key_cols", alternate_key_cols)

  # list cols to coerce to MAP
  data.table::setattr(x, "map_cols", map_cols)

  # partition_cols for hive style partitioning
  data.table::setattr(x, "partition_cols", partition_cols)

  validate(x)
}

#' @export
validate.ducklake_db_t <- function(x, ...) {
  attrs_to_check <- setdiff(
    names(attributes(x)),
    c(
      # data.table internal attrs, not committed to the table metadata
      "names",
      "row.names",
      "class",
      ".internal.selfref",
      "index"
    )
  )

  stopifnot(
    "needs to be a data.table" = inherits(x, "data.table"),
    "all attributes need to be atomic" = all(
      vapply(
        attrs_to_check,
        function(y) is.atomic(attr(x, y)),
        logical(1)
      )
    ),
    "key_cols, alternate_key_cols, map_cols, and partition_cols must be present" = all(
      c(
        attr(x, "key_cols"),
        attr(x, "alternate_key_cols"),
        attr(x, "map_cols"),
        attr(x, "partition_cols")
      ) %in%
        names(x)
    )
  )

  for (col in names(x)) {
    if (
      col %in%
        c(
          attr(x, "key_cols"),
          attr(x, "alternate_key_cols"),
          attr(x, "partition_cols")
        )
    ) {
      if (anyNA(x[[col]])) {
        stop(glue::glue(
          "Column '{col}' cannot contain NA values"
        ))
      }
    } else if (col %in% attr(x, "map_cols")) {
      for (val in x[[col]]) {
        is_valid_mapcol <- function(v) {
          # allow NULLs
          is.null(v) ||
            # allow empty lists
            (is.list(v) && length(v) == 0) ||
            # allow named lists with atomic values
            (is.list(v) && !is.null(names(v)) && all(vapply(v, is.atomic, logical(1))))
        }

        # negation makes the above condition easier to read
        if (Negate(is_valid_mapcol)(val)) {
          stop(glue::glue(
            "Column '{col}' specified as map_cols must be a list of ",
            "named lists with atomic values"
          ))
        }
      }
    } else if (is.list(x[[col]])) {
      distinct_classes <-
        lapply(x[[col]], \(y) if (is.null(y)) NULL else class(y)) |>
        unique() |>
        Filter(Negate(is.null), x = _) # drop NULLs
      if (length(distinct_classes) > 1) {
        stop(glue::glue(
          "All elements of list column '{col}' must have the same class"
        ))
      }
    }
  }

  if (
    length(attr(x, "key_cols")) &&
      anyDuplicated(x, by = attr(x, "key_cols"))
  ) {
    stop(glue::glue(
      "Duplicates found in key_cols\n",
      "  key_cols: {toString(attr(x, 'key_cols'))}"
    ))
  }
  if (
    length(attr(x, "alternate_key_cols")) &&
      anyDuplicated(x, by = attr(x, "alternate_key_cols"))
  ) {
    stop(glue::glue(
      "Duplicates found in alternate_key_cols\n",
      "  alternate_key_cols: {toString(attr(x, 'alternate_key_cols'))}"
    ))
  }

  x
}


#' @describeIn ducklake_db_utils Read one of the package's `inst/*.sql` files, for
#' passing to `$get_query()` -- which interpolates it, so the placeholders in the
#' file are filled from the named arguments given alongside it.
#' @param file Name of a SQL file in the installed package.
read_sql <- function(file) {
  paste(readLines(system.file(file, package = "evoland.plus")), collapse = "\n")
}

#' @describeIn ducklake_db_utils Serialize a named list of atomic vectors into the
#' single string that [ducklake_db] stores as a comment on the catalog table. Each
#' entry becomes one line, `key: "value1", "value2"`; a `"` in a value is written `\"`.
#' @param metadata Named list of atomic vectors; keys match `\w+`, values hold no newline.
#' @keywords internal
serialize_metadata <- function(metadata) {
  if (length(metadata) == 0L) {
    return("")
  }
  stopifnot(
    "metadata keys must match \\w+" = all(grepl("^\\w+$", names(metadata))),
    "metadata values must not contain a newline" = !any(grepl("\n", unlist(metadata)))
  )

  values <- vapply(
    metadata,
    function(value) paste0('"', gsub('"', '\\"', value, fixed = TRUE), '"', collapse = ", "),
    character(1)
  )

  paste0(names(metadata), ": ", values, collapse = "\n")
}

#' @describeIn ducklake_db_utils Recover the named list written by
#' [serialize_metadata()]. Values are type-converted, so that a numeric
#' attribute survives the round-trip as a number rather than as a string.
#' @param comment Character scalar, or NA for a table without metadata.
#' @keywords internal
deserialize_metadata <- function(comment) {
  if (length(comment) == 0L || is.na(comment) || !nzchar(comment)) {
    return(list())
  }

  entries <- strsplit(comment, "\n", fixed = TRUE)[[1]]
  keys <- sub(": .*$", "", entries)
  values <- sub("^[^:]*: ", "", entries)

  parsed <- lapply(values, function(value) {
    strsplit(value, '(?<!\\\\)", "', perl = TRUE)[[1]] |>
      gsub('^"|"$', "", x = _) |>
      gsub('\\"', '"', x = _, fixed = TRUE) |>
      utils::type.convert(as.is = TRUE)
  })

  stats::setNames(parsed, keys)
}

#' @describeIn ducklake_db_utils Convert list columns by applying `fn`
#' @param cols The list columns to convert.
#' @param fn The function to apply to each element of the list columns
#' @keywords internal
convert_list_cols <- function(x, cols, fn) {
  for (col in cols) {
    if (!col %in% names(x)) {
      warning(glue::glue("Column '{col}' not found in data; skipping conversion"))
      next
    }
    x[[col]] <- lapply(x[[col]], fn)
  }
  x
}

#' @describeIn ducklake_db_utils Convert a named list to a data.frame with "key" and
#' "value" columns.
#' @keywords internal
list_to_kv_df <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(data.frame(
      key = character(0),
      value = character(0),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    key = names(x),
    value = as.character(unlist(x)),
    stringsAsFactors = FALSE
  )
}

#' @describeIn ducklake_db_utils Convert a data.frame with "key" and "value" columns to a
#' named list.
#' @keywords internal
kv_df_to_list <- function(x) {
  if (is.null(x) || nrow(x) == 0) {
    return(NULL)
  }

  out <- list()

  for (row in seq_len(nrow(x))) {
    key <- x$key[row]
    val <- utils::type.convert(x$value[row], as.is = TRUE)
    out[[key]] <- val
  }

  out
}

#' @describeIn ducklake_db_utils Create a binding function for a table, which can be used
#' to fetch or commit data to that table. The active binding either returns the table
#' (missing argument), or upserts to it (assignment operation)
#' @param table_name The name of the table to bind to.
#' @param mode The mode of the binding, which determines the behavior when
#' committing data. Options are: "write_once" (default, only allows writing if
#' table doesn't exist), "upsert", "append", and "overwrite".
create_table_binding <- function(
  table_name,
  mode = c("write_once", "upsert", "append", "overwrite")
) {
  force(table_name)
  mode <- match.arg(mode)
  # Use bquote to "bake in" the table_name value directly into the function body,
  # so it doesn't depend on the enclosing environment (which R6 will replace).
  fn <- eval(bquote(
    function(x) {
      tbl <- .(table_name)
      md <- .(mode)

      if (missing(x)) {
        # a binding is an interactive view, so its read is cut off after binding_timeout;
        # DuckDB's limit holds for the whole connection, so it is set for this read only
        if (!is.null(self$binding_timeout)) {
          # at least 1 ms, as 0 means no limit
          self$execute(
            "set max_execution_time = {max(1, round(self$binding_timeout * 1000))}"
          )
          on.exit(self$execute("reset max_execution_time"), add = TRUE)
        }
        fetched <- tryCatch(self$fetch(table_name = tbl), error = function(e) {
          if (!grepl("exceeded maximum execution time", conditionMessage(e), fixed = TRUE)) {
            stop(e)
          }
          stop(
            "Reading `",
            tbl,
            "` through its table binding took longer than ",
            self$binding_timeout,
            " seconds and was interrupted. Read a large table with ",
            "`$fetch()`, narrowed with `where` or `id_run`, or raise the database's ",
            "`binding_timeout`.",
            call. = FALSE
          )
        })
        as_fn <- paste0("as_", tbl) |> get()
        return(as_fn(fetched))
      }

      stopifnot(inherits(x, tbl))

      if (md == "overwrite" && interactive() && tbl %in% self$list_tables()) {
        confirm_overwrite <- utils::menu(
          c("Yes", "No"),
          title = glue::glue(
            "Table '{tbl}' already exists. Are you sure you want to overwrite it?"
          )
        )
        if (confirm_overwrite != 1) {
          return(NULL)
        }
      }

      if (md == "write_once") {
        if (tbl %in% self$list_tables()) {
          warning(
            glue::glue(
              "Table '{tbl}' is write-once! Delete manually and write ",
              "anew ONLY if you know what you are doing (-> foreign keys)!"
            ),

            call. = FALSE
          )
          return(NULL)
        }
        md <- "overwrite"
      }

      self$commit(x, table_name = tbl, method = md)
    }
  ))
  fn
}

#' @describeIn ducklake_db_utils Helper function to bind a function as a method to an R6
#' generator. Simply passes the R6 method's arguments as-is to the underlying "pure"
#' function, passing the R6 reference to `self`, but not to `private`.
#' @keywords internal
#' @param fun The underlying function to bind as an R6 method, which must have a `self` argument
#' @param with_private Whether to also pass the R6 reference to `private` as an argument
#' @param with_super Whether to also pass the R6 reference to `super` as an argument
create_method_binding <- function(fun, with_private = FALSE, with_super = FALSE) {
  # Capture the original call (e.g., obj$method(arg1 = 1)); expands ... arguments
  cl <- match.call(definition = sys.function(-1), call = sys.call(-1))

  # Modify 'obj$method' to 'fun' instead
  cl[[1]] <- fun

  # Inject 'self' as named arg. Use get() to retrieve R6 instance's self from calling env
  cl[["self"]] <- get("self", envir = parent.frame())
  if (with_private) {
    cl[["private"]] <- get("private", envir = parent.frame())
  }
  if (with_super) {
    cl[["super"]] <- get("super", envir = parent.frame())
  }

  # Evaluate in original environment. Preserves lazy evaluation of arguments.
  eval(cl, envir = parent.frame(2))
}

#' @describeIn ducklake_db_utils Declare a table once, as data. Its constructor
#' (`as_table_t()`), its help page's column list (`table_spec_rd()`) and, for the tables of
#' [evoland_db], the ER diagram and the reads through the run lineage are derived from it.
#' @param name The table name, which is also the class of its objects.
#' @param columns Named character vector, column name to type, in the canonical order. Types
#'   are those of `cast_dt_col()` (`"int"`, `"float"`, `"bool"`, `"char"`, `"date"`,
#'   `"factor"`), `"list"`, and `"map"` for a list column stored as a DuckDB `MAP`.
#' @param descriptions Named character vector, the meaning of every column, for the
#'   documentation.
#' @param inherit_by `NULL` for a table that is not run-scoped; otherwise the columns
#'   besides `id_run` that form one slice of inheritance through the run lineage.
#' @param optional_cols Like `columns`, for columns that are cast when present but not added.
#' @param levels Named list, factor column to its levels.
#' @param allow_extra_cols If `FALSE`, columns not declared are rejected.
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
  stopifnot(
    "descriptions must name every column" = setequal(
      names(descriptions),
      names(c(columns, optional_cols))
    ),
    "keys, partitions and inherit_by must be columns" = all(
      c(key_cols, alternate_key_cols, partition_cols, inherit_by) %in% names(columns)
    )
  )
  map_cols <- names(columns)[columns == "map"]
  structure(
    list(
      name = name,
      columns = columns,
      optional_cols = optional_cols,
      descriptions = descriptions[names(c(columns, optional_cols))],
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

#' @describeIn ducklake_db_utils Coerce `x` (or, if missing, an empty table) to the table
#' `spec` declares: cast every column, add the missing ones as `NA` (`NULL` for list
#' columns), apply `derive`, set the column order and finish with [as_ducklake_db_t()].
#' @param spec A `table_spec`.
#' @param derive Optional function computing derived columns of the cast table by reference.
as_table_t <- function(x, spec, derive = NULL) {
  x <- if (missing(x)) data.table::data.table() else data.table::setDT(x)
  all_cols <- c(spec[["columns"]], spec[["optional_cols"]])
  for (col in names(all_cols)) {
    type <- all_cols[[col]]
    if (!col %in% names(x)) {
      if (col %in% names(spec[["optional_cols"]])) {
        next
      }
      # wrapped in a list, or a single NULL would remove the column
      value <- if (type %in% c("list", "map")) {
        list(rep(list(NULL), nrow(x)))
      } else {
        rep(NA, nrow(x))
      }
      data.table::set(x, j = col, value = value)
    }
    if (!type %in% c("list", "map")) {
      cast_dt_col(x, col, type, levels = spec[["levels"]][[col]])
    }
  }
  if (!is.null(derive)) {
    derive(x)
  }
  if (!spec[["allow_extra_cols"]]) {
    stop_on_extra_cols(x, names(all_cols))
  }
  data.table::setcolorder(x, intersect(names(all_cols), names(x)))
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

#' @describeIn ducklake_db_utils The `@return` section listing a table's columns, for
#' `#' @eval table_spec_rd(SPEC)`.
table_spec_rd <- function(spec) {
  optional <- names(spec[["descriptions"]]) %in% names(spec[["optional_cols"]])
  paste(
    c(
      paste0('@return A data.table of class "', spec[["name"]], '" with columns:'),
      paste0(
        "  - `",
        names(spec[["descriptions"]]),
        "`: ",
        ifelse(optional, "Optional. ", ""),
        spec[["descriptions"]]
      )
    ),
    collapse = "\n"
  )
}
