LULC_META_T <- table_spec(
  "lulc_meta_t",
  columns = c(
    id_lulc = "int",
    name = "char",
    pretty_name = "char",
    description = "char",
    src_classes = "list"
  ),
  descriptions = c(
    id_lulc = "Unique ID for each land use class",
    name = "Concise name for use in code and queries, e.g. `forest_dense`",
    pretty_name = "Long name for plots and output, e.g. _Dense/Old Growth Forest_",
    description = "Long description / operationalisation",
    src_classes = "List of the source class IDs that are aggregated into this class"
  ),
  key_cols = "name",
  alternate_key_cols = "id_lulc"
)

#' Create LULC Metadata Table
#'
#' Create a `lulc_meta_t` table or rows thereof that meet the schema requirements. These
#' metadata describe the characteristics of land use classes.
#'
#' @name lulc_meta_t
#' @include table_spec.R
#'
#' @param x List-like coercible to `lulc_meta_t` via [data.table::as.data.table()]
#' @eval table_spec_rd(LULC_META_T)
#' @export
as_lulc_meta_t <- function(x) {
  as_table_t(x, LULC_META_T)
}

#' @describeIn lulc_meta_t Creates a lulc_meta_t table from land use class specifications
#' @param lulc_class_spec A named list of land use class specifications, see examples
#' @examples create_lulc_meta_t(list(
#'   forest = list(
#'     pretty_name = "Dense Forest",
#'     description = "Mature, dense forest cover",
#'     src_classes = c(1L, 2L)
#'   ),
#'   arable = list(
#'     pretty_name = "Arable land",
#'     description = "Non-permanent cropland",
#'     src_classes = 3:6
#'   )
#' ))
#' @export
create_lulc_meta_t <- function(lulc_class_spec) {
  if (!is.list(lulc_class_spec) || length(lulc_class_spec) == 0) {
    stop("lulc_classes must be a non-empty list")
  }

  # Extract class names
  class_names <- names(lulc_class_spec)
  if (is.null(class_names) || any(class_names == "")) {
    stop("All lulc_classes must have names")
  }

  # Create the data.table
  x <- data.table::data.table(
    id_lulc = seq_along(class_names),
    name = class_names,
    pretty_name = vapply(
      lulc_class_spec,
      function(cls) cls[["pretty_name"]] %||% cls[["name"]] %||% NA_character_,
      character(1)
    ),
    description = vapply(
      lulc_class_spec,
      function(cls) cls[["description"]] %||% NA_character_,
      character(1)
    ),
    src_classes = lapply(
      pluck_wildcard(lulc_class_spec, NA, "src_classes"),
      \(y) if (is.null(y)) integer(0) else as.integer(y)
    )
  )

  data.table::setkey(x, "id_lulc")

  as_lulc_meta_t(x)
}

#' @export
validate.lulc_meta_t <- function(x, ...) {
  NextMethod()
  validate_table_spec(x, LULC_META_T)

  stopifnot("name cannot be empty" = !any(x[["name"]] == ""))

  return(x)
}

#' @export
#' @describeIn lulc_meta_t Print an lulc_meta_t object, passing params to data.table print
#' @param nrow see [data.table::print.data.table]
#' @param ... passed to [data.table::print.data.table]
print.lulc_meta_t <- function(x, nrow = 10, ...) {
  if (nrow(x) > 0) {
    cat(glue::glue(
      "LULC Metadata Table\n",
      "Number of classes: {nrow(x)}\n\n"
    ))
  } else {
    cat("LULC Metadata Table (empty)\n")
  }
  NextMethod(nrow = nrow, ...)
  invisible(x)
}
