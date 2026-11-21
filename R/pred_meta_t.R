# leaving out POSIXct / Date for now; can be operationalized as int/float
PRED_DATA_TYPES <- c("int", "float", "bool", "factor", "ordered")

PRED_META_T <- table_spec(
  "pred_meta_t",
  columns = c(
    id_pred = "int",
    name = "char",
    pretty_name = "char",
    description = "char",
    orig_format = "char",
    sources = "list",
    unit = "char",
    factor_levels = "list",
    data_type = "factor",
    fill_value = "char"
  ),
  descriptions = c(
    id_pred = "Unique ID for each predictor",
    name = "Concise name for use in code and queries, e.g. `distance_road`",
    pretty_name = "Name for plots and output, e.g. _Distance To Closest Road_",
    description = "Long description / operationalisation",
    orig_format = "Original format, e.g. \"100m raster, daily\"",
    sources = paste(
      "Where the data were fetched: list column of data.frames with columns `url` and",
      "`md5sum`"
    ),
    unit = paste(
      "SI units for physical properties, or more complex descriptors like \"bed nights/year\"",
      "as a proxy for touristic activity; `NA` for categorical values"
    ),
    factor_levels = paste(
      "List of character vectors, the ordered level labels when `data_type` is",
      "\"factor\" or \"ordered\""
    ),
    data_type = paste(
      "Factor with levels \"int\", \"float\", \"bool\", \"factor\", \"ordered\": the type",
      "the stored values are cast back to"
    ),
    fill_value = paste(
      "Value of the [coords_t] coordinate points that are not explicitly stored in",
      "[pred_data_t]"
    )
  ),
  key_cols = "name",
  alternate_key_cols = "id_pred",
  levels = list(data_type = PRED_DATA_TYPES)
)

#' Create Predictor Metadata Table
#'
#' Construct and validate `pred_meta_t` objects, which are used to store predictor
#' metadata.
#'
#' @name pred_meta_t
#' @include ducklake_db_utils.R
#'
#' @param x An object that is accepted by [data.table::setDT()]
#'
#' @eval table_spec_rd(PRED_META_T)
#' @export
as_pred_meta_t <- function(x) {
  as_table_t(x, PRED_META_T, derive = function(x) {
    x[,
      sources := lapply(
        sources,
        function(src) {
          # coerce to data.table with exactly url & md5sum
          if (
            inherits(src, "data.frame") &&
              all(utils::hasName(src, c("url", "md5sum")))
          ) {
            return(src[, c("url", "md5sum")])
          }
          src_dt <- data.table::rbindlist(src, use.names = TRUE)
          if (length(src_dt) == 0L) {
            # length 0 is a null data frame, e.g. if src is NULL or list()
            src_dt <- data.table::data.table(url = character(), md5sum = character())
          }
          src_dt[, .(url, md5sum)]
        }
      )
    ]
  })
}


#' @describeIn pred_meta_t Creates a `pred_meta_t` table from predictor
#' specifications. Attributes incremental IDs to predictors, starting from
#' `starting_id`. See examples for schema of `pred_spec`.
#' @param pred_spec A list of predictor specifications, schema: see examples
#' @param starting_id Integer, starting ID for predictors; default 1L
#' @examples create_pred_meta_t(list(
#'    noise = list(
#'      unit = "dBa",
#'      pretty_name = "Maximum noise exposure",
#'      orig_format = "10m*10m raster",
#'      description = "daytime & nighttime road & rail noise exposure",
#'      sources = list(
#'        list(
#'          url = "https://example.com/noisedata1.tif",
#'          md5sum = "a4b9f1c04ee63824f18852bfd1eecbdd"
#'        ),
#'        list(
#'          url = "https://example.com/noisedata2.tif",
#'          md5sum = "4b782128495b5af8467e2259bd57def2"
#'        )
#'      ),
#'      data_type = "float"
#'    ),
#'    distance_to_lake = list(
#'      unit = "m",
#'      pretty_name = "Distance to closest lake",
#'      orig_format = "vector",
#'      description = "Derived from swissTLM3D",
#'      sources = list(list(
#'        url = "https://example.com/dist_to_lake.tif",
#'        md5sum = "ecb3bcfbf6316c6e7542e20de24f61b7"
#'      )),
#'      data_type = "float"
#'    )
#'  ))
#' @export
create_pred_meta_t <- function(pred_spec, starting_id = 1L) {
  # Extract predictor names
  pred_names <- names(pred_spec)
  if (is.null(pred_names) || any(pred_names == "")) {
    stop("All pred_data_spec entries must have names")
  }

  x <- data.table::data.table(
    id_pred = seq(from = starting_id, length.out = length(pred_spec)),
    name = pred_names,
    # path is pred_spec > pred_name > leaf_name
    # we pluck each element, then replace potential null using %||%
    pretty_name = unlist(
      mapply(
        function(x, y) x %||% y,
        pluck_wildcard(pred_spec, NA, "pretty_name"),
        pred_names,
        SIMPLIFY = FALSE
      )
    ),
    description = unlist(
      lapply(pluck_wildcard(pred_spec, NA, "description"), function(x) x %||% NA_character_)
    ),
    orig_format = unlist(
      lapply(pluck_wildcard(pred_spec, NA, "orig_format"), function(x) x %||% NA_character_)
    ),
    sources = lapply(
      pluck_wildcard(pred_spec, NA, "sources"),
      function(src) {
        src_dt <- data.table::rbindlist(src, use.names = TRUE)
        if (length(src_dt) == 0L) {
          src_dt <- data.table::data.table(url = character(), md5sum = character())
        }
        src_dt[, .(url, md5sum)]
      }
    ),
    unit = unlist(
      lapply(pluck_wildcard(pred_spec, NA, "unit"), function(x) x %||% NA_character_)
    ),
    data_type = {
      lapply(pluck_wildcard(pred_spec, NA, "data_type"), function(x) {
        x %||% NA_character_
      }) |>
        unlist() |>
        factor(levels = PRED_DATA_TYPES)
    },
    fill_value = unlist(
      lapply(pluck_wildcard(pred_spec, NA, "fill_value"), function(x) x %||% NA_character_)
    ),
    factor_levels = lapply(
      pluck_wildcard(pred_spec, NA, "factor_levels"),
      \(y) if (is.null(y)) character(0) else as.character(y)
    )
  )

  as_pred_meta_t(x)
}

#' @export
validate.pred_meta_t <- function(x, ...) {
  NextMethod()

  sources_dt <- unique(data.table::rbindlist(x[["sources"]]))

  stopifnot(
    "data_type must be set" = !any(is.na(x[["data_type"]])),
    "data_type can only be one of 'int', 'float', 'bool', 'factor', 'ordered'" = setequal(
      levels(x[["data_type"]]),
      PRED_DATA_TYPES
    ),
    "name cannot be empty" = !any(x[["name"]] == ""),
    "single URL with multiple checksums present" = !anyDuplicated(sources_dt[["url"]]),
    "sources must hold exactly a url and an md5sum field" = length(sources_dt) == 0L ||
      setequal(names(sources_dt), c("url", "md5sum"))
  )

  return(x)
}

#' @export
#' @describeIn pred_meta_t Print an `pred_meta_t` object, passing params to data.table print
#' @param ... passed to [data.table::print.data.table()]
print.pred_meta_t <- function(x, ...) {
  if (nrow(x) > 0) {
    cat(glue::glue(
      "Predictor Metadata Table\n",
      "Number of predictors: {nrow(x)}\n\n"
    ))
  } else {
    cat("Predictor Metadata Table (empty)\n")
    return(invisible(x))
  }

  NextMethod(trunc.cols = TRUE, ...)

  invisible(x)
}
