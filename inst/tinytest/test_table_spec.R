schema <- evoland.plus:::EVOLAND_SCHEMA
as_table_t <- evoland.plus:::as_table_t

# the registry holds one spec per exported table constructor, under the table's name
constructors <- setdiff(
  grep("^as_.*_t$", getNamespaceExports("evoland.plus"), value = TRUE),
  "as_ducklake_db_t"
)
expect_equal(sort(names(schema)), sort(sub("^as_", "", constructors)))
expect_equal(names(schema), unname(vapply(schema, \(spec) spec[["name"]], character(1))))

# every spec builds an empty table with the columns, types and attributes the constructor
# gives; runs_t alone needs a row, the base run
for (spec in schema) {
  from_spec <- if (spec[["name"]] == "runs_t") {
    as_table_t(list(id_run = 0L), spec)
  } else {
    as_table_t(spec = spec)
  }
  from_constructor <- get(paste0("as_", spec[["name"]]))()
  expect_equal(names(from_spec), names(spec[["columns"]]), info = spec[["name"]])
  expect_equal(lapply(from_spec, class), lapply(from_constructor, class), info = spec[["name"]])
  expect_equal(lapply(from_spec, levels), lapply(from_constructor, levels), info = spec[["name"]])
  expect_equal(
    attributes(from_spec)[c("class", "key_cols", "alternate_key_cols", "map_cols")],
    attributes(from_constructor)[c("class", "key_cols", "alternate_key_cols", "map_cols")],
    info = spec[["name"]]
  )
  expect_equal(attr(from_spec, "partition_cols"), spec[["partition_cols"]], info = spec[["name"]])
  expect_equal(attr(from_spec, "inherit_by"), spec[["inherit_by"]], info = spec[["name"]])
}

# inherit_by declares what get_evoland_db_read_expr() used to guess from the column names:
# every table with an id_run but runs_t is run-scoped, sliced by the id columns among
# id_period, id_period_post, id_trans and id_pred. For trans_preds_t that is each pair of
# id_pred and id_trans, so that a run deselects an inherited predictor by holding the pair.
for (spec in schema) {
  guessed <- if ("id_run" %in% names(spec[["columns"]]) && spec[["name"]] != "runs_t") {
    intersect(names(spec[["columns"]]), c("id_period", "id_period_post", "id_trans", "id_pred"))
  }
  expect_true(setequal(spec[["inherit_by"]], guessed), info = spec[["name"]])
  expect_equal(is.null(spec[["inherit_by"]]), is.null(guessed), info = spec[["name"]])
}
expect_equal(schema[["trans_preds_t"]][["inherit_by"]], c("id_pred", "id_trans"))

# a missing column is added, NA or, for a list column, NULL; a single fitted model without
# cross-validation can therefore be read back (ethzplus/evoland.plus#152)
one_model <- as_trans_models_t(data.table::data.table(
  id_run = 0L,
  id_trans = 1L,
  learner_id = "classif.rpart",
  learner_params = list(list(cp = 0.01)),
  crossval_score = list(list(no.crossval = 1))
))
expect_equal(names(one_model), names(schema[["trans_models_t"]][["columns"]]))
expect_equal(one_model[["crossval_predictions"]], list(NULL))
expect_equal(as_lulc_meta_t(list(id_lulc = 1L, name = "forest"))[["description"]], NA_character_)

# columns are cast, ad-hoc columns pass through after the declared ones, and tables that
# declare a closed set of columns reject the others
alloc <- as_alloc_params_t(list(id_trans = 1, id_run = 0, note = "manual"))
expect_identical(alloc[["id_trans"]], 1L)
expect_equal(names(alloc), c(names(schema[["alloc_params_t"]][["columns"]]), "note"))
expect_error(
  as_lulc_data_t(list(id_run = 0L, id_coord = 1L, id_lulc = 1L, id_period = 1L, note = "x")),
  "lulc_data_t holds only id_run, id_coord, id_lulc, id_period; drop note"
)

# validate() checks the declared types
runs <- as_runs_t()
runs[, rng_seed := 1.5]
expect_error(validate(runs), "rng_seed is not integer, but numeric")

# optional columns are cast when present, and not added when missing
expect_false("distance_class" %in% names(as_neighbors_t()))
neighbors <- as_neighbors_t(list(
  id_coord_origin = 1L,
  id_coord_neighbor = 2L,
  distance = 1L,
  distance_class = "[0,100)"
))
expect_inherits(neighbors[["distance_class"]], "factor")
expect_identical(neighbors[["distance"]], 1)

# a spec is checked for consistency when it is declared
expect_error(
  evoland.plus:::table_spec(
    "bad_t",
    columns = c(id = "int"),
    descriptions = c(id = "ID"),
    key_cols = "missing"
  ),
  "must be in columns"
)
expect_error(
  evoland.plus:::table_spec("bad_t", columns = c(id = "integer"), descriptions = c(id = "ID")),
  "column types must be one of"
)
expect_error(
  evoland.plus:::table_spec(
    "bad_t",
    columns = c(id = "int"),
    descriptions = c(id = "ID"),
    inherit_by = "id"
  ),
  "needs an id_run column"
)

# the column lists for the help pages and for devdocs/design/database.md
expect_equal(
  evoland.plus:::table_spec_rd(schema[["reporting_t"]]),
  paste(
    '@return A data.table of class "reporting_t" with columns:',
    "  - `key`: Unique character key",
    "  - `value`: Value as character",
    sep = "\n"
  )
)
expect_equal(
  evoland.plus:::table_spec_md(schema[["reporting_t"]]),
  c(
    "| uniqueness | colname | type   | description          |",
    "| ---------- | ------- | ------ | -------------------- |",
    "| \\*         | key     | string | Unique character key |",
    "|            | value   | string | Value as character   |"
  )
)
