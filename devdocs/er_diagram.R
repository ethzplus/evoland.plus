# Regenerates the ER diagram and the column tables in devdocs/design/database.md from the
# table specifications.
#
#   Rscript devdocs/er_diagram.R
#
# Every table is declared once as a `table_spec()` next to its `as_*_t()` constructor, and
# the registry `EVOLAND_SCHEMA` lists them all; the diagram and the column tables are derived
# from it. A column is a foreign key if it is named after another table's id: `id_<entity>`
# for the table `<entity>_t`, `<entity>s_t` or `<entity>_meta_t`, optionally with a role
# prefix or suffix (`parent_id_run`, `id_lulc_anterior`, `id_period_post`). CI fails if the
# committed document differs from the generated one, see .github/workflows/R-CMD-check.yaml.

if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(quiet = TRUE)
}
schema <- utils::getFromNamespace("EVOLAND_SCHEMA", "evoland.plus")
schema <- schema[sort(names(schema))]
column_type <- \(type) utils::getFromNamespace("TABLE_SPEC_TYPES", "evoland.plus")[[type]][["doc"]]
table_spec_md <- utils::getFromNamespace("table_spec_md", "evoland.plus")

target_file <- file.path("devdocs", "design", "database.md")
generated_note <- "run `Rscript devdocs/er_diagram.R`, do not edit"

# id column -> table it identifies, from the naming convention described above
id_owner <- list()
for (spec in schema) {
  base <- sub("_t$", "", spec[["name"]])
  entities <- unique(c(base, sub("_meta$", "", base), sub("s$", "", base)))
  key_cols <- c(spec[["key_cols"]], spec[["alternate_key_cols"]])
  for (id_col in intersect(paste0("id_", entities), key_cols)) {
    id_owner[[id_col]] <- spec[["name"]]
  }
}

referenced_table <- function(column_name, type) {
  if (type != "int") {
    return(NA_character_)
  }
  for (id_col in names(id_owner)) {
    pattern <- paste0("^([a-z]+_)?", id_col, "(_[a-z]+)?$")
    if (grepl(pattern, column_name)) {
      return(id_owner[[id_col]])
    }
  }
  NA_character_
}

entity_lines <- character()
relation_lines <- character()
for (spec in schema) {
  table_name <- spec[["name"]]
  entity_lines <- c(entity_lines, paste0("    ", table_name, " {"))
  for (column_name in names(spec[["columns"]])) {
    type <- spec[["columns"]][[column_name]]
    parent <- referenced_table(column_name, type)
    is_foreign <- !is.na(parent) && !(parent == table_name && column_name %in% names(id_owner))
    markers <- c(
      if (column_name %in% spec[["key_cols"]]) "PK",
      if (column_name %in% spec[["alternate_key_cols"]]) "UK",
      if (is_foreign) "FK"
    )
    entity_lines <- c(
      entity_lines,
      paste("       ", column_type(type), column_name, paste(markers, collapse = ", ")) |>
        trimws(which = "right")
    )
    if (is_foreign) {
      cardinality <- if (startsWith(column_name, "parent_")) "}o--o|" else "}o--||"
      relation_lines <- c(
        relation_lines,
        glue::glue('    {table_name} {cardinality} {parent} : "{column_name}"')
      )
    }
  }
  entity_lines <- c(entity_lines, "    }")
}

# Replace the lines between a begin and an end marker, each present exactly once
replace_between <- function(document, begin_marker, end_marker, lines) {
  begin <- which(document == begin_marker)
  end <- which(document == end_marker)
  if (!(length(begin) == 1L && length(end) == 1L && begin < end)) {
    stop("database.md must contain these markers exactly once, in order: ", begin_marker)
  }
  c(document[seq_len(begin)], lines, document[-seq_len(end - 1L)])
}

document <- replace_between(
  readLines(target_file),
  paste0("<!-- BEGIN GENERATED ER DIAGRAM: ", generated_note, " -->"),
  "<!-- END GENERATED ER DIAGRAM -->",
  c("```mermaid", "erDiagram", entity_lines, "", relation_lines, "```")
)
for (spec in schema) {
  document <- replace_between(
    document,
    glue::glue("<!-- BEGIN GENERATED COLUMNS {spec[['name']]}: {generated_note} -->"),
    glue::glue("<!-- END GENERATED COLUMNS {spec[['name']]} -->"),
    table_spec_md(spec)
  )
}
writeLines(document, target_file)
