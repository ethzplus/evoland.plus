#' Utility functions and methods for `evoland_db`
#'
#' Factored out utility functions and methods specifically used in conjunction with
#' `evoland_db` and its subclasses.
#'
#' @name evoland_db_util
NULL

#' @describeIn evoland_db_util Gets the read expression for a given table, taking into
#' account the run lineage and distinctness columns. If the table has no `id_run` column
#' or is the `runs_t` table itself, falls back to the default read expression.
#' Otherwise, constructs a read expression that filters for `id_run` and its lineage,
#' using the distinctness columns to determine which slices of data to include.
#' @param self The `evoland_db` instance
#' @param super The superclass (parent class) environment.
#' @param table_name The name of the table to read from
#' @param id_run Integer, the run whose lineage is read; `NULL` reads every run as stored
#' @return A SQL expression string to read most specific data for slice
get_evoland_db_read_expr <- function(self, super, table_name, id_run) {
  base_read_expr <- super$get_read_expr(table_name)
  all_cols <- self$get_query(
    "select column_name from (describe {base_read_expr})",
    as_atomic = TRUE
  )

  if (
    is.null(id_run) || # no run to read
      !("id_run" %in% all_cols) || # no id_run in table
      table_name == "runs_t" # runs_t table itself
  ) {
    return(super$get_read_expr(table_name))
  }

  # Hardcode the set of columns to probe for distinctness; i.e. if any value is present
  # for a given distinctness tuple, the corresponding slice of data will be used. This
  # is a tradeoff: we could also choose to fall through on each coordinate point, but
  # that would blow up the query time.
  theoretical_distinctness_cols <- c(
    "id_run",
    "id_period",
    "id_period_post",
    "id_trans",
    "id_pred"
  )
  distinctness_cols <- intersect(all_cols, theoretical_distinctness_cols)
  inheritance_key_cols <- setdiff(distinctness_cols, "id_run")

  # .envir, or the statement would be interpolated in this helper's frame
  sql <- function(...) glue::glue_sql(..., .con = self$connection, .envir = parent.frame())
  run_lineage <- get_lineage(self$runs_t, as.integer(id_run))

  # Single run in lineage: just filter for id_run
  if (length(run_lineage) == 1L) {
    return(sql("(select * from {base_read_expr} where id_run = {id_run})"))
  }

  # map each id_run in lineage to its distance from the requested run; used to
  # find the minimum distance. e.g. if run_lineage is (3, 2, 0) we get
  # case id_run when 3 then 1 when 2 then 2 when 0 then 3 else 999999 end
  run_case <- sql(
    "case b.id_run {when_clauses} else 999999 end",
    when_clauses = glue::glue_sql_collapse(
      sql("when {run_lineage} then {seq_along(run_lineage)}"),
      sep = " "
    )
  )

  ctes <- list()

  # return one row per tuple of available data within lineage
  ctes[["data_present"]] <- sql(
    r"[
    select distinct
      {`distinctness_cols`*}
    from
      {base_read_expr}
    where
      id_run in ({run_lineage*})
    ]"
  )

  # Special case for id_period: need self-join on id_period AND id_period=0;
  # currently only relevant for pred_data_t, but let's generalize just in case
  if ("id_period" %in% inheritance_key_cols) {
    join_conditions <- lapply(
      inheritance_key_cols,
      function(col) {
        if (col == "id_period") {
          # special case for id_period: self-join data_present to allow for
          # id_period=0 fallback
          return(sql("(a.{`col`} = b.{`col`} or b.{`col`} = 0)"))
        }
        sql("a.{`col`} = b.{`col`}")
      }
    )
    # reduce data_present to one row per most specific id_run
    ctes[["best_run"]] <- sql(
      r"[
      select
        {group_cols},
        -- arg_min returns id_run for the single row where run_case is minimal
        arg_min(b.id_run, {run_case}) as id_run
      from
        data_present a,
        data_present b
      where
        {glue::glue_sql_collapse(join_conditions, sep = " and ")}
      group by
        {group_cols}
      ]",
      group_cols = glue::glue_sql_collapse(sql("b.{`inheritance_key_cols`}"), sep = ", ")
    )
  } else if (length(inheritance_key_cols) > 0L) {
    # general case: find minimum distance id_run for each tuple of inheritance key cols
    ctes[["best_run"]] <- sql(
      r"[
      select
        {group_cols},
        arg_min(b.id_run, {run_case}) as id_run
      from
        data_present b
      group by
        {group_cols}
      ]",
      group_cols = glue::glue_sql_collapse(sql("b.{`inheritance_key_cols`}"), sep = ", ")
    )
  } else {
    # no inheritance keys: the nearest run holding any data wins
    ctes[["best_run"]] <- sql(
      "select arg_min(b.id_run, {run_case}) as id_run from data_present b"
    )
  }

  # use a semi join to filter the table using best_run
  sql(
    r"[(
    with
      data_present as (
        {ctes[["data_present"]]}
      ),
      best_run as (
        {ctes[["best_run"]]}
      )
    from
      {base_read_expr} c
    semi join
      best_run b
      using ({`distinctness_cols`*})
    where
      c.id_run in ({run_lineage*})
    )]"
  )
}

#' @describeIn evoland_db_util Apply a worker function to each item, serially or on mirai
#' daemons. Without `compute`, the worker gets `parent_db` itself. With `compute`, each daemon
#' of that compute profile opens one connection to the same catalog and data files as
#' `parent_db` and reuses it for every item it is given. Workers must therefore not rely on
#' the connection's active run: an item carries the run it is for.
#'
#' Set the daemons up beforehand, e.g. `mirai::daemons(4, .compute = "evoland")`, and pass
#' `compute = "evoland"`; evoland.plus must be installed where the daemons run.
#'
#' @param worker_fun A function to apply to each item, with signature `worker_fun(item,
#' db, ...)`.
#' @param items A list of items to iterate over
#' @param parent_db An [evoland_db] instance: passed to the worker in the serial case, and
#' giving the catalog and data path the daemons connect to otherwise.
#' @param compute `NULL` (default) to run serially, or the name of a mirai compute profile
#' whose daemons are set, see [mirai::daemons()].
#' @param ... Additional arguments passed to `worker_fun`.
#'
#' @return A list of results
#' @export
run_parallel_evoland <- function(
  items,
  worker_fun,
  parent_db,
  compute = NULL,
  ...
) {
  if (is.null(compute)) {
    return(lapply(X = items, FUN = worker_fun, db = parent_db, ...))
  }

  stopifnot(
    "compute must be NULL or the name of a mirai compute profile" = is_scalar_character(compute),
    "no mirai daemons are set for this compute profile; see ?mirai::daemons" = {
      mirai::daemons_set(.compute = compute)
    }
  )

  # one connection per daemon, opened once and reused for every item; a daemon already
  # connected to this database from an earlier call keeps its connection
  mirai::everywhere(
    {
      connected <- exists(".evoland_worker_db", envir = globalenv()) &&
        identical(.evoland_worker_db$catalog, catalog)
      if (!connected) {
        assign(
          ".evoland_worker_db",
          evoland.plus::evoland_db$new(catalog = catalog, data_path = data_path),
          envir = globalenv()
        )
        .evoland_worker_db$execute("set threads to 1")
      }
    },
    catalog = parent_db$catalog,
    data_path = parent_db$data_path,
    .compute = compute
  )

  mirai::mirai_map(
    items,
    function(item, worker_fun, args) {
      do.call(
        worker_fun,
        c(list(item = item, db = get(".evoland_worker_db", envir = globalenv())), args)
      )
    },
    .args = list(worker_fun = worker_fun, args = list(...)),
    .compute = compute
  )[.stop]
}
