#' Utility functions and methods for `evoland_db`
#'
#' Factored out utility functions and methods specifically used in conjunction with
#' `evoland_db` and its subclasses.
#'
#' @name evoland_db_util
NULL

#' @describeIn evoland_db_util Gets the read expression for a given table, taking into
#' account the run lineage. Whether a table is run-scoped, and which columns form one
#' inheritance slice, is the `inherit_by` of its `table_spec()` in `EVOLAND_SCHEMA`; a table
#' that is not in the registry is run-scoped if its stored metadata carries `inherit_by`.
#' A table that is not run-scoped, or a read without `id_run`, falls back to the default
#' read expression. Otherwise, the read expression filters for `id_run` and its lineage,
#' taking each slice from the closest run in the lineage that holds it.
#' @param self The `evoland_db` instance
#' @param super The superclass (parent class) environment.
#' @param table_name The name of the table to read from
#' @param id_run Integer, the run whose lineage is read; `NULL` reads every run as stored
#' @return A SQL expression string to read most specific data for slice
get_evoland_db_read_expr <- function(self, super, table_name, id_run) {
  base_read_expr <- super$get_read_expr(table_name)
  if (is.null(id_run)) {
    return(base_read_expr)
  }

  inherit_by <- if (table_name %in% names(EVOLAND_SCHEMA)) {
    EVOLAND_SCHEMA[[table_name]][["inherit_by"]]
  } else {
    self$get_table_metadata(table_name)[["inherit_by"]]
  }
  if (is.null(inherit_by)) {
    return(base_read_expr)
  }

  distinctness_cols <- c("id_run", inherit_by)

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
  if ("id_period" %in% inherit_by) {
    join_conditions <- lapply(
      inherit_by,
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
      group_cols = glue::glue_sql_collapse(sql("b.{`inherit_by`}"), sep = ", ")
    )
  } else if (length(inherit_by) > 0L) {
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
      group_cols = glue::glue_sql_collapse(sql("b.{`inherit_by`}"), sep = ", ")
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

#' @describeIn evoland_db_util Apply a worker function to each item, on mirai daemons when
#' they are set and serially otherwise, and return the results once all items are done.
#'
#' Turn parallelism on with `mirai::daemons()`, e.g. `mirai::daemons(4)`, and off with
#' `mirai::daemons(0)`; nothing else changes. With daemons set for `compute`, each daemon opens
#' one connection to the same catalog and data files as `parent_db` and reuses it for every
#' item it is given; evoland.plus must be installed where the daemons run. DuckDB on each
#' daemon gets an equal share of its machine's cores, so that the daemons together do not
#' oversubscribe it. Without daemons, the worker gets `parent_db` itself. Workers must not rely
#' on anything but the item and the database: an item carries the run it is for.
#'
#' @param worker_fun A function to apply to each item, with signature `worker_fun(item,
#' db, ...)`.
#' @param items A list of items to iterate over
#' @param parent_db An [evoland_db] instance: passed to the worker in the serial case, and
#' giving the catalog and data path the daemons connect to otherwise.
#' @param compute Name of the mirai compute profile to run on, see [mirai::daemons()]. The
#' default profile, `"default"`, runs serially while it has no daemons; any other profile
#' must have daemons set. `NULL` always runs serially.
#' @param ... Additional arguments passed to `worker_fun`.
#'
#' @return A list of results
#' @export
map_evoland <- function(
  items,
  worker_fun,
  parent_db,
  compute = "default",
  ...
) {
  stopifnot(
    "compute must be NULL or the name of a mirai compute profile" = is.null(compute) ||
      is_scalar_character(compute)
  )
  serial <- is.null(compute) ||
    (identical(compute, "default") && !mirai::daemons_set(.compute = compute))
  if (serial) {
    return(lapply(X = items, FUN = worker_fun, db = parent_db, ...))
  }
  stopifnot(
    "no mirai daemons are set for this compute profile; see ?mirai::daemons" = {
      mirai::daemons_set(.compute = compute)
    }
  )

  # one connection per daemon, opened once and reused for every item; a daemon already
  # connected to this database from an earlier call keeps its connection. Its DuckDB threads
  # are set on every call, as the number of daemons may have changed.
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
      }
      # DuckDB's default thread count is the number of cores it can use
      .evoland_worker_db$execute("reset threads")
      cores <- .evoland_worker_db$get_query(
        "select current_setting('threads')",
        as_atomic = TRUE
      )
      threads <- max(1L, as.integer(cores) %/% n_daemons)
      .evoland_worker_db$execute(sprintf("set threads to %d", threads))
    },
    catalog = parent_db$catalog,
    data_path = parent_db$data_path,
    n_daemons = mirai::info(.compute = compute)[["connections"]],
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

# The number of daemons map_evoland() spreads items over on `compute`; 1 when it runs serially
n_evoland_workers <- function(compute = "default") {
  if (is.null(compute) || !mirai::daemons_set(.compute = compute)) {
    return(1L)
  }
  mirai::info(.compute = compute)[["connections"]]
}
