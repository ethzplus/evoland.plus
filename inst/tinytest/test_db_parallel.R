# Setup temporary DB
library(tinytest)
temp_dir <- tempfile()
db <- evoland_db$new(
  temp_dir,
  reporting_nonce = "vorpal sword"
)

# 1. Test Serial Execution: the default profile has no daemons, so map_evoland() runs serially
expect_false(mirai::daemons_set())
expect_silent(
  res <- map_evoland(
    items = list(1, 2, 3),
    worker_fun = function(item, db) {
      paste(db$read_only, item * 2)
    },
    parent_db = db
  )
)
expect_equal(res, as.list(paste(FALSE, c(2, 4, 6))))

# 2. Test Parallel Execution (Requires installed package for workers)
# In development (pkgload load_all()), daemons can't load the package, so we check if we
# can run a minimal task before proceeding. tinytest::exit_file() does not stop the file
# when it is run by run_test_file() or run_test_dir(), so everything that needs daemons
# hangs off can_run_parallel instead.
can_run_parallel <- FALSE

# the address of the connection a worker was given
connection_of <- function(item, db) format(environment(db$fetch))

if (!at_home() || !requireNamespace("evoland.plus", quietly = TRUE)) {
  message(
    "\n  Skipping parallel tests: workers could not initialize (package likely not installed)"
  )
} else {
  tryCatch(
    {
      # one daemon, so that both calls below land on the same one
      mirai::daemons(1, .compute = "evoland_compute")

      res_parallel <- map_evoland(
        items = list(1, 2),
        worker_fun = function(item, db) db$reporting_t,
        parent_db = db,
        compute = "evoland_compute"
      )
      connections <- c(
        map_evoland(list(1), connection_of, parent_db = db, compute = "evoland_compute"),
        map_evoland(list(2), connection_of, parent_db = db, compute = "evoland_compute")
      )
      can_run_parallel <- TRUE
    },
    error = function(e) {
      message(
        "\n  Skipping parallel tests: workers could not initialize (",
        conditionMessage(e),
        ")"
      )
    },
    finally = mirai::daemons(0, .compute = "evoland_compute")
  )
}

expect_error(
  map_evoland(list(1), connection_of, parent_db = db, compute = "no_daemons_here"),
  "no mirai daemons are set for this compute profile"
)

if (can_run_parallel) {
  # should be exactly the same
  expect_equal(db$reporting_t, res_parallel[[1]])
  expect_equal(db$reporting_t, res_parallel[[2]])
  # the daemon opened one connection and reused it for the second call
  expect_equal(connections[[1]], connections[[2]])

  # daemons on the default profile are used without passing compute, each with an equal share
  # of the cores for DuckDB; compute = NULL still runs serially
  pid_of <- function(item, db) Sys.getpid()
  threads_of <- function(item, db) {
    db$get_query("select current_setting('threads')", as_atomic = TRUE)
  }
  mirai::daemons(2)
  default_pids <- unlist(map_evoland(list(1, 2), pid_of, parent_db = db))
  default_threads <- unlist(map_evoland(list(1, 2), threads_of, parent_db = db))
  db$execute("reset threads")
  n_cores <- db$get_query("select current_setting('threads')", as_atomic = TRUE)
  serial_pid <- unlist(map_evoland(list(1), pid_of, parent_db = db, compute = NULL))
  mirai::daemons(0)
  expect_true(all(default_pids != Sys.getpid()))
  expect_equal(
    unique(default_threads),
    as.numeric(max(1L, as.integer(n_cores) %/% 2L))
  )
  expect_equal(serial_pid, Sys.getpid())

  # 3. Uncoordinated concurrent upserts into a single table.
  # Each worker is its own process with its own connection and no knowledge of the
  # others; the DuckLake catalog plus ducklake_db's retry wrapper have to get all of
  # them through without losing rows or duplicating keys. Four writers is enough to
  # lose rows without the retry -- even two do -- and the retry logic itself is
  # covered deterministically in test_db_ducklake.R.
  conc_dir <- tempfile("evoland_conc_")
  conc_seed <- data.table::data.table(id_worker = 0L, value = 0)
  data.table::setattr(conc_seed, "key_cols", "id_worker")
  ducklake_db$new(conc_dir)$commit(conc_seed, "conc_t", method = "overwrite")

  n_workers <- 4L
  mirai::daemons(n_workers, .compute = "evoland_compute")

  conc_results <- mirai::mirai_map(
    seq_len(n_workers),
    function(id_worker, path) {
      db <- evoland.plus::ducklake_db$new(path = path)
      row <- data.table::data.table(id_worker = as.integer(id_worker), value = id_worker * 1.0)
      data.table::setattr(row, "key_cols", "id_worker")
      tryCatch(
        {
          db$commit(row, "conc_t", method = "upsert")
          "ok"
        },
        error = function(e) conditionMessage(e)
      )
    },
    .args = list(path = conc_dir),
    .compute = "evoland_compute"
  )[]

  mirai::daemons(0, .compute = "evoland_compute")

  expect_equal(unique(unlist(conc_results)), "ok")

  conc_final <- ducklake_db$new(conc_dir)$fetch("conc_t")
  # every writer's row survived exactly once
  expect_equal(sort(conc_final[["id_worker"]]), 0:n_workers)
  expect_equal(conc_final[["value"]][order(conc_final[["id_worker"]])], as.numeric(0:n_workers))

  unlink(conc_dir, recursive = TRUE)
}

# Cleanup
unlink(temp_dir, recursive = TRUE)
