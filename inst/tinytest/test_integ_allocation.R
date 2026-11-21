library(tinytest)

if (!at_home()) {
  exit_file("Integration tests skipped (not at_home)")
}

# The allocation loop around the allocators: prediction, the allocation context, committing
# the map and the neighbor predictors, runs and their lineage. It is the same for every
# allocator, so the deterministic greedy allocator stands in for all of them; the allocators
# themselves are tested on hand-built contexts in test_alloc_*.R

source(file.path(system.file("tinytest", package = "evoland.plus"), "helper_testdb.R"))
db <- make_test_db()

db$trans_rates_t <- db$get_obs_trans_rates()
db$trans_rates_t <- extrapolate_trans_rates(db$trans_rates_t, db$periods_t)

test_learner <- mlr3::lrn("classif.featureless", predict_type = "prob")
test_measures <- c("classif.auc")

# test the package's featureless learner fit and append to disk
expect_message(
  db$trans_models_t <- db$fit_partial_models(
    learner = test_learner,
    measures = test_measures,
    seed = 1244244
  ),
  "Fitting partial models for 4 transitions..."
)

# Score-select mode: pick best partial model by classif.auc and retrain on full data
expect_message(
  db$trans_models_t <- db$fit_full_models(
    select_score = "classif.auc",
    select_maximize = TRUE
  ),
  "Fitting full models for"
)

# helpers: the rows of a run's map in a period, read through the run's lineage
lulc_of <- function(id_run, id_period) {
  db$get_query(
    "select id_coord, id_lulc from {db$get_read_expr('lulc_data_t', id_run = id_run)}
     where id_period = {id_period} order by id_coord"
  )
}
viable <- db$trans_meta_t[is_viable == TRUE]
# the cells that changed between an anterior and a posterior map, by transition; the anterior
# map is the observed period 3 unless given
changes_of <- function(posterior, anterior = lulc_of(0L, 3L)) {
  anterior[posterior, on = "id_coord"][
    id_lulc != i.id_lulc,
    .N,
    by = .(id_lulc_anterior = id_lulc, id_lulc_posterior = i.id_lulc)
  ]
}
# the number of potentials predicted from an anterior map: one per viable transition and cell
# of its anterior class
n_potentials <- function(anterior) {
  nrow(anterior[viable, on = .(id_lulc = id_lulc_anterior), nomatch = NULL, allow.cartesian = TRUE])
}
# the changes per viable transition, to compare with the quota of the allocation context
changes_by_trans <- function(posterior, anterior = lulc_of(0L, 3L)) {
  changes_of(posterior, anterior)[
    viable,
    on = .(id_lulc_anterior, id_lulc_posterior),
    .(id_trans, N = data.table::fcoalesce(N, 0L))
  ]
}

# --------------------------------------------------------------------------
# allocate(): run 1 predicts its potentials, allocates and commits period 4
# --------------------------------------------------------------------------
expect_equal(nrow(lulc_of(1L, 4L)), 0L)
expect_message(
  allocate(
    db,
    id_periods = 4L,
    allocator = alloc_greedy_spec(),
    select_score = "classif.auc",
    select_maximize = TRUE,
    id_run = 1L
  ),
  "Allocating with greedy (joint): run 1",
  fixed = TRUE
)
expect_equal(nrow(lulc_of(1L, 4L)), 900L)
# only viable transitions happen
expect_true(all(
  changes_of(lulc_of(1L, 4L))[!viable, on = .(id_lulc_anterior, id_lulc_posterior), .N] == 0L
))
# the potentials were predicted under the allocated run
expect_equal(db$row_count("trans_pot_t"), n_potentials(lulc_of(0L, 3L)))
expect_equal(
  db$get_query("select distinct id_run from dl_db.trans_pot_t", as_atomic = TRUE),
  1L
)
adj_pots <- db$adjusted_trans_pot_v(4L, id_run = 1L)
expect_equal(nrow(adj_pots), n_potentials(lulc_of(0L, 3L)))
expect_true(all(adj_pots$value >= 0 & adj_pots$value <= 1))

# a run with its own models does not take an ancestor's potentials; run 3 predicts its own
expect_message(
  allocate(
    db,
    id_periods = 4L,
    allocator = alloc_greedy_spec(),
    select_score = "classif.auc",
    select_maximize = TRUE,
    id_run = 3L
  ),
  "Predicting transition potential"
)
expect_equal(nrow(lulc_of(3L, 4L)), 900L)

# force_predict predicts again, also with prefetching
options("evoland.plus.use_prefetch_predict" = TRUE)
expect_message(
  allocate(
    db,
    id_periods = 4L,
    allocator = alloc_greedy_spec(),
    select_score = "classif.auc",
    select_maximize = TRUE,
    id_run = 3L,
    force_predict = TRUE
  ),
  "Predicting transition potential"
)
options("evoland.plus.use_prefetch_predict" = NULL)

# --------------------------------------------------------------------------
# several periods: run 2 allocates period 3 from the observed period 2, then period 4 from its
# own period 3, with the neighbor predictors recomputed in between
# --------------------------------------------------------------------------
expect_message(
  allocate(
    db,
    id_periods = 3:4,
    allocator = alloc_greedy_spec(),
    select_score = "classif.auc",
    select_maximize = TRUE,
    id_run = 2L
  ),
  "periods 3 -> 4"
)
expect_equal(nrow(lulc_of(2L, 3L)), 900L)
expect_equal(nrow(lulc_of(2L, 4L)), 900L)
# period 4 starts from the run's own allocation of period 3, not from the observed map
ctx_4 <- evoland.plus:::allocation_context(db, 4L, 2L, 2L)
expect_identical(ctx_4$anterior, lulc_of(2L, 3L))
expect_false(isTRUE(all.equal(lulc_of(2L, 3L), lulc_of(0L, 3L))))
# greedy fills each period's quota exactly
expect_equal(
  changes_by_trans(lulc_of(2L, 3L), anterior = lulc_of(0L, 2L)),
  evoland.plus:::allocation_context(db, 3L, 2L, 2L)$quota[, .(id_trans, N = quota)]
)
expect_equal(
  changes_by_trans(lulc_of(2L, 4L), anterior = lulc_of(2L, 3L)),
  ctx_4$quota[, .(id_trans, N = quota)]
)
# potentials were predicted, and neighbor predictors written, for both periods under run 2
expect_equal(
  db$get_query(
    "select id_period_post, count(*) as n from dl_db.trans_pot_t where id_run = 2
     group by all order by all"
  ),
  data.table::data.table(
    id_period_post = 3:4,
    n = c(n_potentials(lulc_of(0L, 2L)), n_potentials(lulc_of(2L, 3L)))
  )
)
expect_equal(
  db$get_query(
    "select distinct id_period from dl_db.pred_data_t where id_run = 2 order by all",
    as_atomic = TRUE
  ),
  3:4
)

# --------------------------------------------------------------------------
# a branch: run 4 branches off run 3 at period 4, and allocates on its parent's potentials
# without predicting
# --------------------------------------------------------------------------
db$runs_t <- as_runs_t(rbind(
  db$runs_t,
  list(id_run = 4, parent_id_run = 3, description = "stochastic alloc", id_period_eval_start = 4),
  fill = TRUE
))
alloc_messages <- character()
withCallingHandlers(
  allocate(
    db,
    id_periods = 4L,
    allocator = alloc_greedy_spec(),
    select_score = "classif.auc",
    select_maximize = TRUE,
    id_run = 4L
  ),
  message = function(m) {
    alloc_messages <<- c(alloc_messages, conditionMessage(m))
    invokeRestart("muffleMessage")
  }
)
expect_false(any(grepl("Predicting transition", alloc_messages)))
expect_equal(nrow(lulc_of(4L, 4L)), 900L)
expect_equal(
  db$get_query("select count(*) from dl_db.trans_pot_t where id_run = 4", as_atomic = TRUE),
  0
)

# --------------------------------------------------------------------------
# allocate_one_period(): returns the map without committing it
# --------------------------------------------------------------------------
n_lulc_before <- db$row_count("lulc_data_t")
expect_message(
  one_period <- allocate_one_period(
    db,
    id_period_post = 4L,
    allocator = alloc_greedy_spec(),
    select_score = "classif.auc",
    select_maximize = TRUE,
    id_run = 4L,
    potentials_ancestor = 1L
  ),
  "Changed [0-9]+ of 900 cells"
)
expect_inherits(one_period, "lulc_data_t")
expect_equal(nrow(one_period), 900L)
expect_equal(unique(one_period[["id_run"]]), 4L)
expect_equal(db$row_count("lulc_data_t"), n_lulc_before)
# greedy is deterministic and fills the quota the context gives every allocator exactly
expect_equal(
  changes_by_trans(one_period),
  evoland.plus:::allocation_context(db, 4L, 4L, 3L)$quota[, .(id_trans, N = quota)]
)

# --------------------------------------------------------------------------
# update_neighbors = FALSE skips the neighbor predictors after the last requested period
# --------------------------------------------------------------------------
db$runs_t <- as_runs_t(rbind(
  db$runs_t,
  list(
    id_run = 5,
    parent_id_run = 3,
    description = "no neighbour update",
    id_period_eval_start = 4
  ),
  fill = TRUE
))
count_run_5_preds <- function() {
  db$get_query(
    "select count(*) from dl_db.pred_data_t where id_run = 5 and id_period = 4",
    as_atomic = TRUE
  )
}
alloc_run_5 <- function(update_neighbors, checkpoint = FALSE) {
  allocate(
    db,
    id_periods = 4L,
    allocator = alloc_greedy_spec(),
    select_score = "classif.auc",
    select_maximize = TRUE,
    id_run = 5L,
    update_neighbors = update_neighbors,
    checkpoint = checkpoint
  )
}
expect_message(alloc_run_5(update_neighbors = FALSE), "Allocating with greedy")
expect_equal(count_run_5_preds(), 0)
expect_message(alloc_run_5(update_neighbors = TRUE), "Allocating with greedy")
expect_true(count_run_5_preds() > 0)

# checkpoint = TRUE checkpoints once after the last period: with a retention of zero, every
# older snapshot is expired. DuckLake cannot unset the option again; nothing below checkpoints
count_snapshots <- function() {
  db$get_query("select count(*) from ducklake_snapshots(dl_db)", as_atomic = TRUE)
}
db$execute("call dl_db.set_option('expire_older_than', '0 seconds')")
# allocating one period commits once
snapshots_before <- count_snapshots()
expect_message(alloc_run_5(update_neighbors = FALSE), "Allocating with greedy")
expect_equal(count_snapshots(), snapshots_before + 1)
# the allocation commits once more, then the checkpoint expires every snapshot before its own
# two: one flushes the inlined rows to files, one merges adjacent files
expect_message(
  alloc_run_5(update_neighbors = FALSE, checkpoint = TRUE),
  "Allocating with greedy"
)
expect_equal(count_snapshots(), 2)

# --------------------------------------------------------------------------
# several runs: run 7 branches off run 6 at period 4, which allocate() reads from runs_t
# --------------------------------------------------------------------------
db$runs_t <- as_runs_t(rbind(
  db$runs_t,
  data.table::data.table(
    id_run = 6:7,
    parent_id_run = c(0L, 6L),
    description = c("trunk", "branch"),
    id_period_eval_start = c(NA, 4L)
  ),
  fill = TRUE
))
expect_message(
  allocate(
    db,
    id_periods = 3:4,
    allocator = alloc_greedy_spec(),
    select_score = "classif.auc",
    select_maximize = TRUE,
    id_run = 6:7
  ),
  "run 7, periods 4\n",
  fixed = TRUE
)
own_periods <- function(id_run) {
  db$get_query(
    "select distinct id_period from dl_db.lulc_data_t where id_run = {id_run} order by all",
    as_atomic = TRUE
  )
}
expect_equal(own_periods(6L), 3:4)
expect_equal(own_periods(7L), 4L)
# the branch starts from the trunk's period 3; greedy then allocates the same period 4
expect_identical(evoland.plus:::allocation_context(db, 4L, 7L, 7L)$anterior, lulc_of(6L, 3L))
expect_equal(lulc_of(7L, 4L), lulc_of(6L, 4L))
# where it branches off, run 7 shares run 6's potentials and allocation context
expect_equal(
  evoland.plus:::allocation_schedule(db, 3:4, 6:7)[, .(
    id_period,
    id_run,
    potentials_run,
    context_run
  )],
  data.table::data.table(
    id_period = c(3L, 4L, 4L),
    id_run = c(6L, 6L, 7L),
    potentials_run = 6L,
    context_run = 6L
  )
)
expect_equal(
  db$get_query("select count(*) from dl_db.trans_pot_t where id_run = 7", as_atomic = TRUE),
  0
)
expect_error(
  allocate(
    db,
    id_periods = 3L,
    allocator = alloc_greedy_spec(),
    select_score = "classif.auc",
    select_maximize = TRUE,
    id_run = 7L
  ),
  "Runs 7 branch off after the last of id_periods"
)

# --------------------------------------------------------------------------
# rng_seed: a seeded run replays its maps whichever runs it is allocated with, and the
# caller's generator is left alone
# --------------------------------------------------------------------------
db$runs_t <- as_runs_t(rbind(
  db$runs_t,
  data.table::data.table(
    id_run = 8:11,
    parent_id_run = 3L,
    description = "clumpy realisation",
    rng_seed = c(42L, 42L, 43L, 42L),
    id_period_eval_start = 4L
  ),
  fill = TRUE
))
alloc_clumpy <- function(id_run) {
  suppressMessages(allocate(
    db,
    id_periods = 4L,
    allocator = alloc_clumpy_spec(),
    select_score = "classif.auc",
    select_maximize = TRUE,
    id_run = id_run,
    update_neighbors = FALSE
  ))
}
set.seed(1)
caller_state <- .Random.seed
alloc_clumpy(8:10)
expect_identical(.Random.seed, caller_state)
# the same seed gives the same map, a different one another map
expect_equal(lulc_of(9L, 4L), lulc_of(8L, 4L))
expect_false(isTRUE(all.equal(lulc_of(10L, 4L), lulc_of(8L, 4L))))
# on its own, a run replays the map it got within the batch
alloc_clumpy(11L)
expect_equal(lulc_of(11L, 4L), lulc_of(8L, 4L))

# allocate_one_period() allocates several runs on their parent's context, with their seeds
shared <- suppressMessages(allocate_one_period(
  db,
  id_period_post = 4L,
  allocator = alloc_clumpy_spec(),
  select_score = "classif.auc",
  select_maximize = TRUE,
  id_run = c(8L, 10L),
  potentials_ancestor = 1L
))
expect_equal(
  as.data.frame(shared[id_run == 8L, .(id_coord, id_lulc)]),
  as.data.frame(lulc_of(8L, 4L))
)
expect_equal(
  as.data.frame(shared[id_run == 10L, .(id_coord, id_lulc)]),
  as.data.frame(lulc_of(10L, 4L))
)
expect_error(
  allocate_one_period(
    db,
    id_period_post = 4L,
    allocator = alloc_clumpy_spec(),
    select_score = "classif.auc",
    select_maximize = TRUE,
    id_run = c(8L, 1L),
    potentials_ancestor = 1L
  ),
  "if allocating for multiple runs, they must share their ancestor"
)

# on mirai daemons, a seeded run gets the same map as on the calling process. The daemons load
# the installed package, which pkgload::load_all() does not give them, so probe them first, as
# test_db_parallel.R does
mirai::daemons(2)
can_run_parallel <- tryCatch(
  isTRUE(map_evoland(list(1L), function(item, db) TRUE, parent_db = db)[[1L]]),
  error = function(e) FALSE
)
if (can_run_parallel) {
  db$runs_t <- as_runs_t(rbind(
    db$runs_t,
    data.table::data.table(
      id_run = 12:13,
      parent_id_run = 3L,
      description = "clumpy realisation on daemons",
      rng_seed = c(42L, 43L),
      id_period_eval_start = 4L
    ),
    fill = TRUE
  ))
  tryCatch(alloc_clumpy(12:13), finally = mirai::daemons(0))
  expect_equal(lulc_of(12L, 4L), lulc_of(8L, 4L))
  expect_equal(lulc_of(13L, 4L), lulc_of(10L, 4L))
} else {
  mirai::daemons(0)
  message("Skipping the allocation on daemons: they cannot load the installed package")
}

# --------------------------------------------------------------------------
# viability per run (ethzplus/evoland.plus#142): run 20 disallows transition 2 (1 -> 3), which
# its parent allows, and run 21 allows it again under run 20. Each calibrates, predicts and
# allocates with its own set; their ancestors and siblings keep theirs
# --------------------------------------------------------------------------
db$runs_t <- as_runs_t(rbind(
  db$runs_t,
  data.table::data.table(
    id_run = 20:22,
    parent_id_run = c(0L, 20L, 3L),
    description = c("without 1 -> 3", "with 1 -> 3 again", "branch without 1 -> 3"),
    id_period_eval_start = c(NA, NA, 4L)
  ),
  fill = TRUE
))
db$trans_viability_t <- as_trans_viability_t(data.table::data.table(
  id_run = 20:22,
  id_trans = 2L,
  is_viable = c(FALSE, TRUE, FALSE)
))
expect_equal(viable[id_trans == 2L, c(id_lulc_anterior, id_lulc_posterior)], c(1L, 3L))

# calibration: predictor sets and models for the three transitions run 20 allows; the partial
# models of transition 2 it inherits from run 0 are not refitted
db$set_full_trans_preds(id_run = 20L)
expect_equal(
  sort(unique(db$fetch("trans_preds_t", where = "id_run = 20")[["id_trans"]])),
  c(1L, 3L, 4L)
)
# run 20 still inherits run 0's predictor set of transition 2, which scoring skips
run_20_scores <- suppressMessages(db$get_pred_filter_score("correlation", id_run = 20L))
expect_equal(sort(unique(run_20_scores[["id_trans"]])), c(1L, 3L, 4L))
expect_message(
  db$trans_models_t <- db$fit_partial_models(
    learner = test_learner,
    measures = test_measures,
    seed = 1244244,
    id_run = 20L
  ),
  "Fitting partial models for 3 transitions..."
)
expect_message(
  db$trans_models_t <- db$fit_full_models(select_score = "classif.auc", id_run = 20L),
  "Fitting full models for 3 transitions..."
)

# prediction and allocation: run 21 reads the model of transition 2 from run 0
suppressMessages(allocate(
  db,
  id_periods = 4L,
  allocator = alloc_greedy_spec(),
  select_score = "classif.auc",
  select_maximize = TRUE,
  id_run = 20:21
))
own_pot_trans <- function(id_run) {
  db$get_query(
    "select distinct id_trans from dl_db.trans_pot_t where id_run = {id_run} order by all",
    as_atomic = TRUE
  )
}
expect_equal(own_pot_trans(20L), c(1L, 3L, 4L))
expect_equal(own_pot_trans(21L), 1:4)
ctx_20 <- evoland.plus:::allocation_context(db, 4L, 20L, 20L)
expect_equal(ctx_20$viable_trans[["id_trans"]], c(1L, 3L, 4L))
# greedy meets the quota of the allowed transitions, and no cell takes the disallowed one
expect_equal(
  changes_by_trans(lulc_of(20L, 4L)),
  ctx_20$quota[
    viable[, .(id_trans)],
    on = "id_trans",
    .(id_trans, N = data.table::fcoalesce(quota, 0L))
  ]
)
ctx_21 <- evoland.plus:::allocation_context(db, 4L, 21L, 21L)
expect_true(ctx_21$quota[id_trans == 2L, quota] > 0L)
expect_equal(changes_by_trans(lulc_of(21L, 4L)), ctx_21$quota[, .(id_trans, N = quota)])

# the parent and the siblings keep the default set
expect_true(all(db$trans_meta_t[["is_viable"]]))
for (run in c(0L, 1L, 3L)) {
  expect_equal(
    evoland.plus:::allocation_context(db, 4L, run, run)$viable_trans[["id_trans"]],
    viable[["id_trans"]]
  )
}

# a branch with its own overrides predicts its own potentials, where run 4 takes run 3's
expect_equal(
  evoland.plus:::allocation_schedule(db, 4L, c(4L, 22L))[,
    .(id_run, potentials_run, context_run)
  ],
  data.table::data.table(
    id_run = c(4L, 22L),
    potentials_run = c(3L, 22L),
    context_run = c(3L, 22L)
  )
)

# --------------------------------------------------------------------------
# Argument checks
# --------------------------------------------------------------------------
expect_error(
  allocate(
    db,
    id_periods = c(1L, 3L),
    allocator = alloc_greedy_spec(),
    select_score = "classif.auc",
    select_maximize = TRUE
  ),
  "id_periods must be contiguous"
)
expect_error(
  allocate(
    db,
    id_periods = 4L,
    allocator = alloc_greedy_spec(),
    select_score = "classif.auc",
    select_maximize = TRUE,
    id_run = c(1L, 1L)
  ),
  "id_run must be distinct integers"
)
expect_error(
  allocate_one_period(
    db,
    id_period_post = 4L,
    allocator = alloc_greedy_spec(),
    select_score = "classif.auc",
    select_maximize = TRUE,
    id_run = 4L,
    potentials_ancestor = 5L
  ),
  "Run 4 has no ancestor 5 generations up"
)
expect_error(
  allocate(db, id_periods = 4L, allocator = "clumpy"),
  "allocator must be an allocator specification"
)
expect_error(
  evoland.plus:::new_allocator_spec("x", list(), one_period = function() NULL),
  "one_period must be a function of one argument"
)
