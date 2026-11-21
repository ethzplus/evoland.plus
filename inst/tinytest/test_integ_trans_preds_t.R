library(tinytest)

# Gate: skip during R CMD check; run with build_install_test()
if (!at_home()) {
  exit_file("Integration tests skipped (not at_home)")
}

# Fixture
source(file.path(
  system.file("tinytest", package = "evoland.plus"),
  "helper_testdb.R"
))
db <- make_test_db(
  include_neighbors = FALSE,
  include_trans_preds = FALSE,
  include_alloc_params = FALSE
)

# Test empty table
expect_stdout(print(as_trans_preds_t()), "Transition-Predictor Relationships")

# suppress info logs from mlr3 during testing
lgr::get_logger("mlr3")$set_threshold("warn")

set.seed(123)
# Test covariance filter
expect_message(
  perf_results <- db$get_pred_filter_score(
    # the default performance measure for a classification task is classif.ce
    # with minimize TRUE so we expect scores in [-1,0]
    filter = mlr3filters::FilterPerformance$new(resampling = mlr3::rsmp("cv", folds = 2))
  ),
  "Processing 4 transitions"
)

perf_expected <-
  as_trans_preds_t(data.table::rowwiseDT( # nolint start
      id_run=, id_pred=, id_trans=, performance=,
      0,       1,        1,         -0.12       ,
      0,       1,        2,         -0.02428571 ,
      0,       1,        3,         -0.02073171 ,
      0,       1,        4,         -0.03658537 ,
      0,       2,        1,         -0.12       ,
      0,       2,        2,         -0.02428571 ,
      0,       2,        3,         -0.02073171 ,
      0,       2,        4,         -0.03658537 ,
      0,       3,        1,         -0.12       ,
      0,       3,        2,         -0.02428571 ,
      0,       3,        3,         -0.02073171 ,
      0,       3,        4,         -0.03658537 ,
      0,       4,        1,         -0.12       ,
      0,       4,        2,         -0.02428571 ,
      0,       4,        3,         -0.02073171 ,
      0,       4,        4,         -0.03658537
  )) # nolint end
expect_equal(perf_results, perf_expected, tol = 1e-6)

# Test an importance filter, the usual covariate selection step
set.seed(13233)
expect_message(
  importance_results <- db$get_pred_filter_score(
    filter = mlr3filters::FilterImportance$new(learner = mlr3::lrn("classif.rpart"))
  ),
  "Processing 4 transitions"
)

# rpart leaves a predictor unscored where it does not split on it, and every predictor of the
# rare transitions 3 and 4, where it does not split at all
importance_expected <-
  as_trans_preds_t(data.table::rowwiseDT( # nolint start
      id_run=, id_pred=, id_trans=, importance=,
      0,       1,        1,         15.984901  ,
      0,       1,        2,         8.539657   ,
      0,       1,        3,         NA         ,
      0,       1,        4,         NA         ,
      0,       2,        1,         11.711649  ,
      0,       2,        2,         6.169272   ,
      0,       2,        3,         NA         ,
      0,       2,        4,         NA         ,
      0,       3,        1,         4.853707   ,
      0,       3,        2,         NA         ,
      0,       3,        3,         NA         ,
      0,       3,        4,         NA         ,
      0,       4,        1,         1.534070   ,
      0,       4,        2,         1.676494   ,
      0,       4,        3,         NA         ,
      0,       4,        4,         NA
  )) # nolint end
expect_equal(importance_results, importance_expected, tol = 1e-6)

# get_pred_filter_score restricted to id_trans == 1; expect exactly 1 transition processed
db$set_full_trans_preds(reset = TRUE)

set.seed(123)
expect_message(
  perf_results_manual <- db$get_pred_filter_score(
    filter = mlr3filters::FilterPerformance$new(resampling = mlr3::rsmp("cv", folds = 2)),
    id_trans = 1L
  ),
  "Processing 1 transitions"
)

# Result should contain only id_trans == 1 rows, identical to the full run
expect_equal(
  perf_results_manual,
  perf_results[id_trans == 1L],
  tol = 1e-7
)

# a second filter only scores the pairs the first one kept
db$set_full_trans_preds(reset = TRUE)
importance_results[, is_selected := id_trans != 1L | importance > 5]
suppressMessages(db$trans_preds_t <- importance_results)
set.seed(123)
expect_message(
  perf_refined <- db$get_pred_filter_score(
    filter = mlr3filters::FilterPerformance$new(resampling = mlr3::rsmp("cv", folds = 2))
  ),
  "Processing 4 transitions"
)
expect_equal(perf_refined[id_trans == 1L, id_pred], 1:2)
expect_equal(perf_refined[id_trans == 2L, id_pred], 1:4)
