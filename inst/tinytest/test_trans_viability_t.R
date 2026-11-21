library(tinytest)

expect_inherits(as_trans_viability_t(), "trans_viability_t")
expect_equal(nrow(as_trans_viability_t()), 0L)
expect_stdout(print(as_trans_viability_t()), "(empty)")

trans_viability_t <- as_trans_viability_t(list(
  id_trans = c(2, 1),
  id_run = c(1, 1),
  is_viable = c(FALSE, TRUE)
))
expect_equal(names(trans_viability_t), c("id_run", "id_trans", "is_viable"))
expect_identical(trans_viability_t[["id_trans"]], 1:2)
expect_stdout(print(trans_viability_t), "allowed: 1, disallowed: 1")

# an override says either way; NA would fall through to the default unnoticed
expect_error(
  as_trans_viability_t(list(id_run = 1L, id_trans = 1L, is_viable = NA)),
  "is_viable must be TRUE or FALSE"
)
# one override per run and transition
expect_error(
  as_trans_viability_t(list(id_run = 1L, id_trans = c(1L, 1L), is_viable = c(TRUE, FALSE)))
)
