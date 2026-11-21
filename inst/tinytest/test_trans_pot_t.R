library(tinytest)

# --------------------------------------------------------------------------
# Unit tests for trans_pot_t schema
# --------------------------------------------------------------------------

# Values must remain in [0, 1]
expect_error(
  as_trans_pot_t(data.frame(
    id_run = 0L,
    id_trans = c(1L, 2L, 1L, 2L),
    id_period_post = 2L,
    id_coord = c(1L, 1L, 2L, 2L),
    value = c(1.4, 0.4, 0.6, 0.3)
  )),
  'all(x[["value"]] <= 1) is not TRUE',
  fixed = TRUE
)
