test_that("ohoegdm() stores every draw when burnin = 0", {
  set.seed(1)
  y <- matrix(rbinom(60 * 12, 1, 0.5), nrow = 60, ncol = 12)

  fit <- ohoegdm(y, k = 2, m = 2, burnin = 0, chain_length = 20)

  # Every stored -2 log-likelihood is positive; unfilled entries are 0.
  expect_true(all(fit$chain$m2lls > 0))
})

test_that("ohoegdm() rejects responses with missing values", {
  set.seed(1)
  y <- matrix(rbinom(60 * 12, 1, 0.5), nrow = 60, ncol = 12)
  y[1, 1] <- NA

  expect_error(
    ohoegdm(y, k = 2, m = 2, burnin = 5, chain_length = 10),
    "must not contain missing values"
  )
})
