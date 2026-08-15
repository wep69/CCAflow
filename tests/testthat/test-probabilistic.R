test_that("probabilistic CCA yields finite model correlation", {
  d <- cca_simulate(n=80,scenario="linear",seed=18)
  f <- cca_probabilistic(d$X,d$Y,ncomp=2,max_iter=100,tol=1e-5)
  expect_true(all(is.finite(f$correlations)))
  expect_true(length(f$log_likelihood) >= 1)
})
