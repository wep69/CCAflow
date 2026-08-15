test_that("rank prediction uses training empirical maps", {
  d <- cca_simulate(n=90,scenario="monotone",seed=5)
  f <- cca_rank(d$X,d$Y,ncomp=2)
  pr <- predict(f,d$X[1:8,,drop=FALSE],d$Y[1:8,,drop=FALSE])
  expect_equal(nrow(pr$x),8)
  expect_equal(ncol(pr$x),2)
})

test_that("kernel CCA retains original-variable interpretation and predicts", {
  d <- cca_simulate(n=70,scenario="nonmonotone",seed=6)
  f <- cca_kernel(d$X,d$Y,lambda_x=.2,lambda_y=.2,ncomp=2)
  expect_equal(nrow(f$loadings_x),ncol(d$X))
  pr <- predict(f,d$X[1:10,,drop=FALSE],d$Y[1:10,,drop=FALSE])
  expect_equal(dim(pr$x),c(10,2))
})

test_that("bootstrap returns aligned reference component labels", {
  d <- cca_simulate(n=70,scenario="linear",seed=10)
  f <- cca_ridge(d$X,d$Y,.1,.1,ncomp=2)
  b <- cca_bootstrap(f,B=25,seed=10)
  expect_true(all(unique(b$variable_draws$component) %in% c("CV1","CV2")))
})
