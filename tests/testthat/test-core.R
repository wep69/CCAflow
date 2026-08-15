test_that("classical CCA agrees with base cancor correlations", {
  d <- cca_simulate(n=120,p=5,q=4,scenario="linear",seed=2)
  f <- cca_classic(d$X,d$Y,ncomp=4)
  b <- cancor(scale(d$X),scale(d$Y))
  expect_equal(f$correlations, b$cor[seq_along(f$correlations)], tolerance=1e-6)
  expect_equal(nrow(f$loadings_x), ncol(d$X))
  expect_equal(nrow(f$cross_loadings_y), ncol(d$Y))
})

test_that("ridge CCA handles high dimension", {
  d <- cca_simulate(n=35,p=45,q=40,scenario="highdim",seed=3)
  f <- cca_ridge(d$X,d$Y,.2,.2,ncomp=3)
  expect_true(all(is.finite(f$correlations)))
  expect_lte(max(f$correlations),1+1e-8)
})

test_that("categorical inputs are not silently coerced", {
  X <- data.frame(a=1:10,b=factor(rep(c("a","b"),5)))
  Y <- matrix(rnorm(20),10,2)
  expect_error(cca_classic(X,Y),"numeric")
})
