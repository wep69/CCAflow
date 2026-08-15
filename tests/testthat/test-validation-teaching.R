test_that("validation creates held-out correlations", {
  d <- cca_simulate(n=80,scenario="linear",seed=12)
  cv <- cca_validate(d$X,d$Y,"ridge",list(lambda_x=.1,lambda_y=.1,ncomp=2),folds=4,repeats=2,seed=12)
  expect_gt(nrow(cv$results),0)
  expect_true(all(cv$results$fold %in% 1:4))
})

test_that("grouped validation keeps groups together", {
  d <- cca_simulate(n=80,scenario="linear",seed=13)
  g <- rep(1:20,each=4)
  cv <- cca_validate_grouped(d$X,d$Y,g,"ridge",list(lambda_x=.1,lambda_y=.1,ncomp=1),folds=4,repeats=1,seed=13)
  expect_gt(nrow(cv$results),0)
})

test_that("tour keyframes are reproducible matrices", {
  d <- cca_agri_data(n=60,seed=14)
  f <- cca_classic(d$X,d$Y,ncomp=2)
  tr <- cca_tour(f,"x",frames=12)
  k <- cca_tour_keyframes(tr,n=4)
  expect_equal(length(k$bases),4)
  expect_equal(nrow(k$bases[[1]]),ncol(d$X))
})
