library(CCAflow)
d <- cca_agri_data(seed = 42)
X <- as.matrix(d$X)
Y <- as.matrix(d$Y)

audit <- cca_audit(X, Y)
fit <- cca_ridge(X, Y, lambda_x = 0.15, lambda_y = 0.15, ncomp = 3)
inf <- cca_permutation(fit, B = 199, seed = 42)
boot <- cca_bootstrap(fit, B = 199, seed = 42)
cv <- cca_validate(X, Y, method = "ridge",
                   params = list(lambda_x = 0.15, lambda_y = 0.15, ncomp = 3),
                   folds = 5, repeats = 3, seed = 42)
cca_table(fit, "loadings")
cca_plot(fit, "dual_circle")
cca_plot(cv)
