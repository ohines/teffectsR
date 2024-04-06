N <- 50
X <- rnorm(N) # confounder
A <- rbinom(N, 1, plogis(X)) # treatment variable
Y <- X + 0.5 * A # continuous outcome
Z <- rbinom(N, 1, plogis(Y)) # binary outcome
df <- data.frame(X = X, A = A, Y = Y, Z = Z)

test_that("AIPW runs without errors", {
  # Fits AIPW by default
  teffect(A ~ X, Y ~ X, data = df)
  teffect(A ~ X, Z ~ X, data = df, outcome.family = "binomial")
})

test_that("IPW runs without errors", {
  teffect(A ~ X, Y ~ 1, data = df, method = "IPW")
  teffect(A ~ X, Z ~ 1, data = df, method = "IPW")
})

test_that("RA runs without errors", {
  teffect(A ~ 1, Y ~ X, data = df, method = "RA")
  teffect(A ~ 1, Z ~ X, data = df, outcome.family = "binomial", method = "RA")
})

test_that("PSMatch runs without errors", {
  teffect(A ~ X, Y ~ 1, data = df, method = "PSMatch")
  teffect(A ~ X, Z ~ 1, data = df, method = "PSMatch")
})
