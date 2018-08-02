library(RRembo)
context("mappings")

test_that("mappings",{
  d <- 10; D <- 201
  A <- selectA(d, D, type = 'Gaussian')
  
  n <- 1000
  size <- 10
  Y <- size * (2 * matrix(runif(n * d), n) - 1)
  X <- randEmb(Y, A)
  Z <- ortProj(X, t(A))
  Xback <- mapZX(Z, A)
  
  expect_equal(max(X - Xback), 0, tolerance = 1e-3)
})