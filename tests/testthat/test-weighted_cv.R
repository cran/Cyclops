library(testthat)

test_that("Small logistic dense regression", {
    n1 <- 100
    n2 <- 100
    p <- 10
    x_prob <- 0.5
    beta <- rnorm(p, mean = 0, sd = 1)
    X <- matrix(rbinom((n1 + n1) * p, 1, x_prob),
                ncol = p)
    exb <- exp(X %*% beta)
    p_i <- exb / (1 + exb)
    Y <- apply(p_i, 1, function(p) { rbinom(1,1,p) })
})
