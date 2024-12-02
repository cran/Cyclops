library("testthat")
library("survival")

suppressWarnings(RNGversion("3.5.0"))

test_that("Check very small Cox example with centering", {
    test <- read.table(header=T, sep = ",", text = "
start, length, event, x1, x2
0, 4,  1,0,0
0, 3.5,1,2,0
0, 3,  0,0,1
0, 2.5,1,0,1
0, 2,  1,1,1
0, 1.5,0,1,0
0, 1,  1,1,0
")

    test$x1c <- test$x1 - mean(test$x1)
    test$x2c <- test$x2 - mean(test$x2)

    gold <- coxph(Surv(length, event) ~ x1 + x2, test)
    summary(gold)

    cent <- coxph(Surv(length, event) ~ x1c + x2c, test)
    summary(cent)
})

test_that("Small Bernoulli dense regression", {
    binomial_bid <- c(1,5,10,20,30,40,50,75,100,150,200)
    binomial_n <- c(31,29,27,25,23,21,19,17,15,15,15)
    binomial_y <- c(0,3,6,7,9,13,17,12,11,14,13)

    log_bid <- log(c(rep(rep(binomial_bid, binomial_n - binomial_y)), rep(binomial_bid, binomial_y)))
    y <- c(rep(0, sum(binomial_n - binomial_y)), rep(1, sum(binomial_y)))

    tolerance <- 1E-4

    glmFit <- glm(y ~ log_bid, family = binomial()) # gold standard

    c_log_bid <- log_bid - mean(log_bid)

    glmCen <- glm(y ~ c_log_bid, family = binomial())

    coef(glmFit)
    coef(glmCen)

    dptr <- createCyclopsData(y ~ log_bid, modelType = "lr")
    cyclopsFit <- fitCyclopsModel(dptr, prior = createPrior("none"),
                                   control = createControl(noiseLevel = "silent"))

    dptrC <- createCyclopsData(y ~ c_log_bid, modelType = "lr")
    cyclopsFitC <- fitCyclopsModel(dptrC, prior = createPrior("none"),
                                  control = createControl(noiseLevel = "silent"))

    coef(cyclopsFit)
    coef(cyclopsFitC)
})

