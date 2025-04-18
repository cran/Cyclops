library("testthat")
library("survival")
library("gnm")

context("test-smallCLR.R")

#
# Small conditional logistic regression
#

test_that("Small conditional logistic regression", {
    gold <- clogit(case ~ spontaneous + induced + strata(stratum), data=infert)

    dataPtr <- createCyclopsData(case ~ spontaneous + induced + strata(stratum),
                                  data = infert,
                                  modelType = "clr")

    cyclopsFit <- fitCyclopsModel(dataPtr, prior = createPrior("none"))

	tolerance <- 1E-4

	expect_equal(coef(cyclopsFit), coef(gold), tolerance = tolerance)
	expect_equal(cyclopsFit$log_likelihood, logLik(gold)[[1]], tolerance = tolerance)

    expect_equal(vcov(cyclopsFit), vcov(gold), tolerance = tolerance)

    expect_equal(aconfint(cyclopsFit), confint(gold), tolerance = tolerance)

    expect_equal(confint(cyclopsFit, c(1:2), includePenalty = TRUE),
                 confint(cyclopsFit, c(1:2), includePenalty = FALSE))

    # dataPtrR <- createCyclopsData(case ~ spontaneous + induced + strata(stratum),
    #                                    data = infert,
    #                                    modelType = "clr")
    #
    # cyclopsFitR <- fitCyclopsModel(dataPtrR,
    #                                prior = createPrior("laplace", 1, exclude = 1))

#     expect_not_equal(confint(cyclopsFitR, c(1), includePenalty = TRUE),
#                  confint(cyclopsFitR, c(1), includePenalty = FALSE))
#       # How to test for inequality?
# This is broken:
# 	expect_equal(predict(cyclopsFit), predict(gold, type = "risk"), tolerance = tolerance)
})

test_that("Small conditional poisson regression with an offset", {

    gold.cp <- gnm(event ~ exgr + agegr + offset(loginterval),
                   family = poisson, eliminate = indiv,
                   data = Cyclops::oxford)
    vcov(gold.cp)

    dataPtr <- createCyclopsData(event ~ exgr + agegr + strata(indiv) + offset(loginterval),
                                 data = Cyclops::oxford,
                                 modelType = "cpr")
    cyclopsFit <- fitCyclopsModel(dataPtr,
                                  prior = createPrior("none"))

    expect_equal(as.numeric(coef(cyclopsFit)), as.numeric(coef(gold.cp)),
                 tolerance = 1E-6)
    expect_equal(as.numeric(vcov(cyclopsFit)), as.numeric(vcov(gold.cp)),
                 tolerance = 1E-6)
})


