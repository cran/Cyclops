library(testthat)

test_that("Stratified Cox with cross validation", {
    skip("Only manual debugging")
    
    load(file = "crash.rda")
    cyclopsData <- convertToCyclopsData(data, covariates, modelType="cox", quiet=TRUE)
    
    
    prior = createPrior("laplace", useCrossValidation = TRUE, exclude = 1)
    control = createControl(noiseLevel = "quiet", seed = 666, resetCoefficients = TRUE)
    fit <- fitCyclopsModel(cyclopsData, prior = prior, control = control) 
    test1 <- coef(fit)[coef(fit) != 0.0]
       
    prior <- createPrior("laplace", 0.125992104989487, exclude = 1)
    fit0 <- fitCyclopsModel(cyclopsData, forceNewObject = TRUE, prior = prior, 
                            startingCoefficients = c(2.0, -1.0, rep(0, 13716)),
                            control = control)
    test2 <- coef(fit0)[coef(fit0) != 0.0]
    
    tolerance <- 1E-4
    expect_equal(test1, test2, tolerance = tolerance)
        
    w0 <- c(1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1)
    w1 <- c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,0)
    
    prior <- createPrior("laplace",  0.01, exclude = 1)
    control = createControl(noiseLevel = "quiet", seed = 666, resetCoefficients = TRUE)
    fit <- fitCyclopsModel(cyclopsData, forceNewObject = TRUE, weights = w0, prior = prior, control = control)      
    pred0 <- getCyclopsPredictiveLogLikelihood(fit, 1 - w0)
    pred0
    expect_equal(pred0, 0, tolerance = tolerance)
    
    prior <- createPrior("laplace",  0.01, exclude = 1)
    control = createControl(noiseLevel = "quiet", seed = 666, resetCoefficients = TRUE)
    fit <- fitCyclopsModel(cyclopsData, forceNewObject = TRUE, weights = w1, prior = prior, control = control)
    pred1 <- getCyclopsPredictiveLogLikelihood(fit, 1 - w1)    
    pred1
    expect_equal(pred1, 0, tolerance = tolerance)
        
    prior = createPrior("laplace", useCrossValidation = TRUE, exclude = 1)
    control = createControl(noiseLevel = "quiet", seed = 666, 
                            selectorType = "byPid",
                            resetCoefficients = TRUE)
    fit <- fitCyclopsModel(cyclopsData, prior = prior, control = control)      
    
    
})

test_that("Use KKT swindle on Big Data", {
    skip_on_cran() # Do not run on CRAN    
    skip("Only for manual testing")
    
    library(Cyclops)
    library(ffbase)
    load.ffdf("~/Dropbox/Projects/ForMarc")
    
    cyclopsData <- convertToCyclopsData(cohortSubset,
                                        covariateSubset,
                                        modelType="lr",
                                        quiet=TRUE,
                                        checkSorting = FALSE)
    variance <- 0.0138 # Found using CV
    
    prior = createPrior("laplace",
                        variance,
                        useCrossValidation = FALSE,
                        exclude=0)
    
    control = createControl(noiseLevel = "quiet",
                            useKKTSwindle = FALSE,
                            tuneSwindle = 1000
    )
    start <- Sys.time()
    cyclopsFit <- fitCyclopsModel(cyclopsData,
                                  forceNewObject = TRUE,
                                  prior = prior,
                                  control = control)
    delta <- Sys.time() - start
    cyclopsFit$timeFit
    delta
        
    
    library(Cyclops)
    library(ffbase)
    load.ffdf("~/Dropbox/Projects/ForMarc")
    
    cyclopsData <- convertToCyclopsData(cohortSubset,
                                        covariateSubset,
                                        modelType="lr",
                                        quiet=TRUE,
                                        checkSorting = FALSE)
    variance <- 0.0138 # Found using CV
    
    prior = createPrior("laplace",
                        variance,
                        useCrossValidation = FALSE,
                        exclude=0)
    
    
    
    control = createControl(noiseLevel = "quiet",
                            useKKTSwindle = FALSE,
                            maxIterations = 10,
                            tuneSwindle = 1000
    )
    start <- Sys.time()
    cyclopsFit <- fitCyclopsModel(cyclopsData,
                                  forceNewObject = TRUE,
                                  prior = prior,
                                  control = control)
    delta <- Sys.time() - start   
    cyclopsFit$timeFit
    delta
})

test_that("Medium Poisson", {
    skip_on_cran() # Do not run on CRAN    
    skip("Only for manual testing")
    seed <- 666
    tolerance <- 1E-4
    
    data <- simulateCyclopsData(nstrata = 1,
                                nrows = 10000,
                                ncovars = 20000,
                                zeroEffectSizeProp = 0.99,
                                model = "poisson")
    save(data, file = "poissonCyclopsData.RData")
    
    load("poissonCyclopsData.RData")
    
    cyclopsData <- convertToCyclopsData(data$outcomes,
                                        data$covariates,
                                        modelType = "pr",
                                        addIntercept = TRUE) 
    
    system.time(
        slowFit <- fitCyclopsModel(cyclopsData,
                                   forceNewObject = TRUE, # Cold start for fair comparison
                                   prior = createPrior("laplace", variance = 0.01, exclude = c(0)))
    )            
    slowFit$timeFit
    coef(slowFit)[1:10]    
})

test_that("Use KKT swindle", {
    skip_on_cran() # Do not run on CRAN    
    skip("Only for manual testing")
    seed <- 666
    tolerance <- 1E-4
    
    data <- simulateCyclopsData(nstrata = 1,
                                nrows = 10000,
                                ncovars = 20000,
                                zeroEffectSizeProp = 0.99,
                                model = "logistic")
    
    load("someCyclopsData.RData")
    
    cyclopsData <- convertToCyclopsData(data$outcomes,
                                        data$covariates,
                                        modelType = "lr",
                                        addIntercept = TRUE) 
    
    system.time(
        fastFit <- fitCyclopsModel(cyclopsData,
                                   forceNewObject = TRUE, # Cold start for fair comparison
                                   prior = createPrior("laplace", variance = 0.01, exclude = c(0)),
                                   control = createControl(noiseLevel = "quiet",
                                                           useKKTSwindle = TRUE,
                                                           tuneSwindle = 100))    
    )    
    
    system.time(
    slowFit <- fitCyclopsModel(cyclopsData,
                               forceNewObject = TRUE, # Cold start for fair comparison
                               prior = createPrior("laplace", variance = 0.01, exclude = c(0)))
    )
    

    
    slowFit$timeFit
    fastFit$timeFit
    
    system.time(
        slowFit <- fitCyclopsModel(cyclopsData,
                                   forceNewObject = TRUE, # Cold start for fair comparison
                                   prior = createPrior("laplace", 
                                                       useCrossValidation = TRUE, exclude = c(0)),
                                   control = createControl(noiseLevel = "silent",
                                                           seed = seed,
                                                           cvType = "auto"))
    )
    
    system.time(
        fastFit <- fitCyclopsModel(cyclopsData,
                                   forceNewObject = TRUE, # Cold start for fair comparison
                                   prior = createPrior("laplace", 
                                                       useCrossValidation = TRUE, exclude = c(0)),
                                   control = createControl(noiseLevel = "silent",
                                                           seed = seed,
                                                           cvType = "auto",
                                                           useKKTSwindle = TRUE,
                                                           tuneSwindle = 100))    
    )
    
    library(microbenchmark)
    
#    start_profiler("samples.log")
    microbenchmark(
        fastFit <- fitCyclopsModel(cyclopsData,
                                   forceNewObject = TRUE, # Cold start for fair comparison
                                   prior = createPrior("laplace", 
                                                       useCrossValidation = TRUE, exclude = c(0)),
                                   control = createControl(noiseLevel = "silent",
                                                           seed = seed,
                                                           cvType = "auto",
                                                           useKKTSwindle = TRUE,
                                                           tuneSwindle = 100))   
    , times=5L)
#    stop_profiler()
    
    slowFit$timeFit
    fastFit$timeFit    
    
    
    expect_equal(coef(slowFit), coef(fastFit), tolerance = tolerance)
})