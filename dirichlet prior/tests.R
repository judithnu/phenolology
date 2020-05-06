# tests for priorAnalysisGammaFunctions.R

library(testthat)
source('priorAnalysisGammaFunctions.R')

# rtnorm
## test that rtnorm returns returns a vector with the correct length
test_that("correct length", {
    expect_equal(length(rtnorm(n=10, mean=0, sd=1, min=-0.5, max=2)), 20)
})

## test that rtnorm doesn't return values outside of the range

# simulate_data
## test that simulate_data without groups returns a list of length 4 and with groups a list of length 6
## test that simulate_data chooses the right simulation file


# fit_gamma_model
## test that the right file is being chosen

# label_names
#

# posterior_differencer
