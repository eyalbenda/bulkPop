library(bulkPop)

test_that("Individuals are created correctly and can be accessed",{
  expect_identical(as.logical(newIndividual(N2xCB4856.genome,genotype=0,sex = 0)$hap1),
                       rep(F,length(N2xCB4856.genome$markerNames)))
  expect_identical(as.logical(newIndividual(N2xCB4856.genome,genotype=1,sex = 0)$hap1),
                   rep(T,length(N2xCB4856.genome$markerNames)))
  expect_identical(as.logical(newIndividual(N2xCB4856.genome,genotype=0,sex = 0)$hap2),
                   rep(F,sum(N2xCB4856.genome$markerChrom!="X")))
  expect_identical(as.logical(newIndividual(N2xCB4856.genome,genotype=1,sex = 0)$hap2),
                   rep(T,sum(N2xCB4856.genome$markerChrom!="X")))
  expect_error(newIndividual(N2xCB4856.genome,hap1 = rep(2,10),sex = 0))
  expect_error(newIndividual(N2xCB4856.genome,hap1 = rep(1,10),sex = 0))
})

