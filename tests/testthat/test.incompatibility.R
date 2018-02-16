require(bulkPop)

N2male = newIndividual(N2xCB4856.genome,genotype = 1,sex = 0)
N2fem = newIndividual(N2xCB4856.genome,genotype = 1,sex = 1)
Hmale = newIndividual(N2xCB4856.genome,genotype = 0,sex = 0)
Hfem = newIndividual(N2xCB4856.genome,genotype = 0,sex = 1)
F1male = newIndividual(N2xCB4856.genome,hap1 = N2fem$hap1,hap2 = Hmale$hap2,sex = "male")
F1fem = newIndividual(N2xCB4856.genome,hap1 = N2fem$hap1,hap2 = Hfem$hap2,sex = "female")

test_that("zeel-1/peel-1 selection is working correctly",{
expect_false(zeelPeelSelection(N2xCB4856.genome,N2male,N2fem,c(F,T),c(1,1),list(marker=N2xCB4856.genome$markerNames[10],kill=0)))
expect_false(zeelPeelSelection(N2xCB4856.genome,Hmale,N2fem,c(F,T),c(1,1),list(marker=N2xCB4856.genome$markerNames[10],kill=0)))
expect_false(zeelPeelSelection(N2xCB4856.genome,N2male,Hfem,c(F,T),c(1,1),list(marker=N2xCB4856.genome$markerNames[10],kill=0)))
expect_false(zeelPeelSelection(N2xCB4856.genome,Hmale,Hfem,c(F,T),c(1,1),list(marker=N2xCB4856.genome$markerNames[10],kill=0)))
expect_false(zeelPeelSelection(N2xCB4856.genome,F1male,Hfem,c(F,T),c(1,1),list(marker=N2xCB4856.genome$markerNames[10],kill=0)))
expect_false(zeelPeelSelection(N2xCB4856.genome,F1male,N2fem,c(F,T),c(1,1),list(marker=N2xCB4856.genome$markerNames[10],kill=0)))
expect_false(zeelPeelSelection(N2xCB4856.genome,F1male,N2fem,c(F,T),c(2,1),list(marker=N2xCB4856.genome$markerNames[10],kill=0)))
expect_true(zeelPeelSelection(N2xCB4856.genome,F1male,Hfem,c(F,T),c(2,1),list(marker=N2xCB4856.genome$markerNames[10],kill=0)))
expect_true(zeelPeelSelection(N2xCB4856.genome,F1male,F1fem,c(F,T),c(2,2),list(marker=N2xCB4856.genome$markerNames[10],kill=0)))
})

medea = list(markers=list("toxin"=N2xCB4856.genome$markerNames[10],"antidote"=N2xCB4856.genome$markerNames[10]),kill=0)

test_that("medea selection is working correctly",{
  expect_false(medeaSelection(N2xCB4856.genome,N2male,N2fem,c(F,T),c(1,1),medea))
  expect_false(medeaSelection(N2xCB4856.genome,Hmale,N2fem,c(F,T),c(1,1),medea))
  expect_false(medeaSelection(N2xCB4856.genome,N2male,Hfem,c(F,T),c(1,1),medea))
  expect_false(medeaSelection(N2xCB4856.genome,Hmale,Hfem,c(F,T),c(1,1),medea))
  expect_false(medeaSelection(N2xCB4856.genome,F1male,Hfem,c(F,T),c(1,1),medea))
  expect_false(medeaSelection(N2xCB4856.genome,F1male,N2fem,c(F,T),c(1,1),medea))
  expect_false(medeaSelection(N2xCB4856.genome,F1male,N2fem,c(F,T),c(2,1),medea))
  expect_false(medeaSelection(N2xCB4856.genome,F1male,Hfem,c(F,T),c(2,1),medea))
  expect_false(medeaSelection(N2xCB4856.genome,Hmale,F1fem,c(F,T),c(2,1),medea))
  expect_false(medeaSelection(N2xCB4856.genome,Hmale,F1fem,c(F,T),c(1,1),medea))

  expect_true(medeaSelection(N2xCB4856.genome,Hmale,F1fem,c(F,T),c(2,2),medea))
  expect_true(medeaSelection(N2xCB4856.genome,Hmale,F1fem,c(F,T),c(1,2),medea))
  expect_true(medeaSelection(N2xCB4856.genome,F1male,F1fem,c(F,T),c(2,2),medea))
})
