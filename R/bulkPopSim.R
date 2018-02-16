#' @import qtl
#' @import bit
#' @import nnet
#' @import parallel
#' @importFrom stats rnorm


## Internal function to generate offspring from mating a female to a pool of males
## x is the index of the female chosen for mating, a single male is chosen at random
## then offPerMating offsprings are generated, using the parameters from genome,Nrecombs and incompatibility
matePop = function(x,females,males,offPerMating,genome,Nrecombs,incompatibility)
{
  curMale = sample(males,size = 1)[[1]]
  curFemale = females[[x]]
  numOffspring = max(1,round(rnorm(1,offPerMating,1)))
  newPop = list()
  for(k in 1:offPerMating)
  {
    newInd = crossIndividuals(parent1 = curMale, parent2 = curFemale,genome = genome,Nrecombs = Nrecombs,incompatibility = incompatibility,sex = "random")
    if(!is.character(newInd))
      newPop = c(newPop,list(newInd))
  }
  newPop
}

#' Mate a population of males to a population of females.
#' @param males a list of male individuals
#' @param females a list of female individuals
#' @param genome A genome object generated with the package
#' @param offPerMating how many offspring does each female generate
#' @param ncores Number of cores to be used. If more than one core is specified, a cluster object will be created and used.
#' @param Nrecombs the number of recombinations per chromosome per cross (defaults to 1)
#' @param incompatibility a list of incompatibilities. Each incompatibility is a list with two elements: params and func. The func is the function that tests the recombined parents and the offspring and determines if it's dead or not. See the functions zeelPeelSelection and medeaSelection for example of functions included in the package.
#' @export
performMating = function(males,females,genome,offPerMating,ncores = 1,Nrecombs = 1,incompatibility = NULL)
{
  if(ncores>1)
  {
    cl = initiateCluster(ncores = ncores)
    clusterCall(cl,function(x)library(bulkPop))
    clusterExport(cl,c("genome","offPerMating","males","females","incompatibility"),
                  envir = environment())
    bigPop = parLapply(cl,1:length(females),matePop,females,males,offPerMating,genome,Nrecombs,incompatibility)
    stopCluster(cl)
  } else
  {
    bigPop = lapply(1:length(females),matePop,females,males,offPerMating,genome,Nrecombs,incompatibility=incompatibility)
  }
  unlist(bigPop,recursive=F,use.names=F)
}

#' Use the fitness loci to choose the submating population out of the pool of males or females.
#' @param genome A genome object generated with the package
#' @param pop A list of individuals of males and females to pick from
#' @param matingFrac the fraction of the relevant sex population which will be chosen. between 0 and 1 (1 is no selection)
#' @param fitnessLoci A data frame with the following columns:
#' chrom: chromosome
#' marker: name of the marker to use (see findProxyForPosition for easy way to get marker for locus of interest)
#' direction: (0,1) which allele confers an advantage?
#' effect: the effect size of the locus. How much advantage does it confer?
#' @param sex choose a population of males or females?
#' @export
getMatingSubpopulation = function(genome,pop,matingFrac,fitnessLoci = NULL,sex=c("male","female"))
{
  newPop = list()
  sexPop = sapply(pop,function(x)x$sex)
  if(sex == "male")
  {
    choos = which(sexPop==0)
  } else
  {
    choos = which(sexPop==1)
  }
  if(is.null(fitnessLoci))
  {
    mateProb = rep(1,length(choos))
  } else if(any(fitnessLoci$sex==sex))
  {
    mateProb = getIndividualMatingProbabilities(pop[choos],fitnessLoci[fitnessLoci$sex==sex,],genome)
  } else
  {
    mateProb = rep(1,length(choos))
  }
  choos = sample(choos,size=length(choos)*(matingFrac),prob = mateProb)
  return(pop[choos])
}

#' Return a relative measure of mating probability for a list of individuals
#' This function shouldn't usually be called directly, use getMatingSubpopulation instead
#' @param individuals a list of individuals
#' @param lociList see fitnessLoci parameter in getMatingSubpopulation
#' @param genome a genome object
#' @export
getIndividualMatingProbabilities = function(individuals,lociList,genome)
{
  if(ncol(lociList)==0) return(rep(1,length(individuals)))
  indProbs = rep(1,length(individuals))/length(individuals)
  effect = as.numeric(lociList$effect)
  marker = as.character(lociList$marker)
  direction = lociList$direction
  for(ind in 1:length(individuals))
    {
      for(hap in 1:2)
      {
        for(locus in 1:nrow(lociList))
        {
        curGenotype = getGenotypes(individual = individuals[[ind]],genome = genome, haplotype = hap,var = marker[locus])
        if(curGenotype==as.numeric(direction[locus]))
          indProbs[ind] = indProbs[ind] * effect[locus]
        }
      }
    }
  indProbs
}

