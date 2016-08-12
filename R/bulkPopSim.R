#' @import qtl
#' @import bit
#' @import nnet
#' @import parallel
#' @importFrom stats rnorm
#' @export

#' @export
performMating = function(males,females,genome,offPerMating,ncores = 1,...)
{
  cl = initiateCluster(ncores = ncores)
  if(!exists("zeelPeel"))
    zeelPeel = NULL
  clusterCall(cl,function(x)library(bulkPop))
  clusterExport(cl,c("genome","offPerMating","males","females","zeelPeel"),
                envir = environment())
  bigPop = parLapply(cl,1:length(females),function(x)
  {
    curMale = sample(males,size = 1)[[1]]
    curFemale = females[[x]]
    numOffspring = max(1,round(rnorm(1,offPerMating,1)))
    newPop = list()
    for(k in 1:offPerMating)
    {
      newInd = crossIndividuals(parent1 = curMale, parent2 = curFemale,genome = genome,Nrecombs=1,zeelPeel = zeelPeel)
      if(!is.character(newInd))
        newPop = c(newPop,list(newInd))
    }
    newPop
  })
  stopCluster(cl)
  unlist(bigPop,recursive=F,use.names=F)
}

#' @export
getMatingSubpopulation = function(genome,pop,matingFrac,fitnessLoci,sex=c("male","female"))
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
  if(any(fitnessLoci$sex==sex))
  {
    mateProb = getIndividualMatingProbabilities(pop[choos],fitnessLoci[fitnessLoci$sex==sex,],genome)
  } else
  {
    mateProb = rep(1,length(choos))
  }
  choos = sample(choos,size=length(choos)*(matingFrac),prob = mateProb)
  return(pop[choos])
}

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

