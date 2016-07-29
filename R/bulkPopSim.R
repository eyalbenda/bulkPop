#' @import qtl
#' @import bit
#' @export
c.individual = function(...)
{
  list(...)
}
#' @export
is.individual = function(x)
  inherits(x,"individual")

#' @export
newIndividual = function(genome,genotype = NULL, hap1 = NULL,hap2 = NULL,sex)
{
  if(!is.genome(genome))
    stop("function requires a genome object")
  if(!xor(is.null(genotype),(is.null(hap1)|is.null(hap2))))
    stop("Individual should either have a genotype(0 - hom, 1 - hom, 2 - het) or two haplotypes")
  if(is.null(genotype))
  {
    if(!is.bit(hap1))
      if(!all(hap1 %in% c(0,1)))
        stop("haplotypes must be vectors of 0 and 1 or bit objects")
    if(!is.bit(hap2))
      if(!all(hap2 %in% c(0,1)))
        stop("haplotypes must be vectors of 0 and 1 or bit objects")
  } else if (!genotype %in% c(0,1,2) | length(genotype)!=1)
  {
    stop("genotype should be either 0 for homozygote,1 for homozygot, or 2 for heterozygote")
  }
  if(!sex %in% c(0,1,"m","f","male","female",as.bit(0),as.bit(1)))
    stop("sex should be 0, 1, male, female, f or m")
  if(sex %in% c("male","m","0",as.bit(0)))
  {
    sex = as.bit(0)
  } else
  {
    sex = as.bit(1)
  }
  # if(!is.list(hap1)|!is.list(hap2)|!is.list(sex))
  #   stop("hap1, hap2 and sex must all be lists")
  if(is.null(genotype))
  {
    return(structure(list(hap1 = as.bit(hap1), hap2 = as.bit(hap2),sex = sex),class = "individual"))
  } else if (genotype==0)
  {
    ind = structure(list(hap1 = as.bit(rep(0,length(genome$markerChrom))), hap2 = as.bit(rep(0,length(genome$markerChrom)),sex = sex)),class = "individual")
  } else if (genotype==1)
  {
    ind = structure(list(hap1 = as.bit(rep(1,length(genome$markerChrom))), hap2 = as.bit(rep(1,length(genome$markerChrom)),sex = sex)),class = "individual")
  } else
  {
    hap1 = as.bit(rep(0,length(genome$markerChrom)))
    hap2 = as.bit(rep(1,length(genome$markerChrom)))
    ind = sample(c(structure(list(hap1 = hap1, hap2 = hap2,sex = sex),class = "individual"),
                   structure(list(hap1 = hap2, hap2 = hap1,sex = sex),class = "individual")),size=1)[[1]]
  }
  if(!as.logical.bit(sex))
  {
    ind$hap2 = as.bit(rev(rev(ind$hap2)[-c(1:sum(genome$markerChrom=="X"))]))
    ind$sex = as.bit(0)
  } else
  {
    ind$sex = as.bit(1)
  }
  return(ind)
}

#' @export
getGenotypes = function(individual,genome,haplotype,names = NULL,chromosomes=NULL,positions =NULL)
{
  if(!is.individual(individual))
    stop("individual given must be a proper individual")
  if(is.null(names) & (is.null(chromosomes)|is.null(positions))) stop("either give a name of the marker, or BOTH chromosome and position")
  marker = NULL
  if(!is.null(names))
  {
    marker = match(names,genome$markerNames)
    if(any(is.na(marker))) stop("Invalid marker name given")
  } else
  {
    marker = intersect(match(chromosome,genome$markerChrom),match(cbind(positions,genome$markerPos)))
    if(any(is.na(marker))) stop("Invalid marker positions given")
  }
  genotypeVector = individual[[paste("hap",haplotype,sep="")]]
  if(!all(marker %in% 1:length(genotypeVector)))
  {
    return(NULL)
  } else
  {
    return(individual[[paste("hap",haplotype,sep="")]][marker])
  }
}

#' @export
recombineIndividual = function(individual,genome,Nrecombs=1)
{
  if(!is.individual(individual))
    stop("individual given must be a proper individual")
  if(as.logical.bit(individual$sex))
  {
    if(length(individual$hap1)!=length(individual$hap2))
      stop("Specified female, but different haplotype lengths")
    recombChroms = unique(genome$markerChrom)
  } else
  {
    if(length(individual$hap1)==length(individual$hap2))
      stop("Specified male, but different haplotype lengths")
    recombChroms = rev(rev(unique(genome$markerChrom))[-1])
  }
  for(chrom in recombChroms)
  {
    chromIndices = which(genome$markerChrom==chrom)
    recombSpot = sample(chromIndices[-length(chromIndices)],size = Nrecombs,prob = genome$recombProbs[chromIndices[-length(chromIndices)]])
    newHap1 = as.bit(c(individual$hap1[min(chromIndices):recombSpot],individual$hap2[(recombSpot+1):max(chromIndices)]))
    newHap2 = as.bit(c(individual$hap2[min(chromIndices):recombSpot],individual$hap1[(recombSpot+1):max(chromIndices)]))
    newHaps = list(newHap1,newHap2)
    ord = sample(1:2,size=2)
    individual$hap1[chromIndices] = newHaps[[ord[1]]]
    individual$hap2[chromIndices] = newHaps[[ord[2]]]
  }
  return(individual)
}


#' @export
crossIndividuals = function(parent1,parent2,genome,zeelPeelSelection=T,sex = "random",Nrecombs = 1)
{
  if(as.logical(parent1$sex==parent2$sex)) stop("Error: same-sex coupling doesn't exist for c. elegans")
  recombParent1 = recombineIndividual(parent1,genome,Nrecombs)
  recombParent2 = recombineIndividual(parent2,genome,Nrecombs)
  if(sex == "random")
  {
    sex = sample(c("male","female"),size=1)
  }
  sexes = as.logical.bit(c(recombParent1$sex,recombParent2$sex))
  r1l = sapply(recombParent1[1:2],length)
  r2l = sapply(recombParent2[1:2],length)
  if(sex == "female")
  {
    chosenFrom = c(which.max(r1l),which.max(r2l))
  } else
  {
    chosenFrom = c(which.min(r1l),which.min(r2l))
  }
  newHap1 = recombParent1[[chosenFrom[1]]]
  newHap2 = recombParent2[[chosenFrom[2]]]
  if(zeelPeelSelection)
    if(zeelPeel(genome,parent1,parent2,sexes,chosenFrom))
      return("dead")
  return(newIndividual(genome = genome,hap1 = newHap1,hap2 = newHap2,sex=sex))
}

zeelPeel = function(genome,parent1,parent2,sexes,chosenFrom)
{
  parentsCombined = list(parent1,parent2)
  femPar = which(sexes)
  if(!exists("zeelPeelMarker")|!exists("zeelPeelKill")) stop(c("zeelPeelMarker or zeelPeelKill undefined, can't perform selection"))
  zeelPeel = 0
  if(as.logical(getGenotypes(parentsCombined[[femPar]],genome,haplotype = chosenFrom[[femPar]][1], name = zeelPeelMarker)==zeelPeelKill))
    zeelPeel = zeelPeel+1
  if(as.logical(xor(getGenotypes(parentsCombined[[-femPar]],genome,1,name = zeelPeelMarker)==zeelPeelKill,
                    getGenotypes(parentsCombined[[-femPar]],genome,2,name = zeelPeelMarker)==zeelPeelKill)))
    {
      if(as.logical(getGenotypes(parentsCombined[[-femPar]],genome,chosenFrom[[-femPar]],name = zeelPeelMarker)==zeelPeelKill))
        zeelPeel = zeelPeel + 1
    }
  return(zeelPeel == 2)
}


#' @export
selfHermaphrodite = function(parent,genome,zeelPeelSelection=T,sex = "female",Nrecombs = 1)
{
  if(as.logical(parent$sex)!=1) stop("Error: attempting to self a male")
  alleleAmale = newIndividual(genome = genome,hap1 = parent$hap1[which(genome$markerChrom!="X")],
                              hap2 = parent$hap2,sex = "male")
  alleleBmale = newIndividual(genome = genome,hap1 = parent$hap2[which(genome$markerChrom!="X")],
                              hap2 = parent$hap1,sex = "male")
  coin = sample(1:2)[1]
  if(coin==1)
  {
    return(crossIndividuals(parent,alleleAmale,genome = genome,zeelPeelSelection = zeelPeelSelection,sex=sex,Nrecombs = Nrecombs))
  } else
  {
    return(crossIndividuals(parent,alleleBmale,genome = genome,zeelPeelSelection = zeelPeelSelection,sex=sex,Nrecombs = Nrecombs))
  }
}

#' @export
getIndividualMatingProbabilities = function(individuals,lociList,genome)
{
  if(ncol(lociList)==0) return(rep(1,length(individuals)))
  lociList = as.matrix(lociList)
  indProbs = rep(1,length(individuals))
  for(ind in 1:length(individuals))
    {
      for(hap in 1:2)
      {
        for(locus in 1:nrow(lociList))
        {
        curLocus = as.character(lociList[locus,])
        curGenotype = getGenotypes(individual = individuals[[ind]],genome = genome, haplotype = hap,names = curLocus[2])
        if(!is.null(curGenotype))
          if(curGenotype==as.numeric(curLocus[3]))
            indProbs[ind] = indProbs[ind] * fitnessLoci[locus,4]
        }
      }
    }
  indProbs
}

