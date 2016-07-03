#' @import qtl
#' @import bit


.newIndividual = setClass(Class = "individual",slots = c(hap1 = "list", hap2 = "list",sex = "list"))


#' @export
newIndividual = function(hap1,hap2,sex)
{
  if(!is.list(hap1)|!is.list(hap2)|!is.list(sex))
    stop("hap1, hap2 and sex must all be lists")
  .newIndividual(hap1 =  hap1,hap2 = hap2,sex = sex)
}

#' @export
getGenotypes = function(individual,genome,haplotype,names = NULL,chromosomes=NULL,positions =NULL)
{
  if(is.null(names) & (is.null(chromosomes)|is.null(positions))) stop("either give a name of the marker, or BOTH chromosome and position")
  marker = NULL
  if(!is.null(names))
  {
    marker = match(names,genome@markerNames)
    if(any(is.na(marker))) stop("Invalid marker name given")
  } else
  {
    marker = intersect(match(chromosome,genome@markerChrom),match(cbind(positions,genome@markerPos)))
    if(any(is.na(marker))) stop("Invalid marker positions given")
  }
  genotypeVector = slot(individual,paste("hap",haplotype,sep=""))[[1]]
  if(!all(marker %in% 1:length(genotypeVector)))
  {
    return(NULL)
  } else
  {
    return(slot(individual,paste("hap",haplotype,sep=""))[[1]][marker])
  }
}

#' @export
recombineIndividual = function(individual,genome,Nrecombs=1)
{
  if(as.logical(individual@sex[[1]]))
  {
    if(length(individual@hap1[[1]])!=length(individual@hap2[[1]]))
      stop("Specified female, but different haplotype lengths")
    recombChroms = unique(genome@markerChrom)
  } else
  {
    if(length(individual@hap1[[1]])==length(individual@hap2[[1]]))
      stop("Specified male, but different haplotype lengths")
    recombChroms = rev(rev(unique(genome@markerChrom))[-1])
  }
  for(chrom in recombChroms)
  {
    chromIndices = which(genome@markerChrom==chrom)
    recombSpot = sample(chromIndices[-length(chromIndices)],size = Nrecombs,prob = genome@recombProbs[chromIndices[-length(chromIndices)]])
    newHap1 = as.bit(c(individual@hap1[[1]][min(chromIndices):recombSpot],individual@hap2[[1]][(recombSpot+1):max(chromIndices)]))
    newHap2 = as.bit(c(individual@hap2[[1]][min(chromIndices):recombSpot],individual@hap1[[1]][(recombSpot+1):max(chromIndices)]))
    individual@hap1[[1]][chromIndices] = newHap1
    individual@hap2[[1]][chromIndices] = newHap2
  }
  return(individual)
}


#' @export
crossIndividuals = function(parent1,parent2,genome,zeelPeelSelection=T,sex = "random",Nrecombs = 1)
{
  if(as.logical.bit(!xor(parent1@sex[[1]],parent2@sex[[1]]))) stop("Error: same-sex coupling doesn't exist for c. elegans")
  recombParent1 = recombineIndividual(parent1,genome,Nrecombs)
  recombParent2 = recombineIndividual(parent2,genome,Nrecombs)
  newHap1 = bit(length(genome@markerChrom))
  newHap2 = bit(length(genome@markerChrom))
  chosenFrom1 = NULL
  chosenFrom2 = NULL
  autosomes = rev(rev(unique(genome@markerChrom))[-1])
  for(chrom in autosomes)
  {
    chromIndices = which(genome@markerChrom==chrom)
    from1 = sample(1:2,size = 1)
    from2 = sample(1:2,size = 1)
    hapOrder =sample(1:2, size = 1)

    if(hapOrder == 1)
    {
      newHap1[chromIndices] = slot(recombParent1,paste("hap",from1,sep=""))[[1]][chromIndices]
      newHap2[chromIndices] = slot(recombParent2,paste("hap",from2,sep=""))[[1]][chromIndices]
    } else
    {
      newHap1[chromIndices] = slot(recombParent2,paste("hap",from2,sep=""))[[1]][chromIndices]
      newHap2[chromIndices] = slot(recombParent1,paste("hap",from1,sep=""))[[1]][chromIndices]
    }


    chosenFrom1 = c(chosenFrom1,from1)
    chosenFrom2 = c(chosenFrom2,from2)
  }
  chosenFrom = list(chosenFrom1,chosenFrom2)
  sexes = as.logical(c(parent1@sex[[1]],parent2@sex[[1]]))
  parentsCombined = list(recombParent1,recombParent2)
  femPar = which(sexes)
  if(zeelPeelSelection)
  {
    if(!exists("zeelPeelMarker")|!exists("zeelPeelKill")) stop(c("zeelPeelMarker or zeelPeelKill undefined, can't perform selection"))
    zeelPeel = 0
    if(as.logical(getGenotypes(parentsCombined[[femPar]],genome,haplotype = chosenFrom[[femPar]][1], name = zeelPeelMarker)==zeelPeelKill))
      zeelPeel = zeelPeel+1
    if(as.logical(xor(getGenotypes(parentsCombined[[-femPar]],genome,1,name = zeelPeelMarker)==zeelPeelKill,
           getGenotypes(parentsCombined[[-femPar]],genome,2,name = zeelPeelMarker)==zeelPeelKill)))
    {
      if(as.logical(getGenotypes(parentsCombined[[-femPar]],genome,chosenFrom[[-femPar]][1],name = zeelPeelMarker)==zeelPeelKill))
        zeelPeel = zeelPeel + 1
    }
    if(zeelPeel == 2) return("dead")
  }
  if(sex == "random")
  {
    newSex = sample(c("male","female"),size = 1)
  } else
  {
    newSex = sex
  }
  newSexBit = ifelse(newSex=="female",1,0)
  chooseFemHap = sample(c(1,2),size=1)
  sexChromIndices = which(genome@markerChrom=="X")
  XchromMale = list(getGenotypes(parentsCombined[[-femPar]],genome,1,names = genome@markerNames[genome@markerChrom=="X"]),
              getGenotypes(parentsCombined[[-femPar]],genome,2,names = genome@markerNames[genome@markerChrom=="X"]))
  newHap1[sexChromIndices] = getGenotypes(parentsCombined[[femPar]],genome,chooseFemHap,names = genome@markerNames[genome@markerChrom=="X"])
  if(newSex == "female")
  {
    newHap2[sexChromIndices] = as.bit(XchromMale[!sapply(XchromMale,is.null)][[1]])
  } else
  {
    newHap2 = newHap2[-sexChromIndices]
  }
  offspring = sample(c(newIndividual(hap1=list(as.bit(newHap1)),hap2=list(as.bit(newHap2)),sex=list(newSexBit)),
                       newIndividual(hap1=list(as.bit(newHap2)),hap2=list(as.bit(newHap1)),sex=list(newSexBit))),
                     size=1)[[1]]
  return(offspring)
}

#' @export
selfHermaphrodite = function(parent,genome,zeelPeelSelection=T,sex = "female",Nrecombs = 1)
{
  if(parent@sex!=1) stop("Error: attempting to self a male")
  alleleAmale = newIndividual(hap1 = list(as.bit(parent@hap1[[1]][which(genome@markerChrom!="X")])),
                              hap2 = list(as.bit(parent@hap2[[1]])),sex = list(as.bit(0)))
  alleleBmale = newIndividual(hap1 = list(as.bit(parent@hap2[[1]][which(genome@markerChrom!="X")])),
                              hap2 = list(as.bit(parent@hap1[[1]])),sex = list(as.bit(0)))
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

