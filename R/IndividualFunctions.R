#' @import qtl
#' @import bit
#' @import nnet
#' @import hash
#' @importFrom stats approx
#' @export
c.individual = function(...)
{
  list(...)
}
#' @export
is.individual = function(x)
  inherits(x,"individual")


sexParse = function(sex)
{
  if(!sex %in% c(0,1,"m","f","male","female",as.bit(0),as.bit(1)))
    stop("sex should be 0, 1, male, female, f or m")
  if(sex %in% c("male","m","0",as.bit(0)))
  {
    sex = as.bit(0)
  } else
  {
    sex = as.bit(1)
  }
  sex
}

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
  sex = sexParse(sex)
  # if(!is.list(hap1)|!is.list(hap2)|!is.list(sex))
  #   stop("hap1, hap2 and sex must all be lists")
  if(is.null(genotype))
  {
    return(structure(list(hap1 = as.bit(hap1), hap2 = as.bit(hap2),sex = sex),class="individual"))
  } else if (genotype==0)
  {
    ind = structure(list(hap1 = as.bit(rep(0,length(genome$markerChrom))), hap2 = as.bit(rep(0,length(genome$markerChrom))),sex = sex),class = "individual")
  } else if (genotype==1)
  {
    ind = structure(list(hap1 = as.bit(rep(1,length(genome$markerChrom))), hap2 = as.bit(rep(1,length(genome$markerChrom))),sex = sex),class = "individual")
  } else
  {
    hap1 = as.bit(rep(0,length(genome$markerChrom)))
    hap2 = as.bit(rep(1,length(genome$markerChrom)))
    ind = sample(c(structure(list(hap1 = hap1, hap2 = hap2,sex = sex),class = "individual"),
                   structure(list(hap1 = hap2, hap2 = hap1,sex = sex),class = "individual")),size=1)[[1]]
  }
  if(!as.logical.bit(sex))
  {
    ind$hap2 = ind$hap2[-which(genome$markerChrom=="X")]
    ind$sex = as.bit(0)
  } else
  {
    ind$sex = as.bit(1)
  }
  return(ind)
}

#' @export
getGenotypes = function(individual,genome,haplotype = c(1,2),var = NULL,chrom=NULL,pos =NULL)
{
  if(!is.individual(individual))
    stop("individual given must be a proper individual")
  if(is.null(var) & (is.null(chrom)|is.null(pos))) stop("either give a name of the variant, or BOTH chromosome and position")
  if(!is.null(var))
  {
    lookup = var
    if(any(!has.key(key = lookup,hash = genome$.NameToPos))) stop("One or more vars don't exist in genome")
    choos = sapply(lookup,function(x)genome$.NameToPos[[x]])
  } else
  {
    lookup = paste(chrom,pos,sep="_")
    if(any(!has.key(key =lookup ,hash = genome$.LocToPos))) stop("One or more positions aren't associated with a variant")
    choos = sapply(lookup,function(x)genome$.LocToPos[[x]])
  }
  genotypeVector = individual[[paste("hap",haplotype,sep="")]]
  return(individual[[paste("hap",haplotype,sep="")]][choos])
}

#' @export
recombineIndividual = function(individual,genome,Nrecombs=1)
{
  if(!is.individual(individual))
    stop("individual given must be a proper individual")
  chromLast = as.matrix(unlist(values(genome$.chromRanges)))[2,]
  recombChroms = keys(genome$.chromRanges)[chromLast<=length(individual$hap1)&chromLast<=length(individual$hap2)]
  newHap1 = individual$hap1
  newHap2 = individual$hap2
  for(chrom in recombChroms)
  {
    chromIndices = unlist(genome$.chromRanges[[chrom]],use.names = F)
    recombSpot = getAndUpdateRecombSpot(genome,chrom)
    newHap1[chromIndices[1]:recombSpot] = individual$hap2[chromIndices[1]:recombSpot]
    newHap2[chromIndices[1]:recombSpot] = individual$hap1[chromIndices[1]:recombSpot]
    }
  newHaps = list(newHap1,newHap2)
  ord = sample(1:2,size=2)
  individual$hap1 = newHaps[[ord[1]]]
  individual$hap2 = newHaps[[ord[2]]]
  return(individual)
}


#' @export
crossIndividuals = function(parent1,parent2,genome,zeelPeel=NULL,sex = "random",Nrecombs = 1)
{
  if(as.logical(parent1$sex==parent2$sex)) stop("Error: same-sex coupling doesn't exist for c. elegans")
  recombParent1 = recombineIndividual(parent1,genome,Nrecombs)
  recombParent2 = recombineIndividual(parent2,genome,Nrecombs)
  if(sex == "random")
  {
    sex = sample(c("male","female"),size=1)
  }
  sex = sexParse(sex)

  sexes = as.logical.bit(c(recombParent1$sex,recombParent2$sex))
  r1l = sapply(recombParent1[1:2],length)
  r2l = sapply(recombParent2[1:2],length)
  if(as.logical.bit(sex))
  {

    chosenFrom = c(which.is.max(r1l),
                   which.is.max(r2l))
  } else
  {
    chosenFrom = c(which.is.max(-r1l),
                   which.is.max(-r2l))
  }
  newHap1 = recombParent1[[chosenFrom[1]]]
  newHap2 = recombParent2[[chosenFrom[2]]]
  if(!is.null(zeelPeel))
    if(zeelPeelSelection(genome,recombParent1,recombParent2,sexes,chosenFrom,zeelPeel))
      return("dead")
  return(newIndividual(genome = genome,hap1 = newHap1,hap2 = newHap2,sex=sex))
}

zeelPeelSelection = function(genome,parent1,parent2,sexes,chosenFrom,zeelPeel)
{
  parentsCombined = list(parent1,parent2)
  femPar = which(sexes)
  if(is.null(zeelPeel$marker)|is.null(zeelPeel$kill)) stop(c("zeelPeel must be a list with elements marker and kill"))
  zp = 0
  if(as.logical(getGenotypes(parentsCombined[[femPar]],genome,haplotype = chosenFrom[[femPar]], var = zeelPeel$marker)==zeelPeel$kill))
    zp = zp+1
  if(as.logical(xor(getGenotypes(parentsCombined[[-femPar]],genome,1,var = zeelPeel$marker)==zeelPeel$kill,
                    getGenotypes(parentsCombined[[-femPar]],genome,2,var = zeelPeel$marker)==zeelPeel$kill)))
  {
    if(as.logical.bit(getGenotypes(parentsCombined[[-femPar]],genome,chosenFrom[[-femPar]],var = zeelPeel$marker)==zeelPeel$kill))
      zp = zp + 1
  }
  return(zp == 2)
}


#' @export
selfHermaphrodite = function(parent,genome,zeelPeel=NULL,sex = "female",Nrecombs = 1)
{
  if(as.logical(parent$sex)!=1) stop("Error: attempting to self a male")
  alleleAmale = newIndividual(genome = genome,hap1 = parent$hap1,
                              hap2 = parent$hap2,sex = "male")
  alleleBmale = newIndividual(genome = genome,hap1 = parent$hap2,
                              hap2 = parent$hap1,sex = "male")
  coin = sample(1:2)[1]
  if(coin==1)
  {
    return(crossIndividuals(parent,alleleAmale,genome = genome,zeelPeel = zeelPeel,sex=sex,Nrecombs = Nrecombs))
  } else
  {
    return(crossIndividuals(parent,alleleBmale,genome = genome,zeelPeel = zeelPeel,sex=sex,Nrecombs = Nrecombs))
  }
}
