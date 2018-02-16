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

#' Test if x is in an individual object
#' @param x an object to test
#' @export
is.individual = function(x)
  inherits(x,"individual")


### Determine the sex of an individual
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

#' Generate a new individual
#' Flexible function to generate an individual based on a genome. haplotypes and sex can be random or specified.
#' @param genome a genome object
#' @param genotype the genotype of the individual. Relevant when pure individuals are desired. either 0, 1 or 2. 0: individual will be homozygote for 0. 1: individual will be homozygote for 1. 2: individual will be heterozygote genome wide
#' @param hap1 a list including a bit vector which will be the first haplotype of the individual
#' @param hap2 a list including a bit vector which will be the second haplotype of the individual
#' @param sex male or female. male is specified by "male", "m", "0" or 0. female is "female", "f", "1" or 1
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
    ind$hap2 = as.bit(ind$hap2[which(genome$markerChrom!="X")])
    ind$sex = as.bit(0)
  } else
  {
    ind$sex = as.bit(1)
  }
  return(ind)
}

#' Get the genotype of the individual at a variant or position
#' Either an exact chromosome and position or a variant name are required
#' see the function findProxyForPosition for an easy way to find closest variant to any position in the genome
#' @param individual an individual to query
#' @param genome a genome object corresponding to the individual
#' @param haplotype which haplotype to query
#' @param var variant name
#' @param chrom chromosome to query
#' @param pos position to query
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

#' Simulate a meiotic recombination in an individual
#' Workhorse function that takes an individual and performs recombinations based on the genetic map. shouldn't usually be called directly in most cases.
#' @param individual an individual to recombine
#' @param genome a genome object corresponding to the individual
#' @param Nrecombs number of recombinations per chromosome
#' @return a new recombined individual
#' @export
recombineIndividual = function(individual,genome,Nrecombs=1)
{
  if(!is.individual(individual))
    stop("individual given must be a proper individual")
  chromLast = as.matrix(unlist(values(genome$.chromRanges)))[2,]
  whichRecomb = chromLast<=length(individual$hap1)&chromLast<=length(individual$hap2)
  recombChroms = keys(genome$.chromRanges)[whichRecomb]
  newHap1 = individual$hap1
  newHap2 = individual$hap2
  for(chrom in recombChroms)
  {
    chromIndices = unlist(genome$.chromRanges[[chrom]],use.names = F)
    recombSpot = getAndUpdateRecombSpot(genome,chrom)
    recombIndices = sample(list(chromIndices[1]:recombSpot,recombSpot:chromIndices[2]),size = 1)[[1]]
    newHap1[recombIndices] = individual$hap2[recombIndices]
    newHap2[recombIndices] = individual$hap1[recombIndices]
    }
  newHaps = list(newHap1,newHap2)
  ord = sample(1:2,size=2)
  individual$hap1 = newHaps[[ord[1]]]
  individual$hap2 = newHaps[[ord[2]]]
  return(individual)
}


#' Cross two indiviiduals
#' Workhorse function to cross two individuals. Usually isn't called directly but useful for customized simulations
#' Takes two parents, recombines them, than randomly generates a progeny by a combination of one haplotype of each.
#' @param parent1 first parent. an individual object
#' @param parent2 second parent. an indiviudal object
#' @param genome a genome object, should correspond to both individuals
#' @param incompatibility a list of incompatibilities. Each incompatibility is a list with two elements: params and func. The func is the function that tests the recombined parents and the offspring and determines if it's dead or not. See the functions zeelPeelSelection and medeaSelection for example of functions included in the package.
#' @param sex what sex should the offspring have. either "random", or male (0,"0","male","m") or female (1,"1","female","f")
#' @param Nrecombs number of recombination per chromosome
#' @export
crossIndividuals = function(parent1,parent2,genome,incompatibility=NULL,sex = "random",Nrecombs = 1)
{
  if(as.logical(parent1$sex==parent2$sex)) stop("Error: same-sex mating isn't possible in this package")
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
  ord = sample(1:2)
  newHap1 = newHap2 = NULL
  assign(paste("newHap",ord[1],sep=""),recombParent1[[chosenFrom[1]]])
  assign(paste("newHap",ord[2],sep=""),recombParent2[[chosenFrom[2]]])
  if(!is.null(incompatibility))
  {
    for(inc in 1:length(incompatibility))
    {
      if(incompatibility[[inc]]$func(genome,recombParent1,recombParent2,
                                     sexes,chosenFrom,incompatibility[[inc]]$params))
        return("dead")
    }

  }
  return(newIndividual(genome = genome,hap1 = newHap1,hap2 = newHap2,sex=sex))
}


#' Perform a selfing
#' Hasn't been thoroughly tested, basically generates a male version of the same individuals (with random haplotype chosen to be shortened) and crosses them
#' @param parent the individual to self
#' @param genome a genome object. Corresponds to the individual
#' @param incompatibility a list of incompatibilities. Each incompatibility is a list with two elements: params and func. The func is the function that tests the recombined parents and the offspring and determines if it's dead or not. See the functions zeelPeelSelection and medeaSelection for example of functions included in the package.
#' @param sex the sex of the offspring
#' @param Nrecombs number of recombinations per chromosome
#' @export
selfHermaphrodite = function(parent,genome,incompatibility=NULL,sex = "female",Nrecombs = 1)
{
  if(as.logical(parent$sex)!=1) stop("Error: attempting to self a male")
  alleleAmale = newIndividual(genome = genome,hap1 = parent$hap1,
                              hap2 = parent$hap2,sex = "male")
  alleleBmale = newIndividual(genome = genome,hap1 = parent$hap2,
                              hap2 = parent$hap1,sex = "male")
  coin = sample(1:2)[1]
  if(coin==1)
  {
    return(crossIndividuals(parent,alleleAmale,genome = genome,incompatibility = incompatibility,sex=sex,Nrecombs = Nrecombs))
  } else
  {
    return(crossIndividuals(parent,alleleBmale,genome = genome,incompatibility = incompatibility,sex=sex,Nrecombs = Nrecombs))
  }
}
