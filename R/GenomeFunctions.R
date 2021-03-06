
#' Test if a genome
#' @param x an object to test
#' @export
is.genome = function(x)
{
  x$class == "genome"
}

#' Create a new genome object - object including all markers and the recombination rates between them
#' @param name Name of the genome object
#' @param Nchrom Number of chromosomes
#' @param markerNames A vector of the name of the markers
#' @param markerChrom A vector of the chromosome of each marker
#' @param markerPos A vector of the position (in basepairs) of each marker
#' @param recombProbs A vector of the recombination probability between each two successive markers (can be generated from a genetic map using recombProbFromGeneticDistance)
#' @param map The genetic distance between two successive genetic markers
#' @param Nrecs Number of recombination positions to simulate when generating the genome object (recombination spots in crosses will be cycled through this object to save time. Defaults to 1,000,000)
#' @export
newGenome = function(name,Nchrom,markerNames,markerChrom,markerPos,recombProbs, map,Nrecs = 1e6)
{
  lengths = range(length(markerNames),length(markerChrom),length(markerPos),length(recombProbs),length(map))
  if(min(lengths)!=max(lengths)) stop("Names, chromosomes positions, recombination probabilities and map must all be same length")
  .NameToPos = hash(keys = markerNames, values = 1:length(markerNames))
  .LocToPos = hash(keys = paste(markerChrom,markerPos,sep="_"),values = 1:length(markerNames))
  .chromRanges = hash(keys = unique(markerChrom),values = tapply(1:length(markerChrom),markerChrom,function(x)c(min(x),max(x))))
  .recombPos = hash(keys = unique(markerChrom),values = lapply(unique(markerChrom),function(x)
  {
    chromIndices = unlist(.chromRanges[[x]],use.names = F);
    sample(chromIndices[1]:(chromIndices[2]-1),size=Nrecs,replace=T,
           prob=recombProbs[chromIndices[1]:(chromIndices[2]-1)])}))

  .recombInd = hash(keys = unique(markerChrom),values = lapply(unique(markerChrom),function(x)1))
  as.environment(list(name = name,
                      Nchrom = Nchrom,
                      markerNames = markerNames,
                      markerChrom = markerChrom,
                      markerPos = markerPos,
                      recombProbs = recombProbs,
                      map = map,
                      .NameToPos = .NameToPos,
                      .LocToPos = .LocToPos,
                      .chromRanges = .chromRanges,
                      .recombPos = .recombPos,
                      .Nrecs = Nrecs,
                      .recombInd = .recombInd,
                      class="genome"))
}

#' Convert genetic distance to probability of recombination
#' @param distances A vector of genetic distances between the markers
#' @export
#' @note This function assumes a single recombination event per chromosome!!!
recombProbFromGeneticDistance = function(distances){
  distChrom = distances[-1] - distances[-length(distances)]
  if(any(order(distances)!=1:length(distances)))
    stop("Markers must be ordered by their location in the genome")
  finalProbs = c(0,distChrom / sum(distChrom))
  names(finalProbs) = names(distances)
  finalProbs
}


#' The function extends a genome to more variants, by imputing the genetic map in between
#' @param genome a genome object
#' @param variantNames vector of variant names. If not given will be generated as chrom:pos
#' @param variantChrom vector of variant chromsomes. Names must match genome object
#' @param variantPos vector of variant positions.
#' @param matchNames boolean: If true, the function will look for the new variants by position in the old genome and carry over the names to the new genome.
#' @export
expandGenomeToVariantlist = function(genome,variantNames = NULL,variantChrom,variantPos,matchNames = F)
{
  if(!is.character(variantChrom))
    stop("variantChrom must be character vector!")
  if(!is.numeric(variantPos))
    stop("variantPos must be a numeric vector!")
  if(any(!unique(variantChrom) %in% unique(genome$markerChrom)))
    stop("variants in chromosomes that don't exist in supplied genome!")
  if(is.null(variantNames))
    variantNames = paste(variantChrom,variantPos,sep="_")
  if(matchNames)
  {
    variantComb = paste(variantChrom,variantPos)
    genomeComb = paste(genome$markerChrom,genome$markerPos)
    matchedVariants = which(variantComb %in% genomeComb)
    if(length(matchedVariants)>0)
    {
      variantNames[matchedVariants] = genome$markerNames[match(variantComb[matchedVariants],genomeComb)]
    }
  }
  variantMap = NULL;
  for(chrom in unique(variantChrom))
  {
    indexGenome = which(genome$markerChrom==chrom)
    indexVariants = which(variantChrom==chrom)
    variantMap = c(variantMap,approx(x=genome$markerPos[indexGenome],xout = variantPos[indexVariants],y = genome$map[indexGenome],rule=2)$y)
  }
  variantRecomb = as.numeric(unlist(tapply(variantMap,variantChrom,recombProbFromGeneticDistance)))
  newGenome(name = genome$name,Nchrom = length(unique(variantChrom)),markerNames = variantNames,markerChrom = variantChrom,markerPos = variantPos,recombProbs =  variantRecomb,map = variantMap)
}

getAndUpdateRecombSpot = function(genome,chrom)
{
  out = unlist(genome$.recombInd[[chrom]],use.names = F)
  if (genome$.recombInd[[chrom]] == genome$.Nrecs)
  {
    genome$.recombInd[[chrom]] = 1
  } else
  {
    genome$.recombInd[[chrom]] = unlist(genome$.recombInd[[chrom]]) + 1
  }
  return(unlist(genome$.recombPos[[chrom]],use.names = F)[out])
}
