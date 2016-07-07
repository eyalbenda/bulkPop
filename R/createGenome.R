#' @import qtl
#' @import bit

#' @export
is.genome = function(x)
  inherits(x,"genome")
#' @export
newGenome = function(name,Nchrom,markerNames,markerChrom,markerPos,recombProbs, map)
{
  markers
  structure(list(name = name,Nchrom = Nchrom,markerNames = markerNames,markerChrom = markerChrom,markerPos = markerPos,recombProbs = recombProbs, map = map),class="genome")
}

#' This function assumes a single recombination event per chromosome!!!
#' @param distances A vector of genetic distances between the markers
#' @export
recombProbFromGeneticDistance = function(distances){
  distChrom = distances[-1] - distances[-length(distances)]
  finalProbs = c(0,distChrom / sum(distChrom))
  names(finalProbs) = names(distances)
  finalProbs
}


#' The function expands a
#' @param genome a genome object
#' @param variantNames vector of variant names. If not given will be generated as chrom:pos
#' @param variantChrom vector of variant chromsomes. Names must match genome object
#' @param variantPos vector of variant positions.
#' @param matchNames boolean: If true, the function will look for the new variants by position in the old genome and carry over the names to the new genome.
#' @export
expandGenomeToVariantlist = function(genome,variantNames = NULL,variantChrom,variantPos,matchNames = F)
{
  if(any(!unique(variantChrom) %in% unique(N2xCB4856.genome$markerChrom)))
    stop("variants in chromosomes that don't exist in supplied genome!")
  if(is.null(variantNames))
    variantNames = paste(variantChrom,variantPos,sep=":")
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
    variantMap = c(variantMap,approx(x=genome$markerPos[indexGenome],xout = variantPos[indexVariants],y = genome$map[indexGenome]))
  }
  variantRecomb = as.numeric(unlist(tapply(variantMap,variantChrom,recombProbFromGeneticDistance)))
  newGenome(name = genome$name,Nchrom = length(unique(variantChrom)),markerNames = variantNames,markerChrom = variantChrom,markerPos = variantPos,recombProbs =  variantRecomb,map = variantMap)
}

