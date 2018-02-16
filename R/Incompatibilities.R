#' Zeel peel selection
#' A sleection based on the peel-1/zeel-1 incompatibility from Seidel et al 2008.
#' @param genome A genome object. corresponds to both parents
#' @param parent1 the first parent of the individual. The function assumes the parents have already undergone recombination (meiosis)
#' @param parent2 the second parent of the individual. The function assumes the parents have already undergone recombination (meiosis)
#' @param sexes A vector of length 2 with the sexes of the parents (0 is male, 1 is female)
#' @param chosenFrom which haplotype was chosen from each parent.
#' @param zeelPeel a list with two elements. marker: the marker which is being selected. kill: which genotype is killed (0 or 1)
#' @export
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
#' Medea selection
#' A sleection based on the sup-35/pha-1 selfish element.
#' @param genome A genome object. corresponds to both parents
#' @param parent1 the first parent of the individual. The function assumes the parents have already undergone recombination (meiosis)
#' @param parent2 the second parent of the individual. The function assumes the parents have already undergone recombination (meiosis)
#' @param sexes A vector of length 2 with the sexes of the parents (0 is male, 1 is female)
#' @param chosenFrom which haplotype was chosen from each parent.
#' @param medea a list with two elements. marker: the marker which is being selected. kill: which genotype is killed (0 or 1)
#' @export
medeaSelection = function(genome,parent1,parent2,sexes,chosenFrom,medea)
{
  parentsCombined = list(parent1,parent2)
  femPar = which(sexes)
  if(is.null(medea$markers)|is.null(medea$kill)) stop(c("medea must be a list with elements markers and kill. Markers must have two elements - toxin and antidote"))
  if(as.logical(getGenotypes(parentsCombined[[femPar]],genome,chosenFrom[[femPar]],var = medea$markers["antidote"]))==medea$kill &
     as.logical(getGenotypes(parentsCombined[[femPar]],genome,3-chosenFrom[[femPar]],var = medea$marker["toxin"]))!=medea$kill)
  {
    if(as.logical(getGenotypes(parentsCombined[[-femPar]],genome,chosenFrom[[-femPar]],var = medea$marker["antidote"])==medea$kill))
      return (TRUE)
  }
  return (FALSE)
}

