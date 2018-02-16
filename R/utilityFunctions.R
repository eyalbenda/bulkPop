#' @import parallel

initiateCluster<-function(ncores=7)
{
  makeCluster(getOption("cl.cores", ncores))
}


#' Find closest marker to a position
#' @param genome A genome object to search in
#' @param chrom The chromosome to query
#' @param pos position in the chromosome
#' @return the name of the closest variant
#' @export
findProxyForPosition = function(genome,chrom,pos)
{
  if(length(chrom)!=length(pos))
    stop("Chromosome and Position must have same lengths")
  varsOut = rep(NA,length(chrom))
  for(var in 1:length(chrom))
  {
    curRange = unlist(genome$.chromRanges[[chrom[var]]],use.names = F)
    varsOut[var] = genome$markerNames[curRange[1]:curRange[2]][which.min(abs(genome$markerPos[curRange[1]:curRange[2]]-pos))]
  }
  varsOut
}
