#' @import ggplot2
#' @import ggthemes

#' @export
parseAllGenotypes = function(pop,genome,summarize=F)
{
  genotypeAndMarkerTable = data.frame(chrom = genome@markerChrom,
                                      pos = genome@markerPos,
                                      marker = genome@markerNames)
  genotypes = data.frame(matrix(nrow=nrow(genotypeAndMarkerTable),ncol=length(pop)))
  k = 1
  if(summarize)
  {
    genotypes = vector("numeric",length = nrow(genotypeAndMarkerTable))
    Ngenotypes = vector("numeric",length = nrow(genotypeAndMarkerTable))
  }
  for(i in 1:length(pop))
  {
    if(summarize)
    {
      genos = as.numeric(slot(pop[[i]],paste("hap",1,sep=""))[[1]])
      genotypes[1:length(genos)] = genotypes[1:length(genos)] + genos
      genos = as.numeric(slot(pop[[i]],paste("hap",2,sep=""))[[1]])
      genotypes[1:length(genos)] = genotypes[1:length(genos)] + genos
      Ngenotypes[1:length(genos)] = Ngenotypes[1:length(genos)] + 2
    } else
    {
      for(hap in 1:2)
      {
        genos = as.numeric(slot(pop[[i]],paste("hap",hap,sep=""))[[1]])
        if(length(genos)<nrow(genotypes)&pop[[i]]@sex==0)
          genos = c(genos,rep(NA,sum(genome@markerChrom=="X")))
        genotypes[,k] = genos
        k = k+1
      }
    }
  }
  if(summarize)
  {
    genotypeAndMarkerTable = data.frame(genotypeAndMarkerTable,genotypes = genotypes/Ngenotypes)
  } else
  {
    genotypeAndMarkerTable = cbind(genotypeAndMarkerTable,genotypes)
  }
  genotypeAndMarkerTable
}


#' @export
plotGenotypeAndMarkerTable = function(genotypeAndMarkerTable,calcFreqs=T,internal = F)
{
  chroms = unique(genotypeAndMarkerTable$chrom)
  if(calcFreqs)
  {
    freqs = rowMeans(genotypeAndMarkerTable[,-c(1:3)],na.rm = T)
  } else
  {
    freqs = genotypeAndMarkerTable[,4]
  }
  df = data.frame(genotypeAndMarkerTable[,1:3],Frequency = freqs)
  if(internal)
    return(df)
  else
  {
    ggplot(df,aes(x = pos,y=Frequency,color=chrom,fill=chrom)) + geom_point() + ylim(c(0,1)) +
      theme_bw() + scale_color_tableau() + scale_fill_tableau() + facet_wrap(~chrom,nrow=1,scales = "free_x")
  }
}

#' @export
plotGenotypeAndMarkerList = function(genotypeAndMarkerList,calcFreqs=T,chromPos = NULL, markers = NULL)
{
  dfAll = NULL
  for(curGen in genotypeAndMarkerList)
  {
    if(!is.null(chromPos)) curGen$genotypeAndMarkerTable =
        curGen$genotypeAndMarkerTable[paste(curGen$genotypeAndMarkerTable[,1],curGen$genotypeAndMarkerTable[,2]) %in% paste(chromPos[,1],chromPos[,2]),]
    if(!is.null(markers))  curGen$genotypeAndMarkerTable =
        curGen$genotypeAndMarkerTable[curGen$genotypeAndMarkerTable$marker %in% markers,]
    dfAll = rbind(dfAll,
                  data.frame(generation = as.numeric(curGen$generation),plotGenotypeAndMarkerTable(curGen$genotypeAndMarkerTable,calcFreqs = calcFreqs,internal=T)))
  }
  ggplot(dfAll,aes(x = pos,y=Frequency)) + geom_line() + ylim(c(0,1)) +
    theme_bw() + scale_color_gradient2_tableau(palette="Light Red-Green") + facet_wrap(~chrom,nrow=1,scales = "free_x") + theme(axis.text.x = element_text(angle=90,hjust=1))
}
