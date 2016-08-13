#' @import ggplot2
#' @import ggthemes

#' @export
parseAllGenotypes = function(pop,genome,summarize=F)
{
  genotypeAndMarkerTable = data.frame(chrom = genome$markerChrom,
                                      pos = genome$markerPos,
                                      marker = genome$markerNames)
  suppressWarnings(hap1 <- do.call(cbind,lapply(pop,function(x)as.numeric(x$hap1))))
  suppressWarnings(hap2 <- do.call(cbind,lapply(pop,function(x)as.numeric(x$hap2))))
  if(summarize)
  {
    genotypeAndMarkerTable["genotypes"] = (rowMeans(hap1,na.rm = T) + rowMeans(hap2,na.rm=T)) / 2
  } else
  {
    outTable = matrix(nrow = nrow(hap1),ncol=ncol(hap1)+ncol(hap2))
    outTable[,seq(1,ncol(outTable),2)] = hap1
    outTable[,seq(2,ncol(outTable),2)] = hap2
    colnames(outTable)[seq(1,ncol(outTable),2)] = paste(1:ncol(hap1),"hap1",sep="_")
    colnames(outTable)[seq(2,ncol(outTable),2)] = paste(1:ncol(hap2),"hap2",sep="_")
    genotypeAndMarkerTable = cbind(genotypeAndMarkerTable,outTable)
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
plotGenotypeAndMarkerList = function(genotypeAndMarkerList,calcFreqs=T,chromPos = NULL, markers = NULL, internal =F)
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
  if(internal)
    return(dfAll)
  ggplot(dfAll,aes(x = pos,y=Frequency,color = generation)) + geom_line() + ylim(c(0,1)) +
    theme_bw() + scale_color_gradient_tableau(palette="Light Red-Green") + facet_wrap(~chrom,nrow=1,scales = "free_x") + theme(axis.text.x = element_text(angle=90,hjust=1))
}
