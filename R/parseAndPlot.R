#' @import ggplot2
#' @import ggthemes
#' @useDynLib bulkPop


#' @export
parseAllGenotypes = function(pop,genome,summarize=F)
{
  genotypeAndMarkerTable = data.frame(chrom = genome$markerChrom,
                                      pos = genome$markerPos,
                                      marker = genome$markerNames)
  hap1 = lapply(pop,function(x)as.logical(x$hap1))
  hap2 = lapply(pop,function(x)as.logical(x$hap2))
  if(summarize)
  {
    hap1 = summarizePop(hap1,max(sapply(hap1,length)))
    hap2 = summarizePop(hap2,max(sapply(hap2,length)))
    genotypeAndMarkerTable["genotypes"] = (hap1 + hap2) / 2
  } else
  {
    targetLengths = max(c(lengths(hap1),lengths(hap2)))
    outTable = matrix(nrow = targetLengths,ncol=length(hap1)+length(hap2),dimnames = list(1:targetLengths,1:(length(hap1)*2)))
    outTable[,seq(1,ncol(outTable),2)] = do.call(rbind,lapply(hap1,function(x){length(x) = targetLengths;as.numeric(x)}))
    outTable[,seq(2,ncol(outTable),2)] = do.call(rbind,lapply(hap2,function(x){length(x) = targetLengths;as.numeric(x)}))
    colnames(outTable)[seq(1,ncol(outTable),2)] = paste(1:length(pop),"hap1",sep="_")
    colnames(outTable)[seq(2,ncol(outTable),2)] = paste(1:length(pop),"hap2",sep="_")
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
