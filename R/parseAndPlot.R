#' @import ggplot2
#' @import ggthemes
#' @import Rcpp

Rcpp::cppFunction('NumericVector summarizePop(List hapPop,int maxLength) {
  if(maxLength<1)
  {
    throw std::invalid_argument("length of haplotype must be larger than 0");
  }
  NumericVector counts(maxLength);
  NumericVector size(maxLength);
  List::iterator it;
  for(it = hapPop.begin();it!=hapPop.end();++it)
  {
    LogicalVector Hap = as<LogicalVector>(*it);
    for(int i=0;i<Hap.length();i++)
    {
      if(Hap[i])
      {
        counts[i] += 1;
      }
      size[i]+=1;
    }
  }
  return counts/size;
}
')

#' Function to parse genotype table from a population of individuals and genome object
#' Useful for generating data frames to be used for plotting
#' @param pop a population of individuals
#' @param genome a genome object
#' @param summarize If FALSE, output raw genotypes (big matrix). TRUE: Output allele frequencies
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

#' Plot genotype and marker table
#' @param genotypeAndMarkerTable a table generated with parseAllGenotypes
#' @param calcFreqs Do allele frequencies need to be calculated from the table (see summarize in parseAllGenotypes function)
#' @param internal generate plot (FALSE) or just return plotting dataframe (TRUE)?
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
    pos = NULL;Frequency = NULL ;chrom = NULL
    ggplot(df,aes(x = pos,y=Frequency,color=chrom,fill=chrom)) + geom_point() + ylim(c(0,1)) +
      theme_bw() + scale_color_tableau() + scale_fill_tableau() + facet_wrap(~chrom,nrow=1,scales = "free_x")
  }
}

#' Plot a list of genotypes and markers
#' An extension of plotGenotypeAndMarkerTable to plot multiple populations/generations.
#' @param genotypeAndMarkerList a list where each element is a genotypeAndMarkerTable
#' @param calcFreqs Do allele frequencies need to be calculated from the table (see summarize in parseAllGenotypes function)
#' @param chromPos optional matrix with two columns (chromosome and position). If given, only these positions will be plotted
#' @param markers optional vector of variants. if given, only these variants will be plotted.
#' @param internal generate plot (FALSE) or just return plotting dataframe (TRUE)?
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
  pos = NULL;Frequency = NULL ;generation = NULL;
  ggplot(dfAll,aes(x = pos,y=Frequency,color = generation)) + geom_line() + ylim(c(0,1)) +
    theme_bw() + scale_color_gradient_tableau(palette="Light Red-Green") + facet_wrap(~chrom,nrow=1,scales = "free_x") + theme(axis.text.x = element_text(angle=90,hjust=1))
}
