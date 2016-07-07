library(qtl)
require(bit)

load(url('https://github.com/AndersenLab/N2xCB4856-RIAILS/blob/master/data/N2xCB4856_RIAILs_Rqtlfiles.RData?raw=true'))
#this loads an R/QTL object called N2xCB4856.cross
elegansMap = pull.map(N2xCB4856.cross)
elegansMap = unlist(elegansMap)
names(elegansMap) = gsub("^[I,II,III,IV,V,X]+\\.","",names(elegansMap))

newCrossData = read.csv(url("http://www.g3journal.org/content/suppl/2015/03/13/g3.115.017178.DC1/FileS1.csv"),header=T)

elegansPos= newCrossData$WS244.pos[match(names(elegansMap),newCrossData$id)]
elegansChrom= as.character(newCrossData$chr[match(names(elegansMap),newCrossData$id)])
elegansNames = names(elegansMap)
recombProbsHawaiiN2 = as.numeric(unlist(tapply(elegansMap,elegansChrom,recombProbFromGeneticDistance)))
N2xCB4856.genome = newGenome(name = "N2xCB4856",Nchrom = 6, markerNames = elegansNames,markerChrom = elegansChrom,markerPos = elegansPos, recombProbs = recombProbsHawaiiN2, map = elegansMap)

