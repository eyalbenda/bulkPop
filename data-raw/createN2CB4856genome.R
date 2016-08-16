library(qtl)
require(bit)
require(hash)
require(nnet)

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
N2xCB4856.genome = newGenome(name = "N2xCB4856",Nchrom = 6, markerNames = elegansNames,markerChrom = elegansChrom,markerPos = elegansPos, recombProbs = recombProbsHawaiiN2, map = elegansMap,Nrecs = 10000)

genome = N2xCB4856.genome
onlyV = N2xCB4856.genome$markerChrom=="V"
genome = expandGenomeToVariantlist(N2xCB4856.genome, N2xCB4856.genome$markerNames[onlyV],N2xCB4856.genome$markerChrom[onlyV],N2xCB4856.genome$markerPos[onlyV],matchNames = T)

N2 = newIndividual(genome,1,sex="f")
Hawaii = newIndividual(genome,0,sex="m")
system.time(sapply(1:10000,function(x)crossIndividuals(N2,Hawaii,N2xCB4856.genome,zeelPeel = zeelPeel)))
f1n = crossIndividuals(N2,Hawaii,genome)
zeelPeel = list(marker = findProxyForPosition(N2xCB4856.genome,"I",2350468),kill = 0)
