library(BLMA)
# load KEGG pathways
x <- loadKEGGPathways()
# load example data
dataSets <- c("GSE1297", "GSE4757", "GSE18309", "GSE28146")
data(list=dataSets, package="BLMA")
names(dataSets) <- dataSets
dataList <- lapply(dataSets, function(dataset) get(paste0("data_", dataset)))
groupList <- lapply(dataSets, function(dataset) get(paste0("group_", dataset)))
IAComb <- bilevelAnalysisPathway(x$kpg, x$kpn, dataList, groupList)
head(IAComb[, c("Name", "pBLMA", "pBLMA.fdr", "rBLMA")])


gslist <- lapply(x$kpg,FUN=function(y){return (nodes(y));})
gs.names <- x$kpn[names(gslist)]

GSAComb <- bilevelAnalysisGeneset(gslist, gs.names, dataList, groupList, enrichment = "GSA", nperms = 200, random.seed = 1)
head(GSAComb[, c("Name", "pBLMA", "pBLMA.fdr", "rBLMA")])