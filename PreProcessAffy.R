
library(affy)
library(oligo)
library(affyPLM)

getGeneExpression = function(data, annotation) {
  ano=paste(substring(annotation, 1, nchar(annotation)-3),"ENTREZID",sep="")
  require(annotation, character.only = TRUE)
  geneAno=as.data.frame(get(ano))
  geneAnno=geneAno[,2]
  names(geneAnno)=geneAno[,1]
  allProbes=unique(geneAno$probe_id)
  
  data = 2^data
  data <- data[rownames(data)%in%allProbes,]
  rownames(data) <- as.character(geneAnno[rownames(data)])
  
  mydat <- aggregate(data,by=list(rownames(data)),FUN=median)
  rownames(mydat)=mydat[,1]
  mydat <- mydat[,-1]
  mydata <- log(mydat,2)
}

path="~/Desktop/Bioinformatics/IS/datasets/"

# GSE12685 GSE167593885
dataset="GSE36980_FC"

for (datasets in c("GSE36980_FC")) {

  filename=paste(path,dataset,"/",dataset,"_CelFiles.txt",sep="")
  cellFiles=as.matrix(read.table(filename,header=F,sep="\t"))
  
  setwd(paste(path,dataset,"/", dataset, "_RAW/",sep=""))
  
  # read the raw data using ReadAffy
  rawData=ReadAffy(filenames=cellFiles)
  
  # normalize the data
  eset=threestep(rawData,background.method="RMA.2",normalize.method="quantile",summary.method="median.polish")
  mydat.data=exprs(eset)
  View(mydat.data)
  
  filename=paste(path,dataset,"/",dataset,"_Groups.txt",sep="")
  #mydat.groups=read.table(filename,header=T,sep="\t")
  #rownames(mydat.groups)=mydat.groups[,1]
  tmp=read.table(filename,header=T,sep="\t")
  mydat.groups <- tmp$Group
  names(mydat.groups) <- colnames(mydat.data) <- tmp$Sample
  mydat.groups
  
  
  
  # mydat.annotation = "hgu133plus2.db"
  filename=paste(path,dataset,"/",dataset,"_annotation.txt",sep="")
  mydat.annotation=as.character((read.table(filename,header=F,sep="\t"))[1,1])
  
  mydat.gene = getGeneExpression(mydat.data, mydat.annotation)
  rownames(mydat.gene) <- paste0("hsa:",rownames(mydat.gene))
  
  dynamicVariableData <- paste("data_",dataset, sep="")
  dynamicVariableGene <- paste("gene_", dataset, sep="")
  dynamicVariableGroup <- paste("group_",dataset, sep="")
  dynamicVariableAnno <- paste("annotation_",dataset, sep="")
  
  assign(dynamicVariableData, mydat.data)
  assign(dynamicVariableGene, mydat.gene)
  assign(dynamicVariableGroup, mydat.groups)
  assign(dynamicVariableAnno, mydat.annotation)
  
  mylist=c(dynamicVariableData, dynamicVariableGene, dynamicVariableGroup, dynamicVariableAnno)
  save(list=mylist, file=paste(path,dataset,"/",dataset,".RData",sep=""))

}

