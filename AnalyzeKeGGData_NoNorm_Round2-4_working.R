

require(devEMF)


library(gplots)
library(scales)
library(RColorBrewer)
library(corrplot)
library(som)
library(cluster)
library(e1071)

#####library(pcaMethods)

library(factoextra)


##load modules
source("/Users/dybasj/JoeRLib/ProteomicsAnalysisFunctions.R")
source("/Users/dybasj/JoeRLib/PlottingFunctions.R")


##formatted with experiments as columns and samples as rows
KeGGdata_run2 <- read.table("./datafiles/KeGG_GlyGly_Rd2_mbr_NormalizedData.txt", sep="\t", quote="\"", row.names=1, header=TRUE)
#head(KeGGdata_run2)
print("run2 dim")
dim(KeGGdata_run2)
KeGGdata_run3 <- read.table("./datafiles/KeGG_GlyGly_Rd3_mbr_NormalizedData.txt", sep="\t", quote="\"", row.names=1, header=TRUE)
#head(KeGGdata_run3)
print("run3 dim")
dim(KeGGdata_run3)
KeGGdata_run4 <- read.table("./datafiles/KeGG_GlyGly_Rd4_mbr_NormalizedData.txt", sep="\t", quote="\"", row.names=1, header=TRUE)
#head(KeGGdata_run4)
print("run4 dim")
dim(KeGGdata_run4)


##peptides identified in run2, run3, run4
#run2 peptides
KeGGPeptides_Run2_PepName <- row.names(KeGGdata_run2)
KeGGPeptides_Run2_PepName_Num <- length(KeGGPeptides_Run2_PepName)
print("Run2 num peptide id")
print(KeGGPeptides_Run2_PepName_Num)
#run 3 peptides
KeGGPeptides_Run3_PepName <- row.names(KeGGdata_run3)
KeGGPeptides_Run3_PepName_Num <- length(KeGGPeptides_Run3_PepName)
print("Run3 num peptide id")
print(KeGGPeptides_Run3_PepName_Num)
#run 4 peptides
KeGGPeptides_Run4_PepName <- row.names(KeGGdata_run4)
KeGGPeptides_Run4_PepName_Num <- length(KeGGPeptides_Run4_PepName)
print("Run4 num peptide id")
print(KeGGPeptides_Run4_PepName_Num)
#
#run2-run3 intersection
KeGGPeptides_Run2Run3_Int_PepName <- intersect(KeGGPeptides_Run2_PepName, KeGGPeptides_Run3_PepName)
KeGGPeptides_Run2Run3_Int_PepName_Num <- length(KeGGPeptides_Run2Run3_Int_PepName)
print("Run2-Run3 peptide intersection")
print(KeGGPeptides_Run2Run3_Int_PepName_Num)
#
#run2-run4 intersection
KeGGPeptides_Run2Run4_Int_PepName <- intersect(KeGGPeptides_Run2_PepName, KeGGPeptides_Run4_PepName)
KeGGPeptides_Run2Run4_Int_PepName_Num <- length(KeGGPeptides_Run2Run4_Int_PepName)
print("Run2-Run4 peptide intersection")
print(KeGGPeptides_Run2Run4_Int_PepName_Num)
#
#run3-run4 intersection
KeGGPeptides_Run3Run4_Int_PepName <- intersect(KeGGPeptides_Run3_PepName, KeGGPeptides_Run4_PepName)
KeGGPeptides_Run3Run4_Int_PepName_Num <- length(KeGGPeptides_Run3Run4_Int_PepName)
print("Run3-Run4 peptide intersection")
print(KeGGPeptides_Run3Run4_Int_PepName_Num)
#
#intersection of run2-3-4
KeGGPeptides_Run2Run3Run4_Int_PepName <- Reduce(intersect, list(KeGGPeptides_Run2_PepName, KeGGPeptides_Run3_PepName, KeGGPeptides_Run4_PepName))
KeGGPeptides_Run2Run3Run4_Int_PepName_Num <- length(KeGGPeptides_Run2Run3Run4_Int_PepName)
print("Run2-Run3-Run4 peptide intersection")
print(KeGGPeptides_Run2Run3Run4_Int_PepName_Num)
#
#union of run2-3-4
KeGGPeptides_Run2Run3Run4_Union_PepName <- Reduce(union, list(KeGGPeptides_Run2_PepName, KeGGPeptides_Run3_PepName, KeGGPeptides_Run4_PepName))
KeGGPeptides_Run2Run3Run4_Union_PepName_Num <- length(KeGGPeptides_Run2Run3Run4_Union_PepName)
print("Run2-Run3-Run4 peptide union")
print(KeGGPeptides_Run2Run3Run4_Union_PepName_Num)
##end peptides identified in run2, run3, run4


##peptide intersectin triple venn
plottriplevenn(KeGGPeptides_Run2_PepName, KeGGPeptides_Run3_PepName, KeGGPeptides_Run4_PepName, "", "", "", "navy", "dodgerblue4", "steelblue1", "white", "white", "white", "KeGGdata_PepId_Run2Run3Run4_Venn.pdf")


##proteins identified in run2, run3, run4
#run2 prots
KeGGProts_Run2_ProtName <- unique(KeGGdata_run2$Protein.Rd2)
KeGGProts_Run2_ProtName_Num <- length(KeGGProts_Run2_ProtName)
print("Run2 prot id")
print(KeGGProts_Run2_ProtName_Num)
#run3 prots
KeGGProts_Run3_ProtName <- unique(KeGGdata_run3$Protein.Rd3)
KeGGProts_Run3_ProtName_Num <- length(KeGGProts_Run3_ProtName)
print("Run3 prot id")
print(KeGGProts_Run3_ProtName_Num)
#run4 prots
KeGGProts_Run4_ProtName <- unique(KeGGdata_run4$Protein.Rd4)
KeGGProts_Run4_ProtName_Num <- length(KeGGProts_Run3_ProtName)
print("Run4 prot id")
print(KeGGProts_Run4_ProtName_Num)
#
#run2-run3 intersection
KeGGProts_Run2Run3_Int_ProtName <- intersect(KeGGProts_Run2_ProtName, KeGGProts_Run3_ProtName)
KeGGProts_Run2Run3_Int_ProtName_Num <- length(KeGGProts_Run2Run3_Int_ProtName)
print("Run2-Run3 protein intersection")
print(KeGGProts_Run2Run3_Int_ProtName_Num)
#
#run2-run4 intersection
KeGGProts_Run2Run4_Int_ProtName <- intersect(KeGGProts_Run2_ProtName, KeGGProts_Run4_ProtName)
KeGGProts_Run2Run4_Int_ProtName_Num <- length(KeGGProts_Run2Run4_Int_ProtName)
print("Run2-Run4 protein intersection")
print(KeGGProts_Run2Run4_Int_ProtName_Num)
#
#run3-run4 intersection
KeGGProts_Run3Run4_Int_ProtName <- intersect(KeGGProts_Run3_ProtName, KeGGProts_Run4_ProtName)
KeGGProts_Run3Run4_Int_ProtName_Num <- length(KeGGProts_Run3Run4_Int_ProtName)
print("Run3-Run4 protein intersection")
print(KeGGProts_Run3Run4_Int_ProtName_Num)
#
#intersection of run2-3-4
KeGGProts_Run2Run3Run4_Int_ProtName <- Reduce(intersect, list(KeGGProts_Run2_ProtName, KeGGProts_Run3_ProtName, KeGGProts_Run4_ProtName))
KeGGProts_Run2Run3Run4_Int_ProtName_Num <- length(KeGGProts_Run2Run3Run4_Int_ProtName)
print("Run2-Run3-Run4 protein intersection")
print(KeGGProts_Run2Run3Run4_Int_ProtName_Num)
#
#union of run2-3-4
KeGGProts_Run2Run3Run4_Union_ProtName <- Reduce(union, list(KeGGProts_Run2_ProtName, KeGGProts_Run3_ProtName, KeGGProts_Run4_ProtName))
KeGGProts_Run2Run3Run4_Union_ProtName_Num <- length(KeGGProts_Run2Run3Run4_Union_ProtName)
print("Run2-Run3-Run4 protein union")
print(KeGGProts_Run2Run3Run4_Union_ProtName_Num)
##end proteins identified in run2, run3, run4


##peptide intersectin triple venn
plottriplevenn(KeGGProts_Run2_ProtName, KeGGProts_Run3_ProtName, KeGGProts_Run4_ProtName, "", "", "", "navy", "dodgerblue4", "steelblue1", "white", "white", "white", "KeGGdata_ProtId_Run2Run3Run4_Venn.pdf")


##combine dataframes
KeGGdataMerge <- merge(KeGGdata_run2, KeGGdata_run3, by="row.names", all=TRUE)
row.names(KeGGdataMerge)=KeGGdataMerge$Row.names
KeGGdataMerge$Row.names <- NULL
KeGGdataMerge <- merge(KeGGdataMerge, KeGGdata_run4, by="row.names", all=TRUE)
row.names(KeGGdataMerge)=KeGGdataMerge$Row.names
KeGGdataMerge$Row.names <- NULL
dim(KeGGdataMerge)
#
KeGGdataMerge[KeGGdataMerge=="-Inf"] <- NA
#
#function to combine col values
combinecols <- function(x,y,z)
{
	combineval<-ifelse(is.na(x), ifelse(is.na(y), z, y), x)
	return(combineval)
}
#end weight function
#
#combine protein, protein names and gene names columns from each round
KeGGdataMerge$Protein <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) combinecols(x[["Protein.Rd2"]], x[["Protein.Rd3"]], x[["Protein.Rd4"]]))
KeGGdataMerge<-KeGGdataMerge[,c(ncol(KeGGdataMerge),1:(ncol(KeGGdataMerge)-1))]
KeGGdataMerge$Protein.names <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) combinecols(x[["Protein.names.Rd2"]], x[["Protein.names.Rd3"]], x[["Protein.names.Rd4"]]))
KeGGdataMerge<-KeGGdataMerge[,c(1,ncol(KeGGdataMerge),2:(ncol(KeGGdataMerge)-1))]
KeGGdataMerge$Gene.names <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) combinecols(x[["Gene.names.Rd2"]], x[["Gene.names.Rd3"]], x[["Gene.names.Rd4"]]))
KeGGdataMerge<-KeGGdataMerge[,c(1:2,ncol(KeGGdataMerge),3:(ncol(KeGGdataMerge)-1))]
KeGGdataMerge$Protein.Rd2 <- NULL
KeGGdataMerge$Protein.Rd3 <- NULL
KeGGdataMerge$Protein.Rd4 <- NULL
KeGGdataMerge$Protein.names.Rd2 <- NULL
KeGGdataMerge$Protein.names.Rd3 <- NULL
KeGGdataMerge$Protein.names.Rd4 <- NULL
KeGGdataMerge$Gene.names.Rd2 <- NULL
KeGGdataMerge$Gene.names.Rd3 <- NULL
KeGGdataMerge$Gene.names.Rd4 <- NULL
##end combine dataframes
#
#function to combine gene name and KeGG position
genenamekeggpos <- function(x,y,i)
{
	#print(x)
	#print(y)
	#print(i)
	protpos<-x[i]
	protpossplit<-strsplit(protpos, "_")
	keggpos<-ifelse(length(protpossplit[[1]])>1, protpossplit[[1]][2], stop("kegg pos error"))
	genenamesplit<-strsplit(y,";")
	firstgene<-genenamesplit[[1]][1]
	genekeggpos<-paste(firstgene,keggpos,sep='_')
	return(genekeggpos)
}
protpos<-row.names(KeGGdataMerge)
i<-0
KeGGdataMerge$Gene.KeGGPos <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) genenamekeggpos(protpos, x[["Gene.names"]], i<<-i+1))
matchproteinid <- function(x,i)
{
	protpos<-x[i]
	protpossplit<-strsplit(protpos, "_")
	protid<-ifelse(length(protpossplit[[1]])>1, protpossplit[[1]][1], stop("prot id error"))
	return(protid)
}
i<-0
KeGGdataMerge$Protein <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) matchproteinid(protpos, i<<-i+1))
firstgenename <- function(x)
{
	genelist<-strsplit(x, ";")
	firstgene<-genelist[[1]][1]
	return(firstgene)
}
KeGGdataMerge$Gene.names.first<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) firstgenename(x[["Gene.names"]]))
#write.table(KeGGdataMerge, file="KeGGdataMerge_test.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)


##count peptides and proteins found at each time point for each run
##run2 0hr
#pep identified
KeGGPeptides_0hr_Run2_Data <- subset(KeGGdataMerge, !is.na(KeGGdataMerge$Intensity.0hr_Log2_MedNorm.Rd2))
KeGGPeptides_0hr_Run2_PepId <- row.names(KeGGPeptides_0hr_Run2_Data)
KeGGPeptides_0hr_Run2_PepId_Num <- length(KeGGPeptides_0hr_Run2_PepId)
print("KeGGPeptides_0hr_Run2_PepId_Num")
print(KeGGPeptides_0hr_Run2_PepId_Num)
#prot identified
KeGGPeptides_0hr_Run2_ProtId <- unique(KeGGPeptides_0hr_Run2_Data$Protein)
KeGGPeptides_0hr_Run2_ProtId_Num <- length(KeGGPeptides_0hr_Run2_ProtId)
print("KeGGPeptides_0hr_Run2_ProtId_Num")
print(KeGGPeptides_0hr_Run2_ProtId_Num)
#
##run2 6hr
#pep identified
KeGGPeptides_6hr_Run2_Data <- subset(KeGGdataMerge, !is.na(KeGGdataMerge$Intensity.6hr_Log2_MedNorm.Rd2))
KeGGPeptides_6hr_Run2_PepId <- row.names(KeGGPeptides_6hr_Run2_Data)
KeGGPeptides_6hr_Run2_PepId_Num <- length(KeGGPeptides_6hr_Run2_PepId)
print("KeGGPeptides_6hr_Run2_PepId_Num")
print(KeGGPeptides_6hr_Run2_PepId_Num)
#prot identified
KeGGPeptides_6hr_Run2_ProtId <- unique(KeGGPeptides_6hr_Run2_Data$Protein)
KeGGPeptides_6hr_Run2_ProtId_Num <- length(KeGGPeptides_6hr_Run2_ProtId)
print("KeGGPeptides_6hr_Run2_ProtId_Num")
print(KeGGPeptides_6hr_Run2_ProtId_Num)
#
##run2 8hr
#pep identified
KeGGPeptides_8hr_Run2_Data <- subset(KeGGdataMerge, !is.na(KeGGdataMerge$Intensity.8hr_Log2_MedNorm.Rd2))
KeGGPeptides_8hr_Run2_PepId <- row.names(KeGGPeptides_8hr_Run2_Data)
KeGGPeptides_8hr_Run2_PepId_Num <- length(KeGGPeptides_8hr_Run2_PepId)
print("KeGGPeptides_8hr_Run2_PepId_Num")
print(KeGGPeptides_8hr_Run2_PepId_Num)
#prot identifed
KeGGPeptides_8hr_Run2_ProtId <- unique(KeGGPeptides_8hr_Run2_Data$Protein)
KeGGPeptides_8hr_Run2_ProtId_Num <- length(KeGGPeptides_8hr_Run2_ProtId)
print("KeGGPeptides_8hr_Run2_ProtId_Num")
print(KeGGPeptides_8hr_Run2_ProtId_Num)
#
##run2 10hr
#pep identified
KeGGPeptides_10hr_Run2_Data <- subset(KeGGdataMerge, !is.na(KeGGdataMerge$Intensity.10hr_Log2_MedNorm.Rd2))
KeGGPeptides_10hr_Run2_PepId <- row.names(KeGGPeptides_10hr_Run2_Data)
KeGGPeptides_10hr_Run2_PepId_Num <- length(KeGGPeptides_10hr_Run2_PepId)
print("KeGGPeptides_10hr_Run2_PepId_Num")
print(KeGGPeptides_10hr_Run2_PepId_Num)
#prot identified
KeGGPeptides_10hr_Run2_ProtId <- unique(KeGGPeptides_10hr_Run2_Data$Protein)
KeGGPeptides_10hr_Run2_ProtId_Num <- length(KeGGPeptides_10hr_Run2_ProtId)
print("KeGGPeptides_10hr_Run2_ProtId_Num")
print(KeGGPeptides_10hr_Run2_ProtId_Num)
#
##run3 0hr
#pep identified
KeGGPeptides_0hr_Run3_Data <- subset(KeGGdataMerge, !is.na(KeGGdataMerge$Intensity.0hr_Log2_MedNorm.Rd3))
KeGGPeptides_0hr_Run3_PepId <- row.names(KeGGPeptides_0hr_Run3_Data)
KeGGPeptides_0hr_Run3_PepId_Num <- length(KeGGPeptides_0hr_Run3_PepId)
print("KeGGPeptides_0hr_Run3_PepId_Num")
print(KeGGPeptides_0hr_Run3_PepId_Num)
#prot identified
KeGGPeptides_0hr_Run3_ProtId <- unique(KeGGPeptides_0hr_Run3_Data$Protein)
KeGGPeptides_0hr_Run3_ProtId_Num <- length(KeGGPeptides_0hr_Run3_ProtId)
print("KeGGPeptides_0hr_Run3_ProtId_Num")
print(KeGGPeptides_0hr_Run3_ProtId_Num)
#
##run3 6hr
#pep identified
KeGGPeptides_6hr_Run3_Data <- subset(KeGGdataMerge, !is.na(KeGGdataMerge$Intensity.6hr_Log2_MedNorm.Rd3))
KeGGPeptides_6hr_Run3_PepId <- row.names(KeGGPeptides_6hr_Run3_Data)
KeGGPeptides_6hr_Run3_PepId_Num <- length(KeGGPeptides_6hr_Run3_PepId)
print("KeGGPeptides_6hr_Run3_PepId_Num")
print(KeGGPeptides_6hr_Run3_PepId_Num)
#prot identified
KeGGPeptides_6hr_Run3_ProtId <- unique(KeGGPeptides_6hr_Run3_Data$Protein)
KeGGPeptides_6hr_Run3_ProtId_Num <- length(KeGGPeptides_6hr_Run3_ProtId)
print("KeGGPeptides_6hr_Run3_ProtId_Num")
print(KeGGPeptides_6hr_Run3_ProtId_Num)
#
##run3 8hr
#pep identified
KeGGPeptides_8hr_Run3_Data <- subset(KeGGdataMerge, !is.na(KeGGdataMerge$Intensity.8hr_Log2_MedNorm.Rd3))
KeGGPeptides_8hr_Run3_PepId <- row.names(KeGGPeptides_8hr_Run3_Data)
KeGGPeptides_8hr_Run3_PepId_Num <- length(KeGGPeptides_8hr_Run3_PepId)
print("KeGGPeptides_8hr_Run3_PepId_Num")
print(KeGGPeptides_8hr_Run3_PepId_Num)
#prot identified
KeGGPeptides_8hr_Run3_ProtId <- unique(KeGGPeptides_8hr_Run3_Data$Protein)
KeGGPeptides_8hr_Run3_ProtId_Num <- length(KeGGPeptides_8hr_Run3_ProtId)
print("KeGGPeptides_8hr_Run3_ProtId_Num")
print(KeGGPeptides_8hr_Run3_ProtId_Num)
#
##run3 10hr
#pep identified
KeGGPeptides_10hr_Run3_Data <- subset(KeGGdataMerge, !is.na(KeGGdataMerge$Intensity.10hr_Log2_MedNorm.Rd3))
KeGGPeptides_10hr_Run3_PepId <- row.names(KeGGPeptides_10hr_Run3_Data)
KeGGPeptides_10hr_Run3_PepId_Num <- length(KeGGPeptides_10hr_Run3_PepId)
print("KeGGPeptides_10hr_Run3_PepId_Num")
print(KeGGPeptides_10hr_Run3_PepId_Num)
#prot identified
KeGGPeptides_10hr_Run3_ProtId <- unique(KeGGPeptides_10hr_Run3_Data$Protein)
KeGGPeptides_10hr_Run3_ProtId_Num <- length(KeGGPeptides_10hr_Run3_ProtId)
print("KeGGPeptides_10hr_Run3_ProtId_Num")
print(KeGGPeptides_10hr_Run3_ProtId_Num)
#
##run4 0hr
#pep identified
KeGGPeptides_0hr_Run4_Data <- subset(KeGGdataMerge, !is.na(KeGGdataMerge$Intensity.0hr_Log2_MedNorm.Rd4))
KeGGPeptides_0hr_Run4_PepId <- row.names(KeGGPeptides_0hr_Run4_Data)
KeGGPeptides_0hr_Run4_PepId_Num <- length(KeGGPeptides_0hr_Run4_PepId)
print("KeGGPeptides_0hr_Run4_PepId_Num")
print(KeGGPeptides_0hr_Run4_PepId_Num)
#prot identified
KeGGPeptides_0hr_Run4_ProtId <- unique(KeGGPeptides_0hr_Run4_Data$Protein)
KeGGPeptides_0hr_Run4_ProtId_Num <- length(KeGGPeptides_0hr_Run4_ProtId)
print("KeGGPeptides_0hr_Run4_ProtId_Num")
print(KeGGPeptides_0hr_Run4_ProtId_Num)
#
##run4 6hr
#pep identified
KeGGPeptides_6hr_Run4_Data <- subset(KeGGdataMerge, !is.na(KeGGdataMerge$Intensity.6hr_Log2_MedNorm.Rd4))
KeGGPeptides_6hr_Run4_PepId <- row.names(KeGGPeptides_6hr_Run4_Data)
KeGGPeptides_6hr_Run4_PepId_Num <- length(KeGGPeptides_6hr_Run4_PepId)
print("KeGGPeptides_6hr_Run4_PepId_Num")
print(KeGGPeptides_6hr_Run4_PepId_Num)
#prot identified
KeGGPeptides_6hr_Run4_ProtId <- unique(KeGGPeptides_6hr_Run4_Data$Protein)
KeGGPeptides_6hr_Run4_ProtId_Num <- length(KeGGPeptides_6hr_Run4_ProtId)
print("KeGGPeptides_6hr_Run4_ProtId_Num")
print(KeGGPeptides_6hr_Run4_ProtId_Num)
#
##run4 8hr
#pep identified
KeGGPeptides_8hr_Run4_Data <- subset(KeGGdataMerge, !is.na(KeGGdataMerge$Intensity.8hr_Log2_MedNorm.Rd4))
KeGGPeptides_8hr_Run4_PepId <- row.names(KeGGPeptides_8hr_Run4_Data)
KeGGPeptides_8hr_Run4_PepId_Num <- length(KeGGPeptides_8hr_Run4_PepId)
print("KeGGPeptides_8hr_Run4_PepId_Num")
print(KeGGPeptides_8hr_Run4_PepId_Num)
#prot identified
KeGGPeptides_8hr_Run4_ProtId <- unique(KeGGPeptides_8hr_Run4_Data$Protein)
KeGGPeptides_8hr_Run4_ProtId_Num <- length(KeGGPeptides_8hr_Run4_ProtId)
print("KeGGPeptides_8hr_Run4_ProtId_Num")
print(KeGGPeptides_8hr_Run4_ProtId_Num)
#
##run4 10hr
#pep identified
KeGGPeptides_10hr_Run4_Data <- subset(KeGGdataMerge, !is.na(KeGGdataMerge$Intensity.10hr_Log2_MedNorm.Rd4))
KeGGPeptides_10hr_Run4_PepId <- row.names(KeGGPeptides_10hr_Run4_Data)
KeGGPeptides_10hr_Run4_PepId_Num <- length(KeGGPeptides_10hr_Run4_PepId)
print("KeGGPeptides_10hr_Run4_PepId_Num")
print(KeGGPeptides_10hr_Run4_PepId_Num)
#prot identified
KeGGPeptides_10hr_Run4_ProtId <- unique(KeGGPeptides_10hr_Run4_Data$Protein)
KeGGPeptides_10hr_Run4_ProtId_Num <- length(KeGGPeptides_10hr_Run4_ProtId)
print("KeGGPeptides_10hr_Run4_ProtId_Num")
print(KeGGPeptides_10hr_Run4_ProtId_Num)
#
##end count peptides and proteins found at each time point for each run

#plot bar chart of number of peptides identified in unfiltered and filtered replicates
run2pepnum <- c(KeGGPeptides_0hr_Run2_PepId_Num, KeGGPeptides_6hr_Run2_PepId_Num, KeGGPeptides_8hr_Run2_PepId_Num, KeGGPeptides_10hr_Run2_PepId_Num)
run3pepnum <- c(KeGGPeptides_0hr_Run3_PepId_Num, KeGGPeptides_6hr_Run3_PepId_Num, KeGGPeptides_8hr_Run3_PepId_Num, KeGGPeptides_10hr_Run3_PepId_Num)
run4pepnum <- c(KeGGPeptides_0hr_Run4_PepId_Num, KeGGPeptides_6hr_Run4_PepId_Num, KeGGPeptides_8hr_Run4_PepId_Num, KeGGPeptides_10hr_Run4_PepId_Num)
PepNumTable<-matrix(0,3,4)
PepNumTable[1,]<-run2pepnum
PepNumTable[2,]<-run3pepnum
PepNumTable[3,]<-run4pepnum
#pdf("KeGGdata_PepId_Run2Run3Run4_barchart.pdf", height=6, width=7)
pdf("junk.pdf", height=6, width=7)
par(mar=c(5.1, 5.1, 4.1, 1.1))
barplot(PepNumTable, beside=T, names.arg=c("0hr", "6hr", "8hr", "10hr"), xlab="", ylab="", cex.axis=2, cex.names=2, las=2, col=c("navy", "dodgerblue4", "steelblue1", "navy", "dodgerblue4", "steelblue1", "navy", "dodgerblue4", "steelblue1", "navy", "dodgerblue4", "steelblue1"))
box(bty="l", lwd=3)
#mtext("peptides identified", side=2, line=3.25, cex=1)
dev.off()
##end Ub fold change number bar graph


#plot bar chart of number of peptides identified in unfiltered and filtered replicates
run2protnum <- c(KeGGPeptides_0hr_Run2_ProtId_Num, KeGGPeptides_6hr_Run2_ProtId_Num, KeGGPeptides_8hr_Run2_ProtId_Num, KeGGPeptides_10hr_Run2_ProtId_Num)
run3protnum <- c(KeGGPeptides_0hr_Run3_ProtId_Num, KeGGPeptides_6hr_Run3_ProtId_Num, KeGGPeptides_8hr_Run3_ProtId_Num, KeGGPeptides_10hr_Run3_ProtId_Num)
run4protnum <- c(KeGGPeptides_0hr_Run4_ProtId_Num, KeGGPeptides_6hr_Run4_ProtId_Num, KeGGPeptides_8hr_Run4_ProtId_Num, KeGGPeptides_10hr_Run4_ProtId_Num)
ProtNumTable<-matrix(0,3,4)
ProtNumTable[1,]<-run2protnum
ProtNumTable[2,]<-run3protnum
ProtNumTable[3,]<-run4protnum
#pdf("KeGGdata_ProtId_Run2Run3Run4_barchart.pdf", height=6, width=7)
pdf("junk2.pdf", height=6, width=7)
par(mar=c(5.1, 5.1, 4.1, 1.1))
barplot(ProtNumTable, beside=T, names.arg=c("0hr", "6hr", "8hr", "10hr"), xlab="", ylab="", cex.axis=2, cex.names=2, las=2, col=c("navy", "dodgerblue4", "steelblue1", "navy", "dodgerblue4", "steelblue1", "navy", "dodgerblue4", "steelblue1", "navy", "dodgerblue4", "steelblue1"))
box(bty="l", lwd=3)
#mtext("peptides identified", side=2, line=3.25, cex=1)
dev.off()
##end Ub fold change number bar graph


##count all peptides identified at each time point
#0hr
KeGGPeptides_0hr_PepId <- Reduce(union, list(KeGGPeptides_0hr_Run2_PepId, KeGGPeptides_0hr_Run3_PepId, KeGGPeptides_0hr_Run4_PepId))
KeGGPeptides_0hr_PepId_Num <- length(KeGGPeptides_0hr_PepId)
print("KeGGPeptides_0hr_PepId_Num")
print(KeGGPeptides_0hr_PepId_Num)
#6hr
KeGGPeptides_6hr_PepId <- Reduce(union, list(KeGGPeptides_6hr_Run2_PepId, KeGGPeptides_6hr_Run3_PepId, KeGGPeptides_6hr_Run4_PepId))
KeGGPeptides_6hr_PepId_Num <- length(KeGGPeptides_6hr_PepId)
print("KeGGPeptides_6hr_PepId_Num")
print(KeGGPeptides_6hr_PepId_Num)
#8hr
KeGGPeptides_8hr_PepId <- Reduce(union, list(KeGGPeptides_8hr_Run2_PepId, KeGGPeptides_8hr_Run3_PepId, KeGGPeptides_8hr_Run4_PepId))
KeGGPeptides_8hr_PepId_Num <- length(KeGGPeptides_8hr_PepId)
print("KeGGPeptides_8hr_PepId_Num")
print(KeGGPeptides_8hr_PepId_Num)
#10hr
KeGGPeptides_10hr_PepId <- Reduce(union, list(KeGGPeptides_10hr_Run2_PepId, KeGGPeptides_10hr_Run3_PepId, KeGGPeptides_10hr_Run4_PepId))
KeGGPeptides_10hr_PepId_Num <- length(KeGGPeptides_10hr_PepId)
print("KeGGPeptides_10hr_PepId_Num")
print(KeGGPeptides_10hr_PepId_Num)
#
##count all proteins identified at each time point
#0hr
KeGGPeptides_0hr_ProtId <- Reduce(union, list(KeGGPeptides_0hr_Run2_ProtId, KeGGPeptides_0hr_Run3_ProtId, KeGGPeptides_0hr_Run4_ProtId))
KeGGPeptides_0hr_ProtId_Num <- length(KeGGPeptides_0hr_ProtId)
print("KeGGPeptides_0hr_ProtId_Num")
print(KeGGPeptides_0hr_ProtId_Num)
#6hr
KeGGPeptides_6hr_ProtId <- Reduce(union, list(KeGGPeptides_6hr_Run2_ProtId, KeGGPeptides_6hr_Run3_ProtId, KeGGPeptides_6hr_Run4_ProtId))
KeGGPeptides_6hr_ProtId_Num <- length(KeGGPeptides_6hr_ProtId)
print("KeGGPeptides_6hr_ProtId_Num")
print(KeGGPeptides_6hr_ProtId_Num)
#8hr
KeGGPeptides_8hr_ProtId <- Reduce(union, list(KeGGPeptides_8hr_Run2_ProtId, KeGGPeptides_8hr_Run3_ProtId, KeGGPeptides_8hr_Run4_ProtId))
KeGGPeptides_8hr_ProtId_Num <- length(KeGGPeptides_8hr_ProtId)
print("KeGGPeptides_8hr_ProtId_Num")
print(KeGGPeptides_8hr_ProtId_Num)
#10hr
KeGGPeptides_10hr_ProtId <- Reduce(union, list(KeGGPeptides_10hr_Run2_ProtId, KeGGPeptides_10hr_Run3_ProtId, KeGGPeptides_10hr_Run4_ProtId))
KeGGPeptides_10hr_ProtId_Num <- length(KeGGPeptides_10hr_ProtId)
print("KeGGPeptides_10hr_ProtId_Num")
print(KeGGPeptides_10hr_ProtId_Num)
#
##count replicates with identified protein
##use apply function to parse each row of the dataframe
##count replicates by using sum of non-NA intensity values
#0hr
KeGGdataMerge$KeGGInt.0hr_NumReps <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) sum(!is.na(x[["Intensity.0hr_Log2_MedNorm.Rd2"]]),!is.na(x[["Intensity.0hr_Log2_MedNorm.Rd3"]]),!is.na(x[["Intensity.0hr_Log2_MedNorm.Rd4"]])))
#6hr
KeGGdataMerge$KeGGInt.6hr_NumReps <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) sum(!is.na(x[["Intensity.6hr_Log2_MedNorm.Rd2"]]),!is.na(x[["Intensity.6hr_Log2_MedNorm.Rd3"]]),!is.na(x[["Intensity.6hr_Log2_MedNorm.Rd4"]])))
#8hr
KeGGdataMerge$KeGGInt.8hr_NumReps <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) sum(!is.na(x[["Intensity.8hr_Log2_MedNorm.Rd2"]]),!is.na(x[["Intensity.8hr_Log2_MedNorm.Rd3"]]),!is.na(x[["Intensity.8hr_Log2_MedNorm.Rd4"]])))
#10hr
KeGGdataMerge$KeGGInt.10hr_NumReps <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) sum(!is.na(x[["Intensity.10hr_Log2_MedNorm.Rd2"]]),!is.na(x[["Intensity.10hr_Log2_MedNorm.Rd3"]]),!is.na(x[["Intensity.10hr_Log2_MedNorm.Rd4"]])))
#
##count peptides with KeGG quantification in at least two replicates
#0hr
#peptide
KeGGPeptides_0hr_2reps_Data <- subset(KeGGdataMerge, KeGGdataMerge$KeGGInt.0hr_NumReps>=2)
KeGGPeptides_0hr_2reps_PepId <- row.names(KeGGPeptides_0hr_2reps_Data)
KeGGPeptides_0hr_2reps_PepId_Num <- length(KeGGPeptides_0hr_2reps_PepId)
print("KeGGPeptides_0hr_2reps_PepId_Num")
print(KeGGPeptides_0hr_2reps_PepId_Num)
#protein
KeGGPeptides_0hr_2reps_ProtId <- unique(KeGGPeptides_0hr_2reps_Data$Protein)
KeGGPeptides_0hr_2reps_ProtId_Num <- length(KeGGPeptides_0hr_2reps_ProtId)
print("KeGGPeptides_0hr_2reps_ProtId_Num")
print(KeGGPeptides_0hr_2reps_ProtId_Num)
#6hr
#peptide
KeGGPeptides_6hr_2reps_Data <- subset(KeGGdataMerge, KeGGdataMerge$KeGGInt.6hr_NumReps>=2)
KeGGPeptides_6hr_2reps_PepId <- row.names(KeGGPeptides_6hr_2reps_Data)
KeGGPeptides_6hr_2reps_PepId_Num <- length(KeGGPeptides_6hr_2reps_PepId)
print("KeGGPeptides_6hr_2reps_PepId_Num")
print(KeGGPeptides_6hr_2reps_PepId_Num)
#protein
KeGGPeptides_6hr_2reps_ProtId <- unique(KeGGPeptides_6hr_2reps_Data$Protein)
KeGGPeptides_6hr_2reps_ProtId_Num <- length(KeGGPeptides_6hr_2reps_ProtId)
print("KeGGPeptides_6hr_2reps_ProtId_Num")
print(KeGGPeptides_6hr_2reps_ProtId_Num)
#8hr
#peptide
KeGGPeptides_8hr_2reps_Data <- subset(KeGGdataMerge, KeGGdataMerge$KeGGInt.8hr_NumReps>=2)
KeGGPeptides_8hr_2reps_PepId <- row.names(KeGGPeptides_8hr_2reps_Data)
KeGGPeptides_8hr_2reps_PepId_Num <- length(KeGGPeptides_8hr_2reps_PepId)
print("KeGGPeptides_8hr_2reps_PepId_Num")
print(KeGGPeptides_8hr_2reps_PepId_Num)
#protein
KeGGPeptides_8hr_2reps_ProtId <- unique(KeGGPeptides_8hr_2reps_Data$Protein)
KeGGPeptides_8hr_2reps_ProtId_Num <- length(KeGGPeptides_8hr_2reps_ProtId)
print("KeGGPeptides_8hr_2reps_ProtId_Num")
print(KeGGPeptides_8hr_2reps_ProtId_Num)
#10hr
#peptide
KeGGPeptides_10hr_2reps_Data <- subset(KeGGdataMerge, KeGGdataMerge$KeGGInt.10hr_NumReps>=2)
KeGGPeptides_10hr_2reps_PepId <- row.names(KeGGPeptides_10hr_2reps_Data)
KeGGPeptides_10hr_2reps_PepId_Num <- length(KeGGPeptides_10hr_2reps_PepId)
print("KeGGPeptides_10hr_2reps_PepId_Num")
print(KeGGPeptides_10hr_2reps_PepId_Num)
#protein
KeGGPeptides_10hr_2reps_ProtId <- unique(KeGGPeptides_10hr_2reps_Data$Protein)
KeGGPeptides_10hr_2reps_ProtId_Num <- length(KeGGPeptides_10hr_2reps_ProtId)
print("KeGGPeptides_10hr_2reps_ProtId_Num")
print(KeGGPeptides_10hr_2reps_ProtId_Num)
#
##remove proteins with sparce data
#proteins w/ 0 or 1 KeGG replicate identification for all time points
allNAprots <- row.names(subset(KeGGdataMerge, (KeGGdataMerge$KeGGInt.0hr_NumReps<=1 & KeGGdataMerge$KeGGInt.6hr_NumReps<=1 & KeGGdataMerge$KeGGInt.8hr_NumReps<=1 & KeGGdataMerge$KeGGInt.10hr_NumReps<=1)))
allNAprots_num <- length(allNAprots)
print("proteins w/ no quantifiations")
print(allNAprots_num)
dim(KeGGdataMerge)
KeGGdataMerge<-KeGGdataMerge[!row.names(KeGGdataMerge)%in%allNAprots,]
dim(KeGGdataMerge)
##end count replicates with identified protein
#
#plot bar chart of number of peptides identified in unfiltered and filtered replicates
filteredpepnum <- c(KeGGPeptides_0hr_2reps_PepId_Num, KeGGPeptides_6hr_2reps_PepId_Num, KeGGPeptides_8hr_2reps_PepId_Num, KeGGPeptides_10hr_2reps_PepId_Num)
filteredprotnum <- c(KeGGPeptides_0hr_2reps_ProtId_Num, KeGGPeptides_6hr_2reps_ProtId_Num, KeGGPeptides_8hr_2reps_ProtId_Num, KeGGPeptides_10hr_2reps_ProtId_Num)
PepProtNumTable<-matrix(0,2,4)
PepProtNumTable[1,]<-filteredpepnum
PepProtNumTable[2,]<-filteredprotnum
pdf("KeGGdata_PepIdProtId_0hr6hr8hr10hr_barchart.pdf", height=6, width=7)
#pdf("Run2Run3Run4_IdentifiedPepNum_barchart.pdf", height=6, width=7)
par(mar=c(5.1, 5.1, 4.1, 1.1))
barplot(PepProtNumTable, beside=T, names.arg=c("0hr", "6hr", "8hr", "10hr"), xlab="", ylab="", cex.axis=2, cex.names=2, las=2, col=c("grey70", "grey30", "skyblue", "dodgerblue", "skyblue", "dodgerblue", "skyblue", "dodgerblue"))
box(bty="l", lwd=3)
#mtext("peptides identified", side=2, line=3.25, cex=1)
dev.off()
##end Ub fold change number bar graph


####filtered peptide intersectin triple venn
##plottriplevenn(KeGGPeptides_0hr_2reps_PepId, KeGGPeptides_6hr_2reps_PepId, KeGGPeptides_8hr_2reps_PepId, KeGGPeptides_10hr_2reps_PepId, "", "", "", "navy", "dodgerblue4", "steelblue1", "white", "white", "white", "KeGGdata_ProtId_Run2Run3Run4_Venn.pdf")


##average the non-normalized intensity values for three replicates within each time-points
##use apply function to parse each row of the dataframe, use mean function to find averages
#0hr
KeGGdataMerge$Intensity.0hr_Log2MedNorm_Avg<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd2"]],x[["Intensity.0hr_Log2_MedNorm.Rd3"]],x[["Intensity.0hr_Log2_MedNorm.Rd4"]])), na.rm=TRUE))
#6hr
KeGGdataMerge$Intensity.6hr_Log2MedNorm_Avg<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["Intensity.6hr_Log2_MedNorm.Rd2"]],x[["Intensity.6hr_Log2_MedNorm.Rd3"]],x[["Intensity.6hr_Log2_MedNorm.Rd4"]])), na.rm=TRUE))
#8hr
KeGGdataMerge$Intensity.8hr_Log2MedNorm_Avg<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["Intensity.8hr_Log2_MedNorm.Rd2"]],x[["Intensity.8hr_Log2_MedNorm.Rd3"]],x[["Intensity.8hr_Log2_MedNorm.Rd4"]])), na.rm=TRUE))
#10hr
KeGGdataMerge$Intensity.10hr_Log2MedNorm_Avg<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["Intensity.10hr_Log2_MedNorm.Rd2"]],x[["Intensity.10hr_Log2_MedNorm.Rd3"]],x[["Intensity.10hr_Log2_MedNorm.Rd4"]])), na.rm=TRUE))
##end average the intensity values with time-points


##z-scores for average non-normalize intensity values
#0hr
KeGGInt0hr_Mean <- mean(as.numeric(KeGGdataMerge$Intensity.0hr_Log2MedNorm_Avg), na.rm=TRUE)
KeGGInt0hr_Sd <- sd(as.numeric(KeGGdataMerge$Intensity.0hr_Log2MedNorm_Avg), na.rm=TRUE)
print("0hr mean,sd")
print(KeGGInt0hr_Mean)
print(KeGGInt0hr_Sd)
KeGGdataMerge$Intensity.0hr_Log2MedNorm_Zscore <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) zscores(x[["Intensity.0hr_Log2MedNorm_Avg"]],KeGGInt0hr_Mean,KeGGInt0hr_Sd))
#6hr
KeGGInt6hr_Mean <- mean(as.numeric(KeGGdataMerge$Intensity.6hr_Log2MedNorm_Avg), na.rm=TRUE)
KeGGInt6hr_Sd <- sd(as.numeric(KeGGdataMerge$Intensity.6hr_Log2MedNorm_Avg), na.rm=TRUE)
print("6hr mean,sd")
print(KeGGInt6hr_Mean)
print(KeGGInt6hr_Sd)
KeGGdataMerge$Intensity.6hr_Log2MedNorm_Zscore <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) zscores(x[["Intensity.6hr_Log2MedNorm_Avg"]],KeGGInt6hr_Mean,KeGGInt6hr_Sd))
#8hr
KeGGInt8hr_Mean <- mean(as.numeric(KeGGdataMerge$Intensity.8hr_Log2MedNorm_Avg), na.rm=TRUE)
KeGGInt8hr_Sd <- sd(as.numeric(KeGGdataMerge$Intensity.8hr_Log2MedNorm_Avg), na.rm=TRUE)
print("8hr mean,sd")
print(KeGGInt8hr_Mean)
print(KeGGInt8hr_Sd)
KeGGdataMerge$Intensity.8hr_Log2MedNorm_Zscore <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) zscores(x[["Intensity.8hr_Log2MedNorm_Avg"]],KeGGInt8hr_Mean,KeGGInt8hr_Sd))
#10hr
KeGGInt10hr_Mean <- mean(as.numeric(KeGGdataMerge$Intensity.10hr_Log2MedNorm_Avg), na.rm=TRUE)
KeGGInt10hr_Sd <- sd(as.numeric(KeGGdataMerge$Intensity.10hr_Log2MedNorm_Avg), na.rm=TRUE)
print("10hr mean,sd")
print(KeGGInt10hr_Mean)
print(KeGGInt10hr_Sd)
KeGGdataMerge$Intensity.10hr_Log2MedNorm_Zscore <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) zscores(x[["Intensity.10hr_Log2MedNorm_Avg"]],KeGGInt10hr_Mean,KeGGInt10hr_Sd))
##end calculate z-scores


##calculate fold change of the non-normalized average intensity for each time point comparison
##use apply function to parse each row of the dataframe, foldchange function to calc fc
if ( FALSE )
{
	#0hr-6hr comp avg
	KeGGdataMerge$KeGG_Int_0hr6hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["Intensity.0hr_Log2MedNorm_Avg"]])),as.numeric(c(x[["Intensity.6hr_Log2MedNorm_Avg"]]))))
	#0hr-8hr comp avg
	KeGGdataMerge$KeGG_Int_0hr8hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["Intensity.0hr_Log2MedNorm_Avg"]])),as.numeric(c(x[["Intensity.8hr_Log2MedNorm_Avg"]]))))
	#0hr-10hr comp avg
	KeGGdataMerge$KeGG_Int_0hr10hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["Intensity.0hr_Log2MedNorm_Avg"]])),as.numeric(c(x[["Intensity.10hr_Log2MedNorm_Avg"]]))))
	#6hr-8hr comp avg
	KeGGdataMerge$KeGG_Int_6hr8hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["Intensity.6hr_Log2MedNorm_Avg"]])),as.numeric(c(x[["Intensity.8hr_Log2MedNorm_Avg"]]))))
	#8hr-10hr comp avg
	KeGGdataMerge$KeGG_Int_8hr10hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["Intensity.8hr_Log2MedNorm_Avg"]])),as.numeric(c(x[["Intensity.10hr_Log2MedNorm_Avg"]]))))
}
#0hr-6hr comp avg
KeGGdataMerge$KeGG_Int_0hr6hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) ifelse((x[["KeGGInt.0hr_NumReps"]]>=2 & x[["KeGGInt.6hr_NumReps"]]>=2), foldchange(as.numeric(c(x[["Intensity.0hr_Log2MedNorm_Avg"]])),as.numeric(c(x[["Intensity.6hr_Log2MedNorm_Avg"]]))), NA))
#0hr-8hr comp avg
KeGGdataMerge$KeGG_Int_0hr8hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) ifelse((x[["KeGGInt.0hr_NumReps"]]>=2 & x[["KeGGInt.8hr_NumReps"]]>=2), foldchange(as.numeric(c(x[["Intensity.0hr_Log2MedNorm_Avg"]])),as.numeric(c(x[["Intensity.8hr_Log2MedNorm_Avg"]]))), NA))
#0hr-10hr comp avg
KeGGdataMerge$KeGG_Int_0hr10hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) ifelse((x[["KeGGInt.0hr_NumReps"]]>=2 & x[["KeGGInt.10hr_NumReps"]]>=2), foldchange(as.numeric(c(x[["Intensity.0hr_Log2MedNorm_Avg"]])),as.numeric(c(x[["Intensity.10hr_Log2MedNorm_Avg"]]))), NA))
#6hr-8hr comp avg
KeGGdataMerge$KeGG_Int_6hr8hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) ifelse((x[["KeGGInt.6hr_NumReps"]]>=2 & x[["KeGGInt.8hr_NumReps"]]>=2), foldchange(as.numeric(c(x[["Intensity.6hr_Log2MedNorm_Avg"]])),as.numeric(c(x[["Intensity.8hr_Log2MedNorm_Avg"]]))), NA))
#8hr-10hr comp avg
KeGGdataMerge$KeGG_Int_8hr10hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) ifelse((x[["KeGGInt.8hr_NumReps"]]>=2 & x[["KeGGInt.10hr_NumReps"]]>=2), foldchange(as.numeric(c(x[["Intensity.8hr_Log2MedNorm_Avg"]])),as.numeric(c(x[["Intensity.10hr_Log2MedNorm_Avg"]]))), NA))
##calc fold change for individual replicates
#run2
KeGGdataMerge$KeGG_Int_Run2_0hr6hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd2"]])),as.numeric(c(x[["Intensity.6hr_Log2_MedNorm.Rd2"]]))))
KeGGdataMerge$KeGG_Int_Run2_0hr8hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd2"]])),as.numeric(c(x[["Intensity.8hr_Log2_MedNorm.Rd2"]]))))
KeGGdataMerge$KeGG_Int_Run2_0hr10hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd2"]])),as.numeric(c(x[["Intensity.10hr_Log2_MedNorm.Rd2"]]))))
#run3
KeGGdataMerge$KeGG_Int_Run3_0hr6hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd3"]])),as.numeric(c(x[["Intensity.6hr_Log2_MedNorm.Rd3"]]))))
KeGGdataMerge$KeGG_Int_Run3_0hr8hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd3"]])),as.numeric(c(x[["Intensity.8hr_Log2_MedNorm.Rd3"]]))))
KeGGdataMerge$KeGG_Int_Run3_0hr10hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd3"]])),as.numeric(c(x[["Intensity.10hr_Log2_MedNorm.Rd3"]]))))
#run4
KeGGdataMerge$KeGG_Int_Run4_0hr6hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd4"]])),as.numeric(c(x[["Intensity.6hr_Log2_MedNorm.Rd4"]]))))
KeGGdataMerge$KeGG_Int_Run4_0hr8hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd4"]])),as.numeric(c(x[["Intensity.8hr_Log2_MedNorm.Rd4"]]))))
KeGGdataMerge$KeGG_Int_Run4_0hr10hr_Fc <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) foldchange(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd4"]])),as.numeric(c(x[["Intensity.10hr_Log2_MedNorm.Rd4"]]))))
##end calculate fold change for time point comparison


##signficance tests for non-normalized intensity values compared for replicates across each time point comparisons
##use apply to parse each row of dataset, use sigtest function to calc t-test
##unpaired
#0hr-6hr unpaired
sigtestresults<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) sigtest(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd2"]],x[["Intensity.0hr_Log2_MedNorm.Rd3"]],x[["Intensity.0hr_Log2_MedNorm.Rd4"]])),as.numeric(c(x[["Intensity.6hr_Log2_MedNorm.Rd2"]],x[["Intensity.6hr_Log2_MedNorm.Rd3"]],x[["Intensity.6hr_Log2_MedNorm.Rd4"]])),'UNPAIRED'))
KeGGdataMerge$KeGG_Int_0hr6hr_unpairedTtestPval<-sigtestresults[2,]
KeGGdataMerge$KeGG_Int_0hr6hr_unpairedTtestPval_negLog10<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["KeGG_Int_0hr6hr_unpairedTtestPval"]])))
#0hr-8hr unpaired
sigtestresults<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) sigtest(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd2"]],x[["Intensity.0hr_Log2_MedNorm.Rd3"]],x[["Intensity.0hr_Log2_MedNorm.Rd4"]])),as.numeric(c(x[["Intensity.8hr_Log2_MedNorm.Rd2"]],x[["Intensity.8hr_Log2_MedNorm.Rd3"]],x[["Intensity.8hr_Log2_MedNorm.Rd4"]])),'UNPAIRED'))
KeGGdataMerge$KeGG_Int_0hr8hr_unpairedTtestPval<-sigtestresults[2,]
KeGGdataMerge$KeGG_Int_0hr8hr_unpairedTtestPval_negLog10<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["KeGG_Int_0hr8hr_unpairedTtestPval"]])))
#0hr-10hr unpaired
sigtestresults<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) sigtest(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd2"]],x[["Intensity.0hr_Log2_MedNorm.Rd3"]],x[["Intensity.0hr_Log2_MedNorm.Rd4"]])),as.numeric(c(x[["Intensity.10hr_Log2_MedNorm.Rd2"]],x[["Intensity.10hr_Log2_MedNorm.Rd3"]],x[["Intensity.10hr_Log2_MedNorm.Rd4"]])),'UNPAIRED'))
KeGGdataMerge$KeGG_Int_0hr10hr_unpairedTtestPval<-sigtestresults[2,]
KeGGdataMerge$KeGG_Int_0hr10hr_unpairedTtestPval_negLog10<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["KeGG_Int_0hr10hr_unpairedTtestPval"]])))
##paired
#0hr-6hr paired
sigtestresults<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) sigtest(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd2"]],x[["Intensity.0hr_Log2_MedNorm.Rd3"]],x[["Intensity.0hr_Log2_MedNorm.Rd4"]])),as.numeric(c(x[["Intensity.6hr_Log2_MedNorm.Rd2"]],x[["Intensity.6hr_Log2_MedNorm.Rd3"]],x[["Intensity.6hr_Log2_MedNorm.Rd4"]])),'PAIRED'))
KeGGdataMerge$KeGG_Int_0hr6hr_pairedTtestPval<-sigtestresults[2,]
KeGGdataMerge$KeGG_Int_0hr6hr_pairedTtestPval_negLog10<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["KeGG_Int_0hr6hr_pairedTtestPval"]])))
#0hr-8hr paired
sigtestresults<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) sigtest(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd2"]],x[["Intensity.0hr_Log2_MedNorm.Rd3"]],x[["Intensity.0hr_Log2_MedNorm.Rd4"]])),as.numeric(c(x[["Intensity.8hr_Log2_MedNorm.Rd2"]],x[["Intensity.8hr_Log2_MedNorm.Rd3"]],x[["Intensity.8hr_Log2_MedNorm.Rd4"]])),'PAIRED'))
KeGGdataMerge$KeGG_Int_0hr8hr_pairedTtestPval<-sigtestresults[2,]
KeGGdataMerge$KeGG_Int_0hr8hr_pairedTtestPval_negLog10<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["KeGG_Int_0hr8hr_pairedTtestPval"]])))
#0hr-10hr paired
sigtestresults<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) sigtest(as.numeric(c(x[["Intensity.0hr_Log2_MedNorm.Rd2"]],x[["Intensity.0hr_Log2_MedNorm.Rd3"]],x[["Intensity.0hr_Log2_MedNorm.Rd4"]])),as.numeric(c(x[["Intensity.10hr_Log2_MedNorm.Rd2"]],x[["Intensity.10hr_Log2_MedNorm.Rd3"]],x[["Intensity.10hr_Log2_MedNorm.Rd4"]])),'PAIRED'))
KeGGdataMerge$KeGG_Int_0hr10hr_pairedTtestPval<-sigtestresults[2,]
KeGGdataMerge$KeGG_Int_0hr10hr_pairedTtestPval_negLog10<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["KeGG_Int_0hr10hr_pairedTtestPval"]])))
##end signficance tests for time point comparisons


#####head(KeGGdataMerge)


##combine peptide to protein ubiquitination
##unique protein ids in dataset
proteinIDs<-unique(KeGGdataMerge$Protein, incompariables=FALSE)
proteinIDsNum<-length(proteinIDs)
#print(proteinIDsNum)
#head(proteinIDs)

##combine peptide quantification in to protein-based quantification
#initialize data frame to store protein-level information
KeGGdataMerge_ProteinData<-data.frame()

#evaluate all proteins identified in the combined dataset
#parse each unique protein id
#analyze only peptides with NumReps>=2 for 0hr or 10hr
for ( i in 1:length(proteinIDs) )
{

	#current protein
	tempprot<-proteinIDs[i]
	#print(tempprot)
	
	#peptide data for current protein
	tempdata_tmp<-subset(KeGGdataMerge, KeGGdataMerge$Protein==tempprot)
	##########tempdata<-subset(tempdata_tmp, (tempdata_tmp$KeGGInt.10hr_NumReps>=2 | tempdata_tmp$KeGGInt.0hr_NumReps>=2))
	tempdata<-subset(tempdata_tmp, ((tempdata_tmp$KeGGInt.10hr_NumReps>=2 & tempdata_tmp$KeGGInt.0hr_NumReps>=2) | (tempdata_tmp$KeGGInt.10hr_NumReps>=2 & tempdata_tmp$KeGGInt.0hr_NumReps==0) | (tempdata_tmp$KeGGInt.10hr_NumReps==0 & tempdata_tmp$KeGGInt.0hr_NumReps>=2)))
	unique0hr<-NA
	unique10hr<-NA
	if ( sum(tempdata$KeGGInt.0hr_NumReps)==0 ){
		unique10hr<-1
	}
	if ( sum(tempdata$KeGGInt.10hr_NumReps)==0 ){
		unique0hr<-1
	}
	
	#num of identified sites for current protein
	numsites<-dim(tempdata)[1]
	#print(numsites)
	
	#skip "0" sites, these proteins were identified in a different study
	#the number of sites are counted by those that have an HLcount
	if ( numsites>0 )
	{
	
		#populate protein data with protein ID
		KeGGdataMerge_ProteinData[tempprot,"Protein"]<-tempprot

		#populate protein data with number of sites for current protein
		KeGGdataMerge_ProteinData[tempprot,"numsites"]<-numsites
		
		#gene ids in dataset
		geneID<-tempdata[1,"Gene.names"]
		KeGGdataMerge_ProteinData[tempprot,"Gene.names"]<-geneID
		geneIDsplit<-strsplit(geneID, ";")
		firstGeneID<-ifelse(length(geneIDsplit[[1]])==1, geneID, geneIDsplit[[1]][1])
		KeGGdataMerge_ProteinData[tempprot,"GeneID"]<-firstGeneID
		#protein names in dataset
		protName<-tempdata[1,"Protein.names"]
		KeGGdataMerge_ProteinData[tempprot,"Protein.names"]<-protName

		#indicate unique KeGG proteins
		KeGGdataMerge_ProteinData[tempprot,"Unique0hr"]<-unique0hr
		KeGGdataMerge_ProteinData[tempprot,"Unique10hr"]<-unique10hr
				
		#add fold change value for unique peptides (NA fold change)
		#if ( geneID == "RAD50") {
		#	print(geneID)
		#	print(tempdata$KeGG_Int_0hr10hr_Fc)
		#}
		NewFc <- function(fc, num10reps, num0reps)
		{
		      #print(fc)
			  #print(num10reps)
			  #print(num0reps)
			  newfc <- NA
			  if ( is.na(fc) ){
			       if ( num10reps>=2 & num0reps==0 ){
						newfc <- 7
				   }else{
				        if ( num10reps==0 & num0reps>=2 ){
						     newfc <- -7
						}else{
							stop("error")
						}
				   }
			  }else{
			       newfc <- fc
			  }
		      return(newfc)
		}
		tempdata$KeGG_Int_0hr10hr_Fc <- apply(tempdata, MARGIN=1, FUN=function(x) NewFc(x[["KeGG_Int_0hr10hr_Fc"]], x[["KeGGInt.10hr_NumReps"]], x[["KeGGInt.0hr_NumReps"]]))
	
		#average intensity value of 0hr and 10hr for each peptide
		tempdata$avgintensity<-apply(tempdata, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["Intensity.0hr_Log2MedNorm_Avg"]],x[["Intensity.10hr_Log2MedNorm_Avg"]])),na.rm=TRUE))

		#if ( FALSE )
		#{
		##maximum intensity value at 10hr
		#maxinten10hr<-max(tempdata[,"Intensity.10hr_Log2MedNorm_Avg"], na.rm=TRUE)
		#KeGGdataMerge_ProteinData[tempprot,"maxinten10hr"]<-maxinten10hr
		##peptide with the maximum intensity
		#maxintenpep<-row.names(tempdata[which(tempdata[,"Intensity.10hr_Log2MedNorm_Avg"]==maxinten10hr), ])[1]
		##populate protein data frame
		#KeGGdataMerge_ProteinData[tempprot,"maxintenpep"]<-maxintenpep
		##print(maxintenpep)
		#
		#####maxinten0hr<-tempdata[maxintenpep,"Intensity.0hr_Log2MedNorm_Avg"]
		#####KeGGdataMerge_ProteinData[tempprot,"maxinten0hr"]<-maxinten0hr
		#
		##max intensity peptide normLHratio
		#maxintenratio<-tempdata[maxintenpep,"KeGG_Int_0hr10hr_Fc"]
		##print(maxintenratio)
		##populate protein data frame
		#KeGGdataMerge_ProteinData[tempprot,"maxintenpepratio"]<-maxintenratio
		#}

		#intensity weighted mean based on base10 values
		#Log2 intensity and fold change must be converted to base10
		#convert weighted mean fold change back to Log2
		##########tempdata$base10intensity<-apply(tempdata, MARGIN=1, FUN=function(x) InvLog(x[["Intensity.10hr_Log2MedNorm_Avg"]],2))
		tempdata$base10intensity<-apply(tempdata, MARGIN=1, FUN=function(x) InvLog(x[["avgintensity"]],2))
		tempdata$base10foldchange<-apply(tempdata, MARGIN=1, FUN=function(x) InvLog(x[["KeGG_Int_0hr10hr_Fc"]],2))
		wtfcdata<-tempdata[,c("base10foldchange","base10intensity")]
		colnames(wtfcdata) <- c("foldchange", "intensity")
		intensityweightedmeanbase10<-CalcKeGGIntWtFc(wtfcdata)
		intensityweightedmeanbase2<-log(intensityweightedmeanbase10,2)
		#print(intensityweightedmeanbase2)
		#populate protein data frame
		KeGGdataMerge_ProteinData[tempprot,"intensityweightedmeanratio"]<-intensityweightedmeanbase2
	
	}
	
} #end protein site loop

#dim(KeGGdataMerge_ProteinData)
#head(KeGGdataMerge_ProteinData)


KeGG_0hr10hr_ProtFc_Avg <- mean(KeGGdataMerge_ProteinData$intensityweightedmeanratio, na.rm=TRUE)
KeGG_0hr10hr_ProtFc_Sd <- sd(KeGGdataMerge_ProteinData$intensityweightedmeanratio, na.rm=TRUE)
print("KeGG_0hr10hr_ProtFc_Avg")
print(KeGG_0hr10hr_ProtFc_Avg)
print("KeGG_0hr10hr_ProtFc_Sd")
print(KeGG_0hr10hr_ProtFc_Sd)


#read whole cell proteome data
##########WCPdata <- read.table("/Users/dybasj/JoeWork/Ad5_E1B55Kubiquitin/E1B55KE4orf6transfectionWCPData/WCP_0-10hr_mbr_AnalyzedData.txt", sep="\t", quote="\"", row.names=1, header=TRUE)
WCPdata <- read.table("/Users/dybasj/JoeWork/Ad5_E1B55Kubiquitin/E1B55KE4orf6transfectionWCPData/WCP_0-6-8-10hr_mbr_AnalyzedData.txt", sep="\t", quote="\"", row.names=1, header=TRUE)
WCPdata$firstprot<-row.names(WCPdata)
#head(WCPdata)
dim(WCPdata)

##count number of WCP proteins with >= 2 replicates identified
#0hr
WCPprots_0hr_2reps_ProtId <- row.names(subset(WCPdata, WCPdata$iBAQ_0hr_NumReps>=2))
WCPprots_0hr_2reps_ProtId_Num <- length(WCPprots_0hr_2reps_ProtId)
print("WCPprots_0hr_2reps_ProtId_Num")
print(WCPprots_0hr_2reps_ProtId_Num)
#6hr
WCPprots_6hr_2reps_ProtId <- row.names(subset(WCPdata, WCPdata$iBAQ_6hr_NumReps>=2))
WCPprots_6hr_2reps_ProtId_Num <- length(WCPprots_6hr_2reps_ProtId)
print("WCPprots_6hr_2reps_ProtId_Num")
print(WCPprots_6hr_2reps_ProtId_Num)
#8hr
WCPprots_8hr_2reps_ProtId <- row.names(subset(WCPdata, WCPdata$iBAQ_8hr_NumReps>=2))
WCPprots_8hr_2reps_ProtId_Num <- length(WCPprots_8hr_2reps_ProtId)
print("WCPprots_8hr_2reps_ProtId_Num")
print(WCPprots_8hr_2reps_ProtId_Num)
#10hr
WCPprots_10hr_2reps_ProtId <- row.names(subset(WCPdata, WCPdata$iBAQ_10hr_NumReps>=2))
WCPprots_10hr_2reps_ProtId_Num <- length(WCPprots_10hr_2reps_ProtId)
print("WCPprots_10hr_2reps_ProtId_Num")
print(WCPprots_10hr_2reps_ProtId_Num)


#filter WCP data
WCPdata_Filtered <- subset(WCPdata, ((WCPdata$iBAQ_0hr_NumReps>=2 & WCPdata$iBAQ_10hr_NumReps>=2) | (WCPdata$iBAQ_0hr_NumReps>=2 & WCPdata$iBAQ_10hr_NumReps==0) | (WCPdata$iBAQ_0hr_NumReps==0 & WCPdata$iBAQ_10hr_NumReps>=2)))
WCPdata <- WCPdata_Filtered
dim(WCPdata)


##make the WCP first protein to protein group lookup table
WCPlookup<-data.frame()
for ( i in 1:nrow(WCPdata) )
{
	WCPfirstprot <- as.character(WCPdata[i,"firstprot"])
	WCPprotgroupI <- as.character(WCPdata[i,"Protein.IDs.x"])
	WCPprotgroupII <- as.character(WCPdata[i,"Protein.IDs.y"])
	if ( !is.na(WCPprotgroupI) ){
		WCPprotgroupvecI <- unlist(strsplit(WCPprotgroupI, ";"))
	}else{
		WCPprotgroupvecI
	}
	if ( !is.na(WCPprotgroupII) ){
		WCPprotgroupvecII <- unlist(strsplit(WCPprotgroupII, ";"))
	}else{
		WCPprotgroupvecII
	}
	WCPprotgroupvec <- c(WCPprotgroupvecI,WCPprotgroupvecII)
	for ( j in 1:length(WCPprotgroupvec) )
	{
		grpprot <- WCPprotgroupvec[j]
		tmpdata<-data.frame()
		tmpdata[grpprot,"firstprot"]<-WCPfirstprot
		WCPlookup <- rbind(WCPlookup, tmpdata)
	}
}
#head(WCPlookup)


##find KeGG-WCP correspondence
FindWCPProt <- function(KeGGprot, WCPlookup)
{

	WCProwprot<-NA
	WCProwprot<-WCPlookup[KeGGprot,"firstprot"]
	return(WCProwprot)

}
print("get corr prot...")
KeGGdataMerge$WCPcorrProt<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) FindWCPProt(x[["Protein"]], WCPlookup))
#head(KeGGdataMerge)
KeGGdataMerge_ProteinData$WCPcorrProt<-apply(KeGGdataMerge_ProteinData, MARGIN=1, FUN=function(x) FindWCPProt(x[["Protein"]], WCPlookup))
#

KeGGpepWithWcpData <- subset(KeGGdataMerge, !is.na(KeGGdataMerge$WCPcorrProt))
KeGGpepWithWcp <- unique(KeGGpepWithWcpData$Protein)
KeGGpepWithWcp_Num <- length(KeGGpepWithWcp)
print("KeGGpepWithWcp_Num")
print(KeGGpepWithWcp_Num)
#
KeGGpepWithOutWcpData <- subset(KeGGdataMerge, is.na(KeGGdataMerge$WCPcorrProt))
KeGGpepWithOutWcp <- unique(KeGGpepWithOutWcpData$Protein)
KeGGpepWithOutWcp_Num <- length(KeGGpepWithOutWcp)
print("KeGGpepWithOutWcp_Num")
print(KeGGpepWithOutWcp_Num)
#
#
KeGGprotWithWcpData <- subset(KeGGdataMerge_ProteinData, !is.na(KeGGdataMerge_ProteinData$WCPcorrProt))
KeGGprotWithWcp <- unique(KeGGprotWithWcpData$Protein)
KeGGprotWithWcp_Num <- length(KeGGprotWithWcp)
print("KeGGprotWithWcp_Num")
print(KeGGprotWithWcp_Num)
#
KeGGprotWithOutWcpData <- subset(KeGGdataMerge_ProteinData, is.na(KeGGdataMerge_ProteinData$WCPcorrProt))
KeGGprotWithOutWcp <- unique(KeGGprotWithOutWcpData$Protein)
KeGGprotWithOutWcp_Num <- length(KeGGprotWithOutWcp)
print("KeGGprotWithOutWcp_Num")
print(KeGGprotWithOutWcp_Num)
#


##incorporate relevant WCP information for each KeGG peptide
#corresponding WCP fold change
WCPFoldChange <- function(WCPprot, WCPdata, WCPavgfoldchangeId){
	WCPfoldchange <- WCPdata[WCPprot, WCPavgfoldchangeId]
	return(WCPfoldchange)
}
#corresponding WCP significance value
WCPSigVal <- function(WCPprot, WCPdata, WCPpvalId){
	WCPpval <- WCPdata[WCPprot, WCPpvalId]
	return(WCPpval)
}
#normalize KeGG fold change by WCP
NormKeGGFoldChange <- function(KeGGFoldChange, WCPprot, WCPdata, WCPfoldchangeId)
{

	WCPfoldchange <- WCPdata[WCPprot, WCPfoldchangeId]
	NormKeGGfoldchange <- NA
	#if ( is.na(WCPfoldchange) ){
	#	WCPfoldchange <- WCPdata[WCPprot, "iBAQ_0hr10hr_Log2MedNorm_Avg_Fc"]
	#}
	if ( is.na(KeGGFoldChange) | is.na(WCPfoldchange) ){
		NormKeGGfoldchange <- NA
	}else{
		NormKeGGfoldchange <- as.numeric(KeGGFoldChange) - as.numeric(WCPfoldchange)
	}
	return(NormKeGGfoldchange)
}
#populate KeGG data frame with WCP average fold change
KeGGdataMerge$WCP_iBAQ_0hr6hr_Fc_Avg <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) WCPFoldChange(x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr6hr_Log2MedNorm_Avg_Fc"))
KeGGdataMerge$WCP_iBAQ_0hr8hr_Fc_Avg <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) WCPFoldChange(x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr8hr_Log2MedNorm_Avg_Fc"))
KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) WCPFoldChange(x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr10hr_Log2MedNorm_Avg_Fc"))
#populate KeGG data frame with WCP fold change p-value
KeGGdataMerge$WCP_iBAQ_0hr6hr_pval <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) WCPSigVal(x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr6hr_Log2MedNorm_pval"))
KeGGdataMerge$WCP_iBAQ_0hr8hr_pval <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) WCPSigVal(x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr8hr_Log2MedNorm_pval"))
KeGGdataMerge$WCP_iBAQ_0hr10hr_pval <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) WCPSigVal(x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr10hr_Log2MedNorm_pval"))


#populate KeGGprotein data frame with WCP average fold change
KeGGdataMerge_ProteinData$WCP_iBAQ_0hr6hr_Fc_Avg <- apply(KeGGdataMerge_ProteinData, MARGIN=1, FUN=function(x) WCPFoldChange(x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr6hr_Log2MedNorm_Avg_Fc"))
KeGGdataMerge_ProteinData$WCP_iBAQ_0hr8hr_Fc_Avg <- apply(KeGGdataMerge_ProteinData, MARGIN=1, FUN=function(x) WCPFoldChange(x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr8hr_Log2MedNorm_Avg_Fc"))
KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg <- apply(KeGGdataMerge_ProteinData, MARGIN=1, FUN=function(x) WCPFoldChange(x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr10hr_Log2MedNorm_Avg_Fc"))
#populate KeGG data frame with WCP fold change p-value
KeGGdataMerge_ProteinData$WCP_iBAQ_0hr6hr_pval <- apply(KeGGdataMerge_ProteinData, MARGIN=1, FUN=function(x) WCPSigVal(x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr6hr_Log2MedNorm_pval"))
KeGGdataMerge_ProteinData$WCP_iBAQ_0hr8hr_pval <- apply(KeGGdataMerge_ProteinData, MARGIN=1, FUN=function(x) WCPSigVal(x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr8hr_Log2MedNorm_pval"))
KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_pval <- apply(KeGGdataMerge_ProteinData, MARGIN=1, FUN=function(x) WCPSigVal(x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr10hr_Log2MedNorm_pval"))


##normalize KeGG change by WCP change
##normalize for each replicate and for average fold changes
#0hr6hr
KeGGdataMerge$KeGG_Int_Run2_0hr6hr_Fc_Norm <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) NormKeGGFoldChange(x[["KeGG_Int_Run2_0hr6hr_Fc"]],x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr6hr_Log2MedNorm_R1_Fc"))
KeGGdataMerge$KeGG_Int_Run3_0hr6hr_Fc_Norm <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) NormKeGGFoldChange(x[["KeGG_Int_Run3_0hr6hr_Fc"]],x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr6hr_Log2MedNorm_R2_Fc"))
KeGGdataMerge$KeGG_Int_Run4_0hr6hr_Fc_Norm <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) NormKeGGFoldChange(x[["KeGG_Int_Run4_0hr6hr_Fc"]],x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr6hr_Log2MedNorm_R3_Fc"))
KeGGdataMerge$KeGG_Int_0hr6hr_FcNorm_Avg <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["KeGG_Int_Run2_0hr6hr_Fc_Norm"]],x[["KeGG_Int_Run3_0hr6hr_Fc_Norm"]],x[["KeGG_Int_Run4_0hr6hr_Fc_Norm"]])),na.rm=TRUE))
#0hr8hr
KeGGdataMerge$KeGG_Int_Run2_0hr8hr_Fc_Norm <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) NormKeGGFoldChange(x[["KeGG_Int_Run2_0hr8hr_Fc"]],x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr8hr_Log2MedNorm_R1_Fc"))
KeGGdataMerge$KeGG_Int_Run3_0hr8hr_Fc_Norm <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) NormKeGGFoldChange(x[["KeGG_Int_Run3_0hr8hr_Fc"]],x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr8hr_Log2MedNorm_R2_Fc"))
KeGGdataMerge$KeGG_Int_Run4_0hr8hr_Fc_Norm <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) NormKeGGFoldChange(x[["KeGG_Int_Run4_0hr8hr_Fc"]],x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr8hr_Log2MedNorm_R3_Fc"))
KeGGdataMerge$KeGG_Int_0hr8hr_FcNorm_Avg <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["KeGG_Int_Run2_0hr8hr_Fc_Norm"]],x[["KeGG_Int_Run3_0hr8hr_Fc_Norm"]],x[["KeGG_Int_Run4_0hr8hr_Fc_Norm"]])),na.rm=TRUE))
#0hr10hr
KeGGdataMerge$KeGG_Int_Run2_0hr10hr_Fc_Norm <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) NormKeGGFoldChange(x[["KeGG_Int_Run2_0hr10hr_Fc"]],x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr10hr_Log2MedNorm_R1_Fc"))
KeGGdataMerge$KeGG_Int_Run3_0hr10hr_Fc_Norm <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) NormKeGGFoldChange(x[["KeGG_Int_Run3_0hr10hr_Fc"]],x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr10hr_Log2MedNorm_R2_Fc"))
KeGGdataMerge$KeGG_Int_Run4_0hr10hr_Fc_Norm <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) NormKeGGFoldChange(x[["KeGG_Int_Run4_0hr10hr_Fc"]],x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr10hr_Log2MedNorm_R3_Fc"))
KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg <- apply(KeGGdataMerge, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["KeGG_Int_Run2_0hr10hr_Fc_Norm"]],x[["KeGG_Int_Run3_0hr10hr_Fc_Norm"]],x[["KeGG_Int_Run4_0hr10hr_Fc_Norm"]])),na.rm=TRUE))


##normalize KeGG protein change by WCP change
KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm <- apply(KeGGdataMerge_ProteinData, MARGIN=1, FUN=function(x) NormKeGGFoldChange(x[["intensityweightedmeanratio"]],x[["WCPcorrProt"]],WCPdata,"iBAQ_0hr10hr_Log2MedNorm_Avg_Fc"))


##perform significance test on normalized fold change for replicates compared to 0 (no change)
##one-sided test of fold change for deviation from 0 (no change)
#0-6hr
KeGGdataMerge$KeGG_Int_0hr6hr_FcNorm_onesidePval<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) onesetsigtest(as.numeric(c(x[["KeGG_Int_Run2_0hr6hr_Fc_Norm"]],x[["KeGG_Int_Run3_0hr6hr_Fc_Norm"]],x[["KeGG_Int_Run4_0hr6hr_Fc_Norm"]]))))
KeGGdataMerge$KeGG_Int_0hr6hr_FcNorm_onesidePval_negLog10<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["KeGG_Int_0hr6hr_FcNorm_onesidePval"]])))
#0-8hr
KeGGdataMerge$KeGG_Int_0hr8hr_FcNorm_onesidePval<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) onesetsigtest(as.numeric(c(x[["KeGG_Int_Run2_0hr8hr_Fc_Norm"]],x[["KeGG_Int_Run3_0hr8hr_Fc_Norm"]],x[["KeGG_Int_Run4_0hr8hr_Fc_Norm"]]))))
KeGGdataMerge$KeGG_Int_0hr8hr_FcNorm_onesidePval_negLog10<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["KeGG_Int_0hr8hr_FcNorm_onesidePval"]])))
#0-10hr
KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_onesidePval<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) onesetsigtest(as.numeric(c(x[["KeGG_Int_Run2_0hr10hr_Fc_Norm"]],x[["KeGG_Int_Run3_0hr10hr_Fc_Norm"]],x[["KeGG_Int_Run4_0hr10hr_Fc_Norm"]]))))
KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_onesidePval_negLog10<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["KeGG_Int_0hr10hr_FcNorm_onesidePval"]])))
##perform significance test on un-normalized fold change for 0-10hr
#0-6hr
KeGGdataMerge$KeGG_Int_0hr6hr_FcRaw_onesidePval<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) onesetsigtest(as.numeric(c(x[["KeGG_Int_Run2_0hr6hr_Fc"]],x[["KeGG_Int_Run3_0hr6hr_Fc"]],x[["KeGG_Int_Run4_0hr6hr_Fc"]]))))
#0-8hr
KeGGdataMerge$KeGG_Int_0hr8hr_FcRaw_onesidePval<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) onesetsigtest(as.numeric(c(x[["KeGG_Int_Run2_0hr8hr_Fc"]],x[["KeGG_Int_Run3_0hr8hr_Fc"]],x[["KeGG_Int_Run4_0hr8hr_Fc"]]))))
#0-10hr
KeGGdataMerge$KeGG_Int_0hr10hr_FcRaw_onesidePval<-apply(KeGGdataMerge, MARGIN=1, FUN=function(x) onesetsigtest(as.numeric(c(x[["KeGG_Int_Run2_0hr10hr_Fc"]],x[["KeGG_Int_Run3_0hr10hr_Fc"]],x[["KeGG_Int_Run4_0hr10hr_Fc"]]))))


##write table of complete dataset
write.table(KeGGdataMerge, file="KeGGdata_0hr6hr8hr10hr_Run2Run3Run4_AllData_SupTable1.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
write.table(KeGGdataMerge_ProteinData, file="KeGGdata_0hr6hr8hr10hr_ProteinBasedData_AllData_SupTable2.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#head(KeGGdataMerge)


##WCP proteins that are identified or not identified by number of reps quantified
#0hr
WCPprots_0hr_Id_NumReps <- row.names(subset(WCPdata, WCPdata$iBAQ_0hr_NumReps>=2))
WCPprots_0hr_Id_NumReps_Num <- length(WCPprots_0hr_Id_NumReps)
print("WCPprots_0hr_Id_NumReps_Num")
print(WCPprots_0hr_Id_NumReps_Num)
WCPprots_0hr_NotId_NumReps <- row.names(subset(WCPdata, WCPdata$iBAQ_0hr_NumReps==0))
WCPprots_0hr_NotId_NumReps_Num <- length(WCPprots_0hr_NotId_NumReps)
print("WCPprots_0hr_NotId_NumReps_Num")
print(WCPprots_0hr_NotId_NumReps_Num)
#6hr
WCPprots_6hr_Id_NumReps <- row.names(subset(WCPdata, WCPdata$iBAQ_6hr_NumReps>=2))
WCPprots_6hr_Id_NumReps_Num <- length(WCPprots_6hr_Id_NumReps)
print("WCPprots_6hr_Id_NumReps_Num")
print(WCPprots_6hr_Id_NumReps_Num)
WCPprots_6hr_NotId_NumReps <- row.names(subset(WCPdata, WCPdata$iBAQ_6hr_NumReps==0))
WCPprots_6hr_NotId_NumReps_Num <- length(WCPprots_6hr_NotId_NumReps)
print("WCPprots_6hr_NotId_NumReps_Num")
print(WCPprots_6hr_NotId_NumReps_Num)
#8hr
WCPprots_8hr_Id_NumReps <- row.names(subset(WCPdata, WCPdata$iBAQ_8hr_NumReps>=2))
WCPprots_8hr_Id_NumReps_Num <- length(WCPprots_8hr_Id_NumReps)
print("WCPprots_8hr_Id_NumReps_Num")
print(WCPprots_8hr_Id_NumReps_Num)
WCPprots_8hr_NotId_NumReps <- row.names(subset(WCPdata, WCPdata$iBAQ_8hr_NumReps==0))
WCPprots_8hr_NotId_NumReps_Num <- length(WCPprots_8hr_NotId_NumReps)
print("WCPprots_8hr_NotId_NumReps_Num")
print(WCPprots_8hr_NotId_NumReps_Num)
#10hr
WCPprots_10hr_Id_NumReps <- row.names(subset(WCPdata, WCPdata$iBAQ_10hr_NumReps>=2))
WCPprots_10hr_Id_NumReps_Num <- length(WCPprots_10hr_Id_NumReps)
print("WCPprots_10hr_Id_NumReps_Num")
print(WCPprots_10hr_Id_NumReps_Num)
WCPprots_10hr_NotId_NumReps <- row.names(subset(WCPdata, WCPdata$iBAQ_10hr_NumReps==0))
WCPprots_10hr_NotId_NumReps_Num <- length(WCPprots_10hr_NotId_NumReps)
print("WCPprots_10hr_NotId_NumReps_Num")
print(WCPprots_10hr_NotId_NumReps_Num)


##proteins that are increased or decreased in total protein abundance based on wcp fold change mean+sd, fold change threshold, p-value or unique
wcpfcthreshold=1
wcppvalthreshold<-0.05
#6hr mean
WCP_0hr6hr_Data <- subset(WCPdata, (WCPdata$iBAQ_0hr_NumReps>=2 & WCPdata$iBAQ_6hr_NumReps>=2))
WCP_0hr6hr_Fc_Avg <- mean(WCP_0hr6hr_Data$iBAQ_0hr6hr_Log2MedNorm_Avg_Fc, na.rm=TRUE)
WCP_0hr6hr_Fc_Sd <- sd(WCP_0hr6hr_Data$iBAQ_0hr6hr_Log2MedNorm_Avg_Fc, na.rm=TRUE)
print("WCP 6hr mean, sd")
print(WCP_0hr6hr_Fc_Avg)
print(WCP_0hr6hr_Fc_Sd)
WCPprots_6hr_Inc_MeanSd <- row.names(subset(WCP_0hr6hr_Data, WCP_0hr6hr_Data$iBAQ_0hr6hr_Log2MedNorm_Avg_Fc>=WCP_0hr6hr_Fc_Avg+WCP_0hr6hr_Fc_Sd))
WCPprots_6hr_Inc_MeanSd_Num <- length(WCPprots_6hr_Inc_MeanSd)
print("WCPprots_6hr_Inc_MeanSd_Num")
print(WCPprots_6hr_Inc_MeanSd_Num)
WCPprots_6hr_Dec_MeanSd <- row.names(subset(WCP_0hr6hr_Data, WCP_0hr6hr_Data$iBAQ_0hr6hr_Log2MedNorm_Avg_Fc<=WCP_0hr6hr_Fc_Avg-WCP_0hr6hr_Fc_Sd))
WCPprots_6hr_Dec_MeanSd_Num <- length(WCPprots_6hr_Dec_MeanSd)
print("WCPprots_6hr_Dec_MeanSd_Num")
print(WCPprots_6hr_Dec_MeanSd_Num)
WCPprots_6hr_Unch_MeanSd <- row.names(subset(WCP_0hr6hr_Data, (WCP_0hr6hr_Data$iBAQ_0hr6hr_Log2MedNorm_Avg_Fc>WCP_0hr6hr_Fc_Avg-WCP_0hr6hr_Fc_Sd & WCP_0hr6hr_Data$iBAQ_0hr6hr_Log2MedNorm_Avg_Fc<WCP_0hr6hr_Fc_Avg+WCP_0hr6hr_Fc_Sd)))
WCPprots_6hr_Unch_MeanSd_Num <- length(WCPprots_6hr_Unch_MeanSd)
print("WCPprots_6hr_Unch_MeanSd_Num")
print(WCPprots_6hr_Unch_MeanSd_Num)
#6hr thresh
WCPprots_6hr_Inc_Thresh <- row.names(subset(WCP_0hr6hr_Data, WCP_0hr6hr_Data$iBAQ_0hr6hr_Log2MedNorm_Avg_Fc >= wcpfcthreshold))
WCPprots_6hr_Inc_Thresh_Num <- length(WCPprots_6hr_Inc_Thresh)
print("WCPprots_6hr_Inc_Thresh_Num")
print(WCPprots_6hr_Inc_Thresh_Num)
WCPprots_6hr_Dec_Thresh <- row.names(subset(WCP_0hr6hr_Data, WCP_0hr6hr_Data$iBAQ_0hr6hr_Log2MedNorm_Avg_Fc <= (-1*wcpfcthreshold)))
WCPprots_6hr_Dec_Thresh_Num <- length(WCPprots_6hr_Dec_Thresh)
print("WCPprots_6hr_Dec_Thresh_Num")
print(WCPprots_6hr_Dec_Thresh_Num)
WCPprots_6hr_Unch_Thresh <- row.names(subset(WCP_0hr6hr_Data, (WCP_0hr6hr_Data$iBAQ_0hr6hr_Log2MedNorm_Avg_Fc > (-1*wcpfcthreshold) & WCP_0hr6hr_Data$iBAQ_0hr6hr_Log2MedNorm_Avg_Fc < wcpfcthreshold)))
WCPprots_6hr_Unch_Thresh_Num <- length(WCPprots_6hr_Unch_Thresh)
print("WCPprots_6hr_Unch_Thresh_Num")
print(WCPprots_6hr_Unch_Thresh_Num)
#6hr p-value
WCPprots_6hr_Inc_pval <- row.names(subset(WCP_0hr6hr_Data, (WCP_0hr6hr_Data$iBAQ_0hr6hr_Log2MedNorm_pval < wcppvalthreshold & WCP_0hr6hr_Data$iBAQ_0hr6hr_Log2MedNorm_Avg_Fc > 0)))
WCPprots_6hr_Inc_pval_Num <- length(WCPprots_6hr_Inc_pval)
print("WCPprots_6hr_Inc_pval_Num")
print(WCPprots_6hr_Inc_pval_Num)
WCPprots_6hr_Dec_pval <- row.names(subset(WCP_0hr6hr_Data, (WCP_0hr6hr_Data$iBAQ_0hr6hr_Log2MedNorm_pval < wcppvalthreshold & WCP_0hr6hr_Data$iBAQ_0hr6hr_Log2MedNorm_Avg_Fc < 0)))
WCPprots_6hr_Dec_pval_Num <- length(WCPprots_6hr_Dec_pval)
print("WCPprots_6hr_Dec_pval_Num")
print(WCPprots_6hr_Dec_pval_Num)
WCPprots_6hr_Unch_pval <- row.names(subset(WCP_0hr6hr_Data, (WCP_0hr6hr_Data$iBAQ_0hr6hr_Log2MedNorm_pval > wcppvalthreshold)))
WCPprots_6hr_Unch_pval_Num <- length(WCPprots_6hr_Unch_pval)
print("WCPprots_6hr_Unch_pval_Num")
print(WCPprots_6hr_Unch_pval_Num)
#6hr unique
WCPprots_6hr_Inc_unique <- row.names(subset(WCPdata,  (WCPdata$iBAQ_6hr_NumReps>=2 &  WCPdata$iBAQ_0hr_NumReps==0)))
WCPprots_6hr_Inc_unique_Num <- length(WCPprots_6hr_Inc_unique)
print("WCPprots_6hr_Inc_unique_Num")
print(WCPprots_6hr_Inc_unique_Num)
WCPprots_6hr_Dec_unique <- row.names(subset(WCPdata,  (WCPdata$iBAQ_6hr_NumReps==0 &  WCPdata$iBAQ_0hr_NumReps>=2)))
WCPprots_6hr_Dec_unique_Num <- length(WCPprots_6hr_Dec_unique)
print("WCPprots_6hr_Dec_unique_Num")
print(WCPprots_6hr_Dec_unique_Num)
#
#8hr mean
WCP_0hr8hr_Data <- subset(WCPdata, (WCPdata$iBAQ_0hr_NumReps>=2 & WCPdata$iBAQ_8hr_NumReps>=2))
WCP_0hr8hr_Fc_Avg <- mean(WCP_0hr8hr_Data$iBAQ_0hr8hr_Log2MedNorm_Avg_Fc, na.rm=TRUE)
WCP_0hr8hr_Fc_Sd <- sd(WCP_0hr8hr_Data$iBAQ_0hr8hr_Log2MedNorm_Avg_Fc, na.rm=TRUE)
print("WCP 8hr mean, sd")
print(WCP_0hr8hr_Fc_Avg)
print(WCP_0hr8hr_Fc_Sd)
WCPprots_8hr_Inc_MeanSd <- row.names(subset(WCP_0hr8hr_Data, WCP_0hr8hr_Data$iBAQ_0hr8hr_Log2MedNorm_Avg_Fc>=WCP_0hr8hr_Fc_Avg+WCP_0hr8hr_Fc_Sd))
WCPprots_8hr_Inc_MeanSd_Num <- length(WCPprots_8hr_Inc_MeanSd)
print("WCPprots_8hr_Inc_MeanSd_Num")
print(WCPprots_8hr_Inc_MeanSd_Num)
WCPprots_8hr_Dec_MeanSd <- row.names(subset(WCP_0hr8hr_Data, WCP_0hr8hr_Data$iBAQ_0hr8hr_Log2MedNorm_Avg_Fc<=WCP_0hr8hr_Fc_Avg-WCP_0hr8hr_Fc_Sd))
WCPprots_8hr_Dec_MeanSd_Num <- length(WCPprots_8hr_Dec_MeanSd)
print("WCPprots_8hr_Dec_MeanSd_Num")
print(WCPprots_8hr_Dec_MeanSd_Num)
WCPprots_8hr_Unch_MeanSd <- row.names(subset(WCP_0hr8hr_Data, (WCP_0hr8hr_Data$iBAQ_0hr8hr_Log2MedNorm_Avg_Fc>WCP_0hr8hr_Fc_Avg-WCP_0hr8hr_Fc_Sd & WCP_0hr8hr_Data$iBAQ_0hr8hr_Log2MedNorm_Avg_Fc<WCP_0hr8hr_Fc_Avg+WCP_0hr8hr_Fc_Sd)))
WCPprots_8hr_Unch_MeanSd_Num <- length(WCPprots_8hr_Unch_MeanSd)
print("WCPprots_8hr_Unch_MeanSd_Num")
print(WCPprots_8hr_Unch_MeanSd_Num)
#8hr thresh
WCPprots_8hr_Inc_Thresh <- row.names(subset(WCP_0hr8hr_Data, WCP_0hr8hr_Data$iBAQ_0hr8hr_Log2MedNorm_Avg_Fc >= wcpfcthreshold))
WCPprots_8hr_Inc_Thresh_Num <- length(WCPprots_8hr_Inc_Thresh)
print("WCPprots_8hr_Inc_Thresh_Num")
print(WCPprots_8hr_Inc_Thresh_Num)
WCPprots_8hr_Dec_Thresh <- row.names(subset(WCP_0hr8hr_Data, WCP_0hr8hr_Data$iBAQ_0hr8hr_Log2MedNorm_Avg_Fc <= (-1*wcpfcthreshold)))
WCPprots_8hr_Dec_Thresh_Num <- length(WCPprots_8hr_Dec_Thresh)
print("WCPprots_8hr_Dec_Thresh_Num")
print(WCPprots_8hr_Dec_Thresh_Num)
WCPprots_8hr_Unch_Thresh <- row.names(subset(WCP_0hr8hr_Data, (WCP_0hr8hr_Data$iBAQ_0hr8hr_Log2MedNorm_Avg_Fc > (-1*wcpfcthreshold) & WCP_0hr8hr_Data$iBAQ_0hr8hr_Log2MedNorm_Avg_Fc < wcpfcthreshold)))
WCPprots_8hr_Unch_Thresh_Num <- length(WCPprots_8hr_Unch_Thresh)
print("WCPprots_8hr_Unch_Thresh_Num")
print(WCPprots_8hr_Unch_Thresh_Num)
#8hr p-value
WCPprots_8hr_Inc_pval <- row.names(subset(WCP_0hr8hr_Data, (WCP_0hr8hr_Data$iBAQ_0hr8hr_Log2MedNorm_pval < wcppvalthreshold & WCP_0hr8hr_Data$iBAQ_0hr8hr_Log2MedNorm_Avg_Fc > 0)))
WCPprots_8hr_Inc_pval_Num <- length(WCPprots_8hr_Inc_pval)
print("WCPprots_8hr_Inc_pval_Num")
print(WCPprots_8hr_Inc_pval_Num)
WCPprots_8hr_Dec_pval <- row.names(subset(WCP_0hr8hr_Data, (WCP_0hr8hr_Data$iBAQ_0hr8hr_Log2MedNorm_pval < wcppvalthreshold & WCP_0hr8hr_Data$iBAQ_0hr8hr_Log2MedNorm_Avg_Fc < 0)))
WCPprots_8hr_Dec_pval_Num <- length(WCPprots_8hr_Dec_pval)
print("WCPprots_8hr_Dec_pval_Num")
print(WCPprots_8hr_Dec_pval_Num)
WCPprots_8hr_Unch_pval <- row.names(subset(WCP_0hr8hr_Data, (WCP_0hr8hr_Data$iBAQ_0hr8hr_Log2MedNorm_pval > wcppvalthreshold)))
WCPprots_8hr_Unch_pval_Num <- length(WCPprots_8hr_Unch_pval)
print("WCPprots_8hr_Unch_pval_Num")
print(WCPprots_8hr_Unch_pval_Num)
#8hr unique
WCPprots_8hr_Inc_unique <- row.names(subset(WCPdata,  (WCPdata$iBAQ_8hr_NumReps>=2 &  WCPdata$iBAQ_0hr_NumReps==0)))
WCPprots_8hr_Inc_unique_Num <- length(WCPprots_8hr_Inc_unique)
print("WCPprots_8hr_Inc_unique_Num")
print(WCPprots_8hr_Inc_unique_Num)
WCPprots_8hr_Dec_unique <- row.names(subset(WCPdata,  (WCPdata$iBAQ_8hr_NumReps==0 &  WCPdata$iBAQ_0hr_NumReps>=2)))
WCPprots_8hr_Dec_unique_Num <- length(WCPprots_8hr_Dec_unique)
print("WCPprots_8hr_Dec_unique_Num")
print(WCPprots_8hr_Dec_unique_Num)
#
#10hr mean
WCP_0hr10hr_Data <- subset(WCPdata, (WCPdata$iBAQ_0hr_NumReps>=2 & WCPdata$iBAQ_10hr_NumReps>=2))
WCP_0hr10hr_Fc_Avg <- mean(WCP_0hr10hr_Data$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc, na.rm=TRUE)
WCP_0hr10hr_Fc_Sd <- sd(WCP_0hr10hr_Data$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc, na.rm=TRUE)
print("WCP 10hr mean, sd")
print(WCP_0hr10hr_Fc_Avg)
print(WCP_0hr10hr_Fc_Sd)
WCPprots_10hr_Inc_MeanSd <- row.names(subset(WCP_0hr10hr_Data, WCP_0hr10hr_Data$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc>=WCP_0hr10hr_Fc_Avg+WCP_0hr10hr_Fc_Sd))
WCPprots_10hr_Inc_MeanSd_Num <- length(WCPprots_10hr_Inc_MeanSd)
print("WCPprots_10hr_Inc_MeanSd_Num")
print(WCPprots_10hr_Inc_MeanSd_Num)
WCPprots_10hr_Dec_MeanSd <- row.names(subset(WCP_0hr10hr_Data, WCP_0hr10hr_Data$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc<=WCP_0hr10hr_Fc_Avg-WCP_0hr10hr_Fc_Sd))
WCPprots_10hr_Dec_MeanSd_Num <- length(WCPprots_10hr_Dec_MeanSd)
print("WCPprots_10hr_Dec_MeanSd_Num")
print(WCPprots_10hr_Dec_MeanSd_Num)
WCPprots_10hr_Unch_MeanSd <- row.names(subset(WCP_0hr10hr_Data, (WCP_0hr10hr_Data$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc>WCP_0hr10hr_Fc_Avg-WCP_0hr10hr_Fc_Sd & WCP_0hr10hr_Data$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc<WCP_0hr10hr_Fc_Avg+WCP_0hr10hr_Fc_Sd)))
WCPprots_10hr_Unch_MeanSd_Num <- length(WCPprots_10hr_Unch_MeanSd)
print("WCPprots_10hr_Unch_MeanSd_Num")
print(WCPprots_10hr_Unch_MeanSd_Num)
#10hr thresh
WCPprots_10hr_Inc_Thresh <- row.names(subset(WCP_0hr10hr_Data, WCP_0hr10hr_Data$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc >= wcpfcthreshold))
WCPprots_10hr_Inc_Thresh_Num <- length(WCPprots_10hr_Inc_Thresh)
print("WCPprots_10hr_Inc_Thresh_Num")
print(WCPprots_10hr_Inc_Thresh_Num)
WCPprots_10hr_Dec_Thresh <- row.names(subset(WCP_0hr10hr_Data, WCP_0hr10hr_Data$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc <= (-1*wcpfcthreshold)))
WCPprots_10hr_Dec_Thresh_Num <- length(WCPprots_10hr_Dec_Thresh)
print("WCPprots_10hr_Dec_Thresh_Num")
print(WCPprots_10hr_Dec_Thresh_Num)
WCPprots_10hr_Unch_Thresh <- row.names(subset(WCP_0hr10hr_Data, (WCP_0hr10hr_Data$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc > (-1*wcpfcthreshold) & WCP_0hr10hr_Data$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc < wcpfcthreshold)))
WCPprots_10hr_Unch_Thresh_Num <- length(WCPprots_10hr_Unch_Thresh)
print("WCPprots_10hr_Unch_Thresh_Num")
print(WCPprots_10hr_Unch_Thresh_Num)
#10hr p-value
WCPprots_10hr_Inc_pval <- row.names(subset(WCP_0hr10hr_Data, (WCP_0hr10hr_Data$iBAQ_0hr10hr_Log2MedNorm_pval < wcppvalthreshold & WCP_0hr10hr_Data$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc > 0)))
WCPprots_10hr_Inc_pval_Num <- length(WCPprots_10hr_Inc_pval)
print("WCPprots_10hr_Inc_pval_Num")
print(WCPprots_10hr_Inc_pval_Num)
WCPprots_10hr_Dec_pval <- row.names(subset(WCP_0hr10hr_Data, (WCP_0hr10hr_Data$iBAQ_0hr10hr_Log2MedNorm_pval < wcppvalthreshold & WCP_0hr10hr_Data$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc < 0)))
WCPprots_10hr_Dec_pval_Num <- length(WCPprots_10hr_Dec_pval)
print("WCPprots_10hr_Dec_pval_Num")
print(WCPprots_10hr_Dec_pval_Num)
WCPprots_10hr_Unch_pval <- row.names(subset(WCP_0hr10hr_Data, (WCP_0hr10hr_Data$iBAQ_0hr10hr_Log2MedNorm_pval > wcppvalthreshold)))
WCPprots_10hr_Unch_pval_Num <- length(WCPprots_10hr_Unch_pval)
print("WCPprots_10hr_Unch_pval_Num")
print(WCPprots_10hr_Unch_pval_Num)
#10hr unique
WCPprots_10hr_Inc_unique <- row.names(subset(WCPdata,  (WCPdata$iBAQ_10hr_NumReps>=2 &  WCPdata$iBAQ_0hr_NumReps==0)))
WCPprots_10hr_Inc_unique_Num <- length(WCPprots_10hr_Inc_unique)
print("WCPprots_10hr_Inc_unique_Num")
print(WCPprots_10hr_Inc_unique_Num)
WCPprots_10hr_Dec_unique <- row.names(subset(WCPdata,  (WCPdata$iBAQ_10hr_NumReps==0 &  WCPdata$iBAQ_0hr_NumReps>=2)))
WCPprots_10hr_Dec_unique_Num <- length(WCPprots_10hr_Dec_unique)
print("WCPprots_10hr_Dec_unique_Num")
print(WCPprots_10hr_Dec_unique_Num)
#
#WCP decreased at 10hr based on pval<0.05 and fc<mean-sd or unique
WCPprots_10hr_dec_PvalMeanSd <- intersect(WCPprots_10hr_Dec_MeanSd, WCPprots_10hr_Dec_pval)
WCPprots_10hr_dec_PvalMeanSdUnique <- unique(union(WCPprots_10hr_dec_PvalMeanSd, WCPprots_10hr_Dec_unique))
WCPprots_10hr_dec_PvalMeanSdUnique_Num <- length(WCPprots_10hr_dec_PvalMeanSdUnique)
print("WCPprots_10hr_dec_PvalMeanSdUnique_Num")
print(WCPprots_10hr_dec_PvalMeanSdUnique_Num)
WCPprots_10hr_dec_PvalMeanSdUnique_Data <- WCPdata[row.names(WCPdata)%in%WCPprots_10hr_dec_PvalMeanSdUnique,]
write.table(WCPprots_10hr_dec_PvalMeanSdUnique_Data, file="WcpData_PvalMeanSdUnique_Decreased.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#WCP increased at 10hr based on pval<0.05 and fc>mean+sd or unique
WCPprots_10hr_inc_PvalMeanSd <- intersect(WCPprots_10hr_Inc_MeanSd, WCPprots_10hr_Inc_pval)
WCPprots_10hr_inc_PvalMeanSdUnique <- unique(union(WCPprots_10hr_inc_PvalMeanSd, WCPprots_10hr_Inc_unique))
WCPprots_10hr_inc_PvalMeanSdUnique_Num <- length(WCPprots_10hr_inc_PvalMeanSdUnique)
print("WCPprots_10hr_inc_PvalMeanSdUnique_Num")
print(WCPprots_10hr_inc_PvalMeanSdUnique_Num)
#WCP unchanged at 10hr based on pval>0.05 or pval<0.05 and not in threshold
WCPprots_10hr_unch_PvalMeanSd <- union(WCPprots_10hr_Unch_MeanSd, WCPprots_10hr_Unch_pval)
WCPprots_10hr_unch_PvalMeanSd_Num <- length(WCPprots_10hr_unch_PvalMeanSd)
print("WCPprots_10hr_unch_PvalMeanSd_Num")
print(WCPprots_10hr_unch_PvalMeanSd_Num)


##proteins that are increased or decreased in total protein KeGG
#based on the weighted average peptide KeGG fc and norm to WCP fc
normprotkeggthresh <- 1
KeGGProt_10hr_ProtFcNorm_Inc_thresh_Data <- subset(KeGGdataMerge_ProteinData, KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm >= normprotkeggthresh)
KeGGProt_10hr_ProtFcNorm_Inc_thresh <- KeGGProt_10hr_ProtFcNorm_Inc_thresh_Data$Protein
KeGGProt_10hr_ProtFcNorm_Inc_thresh_Num <- length(KeGGProt_10hr_ProtFcNorm_Inc_thresh)
print("KeGGProt_10hr_ProtFcNorm_Inc_thresh_Num")
print(KeGGProt_10hr_ProtFcNorm_Inc_thresh_Num)
KeGGProt_10hr_ProtFcNorm_Dec_thresh_Data <- subset(KeGGdataMerge_ProteinData, KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm <= -1*normprotkeggthresh)
KeGGProt_10hr_ProtFcNorm_Dec_thresh <- KeGGProt_10hr_ProtFcNorm_Dec_thresh_Data$Protein
KeGGProt_10hr_ProtFcNorm_Dec_thresh_Num <- length(KeGGProt_10hr_ProtFcNorm_Dec_thresh)
print("KeGGProt_10hr_ProtFcNorm_Dec_thresh_Num")
print(KeGGProt_10hr_ProtFcNorm_Dec_thresh_Num)
KeGGProt_10hr_ProtFcNorm_Unch_thresh_Data <- subset(KeGGdataMerge_ProteinData, (KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm < normprotkeggthresh & KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm > -1*normprotkeggthresh))
KeGGProt_10hr_ProtFcNorm_Unch_thresh <- KeGGProt_10hr_ProtFcNorm_Unch_thresh_Data$Protein
KeGGProt_10hr_ProtFcNorm_Unch_thresh_Num <- length(KeGGProt_10hr_ProtFcNorm_Unch_thresh)
print("KeGGProt_10hr_ProtFcNorm_Unch_thresh_Num")
print(KeGGProt_10hr_ProtFcNorm_Unch_thresh_Num)


##combine the list of proteins decreased based on WCP with the list of proteins that are increased in total protein Kegg
##the WCP decreased proteins are based on pvalue and mean,sd or unique
##the total protein KeGG is based on a threshold of normalized total protein KeGG fold change > 1
#
#decreased WCP
WCPprots_10hr_dec_MeanSdUnique <- union(WCPprots_10hr_Dec_MeanSd, WCPprots_10hr_Dec_unique)
WCPprots_10hr_dec_MeanSdUnique_Num <- length(WCPprots_10hr_dec_MeanSdUnique)
print("WCPprots_10hr_dec_MeanSdUnique_Num")
print(WCPprots_10hr_dec_MeanSdUnique_Num)
#increased WCP
WCPprots_10hr_inc_MeanSdUnique <- union(WCPprots_10hr_Inc_MeanSd, WCPprots_10hr_Inc_unique)
WCPprots_10hr_inc_MeanSdUnique_Num <- length(WCPprots_10hr_inc_MeanSdUnique)
print("WCPprots_10hr_inc_MeanSdUnique_Num")
print(WCPprots_10hr_inc_MeanSdUnique_Num)
#


########## THESE ARE THE PREDICTED SUBSTRATES ########## 
#
#
#
#
#
#increased KeGG and decreased WCP
KeGGProt_10hr_Inc_WcpProt_10hr_Dec <- intersect(WCPprots_10hr_dec_MeanSdUnique, KeGGProt_10hr_ProtFcNorm_Inc_thresh)
KeGGProt_10hr_Inc_WcpProt_10hr_Dec_Num <- length(KeGGProt_10hr_Inc_WcpProt_10hr_Dec)
print("KeGGProt_10hr_Inc_WcpProt_10hr_Dec_Num")
print(KeGGProt_10hr_Inc_WcpProt_10hr_Dec_Num)
print(KeGGProt_10hr_Inc_WcpProt_10hr_Dec)
KeGGProt_10hr_Inc_WcpProt_10hr_Dec_ProtData <- KeGGdataMerge_ProteinData[row.names(KeGGdataMerge_ProteinData)%in%KeGGProt_10hr_Inc_WcpProt_10hr_Dec,]
dim(KeGGProt_10hr_Inc_WcpProt_10hr_Dec_ProtData)
write.table(KeGGProt_10hr_Inc_WcpProt_10hr_Dec_ProtData, file="KeGGProteinData_KeGGProtIncWcpDec.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
KeGGProt_10hr_Inc_WcpProt_10hr_Dec_PepData <- KeGGdataMerge[KeGGdataMerge$Protein%in%KeGGProt_10hr_Inc_WcpProt_10hr_Dec,]
dim(KeGGProt_10hr_Inc_WcpProt_10hr_Dec_PepData)
write.table(KeGGProt_10hr_Inc_WcpProt_10hr_Dec_PepData, file="KeGGPeptideData_KeGGProtIncWcpDec.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
#inceased KeGG and unchanged WCP
KeGGProt_10hr_Inc_WcpProt_10hr_Unch <- intersect(WCPprots_10hr_Unch_MeanSd, KeGGProt_10hr_ProtFcNorm_Inc_thresh)
KeGGProt_10hr_Inc_WcpProt_10hr_Unch_Num <- length(KeGGProt_10hr_Inc_WcpProt_10hr_Unch)
print("KeGGProt_10hr_Inc_WcpProt_10hr_Unch_Num")
print(KeGGProt_10hr_Inc_WcpProt_10hr_Unch_Num)
print(KeGGProt_10hr_Inc_WcpProt_10hr_Unch)
KeGGProt_10hr_Inc_WcpProt_10hr_Unch_ProtData <- KeGGdataMerge_ProteinData[row.names(KeGGdataMerge_ProteinData)%in%KeGGProt_10hr_Inc_WcpProt_10hr_Unch,]
dim(KeGGProt_10hr_Inc_WcpProt_10hr_Unch_ProtData)
write.table(KeGGProt_10hr_Inc_WcpProt_10hr_Unch_ProtData, file="KeGGProteinData_KeGGProtIncWcpUnch.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
KeGGProt_10hr_Inc_WcpProt_10hr_Unch_PepData <- KeGGdataMerge[KeGGdataMerge$Protein%in%KeGGProt_10hr_Inc_WcpProt_10hr_Unch,]
dim(KeGGProt_10hr_Inc_WcpProt_10hr_Unch_PepData)
write.table(KeGGProt_10hr_Inc_WcpProt_10hr_Unch_PepData, file="KeGGPeptideData_KeGGProtIncWcpUnch.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
#
#increased KeGG and increased WCP
KeGGProt_10hr_Inc_WcpProt_10hr_Inc <- intersect(WCPprots_10hr_inc_MeanSdUnique, KeGGProt_10hr_ProtFcNorm_Inc_thresh)
KeGGProt_10hr_Inc_WcpProt_10hr_Inc_Num <- length(KeGGProt_10hr_Inc_WcpProt_10hr_Inc)
print("KeGGProt_10hr_Inc_WcpProt_10hr_Inc_Num")
print(KeGGProt_10hr_Inc_WcpProt_10hr_Inc_Num)
print(KeGGProt_10hr_Inc_WcpProt_10hr_Inc)
KeGGProt_10hr_Inc_WcpProt_10hr_Inc_ProtData <- KeGGdataMerge_ProteinData[row.names(KeGGdataMerge_ProteinData)%in%KeGGProt_10hr_Inc_WcpProt_10hr_Inc,]
dim(KeGGProt_10hr_Inc_WcpProt_10hr_Inc_ProtData)
write.table(KeGGProt_10hr_Inc_WcpProt_10hr_Inc_ProtData, file="KeGGProteinData_KeGGProtIncWcpInc.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
KeGGProt_10hr_Inc_WcpProt_10hr_Inc_PepData <- KeGGdataMerge[KeGGdataMerge$Protein%in%KeGGProt_10hr_Inc_WcpProt_10hr_Inc,]
dim(KeGGProt_10hr_Inc_WcpProt_10hr_Inc_PepData)
write.table(KeGGProt_10hr_Inc_WcpProt_10hr_Inc_PepData, file="KeGGPeptideData_KeGGProtIncWcpInc.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
#
#
#
#
#
########## THESE ARE THE PREDICTED SUBSTRATES ########## 


#decreased KeGG and increased WCP
KeGGProt_10hr_Dec_WcpProt_10hr_Inc <- intersect(WCPprots_10hr_inc_MeanSdUnique, KeGGProt_10hr_ProtFcNorm_Dec_thresh)
KeGGProt_10hr_Dec_WcpProt_10hr_Inc_Num <- length(KeGGProt_10hr_Dec_WcpProt_10hr_Inc)
print("KeGGProt_10hr_Dec_WcpProt_10hr_Inc_Num")
print(KeGGProt_10hr_Dec_WcpProt_10hr_Inc_Num)
print(KeGGProt_10hr_Dec_WcpProt_10hr_Inc)
KeGGProt_10hr_Dec_WcpProt_10hr_Inc_Data <- KeGGdataMerge_ProteinData[row.names(KeGGdataMerge_ProteinData)%in%KeGGProt_10hr_Dec_WcpProt_10hr_Inc,]
dim(KeGGProt_10hr_Dec_WcpProt_10hr_Inc_Data)
write.table(KeGGProt_10hr_Dec_WcpProt_10hr_Inc_Data, file="KeGGProteinData_KeGGProtDecWcpInc.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
#decreased KeGG and unchanged WCP
KeGGProt_10hr_Dec_WcpProt_10hr_Unch <- intersect(WCPprots_10hr_Unch_MeanSd, KeGGProt_10hr_ProtFcNorm_Dec_thresh)
KeGGProt_10hr_Dec_WcpProt_10hr_Unch_Num <- length(KeGGProt_10hr_Dec_WcpProt_10hr_Unch)
print("KeGGProt_10hr_Dec_WcpProt_10hr_Unch_Num")
print(KeGGProt_10hr_Dec_WcpProt_10hr_Unch_Num)
print(KeGGProt_10hr_Dec_WcpProt_10hr_Unch)
KeGGProt_10hr_Dec_WcpProt_10hr_Unch_Data <- KeGGdataMerge_ProteinData[row.names(KeGGdataMerge_ProteinData)%in%KeGGProt_10hr_Dec_WcpProt_10hr_Unch,]
dim(KeGGProt_10hr_Dec_WcpProt_10hr_Unch_Data)
write.table(KeGGProt_10hr_Dec_WcpProt_10hr_Unch_Data, file="KeGGProteinData_KeGGProtDecWcpUnch.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
#decreased KeGG and decreased WCP
KeGG_10hr_Dec_WcpProt_10hr_Dec <- intersect(WCPprots_10hr_dec_MeanSdUnique, KeGGProt_10hr_ProtFcNorm_Dec_thresh)
KeGG_10hr_Dec_WcpProt_10hr_Dec_Num <- length(KeGG_10hr_Dec_WcpProt_10hr_Dec)
print("KeGG_10hr_Dec_WcpProt_10hr_Dec_Num")
print(KeGG_10hr_Dec_WcpProt_10hr_Dec_Num)
print(KeGG_10hr_Dec_WcpProt_10hr_Dec)
KeGG_10hr_Dec_WcpProt_10hr_Dec_Data <- KeGGdataMerge_ProteinData[row.names(KeGGdataMerge_ProteinData)%in%KeGG_10hr_Dec_WcpProt_10hr_Dec,]
dim(KeGG_10hr_Dec_WcpProt_10hr_Dec_Data)
write.table(KeGG_10hr_Dec_WcpProt_10hr_Dec_Data, file="KeGGProteinData_KeGGProtDecWcpDec.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#

##peptides that are increased or decreased in KeGG abundance based on normalized KeGG fold change mean, p-value or unique
##limit analysis to comparisons in which the peptide was identified in at least 2 replicates for 0 and time point
##the threshold for increase/decrease/unchange is based on mean fold change +/- 1 standard deviation
keggfcthreshold=1
keggpvalthreshold<-0.1
#6hr mean
KeGGNormFc_0hr6hr_Data <- subset(KeGGdataMerge, (KeGGdataMerge$KeGGInt.0hr_NumReps>=2 | KeGGdataMerge$KeGGInt.6hr_NumReps>=2))
KeGGNormFc_0hr6hr_Avg <- mean(KeGGNormFc_0hr6hr_Data$KeGG_Int_0hr6hr_FcNorm_Avg, na.rm=TRUE)
KeGGNormFc_0hr6hr_Sd <- sd(KeGGNormFc_0hr6hr_Data$KeGG_Int_0hr6hr_FcNorm_Avg, na.rm=TRUE)
KeGG_NormFc_6hr_Inc_MeanSd <- row.names(subset(KeGGNormFc_0hr6hr_Data, KeGGNormFc_0hr6hr_Data$KeGG_Int_0hr6hr_FcNorm_Avg>KeGGNormFc_0hr6hr_Avg+KeGGNormFc_0hr6hr_Sd))
KeGG_NormFc_6hr_Inc_MeanSd_Num <- length(KeGG_NormFc_6hr_Inc_MeanSd)
print("KeGG_NormFc_6hr_Inc_MeanSd_Num")
print(KeGG_NormFc_6hr_Inc_MeanSd_Num)
KeGG_NormFc_6hr_Dec_MeanSd <- row.names(subset(KeGGNormFc_0hr6hr_Data, KeGGNormFc_0hr6hr_Data$KeGG_Int_0hr6hr_FcNorm_Avg<KeGGNormFc_0hr6hr_Avg-KeGGNormFc_0hr6hr_Sd))
KeGG_NormFc_6hr_Dec_MeanSd_Num <- length(KeGG_NormFc_6hr_Dec_MeanSd)
print("KeGG_NormFc_6hr_Dec_MeanSd_Num")
print(KeGG_NormFc_6hr_Dec_MeanSd_Num)
KeGG_NormFc_6hr_Unch_MeanSd <- row.names(subset(KeGGNormFc_0hr6hr_Data, (KeGGNormFc_0hr6hr_Data$KeGG_Int_0hr6hr_FcNorm_Avg<KeGGNormFc_0hr6hr_Avg+KeGGNormFc_0hr6hr_Sd & KeGGNormFc_0hr6hr_Data$KeGG_Int_0hr6hr_FcNorm_Avg>KeGGNormFc_0hr6hr_Avg-KeGGNormFc_0hr6hr_Sd)))
KeGG_NormFc_6hr_Unch_MeanSd_Num <- length(KeGG_NormFc_6hr_Unch_MeanSd)
print("KeGG_NormFc_6hr_Unch_MeanSd_Num")
print(KeGG_NormFc_6hr_Unch_MeanSd_Num)
#6hr p-val
KeGG_NormFc_6hr_Inc_pval <- row.names(subset(KeGGNormFc_0hr6hr_Data, (KeGGNormFc_0hr6hr_Data$KeGG_Int_0hr6hr_FcNorm_onesidePval < 0.1 & KeGGNormFc_0hr6hr_Data$KeGG_Int_0hr6hr_FcNorm_Avg > 0)))
KeGG_NormFc_6hr_Inc_pval_Num <- length(KeGG_NormFc_6hr_Inc_pval)
print("KeGG_NormFc_6hr_Inc_pval_Num")
print(KeGG_NormFc_6hr_Inc_pval_Num)
KeGG_NormFc_6hr_Dec_pval <- row.names(subset(KeGGNormFc_0hr6hr_Data, (KeGGNormFc_0hr6hr_Data$KeGG_Int_0hr6hr_FcNorm_onesidePval < 0.1 & KeGGNormFc_0hr6hr_Data$KeGG_Int_0hr6hr_FcNorm_Avg < 0)))
KeGG_NormFc_6hr_Dec_pval_Num <- length(KeGG_NormFc_6hr_Dec_pval)
print("KeGG_NormFc_6hr_Dec_pval_Num")
print(KeGG_NormFc_6hr_Dec_pval_Num)
KeGG_NormFc_6hr_Unch_pval <- row.names(subset(KeGGNormFc_0hr6hr_Data, (KeGGNormFc_0hr6hr_Data$KeGG_Int_0hr6hr_FcNorm_onesidePval > 0.1)))
KeGG_NormFc_6hr_Unch_pval_Num <- length(KeGG_NormFc_6hr_Unch_pval)
print("KeGG_NormFc_6hr_Unch_pval_Num")
print(KeGG_NormFc_6hr_Unch_pval_Num)
#6hr unique
KeGG_NormFc_6hr_Inc_Unique <- row.names(subset(KeGGdataMerge, (KeGGdataMerge$KeGGInt.0hr_NumReps==0 & KeGGdataMerge$KeGGInt.6hr_NumReps>=2)))
KeGG_NormFc_6hr_Inc_Unique_Num <- length(KeGG_NormFc_6hr_Inc_Unique)
print("KeGG_NormFc_6hr_Inc_Unique_Num")
print(KeGG_NormFc_6hr_Inc_Unique_Num)
KeGG_NormFc_6hr_Dec_Unique <- row.names(subset(KeGGdataMerge, (KeGGdataMerge$KeGGInt.0hr_NumReps>=2 & KeGGdataMerge$KeGGInt.6hr_NumReps==0)))
KeGG_NormFc_6hr_Dec_Unique_Num <- length(KeGG_NormFc_6hr_Dec_Unique)
print("KeGG_NormFc_6hr_Dec_Unique_Num")
print(KeGG_NormFc_6hr_Dec_Unique_Num)
#8hr mean
KeGGNormFc_0hr8hr_Data <- subset(KeGGdataMerge, (KeGGdataMerge$KeGGInt.0hr_NumReps>=2 | KeGGdataMerge$KeGGInt.8hr_NumReps>=2))
KeGGNormFc_0hr8hr_Avg <- mean(KeGGNormFc_0hr8hr_Data$KeGG_Int_0hr8hr_FcNorm_Avg, na.rm=TRUE)
KeGGNormFc_0hr8hr_Sd <- sd(KeGGNormFc_0hr8hr_Data$KeGG_Int_0hr8hr_FcNorm_Avg, na.rm=TRUE)
KeGG_NormFc_8hr_Inc_MeanSd <- row.names(subset(KeGGNormFc_0hr8hr_Data, KeGGNormFc_0hr8hr_Data$KeGG_Int_0hr8hr_FcNorm_Avg>KeGGNormFc_0hr8hr_Avg+KeGGNormFc_0hr8hr_Sd))
KeGG_NormFc_8hr_Inc_MeanSd_Num <- length(KeGG_NormFc_8hr_Inc_MeanSd)
print("KeGG_NormFc_8hr_Inc_MeanSd_Num")
print(KeGG_NormFc_8hr_Inc_MeanSd_Num)
KeGG_NormFc_8hr_Dec_MeanSd <- row.names(subset(KeGGNormFc_0hr8hr_Data, KeGGNormFc_0hr8hr_Data$KeGG_Int_0hr8hr_FcNorm_Avg<KeGGNormFc_0hr8hr_Avg-KeGGNormFc_0hr8hr_Sd))
KeGG_NormFc_8hr_Dec_MeanSd_Num <- length(KeGG_NormFc_8hr_Dec_MeanSd)
print("KeGG_NormFc_8hr_Dec_MeanSd_Num")
print(KeGG_NormFc_8hr_Dec_MeanSd_Num)
KeGG_NormFc_8hr_Unch_MeanSd <- row.names(subset(KeGGNormFc_0hr8hr_Data, (KeGGNormFc_0hr8hr_Data$KeGG_Int_0hr8hr_FcNorm_Avg<KeGGNormFc_0hr8hr_Avg+KeGGNormFc_0hr8hr_Sd & KeGGNormFc_0hr8hr_Data$KeGG_Int_0hr8hr_FcNorm_Avg>KeGGNormFc_0hr8hr_Avg-KeGGNormFc_0hr8hr_Sd)))
KeGG_NormFc_8hr_Unch_MeanSd_Num <- length(KeGG_NormFc_8hr_Unch_MeanSd)
print("KeGG_NormFc_8hr_Unch_MeanSd_Num")
print(KeGG_NormFc_8hr_Unch_MeanSd_Num)
#8hr p-val
KeGG_NormFc_8hr_Inc_pval <- row.names(subset(KeGGNormFc_0hr8hr_Data, (KeGGNormFc_0hr8hr_Data$KeGG_Int_0hr8hr_FcNorm_onesidePval < 0.1 & KeGGNormFc_0hr8hr_Data$KeGG_Int_0hr8hr_FcNorm_Avg > 0)))
KeGG_NormFc_8hr_Inc_pval_Num <- length(KeGG_NormFc_8hr_Inc_pval)
print("KeGG_NormFc_8hr_Inc_pval_Num")
print(KeGG_NormFc_8hr_Inc_pval_Num)
KeGG_NormFc_8hr_Dec_pval <- row.names(subset(KeGGNormFc_0hr8hr_Data, (KeGGNormFc_0hr8hr_Data$KeGG_Int_0hr8hr_FcNorm_onesidePval < 0.1 & KeGGNormFc_0hr8hr_Data$KeGG_Int_0hr8hr_FcNorm_Avg < 0)))
KeGG_NormFc_8hr_Dec_pval_Num <- length(KeGG_NormFc_8hr_Dec_pval)
print("KeGG_NormFc_8hr_Dec_pval_Num")
print(KeGG_NormFc_8hr_Dec_pval_Num)
KeGG_NormFc_8hr_Unch_pval <- row.names(subset(KeGGNormFc_0hr8hr_Data, (KeGGNormFc_0hr8hr_Data$KeGG_Int_0hr8hr_FcNorm_onesidePval > 0.1)))
KeGG_NormFc_8hr_Unch_pval_Num <- length(KeGG_NormFc_8hr_Unch_pval)
print("KeGG_NormFc_8hr_Unch_pval_Num")
print(KeGG_NormFc_8hr_Unch_pval_Num)
#8hr unique
KeGG_NormFc_8hr_Inc_Unique <- row.names(subset(KeGGdataMerge, (KeGGdataMerge$KeGGInt.0hr_NumReps==0 & KeGGdataMerge$KeGGInt.8hr_NumReps>=2)))
KeGG_NormFc_8hr_Inc_Unique_Num <- length(KeGG_NormFc_8hr_Inc_Unique)
print("KeGG_NormFc_8hr_Inc_Unique_Num")
print(KeGG_NormFc_8hr_Inc_Unique_Num)
KeGG_NormFc_8hr_Dec_Unique <- row.names(subset(KeGGdataMerge, (KeGGdataMerge$KeGGInt.0hr_NumReps>=2 & KeGGdataMerge$KeGGInt.8hr_NumReps==0)))
KeGG_NormFc_8hr_Dec_Unique_Num <- length(KeGG_NormFc_8hr_Dec_Unique)
print("KeGG_NormFc_8hr_Dec_Unique_Num")
print(KeGG_NormFc_8hr_Dec_Unique_Num)
#10hr mean
KeGGNormFc_0hr10hr_Data <- subset(KeGGdataMerge, (KeGGdataMerge$KeGGInt.0hr_NumReps>=2 | KeGGdataMerge$KeGGInt.10hr_NumReps>=2))
KeGGNormFc_0hr10hr_Avg <- mean(KeGGNormFc_0hr10hr_Data$KeGG_Int_0hr10hr_FcNorm_Avg, na.rm=TRUE)
KeGGNormFc_0hr10hr_Sd <- sd(KeGGNormFc_0hr10hr_Data$KeGG_Int_0hr10hr_FcNorm_Avg, na.rm=TRUE)
KeGG_NormFc_10hr_Inc_MeanSd <- row.names(subset(KeGGNormFc_0hr10hr_Data, KeGGNormFc_0hr10hr_Data$KeGG_Int_0hr10hr_FcNorm_Avg>KeGGNormFc_0hr10hr_Avg+KeGGNormFc_0hr10hr_Sd))
KeGG_NormFc_10hr_Inc_MeanSd_Num <- length(KeGG_NormFc_10hr_Inc_MeanSd)
print("KeGG_NormFc_10hr_Inc_MeanSd_Num")
print(KeGG_NormFc_10hr_Inc_MeanSd_Num)
KeGG_NormFc_10hr_Dec_MeanSd <- row.names(subset(KeGGNormFc_0hr10hr_Data, KeGGNormFc_0hr10hr_Data$KeGG_Int_0hr10hr_FcNorm_Avg<KeGGNormFc_0hr10hr_Avg-KeGGNormFc_0hr10hr_Sd))
KeGG_NormFc_10hr_Dec_MeanSd_Num <- length(KeGG_NormFc_10hr_Dec_MeanSd)
print("KeGG_NormFc_10hr_Dec_MeanSd_Num")
print(KeGG_NormFc_10hr_Dec_MeanSd_Num)
KeGG_NormFc_10hr_Unch_MeanSd <- row.names(subset(KeGGNormFc_0hr10hr_Data, (KeGGNormFc_0hr10hr_Data$KeGG_Int_0hr10hr_FcNorm_Avg<KeGGNormFc_0hr10hr_Avg+KeGGNormFc_0hr10hr_Sd & KeGGNormFc_0hr10hr_Data$KeGG_Int_0hr10hr_FcNorm_Avg>KeGGNormFc_0hr10hr_Avg-KeGGNormFc_0hr10hr_Sd)))
KeGG_NormFc_10hr_Unch_MeanSd_Num <- length(KeGG_NormFc_10hr_Unch_MeanSd)
print("KeGG_NormFc_10hr_Unch_MeanSd_Num")
print(KeGG_NormFc_10hr_Unch_MeanSd_Num)
#10hr p-val
KeGG_NormFc_10hr_Inc_pval <- row.names(subset(KeGGNormFc_0hr10hr_Data, (KeGGNormFc_0hr10hr_Data$KeGG_Int_0hr10hr_FcNorm_onesidePval < keggpvalthreshold & KeGGNormFc_0hr10hr_Data$KeGG_Int_0hr10hr_FcNorm_Avg > keggfcthreshold)))
KeGG_NormFc_10hr_Inc_pval_Num <- length(KeGG_NormFc_10hr_Inc_pval)
print("KeGG_NormFc_10hr_Inc_pval_Num")
print(KeGG_NormFc_10hr_Inc_pval_Num)
KeGG_NormFc_10hr_Dec_pval <- row.names(subset(KeGGNormFc_0hr10hr_Data, (KeGGNormFc_0hr10hr_Data$KeGG_Int_0hr10hr_FcNorm_onesidePval < keggpvalthreshold & KeGGNormFc_0hr10hr_Data$KeGG_Int_0hr10hr_FcNorm_Avg < keggfcthreshold)))
KeGG_NormFc_10hr_Dec_pval_Num <- length(KeGG_NormFc_10hr_Dec_pval)
print("KeGG_NormFc_10hr_Dec_pval_Num")
print(KeGG_NormFc_10hr_Dec_pval_Num)
KeGG_NormFc_10hr_Unch_pval <- row.names(subset(KeGGNormFc_0hr10hr_Data, (KeGGNormFc_0hr10hr_Data$KeGG_Int_0hr10hr_FcNorm_onesidePval > keggpvalthreshold)))
KeGG_NormFc_10hr_Unch_pval_Num <- length(KeGG_NormFc_10hr_Unch_pval)
print("KeGG_NormFc_10hr_Unch_pval_Num")
print(KeGG_NormFc_10hr_Unch_pval_Num)
#10hr unique
KeGG_NormFc_10hr_Inc_Unique <- row.names(subset(KeGGdataMerge, (KeGGdataMerge$KeGGInt.0hr_NumReps==0 & KeGGdataMerge$KeGGInt.10hr_NumReps>=2)))
KeGG_NormFc_10hr_Inc_Unique_Num <- length(KeGG_NormFc_10hr_Inc_Unique)
print("KeGG_NormFc_10hr_Inc_Unique_Num")
print(KeGG_NormFc_10hr_Inc_Unique_Num)
KeGG_NormFc_10hr_Dec_Unique <- row.names(subset(KeGGdataMerge, (KeGGdataMerge$KeGGInt.0hr_NumReps>=2 & KeGGdataMerge$KeGGInt.10hr_NumReps==0)))
KeGG_NormFc_10hr_Dec_Unique_Num <- length(KeGG_NormFc_10hr_Dec_Unique)
print("KeGG_NormFc_10hr_Dec_Unique_Num")
print(KeGG_NormFc_10hr_Dec_Unique_Num)

stop("NNNN")

##select increased peptides from proteins that are in the "ubiquitinated and decreased" or "ubiquitinated and unchanged" categories
##the categories are based on the total protein-based KeGG fold change analysis normalized to WCP
##these categories are based on the integration of the protein-based KeGG and WCP data
#combine increased by p-value and unique peptides
KeGG_NormFc_10hr_Inc_PvalUnique <- c(KeGG_NormFc_10hr_Inc_pval, KeGG_NormFc_10hr_Inc_Unique)
KeGG_NormFc_10hr_Inc_PvalUnique_Num <- length(KeGG_NormFc_10hr_Inc_PvalUnique)
print("KeGG_NormFc_10hr_Inc_PvalUnique_Num")
print(KeGG_NormFc_10hr_Inc_PvalUnique_Num)
#filter KeGG data for increased peptides
KeGG_NormFc_10hr_Inc_PvalUnique_Data <- KeGGdataMerge[row.names(KeGGdataMerge)%in%KeGG_NormFc_10hr_Inc_PvalUnique,]
dim(KeGG_NormFc_10hr_Inc_PvalUnique_Data)
#select peptides from proteins in protein-based KeGG increase and WCP decrease set
KeGGPeps_KeGGProt_10hr_Inc_WcpProt_10hr_Dec_Data <- KeGG_NormFc_10hr_Inc_PvalUnique_Data[KeGG_NormFc_10hr_Inc_PvalUnique_Data$Protein%in%KeGGProt_10hr_Inc_WcpProt_10hr_Dec,]
dim(KeGGPeps_KeGGProt_10hr_Inc_WcpProt_10hr_Dec_Data)
write.table(KeGGPeps_KeGGProt_10hr_Inc_WcpProt_10hr_Dec_Data, file="KeGGPeptideData_KeGGProtIncWcpDec.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#select peptides from proteins in protein-based KeGG increase and WCP unchanged set
KeGGPeps_KeGGProt_10hr_Inc_WcpProt_10hr_Unch_Data <- KeGG_NormFc_10hr_Inc_PvalUnique_Data[KeGG_NormFc_10hr_Inc_PvalUnique_Data$Protein%in%KeGGProt_10hr_Inc_WcpProt_10hr_Unch,]
dim(KeGGPeps_KeGGProt_10hr_Inc_WcpProt_10hr_Unch_Data)
write.table(KeGGPeps_KeGGProt_10hr_Inc_WcpProt_10hr_Unch_Data, file="KeGGPeptideData_KeGGProtIncWcpUnch.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)


if ( FALSE )
{
	 ##select increased KeGG peptides by combining peptides found uniquely for a time point with peptides increased at that timepoint
	 #6hr KeGG increase
	 KeGG_6hr_Increased <- Reduce(union, list(KeGGUnique_0hr6hr, KeGGNormFcWcpFc_6hr_Inc_MeanSd))
	 KeGG_6hr_Increased_Num <- length(KeGG_6hr_Increased)
	 print("KeGG_6hr_Increased_Num")
	 print(KeGG_6hr_Increased_Num)
	 #8hr KeGG increase
	 KeGG_8hr_Increased <- Reduce(union, list(KeGGUnique_0hr8hr, KeGGNormFcWcpFc_8hr_Inc_MeanSd))
	 KeGG_8hr_Increased_Num <- length(KeGG_8hr_Increased)
	 print("KeGG_8hr_Increased_Num")
	 print(KeGG_8hr_Increased_Num)
	 #10hr KeGG increase
	 KeGG_10hr_Increased <- Reduce(union, list(KeGGUnique_0hr10hr, KeGGNormFcWcpFc_10hr_Inc_MeanSd))
	 KeGG_10hr_Increased_Num <- length(KeGG_10hr_Increased)
	 print("KeGG_10hr_Increased_Num")
	 print(KeGG_10hr_Increased_Num)
	 #6hr-8hr-10hr increase
	 KeGG_6hr8hr10hr_Inc_Union <- Reduce(union, list(KeGG_6hr_Increased, KeGG_8hr_Increased, KeGG_10hr_Increased))
	 KeGG_6hr8hr10hr_Inc_Union_Num <- length(KeGG_6hr8hr10hr_Inc_Union)
	 print("KeGG_6hr8hr10hr_Inc_Union_Num")
	 print(KeGG_6hr8hr10hr_Inc_Union_Num)
	 KeGG_6hr8hr10hr_Inc_Union_Data <- KeGGdataMerge[row.names(KeGGdataMerge)%in%KeGG_6hr8hr10hr_Inc_Union,]
	 dim(KeGG_6hr8hr10hr_Inc_Union_Data)
	 write.table(KeGG_6hr8hr10hr_Inc_Union_Data, file="KeGGPepData_AllKeGGUniqueAndInc.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
	 #
	 
	 ##combine decreased wcp proteins with increased KeGG peptides
	 #this analysis is stringent. normalized KeGG abundance must be
	 #apparent before or during measured wcp decrease
	 #6-8-10hr
	 tmpdata_6hrKeGG <- KeGGdataMerge[row.names(KeGGdataMerge)%in%KeGG_6hr_IncOrUnch,]
	 ##########tmpdata_6hrKeGG <- KeGGdataMerge[row.names(KeGGdataMerge)%in%KeGG_6hr_Increased,]
	 KeGG_6hr_PredDeg_Data <- tmpdata_6hrKeGG[tmpdata_6hrKeGG$WCPcorrProt%in%WCP_6hr8hr10hr_Decreased,]
	 write.table(KeGG_6hr_PredDeg_Data, file="KeGG_6hr_PredDeg_Data.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
	 KeGG_6hr_PredDeg_Peps <- row.names(KeGG_6hr_PredDeg_Data)
	 KeGG_6hr_PredDeg_Peps_Num <- length(KeGG_6hr_PredDeg_Peps)
	 print("KeGG_6hr_PredDeg_Peps_Num")
	 print(KeGG_6hr_PredDeg_Peps_Num)
	 #8-10hr
	 KeGG_6hr8hr_IncOrUnch <- union(KeGG_6hr_IncOrUnch, KeGG_8hr_IncOrUnch)
	 ##########KeGG_6hr8hr_IncOrUnch <- intersect(KeGG_6hr_IncOrUnch, KeGG_8hr_IncOrUnch)
	 ##########KeGG_6hr8hr_IncOrUnch <- union(KeGG_6hr_Increased, KeGG_8hr_Increased)
	 tmpdata_6hr8hrKeGG <- KeGGdataMerge[row.names(KeGGdataMerge)%in%KeGG_6hr8hr_IncOrUnch,]
	 KeGG_6hr8hr_PredDeg_Data <- tmpdata_6hr8hrKeGG[tmpdata_6hr8hrKeGG$WCPcorrProt%in%WCP_8hr10hr_Decreased,]
	 write.table(KeGG_6hr8hr_PredDeg_Data, file="KeGG_6hr8hr_PredDeg_Data.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
	 KeGG_6hr8hr_PredDeg_Peps <- row.names(KeGG_6hr8hr_PredDeg_Data)
	 KeGG_6hr8hr_PredDeg_Peps_Num <- length(KeGG_6hr8hr_PredDeg_Peps)
	 print("KeGG_6hr8hr_PredDeg_Peps_Num")
	 print(KeGG_6hr8hr_PredDeg_Peps_Num)
	 #10hr
	 KeGG_6hr8hr10hr_IncOrUnch_Union <- union(KeGG_6hr8hr_IncOrUnch, KeGG_10hr_IncOrUnch)
	 ##########KeGG_6hr8hr10hr_IncOrUnch_Union <- intersect(KeGG_6hr8hr_IncOrUnch, KeGG_10hr_IncOrUnch)
	 ##########KeGG_6hr8hr10hr_IncOrUnch_Union <- union(KeGG_6hr8hr_IncOrUnch, KeGG_10hr_Increased)
	 tmpdata_6hr8hr10hrKeGG <- KeGGdataMerge[row.names(KeGGdataMerge)%in%KeGG_6hr8hr10hr_IncOrUnch_Union,]
	 KeGG_6hr8hr10hr_PredDeg_Data <- tmpdata_6hr8hr10hrKeGG[tmpdata_6hr8hr10hrKeGG$WCPcorrProt%in%WCP_10hr_Decreased,]
	 write.table(KeGG_6hr8hr10hr_PredDeg_Data, file="KeGG_6hr8hr10hr_PredDeg_Data.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
	 KeGG_6hr8hr10hr_PredDeg_Peps <- row.names(KeGG_6hr8hr10hr_PredDeg_Data)
	 KeGG_6hr8hr10hr_PredDeg_Peps_Num <- length(KeGG_6hr8hr10hr_PredDeg_Peps)
	 print("KeGG_6hr8hr10hr_PredDeg_Peps_Num")
	 print(KeGG_6hr8hr10hr_PredDeg_Peps_Num)
	 #all predicted deg
	 KeGG_AllPredDeg_Peps <- Reduce(union, list(KeGG_6hr_PredDeg_Peps, KeGG_6hr8hr_PredDeg_Peps, KeGG_6hr8hr10hr_PredDeg_Peps))
	 KeGG_AllPredDeg_Peps_Num <- length(KeGG_AllPredDeg_Peps)
	 print("KeGG_AllPredDeg_Peps_Num")
	 print(KeGG_AllPredDeg_Peps_Num)
	 KeGG_AllPredDeg_Peps_Data <- KeGGdataMerge[row.names(KeGGdataMerge)%in%KeGG_AllPredDeg_Peps,]
	 write.table(KeGG_AllPredDeg_Peps_Data, file="KeGGPepData_AllDecWithKeGGInc_ConsistentTiming.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
	 
	 
	 ##combine decreased wcp proteins with increased KeGG peptides
	 KeGG_AllPredUnch_Data <- KeGG_6hr8hr10hr_Inc_Union_Data[KeGG_6hr8hr10hr_Inc_Union_Data$WCPcorrProt%in%WCP_AllUnch,]
	 print("KeGG_AllPredUnch_Data dim")
	 dim(KeGG_AllPredUnch_Data)
	 KeGG_AllPredUnch_Peps <- row.names(KeGG_AllPredUnch_Data)
	 KeGG_AllPredUnch_Peps_Num <- length(KeGG_AllPredUnch_Peps)
	 print("KeGG_AllPredUnch_Peps_Num")
	 print(KeGG_AllPredUnch_Peps_Num)
	 write.table(KeGG_AllPredUnch_Data, file="KeGGPepData_AllUnchWithKeGGInc_ConsistentTiming.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
	 
	 
	 #proteins that are decreased in the whole cell proteome and have KeGG at 6, 8 or 10 hr
	 #proteins that are quantified uniquely at 0hr
	 WCPprots_0hr_Unique <- row.names(subset(WCPdata, ((WCPdata$iBAQ_0hr_NumReps>=2) & (WCPdata$iBAQ_10hr_NumReps==0))))
	 WCPprots_0hr_Unique_Num <- length(WCPprots_0hr_Unique)
	 print("WCPprots_0hr_Unique_Num")
	 print(WCPprots_0hr_Unique_Num)
	 #proteins that are significantly decreased at 10hr
	 WCPprots_10hr_SigDec <- row.names(subset(WCPdata, ((WCPdata$iBAQ_0hr10hr_Log2MedNorm_pval<=0.05) & (WCPdata$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc < 0) & (WCPdata$iBAQ_0hr_NumReps>=2) & (WCPdata$iBAQ_10hr_NumReps>=2))))
	 WCPprots_10hr_SigDec_Num <- length(WCPprots_10hr_SigDec)
	 print("WCPprots_10hr_SigDec_Num")
	 print(WCPprots_10hr_SigDec_Num)
	 #combine proteins that are unique at 0hr and sig dec at 10hr
	 WCPprots_10hr_Decreased <- c(WCPprots_0hr_Unique, WCPprots_10hr_SigDec)
	 print(WCPprots_10hr_Decreased)
	 WCPprots_10hr_Decreased_Num <- length(WCPprots_10hr_Decreased)
	 print("WCPprots_10hr_Decreased_Num")
	 print(WCPprots_10hr_Decreased_Num)
	 #
	 #subset the KeGG data to get the peptides from these proteins
	 KeGGdata_WCPDecProts <- KeGGdataMerge[KeGGdataMerge$WCPcorrProt%in%WCPprots_10hr_Decreased,]
	 dim(KeGGdata_WCPDecProts)
	 KeGGdata_KeGGPeps_WCPDecProts <- subset(KeGGdata_WCPDecProts, (KeGGdata_WCPDecProts$KeGGInt.6hr_NumReps>=2 | KeGGdata_WCPDecProts$KeGGInt.8hr_NumReps>=2 | KeGGdata_WCPDecProts$KeGGInt.10hr_NumReps>=2))
	 dim(KeGGdata_KeGGPeps_WCPDecProts)
	 write.table(KeGGdata_KeGGPeps_WCPDecProts, file="KeGGdata_0hr10hrWCPdec_6hr8hr10hrKeGGId_KeGGpeptideData.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
	 #
	 
	 ##proteins that are inceased in normalized KeGG abundance
	 #increased by statistical significance
	 KeGGNormFcSigInc_0hr10hr_Data <- subset(KeGGdataMerge, (KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg>1 & KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_onesidePval<0.1))
	 KeGGNormFcSigInc_0hr10hr_Num <- length(row.names(KeGGNormFcSigInc_0hr10hr_Data))
	 print("KeGGNormFcSigInc_0hr10hr_Num")
	 print(KeGGNormFcSigInc_0hr10hr_Num)
	 write.table(KeGGNormFcSigInc_0hr10hr_Data, file="KeGGdata_0hr10hr_KeGGNormFcSigInc_KeGGpeptideData.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
	 #decreased by statistical significance
	 KeGGNormFcSigDec_0hr10hr_Data <- subset(KeGGdataMerge, (KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg < -1 & KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_onesidePval<0.1))
	 KeGGNormFcSigDec_0hr10hr_Num <- length(row.names(KeGGNormFcSigDec_0hr10hr_Data))
	 print("KeGGNormFcSigDec_0hr10hr_Num")
	 print(KeGGNormFcSigDec_0hr10hr_Num)
	 write.table(KeGGNormFcSigDec_0hr10hr_Data, file="KeGGdata_0hr10hr_KeGGNormFcSigDec_KeGGpeptideData.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
	 
	 
	 #increase or decrease by statistical signficance
	 KeGGNormFcWcpFc_SigDiff_0hr10hr <- row.names(subset(KeGGdataMerge, (KeGGdataMerge$WCP_iBAQ_0hr10hr_pval<0.05 & !is.na(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg))))
	 KeGGNormFcWcpFc_SigDiff_0hr10hr_Num <- length(KeGGNormFcWcpFc_SigDiff_0hr10hr)
	 print("KeGGNormFcWcpFc_SigDiff_0hr10hr_Num")
	 print(KeGGNormFcWcpFc_SigDiff_0hr10hr_Num)
	 #
	 
	 ##classify proteins as degradative or non-degradative targets of ubiquitination based on 6hr WCP and non-norm KeGG change
	 KeGGdataMerge_6hr_DegProt <- subset(KeGGdataMerge, (KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg<WCP_0hr10hr_Fc_Avg-(WCP_0hr10hr_Fc_Sd) & KeGGdataMerge$KeGG_Int_0hr6hr_Fc>KeGGNormFc_0hr6hr_Avg+KeGGNormFc_0hr6hr_Sd))
	 DegProtPeps_6hr <- row.names(KeGGdataMerge_6hr_DegProt)
	 KeGGdataMerge_6hr_NondegProt <- subset(KeGGdataMerge, ((KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg>WCP_0hr10hr_Fc_Avg-(WCP_0hr10hr_Fc_Sd) & KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg<WCP_0hr10hr_Fc_Avg+(WCP_0hr10hr_Fc_Sd)) & KeGGdataMerge$KeGG_Int_0hr6hr_Fc>KeGGNormFc_0hr6hr_Avg+KeGGNormFc_0hr6hr_Sd))
	 NondegProtPeps_6hr <- row.names(KeGGdataMerge_6hr_NondegProt)
	 #
	 ##classify proteins as degradative or non-degradative targets of ubiquitination based on 8hr WCP and non-norm KeGG change
	 KeGGdataMerge_8hr_DegProt <- subset(KeGGdataMerge, (KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg<WCP_0hr10hr_Fc_Avg-(WCP_0hr10hr_Fc_Sd) & KeGGdataMerge$KeGG_Int_0hr8hr_Fc>KeGGNormFc_0hr8hr_Avg+KeGGNormFc_0hr8hr_Sd))
	 DegProtPeps_8hr <- row.names(KeGGdataMerge_8hr_DegProt)
	 KeGGdataMerge_8hr_NondegProt <- subset(KeGGdataMerge, ((KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg>WCP_0hr10hr_Fc_Avg-(WCP_0hr10hr_Fc_Sd) & KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg<WCP_0hr10hr_Fc_Avg+(WCP_0hr10hr_Fc_Sd)) & KeGGdataMerge$KeGG_Int_0hr8hr_Fc>KeGGNormFc_0hr8hr_Avg+KeGGNormFc_0hr8hr_Sd))
	 NondegProtPeps_8hr <- row.names(KeGGdataMerge_8hr_NondegProt)
	 #
	 
	 ##classify proteins as degradative or non-degradative targets of ubituitination based on 10hr WCP and norm KeGG change
	 KeGGdataMerge_10hr_DegProt_Data <- subset(KeGGdataMerge, (KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg<WCP_0hr10hr_Fc_Avg-(WCP_0hr10hr_Fc_Sd) & KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg>KeGGNormFc_0hr10hr_Avg+(KeGGNormFc_0hr10hr_Sd)))
	 DegProtPeps_10hr <- row.names(KeGGdataMerge_10hr_DegProt_Data)
	 DegProtPeps_10hr_Num <- length(DegProtPeps_10hr)
	 print("DegProtPeps_10hr_Num")
	 print(DegProtPeps_10hr_Num)
	 write.table(KeGGdataMerge_10hr_DegProt_Data, file="KeGGdata_normKeGGFcWcpFc_DegProt2sd_KeGGpeptideData.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
	 KeGGdataMerge_10hr_NondegProt_Data <- subset(KeGGdataMerge, ((KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg>WCP_0hr10hr_Fc_Avg-(WCP_0hr10hr_Fc_Sd) & KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg<WCP_0hr10hr_Fc_Avg+(WCP_0hr10hr_Fc_Sd)) & KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg>KeGGNormFc_0hr10hr_Avg+(KeGGNormFc_0hr10hr_Sd)))
	 NondegProtPeps_10hr <- row.names(KeGGdataMerge_10hr_NondegProt_Data)
	 NondegProtPeps_10hr_Num <- length(NondegProtPeps_10hr)
	 print("NondegProtPeps_10hr_Num")
	 print(NondegProtPeps_10hr_Num)
	 write.table(KeGGdataMerge_10hr_NondegProt_Data, file="KeGGdata_normKeGGFcWcpFc_NondegProt2sd_KeGGpeptideData.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
	 
	 #get degraded and non-degraded peptides in KeGG id at 6 or 8 hours that aren't at 10hr
	 DegProtPeps_6hr_Ex10hr <- setdiff(DegProtPeps_6hr, DegProtPeps_10hr)
	 DegProtPeps_6hr_Ex10hr_Num <- length(DegProtPeps_6hr_Ex10hr)
	 print("DegProtPeps_6hr_Ex10hr_Num")
	 print(DegProtPeps_6hr_Ex10hr_Num)
	 DegProtPeps_8hr_Ex10hr <- setdiff(DegProtPeps_8hr, DegProtPeps_10hr)
	 DegProtPeps_8hr_Ex10hr_Num <- length(DegProtPeps_8hr_Ex10hr)
	 print("DegProtPeps_8hr_Ex10hr_Num")
	 print(DegProtPeps_8hr_Ex10hr_Num)
	 NondegProtPeps_6hr_Ex10hr <- setdiff(NondegProtPeps_6hr, NondegProtPeps_10hr)
	 NondegProtPeps_6hr_Ex10hr_Num <- length(NondegProtPeps_6hr_Ex10hr)
	 print("NondegProtPeps_6hr_Ex10hr_Num")
	 print(NondegProtPeps_6hr_Ex10hr_Num)
	 NondegProtPeps_8hr_Ex10hr <- setdiff(NondegProtPeps_8hr, NondegProtPeps_10hr)
	 NondegProtPeps_8hr_Ex10hr_Num <- length(NondegProtPeps_8hr_Ex10hr)
	 print("NondegProtPeps_8hr_Ex10hr_Num")
	 print(NondegProtPeps_8hr_Ex10hr_Num)
	 DegProtPeps_6hr8hr_union <- union(DegProtPeps_6hr_Ex10hr,DegProtPeps_8hr_Ex10hr)
	 DegProtPeps_6hr8hr_union_Num <- length(DegProtPeps_6hr8hr_union)
	 print("DegProtPeps_6hr8hr_union_Num")
	 print(DegProtPeps_6hr8hr_union_Num)
	 NondegProtPeps_6hr8hr_union <- union(NondegProtPeps_6hr_Ex10hr,NondegProtPeps_8hr_Ex10hr)
	 NondegProtPeps_6hr8hr_union_Num <- length(NondegProtPeps_6hr8hr_union)
	 print("NondegProtPeps_6hr8hr_union_Num")
	 print(NondegProtPeps_6hr8hr_union_Num)
	 #get degraded peptides that are increased in KeGG but not identified in 10hr wcp
	 
	 KeGGdataMerge_0hrWcpUnique_Data <- KeGGdataMerge[KeGGdataMerge$WCPcorrProt%in%WCPprots_0hr_Unique,]
	 dim(KeGGdataMerge_0hrWcpUnique_Data)
	 KeGGdataMerge_0hrWcpUnique_DegProt <- subset(KeGGdataMerge_0hrWcpUnique_Data, (KeGGdataMerge_0hrWcpUnique_Data$KeGG_Int_0hr6hr_Fc>0 | KeGGdataMerge_0hrWcpUnique_Data$KeGG_Int_0hr8hr_Fc>0 | KeGGdataMerge_0hrWcpUnique_Data$KeGG_Int_0hr10hr_Fc>0))
	 DegProtPeps_0hrWcpUnique <- row.names(KeGGdataMerge_0hrWcpUnique_DegProt)
	 DegProtPeps_0hrWcpUnique_Num <- length(DegProtPeps_0hrWcpUnique)
	 print("DegProtPeps_0hrWcpUnique_Num")
	 print(DegProtPeps_0hrWcpUnique_Num)
	 
	 ##use subset to find peptides with fold change > 0 and p-value < 0.05
	 ##unpaired
	 #0hr-6hr
	 KeGGSigInc_0hr6hr_unpaired <- row.names(subset(KeGGdataMerge, (KeGGdataMerge$KeGG_Int_0hr6hr_Fc>0 & KeGGdataMerge$KeGG_Int_0hr6hr_unpairedTtestPval<0.05)))
	 KeGGSigInc_0hr6hr_unpaired_Num <- length(KeGGSigInc_0hr6hr_unpaired)
	 print("KeGGSigInc_0hr6hr_unpaired_Num")
	 print(KeGGSigInc_0hr6hr_unpaired_Num)
	 #0hr-8hr
	 KeGGSigInc_0hr8hr_unpaired <- row.names(subset(KeGGdataMerge, (KeGGdataMerge$KeGG_Int_0hr8hr_Fc>0 & KeGGdataMerge$KeGG_Int_0hr8hr_unpairedTtestPval<0.05)))
	 KeGGSigInc_0hr8hr_unpaired_Num <- length(KeGGSigInc_0hr8hr_unpaired)
	 print("KeGGSigInc_0hr8hr_unpaired_Num")
	 print(KeGGSigInc_0hr8hr_unpaired_Num)
	 #0hr-10hr
	 KeGGSigInc_0hr10hr_unpaired <- row.names(subset(KeGGdataMerge, (KeGGdataMerge$KeGG_Int_0hr10hr_Fc>0 & KeGGdataMerge$KeGG_Int_0hr10hr_unpairedTtestPval<0.05)))
	 KeGGSigInc_0hr10hr_unpaired_Num <- length(KeGGSigInc_0hr10hr_unpaired)
	 print("KeGGSigInc_0hr10hr_unpaired_Num")
	 print(KeGGSigInc_0hr10hr_unpaired_Num)
	 #union of significantly increased peptides
	 KeGGSigInc_6hr8hr10hr_unpaired_Union <- Reduce(union, list(KeGGSigInc_0hr6hr_unpaired, KeGGSigInc_0hr8hr_unpaired, KeGGSigInc_0hr10hr_unpaired))
	 KeGGSigInc_6hr8hr10hr_unpaired_Union_Num <- length(KeGGSigInc_6hr8hr10hr_unpaired_Union)
	 print("KeGGSigInc_6hr8hr10hr_unpaired_Union_Num")
	 print(KeGGSigInc_6hr8hr10hr_unpaired_Union_Num)
	 KeGGdataMerge6hr8hr10hrUnpairedSigInc <- KeGGdataMerge[row.names(KeGGdataMerge)%in%KeGGSigInc_6hr8hr10hr_unpaired_Union,]
	 #intersection of significantly increased peptides
	 KeGGSigInc_6hr8hr10hr_unpaired_Int <- Reduce(intersect, list(KeGGSigInc_0hr6hr_unpaired, KeGGSigInc_0hr8hr_unpaired, KeGGSigInc_0hr10hr_unpaired))
	 KeGGSigInc_6hr8hr10hr_unpaired_Int_Num <- length(KeGGSigInc_6hr8hr10hr_unpaired_Int)
	 print("KeGGSigInc_6hr8hr10hr_unpaired_Int_Num")
	 print(KeGGSigInc_6hr8hr10hr_unpaired_Int_Num)
	 ##paired
	 #0hr-6hr
	 KeGGSigInc_0hr6hr_paired <- row.names(subset(KeGGdataMerge, (KeGGdataMerge$KeGG_Int_0hr6hr_Fc>0 & KeGGdataMerge$KeGG_Int_0hr6hr_pairedTtestPval<0.05)))
	 KeGGSigInc_0hr6hr_paired_Num <- length(KeGGSigInc_0hr6hr_paired)
	 print("KeGGSigInc_0hr6hr_paired_Num")
	 print(KeGGSigInc_0hr6hr_paired_Num)
	 #0hr-8hr
	 KeGGSigInc_0hr8hr_paired <- row.names(subset(KeGGdataMerge, (KeGGdataMerge$KeGG_Int_0hr8hr_Fc>0 & KeGGdataMerge$KeGG_Int_0hr8hr_pairedTtestPval<0.05)))
	 KeGGSigInc_0hr8hr_paired_Num <- length(KeGGSigInc_0hr8hr_paired)
	 print("KeGGSigInc_0hr8hr_paired_Num")
	 print(KeGGSigInc_0hr8hr_paired_Num)
	 #0hr-10hr
	 KeGGSigInc_0hr10hr_paired <- row.names(subset(KeGGdataMerge, (KeGGdataMerge$KeGG_Int_0hr10hr_Fc>0 & KeGGdataMerge$KeGG_Int_0hr10hr_pairedTtestPval<0.05)))
	 KeGGSigInc_0hr10hr_paired_Num <- length(KeGGSigInc_0hr10hr_paired)
	 print("KeGGSigInc_0hr10hr_paired_Num")
	 print(KeGGSigInc_0hr10hr_paired_Num)
	 #union of significantly increased peptides
	 KeGGSigInc_6hr8hr10hr_paired_Union <- Reduce(union, list(KeGGSigInc_0hr6hr_paired, KeGGSigInc_0hr8hr_paired, KeGGSigInc_0hr10hr_paired))
	 KeGGSigInc_6hr8hr10hr_paired_Union_Num <- length(KeGGSigInc_6hr8hr10hr_paired_Union)
	 print("KeGGSigInc_6hr8hr10hr_paired_Union_Num")
	 print(KeGGSigInc_6hr8hr10hr_paired_Union_Num)
	 KeGGdataMerge6hr8hr10hrPairedSigInc <- KeGGdataMerge[row.names(KeGGdataMerge)%in%KeGGSigInc_6hr8hr10hr_paired_Union,]
	 #intersection of significantly increased peptides
	 KeGGSigInc_6hr8hr10hr_paired_Int <- Reduce(intersect, list(KeGGSigInc_0hr6hr_paired, KeGGSigInc_0hr8hr_paired, KeGGSigInc_0hr10hr_paired))
	 KeGGSigInc_6hr8hr10hr_paired_Int_Num <- length(KeGGSigInc_6hr8hr10hr_paired_Int)
	 print("KeGGSigInc_6hr8hr10hr_paired_Int_Num")
	 print(KeGGSigInc_6hr8hr10hr_paired_Int_Num)
	 #
	 #write table of significantly increased peptide Ub data
	 write.table(KeGGdataMerge6hr8hr10hrUnpairedSigInc, file="KeGGdata_6hr8hr10hr_UnpairedSigInc_KeGGpeptideData_NotFiltered.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
	 write.table(KeGGdataMerge6hr8hr10hrPairedSigInc, file="KeGGdata_6hr8hr10hr_PairedSigInc_KeGGpeptideData_NotFiltered.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
	 ##end proteins that have significant changes compared to 0 hours
}




##########     PLOTTING     ##########

##plot KeGG peptide fold changes vs WCP protein fold changes
pdf("WCPavgFc_KeGGAvgNormFc_0hr10hr_AllData_0211.pdf")
par(bty="l", lwd=2)
#xmin <- floor(min(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg, na.rm=TRUE))
xmin <- -5
xmax <- ceiling(max(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg, na.rm=TRUE))
ymin <- floor(min(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg, na.rm=TRUE))
ymax <- ceiling(max(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg, na.rm=TRUE))
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
par(new=T)
abline(v=0, col="gray70", lwd=1, lty=2)
par(new=T)
abline(h=0, col="gray70", lwd=1, lty=2)
par(new=T)
plot(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg, KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg,
     pch=21,
	 cex=0.75,
	 #cex=4*(abs(KeGGdataMerge$UbFcAvg_4hr_NormWcpFc)),
	 #col="grey80",
	 #bg="grey80",
	 col=ifelse(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg<WCP_0hr10hr_Fc_Avg-WCP_0hr10hr_Fc_Sd, "red", ifelse(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg>WCP_0hr10hr_Fc_Avg+WCP_0hr10hr_Fc_Sd, "red", "grey70")),
	 bg=ifelse(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg<WCP_0hr10hr_Fc_Avg-WCP_0hr10hr_Fc_Sd, "red", ifelse(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg>WCP_0hr10hr_Fc_Avg+WCP_0hr10hr_Fc_Sd, "red", "grey70")),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
	)
	#make best fit line
	par(new=T)
	abline(fit <- lm(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg ~ KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg), col="blue", lwd=4, lty=3) # regression line (y~x)
	#legend("topleft", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)), text.col="gray")
axis(1,c(seq(xmin,xmax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=1)),col="black",lwd=2,cex.axis=1.5)
mtext("WCP log2 fold change",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)),col="black",lwd=2,cex.axis=1.5)
mtext("normalized KeGG log2 fold change",2,line=2.75,cex=1.5)
dev.off()


if ( FALSE )
{
##plot KeGG peptide fold changes vs WCP protein fold changes
pdf("WCPavgFc_KeGGAvgNormFc_0hr10hr_IncUnchKeGG.pdf")
par(bty="l", lwd=2)
xmin <- -4
xmax <- 1 
ymin <- 0
#ymax <- 7 
#xmin <- floor(min(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg, na.rm=TRUE))
#xmax <- ceiling(max(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg, na.rm=TRUE))
#ymin <- floor(min(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg, na.rm=TRUE))
ymax6hr <- ceiling(max(KeGGdataMerge$KeGG_Int_0hr6hr_FcNorm_Avg, na.rm=TRUE))
ymax8hr <- ceiling(max(KeGGdataMerge$KeGG_Int_0hr6hr_FcNorm_Avg, na.rm=TRUE))
ymax10hr <- ceiling(max(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg, na.rm=TRUE))
ymax <- max(ymax6hr, ymax8hr, ymax10hr)
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
par(new=T)
abline(v=0, col="gray70", lwd=1, lty=2)
par(new=T)
abline(h=0, col="gray70", lwd=1, lty=2)
#par(new=T)
#abline(v=WCP_0hr10hr_Fc_Avg+(WCP_0hr10hr_Fc_Sd), col="gray70", lwd=2, lty=1)
#par(new=t)
#abline(v=WCP_0hr10hr_Fc_Avg-(WCP_0hr10hr_Fc_Sd), col="gray70", lwd=2, lty=1)
#par(new=t)
#abline(h=KeGGNormFc_0hr10hr_Avg+(KeGGNormFc_0hr10hr_Sd), col="gray70", lwd=2, lty=1)
#par(new=t)
#abline(h=KeGGNormFc_0hr10hr_Avg-(KeGGNormFc_0hr10hr_Sd), col="gray70", lwd=2, lty=1)
par(new=T)
plot(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg, KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg,
     pch=21,
	 cex=0.75,
	 #cex=4*(abs(KeGGdataMerge$UbFcAvg_4hr_NormWcpFc)),
	 col="grey80",
	 bg="grey80",
	 #col=ifelse(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_onesidePval<0.05, "red", "grey70"),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
	)
par(new=T)
#text(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg,KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg,labels=row.names(KeGGdataMerge),cex = 0.1)
print("6hour...")
for ( i in 1:length(KeGG_6hr_PredDeg_Peps) ){
	par(new=T)
	diffprot<-KeGG_6hr_PredDeg_Peps[i]
	####################wcp<-KeGGdataMerge[diffprot,"WCP_iBAQ_0hr6hr_Fc_Avg"]
	if ( is.na(KeGGdataMerge[diffprot,"WCP_iBAQ_0hr10hr_Fc_Avg"]) ){
		wcp <- xmin
	}else{
		wcp <- KeGGdataMerge[diffprot,"WCP_iBAQ_0hr10hr_Fc_Avg"]
	}
	wcp<-KeGGdataMerge[diffprot,"WCP_iBAQ_0hr10hr_Fc_Avg"]
	####################if ( is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]) ){
	####################	kegg <- 9
	####################}else{
	####################	kegg <- max(c(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]), na.rm=TRUE)
	####################}
	if ( is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"]) ){
		kegg <- ymax
	}else{
		kegg <- KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"]
	}
	if ( is.na(wcp) | is.na(kegg) ){
		print(wcp)
		print(kegg)
	}
	points(wcp,kegg,
	       pch=21,
		   cex=0.75,
		   col="firebrick",
		   bg="firebrick"
		  )
}
print("8hr...")
for ( i in 1:length(KeGG_6hr8hr_PredDeg_Peps) ){
	par(new=T)
	diffprot<-KeGG_6hr8hr_PredDeg_Peps[i]
	####################if ( is.na(KeGGdataMerge[diffprot,"WCP_iBAQ_0hr8hr_Fc_Avg"]) ){
	####################	wcp <- -7
	####################}else{
	####################	wcp<-KeGGdataMerge[diffprot,"WCP_iBAQ_0hr8hr_Fc_Avg"]
	####################}
	if ( is.na(KeGGdataMerge[diffprot,"WCP_iBAQ_0hr10hr_Fc_Avg"]) ){
		wcp <- xmin
	}else{
		wcp <- KeGGdataMerge[diffprot,"WCP_iBAQ_0hr10hr_Fc_Avg"]
	}
	####################if ( is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]) ){
	####################	kegg <- 9
	####################}else{
	####################	kegg <- max(c(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]), na.rm=TRUE)
	####################}
	if ( is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"]) ){
		kegg <- ymax
	}else{
		kegg <- max(c(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"]), na.rm=TRUE)
	}
	if ( is.na(wcp) | is.na(kegg) ){
		print(wcp)
		print(kegg)
	}
	points(wcp,kegg,
	       pch=21,
		   cex=0.75,
		   col="red",
		   bg="red"
		  )
}
print("10hr...")
for ( i in 1:length(KeGG_6hr8hr10hr_PredDeg_Peps) ){
	par(new=T)
	diffprot<-KeGG_6hr8hr10hr_PredDeg_Peps[i]
	if ( is.na(KeGGdataMerge[diffprot,"WCP_iBAQ_0hr10hr_Fc_Avg"]) ){
		wcp <- xmin
	}else{
		wcp<-KeGGdataMerge[diffprot,"WCP_iBAQ_0hr10hr_Fc_Avg"]
	}
	if ( is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]) ){
		kegg <- ymax
	}else{
		kegg <- max(c(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]), na.rm=TRUE)
	}
	if ( is.na(wcp) | is.na(kegg) ){
		print(wcp)
		print(kegg)
	}
	points(wcp,kegg,
	       pch=21,
		   cex=0.75,
		   col="salmon",
		   bg="salmon"
		  )
}
print("non-deg...")
for ( i in 1:length(KeGG_AllPredUnch_Peps) ){
	par(new=T)
	diffprot<-KeGG_AllPredUnch_Peps[i]
	wcp<-KeGGdataMerge[diffprot,"WCP_iBAQ_0hr10hr_Fc_Avg"]
	if ( is.na(wcp) ){
		print("wcp")
		print(wcp)
	}
	if ( is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]) ){
		kegg <- ymax
	}else{
		kegg <- max(c(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]), na.rm=TRUE)
	}
	if ( is.na(kegg) ){
		print("kegg")
		print(kegg)
	}
	points(wcp,kegg,
	       pch=21,
		   cex=0.75,
		   col="olivedrab",
		   bg="olivedrab"
		  )
}
axis(1,c(seq(xmin,xmax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=1)),col="black",lwd=2,cex.axis=1.5)
mtext("WCP log2 fold change",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)),col="black",lwd=2,cex.axis=1.5)
mtext("normalized KeGG log2 fold change",2,line=2.75,cex=1.5)
dev.off()
}


if ( FALSE )
{
print("scan")
RNAbinding <- scan("NonDegRNAbindingProteins.txt", character(), quote = "")
length(RNAbinding)
#####KeGGdataMerge_ProteinData_RNAbind <- KeGGdataMerge_ProteinData[KeGGdataMerge_ProteinData$Gene.names%in%RNAbinding,]
KeGGdataMerge_ProteinData_RNAbind <- KeGGdataMerge_ProteinData[row.names(KeGGdataMerge_ProteinData)%in%RNAbinding,]
dim(KeGGdataMerge_ProteinData_RNAbind)
Nondeg_RNAbind_prots <- row.names(KeGGdataMerge_ProteinData_RNAbind)
length(Nondeg_RNAbind_prots)

UbIncWcpDecProts <- scan("UbiquitinatedAndDecreasedProts.txt", character(), quote = "")
}

##plot KeGG peptide fold changes vs WCP protein fold changes
#pdf("WCPavgFc_KeGGProteinAvgNormFc_0hr10hr_IncUnchKeGGProt.pdf")
####tiff("WCPavgFc_KeGGProteinAvgNormFc_0hr10hr_IncUnchKeGGProt_1218.tiff", width = 7, height = 7, units = 'in', res = 300, compression = 'none')
tiff("WCPavgFc_KeGGProteinAvgNormFc_0hr10hr_IncUnchKeGGProt_0211.tiff", width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(bty="l", lwd=2)
#xmin <- floor(min(KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg, na.rm=TRUE))
#xmax <- ceiling(max(KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg, na.rm=TRUE))
xmin <- -4
xmax <- 1
#ymin <- floor(min(KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm, na.rm=TRUE))
ymax <- ceiling(max(KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm, na.rm=TRUE))
#ymin <- -9
#ymax <- 0
ymin <- 0
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
#par(new=T)
#abline(v=0, col="gray70", lwd=1, lty=2)
#par(new=T)
#abline(h=0, col="gray70", lwd=1, lty=2)
par(new=T)
abline(v=WCP_0hr10hr_Fc_Avg+(WCP_0hr10hr_Fc_Sd), col="gray70", lwd=2, lty=2)
par(new=t)
abline(v=WCP_0hr10hr_Fc_Avg-(WCP_0hr10hr_Fc_Sd), col="gray70", lwd=2, lty=2)
par(new=t)
#abline(h=KeGG_0hr10hr_ProtFc_Avg+(KeGG_0hr10hr_ProtFc_Sd), col="gray70", lwd=2, lty=2)
abline(h= 1, col="gray70", lwd=2, lty=2)
par(new=t)
#abline(h=KeGG_0hr10hr_ProtFc_Avg-(KeGG_0hr10hr_ProtFc_Sd), col="gray70", lwd=2, lty=2)
abline(h= -1, col="gray70", lwd=2, lty=2)
par(new=T)
wcpthresholdinc<-WCP_0hr10hr_Fc_Avg+WCP_0hr10hr_Fc_Sd
wcpthresholddec<-WCP_0hr10hr_Fc_Avg-WCP_0hr10hr_Fc_Sd
#keggthresholdinc<-KeGG_0hr10hr_ProtFc_Avg+KeGG_0hr10hr_ProtFc_Sd
keggthresholdinc<-1
plot(KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg, KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm,
     pch=21,
	 cex=0.75,
	 #cex=4*(abs(KeGGdataMerge$UbFcAvg_4hr_NormWcpFc)),
	 #col="grey60",
	 #bg="grey60",
	 
	 col=ifelse((KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholddec & KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "red", ifelse(((KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "dodgerblue", ifelse((KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholdinc & KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "salmon", ifelse(((KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "skyblue", "grey60")))),
	 bg=ifelse((KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholddec & KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "red", ifelse(((KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "dodgerblue", ifelse((KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholdinc & KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "salmon", ifelse(((KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "skyblue", "grey60")))),
	 #bg=ifelse((KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholddec & KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "red", ifelse(((KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "dodgerblue", "grey60")),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
	)
#text(KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg, KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm, labels=row.names(KeGGdataMerge_ProteinData), cex = 0.1)
#text(KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg, KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm, labels=KeGGdataMerge_ProteinData$Gene.names, cex = 0.1)

#par(new=T)
#RAD50wcp<-KeGGdataMerge_ProteinData["Q92878","WCP_iBAQ_0hr10hr_Fc_Avg"]
#RAD50kegg<-KeGGdataMerge_ProteinData["Q92878","KeGG_Prot_0hr10hr_Fc_Norm"]
#plot(RAD50wcp,RAD50kegg,pch=21,cex=0.75,col="green",bg="green",xlab="",ylab="",xlim=c(xmin,xmax), ylim=c(ymin,ymax),xaxt='n', yaxt='n')
#
#par(new=T)
#MRE11wcp<-KeGGdataMerge_ProteinData["P49959","WCP_iBAQ_0hr10hr_Fc_Avg"]
#MRE11kegg<-KeGGdataMerge_ProteinData["P49959","KeGG_Prot_0hr10hr_Fc_Norm"]
#plot(MRE11wcp,MRE11kegg,pch=21,cex=0.75,col="yellow",bg="yellow",xlab="",ylab="",xlim=c(xmin,xmax), ylim=c(ymin,ymax),xaxt='n', yaxt='n')
#
#par(new=T)
#RALYwcp<-KeGGdataMerge_ProteinData["Q9UKM9","WCP_iBAQ_0hr10hr_Fc_Avg"]
#RALYkegg<-KeGGdataMerge_ProteinData["Q9UKM9","KeGG_Prot_0hr10hr_Fc_Norm"]
#plot(RALYwcp,RALYkegg,pch=21,cex=0.75,col="purple",bg="purple",xlab="",ylab="",xlim=c(xmin,xmax), ylim=c(ymin,ymax),xaxt='n', yaxt='n')
#
#par(new=T)
#RNPCwcp<-KeGGdataMerge_ProteinData["P07910","WCP_iBAQ_0hr10hr_Fc_Avg"]
#RNPCkegg<-KeGGdataMerge_ProteinData["P07910","KeGG_Prot_0hr10hr_Fc_Norm"]
#plot(RNPCwcp,RNPCkegg,pch=21,cex=0.75,col="orange",bg="orange",xlab="",ylab="",xlim=c(xmin,xmax), ylim=c(ymin,ymax),xaxt='n', yaxt='n')

#text(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg,KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg,labels=row.names(KeGGdataMerge),cex = 0.1)
axis(1,c(seq(xmin,xmax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=1)),col="black",lwd=2,cex.axis=1.5)
mtext("WCP log2 fold change",3,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)),col="black",lwd=2,cex.axis=1.5)
mtext("normalized KeGG log2 fold change",2,line=2.75,cex=1.5)
dev.off()

KeGGdataMerge_ProteinData_Unique0hr <- subset(KeGGdataMerge_ProteinData, KeGGdataMerge_ProteinData$Unique0hr==1)
KeGGdataMerge_ProteinData_Unique10hr <- subset(KeGGdataMerge_ProteinData, KeGGdataMerge_ProteinData$Unique10hr==1)
KeGGdataMerge_ProteinData_NonUnique <- subset(KeGGdataMerge_ProteinData, (is.na(KeGGdataMerge_ProteinData$Unique0hr) & is.na(KeGGdataMerge_ProteinData$Unique10hr)))
##plot KeGG peptide fold changes vs WCP protein fold changes
pdf("WCPavgFc_KeGGProteinAvgNormFc_0hr10hr_IncUnchKeGGProt_adduniques.pdf")
#####tiff("WCPavgFc_KeGGProteinAvgNormFc_0hr10hr_IncUnchKeGGProt_adduniques_0211.tiff", width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(bty="l", lwd=2)
xmin <- floor(min(KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg, na.rm=TRUE))
xmax <- ceiling(max(KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg, na.rm=TRUE))
xmin <- -4
xmax <- 1
#ymin <- floor(min(KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm, na.rm=TRUE))
ymax <- ceiling(max(KeGGdataMerge_ProteinData_NonUnique$KeGG_Prot_0hr10hr_Fc_Norm, na.rm=TRUE))+1
#ymin <- -9
#ymax <- 0
ymin <- 0
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
#par(new=T)
#abline(v=0, col="gray70", lwd=1, lty=2)
#par(new=T)
#abline(h=0, col="gray70", lwd=1, lty=2)
par(new=T)
abline(v=WCP_0hr10hr_Fc_Avg+(WCP_0hr10hr_Fc_Sd), col="gray70", lwd=2, lty=2)
par(new=t)
abline(v=WCP_0hr10hr_Fc_Avg-(WCP_0hr10hr_Fc_Sd), col="gray70", lwd=2, lty=2)
par(new=t)
#abline(h=KeGG_0hr10hr_ProtFc_Avg+(KeGG_0hr10hr_ProtFc_Sd), col="gray70", lwd=2, lty=2)
abline(h= 1, col="gray70", lwd=2, lty=2)
par(new=t)
#abline(h=KeGG_0hr10hr_ProtFc_Avg-(KeGG_0hr10hr_ProtFc_Sd), col="gray70", lwd=2, lty=2)
abline(h= -1, col="gray70", lwd=2, lty=2)
par(new=T)
wcpthresholdinc<-WCP_0hr10hr_Fc_Avg+WCP_0hr10hr_Fc_Sd
wcpthresholddec<-WCP_0hr10hr_Fc_Avg-WCP_0hr10hr_Fc_Sd
#keggthresholdinc<-KeGG_0hr10hr_ProtFc_Avg+KeGG_0hr10hr_ProtFc_Sd
keggthresholdinc<-1
plot(KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg, KeGGdataMerge_ProteinData_NonUnique$KeGG_Prot_0hr10hr_Fc_Norm,
     pch=21,
	 cex=0.75,
	 #cex=4*(abs(KeGGdataMerge$UbFcAvg_4hr_NormWcpFc)),
	 #col="grey60",
	 #bg="grey60",
	 
	 col=ifelse((KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholddec & KeGGdataMerge_ProteinData_NonUnique$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "red", ifelse(((KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData_NonUnique$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "dodgerblue", ifelse((KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholdinc & KeGGdataMerge_ProteinData_NonUnique$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "salmon", ifelse(((KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData_NonUnique$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "skyblue", "grey60")))),
	 bg=ifelse((KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholddec & KeGGdataMerge_ProteinData_NonUnique$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "red", ifelse(((KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData_NonUnique$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "dodgerblue", ifelse((KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholdinc & KeGGdataMerge_ProteinData_NonUnique$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "salmon", ifelse(((KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData_NonUnique$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "skyblue", "grey60")))),
	 #bg=ifelse((KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholddec & KeGGdataMerge_ProteinData_NonUnique$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "red", ifelse(((KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData_NonUnique$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData_NonUnique$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "dodgerblue", "grey60")),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
	)
par(new=T)
KeGGdataMerge_ProteinData_Unique0hr$KeGG_Prot_0hr10hr_Fc_Norm<-floor(min(KeGGdataMerge_ProteinData_NonUnique$KeGG_Prot_0hr10hr_Fc_Norm, na.rm=TRUE))
plot(KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg, KeGGdataMerge_ProteinData_Unique0hr$KeGG_Prot_0hr10hr_Fc_Norm,
     pch=21,
	 cex=0.75,
	 #cex=4*(abs(KeGGdataMerge$UbFcAvg_4hr_NormWcpFc)),
	 #col="grey60",
	 #bg="grey60",
	 
	 col=ifelse((KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholddec & KeGGdataMerge_ProteinData_Unique0hr$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "red", ifelse(((KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData_Unique0hr$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "dodgerblue", ifelse((KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholdinc & KeGGdataMerge_ProteinData_Unique0hr$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "salmon", ifelse(((KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData_Unique0hr$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "skyblue", "grey60")))),
	 bg=ifelse((KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholddec & KeGGdataMerge_ProteinData_Unique0hr$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "red", ifelse(((KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData_Unique0hr$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "dodgerblue", ifelse((KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholdinc & KeGGdataMerge_ProteinData_Unique0hr$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "salmon", ifelse(((KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData_Unique0hr$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "skyblue", "grey60")))),
	 #bg=ifelse((KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholddec & KeGGdataMerge_ProteinData_Unique0hr$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "red", ifelse(((KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData_Unique0hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData_Unique0hr$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "dodgerblue", "grey60")),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
	)
par(new=T)
KeGGdataMerge_ProteinData_Unique10hr$KeGG_Prot_0hr10hr_Fc_Norm<-ceiling(max(KeGGdataMerge_ProteinData_NonUnique$KeGG_Prot_0hr10hr_Fc_Norm, na.rm=TRUE))
plot(KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg, KeGGdataMerge_ProteinData_Unique10hr$KeGG_Prot_0hr10hr_Fc_Norm,
     pch=21,
	 cex=0.75,
	 #cex=4*(abs(KeGGdataMerge$UbFcAvg_4hr_NormWcpFc)),
	 #col="grey60",
	 #bg="grey60",
	 
	 col=ifelse((KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholddec & KeGGdataMerge_ProteinData_Unique10hr$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "red", ifelse(((KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData_Unique10hr$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "dodgerblue", ifelse((KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholdinc & KeGGdataMerge_ProteinData_Unique10hr$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "salmon", ifelse(((KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData_Unique10hr$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "skyblue", "grey60")))),
	 bg=ifelse((KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholddec & KeGGdataMerge_ProteinData_Unique10hr$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "red", ifelse(((KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData_Unique10hr$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "dodgerblue", ifelse((KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholdinc & KeGGdataMerge_ProteinData_Unique10hr$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "salmon", ifelse(((KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData_Unique10hr$KeGG_Prot_0hr10hr_Fc_Norm < -1*keggthresholdinc), "skyblue", "grey60")))),
	 #bg=ifelse((KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholddec & KeGGdataMerge_ProteinData_Unique10hr$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "red", ifelse(((KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg > wcpthresholddec & KeGGdataMerge_ProteinData_Unique10hr$WCP_iBAQ_0hr10hr_Fc_Avg < wcpthresholdinc) & KeGGdataMerge_ProteinData_Unique10hr$KeGG_Prot_0hr10hr_Fc_Norm > keggthresholdinc), "dodgerblue", "grey60")),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
	)

#text(KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg, KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm, labels=row.names(KeGGdataMerge_ProteinData), cex = 0.1)
#text(KeGGdataMerge_ProteinData$WCP_iBAQ_0hr10hr_Fc_Avg, KeGGdataMerge_ProteinData$KeGG_Prot_0hr10hr_Fc_Norm, labels=KeGGdataMerge_ProteinData$Gene.names, cex = 0.1)

#par(new=T)
#RAD50wcp<-KeGGdataMerge_ProteinData["Q92878","WCP_iBAQ_0hr10hr_Fc_Avg"]
#RAD50kegg<-KeGGdataMerge_ProteinData["Q92878","KeGG_Prot_0hr10hr_Fc_Norm"]
#plot(RAD50wcp,RAD50kegg,pch=21,cex=0.75,col="green",bg="green",xlab="",ylab="",xlim=c(xmin,xmax), ylim=c(ymin,ymax),xaxt='n', yaxt='n')
#
#par(new=T)
#MRE11wcp<-KeGGdataMerge_ProteinData["P49959","WCP_iBAQ_0hr10hr_Fc_Avg"]
#MRE11kegg<-KeGGdataMerge_ProteinData["P49959","KeGG_Prot_0hr10hr_Fc_Norm"]
#plot(MRE11wcp,MRE11kegg,pch=21,cex=0.75,col="yellow",bg="yellow",xlab="",ylab="",xlim=c(xmin,xmax), ylim=c(ymin,ymax),xaxt='n', yaxt='n')
#
#par(new=T)
#RALYwcp<-KeGGdataMerge_ProteinData["Q9UKM9","WCP_iBAQ_0hr10hr_Fc_Avg"]
#RALYkegg<-KeGGdataMerge_ProteinData["Q9UKM9","KeGG_Prot_0hr10hr_Fc_Norm"]
#plot(RALYwcp,RALYkegg,pch=21,cex=0.75,col="purple",bg="purple",xlab="",ylab="",xlim=c(xmin,xmax), ylim=c(ymin,ymax),xaxt='n', yaxt='n')
#
#par(new=T)
#RNPCwcp<-KeGGdataMerge_ProteinData["P07910","WCP_iBAQ_0hr10hr_Fc_Avg"]
#RNPCkegg<-KeGGdataMerge_ProteinData["P07910","KeGG_Prot_0hr10hr_Fc_Norm"]
#plot(RNPCwcp,RNPCkegg,pch=21,cex=0.75,col="orange",bg="orange",xlab="",ylab="",xlim=c(xmin,xmax), ylim=c(ymin,ymax),xaxt='n', yaxt='n')

#text(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg,KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg,labels=row.names(KeGGdataMerge),cex = 0.1)
axis(1,c(seq(xmin,xmax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=1)),col="black",lwd=2,cex.axis=1.5)
mtext("WCP log2 fold change",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)),col="black",lwd=2,cex.axis=1.5)
mtext("normalized KeGG log2 fold change",2,line=2.75,cex=1.5)
dev.off()

stop("DONE")

if ( FALSE )
{
##plot KeGG peptide fold changes vs WCP protein fold changes
print("scan")
RNAbinding <- scan("NonDegRNAbindingProteins.txt", character(), quote = "")
length(RNAbinding)
KeGGdataMerge_RNAbind <- KeGGdataMerge[KeGGdataMerge$Gene.names%in%RNAbinding,]
dim(KeGGdataMerge_RNAbind)
Nondeg_RNAbind_peps <- row.names(KeGGdataMerge_RNAbind)
length(Nondeg_RNAbind_peps)
Nondeg_RNAbind_peps <- intersect(KeGG_AllPredUnch_Peps,Nondeg_RNAbind_peps)
print("makefig")
pdf("WCPavgFc_KeGGAvgNormFc_0hr10hr_RNAbinding_0613.pdf")
par(bty="l", lwd=2)
xmin <- -7
xmax <- 1 
ymin <- 0
#ymax <- 7 
#xmin <- floor(min(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg, na.rm=TRUE))
#xmax <- ceiling(max(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg, na.rm=TRUE))
#ymin <- floor(min(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg, na.rm=TRUE))
ymax <- ceiling(max(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg, na.rm=TRUE))
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
par(new=T)
abline(v=0, col="gray70", lwd=1, lty=2)
par(new=T)
abline(h=0, col="gray70", lwd=1, lty=2)
#par(new=T)
#abline(v=WCP_0hr10hr_Fc_Avg+(WCP_0hr10hr_Fc_Sd), col="gray70", lwd=2, lty=1)
#par(new=t)
#abline(v=WCP_0hr10hr_Fc_Avg-(WCP_0hr10hr_Fc_Sd), col="gray70", lwd=2, lty=1)
#par(new=t)
#abline(h=KeGGNormFc_0hr10hr_Avg+(KeGGNormFc_0hr10hr_Sd), col="gray70", lwd=2, lty=1)
#par(new=t)
#abline(h=KeGGNormFc_0hr10hr_Avg-(KeGGNormFc_0hr10hr_Sd), col="gray70", lwd=2, lty=1)
par(new=t)
plot(KeGGdataMerge$WCP_iBAQ_0hr10hr_Fc_Avg, KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg,
     pch=21,
	 cex=0.75,
	 #cex=4*(abs(KeGGdataMerge$UbFcAvg_4hr_NormWcpFc)),
	 col="grey80",
	 bg="grey80",
	 #col=ifelse(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_onesidePval<0.05, "red", "grey70"),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
	)
par(new=T)
#text(KeGGdddataMerge$WCP_iBAQ_0hr10hr_Fc_Avg,KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg,labels=row.names(KeGGdataMerge),cex = 0.1)
print("6hour...")
for ( i in 1:length(KeGG_6hr_PredDeg_Peps) ){
	par(new=T)
	diffprot<-KeGG_6hr_PredDeg_Peps[i]
	wcp<-KeGGdataMerge[diffprot,"WCP_iBAQ_0hr6hr_Fc_Avg"]
	if ( is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]) ){
		kegg <- 9
	}else{
		kegg <- max(c(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]), na.rm=TRUE)
	}
	if ( is.na(wcp) | is.na(kegg) ){
		print(wcp)
		print(kegg)
	}
	points(wcp,kegg,
	       pch=21,
		   cex=0.75,
		   col="firebrick",
		   bg="firebrick"
		  )
}
print("8hr...")
for ( i in 1:length(KeGG_6hr8hr_PredDeg_Peps) ){
	par(new=T)
	diffprot<-KeGG_6hr8hr_PredDeg_Peps[i]
	wcp<-KeGGdataMerge[diffprot,"WCP_iBAQ_0hr8hr_Fc_Avg"]
	if ( is.na(wcp) ){
		#wcp<-KeGGdataMerge[diffprot,"WCP_iBAQ_0hr6hr_Fc_Avg"]
		wcp <- -7
	}
	if ( is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]) ){
		kegg <- 9
	}else{
		kegg <- max(c(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]), na.rm=TRUE)
	}
	if ( is.na(wcp) | is.na(kegg) ){
		print(wcp)
		print(kegg)
	}
	points(wcp,kegg,
	       pch=21,
		   cex=0.75,
		   col="red",
		   bg="red"
		  )
}
print("10hr...")
for ( i in 1:length(KeGG_6hr8hr10hr_PredDeg_Peps) ){
	par(new=T)
	diffprot<-KeGG_6hr8hr10hr_PredDeg_Peps[i]
	wcp<-KeGGdataMerge[diffprot,"WCP_iBAQ_0hr10hr_Fc_Avg"]
	if ( is.na(wcp) ){
		#wcp<-KeGGdataMerge[diffprot,"WCP_iBAQ_0hr8hr_Fc_Avg"]
		wcp <- -7
	}
	if ( is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]) ){
		kegg <- 9
	}else{
		kegg <- max(c(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]), na.rm=TRUE)
	}
	if ( is.na(wcp) | is.na(kegg) ){
		print(wcp)
		print(kegg)
	}
	points(wcp,kegg,
	       pch=21,
		   cex=0.75,
		   col="salmon",
		   bg="salmon"
		  )
}
print("non-deg...")
for ( i in 1:length(KeGG_AllPredUnch_Peps) ){
	par(new=T)
	diffprot<-KeGG_AllPredUnch_Peps[i]
	wcp<-KeGGdataMerge[diffprot,"WCP_iBAQ_0hr10hr_Fc_Avg"]
	if ( is.na(wcp) ){
		print("wcp")
		print(wcp)
	}
	if ( is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]) ){
		kegg <- 9
	}else{
		kegg <- max(c(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]), na.rm=TRUE)
	}
	if ( is.na(kegg) ){
		print("kegg")
		print(kegg)
	}
	points(wcp,kegg,
	       pch=21,
		   cex=0.75,
		   col="olivedrab",
		   bg="olivedrab"
		  )
}
for ( i in 1:length(Nondeg_RNAbind_peps) ){
	par(new=T)
	diffprot<-Nondeg_RNAbind_peps[i]
	wcp<-KeGGdataMerge[diffprot,"WCP_iBAQ_0hr10hr_Fc_Avg"]
	if ( is.na(wcp) ){
		print("wcp")
		print(wcp)
	}
	if ( is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"]) & is.na(KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]) ){
		kegg <- 9
	}else{
		kegg <- max(c(KeGGdataMerge[diffprot,"KeGG_Int_0hr6hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr8hr_FcNorm_Avg"], KeGGdataMerge[diffprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]), na.rm=TRUE)
	}
	if ( is.na(kegg) ){
		print("kegg")
		print(kegg)
	}
	points(wcp,kegg,
	       pch=21,
		   cex=0.75,
		   col="dodgerblue",
		   bg="dodgerblue"
		  )
}
axis(1,c(seq(xmin,xmax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=1)),col="black",lwd=2,cex.axis=1.5)
mtext("WCP log2 fold change",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)),col="black",lwd=2,cex.axis=1.5)
mtext("KeGG log2 fold change",2,line=2.75,cex=1.5)
dev.off()

}

##volcano plot
#xmin <- floor(min(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg, na.rm=TRUE))
#xmax <- ceiling(max(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg, na.rm=TRUE))
ymin <- floor(min(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_onesidePval_negLog10, na.rm=TRUE))
ymax <- ceiling(max(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_onesidePval_negLog10, na.rm=TRUE))
xmin <- -6
xmax <- 8
#ymin <- 0
#ymax <- 4.5
#pdf("KeGGdata_0hr10hr_volcanoplot.pdf")
#####tiff("KeGGdata_0hr10hr_volcanoplot.tiff", width = 7, height = 7, units = 'in', res = 300, compression = 'none')
tiff("KeGGdata_0hr10hr_volcanoplot_0211.tiff", width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(bty="l", lwd=2)
plot(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg,
     KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_onesidePval_negLog10,
	 pch=19,
	 cex=1,
	 #col="gray50",
	 col=ifelse(KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_onesidePval_negLog10<1, "grey70", ifelse((KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg > 1 | KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_Avg < -1), "grey30", "grey70")),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n',yaxt='n'
	)
#targets<-c("MRE11", "RAD50", "NBN", "LIG4")
targets<-c("P49959_434", "P49959_442", "P49959_496", "P49959_510")
for( i in 1:length(targets) ){
     par(new=T)
	 currprot<-targets[i]
	 fc<-KeGGdataMerge[currprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]
	 pval<-KeGGdataMerge[currprot,"KeGG_Int_0hr10hr_FcNorm_onesidePval_negLog10"]
	 plot(fc,pval,
	      pch=19,
		  cex=1.5,
		  col="red",
		  xlab="",ylab="",
		  xlim=c(xmin,xmax), ylim=c(ymin,ymax),
		  xaxt='n',yaxt='n'
		 )
}
targets<-c("Q92878_1103", "Q92878_1126", "Q92878_1178")
for( i in 1:length(targets) ){
     par(new=T)
	 currprot<-targets[i]
	 fc<-KeGGdataMerge[currprot,"KeGG_Int_0hr10hr_FcNorm_Avg"]
	 pval<-KeGGdataMerge[currprot,"KeGG_Int_0hr10hr_FcNorm_onesidePval_negLog10"]
	 plot(fc,pval,
	      pch=19,
		  cex=1.5,
		  col="red",
		  xlab="",ylab="",
		  xlim=c(xmin,xmax), ylim=c(ymin,ymax),
		  xaxt='n',yaxt='n'
		 )
}
#label and format axes
#axis(1, c(seq(xmin,xmax,by=1)), labels=F, col="black",cex.axis=1, tck=-0.01)
axis(1, c(seq(xmin,xmax,by=1)), labels=F, col="black",lwd=2,cex.axis=1.5)
axis(1, c(seq(xmin,xmax,by=2)), labels=T, col="black",lwd=2,cex.axis=1.5)
mtext("normalized log2 fold change",1,line=2.75,cex=1.5)
#axis(2, c(seq(ymin,ymax,by=0.5)), labels=F, col="black",cex.axis=1, tck=-0.01)
axis(2, c(seq(ymin,ymax,by=0.5)), labels=F, col="black",lwd=2,cex.axis=1.5, tck=-0.01)
axis(2, c(seq(ymin,ymax,by=1)), labels=T, col="black",lwd=2,cex.axis=1.5)
mtext("-log10 p-value",2,line=2.75,cex=1.25)
dev.off()


##volcano plot
#xmin <- floor(min(WCPdata$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc, na.rm=TRUE))
#xmax <- ceiling(max(WCPdata$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc, na.rm=TRUE))
#ymin <- floor(min(WCPdata$iBAQ_0hr10hr_Log2MedNorm_negLog10pval, na.rm=TRUE))
ymax <- ceiling(max(WCPdata$iBAQ_0hr10hr_Log2MedNorm_negLog10pval, na.rm=TRUE))
xmin <- -5
xmax <- 5
ymin <- 0
#pdf("WCPdata_0hr10hr_volcanoplot.pdf")
#tiff("WCPdata_0hr10hr_volcanoplot.tiff", width = 7, height = 7, units = 'in', res = 300, compression = 'none')
tiff("WCPdata_0hr10hr_volcanoplot_0211.tiff", width = 7, height = 7, units = 'in', res = 300, compression = 'none')
par(bty="l", lwd=2)
plot(WCPdata$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc,
     WCPdata$iBAQ_0hr10hr_Log2MedNorm_negLog10pval,
	 pch=19,
	 cex=1,
	 #col="gray50",
	 col=ifelse(WCPdata$iBAQ_0hr10hr_Log2MedNorm_negLog10pval<1, "grey70", ifelse((WCPdata$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc > 1 | WCPdata$iBAQ_0hr10hr_Log2MedNorm_Avg_Fc < -1), "grey30", "grey70")),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n',yaxt='n'
	)

#targets<-c("MRE11", "RAD50", "NBN", "LIG4", "BLM")
#targets<-c("P49959", "Q92878", "O60934", "P49917", "P54132")
targets<-c("P49959", "Q92878", "O60934", "P49917")
for( i in 1:length(targets) ){
     par(new=T)
	 currprot<-targets[i]
	 fc<-WCPdata[currprot,"iBAQ_0hr10hr_Log2MedNorm_Avg_Fc"]
	 pval<-WCPdata[currprot,"iBAQ_0hr10hr_Log2MedNorm_negLog10pval"]
	 plot(fc,pval,
	      pch=19,
		  cex=1.5,
		  col="red",
		  xlab="",ylab="",
		  xlim=c(xmin,xmax), ylim=c(ymin,ymax),
		  xaxt='n',yaxt='n'
		 )
}
#sigdec<-WCPprots_10hr_dec_PvalMeanSd
#for( j in 1:length(sigdec) ){
#	par(new=T)
#	currprot<-sigdec[j]
#	currgene<-WCPdata[currprot,"Gene.names.x"]
#	fc<-WCPdata[currprot,"iBAQ_0hr10hr_Log2MedNorm_Avg_Fc"]
#	pval<-WCPdata[currprot,"iBAQ_0hr10hr_Log2MedNorm_negLog10pval"]
#	text(fc,pval,labels=currgene,col="red",cex=0.25)
#}
#keggprot<-KeGGProt_10hr_Inc_WcpProt_10hr_Dec
#for( k in 1:length(keggprot) ){
#	par(new=T)
#	currprot<-keggprot[k]
#	currgene<-WCPdata[currprot,"Gene.names.x"]
#	fc<-WCPdata[currprot,"iBAQ_0hr10hr_Log2MedNorm_Avg_Fc"]
#	pval<-WCPdata[currprot,"iBAQ_0hr10hr_Log2MedNorm_negLog10pval"]
#	plot(fc,pval,
#	     pch=19,
#		 cex=1,
#		 col="salmon",
#		 xlab="",ylab="",
#		 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
#		 xaxt='n',yaxt='n'
#		)
#	text(fc,pval,labels=currgene,col="red",cex=0.25)
#}
#label and format axes
#axis(1, c(seq(xmin,xmax,by=0.5)), labels=F, col="black",cex.axis=1, tck=-0.01)
axis(1, c(seq(xmin,xmax,by=1)), labels=T, col="black",lwd=2,cex.axis=1.5)
mtext("log2 fold change",1,line=2.75,cex=1.5)
#axis(2, c(seq(ymin,ymax,by=0.5)), labels=F, col="black",cex.axis=1, tck=-0.01)
axis(2, c(seq(ymin,ymax,by=0.5)), labels=F, col="black",lwd=2,cex.axis=1.5, tck=-0.01)
axis(2, c(seq(ymin,ymax,by=1)), labels=T, col="black",lwd=2,cex.axis=1.5)
mtext("-log10 p-value",2,line=2.75,cex=1.25)
dev.off()

stop("VOL")

##pairwise comparison of pvalues
#0-6hr
makedensitycoloredscatterplot(KeGGdataMerge$KeGG_Int_0hr6hr_FcRaw_onesidePval,KeGGdataMerge$KeGG_Int_0hr6hr_FcNorm_onesidePval, "NA", "NA", "NA", "NA", "one-side p-value (raw Fc)", "one-side p-value (norm Fc)", "KeGGdata_0hr6hr_RawFcOnesideNormFcOneside_PvalComp_DensityPlot.pdf")
#0-8hr
makedensitycoloredscatterplot(KeGGdataMerge$KeGG_Int_0hr8hr_FcRaw_onesidePval,KeGGdataMerge$KeGG_Int_0hr8hr_FcNorm_onesidePval, "NA", "NA", "NA", "NA", "one-side p-value (raw Fc)", "one-side p-value (norm Fc)", "KeGGdata_0hr8hr_RawFcOnesideNormFcOneside_PvalComp_DensityPlot.pdf")
#0-10hr
makedensitycoloredscatterplot(KeGGdataMerge$KeGG_Int_0hr10hr_FcRaw_onesidePval,KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_onesidePval, "NA", "NA", "NA", "NA", "one-side p-value (raw Fc)", "one-side p-value (norm Fc)", "KeGGdata_0hr10hr_RawFcOnesideNormFcOneside_PvalComp_DensityPlot.pdf")
makedensitycoloredscatterplot(KeGGdataMerge$KeGG_Int_0hr10hr_pairedTtestPval,KeGGdataMerge$KeGG_Int_0hr10hr_FcNorm_onesidePval, "NA", "NA", "NA", "NA", "paired p-value (intensity)", "one-side p-value (norm Fc)", "KeGGdata_0hr10hr_KeGGIntPairedNormFcOneside_PvalComp_DensityPlot.pdf")
makedensitycoloredscatterplot(KeGGdataMerge$KeGG_Int_0hr10hr_pairedTtestPval,KeGGdataMerge$KeGG_Int_0hr10hr_FcRaw_onesidePval, "NA", "NA", "NA", "NA","paired p-value (intensity)", "one-side p-value (raw Fc)", "KeGGdata_0hr10hr_KeGGIntPairedRawFcOneside_PvalComp_DensityPlot.pdf")
makedensitycoloredscatterplot(KeGGdataMerge$KeGG_Int_0hr10hr_pairedTtestPval,KeGGdataMerge$KeGG_Int_0hr10hr_unpairedTtestPval, "NA", "NA", "NA", "NA","paired p-value (intensity)", "unpaired p-value (intensity)", "KeGGdata_0hr10hr_KeGGIntPairedKeGGIntUnpaired_PvalComp_DensityPlot.pdf")


##pairwise comparison of replicate peptide quantifications
makedensitycoloredscatterplot(KeGGdataMerge$Intensity.0hr_Log2_MedNorm.Rd2, KeGGdataMerge$Intensity.0hr_Log2_MedNorm.Rd3, "NA", "NA", "NA", "NA", "0hr Run2", "0hr Run3", "KeGGdata_0hr_Run2Run3_IntensityComp_DensityPlot.pdf")
makedensitycoloredscatterplot(KeGGdataMerge$Intensity.0hr_Log2_MedNorm.Rd2, KeGGdataMerge$Intensity.0hr_Log2_MedNorm.Rd4, "NA", "NA", "NA", "NA", "0hr Run2", "0hr Run4", "KeGGdata_0hr_Run2Run4_IntensityComp_DensityPlot.pdf")
makedensitycoloredscatterplot(KeGGdataMerge$Intensity.0hr_Log2_MedNorm.Rd3, KeGGdataMerge$Intensity.0hr_Log2_MedNorm.Rd4, "NA", "NA", "NA", "NA", "0hr Run3", "0hr Run4", "KeGGdata_0hr_Run3Run4_IntensityComp_DensityPlot.pdf")
makedensitycoloredscatterplot(KeGGdataMerge$Intensity.6hr_Log2_MedNorm.Rd2, KeGGdataMerge$Intensity.6hr_Log2_MedNorm.Rd3, "NA", "NA", "NA", "NA", "6hr Run2", "6hr Run3", "KeGGdata_6hr_Run2Run3_IntensityComp_DensityPlot.pdf")
makedensitycoloredscatterplot(KeGGdataMerge$Intensity.6hr_Log2_MedNorm.Rd2, KeGGdataMerge$Intensity.6hr_Log2_MedNorm.Rd4, "NA", "NA", "NA", "NA", "6hr Run2", "6hr Run4", "KeGGdata_6hr_Run2Run4_IntensityComp_DensityPlot.pdf")
makedensitycoloredscatterplot(KeGGdataMerge$Intensity.6hr_Log2_MedNorm.Rd3, KeGGdataMerge$Intensity.6hr_Log2_MedNorm.Rd4, "NA", "NA", "NA", "NA", "6hr Run3", "6hr Run4", "KeGGdata_6hr_Run3Run4_IntensityComp_DensityPlot.pdf")
makedensitycoloredscatterplot(KeGGdataMerge$Intensity.8hr_Log2_MedNorm.Rd2, KeGGdataMerge$Intensity.8hr_Log2_MedNorm.Rd3, "NA", "NA", "NA", "NA", "8hr Run2", "8hr Run3", "KeGGdata_8hr_Run2Run3_IntensityComp_DensityPlot.pdf")
makedensitycoloredscatterplot(KeGGdataMerge$Intensity.8hr_Log2_MedNorm.Rd2, KeGGdataMerge$Intensity.8hr_Log2_MedNorm.Rd4, "NA", "NA", "NA", "NA", "8hr Run2", "8hr Run4", "KeGGdata_8hr_Run2Run4_IntensityComp_DensityPlot.pdf")
makedensitycoloredscatterplot(KeGGdataMerge$Intensity.8hr_Log2_MedNorm.Rd3, KeGGdataMerge$Intensity.8hr_Log2_MedNorm.Rd4, "NA", "NA", "NA", "NA", "8hr Run3", "8hr Run4", "KeGGdata_8hr_Run3Run4_IntensityComp_DensityPlot.pdf")
makedensitycoloredscatterplot(KeGGdataMerge$Intensity.10hr_Log2_MedNorm.Rd2, KeGGdataMerge$Intensity.10hr_Log2_MedNorm.Rd3, "NA", "NA", "NA", "NA", "10hr Run2", "10hr Run3", "KeGGdata_10hr_Run2Run3_IntensityComp_DensityPlot.pdf")
makedensitycoloredscatterplot(KeGGdataMerge$Intensity.10hr_Log2_MedNorm.Rd2, KeGGdataMerge$Intensity.10hr_Log2_MedNorm.Rd4, "NA", "NA", "NA", "NA", "10hr Run2", "10hr Run4", "KeGGdata_10hr_Run2Run4_IntensityComp_DensityPlot.pdf")
makedensitycoloredscatterplot(KeGGdataMerge$Intensity.10hr_Log2_MedNorm.Rd3, KeGGdataMerge$Intensity.10hr_Log2_MedNorm.Rd4, "NA", "NA", "NA", "NA", "10hr Run3", "10hr Run4", "KeGGdata_10hr_Run3Run4_IntensityComp_DensityPlot.pdf")
##end pairwise comparison of replicate peptide quantifications


if ( FALSE )
{
#wcpdecrease <- WCP_AllDec_Data[,c("iBAQ_0hr6hr_Log2MedNorm_Avg_Fc", "iBAQ_0hr8hr_Log2MedNorm_Avg_Fc", "iBAQ_0hr10hr_Log2MedNorm_Avg_Fc")]
#wcpdecrease <- cbind(a = 0, wcpdecrease)
wcpdecrease <- WCP_AllDec_Data[,c("iBAQ.0hr_Log2MedNorm_Zscore","iBAQ.6hr_Log2MedNorm_Zscore","iBAQ.8hr_Log2MedNorm_Zscore","iBAQ.10hr_Log2MedNorm_Zscore")]
#wcpdecrease <- WCP_10hr_Decreased_Data[,c("iBAQ.0hr_Log2MedNorm_Zscore","iBAQ.6hr_Log2MedNorm_Zscore","iBAQ.8hr_Log2MedNorm_Zscore","iBAQ.10hr_Log2MedNorm_Zscore")]
head(wcpdecrease)
filename <- "decreasedwcpprots.pdf"
pdf(filename)
par(bty="l")
ymin <- floor(min(wcpdecrease,na.rm=TRUE))
ymax <- ceiling(max(wcpdecrease,na.rm=TRUE))
xmin <- 1
xmax <- 4
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
for ( i in 1:length(row.names(wcpdecrease)) ){
#for ( i in 1:5 ){
	currdata<-t(wcpdecrease[i,])
	lines(currdata, col="grey50", lwd=1, xlab="", ylab="", xaxt='n', yaxt='n')
}
axis(1, c(seq(1,4,by=1)), labels=c("0", "6", "8", "10"), col="black", cex.axis=0.75, tck=-0.01)
mtext("transduction time course",1,line=1,cex=0.5)
axis(2, c(seq(ymin,ymax,by=1)), col="black", cex.axis=0.75, tck=-0.01)
mtext("abundance z-score",2,line=1,cex=0.5)
mtext("fold change",3,line=0.25,cex=0.5)
}

##plot results for MRE11
#average intensities
MRE11data<-subset(KeGGdataMerge, KeGGdataMerge$Gene.names=="MRE11A")
MRE11dataTmp<-MRE11data[,c("Intensity.0hr_Log2MedNorm_Avg", "Intensity.6hr_Log2MedNorm_Avg", "Intensity.8hr_Log2MedNorm_Avg", "Intensity.10hr_Log2MedNorm_Avg")]
filename <- "MRE11data_AvgIntensity.pdf"
pdf(filename, width=6, height=5)
par(bty="l")
xmin<-1
xmax<-4
ymin<-floor(min(MRE11dataTmp, na.rm=TRUE))
ymax<-ceiling(max(MRE11dataTmp, na.rm=TRUE))
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
peptides<-row.names(MRE11data)
cols<-c("red", "red", "grey50", "grey50", "grey50", "red", "red", "red", "grey50")
#cols<-c("red", "orange", "yellow", "green", "blue", "purple", "grey70", "black", "grey30")
for ( i in 1:length(row.names(MRE11dataTmp)) ){
	peptide<-peptides[i]
	abundancedata<-c(MRE11dataTmp[peptide,"Intensity.0hr_Log2MedNorm_Avg"], MRE11dataTmp[peptide,"Intensity.6hr_Log2MedNorm_Avg"], MRE11dataTmp[peptide,"Intensity.8hr_Log2MedNorm_Avg"], MRE11dataTmp[peptide,"Intensity.10hr_Log2MedNorm_Avg"])
	par(new=T)
	colid<-cols[i]
	lines(abundancedata, col=colid, lwd=4, xlab="", ylab="", xaxt='n', yaxt='n')
}
legend("topleft", col=cols, lwd=1, legend=peptides, cex = 0.5, bty = "n")
axis(1, c(seq(xmin,xmax,by=1)), labels=c("0hr", "6hr", "8hr", "10hr"), col="black", cex.axis=0.9, tck=-0.01)
mtext("transduction time course",1,line=2.5,cex=1)
axis(2, c(seq(ymin,ymax,by=1)), col="black", cex.axis=0.9, tck=-0.01)
mtext("normalized intensity",2,line=2,cex=1)
dev.off()
#average intensity zscores
MRE11dataTmp<-MRE11data[,c("Intensity.0hr_Log2MedNorm_Zscore", "Intensity.6hr_Log2MedNorm_Zscore", "Intensity.8hr_Log2MedNorm_Zscore", "Intensity.10hr_Log2MedNorm_Zscore")]
filename <- "MRE11data_IntensityZscore.pdf"
pdf(filename, width=6, height=5)
par(bty="l")
xmin<-1
xmax<-4
ymin<-floor(min(MRE11dataTmp, na.rm=TRUE))
ymax<-ceiling(max(MRE11dataTmp, na.rm=TRUE))
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
peptides<-row.names(MRE11dataTmp)
cols<-c("red", "red", "grey50", "grey50", "grey50", "red", "red", "red", "grey50")
for ( i in 1:length(row.names(MRE11dataTmp)) ){
	peptide<-peptides[i]
	abundancedata<-c(MRE11dataTmp[peptide,"Intensity.0hr_Log2MedNorm_Zscore"], MRE11dataTmp[peptide,"Intensity.6hr_Log2MedNorm_Zscore"], MRE11dataTmp[peptide,"Intensity.8hr_Log2MedNorm_Zscore"], MRE11dataTmp[peptide,"Intensity.10hr_Log2MedNorm_Zscore"])
	par(new=T)
	colid<-cols[i]
	lines(abundancedata, col=colid, lwd=4, xlab="", ylab="", xaxt='n', yaxt='n')
}
legend("topleft", col=cols, lwd=1, legend=peptides, cex = 0.5, bty = "n")
axis(1, c(seq(xmin,xmax,by=1)), labels=c("0hr", "6hr", "8hr", "10hr"), col="black", cex.axis=0.9, tck=-0.01)
mtext("transduction time course",1,line=2.5,cex=1)
axis(2, c(seq(ymin,ymax,by=1)), col="black", cex.axis=0.9, tck=-0.01)
mtext("intensity zscore",2,line=2,cex=1)
dev.off()


##plot results for RAD50
#average intensity zscores
RAD50data<-subset(KeGGdataMerge, KeGGdataMerge$Gene.names=="RAD50")
RAD50dataTmp<-RAD50data[,c("Intensity.0hr_Log2MedNorm_Zscore", "Intensity.6hr_Log2MedNorm_Zscore", "Intensity.8hr_Log2MedNorm_Zscore", "Intensity.10hr_Log2MedNorm_Zscore")]
filename <- "RAD50data_IntensityZscore.pdf"
pdf(filename, width=6, height=5)
par(bty="l")
xmin<-1
xmax<-4
ymin<-floor(min(RAD50dataTmp, na.rm=TRUE))
ymax<-ceiling(max(RAD50dataTmp, na.rm=TRUE))
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
peptides<-row.names(RAD50dataTmp)
cols<-c("red", "orange", "yellow", "green", "blue", "purple", "grey50", "gold", "salmon")
for ( i in 1:length(row.names(RAD50dataTmp)) ){
	peptide<-peptides[i]
	abundancedata<-c(RAD50dataTmp[peptide,"Intensity.0hr_Log2MedNorm_Zscore"], RAD50data[peptide,"Intensity.6hr_Log2MedNorm_Zscore"], RAD50data[peptide,"Intensity.8hr_Log2MedNorm_Zscore"], RAD50data[peptide,"Intensity.10hr_Log2MedNorm_Zscore"])
	par(new=T)
	colid<-cols[i]
	lines(abundancedata, col=colid, lwd=4, xlab="", ylab="", xaxt='n', yaxt='n')
}
legend("bottomright", col=cols, lwd=1, legend=peptides, cex = 0.5, bty = "n")
axis(1, c(seq(xmin,xmax,by=1)), labels=c("0hr", "6hr", "8hr", "10hr"), col="black", cex.axis=0.9, tck=-0.01)
mtext("transduction time course",1,line=2.5,cex=1)
axis(2, c(seq(ymin,ymax,by=1)), col="black", cex.axis=0.9, tck=-0.01)




##plot results for MRE11 and RAD50
##matrix of intensities
#all replicate intensities
#####KeGGdata_RepIntensities<-KeGGdataMerge6hr8hr10hrPairedSigInc[,c("Intensity.0hr_Log2_MedNorm.Rd2","Intensity.0hr_Log2_MedNorm.Rd3","Intensity.0hr_Log2_MedNorm.Rd4","Intensity.6hr_Log2_MedNorm.Rd2","Intensity.6hr_Log2_MedNorm.Rd3","Intensity.6hr_Log2_MedNorm.Rd4","Intensity.8hr_Log2_MedNorm.Rd2","Intensity.8hr_Log2_MedNorm.Rd3","Intensity.8hr_Log2_MedNorm.Rd4","Intensity.10hr_Log2_MedNorm.Rd2","Intensity.10hr_Log2_MedNorm.Rd3","Intensity.10hr_Log2_MedNorm.Rd4")]
#####KeGGdata_RepIntensities<-KeGGdataMerge6hr8hr10hrUnpairedSigInc[,c("Intensity.0hr_Log2_MedNorm.Rd2","Intensity.0hr_Log2_MedNorm.Rd3","Intensity.0hr_Log2_MedNorm.Rd4","Intensity.6hr_Log2_MedNorm.Rd2","Intensity.6hr_Log2_MedNorm.Rd3","Intensity.6hr_Log2_MedNorm.Rd4","Intensity.8hr_Log2_MedNorm.Rd2","Intensity.8hr_Log2_MedNorm.Rd3","Intensity.8hr_Log2_MedNorm.Rd4","Intensity.10hr_Log2_MedNorm.Rd2","Intensity.10hr_Log2_MedNorm.Rd3","Intensity.10hr_Log2_MedNorm.Rd4")]
#####KeGGdata_RepIntensities<-KeGGdataMerge6hr8hr10hrUnique[,c("Intensity.0hr_Log2_MedNorm.Rd2","Intensity.0hr_Log2_MedNorm.Rd3","Intensity.0hr_Log2_MedNorm.Rd4","Intensity.6hr_Log2_MedNorm.Rd2","Intensity.6hr_Log2_MedNorm.Rd3","Intensity.6hr_Log2_MedNorm.Rd4","Intensity.8hr_Log2_MedNorm.Rd2","Intensity.8hr_Log2_MedNorm.Rd3","Intensity.8hr_Log2_MedNorm.Rd4","Intensity.10hr_Log2_MedNorm.Rd2","Intensity.10hr_Log2_MedNorm.Rd3","Intensity.10hr_Log2_MedNorm.Rd4")]
#####KeGGdata_RepIntensities<-KeGGdataMerge6hr8hr10hr2foldInc3reps[,c("Intensity.0hr_Log2_MedNorm.Rd2","Intensity.0hr_Log2_MedNorm.Rd3","Intensity.0hr_Log2_MedNorm.Rd4","Intensity.6hr_Log2_MedNorm.Rd2","Intensity.6hr_Log2_MedNorm.Rd3","Intensity.6hr_Log2_MedNorm.Rd4","Intensity.8hr_Log2_MedNorm.Rd2","Intensity.8hr_Log2_MedNorm.Rd3","Intensity.8hr_Log2_MedNorm.Rd4","Intensity.10hr_Log2_MedNorm.Rd2","Intensity.10hr_Log2_MedNorm.Rd3","Intensity.10hr_Log2_MedNorm.Rd4")]
#####KeGGdata_RepIntensities <- KeGGdataMerge[,c("Intensity.0hr_Log2_MedNorm.Rd2","Intensity.0hr_Log2_MedNorm.Rd3","Intensity.0hr_Log2_MedNorm.Rd4","Intensity.6hr_Log2_MedNorm.Rd2","Intensity.6hr_Log2_MedNorm.Rd3","Intensity.6hr_Log2_MedNorm.Rd4","Intensity.8hr_Log2_MedNorm.Rd2","Intensity.8hr_Log2_MedNorm.Rd3","Intensity.8hr_Log2_MedNorm.Rd4","Intensity.10hr_Log2_MedNorm.Rd2","Intensity.10hr_Log2_MedNorm.Rd3","Intensity.10hr_Log2_MedNorm.Rd4")]
#
#avg replicate intensities
#####KeGGdata_RepIntensities<-KeGGdataMerge6hr8hr10hrPairedSigInc[,c("Intensity.0hr_Log2MedNorm_Avg", "Intensity.6hr_Log2MedNorm_Avg", "Intensity.8hr_Log2MedNorm_Avg", "Intensity.10hr_Log2MedNorm_Avg")]
#####KeGGdata_RepIntensities<-KeGGdataMerge6hr8hr10hrUnpairedSigInc[,c("Intensity.0hr_Log2MedNorm_Avg", "Intensity.6hr_Log2MedNorm_Avg", "Intensity.8hr_Log2MedNorm_Avg", "Intensity.10hr_Log2MedNorm_Avg")]
#####KeGGdata_RepIntensities<-KeGGdataMerge6hr8hr10hr2foldInc3reps[,c("Intensity.0hr_Log2MedNorm_Avg", "Intensity.6hr_Log2MedNorm_Avg", "Intensity.8hr_Log2MedNorm_Avg", "Intensity.10hr_Log2MedNorm_Avg")]
#####KeGGdata_RepIntensities<-KeGGdataMerge6hr8hr10hrUnique[,c("Intensity.0hr_Log2MedNorm_Avg", "Intensity.6hr_Log2MedNorm_Avg", "Intensity.8hr_Log2MedNorm_Avg", "Intensity.10hr_Log2MedNorm_Avg")]




##read RNA binding proteins from reactome analysis
#RNAbindingProteins <- scan("./KeggProteinInteractionAnalysis/RNABindingProts_PolyARNABindingProts_KeGGIncProteinInteraction_Network.txt", character(), quote = "")
#RNAbindingProteins_ProtId <- scan("./KeggProteinInteractionAnalysis/RNABindingProts_PolyARNABindingProts_KeGGIncProteinInteraction_ProtId_Network.txt", character(), quote = "")

##read RALY/hnRNPC module proteins from reactome analysis
RNAbindingProteins_ProtId <- scan("./KeggProteinInteractionAnalysis/KeGGdata_0hr6hr8hr10hr_RALYhnRNPCmodule_KeGGIncProteinInteraction_CytoscapeReactomeModules_ProtId.txt", character(), quote = "")

##########ProtData <- KeGGdataMerge_ProteinData[KeGGdataMerge_ProteinData$GeneID%in%RNAbindingProteins,]
ProtData <- KeGGdataMerge_ProteinData[KeGGdataMerge_ProteinData$Protein%in%RNAbindingProteins_ProtId,]
ProtData <- ProtData[order(ProtData$intensityweightedmeanratio, decreasing=TRUE),]

#order the protein data matrix by specified rows
#ProtData <- ProtData[c("Q9UKM9", "O75390", "P38919", "Q01844", "Q14527", "Q7KZF4", "O43148", "P26639", "Q00839", "Q14204", "P61313", "P62424", "P07910", "P78527", "Q8NC51", "P27824", "P04792", "Q12905", "P08865", "Q86U42", "P62241", "P04150", "P39687", "P05204", "P62244", "P63165"),]
ProtData <- ProtData[c("Q9UKM9", "Q7KZF4", "P38919", "Q00839", "P04844", "P61313", "P62424", "P07910", "P28288", "P33897", "P08865", "Q9Y450", "Q86U42", "P67936", "P62241", "P07951", "P62244"),]

#ProtData <- subset(ProtData, ProtData$KeGG_Prot_0hr10hr_Fc_Norm>1)

write.table(ProtData, file="ProtData_0hr6hr8hr10hr_RALYhnRNPCmoduleProteins_CytoscapeNetwork_Data_Unfiltered.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)


#ProtData <- subset(KeGGdataMerge, KeGGdataMerge$Gene.names=="MRE11A")
#ProtData <- subset(KeGGdataMerge, KeGGdataMerge$Gene.names=="RAD50")
#ProtData <- subset(KeGGdataMerge, grepl("P49959", row.names(KeGGdataMerge)))
#ProtData <- subset(KeGG_6hr8hr10hr_Inc_Union_Data, grepl("P49959", row.names(KeGG_6hr8hr10hr_Inc_Union_Data)))
#ProtData <- subset(WCPdata, (WCPdata$Gene.names.x=="MRE11A" | WCPdata$Gene.names.x=="RAD50"))
#dim(ProtData)
#KeGGdata_RepIntensities <- ProtData[,c("Intensity.0hr_Log2_MedNorm.Rd2", "Intensity.0hr_Log2_MedNorm.Rd3", "Intensity.0hr_Log2_MedNorm.Rd4", "Intensity.6hr_Log2_MedNorm.Rd2", "Intensity.6hr_Log2_MedNorm.Rd3", "Intensity.6hr_Log2_MedNorm.Rd4", "Intensity.8hr_Log2_MedNorm.Rd2", "Intensity.8hr_Log2_MedNorm.Rd3", "Intensity.8hr_Log2_MedNorm.Rd4", "Intensity.10hr_Log2_MedNorm.Rd2", "Intensity.10hr_Log2_MedNorm.Rd3", "Intensity.10hr_Log2_MedNorm.Rd4")]
#KeGGdata_RepIntensities <- ProtData[,c("KeGG_Int_0hr6hr_FcNorm_Avg", "KeGG_Int_0hr8hr_FcNorm_Avg", "KeGG_Int_0hr10hr_FcNorm_Avg")]
#KeGGdata_RepIntensities <- ProtData[,c("Intensity.0hr_Log2MedNorm_Avg", "Intensity.6hr_Log2MedNorm_Avg", "Intensity.8hr_Log2MedNorm_Avg", "Intensity.10hr_Log2MedNorm_Avg")]
#KeGGdata_RepIntensities <- ProtData[,c("Intensity.0hr_Log2MedNorm_Zscore", "Intensity.6hr_Log2MedNorm_Zscore", "Intensity.8hr_Log2MedNorm_Zscore", "Intensity.10hr_Log2MedNorm_Zscore")]
#KeGGdata_RepIntensities <- ProtData[,c("iBAQ_0hr6hr_Log2MedNorm_Avg_Fc", "iBAQ_0hr8hr_Log2MedNorm_Avg_Fc", "iBAQ_0hr10hr_Log2MedNorm_Avg_Fc")]
#KeGGdata_RepIntensities <- ProtData[,c("maxinten0hr","maxinten10hr")]

#plot bar chart of number of peptides identified in unfiltered and filtered replicates
pdf("RALYhnRNPCmoduleProteins_WCPFc_barchart_CytoscapeNetwork_0211.pdf", height=4, width=7)
#tiff("RNAbindingProteins_KeGGProtFc_barchart_CytoscapeNetwork_0211.tiff", width = 7, height = 4, units = 'in', res = 300, compression = 'none')
#pdf("RALYhnRNPCmoduleProteins_KeGGProtFc_barchart_CytoscapeNetwork_0211.pdf", height=4, width=7)
#tiff("RALYhnRNPCmoduleProteins_KeGGProtFc_barchart_CytoscapeNetwork_0211.tiff", width = 7, height = 4, units = 'in', res = 300, compression = 'none')
par(mar=c(5.1, 5.1, 4.1, 1.1))
#RNAbindingFoldChanges<-ProtData$intensityweightedmeanratio
RNAbindingFoldChanges<-ProtData$WCP_iBAQ_0hr10hr_Fc_Avg
RNAbindingGenes<-ProtData$GeneID
barplot(RNAbindingFoldChanges, horiz=FALSE, beside=T, xlab="", ylab="", ylim=c(-1,1), lwd=2, cex.axis=1, cex.names=0.9, las=2, col=c("grey70"), names.arg=RNAbindingGenes)
abline(h=0, lwd=2)
abline(h = -0.6092695, lwd=1, lty=2, col="grey")
abline(h = 0.6842065, lwd=1, lty=2, col="grey")
box(bty="n", lwd=2)
#mtext("", side=2, line=2, cex=1)
dev.off()
##end Ub fold change number bar graph



#for ( i in 1:length(RNAbindingProteins) )
#{

#NAbindingProt <- RNAbindingProteins[i]

#KeGGdata_RepIntensities<-KeGGdataMerge[KeGGdataMerge$Gene.names%in%RNAbindingProt,]
#####KeGGdata_RepIntensities<-KeGGdataMerge[KeGGdataMerge$Gene.names.first%in%RNAbindingProteins,]
KeGGdata_RepIntensities<-KeGGdataMerge[KeGGdataMerge$Protein%in%RNAbindingProteins_ProtId,]
write.table(KeGGdata_RepIntensities, file="KeGGdata_0hr6hr8hr10hr_RALYhnRNPCmoduleProteins_IntensitiesZscores_heatmap_CytoscapeNetwork_Data_Unfiltered.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
KeGGdata_RepIntensities<-subset(KeGGdata_RepIntensities, ((KeGGdata_RepIntensities$KeGGInt.10hr_NumReps>=2 & KeGGdata_RepIntensities$KeGGInt.0hr_NumReps>=2) | (KeGGdata_RepIntensities$KeGGInt.10hr_NumReps>=2 & KeGGdata_RepIntensities$KeGGInt.0hr_NumReps==0) | (KeGGdata_RepIntensities$KeGGInt.10hr_NumReps==0 & KeGGdata_RepIntensities$KeGGInt.0hr_NumReps>=2)))
stop("NN")

#use this to filter...
#KeGGdata_RepIntensities<-subset(KeGGdata_RepIntensities, ((KeGGdata_RepIntensities$Intensity.0hr_Log2MedNorm_Zscore> -0.5 | 
#                                                          KeGGdata_RepIntensities$Intensity.6hr_Log2MedNorm_Zscore> -0.5 |
#														  KeGGdata_RepIntensities$Intensity.8hr_Log2MedNorm_Zscore> -0.5 | 
#														  KeGGdata_RepIntensities$Intensity.10hr_Log2MedNorm_Zscore> -0.5) & 
#														  ((is.na(KeGGdata_RepIntensities$KeGG_Int_0hr6hr_Fc) | KeGGdata_RepIntensities$KeGG_Int_0hr6hr_Fc>1) | 
#														   (is.na(KeGGdata_RepIntensities$KeGG_Int_0hr8hr_Fc) | KeGGdata_RepIntensities$KeGG_Int_0hr8hr_Fc>1) |
#														   (is.na(KeGGdata_RepIntensities$KeGG_Int_0hr10hr_Fc) | KeGGdata_RepIntensities$KeGG_Int_0hr10hr_Fc>1))))
#remove<-c("SUMO1_25","PRKDC_810")
#KeGGdata_RepIntensities<-KeGGdata_RepIntensities[!(KeGGdata_RepIntensities$Gene.KeGGPos%in%remove),]






#####KeGGdata_RepIntensities<-KeGGdataMerge[KeGGdataMerge$Protein%in%UbIncWcpDecProts,]
#KeGGdata_RepIntensities<-KeGGdata_RepIntensities[,c("Intensity.0hr_Log2MedNorm_Avg", "Intensity.6hr_Log2MedNorm_Avg", "Intensity.8hr_Log2MedNorm_Avg", "Intensity.10hr_Log2MedNorm_Avg")]
row.names(KeGGdata_RepIntensities)<-KeGGdata_RepIntensities$Gene.KeGGPos


KeGGdata_RepIntensities<-KeGGdata_RepIntensities[,c("Intensity.0hr_Log2MedNorm_Zscore", "Intensity.6hr_Log2MedNorm_Zscore", "Intensity.8hr_Log2MedNorm_Zscore", "Intensity.10hr_Log2MedNorm_Zscore")]
#KeGGdata_RepIntensities<-KeGGdata_RepIntensities[,c("KeGG_Int_0hr6hr_FcNorm_Avg", "KeGG_Int_0hr8hr_FcNorm_Avg", "KeGG_Int_0hr10hr_FcNorm_Avg")]
#head(KeGGdata_RepIntensities)


KeGGdata_RepIntensities[KeGGdata_RepIntensities=="NaN"] <- NA
narows <- apply(KeGGdata_RepIntensities, 1, function(x) all(is.na(x)))
KeGGdata_RepIntensities<- KeGGdata_RepIntensities[ !narows, ]

#row.names(KeGGdata_RepIntensities)<-ProtData$Gene.names.x
#head(KeGGdata_RepIntensities)
if ( dim(KeGGdata_RepIntensities)[1]==1 ){
	KeGGdata_RepIntensities["XX",]<-c(0,0,0,0)
}
KeGGdataRepIntensitiesMat <- as.matrix(KeGGdata_RepIntensities)
KeGGdataRepIntensitiesMat <- apply(KeGGdataRepIntensitiesMat,2,as.numeric)
KeGGdataRepIntensitiesMat <- round(KeGGdataRepIntensitiesMat, 2)

row.names(KeGGdataRepIntensitiesMat) <- row.names(KeGGdata_RepIntensities)

#KeGGdataRepIntensitiesMat <- KeGGdataRepIntensitiesMat[order(KeGGdataRepIntensitiesMat[,"Intensity.10hr_Log2MedNorm_Zscore"], decreasing = TRUE),]

#head(KeGGdataRepIntensitiesMat)

#row.names(KeGGdataRepIntensitiesMat) <- row.names(KeGGdata_RepIntensities)


#
##heatmap of avg intensities
print("heatmap")
#heatmapfile<-paste(RNAbindingProt, '_0hr6hr8hr10hr_RNAbindingProteins_IntensitiesZscores_heatmap.pdf', sep='')
#pdf(heatmapfile)
#emf(file = "KeGGdata_0hr6hr8hr10hr_RNAbindingProteins_IntensitiesZscores_heatmap_CytoscapeNetworkModule_filtered.emf", width = 7, height = 7, bg = "transparent", fg = "black", pointsize = 12, family = "Helvetica", coordDPI = 300, custom.lty=emfPlus, emfPlus=TRUE, emfPlusFont = FALSE, emfPlusRaster = FALSE)
pdf("KeGGdata_0hr6hr8hr10hr_RALYhnRNPCmoduleProteins_IntensitiesZscores_heatmap_CytoscapeNetworkModule_unfiltered_0211.pdf")
par(family="ArialMT")
#tiff("KeGGdata_0hr6hr8hr10hr_RNAbindingProteins_IntensitiesZscores_heatmap_CytoscapeNetwork.tiff", width = 8.5, height = 11, units = 'in', res = 300, compression = 'none')
#win.metafile("KeGGdata_0hr6hr8hr10hr_RNAbindingProteins_IntensitiesZscores_heatmap_CytoscapeNetwork_test.wmf")
col_palette <- colorRampPalette(c("yellow", "green", "blue"))(n = 98)
#
#hmcol<-brewer.pal(9,"YlGnBu")
hmcol<-brewer.pal(9,"YlOrRd")
#hmcol<-brewer.pal(9, "RdYlBu")
#hmcol<-brewer.pal(9, "RdBu")
#hmcol<-brewer.pal(9, "Reds")
heatmap.2(KeGGdataRepIntensitiesMat,
          #cellnote=KeGGdataRepIntensitiesMat,
		  notecex=1.5,
		  notecol="grey80",      # change font color of cell labels to black
		  density.info="none",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  #dendrogram="row",    # only draw a row dendrogram
		  #Rowv=TRUE,
		  dendrogram="none",
		  Rowv="NA",
		  Colv="NA",            # turn off column clustering
		  #col=bluered(100),     #col_palette or col=bluered(100),
		  #col=col_palette,
		  #breaks=col_breaks,
		  col=hmcol,
		  na.color = "grey",
		  scale="none",
		  #labRow=CullinData$GeneId.x,
		  labCol=c("0hr", "6hr", "8hr", "10hr"),
		  margins=c(6,20),
		  cexRow=0.5,
		  cexCol=1,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1.0,
		  sepwidth=c(0.0, 0.01),  # width of the borders
		  sepcolor='white',
		  colsep=c()
		  #colsep=seq(0,ncol(KeGGdataRepIntensitiesMat),by=1),
		  #rowsep=0:nrow(KeGGdataRepIntensitiesMat)
		 )
dev.off()

#}

stop("RNAbinding")

#cmeanobject<-cmeans(KeGGdataRepIntensitiesMat, 20, iter.max=100, verbose=FALSE, dist="euclidean", method="cmeans", m=2, rate.par = NULL)
#print(cmeanobject)

if ( FALSE )
{

##heatmap of wcp z-scores
#####WcpDecrease <- WCP_10hr_Decreased_Data[,c("iBAQ.0hr_Log2MedNorm_Zscore","iBAQ.6hr_Log2MedNorm_Zscore","iBAQ.8hr_Log2MedNorm_Zscore","iBAQ.10hr_Log2MedNorm_Zscore")]
WcpDecrease <- WCP_10hr_Decreased_Data[,c("iBAQ_0hr6hr_Log2MedNorm_Avg_Fc", "iBAQ_0hr8hr_Log2MedNorm_Avg_Fc", "iBAQ_0hr10hr_Log2MedNorm_Avg_Fc")]
WcpDecreaseMat <- as.matrix(WcpDecrease)
WcpDecreaseMat <- apply(WcpDecreaseMat,2,as.numeric)
rownames(WcpDecreaseMat) <- WCP_10hr_Decreased_Data$Gene.names.x
narows <- apply(WcpDecreaseMat, 1, function(x) all(is.na(x)))
WcpDecreaseMat <- WcpDecreaseMat[ !narows, ]
head(WcpDecreaseMat)
pdf("WCPdata_10hr_FoldChange_DecreasedProteins_heatmap.pdf")
#col_palette <- colorRampPalette(c("yellow", "green", "blue"))(n = 98)
#hmcol<-brewer.pal(9,"YlGnBu")
hmcol<-brewer.pal(9,"YlOrRd")
#hmcol<-brewer.pal(9, "RdYlBu")
heatmap.2(WcpDecreaseMat,
          notecol="black",      # change font color of cell labels to black
		  density.info="none",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  dendrogram="row",    # only draw a row dendrogram
		  Rowv=TRUE,
		  #dendrogram="none",
		  #Rowv="NA",
		  Colv="NA",            # turn off column clustering
		  #col=bluered(100),     #col_palette or col=bluered(100),
		  #col=col_palette,
		  #breaks=col_breaks,
		  col=hmcol,
		  na.color = "grey",
		  scale="none",
		  labRow=row.names(WcpDecreaseMat),
		  labCol=c("0hr", "6hr", "8hr", "10hr"),
		  margins=c(6,20),
		  cexRow=0.25,
		  cexCol=1,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1.0,
		  sepwidth=c(0.0, 0.01),  # width of the borders
		  sepcolor='grey',
		  colsep=0:ncol(WcpDecreaseMat),
		  #rowsep=0:nrow(WcpDecreaseMat)
		 )
dev.off()

##IntDataNorm<-normalize(WcpDecreaseMat)
##IntDataSom<-som(WcpDecreaseMat, xdim=5, ydim=5, topol="hexa",neigh="gaussian")
##plot(IntDataSom)

dim(WcpDecreaseMat)
narows <- apply(WcpDecreaseMat, 1, function(x) any(is.na(x)))
WcpDecreaseMat <- WcpDecreaseMat[ !narows, ]
#dim(WcpDecreaseMat)
#head(WcpDecreaseMat)
#WcpDecreaseMat <- apply(WcpDecreaseMat, 2, scale)
#head(WcpDecreaseMat)
set.seed(123)
gap_stat <- clusGap(WcpDecreaseMat, FUN = kmeans, nstart = 25, K.max = 50, B = 500)
print(gap_stat, method = "firstmax")
#plot(gap_stat, xlab = "Number of clusters k")



kmc <- kmeans(WcpDecreaseMat, 10, iter.max = 100, nstart = 1, algorithm = c("Hartigan-Wong"), trace=FALSE)
WcpDecreaseMat<-cbind(WcpDecreaseMat,kmc$cluster)
print(WcpDecreaseMat)
names<-colnames(WcpDecreaseMat)
print(names)
names[length(names)] <- "cluster"
print(names)
colnames(WcpDecreaseMat)<-names
head(WcpDecreaseMat)

WcpDecreaseData <- as.data.frame(WcpDecreaseMat)

write.table(WcpDecreaseData, file="WCPProtData_AllDecWcp_FoldChange_Kmeans.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)

clusnumsid<-unique(WcpDecreaseData$cluster)
cluscolors<-c("salmon", "dodgerblue", "cadetblue", "olivedrab3", "orange", "goldenrod3", "gold", "coral2", "mediumorchid3", "plum3")
#pdf("clusteranalysis.pdf")
#par(bty="l") 
ymin <- floor(min(WcpDecreaseData,na.rm=TRUE))
#ymax <- ceiling(max(WcpDecreaseData,na.rm=TRUE))
ymax <- 6
xmin <- 1
xmax <- 4
#plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
for ( i in 1:length(clusnumsid) ){
	currclus<-clusnumsid[i]
	filename <- paste (currclus, "_wcp_foldchange_cluster.pdf", sep = "", collapse = NULL)
	pdf(filename)
	par(bty="l", lwd=2)
	plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
	tmpI <- subset(WcpDecreaseData, WcpDecreaseData$cluster==currclus)
	tmpI$cluster<-NULL
	tmpI <- cbind(a = 0, tmpI)
	for ( j in 1:length(row.names(tmpI)) ){
		currdata<-t(tmpI[j,])
		lines(currdata, 
		      #col="grey80",
			  col=cluscolors[i],
			  lwd=4, xlab="", ylab="", xaxt='n', yaxt='n')
	}
	axis(1, c(seq(1,4,by=1)), labels=c("0", "6", "8", "10"), col="black", cex.axis=1.5, tck=-0.01)
	mtext("transduction time course (hpi)",1,line=2.5,cex=2)
	axis(2, c(seq(ymin,ymax,by=1)), col="black", cex.axis=1.5, tck=-0.01)
	mtext("fold change",2,line=2.5,cex=2)
	dev.off()
}
#axis(1, c(seq(1,5,by=1)), labels=c("0", "3", "6", "9", "12"), col="black", cex.axis=0.75, tck=-0.01)
#mtext("transduction time course (hpi)",1,line=1,cex=0.5)
#axis(2, c(seq(ymin,ymax,by=1)), col="black", cex.axis=0.75, tck=-0.01)
#mtext("abundance z-score",2,line=1,cex=0.5)

}
print("distmat")
#library(dtw)
#print(KeGGdataAvgIntensitiesMat)
#KeGGdata_AvgIntensities_DistMatrix <- dist(KeGGdataRepIntensitiesMat, method="euclidean")
#print(KeGGdata_AvgIntensities_DistMatrix)
#KeGGdata_AvgIntensities_DistMatrix <- dist(testmat, method="euclidean")
#KeGGdata_AvgIntensities_DistMatrix <- daisy(KeGGdataAvgIntensitiesMat,metric = "gower")
##print(KeGGdata_AvgIntensities_DistMatrix)

# hierarchical clustering
#hc <- hclust(KeGGdata_AvgIntensities_DistMatrix, method="complete")


#KeGGdataRepIntensitiesMatUnique<-KeGGdataRepIntensitiesMat[,c(2:4)]
#print(KeGGdataRepIntensitiesMatUnique)

#set.seed(123)
#gap_stat <- clusGap(KeGGdataRepIntensitiesMat, FUN = kmeans, nstart = 25, K.max = 20, B = 500)
#print(gap_stat, method = "firstmax")
#plot(gap_stat, frame = FALSE, xlab = "Number of clusters k")


#kmc <- kmeans(KeGGdataRepIntensitiesMat, 10, iter.max = 100, nstart = 1, algorithm = c("Hartigan-Wong"), trace=FALSE)
#print(kmc)

#KeGGdataRepIntensitiesMat<-cbind(KeGGdataRepIntensitiesMat,kmc$cluster)
#print(KeGGdataRepIntensitiesMat)

#write.table(KeGGdataRepIntensitiesMat, file="KeGGdata_6hr8hr10hr_AllData_Clustered_KeGGpeptideData.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)

stop("TE")

##groups<-cutree(hc, h=10)
##print(groups)
##KeGGdataRepIntensitiesMat<-cbind(KeGGdataRepIntensitiesMat,groups)
###print(KeGGdataAvgIntensitiesMat)
##
##write.table(KeGGdataRepIntensitiesMat, file="KeGGdata_6hr8hr10hr_2IncReps_Clustered_KeGGpeptideData.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)


###cluster by self-organizing map
#IntData<-KeGGdata[,c("Intensity.0hr_Log2MedNorm_Avg", "Intensity.6hr_Log2MedNorm_Avg", "Intensity.8hr_Log2MedNorm_Avg", "Intensity.10hr_Log2MedNorm_Avg")]
IntDataNorm<-normalize(KeGGdataAvgIntensitiesMat)
IntDataSom<-som(IntDataNorm, xdim=5, ydim=5, topol="hexa",neigh="gaussian")
plot(IntDataSom)


stop("Y")

##signficance tests
#sigtestresults<-apply(KeGGdata, MARGIN=1, FUN=function(x) sigtest(as.numeric(c(x[["Intensity.WT1kgg_Log2_MedNorm"]],x[["Intensity.WT2kgg_Log2_MedNorm"]],x[["Intensity.WT3kgg_Log2_MedNorm"]])),as.numeric(c(x[["Intensity.KO1kgg_Log2_MedNorm"]],x[["Intensity.KO2kgg_Log2_MedNorm"]],x[["Intensity.KO3kgg_Log2_MedNorm"]])),'UNPAIRED'))
#KeGGdata$WtKO_TtestPval<-sigtestresults[2,]
#KeGGdata$WtKO_TtestPval_negLog10<-apply(KeGGdata, MARGIN=1, FUN=function(x) -log10(as.numeric(x[["WtKO_TtestPval"]])))
#head(KeGGdata)

###unique protein ids in dataset
#proteinIDs<-vector()
#proteinIDs<-unique(KeGGdata$ProtId, incompariables=FALSE)
#proteinIDsNum<-length(proteinIDs)
##print(proteinIDsNum)
##print(proteinIDs)


##combine peptide quantification in to protein-based quantification
#initialize data frame to store protein-level information
KeGGProteinData<-data.frame()
#dim(KeGGProteinData)
#head(KeGGProteinData)

#evaluate all proteins identified in the combined dataset
#parse each unique protein id
for ( i in 1:length(proteinIDs) )
{

	#current protein
	tempprot<-as.character(proteinIDs[i])
	#print(tempprot)
	
	#current gene
	tempgene<-as.character(KeGGdata[which(KeGGdata$ProtId==tempprot),"Gene"][1])
	#print(tempgene)
	
	#populate KeGGProteinData dataframe with current gene id
	KeGGProteinData[tempprot,"GeneId"]<-tempgene
	#head(KeGGProteinData)
	
	#peptide data for each identified KeGG site for the current protein
	tempdata<-subset(KeGGdata, KeGGdata$ProtId==tempprot)
	#print("tmpdata")
	#print(tempdata)
	
	#num of identified sites for current protein
	numsites<-dim(tempdata)[1]
	#print(numsites)
	
	#populate protein data with number of sites for current protein
	KeGGProteinData[tempprot,"numsites"]<-numsites
	#print(numsites)
	
	######max intensity peptide value from all identified sites
	#####maxinten<-ifelse(!all(is.na(tempdata[,"WtKOkgg_IntSum"])), max(tempdata[,"WtKOkgg_IntSum"], na.rm=TRUE), NA)
	######print(maxinten)
	######populate protein data frame
	#####KeGGProteinData[tempprot,"maxinten"]<-maxinten
	
	#perform ubiquitin abundance analysis for proteins with at least one site
	#skip "0" sites, these proteins were identified in a different study
	if ( numsites>0 )
	{
		
		#max intensity peptide value from all identified sites
		maxinten<-ifelse(!all(is.na(tempdata[,"WtKOkgg_IntSum"])), max(tempdata[,"WtKOkgg_IntSum"], na.rm=TRUE), NA)
		#print(maxinten)
		#populate protein data frame
		KeGGProteinData[tempprot,"maxinten"]<-maxinten
	
		#peptide that exhibits maximum intensity
		maxintenpep<-row.names(tempdata[which(tempdata[,"WtKOkgg_IntSum"]==maxinten), ])[1]
		#print(maxintenpep)
		#populate protein data frame
		KeGGProteinData[tempprot,"maxintenpep"]<-maxintenpep
		
		#max intensity peptide fold cahge
		maxintenfc<-tempdata[maxintenpep,"WtKOkgg_Fc"]
		#print(maxintenfc)
		#populate protein data frame
		KeGGProteinData[tempprot,"maxintenpepfoldchange"]<-maxintenfc
		
		#intensity weighted mean
		#calc total intensity
		totintensity<-sum(tempdata[,"WtKOkgg_IntSum"], na.rm=TRUE)
		#weight fold changes by by raw intensity
		intensityweighted<-apply(tempdata, MARGIN=1, FUN=function(x) weight(x[["WtKOkgg_Fc"]],x[["WtKOkgg_IntSum"]]))
		#populate temp data frame with weighted intensity
		tempdata[,"intensitywtfoldch"]<-intensityweighted
		#calc intensity weighted mean
		#formula: intensityweightedmean<-sum(tempdata[,"intensitywtfoldch"], na.rm=TRUE)/totintensity
		wtfcsum<-sum(tempdata[,"intensitywtfoldch"], na.rm=TRUE)
		if ( !is.na(wtfcsum) & wtfcsum!=0 & !is.na(totintensity) ){
			intensityweightedmean<-wtfcsum/totintensity
		}else{
			intensityweightedmean<-'NA'
		}
		#print(intensityweightedmean)
		#populate protein data frame
		KeGGProteinData[tempprot,"intensityweightedmeanfc"]<-intensityweightedmean
	
	}else{
	
		KeGGProteinData[tempprot,"maxinten"]<-NA
		KeGGProteinData[tempprot,"maxintenpep"]<-NA
		KeGGProteinData[tempprot,"maxintenpepfoldchange"]<-NA
		KeGGProteinData[tempprot,"intensityweightedmeanfc"]<-NA
	
	} #end numsites>0 loop

} #end protein site loop
#dim(KeGGProteinData)
#head(KeGGProteinData)

##write table of peptide Ub quantification data
write.table(KeGGdata, file="Cul4bKOproteomics_KeGGpeptideData.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
##write table of protein Ub quantification data
write.table(KeGGProteinData, file="Cul4bKOproteomics_KeGGproteinData.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)


##sites for which the peptide has a significance value
PvalPep <- row.names(subset(KeGGdata, !is.na(KeGGdata$WtKO_TtestPval)))
PvalPepNum <- length(PvalPep)
print("num id w/ pval")
print(PvalPepNum)


##sites for which KO has signficantly lower KeGG
KOdecUb <- subset(KeGGdata, (KeGGdata$WtKO_TtestPval<0.05 & KeGGdata$WtKOkgg_Fc <= -0.5849625))
dim(KOdecUb)
KOdecUbProtId<-unique(KOdecUb$ProtId, incompariables=FALSE)
KOdecUbProtIdNum<-length(KOdecUbProtId)
print(KOdecUbProtIdNum)
#protein data for these prots
KeGGProteinDataDecUb<-KeGGProteinData[row.names(KeGGProteinData)%in%KOdecUbProtId,]
#write tables to files
write.table(KOdecUb, file="Cul4bKOproteomics_KeGGdecKOpeptideData.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
write.table(KeGGProteinDataDecUb, file="Cul4bKOproteomics_KeGGdecKOproteinData.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)


##sites for which KO has signficantly higher KeGG
KOincUb <- subset(KeGGdata, (KeGGdata$WtKO_TtestPval<0.05 & KeGGdata$WtKOkgg_Fc >= 0.5849625))
dim(KOincUb)
KOincUbProtId<-unique(KOincUb$ProtId, incompariables=FALSE)
KOincUbProtIdNum<-length(KOincUbProtId)
print(KOincUbProtIdNum)
#protein data for these prots
KeGGProteinDataIncUb<-KeGGProteinData[row.names(KeGGProteinData)%in%KOincUbProtId,]
#write tables to files
write.table(KOincUb, file="Cul4bKOproteomics_KeGGincKOpeptideData.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
write.table(KeGGProteinDataIncUb, file="Cul4bKOproteomics_KeGGincKOproteinData.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)

stop("STOPPED")

KOdecUbNonSig <- subset(KeGGdata, ( !is.na(KeGGdata$WtKOkgg_Fc) & ((!is.na(KeGGdata$Intensity.WT1kgg_Log2_MedNorm) & !is.na(KeGGdata$Intensity.WT2kgg_Log2_MedNorm)) |
                                                                         (!is.na(KeGGdata$Intensity.WT1kgg_Log2_MedNorm) & !is.na(KeGGdata$Intensity.WT3kgg_Log2_MedNorm)) |
																		 (!is.na(KeGGdata$Intensity.WT2kgg_Log2_MedNorm) & !is.na(KeGGdata$Intensity.WT3kgg_Log2_MedNorm))
																	   )))
dim(KOdecUbNonSig)



