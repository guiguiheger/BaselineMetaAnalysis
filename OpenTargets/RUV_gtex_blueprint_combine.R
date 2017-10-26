#!/usr/bin/env Rscript

# Author: Suhaib Mohammed

library(RUVSeq)

##load expressioinSet objects 
load("set.RUVg.RData")
load("set.RUVg.bp.Rdata")

## loading generic functions 
source("generic_funcitons.R")

## extracting normalised counts 
norm.Gtex<-normCounts(set.RUVg)
norm.blueprint<-normCounts(set.RUVg.bp)

## meging two matrices (a,b) of different dimensions. 
## order ensemble-ids from matrix a is maintained, 
## and new ensemble-ids from b are concatenated below
## this function also transforms "NA" -> 0
## more info in generic_funcitons.R
all.norm<-mergeX(norm.Gtex, norm.blueprint)

x<-sapply(strsplit(colnames(all.norm), split="_"),"[",3)
names <- colnames(all.norm)

# experiment detailed expression values
expNames <- t(as.data.frame(sapply(names, function(x) strsplit(x, "_"))))
expNormData <- all.norm
colnames(expNormData) <- paste(expNames[,1], expNames[,3], sep="_")

# summarisation (aggregation) of expression values by median across tissue/cell types
expData <-expNormData
tissue.names<-sapply(strsplit(colnames(expData),split="_"),"[",2)
colnames(expData) <- tissue.names
expData <- t(apply(expData, 1, function(x) tapply(x, colnames(expData), median)))

png(file = paste0(name,".png"), width = 2000, height = 2000, res=180)
plot_heatmap(expData, name="normalised")
dev.off()

png(file = paste0("Distribution_Summary_NormCounts.png"), width = 700, height = 700, res=100)
hist(log(expData),freq=FALSE, col="cornflowerblue", breaks =30, xlab="log normalised counts", main = "")
dev.off()

write.table(expData, file="exp_summary_NormCounts_genes_all_Blueprint2.txt",quote=FALSE, sep="\t")
