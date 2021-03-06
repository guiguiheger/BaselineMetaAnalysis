#!/usr/bin/env Rscript

# Author: Suhaib Mohammed

library( ExpressionAtlas)

## loading generic functions 
source("generic_functions.R")

## get blueprints expriments from exression atlas
atlasData <- getAtlasData(c("E-MTAB-3819", "E-MTAB-3827", "E-MTAB-4754"))

t<-atlasData
all <- vector()

for (i in names(t)) {
  expAcc <- i
  k <- t[[i]]
  exp <- k$rnaseq
  eCounts <- assays(exp)$counts
  samples<-colnames(eCounts)
  average.counts<-technical_replicate_average_bp(exp,expAcc,eCounts)
  all <- cbind(all,average.counts)
}

## sanity check 
head(all)
dim(all)

all_blueprint<-all
save(all_blueprint,file="./output/blueprint/all_blueprint.Rda")

## filtering lower signals
filter <- rowSums(all>10)>=5
filtered <- all[filter,]
x<-as.factor(sapply(strsplit(colnames(filtered),split="_"),"[",3))
dim(filtered)
filterCols <- colSums(filtered == 0) / nrow(filtered) < 0.90
x <- x[filterCols]
filtered <- filtered[,filterCols]
filtered<-na.exclude(filtered)

filtered.genes<-setdiff(rownames(all),rownames(filtered))
write.table(filtered.genes,file="./output/blueprint/filtered_genes_blueprint.txt",sep="\t",col.names='NA')

## tissue names
x<-as.factor(sapply(strsplit(colnames(filtered),split="_"),"[",3))


## coefficient of Variation
co.var <- function(x) ( 100*apply(x,1,sd)/rowMeans(x) ) 

## coefficient of variation across all the sammples
cov.allGenes<-na.omit(co.var(as.matrix(filtered)))

# Using threshold (1% qauntile) for identification of non vriable genes that has least cov.
cov.range<-seq(range(cov.allGenes)[1], range(cov.allGenes)[2], by = 10)

## identyfing number of genes changed acros several cov thresholds range
ngenes<-matrix()
for (i in seq_along(cov.range)){
  ngenes[i]<- sum((cov.allGenes<=cov.range[i])*1)
}

# Using threshold (1% qauntile) for identification of non vriable genes that has least cov.
leastVar.genes<-rownames(as.matrix(sort(cov.allGenes[cov.allGenes < quantile(cov.range, c(.1))[[1]]])))


# plot coefficent of variation against number of genes used to set for negative controls; i.e. genes can assumed not be influened 
## by the set of covariate of interest. 
png(file = paste0("./output/blueprint/covRange_threshold.all.png"), width = 750, height = 750, res=120);
plot(cov.range,ngenes, xlab="Coefficient of Variation",type="b", ylab="Number of genes",pch=16,col="blue",main="")
#abline(v=summary(cov.range)[[2]],lwd=1,col="red")
abline(v=quantile(cov.range, c(.01))[[1]], lwd=1, col="red")
legend("bottomright",c("1% quantile"),lty=1, lwd=1, col="red")
dev.off()

##
library(EDASeq)
library(EDASeq)
library(ggplot2)
library(ggfortify)
library(RUVSeq)

## making new expression set object for batch effect removal
set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(x, row.names=colnames(filtered)))

## color similar tissue types for consistency 
colors.order<-colorOrder(x)

library(RColorBrewer)
colors <- brewer.pal(8, "Set2")
pdf("./output/blueprint/unnormalised_allBluePrint.pdf", width=18, height=18)
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors.order)
plotPCA(set, col=colors.order, cex=0.7)
dev.off()

set <- betweenLaneNormalization(set, which="upper")
save(set, file="./output/blueprint/upperqNorm2K.RData")
pdf("./output/blueprint/upperQnormalisation_allBluePrint.pdf", width=18, height=18)
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors.order)
abline(h=c(2,-2),lty=2, col="blue",lwd=2)
plotPCA(set, col=colors.order, cex=0.7)
dev.off()


# print all object for size check
for ( object in ls() ){
  print(paste0(object," = ",format(object.size(get(object)), units="auto")))
}

#####
set.RUVg.bp <- RUVg(set, leastVar.genes[1:1000] , k=1)
save(set.RUVg.bp,file="./output/blueprint/set.RUVg_bp.Rdata")
pdf(paste0("./output/blueprint/RUVgK1.1000_blueprint.pdf"), width=18, height=18)
plotRLE(set.RUVg.bp, outline=FALSE, ylim=c(-4, 4), col=colors.order)
abline(h=c(2,-2),lty=2, col="blue",lwd=2)
plotPCA(set.RUVg.bp, col=colors.order, cex=0.7)
dev.off()
#####

## normalised expression
dim(normCounts(set.RUVg.bp))
agg_matrix_norm<-summary_tissues(normCounts(set.RUVg.bp))
plot_heatmap(agg_matrix_norm, name="Normalised_blueprint")

# raw expression
dim(counts(set.RUVg.bp))
agg_matrix_raw<-summary_tissues(counts(set.RUVg.bp))
plot_heatmap(agg_matrix_raw, name="Raw_blueprint")