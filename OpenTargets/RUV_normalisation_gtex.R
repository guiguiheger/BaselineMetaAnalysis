#!/usr/bin/env Rscript

library( ExpressionAtlas)

## loading generic functions 
source("generic_functions.R")

## get gtex expreriments from exression atlas
atlasData <- getAtlasData(c("E-MTAB-5214", "E-MTAB-513", "E-MTAB-2836", "E-MTAB-3716", "E-MTAB-4344"))

t<-atlasData
all <- vector()

#################################################
for (i in names(t)) {
  expAcc <- i
  k <- t[[i]]
  exp <- k$rnaseq
  eCounts <- assays(exp)$counts
  samples<-colnames(eCounts)
  average.counts<-technical_replicate_average_gtex(exp,expAcc)
  all <- cbind(all,average.counts)
}

## sanity check 
head(all)
dim(all)

x<-sapply(strsplit(colnames(all),split="_"),"[",3)


save(all,x,file="all.RData")

## filtering low expression signals
  filter <- rowSums(all>10)>=15
  filtered <- all[filter,]
  
  filterCols <- colSums(filtered == 0) / nrow(filtered) < 0.90
  x <- x[filterCols]
  filtered <- filtered[,filterCols]
  
  filtered.genes<-setdiff(rownames(all),rownames(filtered))

  write.table(filtered.genes,file="filtered_genes_GTEx.txt",sep="\t",col.names='NA')

##################################################

x<-as.factor(sapply(strsplit(colnames(filtered),split="_"),"[",3))


## coefficient of Variation
co.var <- function(x) ( 100*apply(x,1,sd)/rowMeans(x) ) 

## coefficient of variation across all the sammples
cov.allGenes<-na.omit(co.var(as.matrix(filtered)))
# Using threshold (1% qauntile) for identification of non vriable genes that has least cov.


## identyfing number of genes changed acros several cov thresholds range
cov.range<-seq(range(cov.allGenes)[1], range(cov.allGenes)[2], by = 10)
ngenes<-matrix()
for (i in seq_along(cov.range)){
  ngenes[i]<- sum((cov.allGenes<=cov.range[i])*1)
}

# Using threshold (1% qauntile) for identification of non vriable genes that has least cov.
leastVar.genes<-rownames(as.matrix(sort(cov.allGenes[cov.allGenes < quantile(cov.range, c(.01))[[1]]])))

## plot coefficent of variation against number of genes used to set for negative controls; i.e. genes can assumed not be influened 
## by the set of covariate of interest. 
png(file = paste0("covRange_threshold.all.png"), width = 750, height = 750, res=120);
plot(cov.range,ngenes, xlab="Coefficient of Variation",type="b", ylab="Number of genes",pch=16,col="blue",main="")
abline(v=quantile(cov.range, c(.01))[[1]], lwd=1, col="red")
legend("bottomright",c("1% quantile"),lty=1, lwd=1, col="red")
dev.off()


library(EDASeq)
library(ggplot2)
library(ggfortify)
library(RUVSeq)

set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(x, row.names=colnames(filtered)))

##############
colors.order<-colorOrder(x)

library(RColorBrewer)
colors <- brewer.pal(8, "Set2")
pdf("unnormalised_all_gtex.pdf", width=18, height=18)
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors.order)
plotPCA(set, col=colors.order, cex=0.7)
dev.off()

## upper quartile normalisation
set <- betweenLaneNormalization(set, which="upper")
save(set, file="upperqNorm2K.RData")
pdf("upperQnormalisation_all_gtex.pdf", width=18, height=18)
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors.order)
abline(h=c(2,-2),lty=2, col="blue",lwd=2)
plotPCA(set, col=colors.order, cex=0.7)
dev.off()

## 
# print all object size check
for ( object in ls() ){
  print(paste0(object," = ",format(object.size(get(object)), units="auto")))
}

#####
set.RUVg <- RUVg(set, leastVar.genes[1:1000] , k=1)
save(set.RUVg,file="set.RUVg.Rdata")
pdf(paste0("RUVgK1.1000.pdf"), width=18, height=18)
plotRLE(set.RUVg, outline=FALSE, ylim=c(-4, 4), col=colors.order)
abline(h=c(2,-2),lty=2, col="blue",lwd=2)
plotPCA(set.RUVg, col=colors.order, cex=0.7)
dev.off()

## normalised expression heatmap
dim(normCounts(set.RUVg))
plot_heatmap(normCounts(set.RUVg), name="Normalised") 

# raw expression
dim(counts(set.RUVg)) heatmap
plot_heatmap(counts(set.RUVg), name="Raw")