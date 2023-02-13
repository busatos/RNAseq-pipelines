setwd("C:/Users/Sebastiano Busato/Desktop/RNAseq INVFA")
x <- (c("readr", "DESeq2", "tximportData", "tximport", "dplyr", "ggplot2", "hexbin", "pheatmap"))
lapply(x, library, character.only = T)
library(tidyr)
library(stringr)
library(edgeR)

#input count data
countData<- as.matrix(read.csv("gene_count_matrix.csv", row.names = "gene_id"))
#read conditions into vector
cols <- colnames(countData)
condition <- str_split(cols, "_", n = 2, simplify = T)[,2]

#generate DGEList 
d <- DGEList(counts=countData,group=factor(condition))
#filtering
dim(d)
d.full <- d # keep whole set 
head(cpm(d))
apply(d$counts, 2, sum)
keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)

d$samples$lib.size <- colSums(d$counts)
d$samples

#normalizing the data
d <- calcNormFactors(d)

#Multidimensional scaling
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("topright", as.character(unique(d$samples$group)), col=1:3, pch=20)

#estimate dispersion
d1 <- estimateCommonDisp(d, verbose=T)
d1 <- estimateTagwiseDisp(d1)
plotBCV(d1)


#GLM estimates of dispersion
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

#Differential expression, Exact Tests 
et12 <- exactTest(d1, pair=c(1,2)) # compare groups 1 and 2
et13 <- exactTest(d1, pair=c(1,3))
et14 <- exactTest(d1, pair=c(1,4))
topTags(et14, n=10)
