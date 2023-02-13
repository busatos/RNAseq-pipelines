#some needed maintenacne stuff, packages and such
setwd("C:/Users/Sebastiano Busato/Desktop/RNAseq INVFA")
x <- (c("readr", "DESeq2", "tximportData", "tximport", "dplyr", "ggplot2", "hexbin", "pheatmap"))
lapply(x, library, character.only = T)
library(tidyr)
library(stringr)
library(data.table)
library(purrr)
library(tidyverse)

#import counts 
countData<- as.matrix(read.csv("gene_count_matrix.csv", row.names = "gene_id"))
#read colnames and generate basic conditions (TRT groups)
cols <- colnames(countData)
condition <- str_split(cols, "_", n = 2, simplify = T)[,2]
#generate specific list that specify only FA type, or...
FA_type <- c("EB", "EB", "EB", "CTR", "FA",
             "FA", "CTR", "FA", "CTR", "FA", "CTR", 
             "EB", "EB", "EB", "FA", "FA", "FA",
             "CTR", "EB", "FA", "EB", "CTR", "EB",
             "FA", "FA", "FA", "FA", "EB", "EB", "EB")
#TEST -- CTR is a FA TYPE (HERE == FA), but concentration == 0 
FA_type_t2 <- c("EB", "EB", "EB", "FA", "FA",
             "FA", "FA", "FA", "FA", "FA", "FA", 
             "EB", "EB", "EB", "FA", "FA", "FA",
             "FA", "EB", "FA", "EB", "FA", "EB",
             "FA", "FA", "FA", "FA", "EB", "EB", "EB")
##generate specific list for FA amount
FA_amount <- c("150", "150", "75", "0", "75",
             "75", "0", "75", "0", "75", "0", 
             "150", "75", "75", "150", "75", "150",
             "0", "75", "150", "75", "0", "150",
             "150", "75", "150", "150", "150", "150", "75") 
##add cycle specifics 
exp_cycle <- c(1,1,1,1,2,
               1,2,1,2,2,
               1,2,2,1,2,
               1,1,1,1,2,
               2,2,2,1,2,
               2,1,1,2,2)
#generate colData matrix
colData <- data.frame(row.names = cols, condition, FA_type_t2, FA_amount, factor(exp_cycle))
#DataSet from Matrix 
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition) 

#filtering for counts
nrow(dds) #initial number of rows
keep <- rowSums(counts(dds)>4) >= 5 #filter out 0 or 1 count across at least one sample
dds <- dds[keep]
nrow(dds) #number of filter rows

#PCA on multidimensional data
rld <- rlog(dds, blind = TRUE) # generate rlog transformation (takes a while )
#head(assay(rld)) #if you want to see it 

df <- as_tibble(assay(rld)[, 1:2])
colnames(df)[1:2] <- c("x", "y")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() #not sure what this is

#generate PCA plots (change returnData to TRUE if you want the numbers instead)
plotPCA(rld, intgroup="condition", returnData=F)  ##additional arguments ntop = number of top genes 
#PCA plot using ggplot (adding sample name)
pca1 <- plotPCA(rld, intgroup="condition", returnData=TRUE)                                    
ggplot(pca1, aes(PC1, PC2, colour=condition, label = exp_cycle)) + geom_point(size=3) + 
  geom_text(aes(label=exp_cycle),hjust=0, vjust=0) + coord_fixed() 


###DEG determination
dds <- DESeq(dds)
resultsNames(dds) #names of contrasts, if needed

#generate results tables for all contrasts

FA075 <- results(dds, contrast=c("condition", "FA075", "CTR"), alpha = 0.2)#independentFiltering=TRUE, alpha=0.1, pAdjustMethod="BH", parallel=TRUE)
rownames(FA075)
FA075 <- data.frame(FA075)
FA150 <- results(dds, contrast=c("condition", "FA150", "CTR"), alpha = 0.2)
FA150 <- data.frame(FA150)
EB075 <- results(dds, contrast=c("condition", "EB075", "CTR"), alpha = 0.2)
EB075 <- data.frame(EB075)
EB150 <- results(dds, contrast=c("condition", "EB150", "CTR"), alpha = 0.2)
EB150 <- data.frame(EB150)
FAvsEB075 <- results(dds, contrast=c("condition", "FA075", "EB075"), alpha = 0.2)
FAvsEB150 <- results(dds, contrast=c("condition", "FA150", "EB150"), alpha = 0.2)
FA150vs075 <- results(dds, contrast=c("condition", "FA150", "FA075"), alpha = 0.2)
EB150vs075 <- results(dds, contrast=c("condition", "EB150", "EB075"), alpha = 0.2)

summary(FAvsEB075)
summary(FAvsEB150)


test_df <- data.frame(FA075$log2FoldChange, FA075$padj, 
                      FA150$log2FoldChange, FA150$padj,
                      EB075$log2FoldChange, EB075$padj,
                      EB150$log2FoldChange, EB150$padj)
rownames(test_df) <- rownames(FA075)
aa <- colnames(test_df)[grep("padj", colnames(test_df))]
sig_only_vsCTR <- na.omit(test_df[test_df$FA075.padj<0.2 | test_df$FA150.padj < 0.2 | test_df$EB075.padj < 0.2 | test_df$EB150.padj < 0.2,])

test_df_2 <- data.frame(FAvsEB075$log2FoldChange, FAvsEB075$padj,
                      FAvsEB150$log2FoldChange, FAvsEB150$padj,
                      FA150vs075$log2FoldChange, FA150vs075$padj,
                      EB150vs075$log2FoldChange, EB150vs075$padj)
rownames(test_df_2) <- rownames(FAvsEB075)
#aa <- colnames(test_df_2)[grep("padj", colnames(test_df_2))]
sig_only_FAs <- na.omit(test_df_2[test_df_2$FAvsEB075.padj<0.2 | test_df_2$FAvsEB150.padj < 0.2 | test_df_2$FA150vs075.padj < 0.2 | test_df_2$EB150vs075.padj < 0.2,])

write.csv(sig_only_vsCTR, file = "sig02_allFAvsCTR.csv")
write.csv(sig_only_FAs, file = "sig02_allFAdosecomp.csv")

aa[1]











###STUFF TO KEEP THAT I DON"T NEED TO RUN YET
#############3

summary(res)
res[order(res$padj),]
res <- lfcShrink(dds, contrast=c("condition", "FA075", "CTR"), res=res, type = "normal")
res[order(res$padj),]
summary(res)
#str(dds)
#resultsNames(dds2)
#res <- results(dds, contrast=c("condition","Lactonase","Untreated"))
#res2 <- results(dds, contrast=c("condition","Lactonase","Heat"), test = "Wald")
res <- results(dds)
res_post <- results(dds_post)
head(res)
summary(res)
res150 <- (results(dds, contrast = c("condition", "CTR", "FA150")))
summary(res150)
res150EB <- (results(dds, contrast = c("condition", "CTR", "EB150")))
summary(res150EB)
ttt <- data.frame(res75[order(res75$padj),])
res75 <- (results(dds, contrast = c("condition", "CTR", "FA075")))
res75EB <- (results(dds, contrast = c("condition", "CTR", "EB075")))
ttt <- data.frame(res75EB[order(res75EB$padj),])

res[order(res$padj),]
#ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
#resultsNames(ddsLRT)
#resLRT <- results(ddsLRT, contrast=c("condition","Lactonase","Heat"))

#resCONTRASTS <- results(dds, contrast=c(0,-1/2,1,-1/2))
#resCONTROL <- results(dds2)
resCONTRAST <- results(dds_pre, contrast = c("condition","CTR","G"))
#rds <- results(dds, contrast = list( c("Lactonase"),c("Heat", "Untreated") ) )


#head(res[order(res$padj),])   
#head(res2[order(res2$padj),]) 
hello <- data.frame(head(res_pre[order(res_pre$padj),],10000))
write.csv(hello, file = "overall_padj.csv")  
hello2 <- data.frame(head(resCONTRAST[order(resCONTRAST$padj),],10000))
write.csv(hello2, file = "CTRvsG.csv")  

#head(resLRT[order(resLRT$padj),],8) 
#head(resCONTRASTS[order(resCONTRASTS$padj),]) 
hello <- data.frame(head(res_post[order(res_post$padj),],100))
write.csv(hello, file = "post_SEvsCTR.csv")  

Pmy.file <- data.frame(resCONTROL[order(resCONTROL$padj),])

#summary(res)
#summary(res2)
summary(res)
#summary(resLRT)
#summary(resCONTRASTS, alpha = 0.05)
summary(res_post)


###heatmapv2__sorted by largest difference in 

pre_top <- rld[head(order(res_pre$padj),27),]
post_top <- rld[head(order(res_post$padj),40),]
#CONTRASTS_top <- rld[head(order(resCONTRASTS$padj),12),]
HeatvsUntreated_top <- rld[head(order(res3$padj),6),]
#rld_top <- rld[head(order(resLRT$padj),20),]
topVarGenes <- head(order(rowVars(assay(pre_top)), decreasing = T), 27)
topVarGenes <- head(order(rowVars(assay(post_top)), decreasing = T),40)
mat <- assay(post_top)[topVarGenes,]
mat <- mat - rowMeans(mat)
#library(pheatmap)
pheatmap(mat, annotation_col = colData_post, 
         cluster_rows = F, 
         cluster_cols = F)


##### 




mcols(res, use.names = T)
summary(res)

plotMA( res, ylim = c(-5, 5))
idx <- identify(res$baseMean, res$log2FoldChange)
idx

rownames(res)[idx]

sigGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = sigGene)

topVarGenes <- head(order(rowVars(assay(rld)), decreasing = T), 20)

# install.packages("pheatmap")
# BiocManager::install("genefilter", version = "3.8")

mat <- assay(rld)[topVarGenes, ]
mat <- mat - rowMeans(mat)
library(pheatmap)
pheatmap(mat, annotation_col = colData, cluster_rows = F, cluster_cols = F)


# BiocManager::install("AnnotationDbi", version = "3.8")
# BiocManager::install("org.Hs.eg.db", version = "3.8")

library("AnnotationDbi") 
library("org.Hs.eg.db")

columns(org.Hs.eg.db)

res$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$genename <-mapIds(org.Hs.eg.db, keys=row.names(res), column="GENENAME", keytype="ENSEMBL", multiVals="first")
print(res)
resOrdered <- res[order(res$pvalue),]
head(resOrdered)

resOrderedDF <-as.data.frame(resOrdered)[1:8, ]
write.csv(resOrderedDF, file = "results_reference.csv")


###IGNORE THIS
colData <- data.frame(row.names = cols, condition, FA_type_t2, FA_amount, factor(exp_cycle))
colData[condition!=]


dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition) 
