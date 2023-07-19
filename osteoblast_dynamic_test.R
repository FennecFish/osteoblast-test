library(devtools)
if (packageVersion("fishpond") < "2.1.36") {
    install_github("mikelove/fishpond")
}

library(dplyr)
suppressPackageStartupMessages(library(SummarizedExperiment))
library(tximeta)
library(fishpond)
library(AnnotationHub)
library(ensembldb)
library(Gviz)
library(UpSetR)
library(pheatmap)

set.seed(1)

#######################################################################################
#################################### Upset Plot ######################################
#######################################################################################

load("CASTxB6.gene.dynamic.rda")

gene_dat <- as.data.frame(mcols(dy))
# select genes with significant AI at 5% FDR
sig_gene <- gene_dat %>% 
    dplyr::filter(qvalue < 0.05) %>% 
    dplyr::select(ends_with("id"), "symbol", "log2FC", "qvalue")
length(unique(sig_gene$group_id)) # 57

load("CASTxB6.fuzzy.tss.dynamic.rda")
# select tss groups with significant AI at 5% FDR
tss_dat <- as.data.frame(mcols(dy))
sig_tss <- tss_dat %>% 
    dplyr::filter(qvalue < 0.05) %>%
    dplyr::select(ends_with("id"), "symbol", "log2FC", "qvalue")
length(unique(sig_tss$gene_id)) # 49

load("CASTxB6.txp.dynamic.rda")
# select txps groups with significant AI at 5% FDR
txp_dat <- as.data.frame(mcols(dy))
sig_txp <- txp_dat %>% 
    dplyr::filter(qvalue < 0.05) %>%
    dplyr::select(ends_with("id"), "symbol", "log2FC", "qvalue")
length(unique(sig_txp$gene_id)) # 23

length(unique(intersect(sig_gene$group_id,intersect(sig_txp$gene_id,sig_tss$gene_id)))) # 19
length(unique(intersect(sig_gene$group_id,sig_tss$gene_id))) #38


data <- list(Gene =unique(sig_gene$group_id), TSS =unique(sig_tss$gene_id), Txp =unique(sig_txp$gene_id))
# pdf("Dynamic_upsetPlot.pdf", height = 10, width = 15)
upset(fromList(data),
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 2, 
      point.size = 4, 
      line.size = 2,
      mainbar.y.label = "Significant Genes Distinct Intersections", 
      sets.x.label = "Significant Genes Set Size")
# dev.off()

#######################################################################################
###################################### Calcoco1 #######################################
#######################################################################################

load("CASTxB6.fuzzy.tss.dynamic.rda")
tss_dat <- mcols(dy) %>% as_tibble

load("CASTxB6.txp.dynamic.rda")
txp_dat <- mcols(dy) %>% as_tibble

# list of genes with dynamic AI and discordant TSS
list <- tss_dat %>% 
    dplyr::filter(pvalue < .01) %>%
    dplyr::group_by(gene_id) %>% 
    dplyr::summarize(count=n(), min=min(log2FC), max=max(log2FC), minQ=min(qvalue)) %>% 
    dplyr::filter(count > 1, sign(min) != sign(max)) %>% 
    as.data.frame()

# Calcoco1 is the 6th of the list
tss_dat %>% dplyr::filter(gene_id == list$gene_id[6]) %>% 
    dplyr::select(gene_id, group_id, symbol, log10mean, stat, log2FC, pvalue, qvalue, meanInfRV)
txp_dat %>% dplyr::filter(gene_id == list$gene_id[6]) %>% 
    dplyr::select(tx_id, gene_id, symbol, log10mean, stat, log2FC, pvalue, qvalue, meanInfRV, fuzzy_tss_group)

load("CASTxB6.fuzzy.tss.dynamic.rda")
levels(dy$allele) <- c("B6","CAST")
# making the plotInfReps plot
# pdf("Dynamic_infRV_Calcoco1.pdf", width = 15, height = 15)
par(mfrow=c(2,2), cex=0.3)
exp <- tss_dat %>% dplyr::filter(gene_id == list$gene_id[6]) 

plotInfReps(dy, idx = exp$group_id[1], x= "day", cov = "allele", shiftX = 0.01)
lines(dy$day[1:9], assay(dy, "mean")[exp$group_id[1],1:9], col="dodgerblue") 
lines(dy$day[10:18], assay(dy, "mean")[exp$group_id[1],10:18], col="goldenrod4")
legend("bottomright", legend=c("B6", "CAST"), pch = 10, col = c("dodgerblue", "goldenrod4"))

plotInfReps(dy, idx = exp$group_id[2], x= "day", cov = "allele", shiftX = 0.01)
lines(dy$day[1:9], assay(dy, "mean")[exp$group_id[2],1:9], col="dodgerblue") 
lines(dy$day[10:18], assay(dy, "mean")[exp$group_id[2],10:18], col="goldenrod4")

png("Dynamic_infRV_Calcoco1_sub.png", res = 300, width = 3000, height = 1500)
par(mfrow=c(1,2), cex=1)
plotInfReps(dy, idx = exp$group_id[3], x= "day", cov = "allele", shiftX = 0.01)
lines(dy$day[1:9], assay(dy, "mean")[exp$group_id[3],1:9], col="dodgerblue") 
lines(dy$day[10:18], assay(dy, "mean")[exp$group_id[3],10:18], col="goldenrod4")

plotInfReps(dy, idx = exp$group_id[4], x= "day", cov = "allele", shiftX = 0.01)
lines(dy$day[1:9], assay(dy, "mean")[exp$group_id[4],1:9], col="dodgerblue") 
lines(dy$day[10:18], assay(dy, "mean")[exp$group_id[4],10:18], col="goldenrod4")
legend("bottomright", legend=c("B6", "CAST"), pch = 10, col = c("dodgerblue", "goldenrod4"))
dev.off()

ah <- AnnotationHub()
edb <- ah[["AH89211"]]

# make plotAllelicGene plot
par(mfrow=c(1,1), cex=0.3)
dy$time_bins <- cut(dy$day,breaks=c(2,6,12,18),include.lowest=TRUE, labels=FALSE)
time_labels <- c("Day02-06","Day08-12","Day14-18")
dy$time_bins <- time_labels[ dy$time_bins ]
plotAllelicGene(dy, list$gene_id[6], edb, cov="time_bins",
                qvalue=FALSE, log2FC=FALSE,  tpmFilter = 0.1, labels=list(a2="B6",a1="CAST"))

par(mfrow=c(1,1), cex=1)
idx <- mcols(dy)$gene_id == list$gene_id[6]
colnames(dy) <- c(paste0("Day ",seq(2,18, by= 2), "-a2"), 
                  paste0("Day ",seq(2,18, by= 2), "-a1"))
row_dat <- data.frame(minusLogQ=-log10(mcols(dy)$qvalue[idx]),
                      row.names=rownames(dy)[idx])
col_dat <- data.frame(time=dy$day[1:9],
                      row.names=paste0("Day ",seq(2,18, by= 2)))
plotAllelicHeatmap(dy, idx=idx,
                   annotation_row=row_dat,
                   annotation_col=col_dat,
                   cluster_rows=FALSE, 
                   labels_row=c("TSS-1","TSS-3","TSS-5","TSS-6"), main = "")

#######################################################################################
#################################### Fig 3 Rasl11b ####################################
#######################################################################################

ah <- AnnotationHub()
edb <- ah[["AH89211"]]

load("CASTxB6.fuzzy.tss.dynamic.rda")
levels(dy$allele) <- c("B6","CAST")
rowMeans(assay(dy, "abundance")[mcols(dy)$gene_id == list$gene_id[15],])

# making infRep plot
par(mfrow=c(2,1), cex=0.5)
plotInfReps(dy, idx = "ENSMUSG00000049907-1", x= "day", cov = "allele", shiftX = 0.01)
lines(dy$day[1:9], assay(dy, "mean")["ENSMUSG00000049907-1",1:9], col="dodgerblue") 
lines(dy$day[10:18], assay(dy, "mean")["ENSMUSG00000049907-1",10:18], col="goldenrod4")
plotInfReps(dy, idx = "ENSMUSG00000049907-3", x= "day", cov = "allele", shiftX = 0.01)
lines(dy$day[1:9], assay(dy, "mean")["ENSMUSG00000049907-3",1:9], col="dodgerblue") 
lines(dy$day[10:18], assay(dy, "mean")["ENSMUSG00000049907-3",10:18], col="goldenrod4")
legend("topright", legend=c("B6", "CAST"), pch = 15, col = c("dodgerblue", "goldenrod4"))

# making allelicGene plot 
par(mfrow=c(1,1), cex=2)
dy$time_bins <- cut(dy$day,breaks=c(2,6,12,18),include.lowest=TRUE, labels=FALSE)
time_labels <- c("Day02-06","Day08-12","Day14-18")
dy$time_bins <- time_labels[ dy$time_bins ]
plotAllelicGene(dy, "ENSMUSG00000049907", edb, cov="time_bins",
                qvalue=FALSE, log2FC=FALSE,  tpmFilter = 0.1, labels=list(a2="B6",a1="CAST"))

# making heatmap plot 
par(mfrow=c(1,1), cex=0.3)
idx <- mcols(dy)$gene_id == list$gene_id[15]
colnames(dy) <- c(paste0("Day ",seq(2,18, by= 2), "-a2"), 
                  paste0("Day ",seq(2,18, by= 2), "-a1"))
row_dat <- data.frame(minusLogQ=-log10(mcols(dy)$qvalue[idx]),
                      row.names=rownames(dy)[idx])
col_dat <- data.frame(time=dy$day[1:9],
                      row.names=paste0("Day ",seq(2,18, by= 2)))
plotAllelicHeatmap(dy, idx=idx,
                   annotation_row=row_dat,
                   annotation_col=col_dat,
                   cluster_rows=FALSE, 
                   labels_row=c("TSS-1","TSS-3"), main = "")


#######################################################################################
###################### Output sig output on all three resolution ######################
#######################################################################################

load("CASTxB6.gene.dynamic.rda")

gene_dat <- as.data.frame(mcols(dy))
sig_gene <- gene_dat %>% 
    dplyr::filter(qvalue < 0.05) %>% 
    dplyr::select("symbol", "log10mean", "stat", "log2FC"
                  , "pvalue", "locfdr", "qvalue", "meanInfRV") %>%
    dplyr::arrange(pvalue)
write.csv(sig_gene, file = "dynamic_sig_gene.csv", row.names = T)

load("CASTxB6.fuzzy.tss.dynamic.rda")
tss_dat <- as.data.frame(mcols(dy))
sig_tss <- tss_dat %>% 
    dplyr::filter(qvalue < 0.05) %>%
    dplyr::select("gene_id", "symbol", "start_position", "log10mean", 
                  "stat", "log2FC", "pvalue", "locfdr", "qvalue", "meanInfRV") %>%
    dplyr::arrange(pvalue)
write.csv(sig_tss, file = "dymamic_sig_fuzzy_tss.csv", row.names = T)

load("CASTxB6.txp.dynamic.rda")
txp_dat <- as.data.frame(mcols(dy))
sig_txp <- txp_dat %>% 
    dplyr::filter(qvalue < 0.05) %>%
    dplyr::select("tx_id", "gene_id", "symbol", "log10mean", 
                  "stat", "log2FC", "pvalue", "locfdr", "qvalue", "meanInfRV", 
                  "fuzzy_tss_group") %>%
    dplyr::arrange(pvalue)
write.csv(sig_txp, file = "dynamic_sig_txp.csv", row.names = T)

