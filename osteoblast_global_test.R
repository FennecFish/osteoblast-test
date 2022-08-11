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

set.seed(1)

#######################################################################################
#################################### Upset Plot ######################################
#######################################################################################

load("CASTxB6.gene.global.rda")

gene_dat <- as.data.frame(mcols(global))
# select genes with significant AI at 5% FDR
sig_gene <- gene_dat %>% 
    dplyr::filter(qvalue < 0.05) %>% 
    dplyr::select(ends_with("id"), "symbol", "log2FC", "qvalue")
length(unique(sig_gene$group_id)) # 5701

load("CASTxB6.fuzzy.tss.global.rda")
# select tss groups with significant AI at 5% FDR
tss_dat <- as.data.frame(mcols(global))
sig_tss <- tss_dat %>% 
    dplyr::filter(qvalue < 0.05) %>%
    dplyr::select(ends_with("id"), "symbol", "log2FC", "qvalue")
length(unique(sig_tss$gene_id)) # 5573

load("CASTxB6.txp.global.rda")
# select txps groups with significant AI at 5% FDR
txp_dat <- as.data.frame(mcols(global))
sig_txp <- txp_dat %>% 
    dplyr::filter(qvalue < 0.05) %>%
    dplyr::select(ends_with("id"), "symbol", "log2FC", "qvalue")
length(unique(sig_txp$gene_id)) # 6116

length(unique(intersect(sig_gene$group_id,sig_tss$gene_id))) # 4873
length(unique(intersect(sig_gene$group_id,sig_txp$gene_id))) # 4880
length(unique(intersect(sig_txp$gene_id,sig_tss$gene_id))) # 5284
length(unique(intersect(sig_gene$group_id,intersect(sig_txp$gene_id,sig_tss$gene_id)))) # 4625

# making upset plot
data <- list(Gene =unique(sig_gene$group_id), TSS =unique(sig_tss$gene_id), Txp =unique(sig_txp$gene_id))
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

#######################################################################################
###################### A complete List of different direction AI ######################
#######################################################################################

load("CASTxB6.gene.global.rda")

gene_dat <- as.data.frame(mcols(global))
sig_gene <- gene_dat %>% 
    dplyr::filter(qvalue < 0.05) %>% 
    dplyr::select(ends_with("id"), "symbol", "log2FC", "qvalue")

load("CASTxB6.fuzzy.tss.global.rda")
tss_dat <- as.data.frame(mcols(global))
sig_tss <- tss_dat %>% 
    dplyr::filter(qvalue < 0.05) %>%
    dplyr::select(ends_with("id"), "symbol", "log2FC", "qvalue")

dat <- inner_join(sig_gene, sig_tss, by = c("group_id" = "gene_id"))
# flag tss groups with opposite direction of AI under the same genes
dat <- dat %>% group_by(group_id) %>% mutate(direction = mean(sign(log2FC.x) == sign(log2FC.y)) != 1)
diff_iso <- dat %>% dplyr::filter(direction == "TRUE") %>% dplyr::select(-c("symbol.y","tx_id.x","direction")) %>%
    dplyr::rename(gene_id = group_id,
                  symbol = symbol.x,
                  log2FC_gene = log2FC.x,
                  qvalue_gene = qvalue.x,
                  tss_group_id = group_id.y,
                  log2FC_tss = log2FC.y,
                  qvalue_tss = qvalue.y)
diff_iso$num_txp <- lengths(diff_iso$tx_id.y)
diff_iso <- diff_iso %>% dplyr::select(-c("tx_id.y"))

diff_gene <- diff_iso %>% group_by(tss_group_id) %>% mutate(direction = mean(sign(log2FC_gene) == sign(log2FC_tss)) != 1)
diff_gene <- diff_gene %>% dplyr::filter(direction == "TRUE") %>% dplyr::select(-c("direction"))
length(unique(diff_gene$gene_id)) #134
diff_gene <- diff_gene %>% dplyr::select(tss_group_id, symbol, log2FC_gene, qvalue_gene, log2FC_tss, qvalue_tss)

write.csv(as.data.frame(diff_gene), file = paste0("Genes_w_Isoforms_Diff_Direction.csv"), row.names = T)

#######################################################################################
###################### Output sig output on all three resolution ######################
#######################################################################################

load("CASTxB6.gene.global.rda")

gene_dat <- as.data.frame(mcols(global))
sig_gene <- gene_dat %>% 
    dplyr::filter(qvalue < 0.05) %>% 
    dplyr::select("symbol", "log10mean", "stat", "log2FC"
                  , "pvalue", "locfdr", "qvalue", "meanInfRV") %>%
    dplyr::arrange(pvalue)
write.csv(sig_gene, file = "global_sig_gene.csv", row.names = T)

load("CASTxB6.fuzzy.tss.global.rda")
tss_dat <- as.data.frame(mcols(global))
sig_tss <- tss_dat %>% 
    dplyr::filter(qvalue < 0.05) %>%
    dplyr::select("gene_id", "symbol", "start_position", "log10mean", 
                  "stat", "log2FC", "pvalue", "locfdr", "qvalue", "meanInfRV") %>%
    dplyr::arrange(pvalue)
write.csv(sig_tss, file = "global_sig_fuzzy_tss.csv", row.names = T)

load("CASTxB6.txp.global.rda")
txp_dat <- as.data.frame(mcols(global))
sig_txp <- txp_dat %>% 
    dplyr::filter(qvalue < 0.05) %>%
    dplyr::select("tx_id", "gene_id", "symbol", "log10mean", 
                  "stat", "log2FC", "pvalue", "locfdr", "qvalue", "meanInfRV", 
                  "fuzzy_tss_group") %>%
    dplyr::arrange(pvalue)
write.csv(sig_txp, file = "global_sig_txp.csv", row.names = T)

#######################################################################################
######################### global significant with TSS percentage ######################
#######################################################################################

load("CASTxB6.gene.global.rda")
gene_dat <- as.data.frame(mcols(global))
sig_gene <- gene_dat %>% 
    dplyr::filter(qvalue < 0.05) %>% 
    dplyr::select(ends_with("id"), "symbol", "log2FC", "qvalue")

load("CASTxB6.fuzzy.tss.global.rda")
tss_dat <- as.data.frame(mcols(global))
sig_tss <- tss_dat %>% 
    dplyr::filter(qvalue < 0.05) %>%
    dplyr::select(ends_with("id"), "symbol", "log2FC", "qvalue")

dat <- inner_join(sig_gene, sig_tss, by = c("group_id" = "gene_id"))
perc <- dat %>% group_by(group_id.y)  %>% 
    summarise(count = mean(sign(log2FC.x) == sign(log2FC.y)))
dat1 <- left_join(perc, sig_tss, by = c("group_id.y" = "group_id")) %>% arrange(gene_id,count)
# select the first row of each gene since count 0 will be the first if they exists
dat1 <- dat1[!duplicated(dat1$gene_id ),]
perc <- dat1 %>% group_by(count) %>% summarise(n = n ())

#######################################################################################
###################################### Fuca2 Plot######################################
#######################################################################################

gene_dat %>% dplyr::filter(symbol == "Fuca2")
tss_dat %>% dplyr::filter(symbol == "Fuca2")

load("CASTxB6.txp.global.rda")
txp_dat <- as.data.frame(mcols(global))
txp_dat %>% dplyr::filter(symbol == "Fuca2")

ah <- AnnotationHub()
edb <- ah[["AH89211"]]

load("CASTxB6.fuzzy.tss.global.rda")
plotAllelicGene(global, symbol="Fuca2", db=edb, labels=list(a2="B6",a1="CAST")) # good example

rowMeans(assay(global, "abundance")[which(mcols(global)$symbol == "Fuca2"),])
