、suppressPackageStartupMessages(library(SummarizedExperiment))
library(tximeta)
if (packageVersion("fishpond") < "1.99.24") {
    install_github("mikelove/fishpond")
}
library(fishpond)

set.seed(1)

############################################
## read in the sample names from 'quants' ##
############################################

dir <- "/proj/milovelab/love/proj/ase/osteoblast_quant/quants"
samps <- list.files(dir)
files <- file.path(dir, samps, "quant.sf")

###############################################################
## define different levels of aggregation                   ##
###############################################################

# gene level analysis
cross <- c("129xB6", "CASTxB6")
for (c in cross) {
    load("/osteoblast-quant/data/gse_filtered.rda")
    y <- gse
    y <- y[,y$cross == c]
    y <- labelKeep(y)
    y <- y[mcols(y)$keep,]
    # remove genes with no information (this code specific to our data)
    n <- ncol(y)/2
    mcols(y)$someInfo <- rowSums(abs(assay(y, "infRep1")[,y$allele == "a2"] -
                                         assay(y, "infRep1")[,y$allele == "a1"]) < 1) < n
    y <- y[mcols(y)$someInfo,]
    # no scaling, because we are comparing alleles within samples
    # global Swish test
    global <- swish(y, x="allele", pair="day") #e.g. performed a pairwise signed rank test within each day
    # access infRV
    global <- computeInfRV(global)
    save(global, file = paste0(c, ".gene.global.rda"))
    
    # dynamics test
    dy <- swish(y, x="allele", pair = "day", cov="day", cor = "pearson")
    dy <- computeInfRV(dy)
    save(dy, file = paste0(c, ".gene.dynamic.rda"))
}

# txp level analysis
cross <- c("129xB6", "CASTxB6")
for (c in cross) {
    load("/osteoblast-quant/data/se_filtered.rda")
    y <- se
    y <- y[,y$cross == c]
    y <- labelKeep(y)
    y <- y[mcols(y)$keep,]
    
    # remove genes with no information (this code specific to our data)
    n <- ncol(y)/2
    mcols(y)$someInfo <- rowSums(abs(assay(y, "infRep1")[,y$allele == "a2"] -
                                         assay(y, "infRep1")[,y$allele == "a1"]) < 1) < n
    y <- y[mcols(y)$someInfo,]

    # global Swish test
    global <- swish(y, x="allele", pair="day") #e.g. performed a pairwise signed rank test within each day
    # assess the InfRV
    global <- computeInfRV(global)
    
    load("/osteoblast-quant/data/fuzzy_50bp_tss_se_filtered.rda")
    index <- pmatch(names(global), unlist(mcols(gse)$tx_id))
    mcols(global)$fuzzy_tss_group <- names(unlist(mcols(gse)$tx_id)[index])
    save(global, file = paste0(c, ".txp.global.rda"))
    
    # dynamics test
    dy <- swish(y, x="allele", pair = "day", cov="day", cor = "pearson")
    dy <- computeInfRV(dy)
    index <- pmatch(names(dy), unlist(mcols(gse)$tx_id))
    mcols(dy)$fuzzy_tss_group <- names(unlist(mcols(gse)$tx_id)[index])
    save(dy, file = paste0(c, ".txp.dynamic.rda"))
}

# fuzzy tss level analysis
cross <- c("129xB6", "CASTxB6")
for (c in cross) {
    load("/osteoblast-quant/data/fuzzy_50bp_tss_se_filtered.rda")
    y <- gse
    y <- y[,y$cross == c]
    y <- labelKeep(y)
    y <- y[mcols(y)$keep,]
    rowData(y)$start_position <- sapply(strsplit(rownames(y), "-"), "[", 2)
    
    # remove genes with no information (this code specific to our data)
    n <- ncol(y)/2
    mcols(y)$someInfo <- rowSums(abs(assay(y, "infRep1")[,y$allele == "a2"] -
                                         assay(y, "infRep1")[,y$allele == "a1"]) < 1) < n
    y <- y[mcols(y)$someInfo,]
    
    # no scaling, because we are comparing alleles within samples
    # global Swish test
    global <- swish(y, x="allele", pair="day") #e.g. performed a pairwise signed rank test within each day
    global <- computeInfRV(global)
    save(global, file = paste0(c, ".fuzzy.tss.global.rda"))
    
    # dynamics test
    dy <- swish(y, x="allele", pair = "day", cov="day", cor = "pearson")
    dy <- computeInfRV(dy)
    save(dy, file = paste0(c, ".fuzzy.tss.dynamic.rda"))
}
