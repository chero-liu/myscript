source("/home/chenglong.liu/RaD/myscript/seurat_fc.r")
source("/home/chenglong.liu/RaD/myscript/Get_colors.R")
library(dplyr)
library(ggplot2)
library("optparse")

option_list <- list(
    make_option(c("--manualfile"), type="character", default="auto", help="Path to manual annotation file"),
    make_option(c("--barcodefile"), type="character", default="auto", help="Path to barcode file"),
    make_option(c("--input"), type="character", help="Path to Seurat input file", metavar="FILE"),
    make_option(c("--outdir"), type="character", help="Output directory"),
    make_option(c("--palette"), type="character", default="customecol2", help="Color palette to use"),
    make_option(c("--reduct"), type="character", default="umap", help="Reduction method to use"),
    make_option(c("--orderbyFreq"), type="character", default="yes", help="Order by frequency [default %default]")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

manualfile <- args$manualfile
barcodefile <- args$barcodefile
input_path <- args$input
outdir_path <- args$outdir
palette <- args$palette
reduct <- args$reduct
orderbyFreq <- args$orderbyFreq

input <- input_path
outdir <- outdir_path
createDir(outdir)
rds = readSeurat(input)

if (manualfile != "auto" && barcodefile != "auto"){
    manualfile <- read.csv(manualfile, sep="\t")
    barcodefile <- read.csv(barcodefile)
    if(colnames(manualfile)[2] != colnames(barcodefile)[2]){
        stop("Please provide the same column name for manualfile and barcodefile")
    }
    useCol = colnames(manualfile)[2]

    refs = unique(rds@meta.data[,colnames(manualfile)[1]])
    rds@meta.data[,colnames(manualfile)[2]] = ""
    for(i in 1:length(refs)){
        ref = refs[i]
        rds@meta.data[,colnames(manualfile)[2]][which(rds@meta.data[,colnames(manualfile)[1]] == ref)] = manualfile[,2][which(manualfile[,1] == ref)]
    }

    for (rb in rds@meta.data$rawbc){
        if (rb %in% barcodefile[,1]){
            rds@meta.data[,colnames(barcodefile)[2]][which(rds@meta.data$rawbc == rb)] = barcodefile[,2][which(barcodefile[,1] == rb)]
        }
    }
}else if (manualfile != "auto" && barcodefile == "auto") {
    manualfile <- read.csv(manualfile, sep="\t")
    useCol = colnames(manualfile)[2]
    refs = unique(rds@meta.data[,colnames(manualfile)[1]])
    rds@meta.data[,colnames(manualfile)[2]] = ""
    for(i in 1:length(refs)){
        ref = refs[i]
        rds@meta.data[,colnames(manualfile)[2]][which(rds@meta.data[,colnames(manualfile)[1]] == ref)] = manualfile[,2][which(manualfile[,1] == ref)]
    }
}else if (manualfile == "auto" && barcodefile != "auto") {
    barcodefile <- read.csv(barcodefile)
    useCol = colnames(barcodefile)[2]
    rds@meta.data[,colnames(barcodefile)[2]] = ""
    for (rb in rds@meta.data$rawbc){
        if (rb %in% barcodefile[,1]){
            rds@meta.data[,colnames(barcodefile)[2]][which(rds@meta.data$rawbc == rb)] = barcodefile[,2][which(barcodefile[,1] == rb)]
        }
    }
}else {
    stop("Please provide manualfile or barcodefile")
}


rds=rds[,rds@meta.data[,useCol] != ""]
rds=rds[,rds@meta.data[,useCol] != "Delete"]
rds=rds[,rds@meta.data[,useCol] != "delete"]
# rds=rds[,rds@meta.data[,useCol] != "Mix"]
if(paste0(useCol,"_col") %in% colnames(rds@meta.data)){
    rds@meta.data <- rds@meta.data[ , !(names(rds@meta.data) %in% c(paste0(useCol,"_col")))]
}

getManualannoResult(rds,outdir,useCol="new_celltype",addFreqCount="no", saveRDS = 'yes',label = FALSE, legend = 'yes',orderbyFreq = orderbyFreq,reduct = reduct,palette=palette,pointsize=0.5)

