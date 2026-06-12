library(linseed)
library(GEOquery)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggpubr)

suppressMessages(library(optparse))


# ------------------------------------------------------------------------------
# 命令行参数传递
# ------------------------------------------------------------------------------

option_list <- list(
  make_option("--SampleFile", type="character", default=NULL, help="config/samples.csv"),
  make_option("--Species", type="character", default=NULL, help="物种"),
  make_option("--ReadsMatrix",type="character",  default=NULL, help="转录组比对原始reads矩阵" ),
  make_option("--OutPath", type = "character", default = "result", help = "结果输出路径")
)

args <- parse_args(OptionParser(option_list=option_list))

SampleFile <- args$SampleFile
target_species <- args$Species
ReadsMatrix <- args$ReadsMatrix
OutPath <- args$OutPath


# ------------------------------------------------------------------------------
# 函数
# ------------------------------------------------------------------------------

# 绘制marker基因热图
plot_signatures_heatmap <- function(signatures, top_n = 20, Species) {
  
  # 取每个细胞类型 top N 个基因
  top_genes <- lapply(1:ncol(signatures), function(i) {
    order(signatures[, i], decreasing = TRUE)[1:top_n]
  })
  
  top_genes_idx <- unique(unlist(top_genes))
  sig_subset <- signatures[top_genes_idx, ]
  
  # 标准化用于热图
  sig_scaled <- t(scale(t(sig_subset)))
  
  # 绘制热图
  p <- pheatmap::pheatmap(sig_scaled,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = F,
           scale = "none",
           main = paste("Top", top_n, "Marker Genes per Cell Type for ", Species),
           color = colorRampPalette(c("blue", "white", "red"))(100),
           fontsize = 8,
           annotation_names_col = TRUE)
  p <- as.ggplot(p)
  
  return(p)
}


# Reads矩阵转换为vst矩阵，并去除线粒体基因
Reads2Vstmat <- function(ReadsMatrix, sampletable){
  ReadsMatrix <- ReadsMatrix
  sampletable <- sampletable
  
  dds <- DESeqDataSetFromMatrix(
    countData = ReadsMatrix,
    colData   = sampletable,
    design    = ~ GroupID
  )
  
  keep <- rowSums(counts(dds) >= 10) >= 5   # 至少5个样本中有 >=10 reads
  dds <- dds[keep, ]
  dds <- dds[ rowSums( counts(dds) ) > 10 , ]
  
  vsd <- vst(dds)
  vst_mat <- assay(vsd)
  
  # 去除线粒体基因
  mt_genes <- grep("_mgp", rownames(vst_mat), value=TRUE)
  vst_mat <- vst_mat[!rownames(vst_mat) %in% mt_genes, ]
  mt_genes <- grep("KEG62", rownames(vst_mat), value=TRUE)
  vst_mat <- vst_mat[!rownames(vst_mat) %in% mt_genes, ]
  
  return(vst_mat)
}


# 解卷积
linseed_lo <- function(data, K){
  
    data <- data
    K <- K
    # 用矩阵创建 LinseedObject
    lo <- LinseedObject$new(data)
    
    # evaluate all pairwise collinearity coefficients
    lo$calculatePairwiseLinearity()
    lo$calculateSpearmanCorrelation()
    lo$calculateSignificanceLevel(100)
    lo$significancePlot(0.01)
    
    # 过滤掉不显著的
    lo$filterDatasetByPval(0.01)
    
    # 奇异值分解方差解释图
    #lo$svdPlot()
    
    lo$setCellTypeNumber(K)
    lo$project("filtered") # projecting full dataset
    lo$projectionPlot(color="filtered")
    
    # 解卷积
    lo$project("filtered")
    lo$smartSearchCorners(dataset="filtered", error="norm")
    lo$deconvolveByEndpoints()
    #plotProportions(lo$proportions)
    
    # lets select 100 genes closest to the simplex corners 
    lo$selectGenes(100)
    #lo$tsnePlot()

    return(lo)
}





# ------------------------------------------------------------------------------
# 读取外部文件
# ------------------------------------------------------------------------------

# 样本信息
sampletable <- read.csv(SampleFile, header = TRUE, stringsAsFactors = FALSE)

# 读取转录组表达矩阵
reads <- read.delim(ReadsMatrix, check.names = F) %>% 
  select(-`transcript_id(s)`) %>% 
  tibble::column_to_rownames(var = "gene_id") 

vst_mat <- Reads2Vstmat(ReadsMatrix = reads, 
                                sampletable = sampletable %>% 
                                filter(Species %in% !!Species) %>% 
                                tibble::column_to_rownames("SampleID")
                        )



# 解卷积，结果可能不稳定，我们需要多次尝试，并提取相对稳定的结果


lo.K8 <- linseed_lo(data = vst_mat, K = 8)
lo.K9 <- linseed_lo(data = vst_mat, K = 9)
lo.K10 <- linseed_lo(data = vst_mat, K = 10)


