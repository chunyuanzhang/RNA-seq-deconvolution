

suppressMessages(library(linseed))
suppressMessages(library(GEOquery))
suppressMessages(library(dplyr))
suppressMessages(library(DESeq2))
suppressMessages(library(DualSimplex))
suppressMessages(library(optparse))


#-------------------------------------------------------------------------------
# 传递外部参数
#-------------------------------------------------------------------------------

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


#-------------------------------------------------------------------------------
# 读取样本信息
#-------------------------------------------------------------------------------

sampletable <- read.csv(SampleFile, header = TRUE, stringsAsFactors = FALSE)
sampletable <- sampletable %>% 
    filter( !is.na(Transcriptome) ) %>% 
    filter(Species %in% target_species) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("SampleID")

#-------------------------------------------------------------------------------
# 读取原始reads矩阵,并标注化
#-------------------------------------------------------------------------------

reads <- read.delim(ReadsMatrix, check.names = F) %>% 
  select(-`transcript_id(s)`) %>% 
  tibble::column_to_rownames(var = "gene_id")

dds <- DESeqDataSetFromMatrix(
  countData = reads,
  colData   = sampletable,
  design    = ~ GroupID
)

keep <- rowSums(counts(dds) >= 10) >= 5   # 至少5个样本中有 >=10 reads
dds <- dds[keep, ]
vsd <- vst(dds)
data_raw <- assay(vsd)

#-------------------------------------------------------------------------------
# 寻找最佳K值
#-------------------------------------------------------------------------------

# 实在算不出来可以调整
hinge_coef <- 2.3
der_coef   <- 0.028     # 固定这一套

set.seed(42)
k_values <- 4:18
n_restarts <- 15        # 推荐10次

best_error <- rep(Inf, length(k_values))
best_sol   <- vector("list", length(k_values))

for (i in seq_along(k_values)) {
  K <- k_values[i]
  cat("测试 K =", K, "\n")
  
  for (r in 1:n_restarts) {
    dso <- DualSimplexSolver$new()
    dso$set_data(data_raw)
    dso$project(K)
    dso$init_solution("random")
    
    dso$optim_solution(5000, optim_config(
      coef_hinge_H = hinge_coef,
      coef_hinge_W = hinge_coef,
      coef_der_X   = der_coef,
      coef_der_Omega = der_coef
    ))
    
    solution <- dso$finalize_solution()
    recon <- solution$W %*% solution$H
    err <- norm(data_raw - recon, type = "F")
    
    if (err < best_error[i] && !is.nan(err)) {
      best_error[i] <- err
      best_sol[[i]] <- solution   # 保存最好的解
    }
  }
  
  cat("K =", K, " | Best Error =", round(best_error[i], 3), "\n")
}

# 查看结果
Best_Error <- data.frame(K = k_values, Best_Error = best_error)
outfile <- paste0(OutPath, "/Best_Error.", target_species, ".csv")
write.csv(x = Best_Error, file = outfile, quote = F, row.names = F )





