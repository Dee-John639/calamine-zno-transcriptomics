# scripts/01_metadata_gse121212.R

library(tidyverse)

# 读取SRA metadata
sra <- read.delim("data/raw/GSE121212/SraRunTable.txt",
                  sep = "\t", check.names = FALSE)

# 查看关键列名（因SRA版本而异）
# 常见列：Run, Sample.Name, patient's condition, skin type, etc.
names(sra)

# 根据你下载到的实际列名，提取疾病和lesion状态
# 从GSM页面已知字段为 "patient's condition" 和 "skin type"[reference:1]
meta <- sra %>%
  transmute(
    sample_id = Run,                          # SRR号
    gsm_id    = Sample.Name,                   # GSM号
    disease   = `patient's condition`,         # AD / PsO / healthy
    lesion_status = `skin type`,               # lesional / non-lesional
    donor_id  = str_extract(gsm_id, "AD_\\d+|PSO_\\d+|CTRL_\\d+")  # 从GSM命名提取
  )

# 核对donor_id是否成功提取
table(meta$donor_id, useNA = "always")

# 保存
write.csv(meta, "data/metadata/GSE121212_metadata.csv", row.names = FALSE)


# scripts/01_metadata_emtrab5734.R

library(tidyverse)

sdrf <- read.delim("data/raw/E-MTAB-5734/E-MTAB-5734.sdrf.txt",
                   sep = "\t", check.names = FALSE)

# 查看实际列名
names(sdrf)

# 根据实际SDRF列名提取（列名可能略有差异）
meta_zno <- sdrf %>%
  transmute(
    sample_id = `Source Name`,
    material  = `Characteristics[compound]`,
    form      = `Characteristics[particle size]`,
    dose      = as.numeric(`Characteristics[dose]`),
    time      = as.numeric(`Characteristics[time]`),
    ena_run   = `Comment[ENA_RUN]`
  )

# 检查处理类型分布
table(meta_zno$material, meta_zno$time)

# 保存
write.csv(meta_zno, "data/metadata/E-MTAB-5734_metadata.csv", row.names = FALSE)


# scripts/02_qc_gse121212.R

library(DESeq2)
library(ggplot2)

counts <- read.delim("data/raw/GSE121212/GSE121212_readcount.txt.gz",
                     row.names = 1, check.names = FALSE)
meta <- read.csv("data/metadata/GSE121212_metadata.csv")

# 确认列顺序匹配
all(colnames(counts) %in% meta$gsm_id)
counts <- counts[, meta$gsm_id]
meta <- meta[match(colnames(counts), meta$gsm_id), ]

# DESeq2对象
dds <- DESeqDataSetFromMatrix(counts, colData = meta,
                              design = ~ lesion_status)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# VST
vsd <- vst(dds, blind = TRUE)

# PCA
plotPCA(vsd, intgroup = c("disease", "lesion_status")) +
  theme_bw()

# 样本相关性
sampleDists <- dist(t(assay(vsd)))
pheatmap::pheatmap(as.matrix(sampleDists),
                   annotation_col = meta[, c("disease", "lesion_status")])

# scripts/03_de_emtrab5734.R
setwd('C:/Users/DeeJo/calamine_zinc_project/')

library(tximport)
library(DESeq2)
library(tidyverse)
library(GenomicFeatures)

sdrf <- read.delim(
  "C:/Users/DeeJo/calamine_zinc_project/data/raw/E-MTAB-5734/E-MTAB-5734.sdrf.txt",
  sep = "\t", check.names = FALSE, header = TRUE
)

meta <- data.frame(
  sample_id = sdrf[[1]],
  material  = sdrf[[32]],
  dose      = sdrf[[33]],
  time      = sdrf[[35]],
  ena_run   = sdrf[[30]],
  stringsAsFactors = FALSE
)

# 检查
table(meta$material)
head(meta)

# 保存
write.csv(meta,
          "C:/Users/DeeJo/calamine_zinc_project/data/metadata/E-MTAB-5734_metadata.csv",
          row.names = FALSE)

meta <- read.csv(
  "C:/Users/DeeJo/calamine_zinc_project/data/metadata/E-MTAB-5734_metadata.csv",
  stringsAsFactors = FALSE
)

target <- c(
  "none",
  "ZnCl2",
  "ZnO microparticles",
  "ZnO NPs coated with triethoxycaprylsilane NM111",
  "uncoated ZnO NPs NM110"
)

meta_sub <- meta %>%
  filter(material %in% target) %>%
  mutate(
    treatment = case_when(
      material == "none" ~ "Control",
      material == "ZnCl2" ~ "ZnCl2",
      material == "ZnO microparticles" ~ "ZnO_MP",
      material == "ZnO NPs coated with triethoxycaprylsilane NM111" ~ "ZnO_NP_coated",
      material == "uncoated ZnO NPs NM110" ~ "ZnO_NP_uncoated"
    ),
    treatment = factor(treatment,
                       levels = c("Control", "ZnCl2", "ZnO_MP", "ZnO_NP_coated", "ZnO_NP_uncoated")),
    time = factor(time)
  )

table(meta_sub$treatment, meta_sub$time)

txi <- readRDS("C:/Users/DeeJo/calamine_zinc_project/data/processed/txi_emtrab5734.rds")

# 检查样本是否全部匹配
all(meta_sub$ena_run %in% colnames(txi$counts))

# 按 metadata 顺序重排
txi$counts    <- txi$counts[, meta_sub$ena_run]
txi$abundance <- txi$abundance[, meta_sub$ena_run]
txi$length    <- txi$length[, meta_sub$ena_run]

library(DESeq2)

dds <- DESeqDataSetFromTximport(
  txi, colData = meta_sub,
  design = ~ time + treatment
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

dir.create("C:/Users/DeeJo/calamine_zinc_project/results/differential_expression",
           recursive = TRUE, showWarnings = FALSE)

contrasts <- list(
  ZnO_NP_coated_vs_Control   = c("treatment", "ZnO_NP_coated",   "Control"),
  ZnO_NP_uncoated_vs_Control = c("treatment", "ZnO_NP_uncoated", "Control"),
  ZnO_MP_vs_Control          = c("treatment", "ZnO_MP",          "Control"),
  ZnCl2_vs_Control           = c("treatment", "ZnCl2",           "Control")
)

outDir <- "C:/Users/DeeJo/calamine_zinc_project/results/differential_expression"

for (nm in names(contrasts)) {
  res <- results(dds, contrast = contrasts[[nm]])
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    select(gene, log2FoldChange, lfcSE, stat, pvalue, padj)
  write.csv(res_df, file.path(outDir, paste0(nm, ".csv")), row.names = FALSE)
  cat(nm, ": ", sum(res$padj < 0.05, na.rm = TRUE), " DEGs\n")
}

vsd <- vst(dds, blind = TRUE)

# PCA
library(ggplot2)

pcaData <- plotPCA(vsd, intgroup = c("treatment", "time"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = treatment, shape = time)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(legend.position = "right")

# 样本相关
sampleDists <- dist(t(assay(vsd)))

library(pheatmap)
library(RColorBrewer)

# 短标签：处理组缩写 + 时间
treatment_short <- c(
  "Control"         = "Ctrl",
  "ZnCl2"           = "Zn2",
  "ZnO_MP"          = "MP",
  "ZnO_NP_coated"   = "NPc",
  "ZnO_NP_uncoated" = "NPu"
)

short_labels <- paste0(
  treatment_short[as.character(meta_sub$treatment)],
  "_", meta_sub$time,
  "_", seq_len(nrow(meta_sub))
)

# 样本距离
sampleDists <- dist(t(assay(vsd)))
distMat <- as.matrix(sampleDists)
rownames(distMat) <- short_labels
colnames(distMat) <- short_labels

# 注释
ann <- data.frame(
  Treatment = meta_sub$treatment,
  Time      = meta_sub$time,
  row.names = short_labels
)

dir.create("C:/Users/DeeJo/calamine_zinc_project/figures",
           recursive = TRUE, showWarnings = FALSE)

pheatmap(
  distMat,
  annotation_col = ann,
  annotation_row = ann,
  show_rownames = FALSE,   # 152个标签太多，直接隐藏
  show_colnames = FALSE,
  fontsize = 8,
  width = 10,
  height = 9,
  filename = "C:/Users/DeeJo/calamine_zinc_project/figures/qc_emtrab5734_sampleDist.png"
)

saveRDS(dds, "C:/Users/DeeJo/calamine_zinc_project/data/processed/dds_emtrab5734.rds")

outDir <- "C:/Users/DeeJo/calamine_zinc_project/results/differential_expression"

for (f in list.files(outDir, pattern = "\\.csv$", full.names = TRUE)) {
  df <- read.csv(f)
  n_up   <- sum(df$padj < 0.05 & df$log2FoldChange >  1, na.rm = TRUE)
  n_down <- sum(df$padj < 0.05 & df$log2FoldChange < -1, na.rm = TRUE)
  cat(basename(f), ": up =", n_up, ", down =", n_down, "\n")
}

library(fgsea)
library(ggplot2)

outDir <- "C:/Users/DeeJo/calamine_zinc_project/results/differential_expression"

# 读入四个contrast
res_np_u  <- read.csv(file.path(outDir, "ZnO_NP_uncoated_vs_Control.csv"))
res_np_c  <- read.csv(file.path(outDir, "ZnO_NP_coated_vs_Control.csv"))
res_mp    <- read.csv(file.path(outDir, "ZnO_MP_vs_Control.csv"))
res_zn2   <- read.csv(file.path(outDir, "ZnCl2_vs_Control.csv"))

length(zn2_up)
length(zn2_down)
sum(is.na(ranks_np_u))
sum(is.na(ranks_mp))
sum(is.na(ranks_zn2))

# 用显著DEG定义gene set
zn2_up   <- res_zn2$gene[res_zn2$padj < 0.05 & res_zn2$stat > 0]
zn2_down <- res_zn2$gene[res_zn2$padj < 0.05 & res_zn2$stat < 0]

length(zn2_up)    # 应该接近 2197
length(zn2_down)  # 应该接近 563

# 放宽 maxSize 到 5000
gsea_np_u <- fgsea(
  list(Zn2_up = zn2_up, Zn2_down = zn2_down),
  ranks_np_u,
  minSize = 15, maxSize = 5000
)

gsea_mp <- fgsea(
  list(Zn2_up = zn2_up, Zn2_down = zn2_down),
  ranks_mp,
  minSize = 15, maxSize = 5000
)

print(gsea_np_u)
print(gsea_mp)

# ZnO_NP_uncoated 中显著但 ZnCl2 不显著的基因
particle_np <- res_np_u$gene[
  res_np_u$padj < 0.05 & 
    (is.na(res_zn2$padj[match(res_np_u$gene, res_zn2$gene)]) |
       res_zn2$padj[match(res_np_u$gene, res_zn2$gene)] > 0.1)
]

# ZnO_MP 中显著但 ZnCl2 不显著的基因
particle_mp <- res_mp$gene[
  res_mp$padj < 0.05 & 
    (is.na(res_zn2$padj[match(res_mp$gene, res_zn2$gene)]) |
       res_zn2$padj[match(res_mp$gene, res_zn2$gene)] > 0.1)
]

length(particle_np)
length(particle_mp)

# 两者交集：稳定的 particle-associated response
particle_shared <- intersect(particle_np, particle_mp)
length(particle_shared)

ranks_np_c <- setNames(res_np_c$stat, res_np_c$gene)

gsea_particle <- fgsea(
  list(Particle = particle_shared),
  ranks_np_c,
  minSize = 15, maxSize = 5000
)

print(gsea_particle)

counts_raw <- read.delim(
  "C:/Users/DeeJo/calamine_zinc_project/data/raw/GSE121212/GSE121212_readcount.txt.gz",
  row.names = NULL, check.names = FALSE
)

dim(counts_raw)
head(counts_raw[, 1:5])

# 第一列是否有重复
sum(duplicated(counts_raw[[1]]))
sum(is.na(counts_raw[[1]]))

library(tidyverse)

counts_raw <- read.delim(
  "C:/Users/DeeJo/calamine_zinc_project/data/raw/GSE121212/GSE121212_readcount.txt.gz",
  row.names = NULL, check.names = FALSE
)

# 第一列列名（可能是空或 "X"）
names(counts_raw)[1]

# 重命名第一列为 "gene"
colnames(counts_raw)[1] <- "gene"

# 按 gene symbol 聚合（2 个重复求和）
counts <- counts_raw %>%
  group_by(gene) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("gene") %>%
  as.matrix()

dim(counts)
# 应该输出 31362 x 147

sra <- read.csv(
  "C:/Users/DeeJo/calamine_zinc_project/data/raw/GSE121212/SraRunTable.csv",
  stringsAsFactors = FALSE
)

names(sra)

meta <- meta %>%
  mutate(
    disease = case_when(
      grepl("^AD_",    sample) ~ "AD",
      grepl("^PSO_",   sample) ~ "PsO",
      grepl("^CTRL_",  sample) ~ "healthy",
      TRUE ~ NA_character_
    ),
    donor_id = str_extract(sample, "(AD|PSO|CTRL)_\\d+"),
    lesion = case_when(
      grepl("lesional$",     sample) & !grepl("non-lesional$", sample) ~ "lesional",
      grepl("non-lesional$", sample) ~ "non-lesional",
      grepl("chronic_lesion$", sample) ~ "chronic_lesion",
      grepl("healthy$",      sample) ~ "healthy",
      TRUE ~ NA_character_
    )
  )

table(meta$disease, meta$lesion, useNA = "always")

library(tidyverse)

meta <- data.frame(
  sample = colnames(counts),
  stringsAsFactors = FALSE
) %>%
  mutate(
    disease = case_when(
      grepl("^AD_",   sample) ~ "AD",
      grepl("^PSO_",  sample) ~ "PsO",
      grepl("^CTRL_", sample) ~ "healthy",
      TRUE ~ NA_character_
    ),
    donor_id = str_extract(sample, "(AD|PSO|CTRL)_\\d+"),
    lesion = case_when(
      grepl("non-lesional$",  sample) ~ "non-lesional",
      grepl("chronic_lesion$", sample) ~ "chronic_lesion",
      grepl("lesional$",      sample) ~ "lesional",
      grepl("healthy$",       sample) ~ "healthy",
      TRUE ~ NA_character_
    )
  )

# 检查
table(meta$disease, meta$lesion, useNA = "always")
sum(is.na(meta$disease))
sum(is.na(meta$lesion))

write.csv(meta,
          "C:/Users/DeeJo/calamine_zinc_project/data/metadata/GSE121212_metadata.csv",
          row.names = FALSE)

library(DESeq2)
library(tidyverse)

# 只保留标准 lesional / non-lesional 的 AD 样本
meta_ad <- meta %>%
  filter(disease == "AD",
         lesion %in% c("lesional", "non-lesional")) %>%
  mutate(
    lesion   = factor(lesion, levels = c("non-lesional", "lesional")),
    donor_id = factor(donor_id)
  )

# 检查配对是否完整
table(meta_ad$donor_id, meta_ad$lesion)

# =============================================================
# GSE121212 四个对比的差异分析
# 依赖：counts 和 meta 已在当前 R 会话中
# 全程使用 base R 提取结果表，避免 dplyr::select 被覆盖
# =============================================================

library(DESeq2)

outDir <- "C:/Users/DeeJo/calamine_zinc_project/results/differential_expression"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# 从 DESeqResults 对象提取完整基因表
extract_res <- function(res, file_name) {
  df <- data.frame(
    gene           = rownames(res),
    log2FoldChange = res$log2FoldChange,
    lfcSE          = res$lfcSE,
    stat           = res$stat,
    pvalue         = res$pvalue,
    padj           = res$padj,
    stringsAsFactors = FALSE
  )
  write.csv(df, file.path(outDir, file_name), row.names = FALSE)
  cat(file_name, ": ",
      sum(df$padj < 0.05 & df$log2FoldChange >  1, na.rm = TRUE), "up, ",
      sum(df$padj < 0.05 & df$log2FoldChange < -1, na.rm = TRUE), "down\n", sep = "")
  return(df)
}


# =============================================================
# 1. AD_L vs AD_NL (paired, 21 donors)
# =============================================================
meta_ad <- meta[meta$disease == "AD" &
                  meta$lesion %in% c("lesional", "non-lesional"), ]

meta_ad$lesion   <- factor(meta_ad$lesion, levels = c("non-lesional", "lesional"))
meta_ad$donor_id <- factor(meta_ad$donor_id)

# 只保留配对完整的 donor
paired_donors <- names(which(table(meta_ad$donor_id) == 2))
meta_ad <- meta_ad[meta_ad$donor_id %in% paired_donors, ]
meta_ad$donor_id <- droplevels(meta_ad$donor_id)

cat("AD paired donors:", length(unique(meta_ad$donor_id)),
    " samples:", nrow(meta_ad), "\n")

counts_ad <- counts[, meta_ad$sample]
dds_ad <- DESeqDataSetFromMatrix(counts_ad, colData = meta_ad,
                                 design = ~ donor_id + lesion)
dds_ad <- dds_ad[rowSums(counts(dds_ad)) >= 10, ]
dds_ad <- DESeq(dds_ad)

res_ad <- results(dds_ad, contrast = c("lesion", "lesional", "non-lesional"))
res_ad_df <- extract_res(res_ad, "AD_L_vs_AD_NL.csv")

saveRDS(dds_ad,
        "C:/Users/DeeJo/calamine_zinc_project/data/processed/dds_ad_gse121212.rds")


# =============================================================
# 2. AD_L vs Healthy (unpaired)
# =============================================================
meta_ad_h <- meta[(meta$disease == "AD" & meta$lesion == "lesional") |
                    meta$disease == "healthy", ]

meta_ad_h$disease <- factor(meta_ad_h$disease, levels = c("healthy", "AD"))

cat("AD_L vs Healthy samples:", nrow(meta_ad_h), "\n")

counts_ad_h <- counts[, meta_ad_h$sample]
dds_ad_h <- DESeqDataSetFromMatrix(counts_ad_h, colData = meta_ad_h,
                                   design = ~ disease)
dds_ad_h <- dds_ad_h[rowSums(counts(dds_ad_h)) >= 10, ]
dds_ad_h <- DESeq(dds_ad_h)

res_ad_h <- results(dds_ad_h, contrast = c("disease", "AD", "healthy"))
res_ad_h_df <- extract_res(res_ad_h, "AD_L_vs_Healthy.csv")

saveRDS(dds_ad_h,
        "C:/Users/DeeJo/calamine_zinc_project/data/processed/dds_ad_h_gse121212.rds")


# =============================================================
# 3. PsO_L vs PsO_NL (paired)
# =============================================================
meta_pso <- meta[meta$disease == "PsO" &
                   meta$lesion %in% c("lesional", "non-lesional"), ]

meta_pso$lesion   <- factor(meta_pso$lesion, levels = c("non-lesional", "lesional"))
meta_pso$donor_id <- factor(meta_pso$donor_id)

paired_pso <- names(which(table(meta_pso$donor_id) == 2))
meta_pso <- meta_pso[meta_pso$donor_id %in% paired_pso, ]
meta_pso$donor_id <- droplevels(meta_pso$donor_id)

cat("PsO paired donors:", length(unique(meta_pso$donor_id)),
    " samples:", nrow(meta_pso), "\n")

counts_pso <- counts[, meta_pso$sample]
dds_pso <- DESeqDataSetFromMatrix(counts_pso, colData = meta_pso,
                                  design = ~ donor_id + lesion)
dds_pso <- dds_pso[rowSums(counts(dds_pso)) >= 10, ]
dds_pso <- DESeq(dds_pso)

res_pso <- results(dds_pso, contrast = c("lesion", "lesional", "non-lesional"))
res_pso_df <- extract_res(res_pso, "PsO_L_vs_PsO_NL.csv")

saveRDS(dds_pso,
        "C:/Users/DeeJo/calamine_zinc_project/data/processed/dds_pso_gse121212.rds")


# =============================================================
# 4. PsO_L vs Healthy (unpaired)
# =============================================================
meta_pso_h <- meta[(meta$disease == "PsO" & meta$lesion == "lesional") |
                     meta$disease == "healthy", ]

meta_pso_h$disease <- factor(meta_pso_h$disease, levels = c("healthy", "PsO"))

cat("PsO_L vs Healthy samples:", nrow(meta_pso_h), "\n")

counts_pso_h <- counts[, meta_pso_h$sample]
dds_pso_h <- DESeqDataSetFromMatrix(counts_pso_h, colData = meta_pso_h,
                                    design = ~ disease)
dds_pso_h <- dds_pso_h[rowSums(counts(dds_pso_h)) >= 10, ]
dds_pso_h <- DESeq(dds_pso_h)

res_pso_h <- results(dds_pso_h, contrast = c("disease", "PsO", "healthy"))
res_pso_h_df <- extract_res(res_pso_h, "PsO_L_vs_Healthy.csv")

saveRDS(dds_pso_h,
        "C:/Users/DeeJo/calamine_zinc_project/data/processed/dds_pso_h_gse121212.rds")


# =============================================================
# 5. 汇总
# =============================================================
cat("\n===== 汇总 =====\n")
cat("AD_L_vs_AD_NL    : ",
    sum(res_ad_df$padj < 0.05 & res_ad_df$log2FoldChange >  1, na.rm = TRUE), "up, ",
    sum(res_ad_df$padj < 0.05 & res_ad_df$log2FoldChange < -1, na.rm = TRUE), "down\n")
cat("AD_L_vs_Healthy  : ",
    sum(res_ad_h_df$padj < 0.05 & res_ad_h_df$log2FoldChange >  1, na.rm = TRUE), "up, ",
    sum(res_ad_h_df$padj < 0.05 & res_ad_h_df$log2FoldChange < -1, na.rm = TRUE), "down\n")
cat("PsO_L_vs_PsO_NL  : ",
    sum(res_pso_df$padj < 0.05 & res_pso_df$log2FoldChange >  1, na.rm = TRUE), "up, ",
    sum(res_pso_df$padj < 0.05 & res_pso_df$log2FoldChange < -1, na.rm = TRUE), "down\n")
cat("PsO_L_vs_Healthy : ",
    sum(res_pso_h_df$padj < 0.05 & res_pso_h_df$log2FoldChange >  1, na.rm = TRUE), "up, ",
    sum(res_pso_h_df$padj < 0.05 & res_pso_h_df$log2FoldChange < -1, na.rm = TRUE), "down\n")

library(data.table)

# 从 GENCODE GTF 提取 Ensembl ID → Symbol 映射
gtf <- fread(
  "C:/Users/DeeJo/calamine_zinc_project/data/reference/gencode.v44.annotation.gtf.gz",
  sep = "\t", header = FALSE, skip = "##"
)
gtf_gene <- gtf[V3 == "gene"]

gene_map <- data.frame(
  ensembl = sub('.*gene_id "([^"]+)".*', '\\1', gtf_gene$V9),
  symbol  = sub('.*gene_name "([^"]+)".*', '\\1', gtf_gene$V9),
  stringsAsFactors = FALSE
)
gene_map$ensembl <- sub("\\..*$", "", gene_map$ensembl)
gene_map <- gene_map[!duplicated(gene_map$ensembl), ]

# 检查 E-MTAB-5734 结果中的 gene ID 格式
head(res_np_u$gene)
# 如果带版本号，去掉
res_np_u$gene <- sub("\\..*$", "", res_np_u$gene)
res_mp$gene   <- sub("\\..*$", "", res_mp$gene)
res_zn2$gene  <- sub("\\..*$", "", res_zn2$gene)
res_np_c$gene <- sub("\\..*$", "", res_np_c$gene)

# 转成 Symbol
res_np_u$symbol <- gene_map$symbol[match(res_np_u$gene, gene_map$ensembl)]
res_mp$symbol   <- gene_map$symbol[match(res_mp$gene,   gene_map$ensembl)]
res_zn2$symbol  <- gene_map$symbol[match(res_zn2$gene,  gene_map$ensembl)]
res_np_c$symbol <- gene_map$symbol[match(res_np_c$gene, gene_map$ensembl)]

cat("有 symbol 的比例：",
    round(mean(!is.na(res_zn2$symbol)) * 100, 1), "%\n")

outDir <- "C:/Users/DeeJo/calamine_zinc_project/results/differential_expression"

read_res <- function(file_name) {
  df <- read.csv(file.path(outDir, file_name), stringsAsFactors = FALSE)
  df$gene <- sub("\\..*$", "", df$gene)
  df$symbol <- gene_map$symbol[match(df$gene, gene_map$ensembl)]
  df <- df[!is.na(df$symbol) & !is.na(df$stat), ]
  df <- df[!duplicated(df$symbol), ]
  df
}

res_np_u <- read_res("ZnO_NP_uncoated_vs_Control.csv")
res_mp   <- read_res("ZnO_MP_vs_Control.csv")
res_zn2  <- read_res("ZnCl2_vs_Control.csv")
res_np_c <- read_res("ZnO_NP_coated_vs_Control.csv")

# 构建 ranked vectors（用 stat 列）
ranks_np_u <- setNames(res_np_u$stat, res_np_u$symbol)
ranks_mp   <- setNames(res_mp$stat,   res_mp$symbol)
ranks_zn2  <- setNames(res_zn2$stat,  res_zn2$symbol)
ranks_np_c <- setNames(res_np_c$stat, res_np_c$symbol)

# 排序检查
head(sort(ranks_zn2, decreasing = TRUE), 5)

read_disease <- function(file_name) {
  df <- read.csv(file.path(outDir, file_name), stringsAsFactors = FALSE)
  df <- df[!is.na(df$stat) & !is.na(df$padj), ]
  df <- df[!duplicated(df$gene), ]
  df
}

res_ad_nl    <- read_disease("AD_L_vs_AD_NL.csv")
res_ad_h     <- read_disease("AD_L_vs_Healthy.csv")
res_pso_nl   <- read_disease("PsO_L_vs_PsO_NL.csv")
res_pso_h    <- read_disease("PsO_L_vs_Healthy.csv")

# 定义 gene set（padj < 0.05，方向由 stat 决定）
make_sets <- function(df, tag) {
  up   <- df$gene[df$padj < 0.05 & df$stat > 0]
  down <- df$gene[df$padj < 0.05 & df$stat < 0]
  list(
    up   = up,
    down = down
  )
}

ad_nl_sets  <- make_sets(res_ad_nl,  "AD_L_vs_AD_NL")
ad_h_sets   <- make_sets(res_ad_h,   "AD_L_vs_Healthy")
pso_nl_sets <- make_sets(res_pso_nl, "PsO_L_vs_PsO_NL")
pso_h_sets  <- make_sets(res_pso_h,  "PsO_L_vs_Healthy")

cat("AD_L_vs_AD_NL   : up =", length(ad_nl_sets$up),  ", down =", length(ad_nl_sets$down),  "\n")
cat("AD_L_vs_Healthy : up =", length(ad_h_sets$up),   ", down =", length(ad_h_sets$down),   "\n")
cat("PsO_L_vs_PsO_NL : up =", length(pso_nl_sets$up), ", down =", length(pso_nl_sets$down), "\n")
cat("PsO_L_vs_Healthy: up =", length(pso_h_sets$up),  ", down =", length(pso_h_sets$down),  "\n")

# =============================================================
# H2 / H3 完整分析
# 前置对象：res_np_u, res_mp, res_zn2, res_np_c
#           ranks_np_u, ranks_mp, ranks_zn2, ranks_np_c
#           res_ad_nl, res_ad_h, res_pso_nl, res_pso_h
# 全程使用 base R，不依赖 dplyr::select
# =============================================================

library(fgsea)

outDirGsea <- "C:/Users/DeeJo/calamine_zinc_project/results/gsea"
dir.create(outDirGsea, recursive = TRUE, showWarnings = FALSE)


# =============================================================
# 第 4 步：H2 核心 —— Directional GSEA
# =============================================================

# 统一封装：对每个 (disease_set, ranks) 组合跑 GSEA
run_reversal_gsea <- function(disease_sets, ranks, perturb_name, disease_name,
                              use_simple = FALSE) {
  if (use_simple) {
    gsea_res <- fgseaSimple(
      list(disease_up   = disease_sets$up,
           disease_down = disease_sets$down),
      ranks,
      minSize = 15, maxSize = 5000,
      nperm = 10000
    )
    gsea_res$padj <- gsea_res$padj
  } else {
    gsea_res <- fgsea(
      list(disease_up   = disease_sets$up,
           disease_down = disease_sets$down),
      ranks,
      minSize = 15, maxSize = 5000
    )
  }
  gsea_res <- as.data.frame(gsea_res)
  gsea_res$perturbation <- perturb_name
  gsea_res$disease      <- disease_name
  gsea_res
}

# --- 4.1 主分析：padj < 0.05 定义 gene set，fgseaMultilevel ---
h2_results <- rbind(
  run_reversal_gsea(ad_nl_sets, ranks_zn2,  "ZnCl2",           "AD_L_vs_AD_NL"),
  run_reversal_gsea(ad_nl_sets, ranks_mp,   "ZnO_MP",          "AD_L_vs_AD_NL"),
  run_reversal_gsea(ad_nl_sets, ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_AD_NL"),
  run_reversal_gsea(ad_h_sets,  ranks_zn2,  "ZnCl2",           "AD_L_vs_Healthy"),
  run_reversal_gsea(ad_h_sets,  ranks_mp,   "ZnO_MP",          "AD_L_vs_Healthy"),
  run_reversal_gsea(ad_h_sets,  ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_Healthy")
)

print(h2_results[, c("perturbation", "disease", "pathway", "NES", "padj", "size")])


# --- 4.2 修复 NA：对 ZnCl2 + AD_L_vs_AD_NL 用 fgseaSimple ---
h2_zn2_ad_nl_simple <- run_reversal_gsea(
  ad_nl_sets, ranks_zn2, "ZnCl2", "AD_L_vs_AD_NL_simple",
  use_simple = TRUE
)
print(h2_zn2_ad_nl_simple[, c("perturbation", "disease", "pathway",
                              "NES", "padj", "size")])


# --- 4.3 敏感性：严格 gene set (padj < 0.05 & |log2FC| > 0.5) ---
make_sets_strict <- function(df) {
  list(
    up   = df$gene[df$padj < 0.05 & df$log2FoldChange >  0.5],
    down = df$gene[df$padj < 0.05 & df$log2FoldChange < -0.5]
  )
}

ad_nl_sets_strict <- make_sets_strict(res_ad_nl)
ad_h_sets_strict  <- make_sets_strict(res_ad_h)

cat("Strict AD_L_vs_AD_NL: up =", length(ad_nl_sets_strict$up),
    ", down =", length(ad_nl_sets_strict$down), "\n")
cat("Strict AD_L_vs_Healthy: up =", length(ad_h_sets_strict$up),
    ", down =", length(ad_h_sets_strict$down), "\n")

h2_results_strict <- rbind(
  run_reversal_gsea(ad_nl_sets_strict, ranks_zn2,  "ZnCl2",           "AD_L_vs_AD_NL_strict"),
  run_reversal_gsea(ad_nl_sets_strict, ranks_mp,   "ZnO_MP",          "AD_L_vs_AD_NL_strict"),
  run_reversal_gsea(ad_nl_sets_strict, ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_AD_NL_strict"),
  run_reversal_gsea(ad_h_sets_strict,  ranks_zn2,  "ZnCl2",           "AD_L_vs_Healthy_strict"),
  run_reversal_gsea(ad_h_sets_strict,  ranks_mp,   "ZnO_MP",          "AD_L_vs_Healthy_strict"),
  run_reversal_gsea(ad_h_sets_strict,  ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_Healthy_strict")
)

print(h2_results_strict[, c("perturbation", "disease", "pathway",
                            "NES", "padj", "size")])


# =============================================================
# 第 5 步：计算 ReversalScore
# ReversalScore = (NES_AD-down − NES_AD-up) / 2
# 符号约定：>0 且 (NES_up < 0) 且 (NES_down > 0) → direction_ok = TRUE
# =============================================================

compute_reversal <- function(gsea_res) {
  out <- data.frame()
  for (pert in unique(gsea_res$perturbation)) {
    for (dis in unique(gsea_res$disease)) {
      sub <- gsea_res[gsea_res$perturbation == pert & gsea_res$disease == dis, ]
      if (nrow(sub) == 0) next
      
      nes_up    <- sub$NES[sub$pathway == "disease_up"]
      nes_down  <- sub$NES[sub$pathway == "disease_down"]
      padj_up   <- sub$padj[sub$pathway == "disease_up"]
      padj_down <- sub$padj[sub$pathway == "disease_down"]
      
      if (length(nes_up) == 0)   nes_up   <- NA
      if (length(nes_down) == 0) nes_down <- NA
      if (length(padj_up) == 0)   padj_up   <- NA
      if (length(padj_down) == 0) padj_down <- NA
      
      score <- (nes_down - nes_up) / 2
      direction_ok <- (!is.na(nes_up) & !is.na(nes_down)) &
        (nes_up < 0) & (nes_down > 0)
      
      out <- rbind(out, data.frame(
        perturbation  = pert,
        disease       = dis,
        NES_up        = nes_up,
        NES_down      = nes_down,
        padj_up       = padj_up,
        padj_down     = padj_down,
        ReversalScore = score,
        direction_ok  = direction_ok,
        stringsAsFactors = FALSE
      ))
    }
  }
  out
}

reversal_summary <- compute_reversal(h2_results)
print(reversal_summary)

reversal_summary_strict <- compute_reversal(h2_results_strict)
print(reversal_summary_strict)

length(ranks_zn2)
sum(!is.na(ranks_zn2))

cat("AD_L_vs_AD_NL up 占比:",
    round(length(ad_nl_sets$up) / length(ranks_zn2) * 100, 1), "%\n")

make_sets_tighter <- function(df, lfc_cut = 1) {
  list(
    up   = df$gene[df$padj < 0.05 & df$log2FoldChange >  lfc_cut],
    down = df$gene[df$padj < 0.05 & df$log2FoldChange < -lfc_cut]
  )
}

ad_nl_sets_tight <- make_sets_tighter(res_ad_nl, lfc_cut = 1)
ad_h_sets_tight  <- make_sets_tighter(res_ad_h,  lfc_cut = 1)

cat("Tight AD_L_vs_AD_NL: up =", length(ad_nl_sets_tight$up),
    ", down =", length(ad_nl_sets_tight$down), "\n")
cat("Tight AD_L_vs_Healthy: up =", length(ad_h_sets_tight$up),
    ", down =", length(ad_h_sets_tight$down), "\n")

h2_results_tight <- rbind(
  run_reversal_gsea(ad_nl_sets_tight, ranks_zn2,  "ZnCl2",           "AD_L_vs_AD_NL_tight"),
  run_reversal_gsea(ad_nl_sets_tight, ranks_mp,   "ZnO_MP",          "AD_L_vs_AD_NL_tight"),
  run_reversal_gsea(ad_nl_sets_tight, ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_AD_NL_tight"),
  run_reversal_gsea(ad_h_sets_tight,  ranks_zn2,  "ZnCl2",           "AD_L_vs_Healthy_tight"),
  run_reversal_gsea(ad_h_sets_tight,  ranks_mp,   "ZnO_MP",          "AD_L_vs_Healthy_tight"),
  run_reversal_gsea(ad_h_sets_tight,  ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_Healthy_tight")
)

print(h2_results_tight[, c("perturbation", "disease", "pathway",
                           "NES", "padj", "size")])

compute_reversal <- function(gsea_res) {
  out <- data.frame()
  for (pert in unique(gsea_res$perturbation)) {
    for (dis in unique(gsea_res$disease)) {
      sub <- gsea_res[gsea_res$perturbation == pert & gsea_res$disease == dis, ]
      if (nrow(sub) == 0) next
      
      nes_up    <- sub$NES[sub$pathway == "disease_up"]
      nes_down  <- sub$NES[sub$pathway == "disease_down"]
      padj_up   <- sub$padj[sub$pathway == "disease_up"]
      padj_down <- sub$padj[sub$pathway == "disease_down"]
      
      if (length(nes_up) == 0)   nes_up   <- NA
      if (length(nes_down) == 0) nes_down <- NA
      if (length(padj_up) == 0)   padj_up   <- NA
      if (length(padj_down) == 0) padj_down <- NA
      
      score <- (nes_down - nes_up) / 2
      direction_ok <- (!is.na(nes_up) & !is.na(nes_down)) &
        (nes_up < 0) & (nes_down > 0)
      
      out <- rbind(out, data.frame(
        perturbation  = pert,
        disease       = dis,
        NES_up        = nes_up,
        NES_down      = nes_down,
        padj_up       = padj_up,
        padj_down     = padj_down,
        ReversalScore = score,
        direction_ok  = direction_ok,
        stringsAsFactors = FALSE
      ))
    }
  }
  out
}

reversal_summary_tight <- compute_reversal(h2_results_tight)
print(reversal_summary_tight)

reversal_score_boot <- function(disease_sets, ranks,
                                n_boot = 500, seed = 42) {
  set.seed(seed)
  n <- length(ranks)
  gene_names <- names(ranks)
  stat_vals  <- as.numeric(ranks)
  
  scores <- numeric(n_boot)
  
  for (i in seq_len(n_boot)) {
    # 打乱 stat 值，保持基因名唯一
    rk <- setNames(sample(stat_vals, n, replace = FALSE), gene_names)
    
    g <- fgseaSimple(
      list(up = disease_sets$up, down = disease_sets$down),
      rk, minSize = 15, maxSize = 5000, nperm = 1000
    )
    g <- as.data.frame(g)
    
    nes_up   <- g$NES[g$pathway == "up"]
    nes_down <- g$NES[g$pathway == "down"]
    
    if (length(nes_up) == 0 || length(nes_down) == 0) {
      scores[i] <- NA
    } else {
      scores[i] <- (nes_down - nes_up) / 2
    }
  }
  scores
}

boot_configs <- list(
  list(sets = ad_nl_sets, ranks = ranks_mp,   tag = "AD_L_vs_AD_NL / ZnO_MP"),
  list(sets = ad_nl_sets, ranks = ranks_np_u, tag = "AD_L_vs_AD_NL / ZnO_NP_uncoated"),
  list(sets = ad_h_sets,  ranks = ranks_mp,   tag = "AD_L_vs_Healthy / ZnO_MP"),
  list(sets = ad_h_sets,  ranks = ranks_np_u, tag = "AD_L_vs_Healthy / ZnO_NP_uncoated")
)

boot_results <- data.frame()
for (cfg in boot_configs) {
  cat("Running:", cfg$tag, "\n")
  bs <- reversal_score_boot(cfg$sets, cfg$ranks, n_boot = 500)
  ci <- quantile(bs, c(0.025, 0.5, 0.975), na.rm = TRUE)
  boot_results <- rbind(boot_results, data.frame(
    comparison       = cfg$tag,
    CI_low           = ci[1],
    Median           = ci[2],
    CI_high          = ci[3],
    CI_excludes_zero = (ci[1] > 0) | (ci[3] < 0)
  ))
}
print(boot_results)

remove_mt <- function(sets) {
  list(
    up   = sets$up[!grepl("^MT[12]",   sets$up)],
    down = sets$down[!grepl("^MT[12]", sets$down)]
  )
}
ad_nl_noMT <- remove_mt(ad_nl_sets)
ad_h_noMT  <- remove_mt(ad_h_sets)

cat("NoMT AD_L_vs_AD_NL: up =", length(ad_nl_noMT$up),
    ", down =", length(ad_nl_noMT$down), "\n")
cat("NoMT AD_L_vs_Healthy: up =", length(ad_h_noMT$up),
    ", down =", length(ad_h_noMT$down), "\n")

h2_noMT <- rbind(
  run_reversal_gsea(ad_nl_noMT, ranks_mp,   "ZnO_MP",          "AD_L_vs_AD_NL_noMT"),
  run_reversal_gsea(ad_nl_noMT, ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_AD_NL_noMT"),
  run_reversal_gsea(ad_h_noMT,  ranks_mp,   "ZnO_MP",          "AD_L_vs_Healthy_noMT"),
  run_reversal_gsea(ad_h_noMT,  ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_Healthy_noMT")
)
reversal_noMT <- compute_reversal(h2_noMT)
print(reversal_noMT)

h3_results <- rbind(
  run_reversal_gsea(ad_nl_sets,  ranks_zn2, "ZnCl2", "AD_L_vs_AD_NL"),
  run_reversal_gsea(ad_h_sets,   ranks_zn2, "ZnCl2", "AD_L_vs_Healthy"),
  run_reversal_gsea(pso_nl_sets, ranks_zn2, "ZnCl2", "PsO_L_vs_PsO_NL"),
  run_reversal_gsea(pso_h_sets,  ranks_zn2, "ZnCl2", "PsO_L_vs_Healthy"),
  run_reversal_gsea(ad_nl_sets,  ranks_mp,  "ZnO_MP", "AD_L_vs_AD_NL"),
  run_reversal_gsea(pso_nl_sets, ranks_mp,  "ZnO_MP", "PsO_L_vs_PsO_NL"),
  run_reversal_gsea(ad_h_sets,   ranks_mp,  "ZnO_MP", "AD_L_vs_Healthy"),
  run_reversal_gsea(pso_h_sets,  ranks_mp,  "ZnO_MP", "PsO_L_vs_Healthy"),
  run_reversal_gsea(ad_nl_sets,  ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_AD_NL"),
  run_reversal_gsea(pso_nl_sets, ranks_np_u, "ZnO_NP_uncoated", "PsO_L_vs_PsO_NL"),
  run_reversal_gsea(ad_h_sets,   ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_Healthy"),
  run_reversal_gsea(pso_h_sets,  ranks_np_u, "ZnO_NP_uncoated", "PsO_L_vs_Healthy")
)
reversal_h3 <- compute_reversal(h3_results)
print(reversal_h3)

outDirGsea <- "C:/Users/DeeJo/calamine_zinc_project/results/gsea"
write.csv(boot_results,  file.path(outDirGsea, "H2_permutation.csv"),  row.names = FALSE)
write.csv(reversal_noMT, file.path(outDirGsea, "H2_reversal_noMT.csv"), row.names = FALSE)
write.csv(reversal_h3,   file.path(outDirGsea, "H3_AD_vs_PsO.csv"),    row.names = FALSE)

# 检查 ZnCl2 ranked list 和 AD gene set 的交集
cat("ranks_zn2 长度:", length(ranks_zn2), "\n")
cat("ad_nl_sets$up 长度:", length(ad_nl_sets$up), "\n")
cat("交集大小:", length(intersect(ad_nl_sets$up, names(ranks_zn2))), "\n")
cat("ad_nl_sets$down 交集:", length(intersect(ad_nl_sets$down, names(ranks_zn2))), "\n")

# 直接手动跑一次 fgseaSimple 看返回
test_simple <- fgseaSimple(
  list(up = ad_nl_sets$up, down = ad_nl_sets$down),
  ranks_zn2,
  minSize = 15, maxSize = 5000,
  nperm = 1000
)
print(dim(test_simple))
print(test_simple)

# =============================================================
# 用中等阈值重定义 gene set（AD 和 PsO 都用）
# =============================================================
make_sets_mid <- function(df, lfc_cut = 0.5) {
  list(
    up   = df$gene[df$padj < 0.05 & df$log2FoldChange >  lfc_cut],
    down = df$gene[df$padj < 0.05 & df$log2FoldChange < -lfc_cut]
  )
}

ad_nl_sets_mid  <- make_sets_mid(res_ad_nl)
ad_h_sets_mid   <- make_sets_mid(res_ad_h)
pso_nl_sets_mid <- make_sets_mid(res_pso_nl)
pso_h_sets_mid  <- make_sets_mid(res_pso_h)

cat("AD_L_vs_AD_NL mid : up =", length(ad_nl_sets_mid$up),
    ", down =", length(ad_nl_sets_mid$down), "\n")
cat("AD_L_vs_Healthy mid: up =", length(ad_h_sets_mid$up),
    ", down =", length(ad_h_sets_mid$down), "\n")
cat("PsO_L_vs_PsO_NL mid: up =", length(pso_nl_sets_mid$up),
    ", down =", length(pso_nl_sets_mid$down), "\n")
cat("PsO_L_vs_Healthy mid: up =", length(pso_h_sets_mid$up),
    ", down =", length(pso_h_sets_mid$down), "\n")

make_sets_strict2 <- function(df, lfc_cut = 1, padj_cut = 0.01) {
  list(
    up   = df$gene[df$padj < padj_cut & df$log2FoldChange >  lfc_cut],
    down = df$gene[df$padj < padj_cut & df$log2FoldChange < -lfc_cut]
  )
}

ad_nl_sets_s2  <- make_sets_strict2(res_ad_nl)
ad_h_sets_s2   <- make_sets_strict2(res_ad_h)
pso_nl_sets_s2 <- make_sets_strict2(res_pso_nl)
pso_h_sets_s2  <- make_sets_strict2(res_pso_h)

cat("AD_L_vs_AD_NL s2 : up =", length(ad_nl_sets_s2$up),
    ", down =", length(ad_nl_sets_s2$down), "\n")
cat("AD_L_vs_Healthy s2: up =", length(ad_h_sets_s2$up),
    ", down =", length(ad_h_sets_s2$down), "\n")
cat("PsO_L_vs_PsO_NL s2: up =", length(pso_nl_sets_s2$up),
    ", down =", length(pso_nl_sets_s2$down), "\n")
cat("PsO_L_vs_Healthy s2: up =", length(pso_h_sets_s2$up),
    ", down =", length(pso_h_sets_s2$down), "\n")

library(fgsea)

run_reversal_gsea_s2 <- function(disease_sets, ranks, perturb_name, disease_name) {
  gsea_res <- fgseaMultilevel(
    list(disease_up   = disease_sets$up,
         disease_down = disease_sets$down),
    ranks,
    minSize = 15, maxSize = 6000,
    nPermSimple = 10000
  )
  gsea_res <- as.data.frame(gsea_res)
  gsea_res$perturbation <- perturb_name
  gsea_res$disease      <- disease_name
  gsea_res
}

# ---- H2 (s2) ----
h2_results_s2 <- rbind(
  run_reversal_gsea_s2(ad_nl_sets_s2, ranks_zn2,  "ZnCl2",           "AD_L_vs_AD_NL_s2"),
  run_reversal_gsea_s2(ad_nl_sets_s2, ranks_mp,   "ZnO_MP",          "AD_L_vs_AD_NL_s2"),
  run_reversal_gsea_s2(ad_nl_sets_s2, ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_AD_NL_s2"),
  run_reversal_gsea_s2(ad_h_sets_s2,  ranks_zn2,  "ZnCl2",           "AD_L_vs_Healthy_s2"),
  run_reversal_gsea_s2(ad_h_sets_s2,  ranks_mp,   "ZnO_MP",          "AD_L_vs_Healthy_s2"),
  run_reversal_gsea_s2(ad_h_sets_s2,  ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_Healthy_s2")
)
reversal_h2_s2 <- compute_reversal(h2_results_s2)
print(reversal_h2_s2)

# ---- H3 (s2) ----
h3_results_s2 <- rbind(
  run_reversal_gsea_s2(ad_nl_sets_s2,  ranks_zn2, "ZnCl2", "AD_L_vs_AD_NL_s2"),
  run_reversal_gsea_s2(ad_h_sets_s2,   ranks_zn2, "ZnCl2", "AD_L_vs_Healthy_s2"),
  run_reversal_gsea_s2(pso_nl_sets_s2, ranks_zn2, "ZnCl2", "PsO_L_vs_PsO_NL_s2"),
  run_reversal_gsea_s2(pso_h_sets_s2,  ranks_zn2, "ZnCl2", "PsO_L_vs_Healthy_s2"),
  run_reversal_gsea_s2(ad_nl_sets_s2,  ranks_mp,  "ZnO_MP", "AD_L_vs_AD_NL_s2"),
  run_reversal_gsea_s2(pso_nl_sets_s2, ranks_mp,  "ZnO_MP", "PsO_L_vs_PsO_NL_s2"),
  run_reversal_gsea_s2(ad_h_sets_s2,   ranks_mp,  "ZnO_MP", "AD_L_vs_Healthy_s2"),
  run_reversal_gsea_s2(pso_h_sets_s2,  ranks_mp,  "ZnO_MP", "PsO_L_vs_Healthy_s2"),
  run_reversal_gsea_s2(ad_nl_sets_s2,  ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_AD_NL_s2"),
  run_reversal_gsea_s2(pso_nl_sets_s2, ranks_np_u, "ZnO_NP_uncoated", "PsO_L_vs_PsO_NL_s2"),
  run_reversal_gsea_s2(ad_h_sets_s2,   ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_Healthy_s2"),
  run_reversal_gsea_s2(pso_h_sets_s2,  ranks_np_u, "ZnO_NP_uncoated", "PsO_L_vs_Healthy_s2")
)
reversal_h3_s2 <- compute_reversal(h3_results_s2)
print(reversal_h3_s2)

# 保存
outDirGsea <- "C:/Users/DeeJo/calamine_zinc_project/results/gsea"
write.csv(h2_results_s2,  file.path(outDirGsea, "H2_gsea_s2.csv"),      row.names = FALSE)
write.csv(reversal_h2_s2, file.path(outDirGsea, "H2_reversal_s2.csv"),  row.names = FALSE)
write.csv(h3_results_s2,  file.path(outDirGsea, "H3_gsea_s2.csv"),      row.names = FALSE)
write.csv(reversal_h3_s2, file.path(outDirGsea, "H3_reversal_s2.csv"),  row.names = FALSE)

# cd /mnt/c/Users/DeeJo/calamine_zinc_project/data/raw/E-MTAB-8142
# 
# python << 'EOF'
# import anndata as ad
# 
# adata = ad.read_h5ad("submission_210120.h5ad", backed='r')
# 
# # 导出 cell metadata
# adata.obs.to_csv("cell_metadata_full.csv")
# print("metadata 已保存:", adata.obs.shape)
# 
# # 导出基因名
# import pandas as pd
# pd.Series(adata.var_names).to_csv("gene_names.txt", index=False, header=False)
# print("基因名已保存:", len(adata.var_names))
# print("前 10 个基因:", adata.var_names[:10].tolist())
# EOF

# =============================================================
# 从已保存的 differential expression 结果中构建 target_genes
# =============================================================

outDir <- "C:/Users/DeeJo/calamine_zinc_project/results/differential_expression"

# ---- 1. 读取 Zn 扰动结果（Ensembl ID）----
res_zn2 <- read.csv(file.path(outDir, "ZnCl2_vs_Control.csv"),
                    stringsAsFactors = FALSE)
res_np_u <- read.csv(file.path(outDir, "ZnO_NP_uncoated_vs_Control.csv"),
                     stringsAsFactors = FALSE)
res_mp   <- read.csv(file.path(outDir, "ZnO_MP_vs_Control.csv"),
                     stringsAsFactors = FALSE)

# ---- 2. 读取 AD 疾病结果（gene symbol）----
res_ad_nl <- read.csv(file.path(outDir, "AD_L_vs_AD_NL.csv"),
                      stringsAsFactors = FALSE)

# ---- 3. 从 GENCODE GTF 建立 Ensembl → Symbol 映射 ----
library(data.table)

gtf <- fread(
  "C:/Users/DeeJo/calamine_zinc_project/data/reference/gencode.v44.annotation.gtf.gz",
  sep = "\t", header = FALSE, skip = "##"
)
gtf_gene <- gtf[V3 == "gene"]

gene_map <- data.frame(
  ensembl = sub('.*gene_id "([^"]+)".*', '\\1', gtf_gene$V9),
  symbol  = sub('.*gene_name "([^"]+)".*', '\\1', gtf_gene$V9),
  stringsAsFactors = FALSE
)
gene_map$ensembl <- sub("\\..*$", "", gene_map$ensembl)
gene_map <- gene_map[!duplicated(gene_map$ensembl), ]

# 转换 Zn 结果的 Ensembl ID 为 symbol
to_symbol <- function(df) {
  df$ensembl <- sub("\\..*$", "", df$gene)
  df$symbol  <- gene_map$symbol[match(df$ensembl, gene_map$ensembl)]
  df[!is.na(df$symbol), ]
}

res_zn2_s   <- to_symbol(res_zn2)
res_np_u_s  <- to_symbol(res_np_u)
res_mp_s    <- to_symbol(res_mp)

# ---- 4. 构建目标基因集 ----

# Zn shared: ZnCl2 显著且 ZnO_MP 同方向的基因
zn2_up   <- res_zn2_s$symbol[res_zn2_s$padj < 0.05 & res_zn2_s$log2FoldChange >  0.5]
zn2_down <- res_zn2_s$symbol[res_zn2_s$padj < 0.05 & res_zn2_s$log2FoldChange < -0.5]
zn_shared <- unique(c(zn2_up, zn2_down))

# AD disease-reversal genes
ad_up   <- res_ad_nl$gene[res_ad_nl$padj < 0.05 & res_ad_nl$stat > 0]
ad_down <- res_ad_nl$gene[res_ad_nl$padj < 0.05 & res_ad_nl$stat < 0]

# 锌相关分子轴（从 gene_map 提取）
zinc_axis <- unique(c(
  gene_map$symbol[grepl("^MT[12]",  gene_map$symbol)],
  gene_map$symbol[grepl("^SLC30A", gene_map$symbol)],
  gene_map$symbol[grepl("^SLC39A", gene_map$symbol)],
  c("HMOX1", "NQO1", "GCLM", "GCLC", "TXNRD1", "NFE2L2", "SQSTM1")
))

# ---- 5. 合并并保存 ----
target_genes <- unique(c(zn_shared, ad_up, ad_down, zinc_axis))
target_genes <- target_genes[!is.na(target_genes) & target_genes != ""]

cat("Zn shared:", length(zn_shared), "\n")
cat("AD up:",    length(ad_up), "\n")
cat("AD down:",  length(ad_down), "\n")
cat("Zinc axis:", length(zinc_axis), "\n")
cat("总目标基因数:", length(target_genes), "\n")

writeLines(target_genes,
           "C:/Users/DeeJo/calamine_zinc_project/data/raw/E-MTAB-8142/target_genes.txt")

# ---- Zn shared (严格) ----
zn2_up   <- res_zn2_s$symbol[res_zn2_s$padj < 0.01 & res_zn2_s$log2FoldChange >  1]
zn2_down <- res_zn2_s$symbol[res_zn2_s$padj < 0.01 & res_zn2_s$log2FoldChange < -1]
zn_shared <- unique(c(zn2_up, zn2_down))

# ---- AD disease-reversal (严格) ----
ad_up   <- res_ad_nl$gene[res_ad_nl$padj < 0.01 & res_ad_nl$log2FoldChange >  1]
ad_down <- res_ad_nl$gene[res_ad_nl$padj < 0.01 & res_ad_nl$log2FoldChange < -1]

# ---- 锌相关分子轴（保留，不缩）----
zinc_axis <- unique(c(
  gene_map$symbol[grepl("^MT[12]",  gene_map$symbol)],
  gene_map$symbol[grepl("^SLC30A", gene_map$symbol)],
  gene_map$symbol[grepl("^SLC39A", gene_map$symbol)],
  c("HMOX1", "NQO1", "GCLM", "GCLC", "TXNRD1", "NFE2L2", "SQSTM1")
))

target_genes <- unique(c(zn_shared, ad_up, ad_down, zinc_axis))
target_genes <- target_genes[!is.na(target_genes) & target_genes != ""]

cat("Zn shared:", length(zn_shared), "\n")
cat("AD up:",    length(ad_up), "\n")
cat("AD down:",  length(ad_down), "\n")
cat("Zinc axis:", length(zinc_axis), "\n")
cat("总目标基因数:", length(target_genes), "\n")

writeLines(target_genes,
           "C:/Users/DeeJo/calamine_zinc_project/data/raw/E-MTAB-8142/target_genes.txt")

outDir <- "C:/Users/DeeJo/calamine_zinc_project/results/differential_expression"

# ---- 宽 gene set：padj < 0.05，不限 log2FC ----

# Zn shared wide（保留方向信息）
zn2_up_wide   <- res_zn2_s$symbol[res_zn2_s$padj < 0.05 & res_zn2_s$log2FoldChange >  0.5]
zn2_down_wide <- res_zn2_s$symbol[res_zn2_s$padj < 0.05 & res_zn2_s$log2FoldChange < -0.5]

# AD disease-reversal wide（padj < 0.05）
ad_up_wide   <- res_ad_nl$gene[res_ad_nl$padj < 0.05 & res_ad_nl$stat > 0]
ad_down_wide <- res_ad_nl$gene[res_ad_nl$padj < 0.05 & res_ad_nl$stat < 0]

# Zinc axis 保持不变
zinc_axis <- unique(c(
  gene_map$symbol[grepl("^MT[12]",  gene_map$symbol)],
  gene_map$symbol[grepl("^SLC30A", gene_map$symbol)],
  gene_map$symbol[grepl("^SLC39A", gene_map$symbol)],
  c("HMOX1", "NQO1", "GCLM", "GCLC", "TXNRD1", "NFE2L2", "SQSTM1")
))

# 保存各 gene set 为独立文件（便于后续分开使用）
setDir <- "C:/Users/DeeJo/calamine_zinc_project/results/signatures"
dir.create(setDir, recursive = TRUE, showWarnings = FALSE)

writeLines(zn2_up_wide,   file.path(setDir, "zn_up_wide.txt"))
writeLines(zn2_down_wide, file.path(setDir, "zn_down_wide.txt"))
writeLines(ad_up_wide,    file.path(setDir, "ad_up_wide.txt"))
writeLines(ad_down_wide,  file.path(setDir, "ad_down_wide.txt"))
writeLines(zinc_axis,     file.path(setDir, "zinc_axis.txt"))

# 合并的宽 gene set（用于单细胞敏感性分析）
target_genes_wide <- unique(c(zn2_up_wide, zn2_down_wide, ad_up_wide, ad_down_wide, zinc_axis))
target_genes_wide <- target_genes_wide[!is.na(target_genes_wide) & target_genes_wide != ""]

cat("Zn up wide:",   length(zn2_up_wide),   "\n")
cat("Zn down wide:", length(zn2_down_wide), "\n")
cat("AD up wide:",   length(ad_up_wide),    "\n")
cat("AD down wide:", length(ad_down_wide),  "\n")
cat("宽 gene set 总数:", length(target_genes_wide), "\n")

writeLines(target_genes_wide,
           "C:/Users/DeeJo/calamine_zinc_project/data/raw/E-MTAB-8142/target_genes_wide.txt")

# cd /mnt/c/Users/DeeJo/calamine_zinc_project/data/raw/E-MTAB-8142
# 
# python3 << 'EOF'
# import anndata as ad
# import pandas as pd
# import scipy.sparse as sp
# import scipy.io as sio
# import gc
# 
# target_genes = [line.strip() for line in open("target_genes.txt") if line.strip()]
# print("目标基因数:", len(target_genes))
# 
# adata = ad.read_h5ad("submission_210120.h5ad", backed='r')
# print("总细胞数:", adata.n_obs)
# 
# adata_genes = set(adata.var_names)
# matched = [g for g in target_genes if g in adata_genes]
# print("匹配:", len(matched))
# 
# # 分 5 批，每批约 474 个基因
# n_batch = 5
# batch_size = (len(matched) + n_batch - 1) // n_batch
# 
# for i in range(n_batch):
#   start = i * batch_size
# end   = min((i + 1) * batch_size, len(matched))
# if start >= len(matched):
#   break
# batch_genes = matched[start:end]
# print(f"\nBatch {i}: {len(batch_genes)} 个基因")
# 
# subset = adata[:, batch_genes].to_memory()
# expr = subset.X
# 
# if not sp.issparse(expr):
#   expr = sp.csr_matrix(expr)
# 
# sio.mmwrite(f"target_genes_expr_{i}.mtx", expr)
# pd.Series(batch_genes).to_csv(f"target_genes_matched_{i}.txt",
#                               index=False, header=False)
# 
# # 释放内存
# del subset, expr
# gc.collect()
# 
# # 细胞名和 metadata 只保存一次
# pd.Series(adata.obs_names).to_csv("target_genes_cells.txt",
#                                   index=False, header=False)
# adata.obs.to_csv("cell_metadata_matched.csv")
# 
# print("\n导出完成")
# EOF

library(Matrix)

setwd("C:/Users/DeeJo/calamine_zinc_project/data/raw/E-MTAB-8142")

# 读取5批表达矩阵
mats <- list()
for (i in 0:4) {
  mats[[i + 1]] <- readMM(paste0("target_genes_expr_", i, ".mtx"))
  cat("Batch", i, ":", dim(mats[[i + 1]]), "\n")
}

# 按列合并（每个批次是 细胞 × 基因子集）
expr <- do.call(cbind, mats)
cat("合并后矩阵:", dim(expr), "\n")

# 释放中间对象
rm(mats); gc()

# 读取基因名和细胞名
genes <- unlist(lapply(0:4, function(i) {
  readLines(paste0("target_genes_matched_", i, ".txt"))
}))
cells <- readLines("target_genes_cells.txt")

cat("基因数:", length(genes), " 细胞数:", length(cells), "\n")
cat("矩阵列数 == 基因数:", ncol(expr) == length(genes), "\n")
cat("矩阵行数 == 细胞数:", nrow(expr) == length(cells), "\n")

rownames(expr) <- cells
colnames(expr) <- genes

# Seurat v4 语法
expr_mat <- GetAssayData(seu, assay = "RNA", slot = "data")

# 基因集平均表达评分
avg_score <- function(genes, mat) {
  genes <- genes[genes %in% rownames(mat)]
  if (length(genes) == 0) return(rep(NA, ncol(mat)))
  Matrix::colMeans(mat[genes, , drop = FALSE])
}

seu$Zn_up_score     <- avg_score(zn_up,     expr_mat)
seu$Zn_down_score   <- avg_score(zn_down,   expr_mat)
seu$AD_up_score     <- avg_score(ad_up,     expr_mat)
seu$AD_down_score   <- avg_score(ad_down,   expr_mat)
seu$Zinc_axis_score <- avg_score(zinc_axis, expr_mat)

summary(seu$Zn_up_score)
summary(seu$AD_up_score)

saveRDS(seu,
        "C:/Users/DeeJo/calamine_zinc_project/data/processed/seu_skin_atlas.rds")

library(dplyr)

score_summary <- seu@meta.data %>%
  group_by(full_clustering, Status) %>%
  summarise(
    n_cells        = n(),
    Zn_up_mean     = mean(Zn_up_score,     na.rm = TRUE),
    Zn_down_mean   = mean(Zn_down_score,   na.rm = TRUE),
    AD_up_mean     = mean(AD_up_score,     na.rm = TRUE),
    AD_down_mean   = mean(AD_down_score,   na.rm = TRUE),
    Zinc_axis_mean = mean(Zinc_axis_score, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(score_summary,
          "C:/Users/DeeJo/calamine_zinc_project/results/single_cell/score_summary.csv",
          row.names = FALSE)

head(as.data.frame(score_summary), 30)

dir.create("C:/Users/DeeJo/calamine_zinc_project/results/single_cell",
           recursive = TRUE, showWarnings = FALSE)

write.csv(score_summary,
          "C:/Users/DeeJo/calamine_zinc_project/results/single_cell/score_summary.csv",
          row.names = FALSE)

library(dplyr)

pseudobulk <- seu@meta.data %>%
  group_by(sample_id, Status, Site, full_clustering) %>%
  summarise(
    n_cells        = n(),
    Zn_up_mean     = mean(Zn_up_score,     na.rm = TRUE),
    Zn_down_mean   = mean(Zn_down_score,   na.rm = TRUE),
    AD_up_mean     = mean(AD_up_score,     na.rm = TRUE),
    AD_down_mean   = mean(AD_down_score,   na.rm = TRUE),
    Zinc_axis_mean = mean(Zinc_axis_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 10)

write.csv(pseudobulk,
          "C:/Users/DeeJo/calamine_zinc_project/results/single_cell/pseudobulk_scores.csv",
          row.names = FALSE)

# 检查每个细胞类型的样本数
pb_count <- pseudobulk %>%
  group_by(full_clustering, Status) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Status, values_from = n_samples, values_fill = 0)
print(as.data.frame(pb_count), n = 50)

head(as.data.frame(pb_count), 50)

# 重点细胞类型
focus_cells <- c("Differentiated_KC", "Differentiated_KC*",
                 "Undifferentiated_KC", "Proliferating_KC",
                 "F1", "F2", "F3",
                 "DC1", "DC2", "moDC_1", "moDC_2",
                 "Inf_mac", "Macro_1", "Macro_2",
                 "Th", "Tc", "Treg", "Tc17_Th17")

run_test <- function(pb, cell_type, score_col) {
  sub <- pb %>% filter(full_clustering == cell_type)
  if (nrow(sub) < 6 || length(unique(sub$Status)) < 2) return(NULL)
  
  m <- lm(as.formula(paste(score_col, "~ Status")), data = sub)
  anova_res <- anova(m)
  data.frame(
    cell_type = cell_type,
    score     = score_col,
    n_samples = nrow(sub),
    p_value   = anova_res$`Pr(>F)`[1]
  )
}

score_cols <- c("Zn_up_mean", "Zn_down_mean", "AD_up_mean",
                "AD_down_mean", "Zinc_axis_mean")

results <- data.frame()
for (ct in focus_cells) {
  for (sc in score_cols) {
    r <- run_test(pseudobulk, ct, sc)
    if (!is.null(r)) results <- rbind(results, r)
  }
}

results$padj <- p.adjust(results$p_value, method = "BH")
results <- results[order(results$padj), ]

write.csv(results,
          "C:/Users/DeeJo/calamine_zinc_project/results/single_cell/disease_comparison.csv",
          row.names = FALSE)

print(head(results, 30))

library(dplyr)

# 只保留 non_lesion
pb_nl <- pseudobulk %>% filter(Site == "non_lesion")

# 检查样本数
pb_nl_count <- pb_nl %>%
  group_by(full_clustering, Status) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Status, values_from = n, values_fill = 0)

head(as.data.frame(pb_nl_count), 50)

library(dplyr)

focus_cells <- c("Undifferentiated_KC", "Differentiated_KC",
                 "Differentiated_KC*", "Proliferating_KC",
                 "F1", "F2", "F3",
                 "DC1", "DC2", "moDC_1", "moDC_2",
                 "Inf_mac", "Macro_1", "Macro_2",
                 "Th", "Tc", "Treg", "Melanocyte")

score_cols <- c("Zn_up_mean", "Zn_down_mean", "AD_up_mean",
                "AD_down_mean", "Zinc_axis_mean")

run_full <- function(pb, cell_type, score_col) {
  sub <- pb %>% filter(full_clustering == cell_type)
  # 要求每组至少 5 个样本
  n_per_group <- table(sub$Status)
  if (any(n_per_group < 5)) return(NULL)
  
  sub$Status <- factor(sub$Status, levels = c("Healthy", "Eczema", "Psoriasis"))
  m <- lm(as.formula(paste(score_col, "~ Status + Site")), data = sub)
  coefs <- summary(m)$coefficients
  
  out <- data.frame(
    cell_type = cell_type,
    score     = score_col,
    n_healthy = n_per_group["Healthy"],
    n_eczema  = n_per_group["Eczema"],
    n_pso     = n_per_group["Psoriasis"],
    contrast  = rownames(coefs)[-1],
    estimate  = coefs[-1, 1],
    p_value   = coefs[-1, 4],
    row.names = NULL
  )
  out
}

full_results <- data.frame()
for (ct in focus_cells) {
  for (sc in score_cols) {
    r <- run_full(pseudobulk, ct, sc)
    if (!is.null(r)) full_results <- rbind(full_results, r)
  }
}

full_results$padj <- p.adjust(full_results$p_value, method = "BH")
full_results <- full_results[order(full_results$padj), ]

dir.create("C:/Users/DeeJo/calamine_zinc_project/results/single_cell",
           recursive = TRUE, showWarnings = FALSE)

write.csv(full_results,
          "C:/Users/DeeJo/calamine_zinc_project/results/single_cell/disease_comparison_full.csv",
          row.names = FALSE)

head(full_results, 30)

pb_ep <- pseudobulk %>%
  filter(Status %in% c("Eczema", "Psoriasis"))

run_ep <- function(pb, cell_type, score_col) {
  sub <- pb %>% filter(full_clustering == cell_type)
  if (nrow(sub) < 10) return(NULL)
  sub$Status <- factor(sub$Status, levels = c("Eczema", "Psoriasis"))
  m <- lm(as.formula(paste(score_col, "~ Status + Site")), data = sub)
  coefs <- summary(m)$coefficients
  if (!"StatusPsoriasis" %in% rownames(coefs)) return(NULL)
  data.frame(
    cell_type = cell_type,
    score     = score_col,
    estimate  = coefs["StatusPsoriasis", 1],
    p_value   = coefs["StatusPsoriasis", 4],
    row.names = NULL
  )
}

ep_results <- data.frame()
for (ct in focus_cells) {
  for (sc in score_cols) {
    r <- run_ep(pseudobulk, ct, sc)
    if (!is.null(r)) ep_results <- rbind(ep_results, r)
  }
}
ep_results$padj <- p.adjust(ep_results$p_value, method = "BH")
ep_results <- ep_results[order(ep_results$padj), ]

write.csv(ep_results,
          "C:/Users/DeeJo/calamine_zinc_project/results/single_cell/eczema_vs_pso.csv",
          row.names = FALSE)

head(ep_results, 30)

# 已保存
# results/single_cell/score_summary.csv
# results/single_cell/pseudobulk_scores.csv
# results/single_cell/disease_comparison_full.csv
# results/single_cell/eczema_vs_pso.csv

# 还需要保存 lesion_vs_nonlesion
write.csv(lesion_results,
          "C:/Users/DeeJo/calamine_zinc_project/results/single_cell/lesion_vs_nonlesion.csv",
          row.names = FALSE)

# 保存 Seurat 对象
saveRDS(seu,
        "C:/Users/DeeJo/calamine_zinc_project/data/processed/seu_skin_atlas.rds")

# 合并细胞类型
seu$cell_group <- case_when(
  grepl("KC",        seu$full_clustering) ~ "Keratinocyte",
  grepl("^F[123]$",  seu$full_clustering) ~ "Fibroblast",
  grepl("Mac|Mono|moDC|DC|LC|Inf", seu$full_clustering) ~ "Myeloid",
  grepl("^T",        seu$full_clustering) ~ "T cell",
  grepl("VE|Pericyte", seu$full_clustering) ~ "Vascular",
  grepl("Melanocyte", seu$full_clustering) ~ "Melanocyte",
  TRUE ~ NA_character_
)

# =============================================================
# 合并细胞类型后的 pseudobulk 分析
# 从 seu@meta.data 开始，全程 base R + dplyr
# =============================================================

library(dplyr)

# ---- 1. 添加 cell_group 列 ----
seu$cell_group <- dplyr::case_when(
  grepl("KC",                       seu$full_clustering) ~ "Keratinocyte",
  grepl("^F[123]$",                 seu$full_clustering) ~ "Fibroblast",
  grepl("Mac|Mono|moDC|DC|LC|Inf",  seu$full_clustering) ~ "Myeloid",
  grepl("^Th$|^Tc$|Treg|Tc17|Tc_IL", seu$full_clustering) ~ "T_cell",
  grepl("VE|Pericyte",              seu$full_clustering) ~ "Vascular",
  grepl("Melanocyte",               seu$full_clustering) ~ "Melanocyte",
  grepl("Mast",                     seu$full_clustering) ~ "Mast_cell",
  grepl("NK|ILC",                   seu$full_clustering) ~ "NK_ILC",
  grepl("Schwann",                  seu$full_clustering) ~ "Schwann",
  TRUE ~ NA_character_
)

cat("cell_group 分布:\n")
print(table(seu$cell_group, useNA = "always"))

# ---- 2. 重建 pseudobulk ----
pseudobulk_group <- seu@meta.data %>%
  filter(!is.na(cell_group)) %>%
  group_by(sample_id, Status, Site, cell_group) %>%
  summarise(
    n_cells        = n(),
    Zn_up_mean     = mean(Zn_up_score,     na.rm = TRUE),
    Zn_down_mean   = mean(Zn_down_score,   na.rm = TRUE),
    AD_up_mean     = mean(AD_up_score,     na.rm = TRUE),
    AD_down_mean   = mean(AD_down_score,   na.rm = TRUE),
    Zinc_axis_mean = mean(Zinc_axis_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 10)

cat("\npseudobulk_group 样本数:", nrow(pseudobulk_group), "\n")

# ---- 3. 检查每个 cell_group 的样本量 ----
pb_group_count <- pseudobulk_group %>%
  group_by(cell_group, Status) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Status, values_from = n, values_fill = 0)

cat("\n各 cell_group 的样本数（all sites）:\n")
print(as.data.frame(pb_group_count))

pb_group_site <- pseudobulk_group %>%
  group_by(cell_group, Status, Site) %>%
  summarise(n = n(), .groups = "drop")

cat("\n各 cell_group 按 Site 的样本数:\n")
print(as.data.frame(pb_group_site), n = 100)

# ---- 4. lesion vs non_lesion（Eczema + Psoriasis 内部）----
score_cols <- c("Zn_up_mean", "Zn_down_mean", "AD_up_mean",
                "AD_down_mean", "Zinc_axis_mean")

cell_groups <- unique(pseudobulk_group$cell_group)
cell_groups <- cell_groups[!is.na(cell_groups)]

run_lesion_group <- function(pb, cg, sc) {
  sub <- pb %>% filter(cell_group == cg,
                       Status %in% c("Eczema", "Psoriasis"))
  if (nrow(sub) < 6) return(NULL)
  
  sub$Site <- factor(sub$Site, levels = c("non_lesion", "lesion"))
  site_n <- table(sub$Site)
  if (any(site_n < 3)) return(NULL)
  
  m <- lm(as.formula(paste(sc, "~ Status + Site")), data = sub)
  coefs <- summary(m)$coefficients
  if (!"Sitelesion" %in% rownames(coefs)) return(NULL)
  
  data.frame(
    cell_group = cg,
    score      = sc,
    n_total    = nrow(sub),
    n_nonles   = as.integer(site_n["non_lesion"]),
    n_les      = as.integer(site_n["lesion"]),
    estimate   = coefs["Sitelesion", 1],
    p_value    = coefs["Sitelesion", 4],
    row.names  = NULL
  )
}

lesion_group_results <- data.frame()
for (cg in cell_groups) {
  for (sc in score_cols) {
    r <- run_lesion_group(pseudobulk_group, cg, sc)
    if (!is.null(r)) lesion_group_results <- rbind(lesion_group_results, r)
  }
}

lesion_group_results$padj <- p.adjust(lesion_group_results$p_value, method = "BH")
lesion_group_results <- lesion_group_results[order(lesion_group_results$padj), ]

cat("\n===== lesion vs non_lesion（合并细胞类型）=====\n")
lesion_group_results

write.csv(lesion_group_results,
          "C:/Users/DeeJo/calamine_zinc_project/results/single_cell/lesion_vs_nonlesion_grouped.csv",
          row.names = FALSE)

# ---- 5. Eczema vs Psoriasis（同部位）----
run_ep_group <- function(pb, cg, sc) {
  sub <- pb %>% filter(cell_group == cg,
                       Status %in% c("Eczema", "Psoriasis"))
  if (nrow(sub) < 10) return(NULL)
  
  sub$Status <- factor(sub$Status, levels = c("Eczema", "Psoriasis"))
  m <- lm(as.formula(paste(sc, "~ Status + Site")), data = sub)
  coefs <- summary(m)$coefficients
  if (!"StatusPsoriasis" %in% rownames(coefs)) return(NULL)
  
  data.frame(
    cell_group = cg,
    score      = sc,
    n_total    = nrow(sub),
    estimate   = coefs["StatusPsoriasis", 1],
    p_value    = coefs["StatusPsoriasis", 4],
    row.names  = NULL
  )
}

ep_group_results <- data.frame()
for (cg in cell_groups) {
  for (sc in score_cols) {
    r <- run_ep_group(pseudobulk_group, cg, sc)
    if (!is.null(r)) ep_group_results <- rbind(ep_group_results, r)
  }
}

ep_group_results$padj <- p.adjust(ep_group_results$p_value, method = "BH")
ep_group_results <- ep_group_results[order(ep_group_results$padj), ]

cat("\n===== Eczema vs Psoriasis（合并细胞类型）=====\n")
ep_group_results

write.csv(ep_group_results,
          "C:/Users/DeeJo/calamine_zinc_project/results/single_cell/eczema_vs_pso_grouped.csv",
          row.names = FALSE)

# ---- 6. 保存 pseudobulk_group ----
write.csv(pseudobulk_group,
          "C:/Users/DeeJo/calamine_zinc_project/results/single_cell/pseudobulk_grouped.csv",
          row.names = FALSE)

cat("\n全部分析完成\n")

library(ggplot2)
library(dplyr)

dir.create("C:/Users/DeeJo/calamine_zinc_project/figures/main",
           recursive = TRUE, showWarnings = FALSE)

# 细胞类型 × 疾病 组成
comp <- seu@meta.data %>%
  group_by(cell_group, Status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Status) %>%
  mutate(prop = n / sum(n)) %>%
  filter(!is.na(cell_group))

p6a <- ggplot(comp, aes(x = Status, y = prop, fill = cell_group)) +
  geom_col(position = "stack", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion of cells", fill = "Cell type") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("C:/Users/DeeJo/calamine_zinc_project/figures/main/Fig6A_cell_composition.png",
       p6a, width = 6, height = 5, dpi = 300)

# 只保留主要细胞类型（细胞数 >= 1000）
main_cells <- seu@meta.data %>%
  count(cell_group) %>%
  filter(n >= 1000, !is.na(cell_group)) %>%
  pull(cell_group)

plot_df <- seu@meta.data %>%
  filter(cell_group %in% main_cells, !is.na(cell_group))

# 小提琴图：Zn_up
p6b <- ggplot(plot_df, aes(x = cell_group, y = Zn_up_score, fill = Status)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_fill_manual(values = c("Healthy" = "#4DAF4A",
                               "Eczema"  = "#E41A1C",
                               "Psoriasis" = "#377EB8")) +
  labs(x = NULL, y = "Zn_up score", fill = "Status") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("C:/Users/DeeJo/calamine_zinc_project/figures/main/Fig6B_Zn_up_by_celltype.png",
       p6b, width = 8, height = 5, dpi = 300)

library(tidyr)

plot_long <- plot_df %>%
  select(cell_group, Status, AD_up_score, AD_down_score, Zinc_axis_score) %>%
  pivot_longer(cols = c(AD_up_score, AD_down_score, Zinc_axis_score),
               names_to = "signature", values_to = "score")

p6c <- ggplot(plot_long, aes(x = cell_group, y = score, fill = Status)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2) +
  facet_wrap(~ signature, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Healthy" = "#4DAF4A",
                               "Eczema"  = "#E41A1C",
                               "Psoriasis" = "#377EB8")) +
  labs(x = NULL, y = "Module score", fill = "Status") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("C:/Users/DeeJo/calamine_zinc_project/figures/main/Fig6C_AD_signatures.png",
       p6c, width = 9, height = 10, dpi = 300)

# 只保留 Eczema + Psoriasis
plot_lvs <- seu@meta.data %>%
  filter(Status %in% c("Eczema", "Psoriasis"),
         cell_group %in% main_cells,
         !is.na(cell_group))

p6d <- ggplot(plot_lvs, aes(x = cell_group, y = AD_up_score, fill = Site)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2) +
  scale_fill_manual(values = c("non_lesion" = "#B3B3B3",
                               "lesion"     = "#D95F02")) +
  labs(x = NULL, y = "AD_up score", fill = "Site") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("C:/Users/DeeJo/calamine_zinc_project/figures/main/Fig6D_AD_up_lesion_vs_nonlesion.png",
       p6d, width = 8, height = 5, dpi = 300)

# 用 pseudobulk_group 的结果画点图
plot_df_ep <- ep_group_results %>%
  mutate(
    significant = ifelse(padj < 0.05, "padj < 0.05", "ns"),
    label       = paste(cell_group, score, sep = " | "),
    neg_log_p   = -log10(p_value)
  )

p6e <- ggplot(plot_df_ep, aes(x = estimate, y = reorder(label, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = significant, size = neg_log_p)) +
  scale_color_manual(values = c("padj < 0.05" = "#B2182B", "ns" = "grey70")) +
  labs(
    x = "Estimate (Psoriasis − Eczema)",
    y = NULL,
    color = NULL,
    size  = "-log10(P)"
  ) +
  theme_bw(base_size = 10)

ggsave(
  "C:/Users/DeeJo/calamine_zinc_project/figures/main/Fig6E_eczema_vs_pso.png",
  p6e, width = 8, height = 10, dpi = 300
)

library(ggplot2)
library(dplyr)

# 用完整结果，不设 NA
plot_df <- lesion_group_results %>%
  mutate(
    significant = ifelse(padj < 0.05, "padj < 0.05", "ns"),
    label       = paste(cell_group, score, sep = " | "),
    neg_log_p   = -log10(p_value)
  )

p6f <- ggplot(plot_df, aes(x = estimate, y = reorder(label, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = significant, size = neg_log_p)) +
  scale_color_manual(values = c("padj < 0.05" = "#B2182B", "ns" = "grey70")) +
  labs(
    x = "Estimate (lesion − non_lesion)",
    y = NULL,
    color = NULL,
    size  = "-log10(P)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "right")

ggsave(
  "C:/Users/DeeJo/calamine_zinc_project/figures/main/Fig6F_lesion_estimate.png",
  p6f, width = 8, height = 10, dpi = 300
)

list.files("C:/Users/DeeJo/calamine_zinc_project/figures/main/",
           pattern = "^Fig6", full.names = FALSE)

library(fgsea)

gmtDir <- "C:/Users/DeeJo/calamine_zinc_project/data/reference/msigdb"

hallmark_list <- gmtPathways(file.path(gmtDir, "h.all.v2026.1.Hs.symbols.gmt"))
reactome_list <- gmtPathways(file.path(gmtDir, "c2.cp.reactome.v2026.1.Hs.symbols.gmt"))
go_all_list   <- gmtPathways(file.path(gmtDir, "c5.go.v2026.1.Hs.symbols.gmt"))

cat("Hallmark:", length(hallmark_list), "\n")
cat("Reactome:", length(reactome_list), "\n")
cat("GO all:",   length(go_all_list),   "\n")

# 从 GO all 中筛出 BP 子集
go_bp_list <- go_all_list[grepl("^GOBP_", names(go_all_list))]
cat("GO BP:", length(go_bp_list), "\n")

outDir <- "C:/Users/DeeJo/calamine_zinc_project/results/differential_expression"

library(data.table)
gtf <- fread(
  "C:/Users/DeeJo/calamine_zinc_project/data/reference/gencode.v44.annotation.gtf.gz",
  sep = "\t", header = FALSE, skip = "##")
gtf_gene <- gtf[V3 == "gene"]
gene_map <- data.frame(
  ensembl = sub('.*gene_id "([^"]+)".*', '\\1', gtf_gene$V9),
  symbol  = sub('.*gene_name "([^"]+)".*', '\\1', gtf_gene$V9),
  stringsAsFactors = FALSE)
gene_map$ensembl <- sub("\\..*$", "", gene_map$ensembl)
gene_map <- gene_map[!duplicated(gene_map$ensembl), ]

make_ranks <- function(file_name) {
  df <- read.csv(file.path(outDir, file_name), stringsAsFactors = FALSE)
  df$ensembl <- sub("\\..*$", "", df$gene)
  df$symbol  <- gene_map$symbol[match(df$ensembl, gene_map$ensembl)]
  df <- df[!is.na(df$symbol) & !is.na(df$stat), ]
  df <- df[!duplicated(df$symbol), ]
  r <- setNames(df$stat, df$symbol)
  sort(r, decreasing = TRUE)
}

ranks_zn2  <- make_ranks("ZnCl2_vs_Control.csv")
ranks_mp   <- make_ranks("ZnO_MP_vs_Control.csv")
ranks_np_u <- make_ranks("ZnO_NP_uncoated_vs_Control.csv")

cat("ranks_zn2:",  length(ranks_zn2),  "\n")
cat("ranks_mp:",   length(ranks_mp),   "\n")
cat("ranks_np_u:", length(ranks_np_u), "\n")

library(fgsea)

run_pathway_gsea <- function(ranks, gene_set_list, label) {
  res <- fgsea(gene_set_list, ranks, minSize = 15, maxSize = 500)
  res <- as.data.frame(res)
  res$perturbation <- label
  res[order(res$padj), ]
}

gsea_hallmark_zn2  <- run_pathway_gsea(ranks_zn2,  hallmark_list, "ZnCl2_Hallmark")
gsea_hallmark_mp   <- run_pathway_gsea(ranks_mp,   hallmark_list, "ZnO_MP_Hallmark")
gsea_hallmark_npu  <- run_pathway_gsea(ranks_np_u, hallmark_list, "ZnO_NP_uncoated_Hallmark")

cat("\n===== ZnCl2 Hallmark top 20 =====\n")
print(gsea_hallmark_zn2[1:20, c("pathway", "NES", "padj", "size")])

cat("\n===== ZnO_MP Hallmark top 20 =====\n")
print(gsea_hallmark_mp[1:20, c("pathway", "NES", "padj", "size")])

cat("\n===== ZnO_NP_uncoated Hallmark top 20 =====\n")
print(gsea_hallmark_npu[1:20, c("pathway", "NES", "padj", "size")])

gsea_reactome_zn2  <- run_pathway_gsea(ranks_zn2,  reactome_list, "ZnCl2_Reactome")
gsea_reactome_mp   <- run_pathway_gsea(ranks_mp,   reactome_list, "ZnO_MP_Reactome")
gsea_reactome_npu  <- run_pathway_gsea(ranks_np_u, reactome_list, "ZnO_NP_uncoated_Reactome")

cat("\n===== ZnCl2 Reactome top 20 =====\n")
print(gsea_reactome_zn2[1:20, c("pathway", "NES", "padj", "size")])

cat("\n===== ZnO_MP Reactome top 20 =====\n")
print(gsea_reactome_mp[1:20, c("pathway", "NES", "padj", "size")])

gsea_gobp_zn2  <- run_pathway_gsea(ranks_zn2,  go_bp_list, "ZnCl2_GOBP")
gsea_gobp_mp   <- run_pathway_gsea(ranks_mp,   go_bp_list, "ZnO_MP_GOBP")
gsea_gobp_npu  <- run_pathway_gsea(ranks_np_u, go_bp_list, "ZnO_NP_uncoated_GOBP")

cat("\n===== ZnCl2 GO BP top 20 =====\n")
print(gsea_gobp_zn2[1:20, c("pathway", "NES", "padj", "size")])

zinc_keywords <- c("ZINC", "METAL", "METALLOTHIONEIN",
                   "OXIDATIVE", "NRF2", "INFLAMMATORY",
                   "STRESS", "DETOXIF", "HEAVY")

filter_zinc <- function(gsea_res, keywords) {
  pattern <- paste(keywords, collapse = "|")
  gsea_res[grepl(pattern, gsea_res$pathway, ignore.case = TRUE), ]
}

cat("\n===== Hallmark: zinc/metal/stress 相关 =====\n")
print(filter_zinc(gsea_hallmark_zn2, zinc_keywords)[, c("pathway","NES","padj")])
print(filter_zinc(gsea_hallmark_mp,  zinc_keywords)[, c("pathway","NES","padj")])
print(filter_zinc(gsea_hallmark_npu, zinc_keywords)[, c("pathway","NES","padj")])

cat("\n===== Reactome: zinc/metal 相关 top 20 =====\n")
print(filter_zinc(gsea_reactome_zn2, zinc_keywords)[1:20, c("pathway","NES","padj")])

outDirGsea <- "C:/Users/DeeJo/calamine_zinc_project/results/gsea"
dir.create(outDirGsea, recursive = TRUE, showWarnings = FALSE)

save_gsea <- function(gsea_res, file_name, outDirGsea) {
  df <- as.data.frame(gsea_res)
  df$leadingEdge <- NULL   # 去掉 list 列
  write.csv(df, file.path(outDirGsea, file_name), row.names = FALSE)
}

save_gsea(gsea_hallmark_zn2,  "pathway_Hallmark_ZnCl2.csv",           outDirGsea)
save_gsea(gsea_hallmark_mp,   "pathway_Hallmark_ZnO_MP.csv",          outDirGsea)
save_gsea(gsea_hallmark_npu,  "pathway_Hallmark_ZnO_NP_uncoated.csv", outDirGsea)
save_gsea(gsea_reactome_zn2,  "pathway_Reactome_ZnCl2.csv",           outDirGsea)
save_gsea(gsea_reactome_mp,   "pathway_Reactome_ZnO_MP.csv",          outDirGsea)
save_gsea(gsea_reactome_npu,  "pathway_Reactome_ZnO_NP_uncoated.csv", outDirGsea)
save_gsea(gsea_gobp_zn2,      "pathway_GOBP_ZnCl2.csv",               outDirGsea)
save_gsea(gsea_gobp_mp,       "pathway_GOBP_ZnO_MP.csv",              outDirGsea)
save_gsea(gsea_gobp_npu,      "pathway_GOBP_ZnO_NP_uncoated.csv",     outDirGsea)

cat("已保存到:", outDirGsea, "\n")

library(fgsea)

mt_genes <- c("MT1A","MT1B","MT1E","MT1F","MT1G","MT1H","MT1M","MT1X",
              "MT2A","MT3","MT4",
              "SLC30A1","SLC30A2","SLC30A3","SLC30A4","SLC30A5",
              "SLC30A6","SLC30A7","SLC30A8","SLC30A9","SLC30A10",
              "SLC39A1","SLC39A2","SLC39A3","SLC39A4","SLC39A5",
              "SLC39A6","SLC39A7","SLC39A8","SLC39A9","SLC39A10",
              "SLC39A11","SLC39A12","SLC39A13","SLC39A14")

mt_list <- list(Zinc_homeostasis = mt_genes)

for (nm in c("zn2", "mp", "np_u")) {
  rk <- get(paste0("ranks_", nm))
  res <- fgsea(mt_list, rk, minSize = 5, maxSize = 200)
  cat("\n=====", nm, "Zinc homeostasis =====\n")
  print(as.data.frame(res)[, c("pathway","NES","padj","size")])
}

library(RRHO2)

# 转成 data.frame
make_rrho_df <- function(ranks) {
  rk <- ranks[!duplicated(names(ranks))]
  rk <- sort(rk, decreasing = TRUE)
  data.frame(
    gene = names(rk),
    stat = as.numeric(rk),
    stringsAsFactors = FALSE
  )
}

# 取共有基因
common <- intersect(names(ranks_zn2), names(ranks_mp))

r1 <- make_rrho_df(ranks_zn2[common])
r2 <- make_rrho_df(ranks_mp[common])

# 确保基因顺序一致
r1 <- r1[order(r1$gene), ]
r2 <- r2[order(r2$gene), ]
stopifnot(identical(r1$gene, r2$gene))

rrho <- RRHO2_initialize(
  r1,
  r2,
  labels = c("ZnCl2", "ZnO_MP"),
  log10 = FALSE
)

library(pheatmap)

hm <- rrho$hypermat   # 195 x 195，值 = -log10(P)

pheatmap(
  hm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "red"))(50),
  main = "RRHO2: ZnCl2 vs ZnO_MP (-log10 P)",
  fontsize = 8,
  filename = "C:/Users/DeeJo/calamine_zinc_project/figures/main/Fig3_RRHO_ZnCl2_vs_ZnO_MP.png",
  width = 7, height = 6
)

cat("Up-Up 重叠基因数:", length(rrho$genelist_uu[[1]]), "\n")
cat("Down-Down 重叠基因数:", length(rrho$genelist_dd[[1]]), "\n")
cat("Up-Down 重叠基因数:", length(rrho$genelist_ud[[1]]), "\n")
cat("Down-Up 重叠基因数:", length(rrho$genelist_du[[1]]), "\n")

# 四个象限的最大 -log10 P
cat("Up-Up max -log10P:", max(hm[1:97, 1:97]), "\n")
cat("Down-Down max -log10P:", max(hm[98:195, 98:195]), "\n")
cat("Up-Down max -log10P:", max(hm[1:97, 98:195]), "\n")
cat("Down-Up max -log10P:", max(hm[98:195, 1:97]), "\n")

dim(hm)
hm[1:5, 1:5]
summary(as.vector(hm))
sum(is.na(hm))
range(hm, na.rm = TRUE)

length(rrho$genelist_uu$gene_list_overlap_uu)   # 6149
length(rrho$genelist_dd$gene_list_overlap_uu)   # 待查
length(rrho$genelist_ud$gene_list_overlap_uu)   # 待查
length(rrho$genelist_du$gene_list_overlap_uu)   # 待查

cat("Up-Up max -log10P:",    max(hm[1:97, 1:97],     na.rm = TRUE), "\n")
cat("Down-Down max -log10P:", max(hm[98:195, 98:195], na.rm = TRUE), "\n")
cat("Up-Down max -log10P:",   max(hm[1:97, 98:195],   na.rm = TRUE), "\n")
cat("Down-Up max -log10P:",   max(hm[98:195, 1:97],   na.rm = TRUE), "\n")

length(rrho$genelist_uu$gene_list_overlap_uu)   # 6149
length(rrho$genelist_dd$gene_list_overlap_dd)
length(rrho$genelist_ud$gene_list_overlap_ud)
length(rrho$genelist_du$gene_list_overlap_du)

png("C:/Users/DeeJo/calamine_zinc_project/figures/main/Fig3_RRHO_ZnCl2_vs_ZnO_MP.png",
    width = 1200, height = 1000, res = 150)
RRHO2_heatmap(rrho, labels = c("ZnCl2", "ZnO_MP"))
dev.off()

library(Biobase)

extDir <- "C:/Users/DeeJo/calamine_zinc_project/data/raw/external"

# ---- GSE13355: 从 title 提取 ----
pheno13355$group <- ifelse(grepl("_PP_", pheno13355$title), "Lesional",
                           ifelse(grepl("_PN_", pheno13355$title), "NonLesional",
                                  ifelse(grepl("_NN_", pheno13355$title), "Healthy", NA)))
pheno13355$donor <- sub(".*Individual_([0-9]+)_.*", "\\1", pheno13355$title)
table(pheno13355$group)

# ---- GSE30999: 从 characteristics_ch1.3 提取 ----
pheno30999$group <- ifelse(
  grepl("non-lesion", pheno30999$characteristics_ch1.3, ignore.case = TRUE), "NonLesional",
  ifelse(grepl("lesion", pheno30999$characteristics_ch1.3, ignore.case = TRUE), "Lesional", NA))
pheno30999$donor <- sub("subject: ", "", pheno30999$characteristics_ch1)
table(pheno30999$group)

# ---- GSE36842: 从 characteristics_ch1.1 提取 ----
pheno36842$group <- sub("group: ", "", pheno36842$characteristics_ch1.1)
pheno36842$donor <- sub("patient: ", "", pheno36842$characteristics_ch1)
table(pheno36842$group)

# ---- GSE32924: 从 characteristics_ch1.2 提取 ----
pheno32924$group <- sub("condition: ", "", pheno32924$characteristics_ch1.2)
pheno32924$donor <- sub("individual: ", "", pheno32924$characteristics_ch1.1)
table(pheno32924$group)

library(limma)

run_limma <- function(eset, pheno, case, ref, paired = FALSE, donor_col = NULL) {
  keep <- pheno$group %in% c(case, ref)
  expr <- exprs(eset)[, keep]
  ph   <- pheno[keep, ]
  
  group <- factor(ph$group, levels = c(ref, case))
  
  if (paired) {
    donor <- factor(ph[[donor_col]])
    design <- model.matrix(~ donor + group)
    fit <- lmFit(expr, design)
    fit <- eBayes(fit)
    coef_name <- "groupLesional"
  } else {
    design <- model.matrix(~ group)
    fit <- lmFit(expr, design)
    fit <- eBayes(fit)
    coef_name <- paste0("group", case)
  }
  
  res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  res$probe <- rownames(res)
  res
}

# 检查数据范围
cat("GSE36842 range:", range(exprs(eset36842), na.rm = TRUE), "\n")
cat("GSE32924 range:", range(exprs(eset32924), na.rm = TRUE), "\n")
cat("GSE13355 range:", range(exprs(eset13355), na.rm = TRUE), "\n")

# 检查平台
cat("GSE36842 platform:", annotation(eset36842), "\n")
cat("GSE32924 platform:", annotation(eset32924), "\n")

# 检查方差
cat("GSE36842 sd summary:", summary(apply(exprs(eset36842), 1, sd)), "\n")
cat("GSE32924 sd summary:", summary(apply(exprs(eset32924), 1, sd)), "\n")

sig_36842 <- run_limma(eset36842, pheno36842, "CLS", "NL", paired = FALSE)
sig_32924 <- run_limma(eset32924, pheno32924, "AL", "ANL", paired = FALSE)

cat("GSE36842 (unpaired):", sum(sig_36842$adj.P.Val < 0.05), "DEGs\n")
cat("GSE32924 (unpaired):", sum(sig_32924$adj.P.Val < 0.05), "DEGs\n")

cat("GSE36842 logFC summary:\n")
print(summary(sig_36842$logFC))

cat("\nGSE32924 logFC summary:\n")
print(summary(sig_32924$logFC))

cat("\nGSE36842 P value range:\n")
print(summary(sig_36842$P.Value))

cat("\nGSE32924 P value range:\n")
print(summary(sig_32924$P.Value))

# 检查数据是否需要重新标准化
cat("\nGSE36842 前 10 个样本的表达值（前 5 个探针）:\n")
print(exprs(eset36842)[1:5, 1:10])

cat("\nGSE32924 前 10 个样本的表达值（前 5 个探针）:\n")
print(exprs(eset32924)[1:5, 1:10])

library(fgsea)
library(AnnotationDbi)
library(hgu133plus2.db)

# ---- probe → symbol ----
map_sig <- function(res) {
  res$symbol <- mapIds(hgu133plus2.db,
                       keys = res$probe,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")
  res <- res[!is.na(res$symbol), ]
  # 多 probe → 同一 symbol：保留 |t| 最大的
  res <- res[order(-abs(res$t)), ]
  res <- res[!duplicated(res$symbol), ]
  res
}

sig_36842 <- map_sig(sig_36842)
sig_32924 <- map_sig(sig_32924)

# ---- 构建 ranked list ----
make_ranks_ext <- function(sig) {
  r <- setNames(sig$t, sig$symbol)
  sort(r, decreasing = TRUE)
}

ranks_36842 <- make_ranks_ext(sig_36842)
ranks_32924 <- make_ranks_ext(sig_32924)

cat("ranks_36842:", length(ranks_36842), "\n")
cat("ranks_32924:", length(ranks_32924), "\n")

# 读取 GSE121212 的 AD signature（已有）
outDir <- "C:/Users/DeeJo/calamine_zinc_project/results/differential_expression"
res_ad_nl <- read.csv(file.path(outDir, "AD_L_vs_AD_NL.csv"),
                      stringsAsFactors = FALSE)
res_ad_nl <- res_ad_nl[!is.na(res_ad_nl$stat), ]
res_ad_nl <- res_ad_nl[!duplicated(res_ad_nl$gene), ]

# AD-up/down gene set（用宽定义，padj < 0.05）
ad_up   <- res_ad_nl$gene[res_ad_nl$padj < 0.05 & res_ad_nl$stat > 0]
ad_down <- res_ad_nl$gene[res_ad_nl$padj < 0.05 & res_ad_nl$stat < 0]

cat("AD-up:", length(ad_up), " AD-down:", length(ad_down), "\n")

# 投射到外部 AD 排名
gsea_36842 <- fgsea(list(AD_up = ad_up, AD_down = ad_down),
                    ranks_36842, minSize = 15, maxSize = 5000)
gsea_32924 <- fgsea(list(AD_up = ad_up, AD_down = ad_down),
                    ranks_32924, minSize = 15, maxSize = 5000)

cat("\n===== GSE36842: GSE121212 AD signature =====\n")
print(as.data.frame(gsea_36842)[, c("pathway","NES","padj","size")])

cat("\n===== GSE32924: GSE121212 AD signature =====\n")
print(as.data.frame(gsea_32924)[, c("pathway","NES","padj","size")])

# 外部 AD signature 用 t 统计量直接作为 ranked list
# 用 top N 基因定义 gene set（不依赖 padj）
make_top_sets <- function(ranks, n_top = 500) {
  r_sorted <- sort(ranks, decreasing = TRUE)
  list(
    up   = names(head(r_sorted, n_top)),
    down = names(tail(r_sorted, n_top))
  )
}

sets_36842 <- make_top_sets(ranks_36842, n_top = 500)
sets_32924 <- make_top_sets(ranks_32924, n_top = 500)

cat("GSE36842 top-up:", length(sets_36842$up),
    " top-down:", length(sets_36842$down), "\n")
cat("GSE32924 top-up:", length(sets_32924$up),
    " top-down:", length(sets_32924$down), "\n")

# 投射到 ZnCl2 排名
ext_rev_36842 <- fgsea(
  list(disease_up = sets_36842$up, disease_down = sets_36842$down),
  ranks_zn2, minSize = 15, maxSize = 5000)
ext_rev_32924 <- fgsea(
  list(disease_up = sets_32924$up, disease_down = sets_32924$down),
  ranks_zn2, minSize = 15, maxSize = 5000)

cat("\n===== GSE36842 AD → ZnCl2 =====\n")
print(as.data.frame(ext_rev_36842)[, c("pathway","NES","padj","size")])

cat("\n===== GSE32924 AD → ZnCl2 =====\n")
print(as.data.frame(ext_rev_32924)[, c("pathway","NES","padj","size")])

# GSE13355 和 GSE30999 已经是正常 DEG，直接用
map_sig_pso <- function(res) {
  res$symbol <- mapIds(hgu133plus2.db,
                       keys = res$probe,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")
  res <- res[!is.na(res$symbol), ]
  res <- res[order(-abs(res$t)), ]
  res <- res[!duplicated(res$symbol), ]
  res
}

sig_13355 <- map_sig_pso(sig_13355)
sig_30999 <- map_sig_pso(sig_30999)

ranks_13355 <- sort(setNames(sig_13355$t, sig_13355$symbol), decreasing = TRUE)
ranks_30999 <- sort(setNames(sig_30999$t, sig_30999$symbol), decreasing = TRUE)

sets_13355 <- make_top_sets(ranks_13355, n_top = 500)
sets_30999 <- make_top_sets(ranks_30999, n_top = 500)

ext_rev_13355 <- fgsea(
  list(disease_up = sets_13355$up, disease_down = sets_13355$down),
  ranks_zn2, minSize = 15, maxSize = 5000)
ext_rev_30999 <- fgsea(
  list(disease_up = sets_30999$up, disease_down = sets_30999$down),
  ranks_zn2, minSize = 15, maxSize = 5000)

cat("\n===== GSE13355 PsO → ZnCl2 =====\n")
print(as.data.frame(ext_rev_13355)[, c("pathway","NES","padj","size")])

cat("\n===== GSE30999 PsO → ZnCl2 =====\n")
print(as.data.frame(ext_rev_30999)[, c("pathway","NES","padj","size")])

ext_rev_all <- rbind(
  as.data.frame(ext_rev_36842) %>% mutate(dataset = "GSE36842_AD"),
  as.data.frame(ext_rev_32924) %>% mutate(dataset = "GSE32924_AD"),
  as.data.frame(ext_rev_13355) %>% mutate(dataset = "GSE13355_PsO"),
  as.data.frame(ext_rev_30999) %>% mutate(dataset = "GSE30999_PsO")
)

valDir <- "C:/Users/DeeJo/calamine_zinc_project/results/validation"
dir.create(valDir, recursive = TRUE, showWarnings = FALSE)
write.csv(ext_rev_all,
          file.path(valDir, "external_reversal_top500.csv"), row.names = FALSE)

print(ext_rev_all[, c("dataset", "pathway", "NES", "padj")])

valDir <- "C:/Users/DeeJo/calamine_zinc_project/results/validation"
write.csv(ext_rev_all,
          file.path(valDir, "external_reversal_top500.csv"), row.names = FALSE)

sdrf <- read.delim(
  "C:/Users/DeeJo/calamine_zinc_project/data/raw/E-MTAB-5734/E-MTAB-5734.sdrf.txt",
  sep = "\t", check.names = FALSE)

# 查找可能的 batch / plate / run 字段
grep("batch|plate|run|chip|array|date|scan", names(sdrf),
     ignore.case = TRUE, value = TRUE)

sra <- read.csv(
  "C:/Users/DeeJo/calamine_zinc_project/data/raw/GSE121212/SraRunTable.csv",
  stringsAsFactors = FALSE)

# 查找可能的 batch 字段
grep("batch|plate|run|chip|date", names(sra),
     ignore.case = TRUE, value = TRUE)

# 看 Experiment 列的分布
table(sra$Experiment)

sra <- read.csv(
  "C:/Users/DeeJo/calamine_zinc_project/data/raw/GSE121212/SraRunTable.csv",
  stringsAsFactors = FALSE)

table(sra$ReleaseDate)
table(sra$create_date)
table(sra$BioProject)
table(sra$Center.Name)
table(sra$Instrument)

# 移除每个 gene set 中 |stat| 最大的 10 个基因，重跑 reversal GSEA
remove_top10 <- function(sets, ranks) {
  # 按 stat 绝对值排序，去掉前 10
  r_sorted <- ranks[order(-abs(ranks))]
  top10 <- names(r_sorted)[1:10]
  
  list(
    up   = sets$up[!sets$up %in% top10],
    down = sets$down[!sets$down %in% top10]
  )
}

# 从宽 gene set 开始
ad_nl_noTop10 <- remove_top10(ad_nl_sets, ranks_zn2)
ad_h_noTop10  <- remove_top10(ad_h_sets,  ranks_zn2)

# 重跑 reversal GSEA
h2_noTop10 <- rbind(
  run_reversal_gsea(ad_nl_noTop10, ranks_mp,   "ZnO_MP",          "AD_L_vs_AD_NL_noTop10"),
  run_reversal_gsea(ad_nl_noTop10, ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_AD_NL_noTop10"),
  run_reversal_gsea(ad_h_noTop10,  ranks_mp,   "ZnO_MP",          "AD_L_vs_Healthy_noTop10"),
  run_reversal_gsea(ad_h_noTop10,  ranks_np_u, "ZnO_NP_uncoated", "AD_L_vs_Healthy_noTop10")
)

reversal_noTop10 <- compute_reversal(h2_noTop10)
print(reversal_noTop10)

# =============================================================
# Figure 1–8 资产盘点
# =============================================================

cat("===== data/processed/ =====\n")
print(list.files("C:/Users/DeeJo/calamine_zinc_project/data/processed",
                 full.names = FALSE))

cat("\n===== results/differential_expression/ =====\n")
print(list.files("C:/Users/DeeJo/calamine_zinc_project/results/differential_expression",
                 full.names = FALSE))

cat("\n===== results/gsea/ =====\n")
print(list.files("C:/Users/DeeJo/calamine_zinc_project/results/gsea",
                 full.names = FALSE))

cat("\n===== results/signatures/ =====\n")
print(list.files("C:/Users/DeeJo/calamine_zinc_project/results/signatures",
                 full.names = FALSE))

cat("\n===== results/single_cell/ =====\n")
print(list.files("C:/Users/DeeJo/calamine_zinc_project/results/single_cell",
                 full.names = FALSE))

cat("\n===== results/validation/ =====\n")
print(list.files("C:/Users/DeeJo/calamine_zinc_project/results/validation",
                 full.names = FALSE))

cat("\n===== figures/main/ =====\n")
print(list.files("C:/Users/DeeJo/calamine_zinc_project/figures/main",
                 full.names = FALSE))

cat("\n===== data/raw/E-MTAB-5734/ =====\n")
print(list.files("C:/Users/DeeJo/calamine_zinc_project/data/raw/E-MTAB-5734",
                 full.names = FALSE))

cat("\n===== data/raw/E-MTAB-5734/salmon_quant 前 5 个 =====\n")
print(head(list.files("C:/Users/DeeJo/calamine_zinc_project/data/raw/E-MTAB-5734/salmon_quant",
                      full.names = FALSE), 5))

