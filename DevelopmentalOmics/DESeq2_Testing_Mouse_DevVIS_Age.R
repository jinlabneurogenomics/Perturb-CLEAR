setwd("/gpfs/group/jin/skim/Mouse_Reference/Allen_DevVIS_Nature/")
library(Seurat)
library(dplyr)
library(DESeq2)

sc_l23 = readRDS("Data/scrna_devvis_L23.rds")
meta = as.data.frame(sc_l23@meta.data)
colnames(sc_l23@meta.data)

table(sc_l23$orig.ident) # SeuratProject
table(sc_l23$donor_name) 
table(sc_l23$donor_name, sc_l23$Age)
table(sc_l23$subclass_label) 

## Subset Data
sc_l23_sub = subset(sc_l23, subset = Age %in% c("P9", "P10", "P11", "P14", "P28"))
sc_l23_sub$Age_group = ifelse(sc_l23_sub$Age %in% c("P9", "P10", "P11"), "P9-11", as.character(sc_l23_sub$Age))
table(sc_l23_sub$Age, sc_l23_sub$Age_group)
meta_l23_sub = as.data.frame(sc_l23_sub@meta.data)
table(sc_l23_sub$donor_name, sc_l23_sub$Age_group)

pseudo_l23 <- AggregateExpression(sc_l23_sub, group.by = c("donor_name", "Age_group", "subclass_label", "sex", "roi"), 
                                  assays = "originalexp", return.seurat = TRUE)
pseudo_l23$subclass_label = "L2/3 IT CTX Glut"
pseudo_l23$sex = "Male" # All male
meta_pseudo_l23 <- pseudo_l23@meta.data


## Run DESeq2
counts <- LayerData(pseudo_l23, assay = "originalexp", layer = "counts")
keep <- rowSums(counts >= 10) >= 2
counts <- counts[keep, , drop = FALSE]

table(rownames(meta_pseudo_l23) == colnames(counts))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta_pseudo_l23,
  design = ~ Age_group
)

dds$Age_group <- relevel(dds$Age_group, ref = "P9-11")
dds <- DESeq(dds)
resultsNames(dds)

res_1 <- results(dds, name = "Age_group_P14_vs_P9.11")
res_2 <- results(dds, name = "Age_group_P28_vs_P9.11")

dds$Age_group <- relevel(dds$Age_group, ref = "P14")
dds <- DESeq(dds)
resultsNames(dds)

res_3 <- results(dds, name="Age_group_P28_vs_P14")

res_1_df <- as.data.frame(res_1)
res_1_df$gene <- rownames(res_1_df)
res_1_df$subclass <- "L2-3 IT CTX"
res_1_df$test = "Age_group_P14_vs_P9.11"

res_2_df <- as.data.frame(res_2)
res_2_df$gene <- rownames(res_2_df)
res_2_df$subclass <- "L2-3 IT CTX"
res_2_df$test = "Age_group_P28_vs_P9.11"

res_3_df <- as.data.frame(res_3)
res_3_df$gene <- rownames(res_3_df)
res_3_df$subclass <- "L2-3 IT CTX"
res_3_df$test = "Age_group_P28_vs_P14"

write.csv(res_1_df, "DEG_Pseudobulk_DESeq2_L2-3_IT_CTX_Age_group_P14_vs_P9-11.csv")
write.csv(res_2_df, "DEG_Pseudobulk_DESeq2_L2-3_IT_CTX_Age_group_P28_vs_P9-11.csv")
write.csv(res_3_df, "DEG_Pseudobulk_DESeq2_L2-3_IT_CTX_Age_group_P28_vs_P14.csv")



### Add: Stage-specific DEG

## P9-11
meta_pseudo_l23$Test = ifelse(meta_pseudo_l23$Age_group == "P9-11", "P9-11", "Others")
table(meta_pseudo_l23$Test, meta_pseudo_l23$Age_group)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta_pseudo_l23,
  design = ~ Test
)

dds <- DESeq(dds)
resultsNames(dds)

res_1 <- results(dds, name="Test_P9.11_vs_Others")
res_1_df <- as.data.frame(res_1)
res_1_df$gene <- rownames(res_1_df)
res_1_df$subclass <- "L2-3 IT CTX"
res_1_df$test = "Test_P9.11_vs_Others"

write.csv(res_1_df, "DEG_Pseudobulk_DESeq2_L2-3_IT_CTX_Test_P9-11_vs_Others.csv")


## P14
meta_pseudo_l23$Test = ifelse(meta_pseudo_l23$Age_group == "P14", "P14", "Others")
table(meta_pseudo_l23$Test, meta_pseudo_l23$Age_group)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta_pseudo_l23,
  design = ~ Test
)

dds <- DESeq(dds)
resultsNames(dds)

res_2 <- results(dds, name="Test_P14_vs_Others")
res_2_df <- as.data.frame(res_2)
res_2_df$gene <- rownames(res_2_df)
res_2_df$subclass <- "L2-3 IT CTX"
res_2_df$test = "Test_P14_vs_Others"

write.csv(res_2_df, "DEG_Pseudobulk_DESeq2_L2-3_IT_CTX_Test_P14_vs_Others.csv")

## P28
meta_pseudo_l23$Test = ifelse(meta_pseudo_l23$Age_group == "P28", "P28", "Others")
table(meta_pseudo_l23$Test, meta_pseudo_l23$Age_group)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta_pseudo_l23,
  design = ~ Test
)

dds <- DESeq(dds)
resultsNames(dds)

res_3 <- results(dds, name="Test_P28_vs_Others")
res_3_df <- as.data.frame(res_3)
res_3_df$gene <- rownames(res_3_df)
res_3_df$subclass <- "L2-3 IT CTX"
res_3_df$test = "Test_P28_vs_Others"

write.csv(res_3_df, "DEG_Pseudobulk_DESeq2_L2-3_IT_CTX_Test_P28_vs_Others.csv")


#############################################



sc_l45 = readRDS("Data/scrna_devvis_L45.rds")
colnames(sc_l45@meta.data)
table(sc_l45$orig.ident) # SeuratProject
table(sc_l45$donor_name) 
table(sc_l45$donor_name, sc_l45$Age)
table(sc_l45$subclass_label) 

## Subset Data
sc_l45_sub = subset(sc_l45, subset = Age %in% c("P9", "P10", "P11", "P14", "P28"))
sc_l45_sub$Age_group = ifelse(sc_l45_sub$Age %in% c("P9", "P10", "P11"), "P9-11", as.character(sc_l45_sub$Age))
table(sc_l45_sub$Age_group)

pseudo_l45 <- AggregateExpression(sc_l45_sub, group.by = c("donor_name", "Age_group", "subclass_label"), 
                                  assays = "originalexp", return.seurat = TRUE)
meta_l45 <- pseudo_l45@meta.data

## Run DESeq2
counts <- LayerData(pseudo_l45, assay = "originalexp", layer = "counts")
keep <- rowSums(counts >= 10) >= 2
counts <- counts[keep, , drop = FALSE]

table(rownames(meta_l45) == colnames(pseudo_l45))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta_l45,
  design = ~ Age_group
)

dds$Age_group <- relevel(dds$Age_group, ref = "P9-11")
dds <- DESeq(dds)
resultsNames(dds)

res_1 <- results(dds, name = "Age_group_P14_vs_P9.11")
res_2 <- results(dds, name = "Age_group_P28_vs_P9.11")

dds$Age_group <- relevel(dds$Age_group, ref = "P14")
dds <- DESeq(dds)
resultsNames(dds)

res_3 <- results(dds, name="Age_group_P28_vs_P14")

res_1_df <- as.data.frame(res_1)
res_1_df$gene <- rownames(res_1_df)
res_1_df$subclass <- "L4-5 IT CTX"
res_1_df$test = "Age_group_P14_vs_P9.11"

res_2_df <- as.data.frame(res_2)
res_2_df$gene <- rownames(res_2_df)
res_2_df$subclass <- "L4-5 IT CTX"
res_2_df$test = "Age_group_P28_vs_P9.11"

res_3_df <- as.data.frame(res_3)
res_3_df$gene <- rownames(res_3_df)
res_3_df$subclass <- "L4-5 IT CTX"
res_3_df$test = "Age_group_P28_vs_P14"


write.csv(res_1_df, "DEG_Pseudobulk_DESeq2_L4-5_IT_CTX_Age_group_P14_vs_P9-11.csv")
write.csv(res_2_df, "DEG_Pseudobulk_DESeq2_L4-5_IT_CTX_Age_group_P28_vs_P9-11.csv")
write.csv(res_3_df, "DEG_Pseudobulk_DESeq2_L4-5_IT_CTX_Age_group_P28_vs_P14.csv")



### Add: Stage-specific DEG

meta_pseudo_l45 <- pseudo_l45@meta.data

## P9-11
meta_pseudo_l45$Test = ifelse(meta_pseudo_l45$Age_group == "P9-11", "P9-11", "Others")
table(meta_pseudo_l45$Test, meta_pseudo_l45$Age_group)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta_pseudo_l45,
  design = ~ Test
)

dds <- DESeq(dds)
resultsNames(dds)

res_1 <- results(dds, name="Test_P9.11_vs_Others")
res_1_df <- as.data.frame(res_1)
res_1_df$gene <- rownames(res_1_df)
res_1_df$subclass <- "L4-5 IT CTX"
res_1_df$test = "Test_P9.11_vs_Others"

write.csv(res_1_df, "DEG_Pseudobulk_DESeq2_L4-5_IT_CTX_Test_P9-11_vs_Others.csv")


## P14
meta_pseudo_l45$Test = ifelse(meta_pseudo_l45$Age_group == "P14", "P14", "Others")
table(meta_pseudo_l45$Test, meta_pseudo_l45$Age_group)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta_pseudo_l45,
  design = ~ Test
)

dds <- DESeq(dds)
resultsNames(dds)

res_2 <- results(dds, name="Test_P14_vs_Others")
res_2_df <- as.data.frame(res_2)
res_2_df$gene <- rownames(res_2_df)
res_2_df$subclass <- "L4-5 IT CTX"
res_2_df$test = "Test_P14_vs_Others"

write.csv(res_2_df, "DEG_Pseudobulk_DESeq2_L4-5_IT_CTX_Test_P14_vs_Others.csv")

## P28
meta_pseudo_l45$Test = ifelse(meta_pseudo_l45$Age_group == "P28", "P28", "Others")
table(meta_pseudo_l45$Test, meta_pseudo_l45$Age_group)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta_pseudo_l45,
  design = ~ Test
)

dds <- DESeq(dds)
resultsNames(dds)

res_3 <- results(dds, name="Test_P28_vs_Others")
res_3_df <- as.data.frame(res_3)
res_3_df$gene <- rownames(res_3_df)
res_3_df$subclass <- "L4-5 IT CTX"
res_3_df$test = "Test_P28_vs_Others"

write.csv(res_3_df, "DEG_Pseudobulk_DESeq2_L4-5_IT_CTX_Test_P28_vs_Others.csv")
