rm(list = ls())
library(dplyr)
library(ggplot2)
library(readxl)
library("pheatmap")

# L2-3
## Load DEG table - Test
deg_l23_p9 = read.csv("Tables/DEG_Pseudobulk_DESeq2_L2-3_IT_CTX_Test_P9-11_vs_Others.csv")
deg_l23_p14 = read.csv("Tables/DEG_Pseudobulk_DESeq2_L2-3_IT_CTX_Test_P14_vs_Others.csv")
deg_l23_p28 = read.csv("Tables/DEG_Pseudobulk_DESeq2_L2-3_IT_CTX_Test_P28_vs_Others.csv")

deg_l23_p9 = deg_l23_p9[!is.na(deg_l23_p9$padj),]
deg_l23_p9$change = ifelse(deg_l23_p9$padj < 0.01 & abs(deg_l23_p9$log2FoldChange) > 0.3,
                           ifelse(deg_l23_p9$log2FoldChange > 0, "Up", "Down"), "None")
table(deg_l23_p9$change) # UP 3315

deg_l23_p14 = deg_l23_p14[!is.na(deg_l23_p14$padj),]
deg_l23_p14$change = ifelse(deg_l23_p14$padj < 0.01 & abs(deg_l23_p14$log2FoldChange) > 0.3,
                            ifelse(deg_l23_p14$log2FoldChange > 0, "Up", "Down"), "None")
table(deg_l23_p14$change) # UP 199

deg_l23_p28 = deg_l23_p28[!is.na(deg_l23_p28$padj),]
deg_l23_p28$change = ifelse(deg_l23_p28$padj < 0.01 & abs(deg_l23_p28$log2FoldChange) > 0.3,
                            ifelse(deg_l23_p28$log2FoldChange > 0, "Up", "Down"), "None")
table(deg_l23_p28$change) # UP 1212

## Make a DEG list
deg_l23_df <- bind_rows(
  deg_l23_p9   %>% filter(change == "Up") %>% slice_max(log2FoldChange, n = 150),
  deg_l23_p14  %>% filter(change == "Up") %>% slice_max(log2FoldChange, n = 150),
  deg_l23_p28  %>% filter(change == "Up") %>% slice_max(log2FoldChange, n = 150)
)
deg_l23_df$DEG_group = ifelse(deg_l23_df$test == "Test_P9.11_vs_Others", "P9-11",
                          ifelse(deg_l23_df$test == "Test_P14_vs_Others", "P14", "P28"))
table(deg_l23_df$DEG_group, deg_l23_df$test)
deg_l23_df = deg_l23_df[,c("gene", "DEG_group")]
rownames(deg_l23_df) = deg_l23_df$gene
deg_l23_df$gene <- NULL


## Load Seurat Object
sc_l23 = readRDS("Data/scrna_devvis_L23.rds") # Contains raw count

sc_l23 <- NormalizeData(sc_l23, normalization.method = "LogNormalize", scale.factor = 10000)
max(sc_l23@assays$originalexp$data) # 5597 -> 7.41296

table(sc_l23$Age) # E17.5 and E18.5 have only 1 and 4 cells each -> remove
sc_l23 = subset(sc_l23, subset = !Age %in% c("E17.5", "E18.5"))
table(sc_l23$Age) 

meta_l23 = as.data.frame(sc_l23@meta.data)
table(sc_l23$donor_name, sc_l23$Age)

## Get Average Expression
sc_l23 <- SetIdent(sc_l23, value= "donor_name")
cluster.averages <- AverageExpression(sc_l23, return.seurat = TRUE) # scale.data layer contains z-score


## Heatmap - Basic Seurat Function DoHeatmap
DoHeatmap(cluster.averages, features = rownames(deg_l23_df), label = TRUE, draw.lines = FALSE, size = 3.5, 
          angle = 0, hjust = 0.5)  + 
  guides(colour=FALSE) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =11, name = "RdBu")), 
                       limits = c(-5, 5)) + 
  theme(axis.text=element_text(size=10, color = "black"))

DoHeatmap(cluster.averages, features = rownames(deg_l23_df), label = TRUE, draw.lines = FALSE, size = 3.5, 
          angle = 0, hjust = 0.5)  + 
  guides(colour=FALSE) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =11, name = "RdBu")), 
                       limits = c(-3, 3)) + # This one gives a better visualization
  theme(axis.text=element_text(size=10, color = "black"))



## For a more flexible visualization, should transit to pheatmap

mat = cluster.averages@assays$originalexp$scale.data # data matrix

df = meta_l23[,c("donor_name", "Age")] %>% unique # sample df
df$Age <- factor(df$Age, levels = intersect(paste0("P", c(0:28, 56)), df$Age))
df <- df[order(df$Age), ]
rownames(df) = df$donor_name
df$donor_name <- NULL

deg_l23_df$DEG_group = factor(deg_l23_df$DEG_group, levels = c("P9-11", "P14", "P28"))


## Make color palette

library(viridis)

age_levels <- levels(df$Age)   
age_colors <- rev(viridis(length(age_levels), option = "magma"))
names(age_colors) <- age_levels

ann_colors = list(
  Age = age_colors
)

library(RColorBrewer)

## Make RdBu color palette
rdBu <- rev(brewer.pal(11, "RdBu"))
neg_breaks <- seq(-4, 0, length.out = 50)
pos_breaks <- seq(0, 6, length.out = 51)
my_breaks <- c(neg_breaks, pos_breaks[-1])
my_colors <- colorRampPalette(rdBu)(length(my_breaks) - 1)


## Visualize heatmap
pheatmap(
  mat[rownames(deg_l23_df), rownames(df)], 
  color = my_colors,
  breaks = my_breaks,
  show_rownames = TRUE, show_colnames = TRUE, cluster_cols = FALSE, cluster_rows = FALSE, 
  annotation_colors = ann_colors,
  annotation_col = df, annotation_row = deg_l23_df
)


## ComplexHeatmap

library(ComplexHeatmap)
library(circlize)
library(viridis)
library(RColorBrewer)

# -----------------------------
# 1) Data matrix
# -----------------------------
mat2 <- mat[rownames(deg_l23_df), rownames(df)]

# -----------------------------
# 2) Annotation colors (Age)
# -----------------------------
age_levels <- levels(df$Age)
age_colors <- rev(viridis(length(age_levels), option = "magma"))
names(age_colors) <- age_levels

column_ha <- HeatmapAnnotation(
  Age_color = df$Age,
  Age_text = anno_text(
    as.character(df$Age),
    rot = 90,
    gp = gpar(fontsize = 8)  
  ),
  col = list(Age_color = age_colors),
  annotation_height = unit.c(unit(3, "mm"), unit(6, "mm")),
  annotation_legend_param = list(
    Age_color = list(title = "Age", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 10))
  )
)

# -----------------------------
# Row annotation (DEG_group)
# -----------------------------
deg_levels <- levels(deg_l23_df$DEG_group)
deg_colors <- rev(viridis(length(deg_levels), option = "viridis"))
names(deg_colors) <- deg_levels

row_ha <- rowAnnotation(
  DEG_group = deg_l23_df$DEG_group,
  col = list(DEG_group = deg_colors),
  annotation_legend_param = list(
    DEG_group = list(title = "DEG group",
                     labels_gp = gpar(fontsize = 10),
                     title_gp = gpar(fontsize = 10))
  ),
  gp = gpar(fontsize = 10),     # ← annotation 내부 텍스트 폰트 10
  annotation_width = unit(1, "mm")
)

# -----------------------------
# 3) Heatmap color palette (RdBu)
# -----------------------------
rdBu <- rev(brewer.pal(11, "RdBu"))
my_colors <- colorRampPalette(rdBu)(length(my_breaks))

col_fun <- colorRamp2(my_breaks, my_colors)

# -----------------------------
# 4) Draw heatmap
# -----------------------------
col_dend <- hclust(dist(t(mat2)))
col_dend_rev <- as.dendrogram(col_dend)

pdf(file="Heatmap_L2-3_IT_DEG_by_group_251126.pdf")

Heatmap(
  mat2,
  name = "expr",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_row_slices = TRUE,
  # cluster_columns = col_dend_rev,
  cluster_columns = FALSE,
  row_split = deg_l23_df$DEG_group,
  
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 8),
  
  top_annotation = column_ha,
  left_annotation = row_ha,
  
  heatmap_legend_param = list(
    title = "Expression (scaled)",
    at = c(-6, -3, 0, 3, 6),
    labels_gp = gpar(fontsize = 10),   # legend label size
    title_gp = gpar(fontsize = 10)     # legend title size
  )
)

dev.off()




# L4-5
## Load DEG table - Test
deg_l45_p9 = read.csv("Tables/DEG_Pseudobulk_DESeq2_L4-5_IT_CTX_Test_P9-11_vs_Others.csv")
deg_l45_p14 = read.csv("Tables/DEG_Pseudobulk_DESeq2_L4-5_IT_CTX_Test_P14_vs_Others.csv")
deg_l45_p28 = read.csv("Tables/DEG_Pseudobulk_DESeq2_L4-5_IT_CTX_Test_P28_vs_Others.csv")

deg_l45_p9 = deg_l45_p9[!is.na(deg_l45_p9$padj),]
deg_l45_p9$change = ifelse(deg_l45_p9$padj < 0.01 & abs(deg_l45_p9$log2FoldChange) > 0.3,
                           ifelse(deg_l45_p9$log2FoldChange > 0, "Up", "Down"), "None")
table(deg_l45_p9$change) # UP 3480

deg_l45_p14 = deg_l45_p14[!is.na(deg_l45_p14$padj),]
deg_l45_p14$change = ifelse(deg_l45_p14$padj < 0.01 & abs(deg_l45_p14$log2FoldChange) > 0.3,
                            ifelse(deg_l45_p14$log2FoldChange > 0, "Up", "Down"), "None")
table(deg_l45_p14$change) # UP 177

deg_l45_p28 = deg_l45_p28[!is.na(deg_l45_p28$padj),]
deg_l45_p28$change = ifelse(deg_l45_p28$padj < 0.01 & abs(deg_l45_p28$log2FoldChange) > 0.3,
                            ifelse(deg_l45_p28$log2FoldChange > 0, "Up", "Down"), "None")
table(deg_l45_p28$change) # UP 1099

## Make a DEG list
deg_l45_df <- bind_rows(
  deg_l45_p9   %>% filter(change == "Up") %>% slice_max(log2FoldChange, n = 150),
  deg_l45_p14  %>% filter(change == "Up") %>% slice_max(log2FoldChange, n = 150),
  deg_l45_p28  %>% filter(change == "Up") %>% slice_max(log2FoldChange, n = 150)
)
deg_l45_df$DEG_group = ifelse(deg_l45_df$test == "Test_P9.11_vs_Others", "P9-11",
                              ifelse(deg_l45_df$test == "Test_P14_vs_Others", "P14", "P28"))
table(deg_l45_df$DEG_group, deg_l45_df$test)
deg_l45_df = deg_l45_df[,c("gene", "DEG_group")]
rownames(deg_l45_df) = deg_l45_df$gene
deg_l45_df$gene <- NULL


## Load Seurat Object
sc_l45 = readRDS("Data/scrna_devvis_L45.rds") # Contains raw count

sc_l45 <- NormalizeData(sc_l45, normalization.method = "LogNormalize", scale.factor = 10000)
max(sc_l45@assays$originalexp$data) # 6.88049

table(sc_l45$Age) # Embryonic stage samples have less than 100 cells -> exclude
sc_l45 = subset(sc_l45, subset = !Age %in% c("E17.5", "E18", "E18.5"))
table(sc_l45$Age) 

meta_l45 = as.data.frame(sc_l45@meta.data)
table(sc_l45$donor_name, sc_l45$Age)

## Get Average Expression
sc_l45 <- SetIdent(sc_l45, value= "donor_name")
cluster.averages <- AverageExpression(sc_l45, return.seurat = TRUE) # scale.data layer contains z-score


## Heatmap - Basic Seurat Function DoHeatmap
DoHeatmap(cluster.averages, features = rownames(deg_l45_df), label = TRUE, draw.lines = FALSE, size = 3.5, 
          angle = 0, hjust = 0.5)  + 
  guides(colour=FALSE) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =11, name = "RdBu")), 
                       limits = c(-5, 5)) + 
  theme(axis.text=element_text(size=10, color = "black"))

DoHeatmap(cluster.averages, features = rownames(deg_l45_df), label = TRUE, draw.lines = FALSE, size = 3.5, 
          angle = 0, hjust = 0.5)  + 
  guides(colour=FALSE) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =11, name = "RdBu")), 
                       limits = c(-3, 3)) + # This one gives a better visualization
  theme(axis.text=element_text(size=10, color = "black"))


## For a more flexible visualization, should transit to pheatmap

mat = cluster.averages@assays$originalexp$scale.data # data matrix

df = meta_l45[,c("donor_name", "Age")] %>% unique # sample df
df$Age <- factor(df$Age, levels = intersect(paste0("P", c(0:28, 56)), df$Age))
df <- df[order(df$Age), ]
rownames(df) = df$donor_name
df$donor_name <- NULL

deg_l45_df$DEG_group = factor(deg_l45_df$DEG_group, levels = c("P9-11", "P14", "P28"))


## Make color palette

library(viridis)

age_levels <- levels(df$Age)   
age_colors <- rev(viridis(length(age_levels), option = "magma"))
names(age_colors) <- age_levels

ann_colors = list(
  Age = age_colors
)

library(RColorBrewer)

## Make RdBu color palette
rdBu <- rev(brewer.pal(11, "RdBu"))
neg_breaks <- seq(-4, 0, length.out = 50)
pos_breaks <- seq(0, 6, length.out = 51)
my_breaks <- c(neg_breaks, pos_breaks[-1])
my_colors <- colorRampPalette(rdBu)(length(my_breaks) - 1)


## Visualize heatmap
pheatmap(
  mat[rownames(deg_l45_df), rownames(df)], 
  color = my_colors,
  breaks = my_breaks,
  show_rownames = TRUE, show_colnames = TRUE, cluster_cols = FALSE, cluster_rows = FALSE, 
  annotation_colors = ann_colors,
  annotation_col = df, annotation_row = deg_l45_df
)


## ComplexHeatmap

library(ComplexHeatmap)
library(circlize)
library(viridis)
library(RColorBrewer)

# -----------------------------
# 1) Data matrix
# -----------------------------
mat2 <- mat[rownames(deg_l45_df), rownames(df)]

# -----------------------------
# 2) Annotation colors (Age)
# -----------------------------
age_levels <- levels(df$Age)
age_colors <- rev(viridis(length(age_levels), option = "magma"))
names(age_colors) <- age_levels

column_ha <- HeatmapAnnotation(
  Age_color = df$Age,
  Age_text = anno_text(
    as.character(df$Age),
    rot = 90,
    gp = gpar(fontsize = 8)  
  ),
  col = list(Age_color = age_colors),
  annotation_height = unit.c(unit(3, "mm"), unit(6, "mm")),
  annotation_legend_param = list(
    Age_color = list(title = "Age", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 10))
  )
)

# -----------------------------
# Row annotation (DEG_group)
# -----------------------------
deg_levels <- levels(deg_l45_df$DEG_group)
deg_colors <- rev(viridis(length(deg_levels), option = "viridis"))
names(deg_colors) <- deg_levels

row_ha <- rowAnnotation(
  DEG_group = deg_l45_df$DEG_group,
  col = list(DEG_group = deg_colors),
  annotation_legend_param = list(
    DEG_group = list(title = "DEG group",
                     labels_gp = gpar(fontsize = 10),
                     title_gp = gpar(fontsize = 10))
  ),
  gp = gpar(fontsize = 10),     # ← annotation 내부 텍스트 폰트 10
  annotation_width = unit(1, "mm")
)

# -----------------------------
# 3) Heatmap color palette (RdBu)
# -----------------------------
rdBu <- rev(brewer.pal(11, "RdBu"))
my_colors <- colorRampPalette(rdBu)(length(my_breaks))

col_fun <- colorRamp2(my_breaks, my_colors)

# -----------------------------
# 4) Draw heatmap
# -----------------------------
col_dend <- hclust(dist(t(mat2)))
col_dend_rev <- as.dendrogram(col_dend)

pdf(file="Heatmap_L4-5_IT_DEG_by_group_251126.pdf")

Heatmap(
  mat2,
  name = "expr",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_row_slices = TRUE,
  # cluster_columns = col_dend_rev,
  cluster_columns = FALSE,
  row_split = deg_l45_df$DEG_group,
  
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 8),
  
  top_annotation = column_ha,
  left_annotation = row_ha,
  
  heatmap_legend_param = list(
    title = "Expression (scaled)",
    at = c(-6, -3, 0, 3, 6),
    labels_gp = gpar(fontsize = 10),   # legend label size
    title_gp = gpar(fontsize = 10)     # legend title size
  )
)

dev.off()

