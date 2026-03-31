# rm(list = ls())
library(dplyr)
library(tidyr)
library(stringr)
library(Seurat)
library(ggplot2)
library(cowplot)
setwd("/gpfs/group/jin/skim/Mouse_Reference/Allen_DevVIS_Nature/")

p14_p9_cut = 0.3 ## fixed
p28_p14_inc_cut = 0.3
p28_p14_dec_cut = -0.3
p28_p14_con_cut = 0.1
date = 251212

# Layer 2/3 IT CTX

## Load DEG result
deg_l23_1 = read.csv("Tables/DEG_Pseudobulk_DESeq2_L2-3_IT_CTX_Age_group_P14_vs_P9-11.csv")
deg_l23_1$X <- NULL
deg_l23_1 = deg_l23_1[deg_l23_1$padj < 0.05 & !is.na(deg_l23_1$padj),]
hist(deg_l23_1$log2FoldChange, breaks = 30)

deg_l23_1 = deg_l23_1[deg_l23_1$log2FoldChange > p14_p9_cut,] ## Select genes with increased expression in P14

deg_l23_2 = read.csv("Tables/DEG_Pseudobulk_DESeq2_L2-3_IT_CTX_Age_group_P28_vs_P14.csv")
deg_l23_2$X <- NULL

deg_l23_merged = merge(deg_l23_1, deg_l23_2, by = "gene", all = F)

pdf(paste("Figures/Histogram_L2-3_P28_vs_P14_LFC_Distribution",
          "Inc", p28_p14_inc_cut,
          "Dec", p28_p14_dec_cut,
          "Con", p28_p14_con_cut,
          date, ".pdf", sep = "_"),
    width = 7, height = 5)

hist(deg_l23_merged$log2FoldChange.y,
     breaks = 15,
     xlim = c(-max(deg_l23_merged$log2FoldChange.y), 
              max(deg_l23_merged$log2FoldChange.y)))

abline(v = p28_p14_inc_cut, col = 'red', lwd = 2, lty = 'dashed')
abline(v = p28_p14_dec_cut, col = 'blue', lwd = 2, lty = 'dashed')
abline(v = p28_p14_con_cut, col = 'black', lwd = 2, lty = 'dashed')
abline(v = -p28_p14_con_cut, col = 'black', lwd = 2, lty = 'dashed')

dev.off()


nrow(deg_l23_merged[deg_l23_merged$log2FoldChange.y > p28_p14_inc_cut,]) # 121
deg_l23_merged_up = deg_l23_merged[deg_l23_merged$log2FoldChange.y > p28_p14_inc_cut,]

nrow(deg_l23_merged[deg_l23_merged$log2FoldChange.y < p28_p14_dec_cut,]) # 96
deg_l23_merged_down = deg_l23_merged[deg_l23_merged$log2FoldChange.y < p28_p14_dec_cut,]

nrow(deg_l23_merged[deg_l23_merged$log2FoldChange.y < p28_p14_con_cut & deg_l23_merged$log2FoldChange.y > -p28_p14_con_cut,]) # 122
deg_l23_merged_flat = deg_l23_merged[deg_l23_merged$log2FoldChange.y < p28_p14_con_cut & deg_l23_merged$log2FoldChange.y > -p28_p14_con_cut,]

## Save grouped genes
deg_l23_merged_all = merge(deg_l23_1, deg_l23_2, by = "gene", all = T)
deg_l23_merged_all$group = ifelse(deg_l23_merged_all$gene %in% deg_l23_merged_up$gene, "Increasing", 
                                  ifelse(deg_l23_merged_all$gene %in% deg_l23_merged_down$gene, "Decreasing", 
                                         ifelse(deg_l23_merged_all$gene %in% deg_l23_merged_flat$gene, "Converging", "None")))
table(deg_l23_merged_all$group)
# write.csv(deg_l23_merged_all, paste("Tables/TrendPlot_L2-3_Grouped_Genes", "Inc", p28_p14_inc_cut, "Dec", p28_p14_dec_cut, "Con", p28_p14_con_cut, 
#                                      date, ".csv", sep = "_"))


# ## Load Seurat Object
# sc_l23 = readRDS("Data/scrna_devvis_L23.rds")
# sc_l23 <- subset(sc_l23, subset = !(Age %in% c("E17.5", "E18.5", "P0", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P56"))) # P9-P28 only
# 
# sc_l23 <- NormalizeData(sc_l23, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# ## Pseudobulk average
# pseudo_l23 <- AverageExpression(sc_l23, group.by = c("donor_name", "Age"), assays = "originalexp", layer = "data", return.seurat = TRUE)
# pseudo_l23_mtx = as.data.frame(pseudo_l23@assays[["originalexp"]]@layers[["data"]])
# colnames(pseudo_l23_mtx) = pseudo_l23$orig.ident
# rownames(pseudo_l23_mtx) = rownames(pseudo_l23)

## Get df with psbulk averaged expr
gene_traj_df = pseudo_l23_mtx[,] %>% t() %>% data.frame()

gene_traj_df$Age = sub(".*_(P[0-9]+(-[0-9]+)?)$", "\\1", rownames(gene_traj_df))
gene_traj_df$Age_num <- dplyr::recode(gene_traj_df$Age, 
                                      "P9" = 9, "P10" = 10, "P11" = 11, "P12" = 12, "P13" = 13, "P14" = 14,
                                      "P15" = 15, "P16" = 16, "P17" = 17, "P19" = 19, "P20" = 20,
                                      "P21" = 21, "P23" = 23, "P25" = 25, "P28" = 28)
unique_ages <- sort(unique(gene_traj_df$Age_num))

colnames(gene_traj_df) <- ifelse(grepl("^[0-9]", colnames(gene_traj_df)), paste0("X", colnames(gene_traj_df)), colnames(gene_traj_df))
colnames(gene_traj_df) <- ifelse(grepl("-", colnames(gene_traj_df)), gsub("-", ".", colnames(gene_traj_df)), colnames(gene_traj_df))

long_df <- gene_traj_df %>%
  dplyr::select(Age_num,colnames(gene_traj_df)[1:32285]) %>%
  pivot_longer(cols = colnames(gene_traj_df)[1:32285],
               names_to = "gene",
               values_to = "expr")

# scaled_long <- long_df %>%
#   group_by(gene) %>%
#   mutate(expr_scaled = (expr - min(expr)) / (max(expr) - min(expr))) %>%
#   ungroup()

## Z-scoring based on P14
### Gene expression trajectories were z-scored per gene relative to the P14 expression level, 
### using the expression variance of individual genes across developmental stages.

scaled_long <- long_df %>%
  group_by(gene) %>%
  mutate(
    expr_P14 = expr[Age_num == 14][1],
    expr_sd = sd(expr),
    expr_scaled = (expr - expr_P14) / expr_sd
  ) %>%
  ungroup()

## Visualization with geom_smooth
up_gene_list = deg_l23_merged_up %>% pull(gene)

scaled_long_up = scaled_long[scaled_long$gene %in% up_gene_list,]

mean_trend_up <- scaled_long_up %>%
  group_by(Age_num) %>%
  summarise(mean_expr = mean(expr_scaled, na.rm = TRUE))

p1 <- ggplot() + 
  geom_smooth(
    data = scaled_long_up,
    aes(x = Age_num, y = expr_scaled, group = gene),
    se = FALSE,
    span = 0.3,
    linewidth = 0.5,
    color = "grey75"
  ) +
  geom_smooth(
    data = mean_trend_up,
    aes(x = Age_num, y = mean_expr),
    se = FALSE,
    span = 0.3,
    linewidth = 1.3,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = unique_ages,
    labels = paste0("P", unique_ages)
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey30"
  ) +
  geom_vline(
    xintercept = 14,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey30"
  ) +
  labs(
    title = paste0("Increasing (n=", length(up_gene_list), ")"),
    x = "Age",
    y = "Z-scored expression (relative to P14)"
  ) +
  ylim(-5, 5) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

p1

flat_gene_list = deg_l23_merged_flat %>% pull(gene)

scaled_long_flat = scaled_long[scaled_long$gene %in% flat_gene_list,]

mean_trend_flat <- scaled_long_flat %>%
  group_by(Age_num) %>%
  summarise(mean_expr = mean(expr_scaled, na.rm = TRUE))

p2 <- ggplot() + 
  geom_smooth(
    data = scaled_long_flat,
    aes(x = Age_num, y = expr_scaled, group = gene),
    se = FALSE,
    span = 0.3,
    linewidth = 0.5,
    color = "grey75"
  ) +
  geom_smooth(
    data = mean_trend_flat,
    aes(x = Age_num, y = mean_expr),
    se = FALSE,
    span = 0.3,
    linewidth = 1.3,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = unique_ages,
    labels = paste0("P", unique_ages)
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey30"
  ) +
  geom_vline(
    xintercept = 14,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey30"
  ) +
  ylim(-5, 5) +
  labs(
    title = paste0("Converging (n=", length(flat_gene_list), ")"),
    x = "Age",
    y = "Z-scored expression (relative to P14)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

p2

down_gene_list = deg_l23_merged_down %>% pull(gene)

scaled_long_down = scaled_long[scaled_long$gene %in% down_gene_list,]

mean_trend_down <- scaled_long_down %>%
  group_by(Age_num) %>%
  summarise(mean_expr = mean(expr_scaled, na.rm = TRUE))

p3 <- ggplot() + 
  geom_smooth(
    data = scaled_long_down,
    aes(x = Age_num, y = expr_scaled, group = gene),
    se = FALSE,
    span = 0.3,
    linewidth = 0.5,
    color = "grey75"
  ) +
  geom_smooth(
    data = mean_trend_down,
    aes(x = Age_num, y = mean_expr),
    se = FALSE,
    span = 0.3,
    linewidth = 1.3,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = unique_ages,
    labels = paste0("P", unique_ages)
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey30"
  ) +
  geom_vline(
    xintercept = 14,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey30"
  ) +
  ylim(-5, 5) +
  labs(
    title = paste0("Decreasing (n=", length(down_gene_list), ")"),
    x = "Age",
    y = "Z-scored expression (relative to P14)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

p3

final_plot_L23 <- ggdraw() +
  draw_label("L2/3 IT Glut", 
             x = 0.5, y = 0.98, 
             fontface = "bold", size = 14, hjust = 0.5) +
  draw_plot(
    plot_grid(p1, p2, p3, nrow = 1),
    x = 0, y = 0, width = 1, height = 0.9
  )

final_plot_L23

ggsave(final_plot_L23, filename = paste("Figures/TrendPlot_P14_norm_L2-3_IT_Glut", "Inc", p28_p14_inc_cut, "Dec", p28_p14_dec_cut, "Con", p28_p14_con_cut, 
                             date, "ylim_5.pdf", sep = "_"), width = 12, height = 4)


## Layer 4/5 IT CTX

## Load DEG result
deg_l45_1 = read.csv("Tables/DEG_Pseudobulk_DESeq2_L4-5_IT_CTX_Age_group_P14_vs_P9-11.csv")
deg_l45_1$X <- NULL
deg_l45_1 = deg_l45_1[deg_l45_1$padj < 0.05 & !is.na(deg_l45_1$padj),]
deg_l45_1 = deg_l45_1[deg_l45_1$log2FoldChange > p14_p9_cut,] ## Select genes with increased expression in P14

deg_l45_2 = read.csv("Tables/DEG_Pseudobulk_DESeq2_L4-5_IT_CTX_Age_group_P28_vs_P14.csv")
deg_l45_2$X <- NULL

deg_l45_merged = merge(deg_l45_1, deg_l45_2, by = "gene", all = F)

pdf(paste("Figures/Histogram_L4-5_P28_vs_P14_LFC_Distribution",
          "Inc", p28_p14_inc_cut,
          "Dec", p28_p14_dec_cut,
          "Con", p28_p14_con_cut,
          date, ".pdf", sep = "_"),
    width = 7, height = 5)

hist(deg_l45_merged$log2FoldChange.y,
     breaks = 15,
     xlim = c(-max(deg_l45_merged$log2FoldChange.y), 
              max(deg_l45_merged$log2FoldChange.y)))

abline(v = p28_p14_inc_cut, col = 'red', lwd = 2, lty = 'dashed')
abline(v = p28_p14_dec_cut, col = 'blue', lwd = 2, lty = 'dashed')
abline(v = p28_p14_con_cut, col = 'black', lwd = 2, lty = 'dashed')
abline(v = -p28_p14_con_cut, col = 'black', lwd = 2, lty = 'dashed')

dev.off()


nrow(deg_l45_merged[deg_l45_merged$log2FoldChange.y > p28_p14_inc_cut,]) # 103
deg_l45_merged_up = deg_l45_merged[deg_l45_merged$log2FoldChange.y > p28_p14_inc_cut,]

nrow(deg_l45_merged[deg_l45_merged$log2FoldChange.y < p28_p14_dec_cut,]) # 105
deg_l45_merged_down = deg_l45_merged[deg_l45_merged$log2FoldChange.y < p28_p14_dec_cut,]

nrow(deg_l45_merged[deg_l45_merged$log2FoldChange.y < p28_p14_con_cut & deg_l45_merged$log2FoldChange.y > -p28_p14_con_cut,]) # 133
deg_l45_merged_flat = deg_l45_merged[deg_l45_merged$log2FoldChange.y < p28_p14_con_cut & deg_l45_merged$log2FoldChange.y > -p28_p14_con_cut,]

## Save grouped genes
deg_l45_merged_all = merge(deg_l45_1, deg_l45_2, by = "gene", all = T)
deg_l45_merged_all$group = ifelse(deg_l45_merged_all$gene %in% deg_l45_merged_up$gene, "Increasing", 
                                  ifelse(deg_l45_merged_all$gene %in% deg_l45_merged_down$gene, "Decreasing", 
                                         ifelse(deg_l45_merged_all$gene %in% deg_l45_merged_flat$gene, "Converging", "None")))
table(deg_l45_merged_all$group)
# write.csv(deg_l45_merged_all, paste("Tables/TrendPlot_L4-5_Grouped_Genes", "Inc", p28_p14_inc_cut, "Dec", p28_p14_dec_cut, "Con", p28_p14_con_cut, 
#                                      date, ".csv", sep = "_"))


# ## Load Seurat Object
# sc_l45 = readRDS("Data/scrna_devvis_L45.rds")
# sc_l45 <- subset(sc_l45, subset = !(Age %in% c("E17.5", "E18", "E18.5", "P0", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P56"))) # P9-P28 only
# 
# sc_l45 <- NormalizeData(sc_l45, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# ## Pseudobulk average
# pseudo_l45 <- AverageExpression(sc_l45, group.by = c("donor_name", "Age"), assays = "originalexp", layer = "data", return.seurat = TRUE)
# pseudo_l45_mtx = as.data.frame(pseudo_l45@assays[["originalexp"]]@layers[["data"]])
# colnames(pseudo_l45_mtx) = pseudo_l45$orig.ident
# rownames(pseudo_l45_mtx) = rownames(pseudo_l45)

## Get df with psbulk averaged expr
gene_traj_df = pseudo_l45_mtx[,] %>% t() %>% data.frame()

gene_traj_df$Age = sub(".*_(P[0-9]+(-[0-9]+)?)$", "\\1", rownames(gene_traj_df))
gene_traj_df$Age_num <- dplyr::recode(gene_traj_df$Age, 
                                      "P9" = 9, "P10" = 10, "P11" = 11, "P12" = 12, "P13" = 13, "P14" = 14,
                                      "P15" = 15, "P16" = 16, "P17" = 17, "P19" = 19, "P20" = 20,
                                      "P21" = 21, "P23" = 23, "P25" = 25, "P28" = 28)
unique_ages <- sort(unique(gene_traj_df$Age_num))

colnames(gene_traj_df) <- ifelse(grepl("^[0-9]", colnames(gene_traj_df)), paste0("X", colnames(gene_traj_df)), colnames(gene_traj_df))
colnames(gene_traj_df) <- ifelse(grepl("-", colnames(gene_traj_df)), gsub("-", ".", colnames(gene_traj_df)), colnames(gene_traj_df))

long_df <- gene_traj_df %>%
  dplyr::select(Age_num,colnames(gene_traj_df)[1:32285]) %>%
  pivot_longer(cols = colnames(gene_traj_df)[1:32285],
               names_to = "gene",
               values_to = "expr")

# scaled_long <- long_df %>%
#   group_by(gene) %>%
#   mutate(expr_scaled = (expr - min(expr)) / (max(expr) - min(expr))) %>%
#   ungroup()

## Z-scoring based on P14
### Gene expression trajectories were z-scored per gene relative to the P14 expression level, 
### using the expression variance of individual genes across developmental stages.

scaled_long <- long_df %>%
  group_by(gene) %>%
  mutate(
    expr_P14 = expr[Age_num == 14][1],
    expr_sd = sd(expr),
    expr_scaled = (expr - expr_P14) / expr_sd
  ) %>%
  ungroup()

## Visualization with geom_smooth
up_gene_list = deg_l45_merged_up %>% pull(gene)

scaled_long_up = scaled_long[scaled_long$gene %in% up_gene_list,]

mean_trend_up <- scaled_long_up %>%
  group_by(Age_num) %>%
  summarise(mean_expr = mean(expr_scaled, na.rm = TRUE))

p1 <- ggplot() + 
  geom_smooth(
    data = scaled_long_up,
    aes(x = Age_num, y = expr_scaled, group = gene),
    se = FALSE,
    span = 0.3,
    linewidth = 0.5,
    color = "grey75"
  ) +
  geom_smooth(
    data = mean_trend_up,
    aes(x = Age_num, y = mean_expr),
    se = FALSE,
    span = 0.3,
    linewidth = 1.3,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = unique_ages,
    labels = paste0("P", unique_ages)
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey30"
  ) +
  geom_vline(
    xintercept = 14,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey30"
  ) +
  labs(
    title = paste0("Increasing (n=", length(up_gene_list), ")"),
    x = "Age",
    y = "Z-scored expression (relative to P14)"
  ) +
  ylim(-5, 5) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

p1

flat_gene_list = deg_l45_merged_flat %>% pull(gene)

scaled_long_flat = scaled_long[scaled_long$gene %in% flat_gene_list,]

mean_trend_flat <- scaled_long_flat %>%
  group_by(Age_num) %>%
  summarise(mean_expr = mean(expr_scaled, na.rm = TRUE))

p2 <- ggplot() + 
  geom_smooth(
    data = scaled_long_flat,
    aes(x = Age_num, y = expr_scaled, group = gene),
    se = FALSE,
    span = 0.3,
    linewidth = 0.5,
    color = "grey75"
  ) +
  geom_smooth(
    data = mean_trend_flat,
    aes(x = Age_num, y = mean_expr),
    se = FALSE,
    span = 0.3,
    linewidth = 1.3,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = unique_ages,
    labels = paste0("P", unique_ages)
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey30"
  ) +
  geom_vline(
    xintercept = 14,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey30"
  ) +
  labs(
    title = paste0("Converging (n=", length(flat_gene_list), ")"),
    x = "Age",
    y = "Z-scored expression (relative to P14)"
  ) +
  ylim(-5, 5) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

p2

down_gene_list = deg_l45_merged_down %>% pull(gene)

scaled_long_down = scaled_long[scaled_long$gene %in% down_gene_list,]

mean_trend_down <- scaled_long_down %>%
  group_by(Age_num) %>%
  summarise(mean_expr = mean(expr_scaled, na.rm = TRUE))

p3 <- ggplot() + 
  geom_smooth(
    data = scaled_long_down,
    aes(x = Age_num, y = expr_scaled, group = gene),
    se = FALSE,
    span = 0.3,
    linewidth = 0.5,
    color = "grey75"
  ) +
  geom_smooth(
    data = mean_trend_down,
    aes(x = Age_num, y = mean_expr),
    se = FALSE,
    span = 0.3,
    linewidth = 1.3,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = unique_ages,
    labels = paste0("P", unique_ages)
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey30"
  ) +
  geom_vline(
    xintercept = 14,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey30"
  ) +
  ylim(-5, 5) +
  labs(
    title = paste0("Decreasing (n=", length(down_gene_list), ")"),
    x = "Age",
    y = "Z-scored expression (relative to P14)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

p3

final_plot_l45 <- ggdraw() +
  draw_label("L4/5 IT Glut", 
             x = 0.5, y = 0.98, 
             fontface = "bold", size = 14, hjust = 0.5) +
  draw_plot(
    plot_grid(p1, p2, p3, nrow = 1),
    x = 0, y = 0, width = 1, height = 0.9
  )

final_plot_l45

ggsave(final_plot_l45, filename = paste("Figures/TrendPlot_P14_norm_L4-5_IT_Glut", "Inc", p28_p14_inc_cut, "Dec", p28_p14_dec_cut, "Con", p28_p14_con_cut, 
                                        date, "ylim_5.pdf", sep = "_"), width = 12, height = 4)
