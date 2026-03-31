rm(list = ls())
library(tidyverse)
library(readxl)
library(xml2)
library(clusterProfiler)
setwd("/gpfs/group/jin/skim/Mouse_Reference/Allen_DevVIS_Nature/")

plotHeatmap <- function(gsA, gsB, orderA, orderB, title){
  # Run fisher test
  df = run_geneEnrich(gsA, gsB)
  
  # Update the plot setting
  df = df[df$setA != 'g_bg' & df$setB != 'g_bg',] # Remove tests over grey and background
  df$fisher_padj <- p.adjust(df$fisher_p, method = 'bonferroni', n = nrow(df))
  
  df1 = df
  df1$setA_names <- factor(df1$setA, levels=orderA)
  df1$setB_names <- factor(df1$setB, levels=rev(orderB))
  
  df1$OR <- ifelse(df1$fisher_OR<1/8, 1/8, ifelse(df1$fisher_OR >=8, 8, df1$fisher_OR))
  
  ## Plot

  p <- ggplot(df1, aes(setA_names, setB_names)) + 
    geom_tile(fill='white', size=0.1, color= NA, show.legend = F) +
    geom_point(aes(fill = log2(OR)), 
               shape=22, size=3, color= 'white', show.legend = T, stroke=0.1) +
    geom_point(data= df1[df1$fisher_padj > 0.05 & df1$fisher_p <= 0.05,],
               aes(fill = log2(OR)), 
               shape=22, size=6, color= 'white', show.legend = F, stroke=0.1) +
    geom_tile(data= df1[df1$fisher_padj<= 0.05,],
              aes(fill = log2(OR)),
              size=0.1, color= NA, show.legend = F) +
    geom_point(data= df1[df1$fisher_padj<= 0.05,], 
               shape=5, size=2, color='black', show.legend = F) +
    scale_fill_gradient2(low="#6A90CA", mid = 'white', high="#CD2836", 
                         limits=c(log2(1/8),log2(8))
                         ) +
    labs(title = title, x='', y='') +
    coord_equal() +
    theme_minimal(base_size=10) + 
    theme(axis.ticks=element_blank(), 
          panel.border= element_rect(size = 1.5, fill = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x=element_text(size=10, angle=320, hjust = 0, colour="black"),
          axis.text.y = element_text(size=10, colour="black"),
          legend.key.size = unit(0.4, "cm")
          #plot.margin = unit(c(0.2, 1.5, 0.2, 0.2), "cm")
          )
  out = list(df,p)
  return(out)
}

# Function: Run Fisher
run_geneEnrich <- function(genesetA, genesetB){
  # genesetA: Gene list of interest (e.g. single cell cluster genes, eGenes)
  # genesetB: Gene list for testing (e.g. gene set for autism neurobiology)
  # genesetA_bg: Background gene list of interest
  # genesetB_bg: Background gene list of a cell type
  
  # Find the background gene list
  genesetA_bg = sort(unique(unlist(genesetA)))
  genesetB_bg = sort(unique(unlist(genesetB)))
  
  print (length(genesetA_bg))
  print (length(genesetB_bg))
  
  # Run fisher exact test
  cols = c('setA', 'setB', 
           'setA_size', 'setB_size', 
           'setA_size_background', 'setB_size_background', 
           'setA_size_hits', 
           'fisher_OR', 'fisher_p', 'overlaps')
  res = as.data.frame(matrix(nrow=0, ncol=length(cols)))
  for (genesetB_1 in names(genesetB)){
    print (genesetB_1)
    genesetB1 = genesetB[[genesetB_1]]
    genesetB_bg1 = genesetB_bg[!(genesetB_bg %in% genesetB1)]
    for (genesetA_1 in names(genesetA)){
      genesetA1 = genesetA[[genesetA_1]]
      genesetA_bg1 = genesetA_bg[!(genesetA_bg %in% genesetA1)]
      
      # create a matrix
      mat = matrix( c( sum(genesetB1 %in% genesetA1), sum(genesetA_bg1 %in% genesetB1),
                       sum(genesetB_bg1 %in% genesetA1), sum(genesetA_bg1 %in% genesetB_bg1)), 
                    ncol=2   )
      
      fis = fisher.test(mat, alternative = 'greater')
      
      out = as.data.frame(matrix(c(genesetA_1, genesetB_1,
                                   length(genesetA1), length(genesetB1), 
                                   length(genesetA_bg), length(genesetB_bg), 
                                   sum(genesetB1 %in% genesetA1), 
                                   fis$estimate, fis$p.value, 
                                   paste(sort(intersect(genesetB1, genesetA1)), collapse=', ' )) , ncol=length(cols)))
      colnames(out) = cols
      res = rbind.data.frame(res, out)
    }
  }
  
  # update the type
  res$fisher_OR <- as.numeric(res$fisher_OR)
  res$fisher_p <- as.numeric(res$fisher_p)
  
  return(res)
}

## our dataset
genes_l23 = read.csv('Tables/TrendPlot_L2-3_Grouped_Genes_Inc_0.3_Dec_-0.3_Con_0.3_251211_.csv')
genes_l45 = read.csv('Tables/TrendPlot_L4-5_Grouped_Genes_Inc_0.3_Dec_-0.3_Con_0.3_251211_.csv')

# genes_l23 = read.csv('Tables/TrendPlot_L2-3_Grouped_Genes_Inc_0.5_Dec_-0.5_Con_0.1_251211_.csv')
# genes_l45 = read.csv('Tables/TrendPlot_L4-5_Grouped_Genes_Inc_0.5_Dec_-0.5_Con_0.1_251211_.csv')

# genes_l23 = read.csv('Tables/TrendPlot_L2-3_Grouped_Genes_Inc_0.3_Dec_-0.3_Con_0.1_251211_.csv')
# genes_l45 = read.csv('Tables/TrendPlot_L4-5_Grouped_Genes_Inc_0.3_Dec_-0.3_Con_0.1_251211_.csv')


genes_list = list()
genes_list[['L2-3_Increasing']] = genes_l23[genes_l23$group == "Increasing",]$gene
genes_list[['L2-3_Stabilizing']] = genes_l23[genes_l23$group == "Converging",]$gene
genes_list[['L2-3_Decreasing']] = genes_l23[genes_l23$group == "Decreasing",]$gene

genes_list[['L4-5_Increasing']] = genes_l45[genes_l45$group == "Increasing",]$gene 
genes_list[['L4-5_Stabilizing']] = genes_l45[genes_l45$group == "Converging",]$gene
genes_list[['L4-5_Decreasing']] = genes_l45[genes_l45$group == "Decreasing",]$gene

# genes_list[['g_bg']] <- unique(c(genes_l23$gene, genes_l45$gene))
# genes_list[['g_bg']] <- unique(c(genes_l23[genes_l23$group %in% c("Increasing", "Converging", "Decreasing"),]$gene, 
#                                  genes_l45[genes_l45$group %in% c("Increasing", "Converging", "Decreasing"),]$gene))

library(tidyr)
library(dplyr)

df <- enframe(genes_list, name = "group", value = "gene") %>%
  unnest(gene)

df

write.csv(df, "GeneList_Con_Stab_Dec_260320.csv", row.names = FALSE)

## Disorder geneset
ref_ndd = read.csv("Resources/Fu_NeuroDisorder_RiskGene.csv")
colnames(ref_ndd)
homolog_db = read.csv("Resources/mart_export.txt")
table(ref_ndd$gene %in% unique(homolog_db$Gene.name))
homolog_db$gene = homolog_db$Gene.name

ref_ndd_homolog = merge(ref_ndd, homolog_db, by = "gene")
ref_ndd_homolog = ref_ndd_homolog[ref_ndd_homolog$Mouse.gene.name != "",] # 21805 -> 17867

ref_ndd_list = list()
ref_ndd_list[['ASD72']] = ref_ndd_homolog[ref_ndd_homolog$ASD72 == "TRUE",]$Mouse.gene.name
ref_ndd_list[['ASD185']] = ref_ndd_homolog[ref_ndd_homolog$ASD185 == "TRUE",]$Mouse.gene.name
ref_ndd_list[['ASD255']] = ref_ndd_homolog[ref_ndd_homolog$FDR_TADA_ASD < 0.1,]$Mouse.gene.name
ref_ndd_list[['DD309']] = ref_ndd_homolog[ref_ndd_homolog$DD309 == "TRUE",]$Mouse.gene.name
ref_ndd_list[['DD477']] = ref_ndd_homolog[ref_ndd_homolog$DD477 == "TRUE",]$Mouse.gene.name
ref_ndd_list[['NDD373']] = ref_ndd_homolog[ref_ndd_homolog$NDD373 == "TRUE",]$Mouse.gene.name
ref_ndd_list[['NDD664']] = ref_ndd_homolog[ref_ndd_homolog$NDD664 == "TRUE",]$Mouse.gene.name
ref_ndd_list[['Satterstrom102']] = ref_ndd_homolog[ref_ndd_homolog$Satterstrom102 == "TRUE",]$Mouse.gene.name
ref_ndd_list[['SCZ10']] = ref_ndd_homolog[ref_ndd_homolog$SCZ10 == "TRUE",]$Mouse.gene.name
ref_ndd_list[['SCZ244']] = ref_ndd_homolog[ref_ndd_homolog$SCZ244 == "TRUE",]$Mouse.gene.name

# ref_ndd_list[['g_bg']] = ref_ndd_homolog %>% pull(Mouse.gene.name)
# ref_ndd_list[['g_bg']] = unlist(ref_ndd_list)


## Trimming for overlapping genes

rna_universe  <- unique(c(genes_l23$gene, genes_l45$gene)) # 23158
ndd_universe  <- unique(ref_ndd_homolog$Mouse.gene.name) # 16210

background_genes <- intersect(rna_universe, ndd_universe)
length(background_genes) # 13389

genes_list_bg <- lapply(genes_list, function(g) intersect(g, background_genes))
genes_list_bg[["g_bg"]] = background_genes
ref_ndd_list_bg <- lapply(ref_ndd_list, function(g) intersect(g, background_genes))
ref_ndd_list_bg[["g_bg"]] = background_genes

## Run Enrichment Test

out_dis <- plotHeatmap(gsA = genes_list_bg,
                   gsB = ref_ndd_list_bg, 
                   orderA = names(genes_list_bg),
                   orderB = names(ref_ndd_list_bg),
                   title = '')
res_table_dis = out_dis[[1]]
out_dis[[2]]

ggsave(out_dis[[2]], filename = "Figures/FisherTest_Disorder_TrendPlot_GroupGenes_Inc_0.3_Dec_-0.3_Con_0.3_251212.pdf",
      width = 4, height = 5)
write.csv(out_dis[[1]], file = "Tables/FisherTest_Disorder_TrendPlot_GroupGenes_Inc_0.3_Dec_-0.3_Con_0.3_251212.csv",
          row.names = FALSE)


## Gene Ontology enrichment
syn_term1 = read.gmt("Resources/GOBP_SYNAPTIC_Mm/GOBP_SYNAPSE_ASSEMBLY.v2025.1.Mm.gmt")
syn_term2 = read.gmt("Resources/GOBP_SYNAPTIC_Mm/GOBP_SYNAPSE_MATURATION.v2025.1.Mm.gmt")
syn_term3 = read.gmt("Resources/GOBP_SYNAPTIC_Mm/GOBP_SYNAPSE_ORGANIZATION.v2025.1.Mm.gmt")
syn_term4 = read.gmt("Resources/GOBP_SYNAPTIC_Mm/GOBP_SYNAPSE_PRUNING.v2025.1.Mm.gmt")
syn_term5 = read.gmt("Resources/GOBP_SYNAPTIC_Mm/GOBP_SYNAPTIC_SIGNALING.v2025.1.Mm.gmt")

term_list <- list(
  GOBP_SYNAPSE_ASSEMBLY     = syn_term1$gene,
  GOBP_SYNAPSE_MATURATION   = syn_term2$gene,
  GOBP_SYNAPSE_ORGANIZATION = syn_term3$gene,
  GOBP_SYNAPSE_PRUNING      = syn_term4$gene,
  GOBP_SYNAPTIC_SIGNALING   = syn_term5$gene
)

term_list[['g_bg']] <- unique(c(genes_l23$gene, genes_l45$gene))
genes_list[['g_bg']] <- unique(c(genes_l23$gene, genes_l45$gene))

## Run Enrichment test
out_go <- plotHeatmap(gsA = genes_list,
                   gsB = term_list, 
                   orderA = names(genes_list),
                   orderB = names(term_list),
                   title = '')
res_table_go = out_go[[1]]
out_go[[2]]

ggsave(out_go[[2]], filename = "Figures/FisherTest_GO_TrendPlot_GroupGenes_Inc_0.3_Dec_-0.3_Con_0.3_251212.pdf",
       width = 4, height = 5)
write.csv(out_go[[1]], file = "Tables/FisherTest_GO_TrendPlot_GroupGenes_Inc_0.3_Dec_-0.3_Con_0.3_251212.csv",
          row.names = FALSE)
