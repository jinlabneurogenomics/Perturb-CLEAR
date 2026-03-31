rm(list = ls())
library(gprofiler2)
library(dplyr)
library(purrr)

### Cutoff: 0.5 / 0.5 / 0.1
genes_l23 = read.csv('Tables/TrendPlot_L2-3_Grouped_Genes_Inc_0.5_Dec_-0.5_Con_0.1_251211_.csv')
genes_l45 = read.csv('Tables/TrendPlot_L4-5_Grouped_Genes_Inc_0.5_Dec_-0.5_Con_0.1_251211_.csv')

genes_list = list()
genes_list[['L2-3_Increasing']] = genes_l23[genes_l23$group == "Increasing",]$gene
genes_list[['L2-3_Converging']] = genes_l23[genes_l23$group == "Converging",]$gene
genes_list[['L2-3_Decreasing']] = genes_l23[genes_l23$group == "Decreasing",]$gene

genes_list[['L4-5_Increasing']] = genes_l45[genes_l45$group == "Increasing",]$gene 
genes_list[['L4-5_Converging']] = genes_l45[genes_l45$group == "Converging",]$gene
genes_list[['L4-5_Decreasing']] = genes_l45[genes_l45$group == "Decreasing",]$gene


gostres = gost(
  genes_list[['L2-3_Increasing']],
  organism = "mmusculus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = TRUE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  sources = "GO",
  as_short_link = FALSE,
  highlight = FALSE
)

gostres_df = as.data.frame(gostres[["result"]])
table(gostres_df$significant)

## Run gost for all list elements

gost_results <- lapply(names(genes_list), function(name) {
  message("Running g:Profiler for ", name)
  
  res <- gost(
    genes_list[[name]],
    organism = "mmusculus",
    ordered_query = FALSE,
    multi_query = FALSE,
    significant = FALSE,
    exclude_iea = TRUE,
    measure_underrepresentation = FALSE,
    evcodes = TRUE,
    user_threshold = 0.05,
    correction_method = "g_SCS",
    sources = "GO",
    as_short_link = FALSE,
    highlight = FALSE
  )
  
  if (!is.null(res$result)) {
    df <- as.data.frame(res$result)
    df$gene_set <- name
    return(df)
  } else {
    return(NULL)
  }
})

gost_results_df <- dplyr::bind_rows(gost_results_with_group)
gost_results_df_clean <- gost_results_df %>%
  mutate(across(
    where(is.list),
    ~ sapply(.x, function(y) paste(y, collapse = ";"))
  ))

# write.csv(gost_results_df_clean, file = 'Tables/gprofileR_TrendPlot_Grouped_Genes_Inc_0.5_Dec_-0.5_Con_0.1_251211_.csv', row.names = FALSE)



### Cutoff: 0.3 / 0.3 / 0.3
genes_l23 = read.csv('Tables/TrendPlot_L2-3_Grouped_Genes_Inc_0.3_Dec_-0.3_Con_0.3_251211_.csv')
genes_l45 = read.csv('Tables/TrendPlot_L4-5_Grouped_Genes_Inc_0.3_Dec_-0.3_Con_0.3_251211_.csv')

genes_list = list()
genes_list[['L2-3_Increasing']] = genes_l23[genes_l23$group == "Increasing",]$gene
genes_list[['L2-3_Converging']] = genes_l23[genes_l23$group == "Converging",]$gene
genes_list[['L2-3_Decreasing']] = genes_l23[genes_l23$group == "Decreasing",]$gene

genes_list[['L4-5_Increasing']] = genes_l45[genes_l45$group == "Increasing",]$gene 
genes_list[['L4-5_Converging']] = genes_l45[genes_l45$group == "Converging",]$gene
genes_list[['L4-5_Decreasing']] = genes_l45[genes_l45$group == "Decreasing",]$gene


gostres = gost(
  genes_list[['L2-3_Increasing']],
  organism = "mmusculus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = TRUE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  sources = "GO",
  as_short_link = FALSE,
  highlight = FALSE
)

gostres_df = as.data.frame(gostres[["result"]])
table(gostres_df$significant)

## Run gost for all list elements

gost_results <- lapply(names(genes_list), function(name) {
  message("Running g:Profiler for ", name)
  
  res <- gost(
    genes_list[[name]],
    organism = "mmusculus",
    ordered_query = FALSE,
    multi_query = FALSE,
    significant = FALSE,
    exclude_iea = TRUE,
    measure_underrepresentation = FALSE,
    evcodes = TRUE,
    user_threshold = 0.05,
    correction_method = "g_SCS",
    sources = "GO",
    as_short_link = FALSE,
    highlight = FALSE
  )
  
  if (!is.null(res$result)) {
    df <- as.data.frame(res$result)
    df$gene_set <- name
    return(df)
  } else {
    return(NULL)
  }
})

# names(gost_results) = c("L2-3_Converging", "L2-3_Increasing", "L2-3_Decreasing", "L4-5_Increasing", "L4-5_Converging", "L4-5_Decreasing")

gost_results_df <- dplyr::bind_rows(gost_results)
gost_results_df_clean <- gost_results_df %>%
  mutate(across(
    where(is.list),
    ~ sapply(.x, function(y) paste(y, collapse = ";"))
  ))

write.csv(gost_results_df_clean, file = 'Tables/gprofileR_TrendPlot_Grouped_Genes_Inc_0.3_Dec_-0.3_Con_0.3_260112.csv', row.names = FALSE)


