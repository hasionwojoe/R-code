## Environment cleanup and load libraries
rm(list = ls())

library(tidyverse)
library(ggpubr)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(GSEABase)
library(enrichplot)
library(ggbreak)
library(fgsea)
library(ggforce)
library(KEGGREST)
library(reshape2)
library(readxl)
library(DESeq2)
library(sva)       # for ComBat_seq

## ==========================
## 1. Simplify group information
## ==========================
group1 <- group %>%
  rownames_to_column("sample") %>%
  filter(group == "Retina") %>%
  mutate(
    subtype = case_when(
      str_detect(subgroup, "normal") ~ "normal",
      str_detect(subgroup, "CNV") ~ "CNV",
      str_detect(subgroup, "AMD|GA") ~ "AMD",
      TRUE ~ "pre-AMD"
    ),
    subtype = factor(subtype, levels = c("normal", "pre-AMD", "AMD", "CNV"))
  )

## ==========================
## 2. Extract expression of genes of interest for plotting
## ==========================
gene_of_interest <- c('S1PR5', 'AIFM2', 'ALDH3B1')

filtered_exp_df <- as.data.frame(exp) %>%
  rownames_to_column("gene") %>%
  filter(gene %in% gene_of_interest) %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
  inner_join(group1, by = "sample")

## ==========================
## 3. Single-gene boxplot
## ==========================
target_gene <- 'ALDH3B1'  # alternatives: 'S1PR5', 'AIFM2'

p <- filtered_exp_df %>%
  filter(gene == target_gene) %>%
  ggplot(aes(x = subtype, y = expression, color = subtype)) +
  geom_boxplot(position = position_dodge(width = 0.2), width = 0.6, alpha = 0.7, outlier.size = 0) +
  geom_point(aes(fill = subtype), size = 1, alpha = 0.4,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6)) +
  scale_fill_manual(values = c("normal" = "#85A0CC", "pre-AMD" = "#AB92BF", "AMD" = "#92B65F", "CNV" = "#E48D83")) +
  scale_color_manual(values = c("normal" = "#85A0CC", "pre-AMD" = "#AB92BF", "AMD" = "#92B65F", "CNV" = "#E48D83")) +
  labs(
    title = paste0(target_gene, " in AMD"),
    y = "Expression Level",
    x = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.ticks = element_line(),
    axis.ticks.length = unit(0.1, "cm"),
    legend.position = "none"
  ) +
  stat_compare_means(
    aes(group = subtype),
    comparisons = list(c("normal", "pre-AMD"), c("normal", "AMD"), c("normal", "CNV")),
    label = "p.signif",
    method = "t.test"
  )

print(p)
ggsave(paste0(target_gene, "_expression_boxplot.pdf"), p, width = 4, height = 5)

## ==========================
## 4. Differential expression analysis (CNV vs normal)
## ==========================
exp_limma <- exp[, group1$sample]
design <- model.matrix(~ 0 + subtype, data = group1)
colnames(design) <- make.names(colnames(design))

contrast <- makeContrasts(CNV_vs_normal = subtypeCNV - subtypenormal, levels = design)

fit <- lmFit(exp_limma, design) %>%
  contrasts.fit(contrast) %>%
  eBayes()

results <- topTable(fit, coef = "CNV_vs_normal", adjust.method = "fdr", number = Inf)
significant_genes <- results %>%
  filter(P.Value < 0.05, abs(logFC) > 1)

## ==========================
## 5. KEGG enrichment analysis
## ==========================
diff_gene <- significant_genes %>%
  rownames_to_column('Symbol') %>%
  dplyr::select(Symbol, logFC)

s2e <- bitr(unique(diff_gene$Symbol), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
diff_gene <- inner_join(diff_gene, s2e, by = c("Symbol" = "SYMBOL"))

KEGG <- enrichKEGG(gene = diff_gene$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
kegg_result <- as.data.frame(KEGG@result) %>% filter(pvalue < 0.05)
write.csv(kegg_result, "kegg_result.csv", row.names = FALSE)

## ==========================
## 6. GSEA analysis (KEGG)
## ==========================
gene_df <- results %>%
  rownames_to_column("gene_id")

gene <- gene_df$gene_id %>%
  bitr(fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::rename(gene_id = SYMBOL)

gene_df <- inner_join(gene_df, gene, by = "gene_id")

geneList <- gene_df$logFC
names(geneList) <- gene_df$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

GSEA <- gseKEGG(
  geneList,
  organism = "hsa",
  pvalueCutoff = 1,
  pAdjustMethod = "BH"
)

gsea <- GSEA@result %>% filter(pvalue < 0.05)

## Collect KEGG pathway class information
pathway_ids <- gsea$ID
pathway_info <- lapply(pathway_ids, keggGet)
gsea$class_joined <- sapply(pathway_info, function(info) {
  class_info <- info[[1]]$CLASS
  if (!is.null(class_info)) paste(class_info, collapse = ", ") else NA
})

write.csv(gsea, "amd-gsea.csv", row.names = FALSE)

save(list = ls(all.names = TRUE), file = "kegg_and_gsea.RData")

## ==========================
## 7. Donut (ring) plot of pathway categories (example)
## ==========================
pathway_data <- data.frame(
  category = c("Amino acid metabolism", "Cardiovascular disease", "Nervous system",
               "Signal transduction", "Transport and catabolism", "others"),
  value = c(6, 5, 4, 4, 3, 13)
) %>%
  mutate(
    prop = value / sum(value),
    end_angle = cumsum(prop) * 2 * pi,
    start_angle = lag(end_angle, default = 0),
    mid_angle = (start_angle + end_angle) / 2
  )

color_palette <- c("#7AA6D9", "#EFA86E", "#88C4B5", "#C989B6", "#D9CA8F", "#A699A6")

ggplot(pathway_data) +
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.5, r = 1,
                   start = start_angle, end = end_angle,
                   fill = category),
               color = "white", size = 0.5) +
  geom_text(
    aes(x = 0.7 * sin(mid_angle), y = 0.7 * cos(mid_angle),
        label = sprintf("%.2f%%", prop * 100),
        angle = ifelse(mid_angle > pi, 180 * mid_angle/pi + 180, 180 * mid_angle/pi)),
    color = "white", size = 3.5, fontface = "bold"
  ) +
  annotate("text", x = 0, y = 0, label = "Functional\nCategories",
           size = 5, fontface = "bold", color = "#333333") +
  coord_fixed() +
  scale_fill_manual(values = color_palette) +
  theme_void() +
  labs(title = "Pathway Ratio in AMD")

ggsave("Pathway_ratio_in_AMD.pdf", width = 8, height = 8)

## ==========================
## Batch correction for two datasets and QC (DE/Enrichment same as above)
## ==========================
common_genes <- intersect(rownames(dr1_exp), rownames(dr2_exp))
expr1 <- dr1_exp[common_genes, ]
expr2 <- dr2_exp[common_genes, ]
expr_combined <- cbind(expr1, expr2)
batch <- c(rep("GSE102485", ncol(expr1)), rep("GSE179568", ncol(expr2)))

# Ensure integer counts
expr_combined <- round(expr_combined)

# Batch correction using ComBat_seq
expr_corrected <- ComBat_seq(counts = as.matrix(expr_combined), batch = as.factor(batch))

## ==========================
## PCA before and after batch correction
## ==========================
# Before correction
dds <- DESeqDataSetFromMatrix(countData = expr_combined, colData = data.frame(batch), design = ~1)
vsd <- vst(dds)
pcaData <- plotPCA(vsd, intgroup = "batch", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = batch)) +
  geom_point(size = 3) +
  labs(title = "PCA before batch correction",
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) +
  theme_bw()

# After correction
dds_corrected <- DESeqDataSetFromMatrix(countData = expr_corrected, colData = data.frame(batch), design = ~1)
vsd_corrected <- vst(dds_corrected)
pcaData_corrected <- plotPCA(vsd_corrected, intgroup = "batch", returnData = TRUE)
percentVar_corrected <- round(100 * attr(pcaData_corrected, "percentVar"))

ggplot(pcaData_corrected, aes(PC1, PC2, color = batch)) +
  geom_point(size = 3) +
  labs(title = "PCA after batch correction",
       x = paste0("PC1: ", percentVar_corrected[1], "%"),
       y = paste0("PC2: ", percentVar_corrected[2], "%")) +
  theme_bw()

## ==========================
## Boxplots before and after batch correction
## ==========================
# Before correction
log_expr_before <- log2(expr_combined + 1)
df_before <- log_expr_before %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
  inner_join(pcaData, by = c("Sample" = "name"))

ggplot(df_before, aes(x = Sample, y = Expression, fill = batch)) +
  geom_boxplot(outlier.size = 0.2) +
  theme_bw(base_size = 12) +
  labs(title = "Boxplot before batch correction", x = "Sample", y = "log2(Expression+1)") +
  theme(axis.text.x = element_blank())

# After correction
log_expr_after <- log2(expr_corrected + 1)
df_after <- melt(log_expr_after)
colnames(df_after) <- c("Gene", "Sample", "Expression")
df_after <- inner_join(df_after, pcaData, by = c("Sample" = "name"))

ggplot(df_after, aes(x = Sample, y = Expression, fill = batch)) +
  geom_boxplot(outlier.size = 0.2) +
  theme_bw(base_size = 12) +
  labs(title = "Boxplot after batch correction", x = "Sample", y = "log2(Expression+1)") +
  theme(axis.text.x = element_blank())

## ==========================
## Correlation calculations (e.g., vs S1PR5)
## ==========================
# Ensure gene S1PR5 exists
if (!"S1PR5" %in% rownames(expr_corrected)) {
  stop("S1PR5 is not present in the expression matrix. Please check gene naming consistency.")
}

# Transpose expression matrix to compute correlations per sample
exp_t <- t(expr_corrected)

# Extract S1PR5 expression vector
s1pr5_expr <- exp_t[, "S1PR5"]

# Compute Pearson correlation and p-value for each gene vs S1PR5
cor_results <- apply(exp_t, 2, function(gene_expr) {
  test <- cor.test(gene_expr, s1pr5_expr, use = "complete.obs")
  c(correlation = test$estimate, p_value = test$p.value)
})

# Organize results
cor_df <- as.data.frame(t(cor_results)) %>%
  rownames_to_column("gene") %>%
  filter(gene != "S1PR5") %>%
  arrange(desc(abs(correlation.cor)))
colnames(cor_df)[2] <- "correlation"

# Export correlations
write_csv(cor_df, "cor_df.csv")

## ==========================
## Ferroptosis-related gene correlation forest plot
## ==========================
ferroptosis_genes <- read_excel("ferroptosis geneset.xlsx")$symbol
ferro_cor <- cor_df %>% filter(gene %in% ferroptosis_genes)
write.csv(ferro_cor, "ferroptosis_correlations.csv", row.names = FALSE)

ggplot(ferro_cor, aes(x = fct_reorder(gene, correlation), y = correlation)) +
  geom_segment(aes(xend = gene, y = 0, yend = correlation), color = "grey") +
  geom_point() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Ferroptosis Gene Correlations with S1PR5", x = "Gene", y = "Correlation")

## ==========================
## S1PR5 vs AIFM2 scatter plot and correlation
## ==========================
exp_pair <- exp_long %>%
  filter(gene %in% c("S1PR5", "AIFM2")) %>%
  pivot_wider(names_from = gene, values_from = expression)

ggplot(exp_pair, aes(x = S1PR5, y = AIFM2)) +
  geom_point() +
  theme_minimal() +
  geom_smooth() +
  labs(title = "S1PR5 vs AIFM2 Expression", x = "S1PR5 Expression", y = "AIFM2 Expression")

print(cor(exp_pair$S1PR5, exp_pair$AIFM2, use = "complete.obs"))

## ==========================
## Volcano plot for ferroptosis-related correlations
## ==========================
volcano_data <- read_csv("ferroptosis_correlations.csv")

gene_of_interest <- c('MAP1LC3A','AKR1C3','FDFT1','ACSL4','NOX4','DPP4','GCH1','ACSL5')

volcano_data <- volcano_data %>%
  mutate(
    label = if_else(gene %in% gene_of_interest, gene, NA_character_),
    negLog10p = -log10(p_value + 1e-300),
    cor_group = case_when(
      correlation > 0.3 ~ "Pos cor (>0.3)",
      correlation < -0.3 ~ "Neg cor (<-0.3)",
      TRUE ~ "Non-sig"
    )
  )

y_max <- max(volcano_data$negLog10p, na.rm = TRUE) * 1.05

p <- ggplot(volcano_data, aes(x = correlation, y = negLog10p)) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = 'dashed', color = 'gray60', linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', color = 'gray40', alpha = 0.8) +
  geom_point(aes(fill = cor_group), shape = 21, color = "gray30", size = 2.5, stroke = 0.2) +
  scale_fill_manual(
    values = c("Neg cor (<-0.3)" = "#6F93CB",
               "Non-sig" = "gray80",
               "Pos cor (>0.3)" = "#F1756D"),
    name = "Correlation"
  ) +
  geom_label_repel(
    data = filter(volcano_data, !is.na(label)),
    aes(label = label, fill = cor_group),
    color = "white",
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.4,
    min.segment.length = 0.2,
    segment.color = "gray30",
    segment.size = 0.2,
    size = 4,
    show.legend = FALSE
  ) +
  coord_cartesian(xlim = c(-1, 1.5), ylim = c(0, y_max)) +
  labs(
    title = "Ferroptosis-Related Gene Correlations",
    subtitle = "Significant genes labeled (|r| > 0.3, p < 0.05)",
    x = "Correlation Coefficient",
    y = expression(-Log[10](italic(p)-value))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
    panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.4),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 10)
  )

print(p)
ggsave("S1PR5_Ferroptosis_Correlation_Volcano.pdf", p, width = 7, height = 7)

## ==========================
## Pathway correlation with genes (example based on GSEA results)
## ==========================
# Load correlation table
cor_gene <- read_csv("cor_df.csv")

# Create multiple threshold subsets for downstream analysis
thresholds <- c(0.3, 0.4, 0.5, 0.6, 0.7)
cor_gene_list <- map(thresholds, ~ {
  cor_gene %>% filter(abs(correlation) > .x)
}) %>% set_names(paste0("cor_gene", thresholds))

# Example: use threshold 0.6
cor_gene0.6 <- cor_gene_list$cor_gene0.6

# ID conversion and deduplication
GSEA_cor <- bitr(
  cor_gene0.6$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  rename(gene = SYMBOL)

cor_gene0.6 <- cor_gene0.6 %>%
  inner_join(GSEA_cor, by = "gene")

# Create named vector for GSEA (names are ENTREZ IDs)
geneList_cor <- cor_gene0.6$correlation
names(geneList_cor) <- cor_gene0.6$ENTREZID
geneList_cor <- sort(geneList_cor, decreasing = TRUE)

# Run GSEA (KEGG)
GSEA_0.6 <- gseKEGG(
  geneList_cor,
  organism = "hsa",
  pvalueCutoff = 1,
  pAdjustMethod = "BH"
)

gsea_0.6 <- GSEA_0.6@result
write_csv(gsea_0.6, "GSEA_0.6_KEGG_results.csv")

# Example GSEA plot for selected pathways (replace IDs with significant ones as needed)
paths <- c("hsa04660", "hsa04650", "hsa04512")
p <- gseaplot2(GSEA_0.6, geneSetID = paths, pvalue_table = TRUE)

print(p)
ggsave("GSEA_0.6_selected_paths.pdf", p, width = 10, height = 8)