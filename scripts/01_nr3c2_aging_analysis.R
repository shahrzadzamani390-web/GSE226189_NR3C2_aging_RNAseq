# NR3C2 aging analysis in GSE226189 primary skin fibroblast RNA-seq
# Data source: GEO GSE226189
# Counts: featureCounts gene-level counts, hg19 Ensembl v82

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(GEOquery)
library(DESeq2)
library(ggplot2)

dir.create("data_raw/gene_counts", recursive = TRUE, showWarnings = FALSE)
dir.create("data_processed", showWarnings = FALSE)
dir.create("metadata", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

nr3c2_ensembl <- "ENSG00000151623"

gse <- getGEO("GSE226189", GSEMatrix = TRUE)
pheno <- pData(gse[[1]])

sample_metadata <- data.frame(
  sample_id = pheno$geo_accession,
  title = pheno$title,
  age = as.numeric(pheno$`age (years):ch1`),
  sex = pheno$`Sex:ch1`,
  cell_type = pheno$`cell type:ch1`
)

sample_metadata$age_group <- ifelse(
  sample_metadata$age < 40, "Young",
  ifelse(sample_metadata$age > 60, "Old", "Middle")
)

gene_count_urls <- pheno$supplementary_file_2
gene_count_files <- file.path("data_raw/gene_counts", basename(gene_count_urls))

for (i in seq_along(gene_count_urls)) {
  if (!file.exists(gene_count_files[i])) {
    message("Downloading ", i, " of ", length(gene_count_urls), ": ", basename(gene_count_files[i]))
    download.file(gene_count_urls[i], destfile = gene_count_files[i], mode = "wb")
  }
}

read_count_file <- function(file) {
  x <- read.delim(gzfile(file), header = TRUE, check.names = FALSE)
  sample_col <- colnames(x)[2]
  colnames(x) <- c("Tracking_ID", sample_col)
  x
}

count_list <- lapply(gene_count_files, read_count_file)

same_genes <- all(sapply(count_list, function(x) identical(x$Tracking_ID, count_list[[1]]$Tracking_ID)))
stopifnot(same_genes)

sample_names_from_files <- sapply(count_list, function(x) colnames(x)[2])
sample_names_from_files <- sub("_COUNT$", "", sample_names_from_files)

count_matrix <- data.frame(
  Tracking_ID = count_list[[1]]$Tracking_ID,
  do.call(cbind, lapply(count_list, function(x) x[[2]])),
  check.names = FALSE
)

colnames(count_matrix)[-1] <- sample_names_from_files

stopifnot(all(sample_metadata$title == colnames(count_matrix)[-1]))

write.csv(sample_metadata, "metadata/sample_metadata.csv", row.names = FALSE)
write.csv(count_matrix, "data_processed/gene_count_matrix.csv", row.names = FALSE)

# Young vs Old analysis, adjusted for sex
analysis_metadata <- sample_metadata[sample_metadata$age_group %in% c("Young", "Old"), ]
analysis_metadata$age_group <- factor(analysis_metadata$age_group, levels = c("Young", "Old"))
analysis_metadata$sex <- factor(analysis_metadata$sex)

counts_for_deseq <- count_matrix[, analysis_metadata$title]
rownames(counts_for_deseq) <- count_matrix$Tracking_ID
counts_for_deseq <- as.matrix(counts_for_deseq)
storage.mode(counts_for_deseq) <- "integer"

stopifnot(all(colnames(counts_for_deseq) == analysis_metadata$title))

dds <- DESeqDataSetFromMatrix(
  countData = counts_for_deseq,
  colData = analysis_metadata,
  design = ~ sex + age_group
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res_old_vs_young <- results(dds, name = "age_group_Old_vs_Young")
nr3c2_result <- res_old_vs_young[nr3c2_ensembl, ]

nr3c2_normalized <- counts(dds, normalized = TRUE)[nr3c2_ensembl, ]

nr3c2_table <- data.frame(
  sample_id = analysis_metadata$sample_id,
  title = analysis_metadata$title,
  age = analysis_metadata$age,
  sex = analysis_metadata$sex,
  age_group = analysis_metadata$age_group,
  nr3c2_normalized_count = as.numeric(nr3c2_normalized)
)

write.csv(as.data.frame(nr3c2_result), "results/NR3C2_DESeq2_old_vs_young_adjusted_for_sex.csv")
write.csv(nr3c2_table, "results/NR3C2_normalized_counts_young_old.csv", row.names = FALSE)

p1 <- ggplot(nr3c2_table, aes(x = age_group, y = nr3c2_normalized_count, color = sex)) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.85) +
  labs(
    title = "NR3C2 expression in young vs old fibroblast samples",
    subtitle = "DESeq2 normalized counts; model adjusted for sex",
    x = "Age group",
    y = "NR3C2 normalized count",
    color = "Sex"
  ) +
  theme_classic(base_size = 13)

ggsave("figures/NR3C2_young_vs_old_boxplot.png", p1, width = 6, height = 4, dpi = 300)

# Continuous age analysis, adjusted for sex
all_metadata <- sample_metadata
all_metadata$sex <- factor(all_metadata$sex)
all_metadata$age_centered <- all_metadata$age - mean(all_metadata$age)

all_counts_for_deseq <- count_matrix[, all_metadata$title]
rownames(all_counts_for_deseq) <- count_matrix$Tracking_ID
all_counts_for_deseq <- as.matrix(all_counts_for_deseq)
storage.mode(all_counts_for_deseq) <- "integer"

stopifnot(all(colnames(all_counts_for_deseq) == all_metadata$title))

dds_age_centered <- DESeqDataSetFromMatrix(
  countData = all_counts_for_deseq,
  colData = all_metadata,
  design = ~ sex + age_centered
)

dds_age_centered <- dds_age_centered[rowSums(counts(dds_age_centered)) >= 10, ]
dds_age_centered <- DESeq(dds_age_centered)

res_age_centered <- results(dds_age_centered, name = "age_centered")
nr3c2_age_centered_result <- res_age_centered[nr3c2_ensembl, ]

nr3c2_all_normalized <- counts(dds_age_centered, normalized = TRUE)[nr3c2_ensembl, ]

nr3c2_all_table <- data.frame(
  sample_id = all_metadata$sample_id,
  title = all_metadata$title,
  age = all_metadata$age,
  age_centered = all_metadata$age_centered,
  sex = all_metadata$sex,
  age_group = all_metadata$age_group,
  nr3c2_normalized_count = as.numeric(nr3c2_all_normalized)
)

write.csv(as.data.frame(nr3c2_age_centered_result), "results/NR3C2_DESeq2_continuous_age_adjusted_for_sex.csv")
write.csv(nr3c2_all_table, "results/NR3C2_normalized_counts_all_samples.csv", row.names = FALSE)

p2 <- ggplot(nr3c2_all_table, aes(x = age, y = nr3c2_normalized_count, color = sex)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  labs(
    title = "NR3C2 expression increases with age",
    subtitle = "DESeq2 normalized counts; model adjusted for sex",
    x = "Age (years)",
    y = "NR3C2 normalized count",
    color = "Sex"
  ) +
  theme_classic(base_size = 13)

ggsave("figures/NR3C2_expression_vs_age.png", p2, width = 6, height = 4, dpi = 300)

print(nr3c2_result)
print(nr3c2_age_centered_result)
# Quality control: PCA, library size, and sample distances
vsd <- vst(dds_age_centered, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = c("age_group", "sex"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

write.csv(pca_data, "results/PCA_sample_coordinates.csv", row.names = FALSE)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = age_group, shape = sex)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of GSE226189 fibroblast RNA-seq samples",
    subtitle = "Variance-stabilized gene counts",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Age group",
    shape = "Sex"
  ) +
  theme_classic(base_size = 13)

ggsave("figures/PCA_age_group_sex.png", p_pca, width = 6, height = 4.5, dpi = 300)

library_size <- colSums(count_matrix[, sample_metadata$title])

library_qc <- data.frame(
  sample_id = sample_metadata$sample_id,
  title = sample_metadata$title,
  age = sample_metadata$age,
  sex = sample_metadata$sex,
  age_group = sample_metadata$age_group,
  library_size = as.numeric(library_size)
)

write.csv(library_qc, "results/library_size_qc.csv", row.names = FALSE)

p_library <- ggplot(library_qc, aes(x = age_group, y = library_size / 1e6, color = sex)) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.85) +
  labs(
    title = "RNA-seq library size by age group",
    subtitle = "Total gene counts per sample",
    x = "Age group",
    y = "Library size (million counts)",
    color = "Sex"
  ) +
  theme_classic(base_size = 13)

ggsave("figures/library_size_qc.png", p_library, width = 6, height = 4, dpi = 300)

if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(pheatmap)

sample_dist <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dist)

rownames(sample_dist_matrix) <- colnames(vsd)
colnames(sample_dist_matrix) <- colnames(vsd)

annotation_col <- data.frame(
  age = sample_metadata$age,
  sex = sample_metadata$sex,
  age_group = sample_metadata$age_group
)

rownames(annotation_col) <- sample_metadata$title

pheatmap(
  sample_dist_matrix,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 8,
  filename = "figures/sample_distance_heatmap.png",
  width = 7,
  height = 6
)

# Genome-wide differential expression: Old vs Young adjusted for sex
res_all_old_vs_young <- results(dds, name = "age_group_Old_vs_Young")

res_all_table <- as.data.frame(res_all_old_vs_young)
res_all_table$ensembl_id <- rownames(res_all_table)

res_all_table <- res_all_table[, c(
  "ensembl_id",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

res_all_table <- res_all_table[order(res_all_table$pvalue), ]

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(AnnotationDbi)

res_all_table$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = res_all_table$ensembl_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res_all_table$gene_name <- mapIds(
  org.Hs.eg.db,
  keys = res_all_table$ensembl_id,
  column = "GENENAME",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res_all_table <- res_all_table[, c(
  "ensembl_id",
  "gene_symbol",
  "gene_name",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

write.csv(res_all_table, "results/DESeq2_all_genes_old_vs_young_adjusted_for_sex_annotated.csv", row.names = FALSE)

significant_genes <- res_all_table[
  !is.na(res_all_table$padj) & res_all_table$padj < 0.05,
]

significant_genes <- significant_genes[order(significant_genes$padj), ]

write.csv(significant_genes, "results/DESeq2_significant_genes_FDR05_old_vs_young_adjusted_for_sex.csv", row.names = FALSE)

volcano_table <- res_all_table
volcano_table$significance <- "Not significant"

volcano_table$significance[
  !is.na(volcano_table$padj) & volcano_table$padj < 0.05 & volcano_table$log2FoldChange > 0
] <- "Higher in Old"

volcano_table$significance[
  !is.na(volcano_table$padj) & volcano_table$padj < 0.05 & volcano_table$log2FoldChange < 0
] <- "Lower in Old"

p_volcano <- ggplot(volcano_table, aes(x = log2FoldChange, y = -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.6, size = 1.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(
    values = c(
      "Higher in Old" = "#D95F02",
      "Lower in Old" = "#1B9E77",
      "Not significant" = "gray70"
    )
  ) +
  labs(
    title = "Differential expression: Old vs Young fibroblasts",
    subtitle = "DESeq2 model adjusted for sex",
    x = "log2 fold-change (Old vs Young)",
    y = "-log10 nominal p-value",
    color = "FDR < 0.05"
  ) +
  theme_classic(base_size = 13)

ggsave("figures/volcano_old_vs_young_adjusted_for_sex.png", p_volcano, width = 6.5, height = 5, dpi = 300)

# Targeted MR-related gene panel analysis
mr_gene_panel <- c(
  "NR3C2",
  "HSD11B1",
  "HSD11B2",
  "FKBP4",
  "FKBP5",
  "HSP90AA1",
  "HSP90AB1",
  "NCOA1",
  "NCOA2",
  "NCOR1",
  "NCOR2",
  "NR3C1",
  "SGK1",
  "SCNN1A",
  "SCNN1B",
  "SCNN1G"
)

mr_panel_results <- res_all_table[
  !is.na(res_all_table$gene_symbol) & res_all_table$gene_symbol %in% mr_gene_panel,
]

mr_panel_results <- mr_panel_results[match(mr_gene_panel, mr_panel_results$gene_symbol), ]
mr_panel_results <- mr_panel_results[!is.na(mr_panel_results$gene_symbol), ]
mr_panel_results$panel_padj <- p.adjust(mr_panel_results$pvalue, method = "BH")

mr_panel_results <- mr_panel_results[, c(
  "ensembl_id",
  "gene_symbol",
  "gene_name",
  "baseMean",
  "log2FoldChange",
  "pvalue",
  "padj",
  "panel_padj"
)]

write.csv(mr_panel_results, "results/MR_gene_panel_old_vs_young_adjusted_for_sex.csv", row.names = FALSE)

mr_panel_plot_table <- mr_panel_results
mr_panel_plot_table$direction <- ifelse(
  mr_panel_plot_table$log2FoldChange > 0,
  "Higher in Old",
  "Lower in Old"
)

mr_panel_plot_table$gene_symbol <- factor(
  mr_panel_plot_table$gene_symbol,
  levels = mr_panel_plot_table$gene_symbol[order(mr_panel_plot_table$log2FoldChange)]
)

p_mr_panel <- ggplot(
  mr_panel_plot_table,
  aes(x = log2FoldChange, y = gene_symbol, color = direction, size = -log10(pvalue))
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(alpha = 0.9) +
  scale_color_manual(
    values = c(
      "Higher in Old" = "#D95F02",
      "Lower in Old" = "#1B9E77"
    )
  ) +
  labs(
    title = "MR-related gene panel: Old vs Young fibroblasts",
    subtitle = "DESeq2 model adjusted for sex",
    x = "log2 fold-change (Old vs Young)",
    y = "Gene",
    color = "Direction",
    size = "-log10 p-value"
  ) +
  theme_classic(base_size = 13)

ggsave("figures/MR_gene_panel_old_vs_young_dotplot.png", p_mr_panel, width = 7, height = 5, dpi = 300)
