# =============================================================================
# CBIO310 Multi-Omics Meta-Analysis Pipeline
# SARS-CoV-2 and Alzheimer's Disease
# =============================================================================
#
# Authors:
#   Tamim Alaa          -- 211002002
#   Hana Mohamed Hikal  -- 231000300
#   Asmaa Ragab         -- 221001773
#   Habiba Farag        -- 221001888
#
# =============================================================================

# Set working directory to the location of this script
# setwd("path/to/your/project")  # Uncomment and modify if needed

# -----------------------------------------------------------------------------
# PACKAGE INSTALLATION
# -----------------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_packages <- c("GEOquery", "DESeq2", "limma", "edgeR", "sva", "WGCNA",
                   "clusterProfiler", "enrichplot", "org.Hs.eg.db", "pathview",
                   "biomaRt", "ComplexHeatmap", "hgu133plus2.db")

cran_packages <- c("tidyverse", "ggplot2", "pheatmap", "VennDiagram", 
                   "ggVennDiagram", "RColorBrewer", "scales", "openxlsx",
                   "data.table", "ggrepel", "cowplot", "viridis")

BiocManager::install(bioc_packages, ask = FALSE, update = FALSE)
install.packages(cran_packages, dependencies = TRUE)

# -----------------------------------------------------------------------------
# LOAD LIBRARIES
# -----------------------------------------------------------------------------

library(GEOquery)
library(tidyverse)
library(DESeq2)
library(limma)
library(sva)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(edgeR)
library(hgu133plus2.db)
library(VennDiagram)
library(ggVennDiagram)
library(openxlsx)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(DOSE)
library(WGCNA)

enableWGCNAThreads()
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Create directories
dir.create("data/raw/sars-cov2", recursive = TRUE, showWarnings = FALSE)
dir.create("data/raw/alzheimers", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/degs", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("results/wgcna", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# DOWNLOAD GEO DATASETS
# -----------------------------------------------------------------------------

download_results <- list()

# Download GSE157103 (COVID-19 Blood)
gse_id <- "GSE157103"
output_path <- "data/raw/sars-cov2"
cat("Downloading:", gse_id, "\n")

tryCatch({
    gse <- getGEO(gse_id, GSEMatrix = TRUE, destdir = output_path)
    if (length(gse) > 0) {
        eset <- gse[[1]]
        write.csv(exprs(eset), file.path(output_path, paste0(gse_id, "_expression.csv")))
        write.csv(pData(eset), file.path(output_path, paste0(gse_id, "_phenotype.csv")))
        write.csv(fData(eset), file.path(output_path, paste0(gse_id, "_features.csv")))
        download_results[[gse_id]] <- list(status = "success", n_samples = ncol(exprs(eset)))
        cat("  Success:", ncol(exprs(eset)), "samples\n")
    }
}, error = function(e) cat("  Failed:", e$message, "\n"))

tryCatch({ getGEOSuppFiles(gse_id, baseDir = output_path) }, error = function(e) {})

# Download GSE164805 (COVID-19 Lung)
gse_id <- "GSE164805"
cat("Downloading:", gse_id, "\n")

tryCatch({
    gse <- getGEO(gse_id, GSEMatrix = TRUE, destdir = output_path)
    if (length(gse) > 0) {
        eset <- gse[[1]]
        write.csv(exprs(eset), file.path(output_path, paste0(gse_id, "_expression.csv")))
        write.csv(pData(eset), file.path(output_path, paste0(gse_id, "_phenotype.csv")))
        write.csv(fData(eset), file.path(output_path, paste0(gse_id, "_features.csv")))
        download_results[[gse_id]] <- list(status = "success", n_samples = ncol(exprs(eset)))
        cat("  Success:", ncol(exprs(eset)), "samples\n")
    }
}, error = function(e) cat("  Failed:", e$message, "\n"))

tryCatch({ getGEOSuppFiles(gse_id, baseDir = output_path) }, error = function(e) {})

# Download GSE5281 (Alzheimer's Brain)
gse_id <- "GSE5281"
output_path <- "data/raw/alzheimers"
cat("Downloading:", gse_id, "\n")

tryCatch({
    gse <- getGEO(gse_id, GSEMatrix = TRUE, destdir = output_path)
    if (length(gse) > 0) {
        eset <- gse[[1]]
        write.csv(exprs(eset), file.path(output_path, paste0(gse_id, "_expression.csv")))
        write.csv(pData(eset), file.path(output_path, paste0(gse_id, "_phenotype.csv")))
        write.csv(fData(eset), file.path(output_path, paste0(gse_id, "_features.csv")))
        download_results[[gse_id]] <- list(status = "success", n_samples = ncol(exprs(eset)))
        cat("  Success:", ncol(exprs(eset)), "samples\n")
    }
}, error = function(e) cat("  Failed:", e$message, "\n"))

tryCatch({ getGEOSuppFiles(gse_id, baseDir = output_path) }, error = function(e) {})

saveRDS(download_results, "data/download_log.rds")

# -----------------------------------------------------------------------------
# DATA PREPROCESSING
# -----------------------------------------------------------------------------

sars_data <- list()
ad_data <- list()

# Process SARS-CoV-2 datasets
cat("\nProcessing SARS-CoV-2 datasets...\n")
sars_dir <- "data/raw/sars-cov2"
expr_files <- list.files(sars_dir, pattern = "_expression.csv", full.names = TRUE)

for (expr_file in expr_files) {
    gse_id <- gsub("_expression.csv", "", basename(expr_file))
    cat("  Processing:", gse_id, "\n")
    
    expr_matrix <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))
    pheno_file <- gsub("expression", "phenotype", expr_file)
    pheno_data <- read.csv(pheno_file, row.names = 1, stringsAsFactors = FALSE)
    
    # Log transform if needed
    if (max(expr_matrix, na.rm = TRUE) > 100) {
        expr_matrix <- log2(expr_matrix + 1)
    }
    
    # Quantile normalization
    expr_normalized <- normalizeBetweenArrays(expr_matrix, method = "quantile")
    
    # QC plots
    cor_matrix <- cor(expr_normalized, use = "pairwise.complete.obs")
    pdf(file.path("results/figures", paste0(gse_id, "_correlation.pdf")), width = 10, height = 10)
    pheatmap(cor_matrix, main = paste(gse_id, "- Sample Correlation"),
             color = colorRampPalette(brewer.pal(9, "Blues"))(100),
             show_rownames = FALSE, show_colnames = FALSE)
    dev.off()
    
    expr_long <- as.data.frame(expr_normalized) %>%
        pivot_longer(everything(), names_to = "Sample", values_to = "Expression")
    p <- ggplot(expr_long, aes(x = Expression, color = Sample)) +
        geom_density(show.legend = FALSE) + theme_bw() +
        labs(title = paste(gse_id, "- Expression Density"))
    ggsave(file.path("results/figures", paste0(gse_id, "_density.pdf")), p, width = 10, height = 6)
    
    sars_data[[gse_id]] <- list(expression = expr_normalized, phenotype = pheno_data)
    saveRDS(sars_data[[gse_id]], file.path("data/processed", paste0(gse_id, "_processed.rds")))
}

# Process Alzheimer's datasets
cat("Processing Alzheimer's datasets...\n")
ad_dir <- "data/raw/alzheimers"
expr_files <- list.files(ad_dir, pattern = "_expression.csv", full.names = TRUE)

for (expr_file in expr_files) {
    gse_id <- gsub("_expression.csv", "", basename(expr_file))
    cat("  Processing:", gse_id, "\n")
    
    expr_matrix <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))
    pheno_file <- gsub("expression", "phenotype", expr_file)
    pheno_data <- read.csv(pheno_file, row.names = 1, stringsAsFactors = FALSE)
    
    if (max(expr_matrix, na.rm = TRUE) > 100) {
        expr_matrix <- log2(expr_matrix + 1)
    }
    
    expr_normalized <- normalizeBetweenArrays(expr_matrix, method = "quantile")
    
    cor_matrix <- cor(expr_normalized, use = "pairwise.complete.obs")
    pdf(file.path("results/figures", paste0(gse_id, "_correlation.pdf")), width = 10, height = 10)
    pheatmap(cor_matrix, main = paste(gse_id, "- Sample Correlation"),
             color = colorRampPalette(brewer.pal(9, "Blues"))(100),
             show_rownames = FALSE, show_colnames = FALSE)
    dev.off()
    
    expr_long <- as.data.frame(expr_normalized) %>%
        pivot_longer(everything(), names_to = "Sample", values_to = "Expression")
    p <- ggplot(expr_long, aes(x = Expression, color = Sample)) +
        geom_density(show.legend = FALSE) + theme_bw() +
        labs(title = paste(gse_id, "- Expression Density"))
    ggsave(file.path("results/figures", paste0(gse_id, "_density.pdf")), p, width = 10, height = 6)
    
    ad_data[[gse_id]] <- list(expression = expr_normalized, phenotype = pheno_data)
    saveRDS(ad_data[[gse_id]], file.path("data/processed", paste0(gse_id, "_processed.rds")))
}

# -----------------------------------------------------------------------------
# DEG ANALYSIS
# -----------------------------------------------------------------------------

# DEG parameters
DEG_log2fc_threshold <- 1
DEG_pvalue_threshold <- 0.05
all_deg_results <- list()

# Process SARS-CoV-2 DEGs
cat("\nDEG Analysis: SARS-CoV-2...\n")
processed_files <- list.files("data/processed", pattern = "GSE.*_processed.rds", full.names = TRUE)
sars_files <- processed_files[grep("157103|164805|188847", processed_files)]

for (rds_file in sars_files) {
    gse_id <- gsub("_processed.rds", "", basename(rds_file))
    cat("  Analyzing:", gse_id, "\n")
    
    data <- readRDS(rds_file)
    expr_matrix <- data$expression
    pheno_data <- data$phenotype
    
    # Find group column
    group_cols <- c("disease state:ch1", "condition:ch1", "group:ch1", 
                    "disease_state", "condition", "characteristics_ch1")
    group_col <- NA
    for (col in group_cols) {
        if (col %in% colnames(pheno_data) && length(unique(pheno_data[[col]])) >= 2) {
            group_col <- col
            break
        }
    }
    
    if (is.na(group_col)) { cat("    No valid group column found\n"); next }
    
    group_vector <- pheno_data[[group_col]]
    group_clean <- make.names(group_vector)
    
    # Limma analysis
    design <- model.matrix(~ 0 + factor(group_clean))
    colnames(design) <- levels(factor(group_clean))
    fit <- lmFit(expr_matrix, design)
    
    groups <- levels(factor(group_clean))
    contrast_formula <- paste(groups[1], "-", groups[2])
    contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
    
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    deg_results <- topTable(fit2, adjust.method = "BH", number = Inf, sort.by = "P")
    deg_results$Gene <- rownames(deg_results)
    deg_results <- deg_results %>%
        mutate(DEG_status = case_when(
            logFC >= DEG_log2fc_threshold & adj.P.Val < DEG_pvalue_threshold ~ "Upregulated",
            logFC <= -DEG_log2fc_threshold & adj.P.Val < DEG_pvalue_threshold ~ "Downregulated",
            TRUE ~ "Not significant"
        ))
    
    n_up <- sum(deg_results$DEG_status == "Upregulated")
    n_down <- sum(deg_results$DEG_status == "Downregulated")
    cat("    DEGs: Up=", n_up, ", Down=", n_down, "\n")
    
    write.csv(deg_results, file.path("results/degs", paste0(gse_id, "_DEG_results.csv")), row.names = FALSE)
    
    # Volcano plot
    plot_data <- deg_results %>% mutate(neg_log10_pval = -log10(adj.P.Val + 1e-300))
    p <- ggplot(plot_data, aes(x = logFC, y = neg_log10_pval, color = DEG_status)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_hline(yintercept = -log10(DEG_pvalue_threshold), linetype = "dashed") +
        geom_vline(xintercept = c(-DEG_log2fc_threshold, DEG_log2fc_threshold), linetype = "dashed") +
        scale_color_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "Not significant" = "gray70")) +
        labs(title = paste(gse_id, "- Volcano Plot"), x = "log2 Fold Change", y = "-log10(adj. p-value)") +
        theme_bw()
    ggsave(file.path("results/figures", paste0(gse_id, "_volcano.pdf")), p, width = 10, height = 8)
    
    # Heatmap
    top_degs <- deg_results %>% filter(DEG_status != "Not significant") %>% arrange(adj.P.Val) %>% head(50) %>% pull(Gene)
    if (length(top_degs) >= 2) {
        expr_subset <- expr_matrix[rownames(expr_matrix) %in% top_degs, ]
        if (nrow(expr_subset) >= 2) {
            tryCatch({
                pdf(file.path("results/figures", paste0(gse_id, "_heatmap.pdf")), width = 12, height = 10)
                pheatmap(t(scale(t(expr_subset))), main = paste(gse_id, "- Top DEGs"),
                         annotation_col = data.frame(Group = group_vector, row.names = colnames(expr_subset)),
                         show_colnames = FALSE, color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
                dev.off()
            }, error = function(e) {})
        }
    }
    
    all_deg_results[[gse_id]] <- deg_results
}

# Process Alzheimer's DEGs
cat("DEG Analysis: Alzheimer's...\n")
ad_files <- processed_files[grep("5281", processed_files)]

for (rds_file in ad_files) {
    gse_id <- gsub("_processed.rds", "", basename(rds_file))
    cat("  Analyzing:", gse_id, "\n")
    
    data <- readRDS(rds_file)
    expr_matrix <- data$expression
    pheno_data <- data$phenotype
    
    group_cols <- c("Disease State:ch1", "disease.state.ch1", "disease state:ch1", 
                    "diagnosis:ch1", "condition:ch1", "disease_state")
    group_col <- NA
    for (col in group_cols) {
        if (col %in% colnames(pheno_data) && length(unique(pheno_data[[col]])) >= 2) {
            group_col <- col
            break
        }
    }
    
    if (is.na(group_col)) { cat("    No valid group column found\n"); next }
    
    group_vector <- pheno_data[[group_col]]
    
    # Filter NA values
    valid_idx <- !is.na(group_vector) & group_vector != "" & !grepl("^\\s*$", group_vector)
    if (sum(!valid_idx) > 0) {
        expr_matrix <- expr_matrix[, valid_idx]
        group_vector <- group_vector[valid_idx]
        pheno_data <- pheno_data[valid_idx, ]
    }
    
    if (length(unique(group_vector)) != 2) { cat("    Need exactly 2 groups\n"); next }
    
    group_clean <- make.names(group_vector)
    design <- model.matrix(~ 0 + factor(group_clean))
    colnames(design) <- levels(factor(group_clean))
    fit <- lmFit(expr_matrix, design)
    
    groups <- levels(factor(group_clean))
    contrast_formula <- paste(groups[1], "-", groups[2])
    contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
    
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    deg_results <- topTable(fit2, adjust.method = "BH", number = Inf, sort.by = "P")
    deg_results$Gene <- rownames(deg_results)
    deg_results <- deg_results %>%
        mutate(DEG_status = case_when(
            logFC >= DEG_log2fc_threshold & adj.P.Val < DEG_pvalue_threshold ~ "Upregulated",
            logFC <= -DEG_log2fc_threshold & adj.P.Val < DEG_pvalue_threshold ~ "Downregulated",
            TRUE ~ "Not significant"
        ))
    
    # Map probes to gene symbols
    if (grepl("5281", gse_id)) {
        gene_symbols <- mapIds(hgu133plus2.db, keys = deg_results$Gene, column = "SYMBOL",
                               keytype = "PROBEID", multiVals = "first")
        deg_results$GeneSymbol <- gene_symbols[deg_results$Gene]
        deg_results <- deg_results[!is.na(deg_results$GeneSymbol), ]
        deg_results$ProbeID <- deg_results$Gene
        deg_results$Gene <- deg_results$GeneSymbol
    }
    
    n_up <- sum(deg_results$DEG_status == "Upregulated")
    n_down <- sum(deg_results$DEG_status == "Downregulated")
    cat("    DEGs: Up=", n_up, ", Down=", n_down, "\n")
    
    write.csv(deg_results, file.path("results/degs", paste0(gse_id, "_DEG_results.csv")), row.names = FALSE)
    
    # Volcano plot
    plot_data <- deg_results %>% mutate(neg_log10_pval = -log10(adj.P.Val + 1e-300))
    p <- ggplot(plot_data, aes(x = logFC, y = neg_log10_pval, color = DEG_status)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_hline(yintercept = -log10(DEG_pvalue_threshold), linetype = "dashed") +
        geom_vline(xintercept = c(-DEG_log2fc_threshold, DEG_log2fc_threshold), linetype = "dashed") +
        scale_color_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "Not significant" = "gray70")) +
        labs(title = paste(gse_id, "- Volcano Plot"), x = "log2 Fold Change", y = "-log10(adj. p-value)") +
        theme_bw()
    ggsave(file.path("results/figures", paste0(gse_id, "_volcano.pdf")), p, width = 10, height = 8)
    
    # Heatmap
    top_degs <- deg_results %>% filter(DEG_status != "Not significant") %>% arrange(adj.P.Val) %>% head(50)
    if ("ProbeID" %in% colnames(top_degs)) top_degs <- top_degs$ProbeID else top_degs <- top_degs$Gene
    if (length(top_degs) >= 2) {
        expr_subset <- expr_matrix[rownames(expr_matrix) %in% top_degs, ]
        if (nrow(expr_subset) >= 2) {
            tryCatch({
                pdf(file.path("results/figures", paste0(gse_id, "_heatmap.pdf")), width = 12, height = 10)
                pheatmap(t(scale(t(expr_subset))), main = paste(gse_id, "- Top DEGs"),
                         annotation_col = data.frame(Group = group_vector, row.names = colnames(expr_subset)),
                         show_colnames = FALSE, color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
                dev.off()
            }, error = function(e) {})
        }
    }
    
    all_deg_results[[gse_id]] <- deg_results
}

saveRDS(all_deg_results, "results/degs/all_DEG_results.rds")

# -----------------------------------------------------------------------------
# INTERSECTION ANALYSIS
# -----------------------------------------------------------------------------

cat("\nIntersection Analysis...\n")

# Load DEG results
deg_files <- list.files("results/degs", pattern = "_DEG_results.csv", full.names = TRUE)
deg_list <- list()
for (f in deg_files) {
    gse_id <- gsub("_DEG_results.csv", "", basename(f))
    deg_list[[gse_id]] <- read.csv(f)
}

# Separate by disease
sars_ids <- names(deg_list)[grep("157103|164805|188847", names(deg_list))]
ad_ids <- names(deg_list)[grep("5281", names(deg_list))]

# Get DEGs
sars_degs_all <- unique(unlist(lapply(sars_ids, function(id) 
    deg_list[[id]] %>% filter(DEG_status != "Not significant") %>% pull(Gene))))
sars_degs_up <- unique(unlist(lapply(sars_ids, function(id) 
    deg_list[[id]] %>% filter(DEG_status == "Upregulated") %>% pull(Gene))))
sars_degs_down <- unique(unlist(lapply(sars_ids, function(id) 
    deg_list[[id]] %>% filter(DEG_status == "Downregulated") %>% pull(Gene))))

ad_degs_all <- unique(unlist(lapply(ad_ids, function(id) 
    deg_list[[id]] %>% filter(DEG_status != "Not significant") %>% pull(Gene))))
ad_degs_up <- unique(unlist(lapply(ad_ids, function(id) 
    deg_list[[id]] %>% filter(DEG_status == "Upregulated") %>% pull(Gene))))
ad_degs_down <- unique(unlist(lapply(ad_ids, function(id) 
    deg_list[[id]] %>% filter(DEG_status == "Downregulated") %>% pull(Gene))))

cat("  SARS-CoV-2 DEGs:", length(sars_degs_all), "\n")
cat("  Alzheimer's DEGs:", length(ad_degs_all), "\n")

# Shared genes
shared_all <- intersect(sars_degs_all, ad_degs_all)
shared_up <- intersect(sars_degs_up, ad_degs_up)
shared_down <- intersect(sars_degs_down, ad_degs_down)
concordant_up <- intersect(sars_degs_up, ad_degs_up)
concordant_down <- intersect(sars_degs_down, ad_degs_down)
discordant_1 <- intersect(sars_degs_up, ad_degs_down)
discordant_2 <- intersect(sars_degs_down, ad_degs_up)

cat("  Shared DEGs:", length(shared_all), "\n")
cat("  Concordant:", length(concordant_up) + length(concordant_down), "\n")
cat("  Discordant:", length(discordant_1) + length(discordant_2), "\n")

# Venn diagrams
venn.plot <- venn.diagram(x = list(sars_degs_all, ad_degs_all),
    category.names = c("SARS-CoV-2", "Alzheimer's"), filename = NULL,
    fill = c("#E41A1C", "#377EB8"), alpha = 0.5, cex = 1.5, cat.cex = 1.2)
pdf("results/figures/venn_all_degs.pdf", width = 8, height = 8)
grid.draw(venn.plot)
dev.off()

venn.plot <- venn.diagram(x = list(sars_degs_up, ad_degs_up),
    category.names = c("COVID Up", "AD Up"), filename = NULL,
    fill = c("#E41A1C", "#377EB8"), alpha = 0.5)
pdf("results/figures/venn_upregulated.pdf", width = 8, height = 8)
grid.draw(venn.plot)
dev.off()

venn.plot <- venn.diagram(x = list(sars_degs_down, ad_degs_down),
    category.names = c("COVID Down", "AD Down"), filename = NULL,
    fill = c("#E41A1C", "#377EB8"), alpha = 0.5)
pdf("results/figures/venn_downregulated.pdf", width = 8, height = 8)
grid.draw(venn.plot)
dev.off()

p <- ggVennDiagram(list("SARS-CoV-2" = sars_degs_all, "Alzheimer's" = ad_degs_all), label_alpha = 0) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    theme(legend.position = "none")
ggsave("results/figures/venn_publication.pdf", p, width = 8, height = 8)

# Save results
shared_results <- list(shared_all = shared_all, shared_up = shared_up, shared_down = shared_down,
    concordant_up = concordant_up, concordant_down = concordant_down,
    discordant_covid_up_ad_down = discordant_1, discordant_covid_down_ad_up = discordant_2,
    sars_degs = sars_degs_all, ad_degs = ad_degs_all)
saveRDS(shared_results, "results/degs/shared_degs_results.rds")

# Excel export
wb <- createWorkbook()
addWorksheet(wb, "Shared_All"); writeData(wb, "Shared_All", data.frame(Gene = shared_all))
addWorksheet(wb, "Concordant_Up"); writeData(wb, "Concordant_Up", data.frame(Gene = concordant_up))
addWorksheet(wb, "Concordant_Down"); writeData(wb, "Concordant_Down", data.frame(Gene = concordant_down))
addWorksheet(wb, "Summary"); writeData(wb, "Summary", data.frame(
    Category = c("SARS-CoV-2 DEGs", "AD DEGs", "Shared", "Concordant Up", "Concordant Down", "Discordant"),
    Count = c(length(sars_degs_all), length(ad_degs_all), length(shared_all),
              length(concordant_up), length(concordant_down), length(discordant_1) + length(discordant_2))))
saveWorkbook(wb, "results/degs/shared_degs_analysis.xlsx", overwrite = TRUE)

# -----------------------------------------------------------------------------
# ENRICHMENT ANALYSIS
# -----------------------------------------------------------------------------

cat("\nEnrichment Analysis...\n")

all_deg_results <- readRDS("results/degs/all_DEG_results.rds")
if (!"GSE5281" %in% names(all_deg_results) && file.exists("results/degs/GSE5281_DEG_results.csv")) {
    all_deg_results$GSE5281 <- read.csv("results/degs/GSE5281_DEG_results.csv")
    saveRDS(all_deg_results, "results/degs/all_DEG_results.rds")
}

go_covid <- NULL; kegg_covid <- NULL
go_ad <- NULL; kegg_ad <- NULL

# COVID-19 enrichment
if ("GSE157103" %in% names(all_deg_results)) {
    covid_degs <- all_deg_results$GSE157103 %>% filter(DEG_status != "Not significant") %>% pull(Gene)
    cat("  COVID-19 DEGs:", length(covid_degs), "\n")
    
    if (length(covid_degs) >= 20) {
        gene_ids <- tryCatch({ bitr(covid_degs, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) },
                             error = function(e) data.frame())
        
        if (nrow(gene_ids) >= 10) {
            go_bp <- tryCatch({ enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE) }, error = function(e) NULL)
            go_mf <- tryCatch({ enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE) }, error = function(e) NULL)
            go_cc <- tryCatch({ enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE) }, error = function(e) NULL)
            kegg_covid <- tryCatch({ enrichKEGG(gene = gene_ids$ENTREZID, organism = "hsa",
                pAdjustMethod = "BH", pvalueCutoff = 0.05) }, error = function(e) NULL)
            
            go_covid <- list(BP = go_bp, MF = go_mf, CC = go_cc)
            
            if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
                ggsave("results/enrichment/COVID_GSE157103_GO_BP_dotplot.pdf",
                       dotplot(go_bp, showCategory = 20), width = 10, height = 8)
            }
            if (!is.null(kegg_covid) && nrow(kegg_covid@result) > 0) {
                ggsave("results/enrichment/COVID_GSE157103_KEGG_dotplot.pdf",
                       dotplot(kegg_covid, showCategory = 20), width = 10, height = 8)
            }
        }
    }
}

# Alzheimer's enrichment
if ("GSE5281" %in% names(all_deg_results)) {
    ad_degs <- all_deg_results$GSE5281 %>% filter(DEG_status != "Not significant") %>% pull(Gene)
    cat("  Alzheimer's DEGs:", length(ad_degs), "\n")
    
    if (length(ad_degs) >= 20) {
        gene_ids <- tryCatch({ bitr(ad_degs, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) },
                             error = function(e) data.frame())
        
        if (nrow(gene_ids) >= 10) {
            go_bp <- tryCatch({ enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE) }, error = function(e) NULL)
            go_mf <- tryCatch({ enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE) }, error = function(e) NULL)
            go_cc <- tryCatch({ enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE) }, error = function(e) NULL)
            kegg_ad <- tryCatch({ enrichKEGG(gene = gene_ids$ENTREZID, organism = "hsa",
                pAdjustMethod = "BH", pvalueCutoff = 0.05) }, error = function(e) NULL)
            
            go_ad <- list(BP = go_bp, MF = go_mf, CC = go_cc)
            
            if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
                ggsave("results/enrichment/AD_GSE5281_GO_BP_dotplot.pdf",
                       dotplot(go_bp, showCategory = 20), width = 10, height = 8)
            }
            if (!is.null(kegg_ad) && nrow(kegg_ad@result) > 0) {
                ggsave("results/enrichment/AD_GSE5281_KEGG_dotplot.pdf",
                       dotplot(kegg_ad, showCategory = 20), width = 10, height = 8)
            }
        }
    }
}

# Save results
enrichment_results <- list(COVID = list(GO = go_covid, KEGG = kegg_covid),
                           AD = list(GO = go_ad, KEGG = kegg_ad))
saveRDS(enrichment_results, "results/enrichment/enrichment_results.rds")

if (!is.null(go_covid$BP)) write.csv(go_covid$BP@result, "results/enrichment/COVID_GO_BP.csv", row.names = FALSE)
if (!is.null(go_ad$BP)) write.csv(go_ad$BP@result, "results/enrichment/AD_GO_BP.csv", row.names = FALSE)
if (!is.null(kegg_covid)) write.csv(kegg_covid@result, "results/enrichment/COVID_KEGG.csv", row.names = FALSE)
if (!is.null(kegg_ad)) write.csv(kegg_ad@result, "results/enrichment/AD_KEGG.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# WGCNA ANALYSIS
# -----------------------------------------------------------------------------

cat("\nWGCNA Analysis...\n")

WGCNA_min_module_size <- 30
WGCNA_merge_cut_height <- 0.25
wgcna_results <- list()

processed_files <- list.files("data/processed", pattern = "_processed.rds", full.names = TRUE)
key_files <- processed_files[sapply(processed_files, function(f) grepl("GSE157103|GSE5281", f))]

for (rds_file in key_files) {
    gse_id <- gsub("_processed.rds", "", basename(rds_file))
    cat("  WGCNA:", gse_id, "\n")
    
    data <- readRDS(rds_file)
    expr_matrix <- data$expression
    pheno_data <- data$phenotype
    
    # Prepare data
    good_genes <- goodSamplesGenes(t(expr_matrix), verbose = 0)
    expr_matrix <- expr_matrix[good_genes$goodGenes, good_genes$goodSamples]
    
    gene_vars <- apply(expr_matrix, 1, var, na.rm = TRUE)
    top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(5000, nrow(expr_matrix))]
    data_expr <- t(expr_matrix[top_genes, ])
    
    # Soft threshold
    powers <- c(1:20)
    sft <- pickSoftThreshold(data_expr, powerVector = powers, verbose = 0, networkType = "signed")
    
    pdf(file.path("results/wgcna", paste0(gse_id, "_soft_threshold.pdf")), width = 12, height = 5)
    par(mfrow = c(1, 2))
    plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
         xlab = "Power", ylab = "Scale Free R^2", main = "Scale independence", type = "n")
    text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, col = "red")
    abline(h = 0.85, col = "red")
    plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Power", ylab = "Mean Connectivity",
         main = "Mean connectivity", type = "n")
    text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
    dev.off()
    
    soft_power <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)
    
    # Build network
    net <- blockwiseModules(data_expr, power = soft_power, networkType = "signed", TOMType = "signed",
        minModuleSize = WGCNA_min_module_size, mergeCutHeight = WGCNA_merge_cut_height,
        numericLabels = TRUE, saveTOMs = TRUE,
        saveTOMFileBase = file.path("results/wgcna", paste0(gse_id, "_TOM")), verbose = 0)
    
    module_colors <- labels2colors(net$colors)
    cat("    Modules:", length(unique(module_colors)) - 1, "\n")
    
    # Dendrogram
    pdf(file.path("results/wgcna", paste0(gse_id, "_dendrogram.pdf")), width = 12, height = 9)
    plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
                        "Module colors", dendroLabels = FALSE, hang = 0.03)
    dev.off()
    
    # Hub genes
    hub_genes_list <- list()
    for (mod in unique(module_colors)[unique(module_colors) != "grey"]) {
        mod_genes <- colnames(data_expr)[module_colors == mod]
        if (length(mod_genes) > 0) {
            connectivity <- softConnectivity(data_expr[, mod_genes], power = 6)
            names(connectivity) <- mod_genes
            top_hubs <- names(sort(connectivity, decreasing = TRUE))[1:min(30, length(mod_genes))]
            hub_genes_list[[mod]] <- data.frame(Module = mod, Gene = top_hubs, Connectivity = connectivity[top_hubs])
        }
    }
    hub_genes_df <- bind_rows(hub_genes_list)
    write.csv(hub_genes_df, file.path("results/wgcna", paste0(gse_id, "_hub_genes.csv")), row.names = FALSE)
    
    # Module-trait correlation
    group_cols <- c("disease.state.ch1", "Disease State:ch1", "disease state:ch1", "condition:ch1")
    group_col <- NA
    for (col in group_cols) if (col %in% colnames(pheno_data)) { group_col <- col; break }
    
    if (!is.na(group_col)) {
        group_values <- pheno_data[[group_col]]
        valid_idx <- !is.na(group_values) & group_values != ""
        common_samples <- intersect(rownames(data_expr), rownames(pheno_data)[valid_idx])
        
        if (length(common_samples) > 10) {
            data_expr_matched <- data_expr[common_samples, ]
            trait_data <- data.frame(Disease = as.numeric(factor(pheno_data[common_samples, group_col])) - 1)
            rownames(trait_data) <- common_samples
            
            MEs <- moduleEigengenes(data_expr_matched, colors = module_colors)$eigengenes
            MEs <- orderMEs(MEs)
            
            module_trait_cor <- cor(MEs, trait_data, use = "p")
            module_trait_pvalue <- corPvalueStudent(module_trait_cor, nrow(data_expr_matched))
            
            text_matrix <- paste(signif(module_trait_cor, 2), "\n(", signif(module_trait_pvalue, 1), ")", sep = "")
            dim(text_matrix) <- dim(module_trait_cor)
            
            pdf(file.path("results/wgcna", paste0(gse_id, "_module_trait_heatmap.pdf")), width = 10, height = 12)
            par(mar = c(6, 8.5, 3, 3))
            labeledHeatmap(Matrix = module_trait_cor, xLabels = colnames(trait_data),
                           yLabels = names(MEs), colors = blueWhiteRed(50), textMatrix = text_matrix,
                           setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1, 1), main = "Module-trait relationships")
            dev.off()
        }
    }
    
    gene_modules <- data.frame(Gene = colnames(data_expr), Module = module_colors)
    write.csv(gene_modules, file.path("results/wgcna", paste0(gse_id, "_gene_modules.csv")), row.names = FALSE)
    
    wgcna_results[[gse_id]] <- list(colors = module_colors, hub_genes = hub_genes_df)
}

saveRDS(wgcna_results, "results/wgcna/wgcna_all_results.rds")

# Cross-disease comparison
shared_results <- readRDS("results/degs/shared_degs_results.rds")
cat("\nCross-disease hub gene overlap:\n")
for (gse_id in names(wgcna_results)) {
    overlap <- intersect(shared_results$shared_all, wgcna_results[[gse_id]]$hub_genes$Gene)
    cat("  ", gse_id, ":", length(overlap), "/", length(shared_results$shared_all), "shared DEGs are hub genes\n")
}

# -----------------------------------------------------------------------------
# COMPLETE
# -----------------------------------------------------------------------------