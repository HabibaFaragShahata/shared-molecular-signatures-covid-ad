# Final Report Files

## Multi-Omics Meta-Analysis: COVID-19 and Alzheimer's Disease Molecular Signatures

### Authors
| Name | Student ID |
|------|------------|
| Tamim Alaa | 211002002 |
| Hana Mohamed Hikal | 231000300 |
| Asmaa Ragab | 221001773 |
| Habiba Farag | 221001888 |

---

## 📊 Main Figures (figures/)

| Figure | Description |
|--------|-------------|
| **Figure1A_COVID_volcano.pdf** | Volcano plot of COVID-19 DEGs (GSE157103) |
| **Figure1B_AD_volcano.pdf** | Volcano plot of Alzheimer's DEGs (GSE5281) |
| **Figure2_Venn_diagram.pdf** | Venn diagram showing overlap between COVID-19 and AD DEGs |
| **Figure3A_COVID_heatmap.pdf** | Heatmap of top COVID-19 DEGs |
| **Figure3B_AD_heatmap.pdf** | Heatmap of top AD DEGs |
| **Figure4A_COVID_GO_BP.pdf** | GO Biological Process enrichment - COVID-19 |
| **Figure4B_AD_GO_BP.pdf** | GO Biological Process enrichment - AD |
| **Figure5A_COVID_KEGG.pdf** | KEGG pathway enrichment - COVID-19 |
| **Figure5B_AD_KEGG.pdf** | KEGG pathway enrichment - AD |
| **Figure6A_COVID_WGCNA_dendrogram.pdf** | WGCNA gene clustering dendrogram - COVID-19 |
| **Figure6B_AD_WGCNA_dendrogram.pdf** | WGCNA gene clustering dendrogram - AD |
| **Figure7A_COVID_module_trait.pdf** | Module-trait correlation heatmap - COVID-19 |
| **Figure7B_AD_module_trait.pdf** | Module-trait correlation heatmap - AD |

---

## 📋 Main Tables (tables/)

| Table | Description |
|-------|-------------|
| **Table1_Shared_DEGs.xlsx** | 20 shared DEGs between COVID-19 and AD with direction analysis |
| **Table2_COVID_Hub_Genes.csv** | 420 hub genes from COVID-19 WGCNA network |
| **Table3_AD_Hub_Genes.csv** | 90 hub genes from AD WGCNA network |
| **Table4_COVID_GO_BP.csv** | Full GO Biological Process results - COVID-19 |
| **Table5_AD_GO_BP.csv** | Full GO Biological Process results - AD (6,175 terms) |
| **Table6_COVID_KEGG.csv** | Full KEGG pathway results - COVID-19 |
| **Table7_AD_KEGG.csv** | Full KEGG pathway results - AD |

---

## 📎 Supplementary Materials (supplementary/)

### Supplementary Figures
| File | Description |
|------|-------------|
| **Supp_Fig_S1_Venn_Upregulated.pdf** | Venn diagram of upregulated genes |
| **Supp_Fig_S2_Venn_Downregulated.pdf** | Venn diagram of downregulated genes |
| **Supp_Fig_S3_COVID_Soft_Threshold.pdf** | WGCNA soft threshold selection - COVID-19 |
| **Supp_Fig_S4_AD_Soft_Threshold.pdf** | WGCNA soft threshold selection - AD |

### Supplementary Tables
| File | Description |
|------|-------------|
| **Supp_Table_S1_COVID_All_DEGs.csv** | Complete list of 231 COVID-19 DEGs |
| **Supp_Table_S2_AD_All_DEGs.csv** | Complete list of 3,654 AD DEGs |
| **Supp_Table_S3_COVID_Gene_Modules.csv** | Gene module assignments - COVID-19 network |
| **Supp_Table_S4_AD_Gene_Modules.csv** | Gene module assignments - AD network |

---

## Key Results Summary

- **COVID-19 DEGs**: 231 genes (118 upregulated, 113 downregulated)
- **AD DEGs**: 3,654 genes (1,411 upregulated, 2,243 downregulated)
- **Shared DEGs**: 20 genes
  - 7 concordantly regulated (same direction)
  - 13 discordantly regulated (opposite direction)
- **COVID-19 WGCNA**: 14 modules, 420 hub genes
- **AD WGCNA**: 3 modules, 90 hub genes
- **Key Finding**: 7/20 shared DEGs are hub genes in the COVID-19 network

---

## Datasets Used

| Dataset | Disease | Platform | Samples |
|---------|---------|----------|---------|
| GSE157103 | COVID-19 | RNA-seq | Blood samples |
| GSE5281 | Alzheimer's | Affymetrix GPL570 | Brain tissue (6 regions) |

---

## Analysis Criteria

- **DEG Threshold**: |log2FC| ≥ 1, adjusted p-value < 0.05
- **WGCNA**: Signed network, minimum module size 30
- **Enrichment**: Benjamini-Hochberg FDR correction

---
