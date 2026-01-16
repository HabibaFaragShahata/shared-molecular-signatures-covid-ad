# Multi-Omics Meta-Analysis of Shared Molecular Signatures Between SARS-CoV-2 Infection and Alzheimer's Disease
  PROJECT RECORD LINK
     https://drive.google.com/file/d/17wEaKsDUtfJ6_FJS9el6yK4wIvAAiJKg/view?usp=drive_link 

## CBIO310 Bioinformatics Research Report

**Analysis Pipeline:** R/Bioconductor  

### Authors
| Name | Student ID |
|------|------------|
| Tamim Alaa | 211002002 |
| Hana Mohamed Hikal | 231000300 |
| Asmaa Ragab | 221001773 |
| Habiba Farag | 221001888 |

---

## Abstract

This study investigates the molecular relationship between SARS-CoV-2 infection and Alzheimer's Disease (AD) through comprehensive bioinformatics analysis of publicly available transcriptomic data. Using differential gene expression analysis, functional enrichment, and weighted gene co-expression network analysis (WGCNA), we identified **20 shared differentially expressed genes (DEGs)** between COVID-19 and AD, with **7 genes showing concordant regulation** (same direction of change in both diseases). Our findings reveal shared biological pathways centered on mitochondrial dysfunction, oxidative phosphorylation, and neuroinflammation, providing molecular evidence for the neurological sequelae observed in COVID-19 patients.

---

## 1. Introduction

### 1.1 Background

The COVID-19 pandemic, caused by SARS-CoV-2, has revealed significant neurological complications in infected patients, including cognitive impairment, memory deficits, and accelerated neurodegeneration. Emerging clinical evidence suggests a potential link between SARS-CoV-2 infection and increased risk of Alzheimer's Disease (AD), though the underlying molecular mechanisms remain poorly understood.

### 1.2 Research Objectives

1. Identify differentially expressed genes (DEGs) in COVID-19 and Alzheimer's Disease datasets
2. Discover shared molecular signatures between the two conditions
3. Characterize common biological pathways through functional enrichment analysis
4. Identify hub genes and co-expression modules using WGCNA

### 1.3 Clinical Significance

Understanding the molecular connections between viral infection and neurodegeneration has important implications for:
- Identifying patients at risk for post-COVID cognitive decline
- Developing therapeutic strategies targeting shared pathways
- Elucidating mechanisms of infection-induced neurodegeneration

---

## 2. Materials and Methods

### 2.1 Data Sources

Gene expression data was obtained from the NCBI Gene Expression Omnibus (GEO):

| Dataset | Disease | Platform | Tissue | Samples |
|---------|---------|----------|--------|---------|
| GSE157103 | COVID-19 | RNA-seq | Blood | 126 (COVID+ vs Healthy) |
| GSE5281 | Alzheimer's | Affymetrix HG-U133 Plus 2.0 | Brain (6 regions) | 161 (AD vs Normal) |

### 2.2 Analysis Pipeline

```
Pipeline.R     → All Analysis
```

### 2.3 DEG Identification Criteria

- **Fold Change Threshold:** |log₂FC| ≥ 1
- **Statistical Significance:** Adjusted p-value < 0.05 (Benjamini-Hochberg correction)
- **Method:** limma (linear models for microarray data)

### 2.4 Functional Enrichment Analysis

- **Gene Ontology (GO):** Biological Process, Molecular Function, Cellular Component
- **KEGG Pathway Analysis:** Disease-related pathways
- **Tools:** clusterProfiler, org.Hs.eg.db

### 2.5 WGCNA Parameters

- **Network Type:** Signed
- **Minimum Module Size:** 30 genes
- **Merge Cut Height:** 0.25
- **Top Variable Genes:** 5,000

---

## 3. Results

### 3.1 Differential Expression Analysis

#### 3.1.1 COVID-19 DEGs (GSE157103)

| Category | Count |
|----------|-------|
| Total DEGs | 231 |
| Upregulated | 118 |
| Downregulated | 113 |

**Top Upregulated Genes:** Immune response, cytokine signaling, inflammatory pathways  
**Top Downregulated Genes:** Metabolic processes, cellular homeostasis

#### 3.1.2 Alzheimer's Disease DEGs (GSE5281)

| Category | Count |
|----------|-------|
| Total DEGs | 3,654 |
| Upregulated | 1,411 |
| Downregulated | 2,243 |

**Top Upregulated Genes:** Inflammatory response, oxidative stress  
**Top Downregulated Genes:** Synaptic transmission, mitochondrial function, energy metabolism

### 3.2 Shared DEGs Between COVID-19 and AD

**Total Shared DEGs: 20**

| Regulation Pattern | Count | Interpretation |
|-------------------|-------|----------------|
| Concordant Upregulated | 3 | Same activation in both diseases |
| Concordant Downregulated | 4 | Same suppression in both diseases |
| **Total Concordant** | **7** | **Shared molecular mechanisms** |
| Discordant (Up COVID, Down AD) | 11 | Disease-specific responses |
| Discordant (Down COVID, Up AD) | 2 | Compensatory mechanisms |

#### Key Concordant Genes (Both Diseases, Same Direction)

The 7 concordantly regulated genes represent the strongest evidence for shared molecular mechanisms between COVID-19 and Alzheimer's Disease. These genes are dysregulated in the same direction in both conditions, suggesting common pathological pathways.

### 3.3 Functional Enrichment Analysis

#### 3.3.1 COVID-19 Pathway Enrichment

**GO Biological Process - Top Terms:**
1. Inflammatory response
2. Immune system activation
3. Cytokine-mediated signaling
4. Defense response to virus
5. Neutrophil activation

**KEGG Pathways Enriched:**
- Cytokine-cytokine receptor interaction
- TNF signaling pathway
- NF-κB signaling pathway
- COVID-19 (hsa05171)

#### 3.3.2 Alzheimer's Disease Pathway Enrichment

**GO Biological Process - Top Terms (6,175 total):**
1. Proton transmembrane transport (p = 3.47×10⁻¹⁴)
2. Aerobic respiration (p = 1.50×10⁻¹¹)
3. Oxidative phosphorylation (p = 6.45×10⁻¹¹)
4. Regulation of neuron projection development (p = 9.01×10⁻¹¹)
5. Proton motive force-driven ATP synthesis (p = 1.09×10⁻⁹)

**KEGG Pathways Enriched (Key Pathways):**
- Oxidative phosphorylation
- Pathways of neurodegeneration (hsa05022)
- Alzheimer's disease (hsa05010)
- Parkinson's disease
- Huntington's disease

#### 3.3.3 Shared Pathway Themes

The enrichment analysis reveals convergent biology:

| Shared Theme | COVID-19 | Alzheimer's |
|--------------|----------|-------------|
| Mitochondrial Dysfunction | Inflammatory stress | Energy metabolism failure |
| Oxidative Stress | Immune-mediated ROS | Neuronal oxidative damage |
| Inflammatory Signaling | Acute inflammation | Chronic neuroinflammation |
| Cellular Stress Response | Cytokine storm | ER stress, protein misfolding |

### 3.4 WGCNA Co-expression Network Analysis

#### 3.4.1 Network Statistics

| Dataset | Modules | Hub Genes | Soft Power |
|---------|---------|-----------|------------|
| GSE157103 (COVID) | 14 | 420 | Auto-detected |
| GSE5281 (AD) | 3 | 90 | Auto-detected |

#### 3.4.2 Shared DEGs as Hub Genes

**Critical Finding:** 7 out of 20 shared DEGs are hub genes in the COVID-19 network

| Module | Hub Genes (Shared DEGs) |
|--------|------------------------|
| Black | 3 |
| Blue | 1 |
| Magenta | 2 |
| Tan | 1 |

This indicates that shared DEGs occupy central positions in the COVID-19 co-expression network, suggesting they may be key drivers of the molecular connection between the two diseases.

#### 3.4.3 Module-Trait Correlations

Modules showing significant correlation with disease status were identified in both datasets, providing network-level validation of the differential expression findings.

---

## 4. Discussion

### 4.1 Key Findings

1. **Molecular Link Established:** We identified 20 genes that are differentially expressed in both COVID-19 and Alzheimer's Disease, with 7 showing concordant regulation.

2. **Shared Pathological Mechanisms:** Both diseases converge on:
   - Mitochondrial dysfunction and energy metabolism impairment
   - Oxidative stress and reactive oxygen species generation
   - Inflammatory signaling cascades
   - Cellular stress responses

3. **Network Centrality:** Shared DEGs occupy hub positions in disease networks, suggesting they are not peripheral players but central to disease pathophysiology.

### 4.2 Biological Interpretation

The convergence on mitochondrial and oxidative stress pathways is particularly significant:

- **COVID-19:** SARS-CoV-2 infection triggers massive inflammatory responses that generate oxidative stress and mitochondrial dysfunction in multiple organs, including the brain.

- **Alzheimer's Disease:** AD is characterized by progressive mitochondrial dysfunction, energy metabolism failure, and accumulation of oxidative damage in neurons.

The shared genes may represent a common vulnerability pathway where:
1. Viral infection induces acute mitochondrial stress
2. This stress accelerates or triggers neurodegenerative processes
3. Pre-existing AD pathology may increase susceptibility to severe COVID-19

### 4.3 Clinical Implications

1. **Risk Stratification:** Patients with elevated expression of shared genes may be at higher risk for post-COVID cognitive decline.

2. **Therapeutic Targets:** Shared pathways (particularly mitochondrial function and oxidative stress) represent potential intervention points.

3. **Monitoring Biomarkers:** The concordant genes could serve as biomarkers for tracking neurological involvement in COVID-19.

### 4.4 Limitations

1. **Dataset Heterogeneity:** Different tissue sources (blood vs. brain) may limit direct comparisons.
2. **Sample Size:** Limited number of shared DEGs (20) restricts statistical power for enrichment of shared gene set.
3. **Cross-sectional Design:** Cannot establish causality or temporal relationships.
4. **Platform Differences:** RNA-seq vs. microarray may introduce technical variability.

### 4.5 Future Directions

1. Validate shared genes in independent cohorts
2. Perform longitudinal studies of COVID-19 patients with cognitive follow-up
3. Investigate protein-level changes corresponding to transcriptomic findings
4. Explore therapeutic targeting of shared pathways

---

## 5. Conclusions

This multi-omics meta-analysis provides molecular evidence for the clinical observations linking SARS-CoV-2 infection to neurological complications and potential neurodegeneration. The identification of 20 shared DEGs between COVID-19 and Alzheimer's Disease, with 7 showing concordant regulation, establishes a molecular foundation for understanding post-COVID cognitive sequelae.

The convergence on mitochondrial dysfunction and oxidative stress pathways suggests that:
1. COVID-19 may accelerate neurodegenerative processes in susceptible individuals
2. Shared molecular mechanisms exist between acute viral infection and chronic neurodegeneration
3. Therapeutic strategies targeting mitochondrial function may benefit both conditions

These findings contribute to our understanding of the long-term neurological impacts of COVID-19 and highlight the importance of monitoring cognitive function in recovered patients.

---

## 6. References

1. **Primary Methodology:** Rahman MR, et al. (2020). Identification of molecular signatures and pathways to identify novel therapeutic targets in Alzheimer's disease: Insights from a systems biomedicine perspective. *Translational Psychiatry*. DOI: 10.1038/s41398-020-01151-3

2. **GEO Datasets:**
   - GSE157103: Overmyer KA, et al. (2020). Large-Scale Multi-omic Analysis of COVID-19 Severity.
   - GSE5281: Liang WS, et al. (2008). Alzheimer's disease is associated with reduced expression of energy metabolism genes in posterior cingulate neurons.

3. **COVID-19 Neurological Impact:**
   - Taquet M, et al. (2021). 6-month neurological and psychiatric outcomes in 236,379 survivors of COVID-19. *Lancet Psychiatry*.
   - Boldrini M, et al. (2021). How COVID-19 Affects the Brain. *JAMA Psychiatry*.

4. **Bioinformatics Tools:**
   - Ritchie ME, et al. (2015). limma powers differential expression analyses. *Nucleic Acids Research*.
   - Yu G, et al. (2012). clusterProfiler: an R package for comparing biological themes. *OMICS*.
   - Langfelder P & Horvath S (2008). WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*.
  
 

---

**Pipeline:** R 4.4.2 / Bioconductor  
**Packages:** limma, WGCNA, clusterProfiler, GEOquery, org.Hs.eg.db

