# MicroRNAs regulation role in corneal endothelium affected with Fuchs endothelial corneal dystrophy

Analysis code for miRNA and RNA-seq studies in Fuchs endothelial corneal dystrophy (FECD).

## Project Structure
You can download the HTML files from the directory to explore interactive visualizations and tables.

### 01_QC

This module performs quality control and differential expression analysis of miRNA NanoString data from FECD patients and healthy controls.

**Key steps:**

- **Data preprocessing** – import from nSolver, outlier removal, metadata integration 
- **Quality control** – PCA. Batch effects and age imbalance were detected
- **Normalization comparison** – heatmaps for different normalisation from nSolver methods. Ligation normalization was selected
- **Differential expression** – t-test, Mann–Whitney U-test, and limma. Low-expression miRNAs were filtered. Identified 98 DEGs (mostly up-regulated in FECD)
- **Comparison with published data** – overlap with Matthaei et al. DEGs, logFC correlation analysis
- **CPM validation** – cross-validation using RNA-seq data (our dataset and Chu et al.). Host gene and mirtrone expression correlations
- **miR-29 analysis** – downregulation of miR-29 family in FECD, upregulation of its targets (collagens, LAMC1). Also examined DICER1/DROSHA/DGCR8
- **Annotation export** – tables with HGNC symbols, MIMAT IDs, and 5p/3p strand information for downstream analyses
