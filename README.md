# MicroRNAs regulation role in corneal endothelium affected with Fuchs endothelial corneal dystrophy

Analysis code for miRNA and RNA-seq studies in Fuchs endothelial corneal dystrophy (FECD).

## Project Structure
You can download the HTML files from the `html_reports` directory to explore interactive visualizations and tables or explore `.md` files on GitHub.

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

### 02_HGNC_formatting

Data preparation script. It performs identifier mapping and annotation formatting that is necessary for downstream analyses but is not intended for direct review.

**What happens here:**

- Standardisation of miRNA names to HGNC symbols (e.g., `LET7A` → `MIRLET7A`)
- Merging chip probe sets with Gencode GTF annotation (genes, coordinates, strands)
- Creating BED files for genomic intersections
- Mapping MIMAT IDs to miRBase GFF3 entries
- Identifying mirtrons (miRNAs located within host genes) and nearby genomic regions

### 03_target_correlation

- Compute Spearman correlations (miRNA vs gene) separately for Control and FECD
- Keep negatively correlated pairs (≤ –0.7) as candidate miRNA–target interactions
- Filter by miRTarBase v10 (experimentally supported targets)
- Annotate with Gencode (HGNC, strands, coordinates) and MIMAT (mirBase)
- Venn diagrams + boxplots for expression validation

### 04_KEGG pathway search

KEGG enrichment analysis for miRNA targets (DEGs and miR-29). Outputs top pathways, frequency of gene occurrence, plots with `pathview`.
