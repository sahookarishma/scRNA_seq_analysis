# scRNA_seq_analysis
📂 Dataset

GEO Accession: GSE221575
Platform: Illumina NovaSeq 6000
Samples: Primary tumor, metastatic, and normal tissues from CRC patients.
🧪 Objective

To validate cell-specific expression of candidate hub genes in colorectal cancer using scRNA-seq datasets and identify cluster-specific differentially expressed genes (DEGs).
🔧 Pipeline Overview
1. Data Loading and Seurat Object Creation

    Each dataset folder (tumor, normal, tumor2, normal2, metastasis) contains compressed matrix files:

        matrix.mtx.gz

        features.tsv.gz

        barcodes.tsv.gz

    These are read and converted to Seurat objects using ReadMtx() and CreateSeuratObject().

2. Merging Samples

    All five Seurat objects are merged using merge() with cell ID prefixes to track sample origin.

    Metadata is updated to extract patient/sample information.

3. Quality Control

    Filtering based on:

        nCount_RNA > 800

        nFeature_RNA > 500

        mitoPercent < 10

    Visualization: VlnPlot(), FeatureScatter() plots for QC metrics.

4. Normalization & Feature Selection

    Data is normalized using LogNormalize.

    Highly variable features selected using FindVariableFeatures().

    Scaling performed on all genes (ScaleData()).

5. Dimensionality Reduction & Clustering

    PCA and UMAP used for dimensionality reduction.

    Clustering performed using:

        FindNeighbors()

        FindClusters() (resolution = 0.5)

6. Marker Gene Identification

    FindAllMarkers() used to identify cluster-specific marker genes.

    Violin plots of hub gene expression are generated.

📊 Results
✅ Clustering

    Cells from five patients clustered into 22 distinct populations using UMAP.

🧬 Hub Gene Expression

    Violin plots demonstrate expression of top 11 hub genes across clusters.

    High expression observed for:

        ACTB, RPS27A, H3F3B

        Moderate: STAT3, APP, HIF1A, CTNNB1

    Out of ~20 hub genes examined, 14 were validated across single-cell clusters, excluding IL6, CSF2, IL10, ALB, H6PD, and PXDN.

📌 Visualization

UMAP clusters 

📌 References

    Butler et al. Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nature Biotechnology (2018).

    Stuart et al. Comprehensive Integration of Single-Cell Data. Cell (2019).

    GEO Dataset GSE221575: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221575
   
    
📌 Publication

    Results of this analysis is published in the following paper 
    Sahoo K, Sundararajan V. IL-1β and associated molecules as prognostic biomarkers linked with immune cell infiltration in colorectal cancer: an integrated statistical and machine learning approach. Discov Oncol. 2025 Feb 28;16(1):252. doi: 10.1007/s12672-025-01989-3. PMID: 40019680; PMCID: PMC11871282.

