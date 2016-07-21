# awesome-single-cell

List of software packages for single-cell data analysis, including RNA-seq, ATAC-seq, etc. [Contributions welcome](https://github.com/seandavi/awesome-single-cell/blob/master/CONTRIBUTE.md)...

## Software packages

### RNA-seq

- [BASiCs](https://github.com/catavallejos/BASiCS) - [R] - Bayesian Analysis of single-cell RNA-seq data. Estimates cell-specific normalization constants. Technical variability is quantified based on spike-in genes. The total variability of the expression counts is decomposed into technical and biological components.
- [BPSC](https://github.com/nghiavtr/BPSC) - [R] - Beta-Poisson model for single-cell RNA-seq data analyses
- [Cellity](https://github.com/teichlab/cellity) - [R] - Classification of low quality cells in scRNA-seq data using R
- [clusterExperiment](https://github.com/epurdom/clusterExperiment) - [R] - Functions for running and comparing many different clusterings of single-cell sequencing data. Meant to work with SCONE and slingshot.
- [ECLAIR](https://github.com/GGiecold/ECLAIR) - [python] - ECLAIR stands for Ensemble Clustering for Lineage Analysis, Inference and Robustness. Robust and scalable inference of cell lineages from gene expression data.
- [Falco](https://github.com/VCCRI/Falco/) - [AWS cloud] - [Falco: A quick and flexible single-cell RNA-seq processing framework on the cloud](http://www.biorxiv.org/content/early/2016/07/15/064006.abstract).
- [HocusPocus](https://github.com/joeburns06/hocuspocus) - [R] - Basic PCA-based workflow for analysis and plotting of single cell RNA-seq data.
- [MAST](https://github.com/RGLab/MAST) - [R] - Model-based Analysis of Single-cell Transcriptomics
(MAST) fits a two-part, generalized linear models that are specially adapted for bimodal and/or zero-inflated single cell gene expression data.
- [Monocle](http://cole-trapnell-lab.github.io/monocle-release/) - [R] - Differential expression and time-series analysis for single-cell RNA-Seq.
- [OEFinder](https://github.com/lengning/OEFinder) - [R] - Identify ordering effect genes in single cell RNA-seq data. OEFinder shiny impelemention depends on packages shiny, shinyFiles, gdata, and EBSeq.
- [SC3](https://github.com/hemberg-lab/sc3) - [R] - SC3 is an interactive tool for the unsupervised clustering of cells from single cell RNA-Seq experiments.
- [scater](bioconductor.org/packages/scater) - [R] - Scater places an emphasis on tools for quality control, visualisation and pre-processing of data before further downstream analysis, filling a useful niche between raw RNA-sequencing count or transcripts-per-million data and more focused downstream modelling tools such as monocle, scLVM, SCDE, edgeR, limma and so on.
- [SCDE](https://github.com/hms-dbmi/scde) - [R] - Differential expression using error models and overdispersion-based identification of important gene sets.
- [SCell](https://github.com/diazlab/SCell) - [matlab] - SCell is an integrated software tool for quality filtering, normalization, feature selection, iterative dimensionality reduction, clustering and the estimation of gene-expression gradients from large ensembles of single-cell RNA-seq datasets. SCell is open source, and implemented with an intuitive graphical interface.
- [scLVM](https://github.com/PMBio/scLVM) - [R] - scLVM is a modelling framework for single-cell RNA-seq data that can be used to dissect the observed heterogeneity into different sources, thereby allowing for the correction of confounding sources of variation.  scLVM was primarily designed to account for cell-cycle induced variations in single-cell RNA-seq data where cell cycle is the primary soure of variability. 
- [SCONE](https://github.com/YosefLab/scone) - [R] - SCONE (Single-Cell Overview of Normalized Expression), a package for single-cell RNA-seq data quality control (QC) and normalization. This data-driven framework uses summaries of expression data to assess the efficacy of normalization workflows.
- [SCOUP](https://github.com/hmatsu1226/SCOUP) - [C++] - Uses probabilistic model based on the Ornstein-Uhlenbeck process to analyze single-cell expression data during differentiation. 
- [scran](http://bioconductor.org/packages/scran) - [R] - This package implements a variety of low-level analyses of single-cell RNA-seq data. Methods are provided for normalization of cell-specific biases, pool-based norms to estimate size factors, assignment of cell cycle phase, and detection of highly variable and significantly correlated genes.
- [SCRAT](https://github.com/zji90/SCRAT) - [R] - SCRAT provides essential tools for users to read in single-cell regolome data (ChIP-seq, ATAC-seq, DNase-seq) and summarize into different types of features. It also allows users to visualize the features, cluster samples and identify key features.
- [SEPA](https://github.com/zji90/SEPA) - [R] - SEPA provides convenient functions for users to assign genes into different gene expression patterns such as constant, monotone increasing and increasing then decreasing. SEPA then performs GO enrichment analysis to analysis the functional roles of genes with same or similar patterns.
- [scTCRseq](https://github.com/ElementoLab/scTCRseq) - [python] - Map T-cell receptor (TCR) repertoires from single cell RNAseq.
- [Seurat](http://www.satijalab.org/seurat.html) - [R] - It contains easy-to-use implementations of commonly used analytical techniques, including the identification of highly variable genes, dimensionality reduction (PCA, ICA, t-SNE), standard unsupervised clustering algorithms (density clustering, hierarchical clustering, k-means), and the discovery of differentially expressed genes and markers.
- [sincell](http://bioconductor.org/packages/sincell) - [R] - Existing computational approaches for the assessment of cell-state hierarchies from single-cell data might be formalized under a general workflow composed of i) a metric to assess cell-to-cell similarities (combined or not with a dimensionality reduction step), and ii) a graph-building algorithm (optionally making use of a cells-clustering step). Sincell R package implements a methodological toolbox allowing flexible workflows under such framework.
- [sincera](https://research.cchmc.org/pbge/sincera.html) - [R] - R-based pipeline for single-cell analysis including clustering and visualization.
- [SinQC](http://www.morgridge.net/SinQC.html) - [R] - A Method and Tool to Control Single-cell RNA-seq Data Quality.
- [SLICER](https://github.com/jw156605/SLICER) - [R] - Selective Locally linear Inference of Cellular Expression Relationships (SLICER) algorithm for inferring cell trajectories.
- [slingshot](https://github.com/kstreet13/slingshot) - [R] - Functions for identifying and characterizing continuous developmental trajectories in single-cell sequencing data.
- [SPADE](http://www.nature.com/nprot/journal/v11/n7/full/nprot.2016.066.html) - [R] - Visualization and cellular hierarchy inference of single-cell data using SPADE.
- [TraCeR](http://github.com/teichlab/tracer) - [python] - Reconstruction of T-Cell receptor sequences from single-cell RNA-seq data.
- [TSCAN](https://github.com/zji90/TSCAN) - [R] - Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis.

### Copy number analysis

- [Gingko](https://github.com/robertaboukhalil/ginkgo) - [R, C] - Gingko is a cloud-based web app single-cell copy-number variation analysis tool.

## Tutorials and workflows
 
- [Aaron Lun's Single Cell workflow on Bioconductor](http://bioconductor.org/help/workflows/simpleSingleCell/) - [R] - This article describes a computational workflow for basic analysis of scRNA-seq data using software packages from the open-source Bioconductor project.
- [Bioconductor2016 Single-cell-RNA-sequencing workshop by Sandrine Dudoit lab](https://github.com/drisso/bioc2016singlecell) - [R] - SCONE, clusterExperiment, and slingshot tutorial.
- [BiomedCentral Single Cell Omics collectin](http://www.biomedcentral.com/collections/singlecellomics) - collection of papers describing techniques for single-cell analysis and protocols.
- [Gilad Lab Single Cell Data Exploration](http://jdblischak.github.io/singleCellSeq/analysis/) - R-based exploration of single cell sequence data. Lots of experimentation. 
- [Harvard STEM Cell Institute Single Cell Workshop 2015](http://hms-dbmi.github.io/scw/) - workshop on common computational analysis techniques for scRNA-seq data from differential expression to subpopulation identification and network analysis. [See course description for more information](http://scholar.harvard.edu/jeanfan/classes/single-cell-workshop-2015) 
- [Hemberg Lab scRNA-seq course materials](http://hemberg-lab.github.io/scRNA.seq.course/index.html) 
- [Using Seurat for unsupervised clustering and biomarker discovery](http://www.satijalab.org/seurat-intro.html) - 301 single cells across diverse tissues from (Pollen et al., Nature Biotechnology, 2014)
- [Using Seurat for spatial inference in single-cell data](http://www.satijalab.org/seurat-intro.html) - 851 single cells from Zebrafish embryogenesis (Satija*, Farrell* et al., Nature Biotechnology, 2015)

## Experimental design

- [Design and computational analysis of single-cell RNA-sequencing experiments](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y)

## Similar lists and collections

- [CrazyHotTommy's RNA-seq analysis list](https://github.com/crazyhottommy/RNA-seq-analysis#single-cell-rna-seq) - Very broad list that includes some single cell RNA-seq packages and papers.
