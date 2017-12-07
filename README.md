# awesome-single-cell

List of software packages (and the people developing these methods) for single-cell data analysis, including RNA-seq, ATAC-seq, etc. [Contributions welcome](https://github.com/seandavi/awesome-single-cell/blob/master/CONTRIBUTING.md)...


## Software packages

### RNA-seq

- [anchor](https://github.com/yeolab/anchor) - [Python] - ⚓ Find bimodal, unimodal, and multimodal features in your data
- [ascend](https://github.com/IMB-Computational-Genomics-Lab/ascend) - [R] - ascend is an R package comprised of fast, streamlined analysis functions optimized to address the statistical challenges of single cell RNA-seq. The package incorporates novel and established methods to provide a flexible framework to perform filtering, quality control, normalization, dimension reduction, clustering, differential expression and a wide-range of plotting. 
- [BackSPIN](https://github.com/linnarsson-lab/BackSPIN) - [Python] - Biclustering algorithm developed taking into account intrinsic features of single-cell RNA-seq experiments.
- [BASiCS](https://github.com/catavallejos/BASiCS) - [R] - Bayesian Analysis of single-cell RNA-seq data. Estimates cell-specific normalization constants. Technical variability is quantified based on spike-in genes. The total variability of the expression counts is decomposed into technical and biological components. BASiCS can also identify genes with differential expression/over-dispersion between two or more groups of cells.
- [BatchEffectRemoval](https://github.com/ushaham/BatchEffectRemoval) - [Python] - [Removal of Batch Effects using Distribution-Matching Residual Networks](https://doi.org/10.1093/bioinformatics/btx196)
- [BEARscc](https://bitbucket.org/bsblabludwig/bearscc) - [R] - BEARscc makes use of ERCC spike-in measurements to model technical variance as a function of gene expression and technical dropout effects on lowly expressed genes.
- [bonvoyage](https://github.com/yeolab/bonvoyage) - [Python] - 📐 Transform percentage-based units into a 2d space to evaluate changes in distribution with both magnitude and direction.
- [BPSC](https://github.com/nghiavtr/BPSC) - [R] - Beta-Poisson model for single-cell RNA-seq data analyses
- [CellCNN](https://github.com/eiriniar/CellCnn) - [Python] - Representation Learning for detection of phenotype-associated cell subsets
- [Cellity](https://github.com/teichlab/cellity) - [R] - Classification of low quality cells in scRNA-seq data using R
- [cellTree](https://www.bioconductor.org/packages/3.3/bioc/html/cellTree.html) - [R] - Cell population analysis and visualization from single cell RNA-seq data using a Latent Dirichlet Allocation model.
- [clusterExperiment](https://github.com/epurdom/clusterExperiment) - [R] - Functions for running and comparing many different clusterings of single-cell sequencing data. Meant to work with SCONE and slingshot.
- [CytoGuide](https://cyteguide.cytosplore.org/) - [C++,D3] - [CyteGuide: Visual Guidance for Hierarchical Single-Cell Analysis](http://ieeexplore.ieee.org/document/8017575/)
- [DECENT](https://github.com/cz-ye/DECENT) - [R] - The unique features of scRNA-seq data have led to the development of novel methods for differential expression (DE) analysis. However, few of the existing DE methods for scRNA-seq data estimate the number of molecules pre-dropout and therefore do not explicitly distinguish technical and biological zeroes. We develop DECENT, a DE method for scRNA-seq data that adjusts for the imperfect capture efficiency by estimating the number of molecules pre-dropout. 
- [DESCEND](https://github.com/jingshuw/descend) - [R] - DESCEND deconvolves the true gene expression distribution across cells for UMI scRNA-seq counts. It provides estimates of several distribution based statistics (five distribution measurements and the coefficients of covariates (such as batches or cell size)).
- [destiny](http://bioconductor.org/packages/destiny/) - [R] - Diffusion maps are spectral method for non-linear dimension reduction introduced by Coifman et al.(2005). Diffusion maps are based on a distance metric (diffusion distance) which is conceptually relevant to how differentiating cells follow noisy diffusion-like dynamics, moving from a pluripotent state towards more differentiated states.
- [DeLorean](https://cran.r-project.org/web/packages/DeLorean/index.html) - [R] - Bayesian pseudotime estimation algorithm that uses Gaussian processes to model gene expression profiles and provides a full posterior for the pseudotimes.
- [dropClust](https://github.com/debsin/dropClust) - [R/Python] - Efficient clustering of ultra-large scRNA-seq data.
- [ECLAIR](https://github.com/GGiecold/ECLAIR) - [python] - ECLAIR stands for Ensemble Clustering for Lineage Analysis, Inference and Robustness. Robust and scalable inference of cell lineages from gene expression data.
- [embeddr](https://github.com/kieranrcampbell/embeddr) - [R] - Embeddr creates a reduced dimensional representation of the gene space using a high-variance gene correlation graph and laplacian eigenmaps. It then fits a smooth pseudotime trajectory using principal curves.
- [Falco](https://github.com/VCCRI/Falco/) - [AWS cloud] - [Falco: A quick and flexible single-cell RNA-seq processing framework on the cloud](http://www.biorxiv.org/content/early/2016/07/15/064006.abstract).
- [FastProject](https://github.com/yoseflab/fastproject) - [Python] - Signature analysis on low-dimensional projections of single-cell expression data.
- [flotilla](https://github.com/yeolab/flotilla) - [Python] - Reproducible machine learning analysis of gene expression and alternative splicing data
- [GPfates](https://github.com/Teichlab/GPfates) - [Python] - Model transcriptional cell fates as mixtures of Gaussian Processes
- [GiniClust](https://github.com/lanjiangboston/GiniClust) - [Python/R] - GiniClust is a clustering method implemented in Python and R for detecting rare cell-types from large-scale single-cell gene expression data.  GiniClust can be applied to datasets originating from different platforms, such as multiplex qPCR data, traditional single-cell RNAseq or newly emerging UMI-based single-cell RNAseq, e.g. inDrops and Drop-seq.
- [HocusPocus](https://github.com/joeburns06/hocuspocus) - [R] - Basic PCA-based workflow for analysis and plotting of single cell RNA-seq data.
- [ICGS](https://github.com/nsalomonis/altanalyze) - [Python] - Iterative Clustering and Guide-gene Selection (Olsson et al. Nature 2016). Identify discrete, transitional and mixed-lineage states from diverse single-cell transcriptomics platforms. Integrated FASTQ pseudoalignment /quantification (Kallisto), differential expression, cell-type prediction and optional cell cycle exclusion analyses. Specialized methods for processing BAM and 10X Genomics spares matrix files. Associated single-cell splicing PSI methods (MultIPath-PSI). Apart of the AltAnalyze toolkit along with accompanying visualization methods (e.g., heatmap, t-SNE, SashimiPlots, network graphs). Easy-to-use graphical user and commandline interfaces.
- [MAGIC](https://github.com/pkathail/magic) - [python or matlab] - Markov Affinity-based Graph Imputation of Cells (MAGIC).
- [MAST](https://github.com/RGLab/MAST) - [R] - Model-based Analysis of Single-cell Transcriptomics (MAST) fits a two-part, generalized linear models that are specially adapted for bimodal and/or zero-inflated single cell gene expression data.
- [mfa](https://github.com/kieranrcampbell/mfa) - [R] - Bayesian modelling of bifurcations using a mixture of factor analysers
- [K-Branches](https://github.com/theislab/kbranches) - [R] - The main idea behind the K-Branches method is to identify regions of interest (branching regions and tips) in differentiation trajectories of single cells. So far, K-Branches is intended to be used on the diffusion map representation of the data, so the user should either provide the data in diffusion map space or use the destiny package perform diffusion map dimensionality reduction.
- [M3Drop](https://github.com/tallulandrews/M3Drop) - [R] - Michaelis-Menten Modelling of Dropouts for scRNASeq.
- [MAST](https://github.com/RGLab/MAST) - [R] - Model-based Analysis of Single-cell Transcriptomics (MAST) fits a two-part, generalized linear models that are specially adapted for bimodal and/or zero-inflated single cell gene expression data
- [MIMOSCA](https://github.com/asncd/MIMOSCA) - [python] - A repository for the design and analysis of pooled single cell RNA-seq perturbation experiments (Perturb-seq).
- [Monocle](http://cole-trapnell-lab.github.io/monocle-release/) - [R] - Differential expression and time-series analysis for single-cell RNA-Seq.
- [NetworkInference](https://github.com/Tchanders/NetworkInference.jl) - [Julia] - Fast implementation of single-cell network inference algorithms: [Gene Regulatory Network Inference from Single-Cell Data Using Multivariate Information Measures]( http://www.cell.com/cell-systems/fulltext/S2405-4712(17)30386-1 )
- [nimfa](https://github.com/ccshao/nimfa) - [Python] - Nimfa is a Python scripting library which includes a number of published matrix factorization algorithms, initialization methods, quality and performance measures and facilitates the combination of these to produce new strategies. The library represents a unified and efficient interface to matrix factorization algorithms and methods.
- [OEFinder](https://github.com/lengning/OEFinder) - [R] - Identify ordering effect genes in single cell RNA-seq data. OEFinder shiny impelemention depends on packages shiny, shinyFiles, gdata, and EBSeq.
- [OncoNEM](https://bitbucket.org/edith_ross/onconem/src) - [R] -  OncoNEM is a probabilistic method for inferring intra-tumor evolutionarylineage trees from somatic single nucleotide variants of single cells. OncoNEM identifies homogeneous cellularsubpopulations and infers their genotypes as well as a tree describing their evolutionary relationships.
- [Ouija](https://github.com/kieranrcampbell/ouija) - [R] - Incorporate prior information into single-cell trajectory (pseudotime) analyses using Bayesian nonlinear factor analysis.
- [outrigger](https://github.com/YeoLab/outrigger) - [Python] - Outrigger is a program to calculate alternative splicing scores of RNA-Seq data based on junction reads and a *de novo*, custom annotation created with a graph database, especially made for single-cell analyses.
- [pcaReduce](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0984-y) - [R] - hierarchical clustering of single cell transcriptional profiles.
- [PhenoPath](https://github.com/kieranrcampbell/phenopath) - [R] - Single-cell pseudotime with heterogeneous genetic and environmental backgrounds, including Bayesian significance testing of iteractions.
- [PoissonUMIs](https://github.com/tallulandrews/PoissonUMIs) - [R] - Poisson Modelling of scRNASeq UMI counts.
- [SAVER](https://github.com/mohuangx/SAVER) - [R] - SAVER (Single-cell Analysis Via Expression Recovery) implements a regularized regression prediction and empirical Bayes method to recover the true gene expression profile in noisy and sparse single-cell RNA-seq data.
- [SAKE](https://github.com/naikai/sake) - [R] - Single-cell RNA-Seq Analysis and Clustering Evaluation.
- [SC3](https://github.com/hemberg-lab/sc3) - [R] - SC3 is a tool for the unsupervised clustering of cells from single cell RNA-Seq experiments.
- [Scanpy](https://github.com/theislab/scanpy) - [Py] - Scanpy provides computationally efficient tools that scale up to very large data sets and enables simple integraton of advanced machine learning algorithms.
- [scater](bioconductor.org/packages/scater) - [R] - Scater places an emphasis on tools for quality control, visualisation and pre-processing of data before further downstream analysis, filling a useful niche between raw RNA-sequencing count or transcripts-per-million data and more focused downstream modelling tools such as monocle, scLVM, SCDE, edgeR, limma and so on.
- [scDD](https://github.com/kdkorthauer/scDD) - [R] - scDD (Single-Cell Differential Distributions) is a framework to identify genes with different expression patterns between biological groups of interest.  In addition to traditional differential expression, it can detect differences that are more complex and subtle than a mean shift.
- [SCDE](https://github.com/hms-dbmi/scde) - [R] - Differential expression using error models and overdispersion-based identification of important gene sets.
- [SCell](https://github.com/diazlab/SCell) - [matlab] - SCell is an integrated software tool for quality filtering, normalization, feature selection, iterative dimensionality reduction, clustering and the estimation of gene-expression gradients from large ensembles of single-cell RNA-seq datasets. SCell is open source, and implemented with an intuitive graphical interface.
- [SCIMITAR](https://github.com/dimenwarper/scimitar) - [Python] - Single Cell Inference of Morphing Trajectories and their Associated Regulation module (SCIMITAR) is a method for inferring biological properties from a pseudotemporal ordering. It can also be used to obtain progression-associated genes that vary along the trajectory, and genes that change their correlation structure over the trajectory; progression co-associated genes.
- [scImput](https://github.com/Vivianstats/scImpute) - [R] - [scImpute: Accurate And Robust Imputation For Single Cell RNA-Seq Data](https://doi.org/10.1101/141598)
- [SCINIC](https://gbiomed.kuleuven.be/english/research/50000622/lcb/tools/scenic) - [R] - [SCENIC: single-cell regulatory network inference and clustering](https://www.biorxiv.org/content/early/2017/05/31/144501)
- [scvis](https://bitbucket.org/jerry00/scvis-dev) - [python] - [Interpretable dimensionality reduction of single cell transcriptome data with deep generative models](https://doi.org/10.1101/178624)
- [scLVM](https://github.com/PMBio/scLVM) - [R] - scLVM is a modelling framework for single-cell RNA-seq data that can be used to dissect the observed heterogeneity into different sources, thereby allowing for the correction of confounding sources of variation.  scLVM was primarily designed to account for cell-cycle induced variations in single-cell RNA-seq data where cell cycle is the primary soure of variability.
- [scTDA](https://github.com/RabadanLab/scTDA) - [Python] - scTDA is an object oriented python library for topological data analysis of high-throughput single-cell RNA-seq data. It includes tools for the preprocessing, analysis, and exploration of single-cell RNA-seq data based on topological representations.
- [scmap](http://bioconductor.org/packages/scmap) - [R] - scmap is a method for projecting cells from a scRNA-seq experiment on to the cell-types identified in a different experiment.
- [SCnorm](https://github.com/rhondabacher/SCnorm) - [R] - A quantile regression based approach for robust normalization of single cell RNA-seq data.
- [SCONE](https://github.com/YosefLab/scone) - [R] - SCONE (Single-Cell Overview of Normalized Expression), a package for single-cell RNA-seq data quality control (QC) and normalization. This data-driven framework uses summaries of expression data to assess the efficacy of normalization workflows.
- [SCORPIUS](https://github.com/rcannood/SCORPIUS) - [R] - SCORPIUS an unsupervised approach for inferring developmental chronologies from single-cell RNA sequencing data. It accurately reconstructs trajectories for a wide variety of dynamic cellular processes. The performance was evaluated using a new, quantitative evaluation pipeline, comparing the performance of current state-of-the-art techniques on 10 publicly available single-cell RNA sequencing datasets. It automatically identifies marker genes, speeding up knowledge discovery.
- [SCOUP](https://github.com/hmatsu1226/SCOUP) - [C++] - Uses probabilistic model based on the Ornstein-Uhlenbeck process to analyze single-cell expression data during differentiation.
- [scran](http://bioconductor.org/packages/scran) - [R] - This package implements a variety of low-level analyses of single-cell RNA-seq data. Methods are provided for normalization of cell-specific biases, pool-based norms to estimate size factors, assignment of cell cycle phase, and detection of highly variable and significantly correlated genes.
- [SCRAT](https://github.com/zji90/SCRAT) - [R] - SCRAT provides essential tools for users to read in single-cell regolome data (ChIP-seq, ATAC-seq, DNase-seq) and summarize into different types of features. It also allows users to visualize the features, cluster samples and identify key features.
- [scTCRseq](https://github.com/ElementoLab/scTCRseq) - [python] - Map T-cell receptor (TCR) repertoires from single cell RNAseq.
- [SCUBA](https://github.com/gcyuan/SCUBA) - [matlab/R] - SCUBA stands for "Single-cell Clustering Using Bifurcation Analysis." SCUBA is a novel computational method for extracting lineage relationships from single-cell gene expression data, and modeling the dynamic changes associated with cell differentiation.
- [SEPA](https://github.com/zji90/SEPA) - [R] - SEPA provides convenient functions for users to assign genes into different gene expression patterns such as constant, monotone increasing and increasing then decreasing. SEPA then performs GO enrichment analysis to analysis the functional roles of genes with same or similar patterns.
- [Seurat](http://www.satijalab.org/seurat.html) - [R] - It contains easy-to-use implementations of commonly used analytical techniques, including the identification of highly variable genes, dimensionality reduction (PCA, ICA, t-SNE), standard unsupervised clustering algorithms (density clustering, hierarchical clustering, k-means), and the discovery of differentially expressed genes and markers.
- [SIMLR](https://github.com/BatzoglouLabSU/SIMLR) - [R, matlab] - SIMLR (Single-cell Interpretation via Multi-kernel LeaRning) learns an appropriate distance metric from the data for dimension reduction, clustering and visualization. SIMLR is capable of separating known subpopulations more accurately in single-cell data sets than do existing dimension reduction methods.
- [sincell](http://bioconductor.org/packages/sincell) - [R] - Existing computational approaches for the assessment of cell-state hierarchies from single-cell data might be formalized under a general workflow composed of i) a metric to assess cell-to-cell similarities (combined or not with a dimensionality reduction step), and ii) a graph-building algorithm (optionally making use of a cells-clustering step). Sincell R package implements a methodological toolbox allowing flexible workflows under such framework.
- [sincera](https://research.cchmc.org/pbge/sincera.html) - [R] - R-based pipeline for single-cell analysis including clustering and visualization.
- [SingleSplice](https://github.com/jw156605/SingleSplice) - [R, perl, C++] - A tool for detecting biological variation in alternative splicing within a population of single cells. See [Welch et al. 2016](https://academic.oup.com/nar/article/44/8/e73/2465993/Robust-detection-of-alternative-splicing-in-a).
- [singlet](https://github.com/iosonofabio/singlet) - [Python] - Single cell RNA-Seq analysis with phenotypes.
- [SinQC](http://www.morgridge.net/SinQC.html) - [R] - A Method and Tool to Control Single-cell RNA-seq Data Quality.
- [SLICER](https://github.com/jw156605/SLICER) - [R] - Selective Locally linear Inference of Cellular Expression Relationships (SLICER) algorithm for inferring cell trajectories.
- [slingshot](https://github.com/kstreet13/slingshot) - [R] - Functions for identifying and characterizing continuous developmental trajectories in single-cell sequencing data.
- [SPADE](http://www.nature.com/nprot/journal/v11/n7/full/nprot.2016.066.html) - [R] - Visualization and cellular hierarchy inference of single-cell data using SPADE.
- [splatter](http://bioconductor.org/packages/splatter/) - [R] - Splatter is a package for the simulation of single-cell RNA sequencing count data. It provides a simple interface for creating complex simulations that are reproducible and well-documented.
- [SPRING](https://github.com/AllonKleinLab/SPRING) - [matlab, javascript, python] - SPRING is a collection of pre-processing scripts and a web browser-based tool for visualizing and interacting with high dimensional data. SPRING was developed for single cell RNA-Seq data but can be applied more generally.
- [switchde](http://github.com/kieranrcampbell/switchde) - [R] - Differential expression analysis across pseudotime. Identify genes that exhibit switch-like up or down regulation along single-cell trajectories along with where in the trajectory the regulation occurs.
- [TASC](https://github.com/scrna-seq/TASC) - [C++, python] - To account for cell-to-cell technical differences, we propose a statistical framework, TASC (Toolkit for Analysis of Single Cell RNA-seq), an empirical Bayes approach to reliably model the cell-specific dropout rates and amplification bias by use of external RNA spike-ins. TASC incorporates the technical parameters, which reflect cell-to-cell batch effects, into a hierarchical mixture model to estimate the biological variance of a gene and detect differentially expressed genes. More importantly, TASC is able to adjust for covariates to further eliminate confounding that may originate from cell size and cell cycle differences.
- [TASIC](https://www.andrew.cmu.edu/user/sabrinar/TASIC) - [matlab] - TASIC is a new method for determining temporal trajectories, branching and cell assignments in single cell time series experiments. Unlike prior approaches TASIC uses on a probabilistic graphical model to integrate expression and time information making it more robust to noise and stochastic variations.
- [TopSLAM](https://github.com/mzwiessele/topslam) - [python] - Extracting and using probabilistic Waddington's landscape recreation from single cell gene expression measurements.
- [TraCeR](http://github.com/teichlab/tracer) - [python] - Reconstruction of T-Cell receptor sequences from single-cell RNA-seq data.
- [TRAPeS](https://github.com/yoseflab/trapes) - [python, C++] - TRAPeS (TCR Reconstruction Algorithm for Paired-End Single-cell), a software for reconstruction of T cell receptors (TCR) using short, paired-end single-cell RNA-sequencing.
- [TSCAN](https://github.com/zji90/TSCAN) - [R] - Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis.
- [ZIFA](https://github.com/epierson9/ZIFA) - [Python] - Zero-inflated dimensionality reduction algorithm for single-cell data.
- [zUMIs](https://github.com/sdparekh/zUMIs) - [R, perl, shell] - [zUMIs: A fast and flexible pipeline to process RNA-seq data with UMIs.](https://www.biorxiv.org/content/early/2017/10/18/153940)


### Copy number analysis

- [Ginkgo](https://github.com/robertaboukhalil/ginkgo) - [R, C] - Ginkgo is a web application for single-cell copy-number variation analysis.

### Variant calling

- [monovar](https://bitbucket.org/hamimzafar/monovar) - [python] - Monovar is a single nucleotide variant (SNV) detection and genotyping algorithm for single-cell DNA sequencing data. It takes a list of bam files as input and outputs a vcf file containing the detected SNVs.
- [SSrGE](https://github.com/lanagarmire/SSrGE) - [python] - SSrGE is an approach to identify SNVs correlated with Gene Expression using multiple regularized linear regressions. It contains its own pipeline to infer SNVs from scRNA-seq reads and is able to identify and sort genes and SNVs for a given cell subgroup. Deposited in [BioRvix in December 2016](http://biorxiv.org/content/early/2017/03/01/095810).

### Other applications

- [BASIC](http://ttic.uchicago.edu/~aakhan/BASIC/) - [python] - BASIC is a semi-de novo assembly method to determine the full-length sequence of the BCR in single B cells from scRNA-seq data.
- [dropSeqPipe](https://github.com/Hoohm/dropSeqPipe) - [python, R, snakemake] - An automatic data handling pipeline for drop-seq/scrb-seq data. It runs from raw fastq.gz data until the final count matrix with QC plots along the way.
- [MAGIC]( https://github.com/pkathail/magic) - [python] - [MAGIC: A diffusion-based imputation method reveals gene-gene interactions in single-cell RNA-sequencing data](http://www.biorxiv.org/content/early/2017/02/25/111591)
- [MATCHER]( https://github.com/jw156605/MATCHER) - [python] - [MATCHER: An algorithm for integrating single cell transcriptomic and epigenomic data using manifold alignment. MATCHER takes multiple types of single cell measurements performed on distinct single cells and infers single cell multi-omic profiles.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1269-0)
- [MetaNeighbor](github.com/maggiecrow/MetaNeighbor) - [R] - [MetaNeighbor: a method to rapidly assess cell type identity using both functional and random gene sets](https://doi.org/10.1101/150524)
- [powsimR](https://github.com/bvieth/powsimR) - [R] - [Power analysis for bulk and single cell RNA-seq experiments.](https://doi.org/10.1101/117150)
- [SiFit](https://bitbucket.org/hamimzafar/sifit) - [Java] - [SiFit: A Method for Inferring Tumor Trees from Single-Cell Sequencing Data under Finite-site Models](http://biorxiv.org/content/early/2016/12/04/091595)
- [sircel](https://github.com/pachterlab/Sircel) - [python] - sircel (pronounced "circle") separates reads in a fastq file based on barcode sequences that occur at known positions of reads. This is an essential first step in analyzing single-cell genomics data from experiments such as Drop-Seq. Barcode sequences often contain deletion and/or mismatch errors that arise during barcode synthesis and sequencing, and we have designed our barcode recovery approach with these issues in mind. In addition to identifying barcodes in an unbiased manner, sircel also quantifies their abundances. [doi](https://doi.org/10.1101/136242)
- [Wishbone](https://github.com/ManuSetty/wishbone) - [python] - [Wishbone is an algorithm to identify bifurcating developmental trajectories from single cell data. Wishbone can applied to both single cell RNA-seq and mass cytometry datasets.](http://www.nature.com/nbt/journal/v34/n6/full/nbt.3569.html)


## Tutorials and workflows

- [Aaron Lun's Single Cell workflow on Bioconductor](http://bioconductor.org/help/workflows/simpleSingleCell/) - [R] - This article describes a computational workflow for basic analysis of scRNA-seq data using software packages from the open-source Bioconductor project.
- [Bioconductor2016 Single-cell-RNA-sequencing workshop by Sandrine Dudoit lab](https://github.com/drisso/bioc2016singlecell) - [R] - SCONE, clusterExperiment, and slingshot tutorial.
- [BiomedCentral Single Cell Omics collectin](http://www.biomedcentral.com/collections/singlecellomics) - collection of papers describing techniques for single-cell analysis and protocols.
- [CSHL Single Cell Analysis - Bioinformatics](https://github.com/YeoLab/single-cell-bioinformatics/) course materials - Uses Shalek 2013 and Macaulay 2016 datasets to teach machine learning to biologists
- [Festival of Genomics California Single Cell Workshop](https://kdkorthauer.github.io/FestivalWorkshopVignettes/) - [R] - Explores basic workflow from exploratory data analysis to normalization and downstream analyses using a dataset of 1679 cells from the Allen Brain Atlas.
- [Gilad Lab Single Cell Data Exploration](http://jdblischak.github.io/singleCellSeq/analysis/) - R-based exploration of single cell sequence data. Lots of experimentation.
- [Harvard STEM Cell Institute Single Cell Workshop 2015](http://hms-dbmi.github.io/scw/) - workshop on common computational analysis techniques for scRNA-seq data from differential expression to subpopulation identification and network analysis. [See course description for more information](http://scholar.harvard.edu/jeanfan/classes/single-cell-workshop-2015)
- [Hemberg Lab scRNA-seq course materials](http://hemberg-lab.github.io/scRNA.seq.course/index.html)
- [Using Seurat (v1.2) for unsupervised clustering and biomarker discovery](http://www.satijalab.org/seurat/get_started_v1_2.html) - 301 single cells across diverse tissues from (Pollen et al., Nature Biotechnology, 2014). Original tutorial using Seurat 1.2
- [Using Seurat (v1.2) for spatial inference in single-cell data](http://www.satijalab.org/seurat/get_started_v1_2.html) - 851 single cells from Zebrafish embryogenesis (Satija*, Farrell* et al., Nature Biotechnology, 2015). Original tutorial using Seurat 1.2
- [Seurat (v2.0) - Guided Clustering Tutorial](http://satijalab.org/seurat/pbmc3k_tutorial.html) - new tutorial using Seurat 2.0

## Web portals and apps

- [10X Genomics datasets](https://support.10xgenomics.com/single-cell/datasets) - 10x genomics public datasets, including 1.3M cell mouse brain dataset.
- [ASAP](http://asap.epfl.ch/) - Automated Single-cell Analysis Pipeline (deposited in [BioRXiv](http://biorxiv.org/content/early/2016/12/22/096222) on December 22, 2016).
- [conquer](http://imlspenticton.uzh.ch:3838/conquer/) - A repository of consistently processed, analysis-ready single-cell RNA-seq data sets.
- [D<sup>3</sup>E](http://www.sanger.ac.uk/sanger/GeneRegulation_D3E/) - Discrete Distributional Differential Expression (D<sup>3</sup>E) is a tool for identifying differentially-expressed genes, based on single-cell RNA-seq data.
- [Ginkgo](http://qb.cshl.edu/ginkgo) - [R, C] - Ginkgo is a web application for single-cell copy-number variation analysis and visualization.
- [Granatum](http://garmiregroup.org/granatum/app) - Granatum 🍇 is a graphical single-cell RNA-seq (scRNA-seq) analysis pipeline for genomics scientists. [Deposited in Feb. 2017](http://biorxiv.org/content/early/2017/02/22/110759).
- [JingleBells](http://jinglebells.bgu.ac.il/) - A repository of standardized single cell RNA-Seq datasets for analysis and visualization in IGV at the single cell level. Currently focused on immune cells (http://www.jimmunol.org/content/198/9/3375.long).
- [Single Cell Portal](https://portals.broadinstitute.org/single_cell) - The Single-Cell Portal was developed to facilitate open data and open science in Single-cell Genomics. The portal currently focuses on sharing scientific results interactively, and sharing associated datasets.
- [SAKE](http://sake.mhammell.tools/) - Single-cell RNA-Seq Analysis and Clustering Evaluation.
- [scmap](http://www.hemberg-lab.cloud/scmap/) - A web tool for fast and accurate mapping of cells to a reference database using scRNA-seq data
- [scRNA.seq.datasets](https://hemberg-lab.github.io/scRNA.seq.datasets) - Collection of public scRNA-Seq datasets used by [Hemberg Lab](http://www.sanger.ac.uk/science/groups/hemberg-group)
- [scRNASeqDB](https://bioinfo.uth.edu/scrnaseqdb/) - A database aggregating human single-cell RNA-seq datasets. [ref](http://biorxiv.org/content/early/2017/01/31/104810)
- [SCPortalen](http://single-cell.clst.riken.jp/) - SCPortalen: human and mouse single-cell centric database. [ref](https://doi.org/10.1093/nar/gkx949)

## Journal articles of general interest

### Paper collections

- [Mendeley Single Cell Sequencing Analysis](https://www.mendeley.com/groups/9329461/single-cell-sequencing-analysis/papers/)
- [BioMedCentral Single-Cell -omics collection](http://www.biomedcentral.com/collections/singlecellomics)
- [Single-Cell Genomics in the Journal Science](http://science.sciencemag.org/content/358/6359) - Special issue on Single-Cell Genomics
 
### Big data approach overview
- [Single-cell Transcriptome Study as Big Data](http://www.sciencedirect.com/science/article/pii/S1672022916000437)

### Experimental design

- [Design and computational analysis of single-cell RNA-sequencing experiments](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y)

### Methods comparisons

- [Comparative analysis of single-cell RNA sequencing methods](https://doi.org/10.1016/j.molcel.2017.01.023) - a comparison of wet lab protocols for scRNA sequencing.
- [Comparison of methods to detect differentially expressed genes between single-cell populations](https://doi.org/10.1093/bib/bbw057) - comparison of five statistical methods to detect differentially expressed genes between two distinct single-cell populations.
- [Single-Cell RNA-Sequencing: Assessment of Differential Expression Analysis Methods](https://www.frontiersin.org/articles/10.3389/fgene.2017.00062/full) - an assessment of main bulk and single-cell differential analysis methods used to analyze scRNA-seq data.

## Similar lists and collections

- [CrazyHotTommy's RNA-seq analysis list](https://github.com/crazyhottommy/RNA-seq-analysis#single-cell-rna-seq) - Very broad list that includes some single cell RNA-seq packages and papers.
- [scRNA-tools.org](https://www.scrna-tools.org) - Database of scRNA-seq analysis tools and their functions. Managed through this [Github repository](https://github.com/Oshlack/scRNA-tools).
- [agitter's Pseudotime estimation list](https://github.com/agitter/single-cell-pseudotime) - An overview of algorithms for estimating pseudotime in single-cell RNA-seq data.

## People

Gender bias at conferences is a well known problem ([http://www.sciencemag.org/careers/2015/07/countering-gender-bias-conferences](http://www.sciencemag.org/careers/2015/07/countering-gender-bias-conferences)). Creating a list of potential speakers can help mitigate this bias and a community of people developing and maintaining helps to further diversify this list beyond smaller networks.

### Female

- [Christina Kendziorski (University of Wisconsin–Madison, USA)](https://www.biostat.wisc.edu/~kendzior/)
- [Sandrine Dudoit (UC Berkeley, USA)](http://www.stat.berkeley.edu/~sandrine/)
- [Keegan Korthauer (Dana Farber Cancer Institute, USA)](http://bcb.dfci.harvard.edu/~keegan/)
- [Stephanie Hicks (Dana Farber Cancer Institute, USA)](http://www.stephaniehicks.com/)
- [Dana Pe'er (Columbia University, USA)](http://www.c2b2.columbia.edu/danapeerlab/html/)
- [Alicia Oshlack (Murdock Children's Research Institute, Australia)](https://www.mcri.edu.au/users/dr-alicia-oshlack)
- [Aviv Regev (Broad Institute, USA)](https://www.broadinstitute.org/scientific-community/science/core-faculty-labs/regev-lab/regev-lab-home)
- [Catalina Vallejos (The Alan Turing Institute & UCL, UK)](https://sites.google.com/view/catalinavallejos)
- [Sarah Teichmann (Wellcome Trust Sanger Institute, UK)](http://www.teichlab.org/)
- [Emma Pierson (Stanford University, USA)](http://cs.stanford.edu/people/emmap1/)
- [Ning Leng (Morgridge Institute for Research, USA)](https://www.biostat.wisc.edu/~ningleng/)
- [Laleh Haghverdi (Institute of Computational Biology, Germany)](https://www.helmholtz-muenchen.de/icb/institute/staff/staff/ma/2453/-Haghverdi/index.html)
- [Rhonda Bacher (University of Wisconsin-Madison, USA)](https://twitter.com/rbacher)
- [Lana X. Garmire, (University oh Hawaii, Cancer Center, USA)](http://garmiregroup.org/)
- [Barbara Treutlein (Max Planck Institute for Evolutionary Anthropology, Germany)](http://www.treutleinlab.org/)
- [Samantha Morris (Depts of Dev. Bio. and Genetics, Washington University, St. Louis)](http://morrislab.wustl.edu/)

### Male

- [Raphael Gottardo (Fred Hutchinson Cancer Research Center, USA)](https://www.fredhutch.org/en/labs/profiles/gottardo-raphael.html)
- [Martin Hemberg (Sanger Institute, UK)](http://www.sanger.ac.uk/science/groups/hemberg-group)
- [Aaron Lun (Cancer Research UK, UK)](http://www.cruk.cam.ac.uk/users/aaron-lun)
- [John Marioni (EBI, UK)](http://www.ebi.ac.uk/research/marioni)
- [Davis McCarthy (EBI, UK)](https://sites.google.com/site/davismcc/)
- [John Reid (MRC Biostatistics Unit, Cambridge University, UK)](http://johnreid.github.io/)
- [Rahul Satija (New York Genome Center)](http://satijalab.org/)
- [Peter Sims (Columbia University, Department of Systems Biology)](http://www.columbia.edu/~pas2182/index.php/home-top.html)
- [Oliver Stegle (EBI, UK)](http://www.ebi.ac.uk/research/stegle)
- [Cole Trapnell (University of Washington, Department of Genome Sciences)](http://cole-trapnell-lab.github.io/)
- [Itai Yanai (New York University, School of Medicine, Institute for Computational Medicine, USA)](https://yanailab.org/about/)
