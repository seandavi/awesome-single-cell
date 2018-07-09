# awesome-single-cell

List of software packages (and the people developing these methods) for single-cell data analysis, including RNA-seq, ATAC-seq, etc. [Contributions welcome](https://github.com/seandavi/awesome-single-cell/blob/master/CONTRIBUTING.md)...

## Citation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1117762.svg)](https://doi.org/10.5281/zenodo.1117762)


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
- [CALISTA](https://github.com/CABSEL/CALISTA) - [R] - CALISTA provides a user-friendly toolbox for the analysis of single cell expression data. CALISTA accomplishes three major tasks: 1)	Identification of cell clusters in a cell population based on single-cell gene expression data, 2)	Reconstruction of lineage progression and produce transition genes, and 3)	Pseudotemporal ordering of cells along any given developmental paths in the lineage progression.
- [ccRemover](https://CRAN.R-project.org/package=ccRemover) - [R] - Removes the Cell-Cycle Effect from Single-Cell RNA-Sequencing Data. [Identifying and removing the cell-cycle effect from single-cell RNA-Sequencing data](https://www.nature.com/articles/srep33892).
- [CellCNN](https://github.com/eiriniar/CellCnn) - [Python] - Representation Learning for detection of phenotype-associated cell subsets
- [Cellity](https://github.com/teichlab/cellity) - [R] - Classification of low quality cells in scRNA-seq data using R
- [CellRanger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) - [Linux Binary] - Cell Ranger is a set of analysis pipelines that process Chromium single-cell RNA-seq output to align reads, generate gene-cell matrices and perform clustering and gene expression analysis. *Software requires registration with 10xgenomics.*
- [cellTree](https://www.bioconductor.org/packages/3.3/bioc/html/cellTree.html) - [R] - Cell population analysis and visualization from single cell RNA-seq data using a Latent Dirichlet Allocation model.
- [clusterExperiment](https://github.com/epurdom/clusterExperiment) - [R] - Functions for running and comparing many different clusterings of single-cell sequencing data. Meant to work with SCONE and slingshot.
- [Clustergrammer](https://github.com/maayanlab/clustergrammer) - [Python, JavaScript] - Interative web-based heatmap for visualizing and analyzing high dimensional biological data, including single-cell RNA-seq. Clustergrammer can be used within a Jupyter notebook as an interative widget that can be shared using GitHub and NBviewer, see [example notebook](http://nbviewer.jupyter.org/github/MaayanLab/CCLE_Clustergrammer/blob/master/notebooks/Clustergrammer_CCLE_Notebook.ipynb).
- [CytoGuide](https://cyteguide.cytosplore.org/) - [C++,D3] - [CyteGuide: Visual Guidance for Hierarchical Single-Cell Analysis](http://ieeexplore.ieee.org/document/8017575/)
- [DECENT](https://github.com/cz-ye/DECENT) - [R] - The unique features of scRNA-seq data have led to the development of novel methods for differential expression (DE) analysis. However, few of the existing DE methods for scRNA-seq data estimate the number of molecules pre-dropout and therefore do not explicitly distinguish technical and biological zeroes. We develop DECENT, a DE method for scRNA-seq data that adjusts for the imperfect capture efficiency by estimating the number of molecules pre-dropout.
- [DECODE](https://github.com/shmohammadi86/DECODE) - [R] - We develop an algorithm, called DECODE, to assess the extent of joint presence/absence of genes across different cells. We show that this network captures biologically-meaningful pathways, cell-type specific modules, and connectivity patterns characteristic of complex networks. We develop a model that uses this network to discriminate biological vs. technical zeros, by exploiting each gene's local neighborhood. For non-biological zeros, we build a predictive model to impute the missing value using their most informative neighbors.
- [DESCEND](https://github.com/jingshuw/descend) - [R] - DESCEND deconvolves the true gene expression distribution across cells for UMI scRNA-seq counts. It provides estimates of several distribution based statistics (five distribution measurements and the coefficients of covariates (such as batches or cell size)).
- [destiny](http://bioconductor.org/packages/destiny/) - [R] - Diffusion maps are spectral method for non-linear dimension reduction introduced by Coifman et al.(2005). Diffusion maps are based on a distance metric (diffusion distance) which is conceptually relevant to how differentiating cells follow noisy diffusion-like dynamics, moving from a pluripotent state towards more differentiated states.
- [DensityPath](https://doi.org/10.1101/276311) - [.] - DensityPath: a level-set algorithm to visualize and reconstruct cell developmental trajectories for large-scale single-cell RNAseq data
- [DeLorean](https://cran.r-project.org/web/packages/DeLorean/index.html) - [R] - Bayesian pseudotime estimation algorithm that uses Gaussian processes to model gene expression profiles and provides a full posterior for the pseudotimes.
- [dropClust](https://github.com/debsin/dropClust) - [R/Python] - Efficient clustering of ultra-large scRNA-seq data.
- [dropsim](https://github.com/marchinilab/dropsim) - [R] - Simulating droplet based scRNA-seq data.
- [dynverse](https://github.com/dynverse/dynverse/) - [R] - [A comparison of single-cell trajectory inference methods: towards more accurate and robust tools](https://doi.org/10.1101/276907)
- [ECLAIR](https://github.com/GGiecold/ECLAIR) - [python] - ECLAIR stands for Ensemble Clustering for Lineage Analysis, Inference and Robustness. Robust and scalable inference of cell lineages from gene expression data.
- [embeddr](https://github.com/kieranrcampbell/embeddr) - [R] - Embeddr creates a reduced dimensional representation of the gene space using a high-variance gene correlation graph and laplacian eigenmaps. It then fits a smooth pseudotime trajectory using principal curves.
- [Falco](https://github.com/VCCRI/Falco/) - [AWS cloud] - [Falco: A quick and flexible single-cell RNA-seq processing framework on the cloud](http://www.biorxiv.org/content/early/2016/07/15/064006.abstract).
- [FastProject](https://github.com/yoseflab/fastproject) - [Python] - Signature analysis on low-dimensional projections of single-cell expression data.
- [flotilla](https://github.com/yeolab/flotilla) - [Python] - Reproducible machine learning analysis of gene expression and alternative splicing data
- [GPfates](https://github.com/Teichlab/GPfates) - [Python] - Model transcriptional cell fates as mixtures of Gaussian Processes
- [GiniClust](https://github.com/lanjiangboston/GiniClust) - [Python/R] - GiniClust is a clustering method implemented in Python and R for detecting rare cell-types from large-scale single-cell gene expression data.  GiniClust can be applied to datasets originating from different platforms, such as multiplex qPCR data, traditional single-cell RNAseq or newly emerging UMI-based single-cell RNAseq, e.g. inDrops and Drop-seq.
- [HocusPocus](https://github.com/joeburns06/hocuspocus) - [R] - Basic PCA-based workflow for analysis and plotting of single cell RNA-seq data.
- [ICGS](https://github.com/nsalomonis/altanalyze) - [Python] - Iterative Clustering and Guide-gene Selection (Olsson et al. Nature 2016). Identify discrete, transitional and mixed-lineage states from diverse single-cell transcriptomics platforms. Integrated FASTQ pseudoalignment /quantification (Kallisto), differential expression, cell-type prediction and optional cell cycle exclusion analyses. Specialized methods for processing BAM and 10X Genomics spares matrix files. Associated single-cell splicing PSI methods (MultIPath-PSI). Apart of the AltAnalyze toolkit along with accompanying visualization methods (e.g., heatmap, t-SNE, SashimiPlots, network graphs). Easy-to-use graphical user and commandline interfaces.
- [inferCNV](https://github.com/broadinstitute/inferCNV) - [R] - Part of the TrinityCTAT (Trinity Cancer Transcriptome Analysis Toolkit).  Provides tools for copy-number inference from single-cell RNA-seq data.
- [iSEE](https://bioconductor.org/packages/iSEE/) - [R] - iSEE, interactive SummarizedExperiment Explorer. The iSEE package aims to provide an interactive user interface for exploring data in objects derived from the SummarizedExperiment class. Particular focus will be given to single-cell data in the SingleCellExperiment derived class. The interface is implemented with RStudio's Shiny, with a multi-panel setup for ease of navigation. Features include: dynamically linked charts, support for reproducibility by recording the exact code for every output, as well as guided tours to learn step-by-step the salient features of the user interface and of the data. A demo instance of the app is available at this address: http://shiny.imbei.uni-mainz.de:3838/iSEE.
- [knn-smoothing](https://github.com/yanailab/knn-smoothing) - [python or R or matlab] - The algorithm is based on the observation that across protocols, the technical noise exhibited by UMI-filtered scRNA-Seq data closely follows Poisson statistics. Smoothing is performed by first identifying the nearest neighbors of each cell in a step-wise fashion, based on variance-stabilized and partially smoothed expression profiles, and then aggregating their transcript counts.
- [MAGIC](https://github.com/pkathail/magic) - [python or matlab] - Markov Affinity-based Graph Imputation of Cells (MAGIC).
- [MAST](https://github.com/RGLab/MAST) - [R] - Model-based Analysis of Single-cell Transcriptomics (MAST) fits a two-part, generalized linear models that are specially adapted for bimodal and/or zero-inflated single cell gene expression data.
- [MERLoT](https://github.com/soedinglab/merlot) - [R/python] - Reconstructing complex lineage trees from scRNA-seq data using MERLoT.
- [mfa](https://github.com/kieranrcampbell/mfa) - [R] - [Probabilistic modeling of bifurcations in single-cell gene expression data using a Bayesian mixture of factor analyzers](https://wellcomeopenresearch.org/articles/2-19/v1)
- [K-Branches](https://github.com/theislab/kbranches) - [R] - The main idea behind the K-Branches method is to identify regions of interest (branching regions and tips) in differentiation trajectories of single cells. So far, K-Branches is intended to be used on the diffusion map representation of the data, so the user should either provide the data in diffusion map space or use the destiny package perform diffusion map dimensionality reduction.
- [M3Drop](https://github.com/tallulandrews/M3Drop) - [R] - Michaelis-Menten Modelling of Dropouts for scRNASeq.
- [MAST](https://github.com/RGLab/MAST) - [R] - Model-based Analysis of Single-cell Transcriptomics (MAST) fits a two-part, generalized linear models that are specially adapted for bimodal and/or zero-inflated single cell gene expression data
- [MIMOSCA](https://github.com/asncd/MIMOSCA) - [python] - A repository for the design and analysis of pooled single cell RNA-seq perturbation experiments (Perturb-seq).
- [Monocle](http://cole-trapnell-lab.github.io/monocle-release/) - [R] - Differential expression and time-series analysis for single-cell RNA-Seq.
- [netSmooth](https://github.com/BIMSBbioinfo/netSmooth) - [R] - netSmooth is a network-diffusion based method that uses priors for the covariance structure of gene expression profiles on scRNA-seq experiments in order to smooth expression values. We demonstrate that netSmooth improves clustering results of scRNA-seq experiments from distinct cell populations, time-course experiments, and cancer genomics.
- [NetworkInference](https://github.com/Tchanders/NetworkInference.jl) - [Julia] - Fast implementation of single-cell network inference algorithms: [Gene Regulatory Network Inference from Single-Cell Data Using Multivariate Information Measures]( http://www.cell.com/cell-systems/fulltext/S2405-4712(17)30386-1 )
- [nimfa](https://github.com/ccshao/nimfa) - [Python] - Nimfa is a Python scripting library which includes a number of published matrix factorization algorithms, initialization methods, quality and performance measures and facilitates the combination of these to produce new strategies. The library represents a unified and efficient interface to matrix factorization algorithms and methods.
- [OEFinder](https://github.com/lengning/OEFinder) - [R] - Identify ordering effect genes in single cell RNA-seq data. OEFinder shiny impelemention depends on packages shiny, shinyFiles, gdata, and EBSeq.
- [OncoNEM](https://bitbucket.org/edith_ross/onconem/src) - [R] -  OncoNEM is a probabilistic method for inferring intra-tumor evolutionarylineage trees from somatic single nucleotide variants of single cells. OncoNEM identifies homogeneous cellularsubpopulations and infers their genotypes as well as a tree describing their evolutionary relationships.
- [ouija](https://github.com/kieranrcampbell/ouija) - [R] - [A descriptive marker gene approach to single-cell pseudotime inference](https://doi.org/10.1101/060442)
- [ouijaflow](http://www.github.com/kieranrcampbell/ouijaflow) - [python] - [A descriptive marker gene approach to single-cell pseudotime inference](https://doi.org/10.1101/060442)
- [outrigger](https://github.com/YeoLab/outrigger) - [Python] - Outrigger is a program to calculate alternative splicing scores of RNA-Seq data based on junction reads and a *de novo*, custom annotation created with a graph database, especially made for single-cell analyses.
- [pcaReduce](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0984-y) - [R] - hierarchical clustering of single cell transcriptional profiles.
- [PHATE - Potential of Heat-diffusion for Affinity-based Transition Embedding](https://github.com/KrishnaswamyLab/PHATE) - [Python, R, Matlab] - PHATE is a tool for visualizing high dimensional single-cell data with natural progressions or trajectories. PHATE uses a novel conceptual framework for learning and visualizing the manifold inherent to biological systems in which smooth transitions mark the progressions of cells from one state to another.
- [PhenoPath](https://github.com/kieranrcampbell/phenopath) - [R] - Single-cell pseudotime with heterogeneous genetic and environmental backgrounds, including Bayesian significance testing of iteractions.
- [PoissonUMIs](https://github.com/tallulandrews/PoissonUMIs) - [R] - Poisson Modelling of scRNASeq UMI counts.
- [powsimR](https://bvieth.github.io/powsimR/) - [R] - Power analysis is essential to optimize the design of RNA-seq experiments and to assess and compare the power to detect differentially expressed genes. [PowsimR](https://academic.oup.com/bioinformatics/article/33/21/3486/3952669) is a flexible tool to simulate and evaluate differential expression from bulk and especially single-cell RNA-seq data making it suitable for a priori and posterior power analyses.
- [PyGMNormalize](https://github.com/ficusss/PyGMNormalize) - [Python] - Python implementation of [edgeR](https://www.ncbi.nlm.nih.gov/pubmed/19910308) normalization method for count matrices.
- [rMATS](http://rnaseq-mats.sourceforge.net/) - [Python] - RNA-Seq Multavariate Analysis of Transcript Splicing.
- [SAVER](https://github.com/mohuangx/SAVER) - [R] - SAVER (Single-cell Analysis Via Expression Recovery) implements a regularized regression prediction and empirical Bayes method to recover the true gene expression profile in noisy and sparse single-cell RNA-seq data.
- [SAKE](https://github.com/naikai/sake) - [R] - Single-cell RNA-Seq Analysis and Clustering Evaluation.
- [SC3](https://github.com/hemberg-lab/sc3) - [R] - SC3 is a tool for the unsupervised clustering of cells from single cell RNA-Seq experiments.
- [Scanpy](https://github.com/theislab/scanpy) - [Python] - Scanpy provides computationally efficient tools that scale up to very large data sets and enables simple integraton of advanced machine learning algorithms.
- [scater](https://bioconductor.org/packages/release/bioc/html/scater.html) - [R] - Scater places an emphasis on tools for quality control, visualisation and pre-processing of data before further downstream analysis, filling a useful niche between raw RNA-sequencing count or transcripts-per-million data and more focused downstream modelling tools such as monocle, scLVM, SCDE, edgeR, limma and so on.
- [scDD](https://github.com/kdkorthauer/scDD) - [R] - scDD (Single-Cell Differential Distributions) is a framework to identify genes with different expression patterns between biological groups of interest.  In addition to traditional differential expression, it can detect differences that are more complex and subtle than a mean shift.
- [SCDE](https://github.com/hms-dbmi/scde) - [R] - Differential expression using error models and overdispersion-based identification of important gene sets.
- [SCell](https://github.com/diazlab/SCell) - [matlab] - SCell is an integrated software tool for quality filtering, normalization, feature selection, iterative dimensionality reduction, clustering and the estimation of gene-expression gradients from large ensembles of single-cell RNA-seq datasets. SCell is open source, and implemented with an intuitive graphical interface.
- [SCIMITAR](https://github.com/dimenwarper/scimitar) - [Python] - Single Cell Inference of Morphing Trajectories and their Associated Regulation module (SCIMITAR) is a method for inferring biological properties from a pseudotemporal ordering. It can also be used to obtain progression-associated genes that vary along the trajectory, and genes that change their correlation structure over the trajectory; progression co-associated genes.
- [scImpute](https://github.com/Vivianstats/scImpute) - [R] - [scImpute: Accurate And Robust Imputation For Single Cell RNA-Seq Data](https://doi.org/10.1101/141598)
- [scLVM](https://github.com/PMBio/scLVM) - [R] - scLVM is a modelling framework for single-cell RNA-seq data that can be used to dissect the observed heterogeneity into different sources, thereby allowing for the correction of confounding sources of variation.  scLVM was primarily designed to account for cell-cycle induced variations in single-cell RNA-seq data where cell cycle is the primary soure of variability. 
- [scImpute](https://github.com/Vivianstats/scImpute) - [R] - [scImpute: Accurate And Robust Imputation For Single Cell RNA-Seq Data](doi:10.1038/s41467-018-03405-7)
- [SCENIC](https://gbiomed.kuleuven.be/english/research/50000622/lcb/tools/scenic) - [R] - [SCENIC: single-cell regulatory network inference and clustering](https://www.nature.com/articles/nmeth.4463)
- [scvis](https://bitbucket.org/jerry00/scvis-dev) - [python] - [Interpretable dimensionality reduction of single cell transcriptome data with deep generative models](https://doi.org/10.1101/178624)
- [scLVM](https://github.com/PMBio/scLVM) - [R] - scLVM is a modelling framework for single-cell RNA-seq data that can be used to dissect the observed heterogeneity into different sources, thereby allowing for the correction of confounding sources of variation.  scLVM was primarily designed to account for cell-cycle induced variations in single-cell RNA-seq data where cell cycle is the primary soure of variability.
- [scTDA](https://github.com/RabadanLab/scTDA) - [Python] - scTDA is an object oriented python library for topological data analysis of high-throughput single-cell RNA-seq data. It includes tools for the preprocessing, analysis, and exploration of single-cell RNA-seq data based on topological representations.
- [scmap](http://bioconductor.org/packages/scmap) - [R] - scmap is a method for projecting cells from a scRNA-seq experiment on to the cell-types identified in a different experiment.
- [SCMarker](https://github.com/KChen-lab/SCMarker) - [R] - SCMarker is a method performing ab initial marker gene set selection from scRNA-seq data to achieve improved clustering/cell-typing results. [SCMarker: ab initio marker selection for single cell transcriptome profiling](https://www.biorxiv.org/content/early/2018/07/04/356634).
- [SCnorm](https://github.com/rhondabacher/SCnorm) - [R] - A quantile regression based approach for robust normalization of single cell RNA-seq data.
- [SCODE](https://github.com/hmatsu1226/SCODE) - [R/Julia]- an efficient regulatory network inference algorithm from single-cell RNA-Seq during differentiation
- [SCONE](https://github.com/YosefLab/scone) - [R] - SCONE (Single-Cell Overview of Normalized Expression), a package for single-cell RNA-seq data quality control (QC) and normalization. This data-driven framework uses summaries of expression data to assess the efficacy of normalization workflows.
- [SCORPIUS](https://github.com/rcannood/SCORPIUS) - [R] - SCORPIUS an unsupervised approach for inferring developmental chronologies from single-cell RNA sequencing data. It accurately reconstructs trajectories for a wide variety of dynamic cellular processes. The performance was evaluated using a new, quantitative evaluation pipeline, comparing the performance of current state-of-the-art techniques on 10 publicly available single-cell RNA sequencing datasets. It automatically identifies marker genes, speeding up knowledge discovery.
- [SCOUP](https://github.com/hmatsu1226/SCOUP) - [C++] - Uses probabilistic model based on the Ornstein-Uhlenbeck process to analyze single-cell expression data during differentiation.
- [scran](http://bioconductor.org/packages/scran) - [R] - This package implements a variety of low-level analyses of single-cell RNA-seq data. Methods are provided for normalization of cell-specific biases, pool-based norms to estimate size factors, assignment of cell cycle phase, and detection of highly variable and significantly correlated genes.
- [SCRL](https://github.com/SuntreeLi/SCRL) - [C++] - [Network embedding-based representation learning for single cell RNA-seq data](https://doi.org/10.1093/nar/gkx750)
- [scruff](https://github.com/compbiomed/scruff) - [R] - An R package for preprocessing single cell RNA-seq (scRNA-seq) FASTQ reads generated by CEL-Seq and CEL-Seq2 protocols. It demultiplexes reads according to predetermined cell barcodes, aligns reads to reference genome using [Rsubread aligner](https://bioconductor.org/packages/release/bioc/html/Rsubread.html), and reports UMI (Unique Molecular Identifier) filtered count matrix ready for downstream analysis. It also provides functions to visualize the quality of data and the alignments of reads for individual cells. See [vignette](https://github.com/compbiomed/scruff/blob/master/vignettes/scruff.pdf).
- [scTCRseq](https://github.com/ElementoLab/scTCRseq) - [python] - Map T-cell receptor (TCR) repertoires from single cell RNAseq.
- [SCUBA](https://github.com/gcyuan/SCUBA) - [matlab/R] - SCUBA stands for "Single-cell Clustering Using Bifurcation Analysis." SCUBA is a novel computational method for extracting lineage relationships from single-cell gene expression data, and modeling the dynamic changes associated with cell differentiation.
- [SEPA](https://github.com/zji90/SEPA) - [R] - SEPA provides convenient functions for users to assign genes into different gene expression patterns such as constant, monotone increasing and increasing then decreasing. SEPA then performs GO enrichment analysis to analysis the functional roles of genes with same or similar patterns.
- [Seurat](http://www.satijalab.org/seurat.html) - [R] - It contains easy-to-use implementations of commonly used analytical techniques, including the identification of highly variable genes, dimensionality reduction (PCA, ICA, t-SNE), standard unsupervised clustering algorithms (density clustering, hierarchical clustering, k-means), and the discovery of differentially expressed genes and markers.
- [SIMLR](https://github.com/BatzoglouLabSU/SIMLR) - [R, matlab] - SIMLR (Single-cell Interpretation via Multi-kernel LeaRning) learns an appropriate distance metric from the data for dimension reduction, clustering and visualization. SIMLR is capable of separating known subpopulations more accurately in single-cell data sets than do existing dimension reduction methods.
- [sincell](http://bioconductor.org/packages/sincell) - [R] - Existing computational approaches for the assessment of cell-state hierarchies from single-cell data might be formalized under a general workflow composed of i) a metric to assess cell-to-cell similarities (combined or not with a dimensionality reduction step), and ii) a graph-building algorithm (optionally making use of a cells-clustering step). Sincell R package implements a methodological toolbox allowing flexible workflows under such framework.
- [sincera](https://research.cchmc.org/pbge/sincera.html) - [R] - R-based pipeline for single-cell analysis including clustering and visualization.
- [SINCERITIES](https://github.com/CABSEL/SINCERITIES) - [R/Matlab] - [Inferring gene regulatory networks from time-stamped single cell transcriptional expression profiles](https://academic.oup.com/bioinformatics/article/34/2/258/4158033)
- [SingleSplice](https://github.com/jw156605/SingleSplice) - [R, perl, C++] - A tool for detecting biological variation in alternative splicing within a population of single cells. See [Welch et al. 2016](https://academic.oup.com/nar/article/44/8/e73/2465993/Robust-detection-of-alternative-splicing-in-a).
- [singlet](https://github.com/iosonofabio/singlet) - [Python] - Single cell RNA-Seq analysis with phenotypes.
- [SinQC](http://www.morgridge.net/SinQC.html) - [R] - A Method and Tool to Control Single-cell RNA-seq Data Quality.
- [SLICER](https://github.com/jw156605/SLICER) - [R] - Selective Locally linear Inference of Cellular Expression Relationships (SLICER) algorithm for inferring cell trajectories.
- [slingshot](https://github.com/kstreet13/slingshot) - [R] - Functions for identifying and characterizing continuous developmental trajectories in single-cell sequencing data.
- [soupX](https://github.com/constantAmateur/SoupX) - [R] - An R package for the estimation and removal of cell free mRNA contamination in droplet based single cell RNA-seq data. The problem this package attempts to solve is that all droplet based single cell RNA-seq experiments also capture ambient mRNAs present in the input solution along with cell specific mRNAs of interest.
- [SPADE](http://www.nature.com/nprot/journal/v11/n7/full/nprot.2016.066.html) - [R] - Visualization and cellular hierarchy inference of single-cell data using SPADE.
- [SpatialIDE](https://media.nature.com/original/nature-assets/nmeth/journal/vaop/ncurrent/extref/nmeth.4636-S4.zip) - [Python, R] - [SpatialDE: identification of spatially variable genes](https://www.nature.com/articles/nmeth.4636)
- [splatter](http://bioconductor.org/packages/splatter/) - [R] - Splatter is a package for the simulation of single-cell RNA sequencing count data. It provides a simple interface for creating complex simulations that are reproducible and well-documented.
- [SPRING](https://github.com/AllonKleinLab/SPRING) - [matlab, javascript, python] - SPRING is a collection of pre-processing scripts and a web browser-based tool for visualizing and interacting with high dimensional data. SPRING was developed for single cell RNA-Seq data but can be applied more generally.
- [switchde](http://github.com/kieranrcampbell/switchde) - [R] - Differential expression analysis across pseudotime. Identify genes that exhibit switch-like up or down regulation along single-cell trajectories along with where in the trajectory the regulation occurs.
- [SWNE](https://github.com/yanwu2014/swne) - [R] - [Visualizing single-cell RNA-seq datasets with Similarity Weighted Nonnegative Embedding (SWNE)](https://www.biorxiv.org/content/early/2018/03/05/276261)
- [TASC](https://github.com/scrna-seq/TASC) - [C++, python] - To account for cell-to-cell technical differences, we propose a statistical framework, TASC (Toolkit for Analysis of Single Cell RNA-seq), an empirical Bayes approach to reliably model the cell-specific dropout rates and amplification bias by use of external RNA spike-ins. TASC incorporates the technical parameters, which reflect cell-to-cell batch effects, into a hierarchical mixture model to estimate the biological variance of a gene and detect differentially expressed genes. More importantly, TASC is able to adjust for covariates to further eliminate confounding that may originate from cell size and cell cycle differences.
- [TASIC](https://www.andrew.cmu.edu/user/sabrinar/TASIC) - [matlab] - TASIC is a new method for determining temporal trajectories, branching and cell assignments in single cell time series experiments. Unlike prior approaches TASIC uses on a probabilistic graphical model to integrate expression and time information making it more robust to noise and stochastic variations.
- [TopSLAM](https://github.com/mzwiessele/topslam) - [python] - Extracting and using probabilistic Waddington's landscape recreation from single cell gene expression measurements.
- [TraCeR](http://github.com/teichlab/tracer) - [python] - Reconstruction of T-Cell receptor sequences from single-cell RNA-seq data.
- [TRAPeS](https://github.com/yoseflab/trapes) - [python, C++] - TRAPeS (TCR Reconstruction Algorithm for Paired-End Single-cell), a software for reconstruction of T cell receptors (TCR) using short, paired-end single-cell RNA-sequencing.
- [trendsceek](https://github.com/edsgard/trendsceek) - [R] - [Identification of spatial expression trends in single-cell gene expression data](https://www.nature.com/articles/nmeth.4634)
- [TSCAN](https://github.com/zji90/TSCAN) - [R] - Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis.
- [UNCURL](https://github.com/yjzhang/uncurl_python) - [Python] - Unsupervised and semi-supervised sampling effect removal for single-cell RNA-seq data.
- [ZIFA](https://github.com/epierson9/ZIFA) - [Python] - Zero-inflated dimensionality reduction algorithm for single-cell data.
- [zinbwaveZinger](https://github.com/statOmics/zinbwaveZinger) - [R] - We introduce a weighting strategy, based on a zero-inflated negative binomial model, that identifies excess zero counts and generates gene- and cell-specific weights to unlock bulk RNA-seq DE pipelines for zero-inflated data, boosting performance for scRNA-seq. https://doi.org/10.1186/s13059-018-1406-4
- [zUMIs](https://github.com/sdparekh/zUMIs) - [R, perl, shell] - [zUMIs: A fast and flexible pipeline to process RNA-seq data with UMIs.](https://www.biorxiv.org/content/early/2017/10/18/153940)


### Doublet Identification

- [demuxlet](https://github.com/statgen/demuxlet) - [shell] - [Multiplexed droplet single-cell RNA-sequencing using natural genetic variation](https://www.nature.com/articles/nbt.4042)
- DoubletFinder - [R] - Doublet detection in single-cell RNA sequencing data using artificial nearest neighbors[BioRxiv](https://www.biorxiv.org/content/early/2018/06/20/352484)
- [DoubletDetection](https://github.com/JonathanShor/DoubletDetection) - [R, Python] - A Python3 package to detect doublets (technical errors) in single-cell RNA-seq count matrices. An [R implementation](https://github.com/TomKellyGenetics/DoubletDetection) is in development.


### Copy number analysis

- [aneufinder](https://github.com/ataudt/aneufinder) - [R] - Bioconductor module for copy-number detection in single-cell whole genome sequencing (scWGS) and strand-seq data using a Hidden Markov Model or binary bisection method.
- [Ginkgo](https://github.com/robertaboukhalil/ginkgo) - [R, C] - Ginkgo is a web application for single-cell copy-number variation analysis.
- [HoneyBADGER](https://github.com/JEFworks/HoneyBADGER) - [R] - HoneyBADGER identifies and infers the presence of CNV and LOH events in single cells and reconstructs subclonal architecture using allele and expression information from single-cell RNA-sequencing data.

### Variant calling

- [monovar](https://bitbucket.org/hamimzafar/monovar) - [python] - Monovar is a single nucleotide variant (SNV) detection and genotyping algorithm for single-cell DNA sequencing data. It takes a list of bam files as input and outputs a vcf file containing the detected SNVs.
- [SCIPhi](https://github.com/cbg-ethz/SCIPhI) - [python] - [Single-cell mutation identification via phylogenetic inference (SCIPhI) is a new approach to mutation detection in individual tumor cells by leveraging the evolutionary relationship among cells.](https://www.biorxiv.org/content/early/2018/03/28/290908)
- [SSrGE](https://github.com/lanagarmire/SSrGE) - [python] - SSrGE is an approach to identify SNVs correlated with Gene Expression using multiple regularized linear regressions. It contains its own pipeline to infer SNVs from scRNA-seq reads and is able to identify and sort genes and SNVs for a given cell subgroup. Deposited in [BioRxiv in December 2016](http://biorxiv.org/content/early/2017/03/01/095810).

### Epigenomics

- [ChromVAR](https://bioconductor.org/packages/release/bioc/html/chromVAR.html) - [R] - Determine variations in chromatin accessibility across sets of annotations or peaks. Designed primarily for single-cell or sparse chromatin accessibility data, e.g. from scATAC-seq or sparse bulk ATAC or DNAse-seq experiments. [BioRxiv](https://www.biorxiv.org/content/early/2017/02/21/110346)
- [DeepCpg](https://github.com/cangermueller/deepcpg) - [python] - DeepCpG is a deep neural network for predicting the methylation state of CpG dinucleotides in multiple cells. It allows to accurately impute incomplete DNA methylation profiles, to discover predictive sequence motifs, and to quantify the effect of sequence mutations.
- [Melissa](https://github.com/andreaskapou/Melissa) - [R] - Melissa (MEthyLation Inference for Single cell Analysis), a Bayesian hierarchical method to quantify spatially-varying methylation profiles across genomic regions from single-cell bisulfite sequencing data (scBS-seq). Melissa clusters individual cells based on local methylation patterns, enabling the discovery of epigenetic differences and similarities among individual cells. The clustering also acts as an effective regularisation method for imputation of methylation on unassayed CpG sites, enabling transfer of information between individual cells. [BioRxiv](https://doi.org/10.1101/312025)
- [SCRAT](https://github.com/zji90/SCRAT) - [R] - SCRAT provides essential tools for users to read in single-cell regolome data (ChIP-seq, ATAC-seq, DNase-seq) and summarize into different types of features. It also allows users to visualize the features, cluster samples and identify key features.

### Multi-assay data integration

- [MATCHER]( https://github.com/jw156605/MATCHER) - [python] - [MATCHER: An algorithm for integrating single cell transcriptomic and epigenomic data using manifold alignment. MATCHER takes multiple types of single cell measurements performed on distinct single cells and infers single cell multi-omic profiles.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1269-0)
- [MOFA]( https://github.com/bioFAM/MOFA/tree/master/mofa/core) - [python, R] - Multi‐Omics Factor Analysis, a framework for unsupervised integration of multi‐omics data sets. MOFA is a method for disentangling the different sources of heterogeneity in bulk and single-cell multi-omics data sets. It identifies the latent factors that drive unique and shared variability in the different assays. The factors can be used for visualisation, pseudotime reconstruction, imputation, among other functionalities. [Paper](http://msb.embopress.org/content/14/6/e8124)
 
### Other applications

- [BASIC](http://ttic.uchicago.edu/~aakhan/BASIC/) - [python] - BASIC is a semi-de novo assembly method to determine the full-length sequence of the BCR in single B cells from scRNA-seq data.
- [dropEst](https://github.com/hms-dbmi/dropEst) - [C++, R] - High-performance pipeline for initial analysis of droplet-based single-cell RNA-seq data (Drop-seq, inDrop, 10x and some others). Allows to estimate gene count matrix as well as diagnostic stats from fastq files with raw reads. Implements corrections for different noise sources.
- [dropSeqPipe](https://github.com/Hoohm/dropSeqPipe) - [python, R, snakemake] - An automatic data handling pipeline for drop-seq/scrb-seq data. It runs from raw fastq.gz data until the final count matrix with QC plots along the way.
- [MAGIC]( https://github.com/pkathail/magic) - [python] - [MAGIC: A diffusion-based imputation method reveals gene-gene interactions in single-cell RNA-sequencing data](http://www.biorxiv.org/content/early/2017/02/25/111591)
- [MetaNeighbor](github.com/maggiecrow/MetaNeighbor) - [R] - [MetaNeighbor: a method to rapidly assess cell type identity using both functional and random gene sets](https://doi.org/10.1101/150524)
- [sasc](https://github.com/sciccolella/sasc) - [C] - sasc stands for Simulated Annealing Single-Cell, an algorithm for performing phylogenetic analysis of single-cell cancer samples.  Manuscript [here](https://www.biorxiv.org/content/early/2018/02/20/268243).
- [SCope](https://github.com/aertslab/SCope) - [python] - SCope is a fast visualization tool for large-scale and high dimensional scRNA-seq datasets. Publication [here](https://doi.org/10.1016/j.cell.2018.05.057).
- [SiFit](https://bitbucket.org/hamimzafar/sifit) - [Java] - [SiFit: A Method for Inferring Tumor Trees from Single-Cell Sequencing Data under Finite-site Models](http://biorxiv.org/content/early/2016/12/04/091595)
- [sircel](https://github.com/pachterlab/Sircel) - [python] - sircel (pronounced "circle") separates reads in a fastq file based on barcode sequences that occur at known positions of reads. This is an essential first step in analyzing single-cell genomics data from experiments such as Drop-Seq. Barcode sequences often contain deletion and/or mismatch errors that arise during barcode synthesis and sequencing, and we have designed our barcode recovery approach with these issues in mind. In addition to identifying barcodes in an unbiased manner, sircel also quantifies their abundances. [doi](https://doi.org/10.1101/136242)
- [Wishbone](https://github.com/ManuSetty/wishbone) - [python] - [Wishbone is an algorithm to identify bifurcating developmental trajectories from single cell data. Wishbone can applied to both single cell RNA-seq and mass cytometry datasets.](http://www.nature.com/nbt/journal/v34/n6/full/nbt.3569.html)
- [Snakemake single-cell-rna-seq workflow](https://github.com/snakemake-workflows/single-cell-rna-seq) - [python, R, snakemake] - An automated pipeline for single cell RNA-seq analysis.


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
- [CellView](https://github.com/mohanbolisetty/CellView) - CellView is an R Shiny web application that allows knowledge-based and hypothesis-driven exploration of processed single cell transcriptomic data. [ref](https://www.biorxiv.org/content/early/2017/04/04/123810).
- [conquer](http://imlspenticton.uzh.ch:3838/conquer/) - A repository of consistently processed, analysis-ready single-cell RNA-seq data sets.
- [D<sup>3</sup>E](http://www.sanger.ac.uk/sanger/GeneRegulation_D3E/) - Discrete Distributional Differential Expression (D<sup>3</sup>E) is a tool for identifying differentially-expressed genes, based on single-cell RNA-seq data.
- [Ginkgo](http://qb.cshl.edu/ginkgo) - [R, C] - Ginkgo is a web application for single-cell copy-number variation analysis and visualization.
- [Granatum](http://garmiregroup.org/granatum/app) - Granatum 🍇 is a graphical single-cell RNA-seq (scRNA-seq) analysis pipeline for genomics scientists. [Published in December 2017](https://doi.org/10.1186/s13073-017-0492-3).
- [JingleBells](http://jinglebells.bgu.ac.il/) - A repository of standardized single cell RNA-Seq datasets for analysis and visualization in IGV at the single cell level. Currently focused on immune cells (http://www.jimmunol.org/content/198/9/3375.long).
- [SAKE](http://sake.mhammell.tools/) - Single-cell RNA-Seq Analysis and Clustering Evaluation.
- [scmap](http://www.hemberg-lab.cloud/scmap/) - A web tool for fast and accurate mapping of cells to a reference database using scRNA-seq data
- [SCPortalen](http://single-cell.clst.riken.jp/) - SCPortalen: human and mouse single-cell centric database. [ref](https://doi.org/10.1093/nar/gkx949)
- [scRNA.seq.datasets](https://hemberg-lab.github.io/scRNA.seq.datasets) - Collection of public scRNA-Seq datasets used by [Hemberg Lab](http://www.sanger.ac.uk/science/groups/hemberg-group)
- [scRNASeqDB](https://bioinfo.uth.edu/scrnaseqdb/) - A database aggregating human single-cell RNA-seq datasets. [ref](http://biorxiv.org/content/early/2017/01/31/104810)
- [ShinyCortex](https://bioinf.eva.mpg.de/shiny/sample-apps/ShinyCortex/) - a resource that brings together data from recent scRNA-seq studies of the developing cortex for further analysis. ShinyCortex is based in R and displays recently published scRNA-seq data from the human and mouse cortex in a comprehensible, dynamic and accessible way, suitable for data exploration by biologists. [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5962798/)
- [Single Cell Portal](https://portals.broadinstitute.org/single_cell) - The Single-Cell Portal was developed to facilitate open data and open science in Single-cell Genomics. The portal currently focuses on sharing scientific results interactively, and sharing associated datasets.
- [singleCellTK](https://bioconductor.org/packages/release/bioc/html/singleCellTK.html) - The singleCellTK is an R/Shiny package and GUI for analyzing and visualizing scRNA-Seq through a web interface. Analysis modules include data summary and filtering, dimensionality reduction and clustering, batch correction, differential expression analysis, pathway activity analysis, and power analysis.

## Journal articles of general interest

### Paper collections

- [Mendeley Single Cell Sequencing Analysis](https://www.mendeley.com/groups/9329461/single-cell-sequencing-analysis/papers/)
- [BioMedCentral Single-Cell -omics collection](http://www.biomedcentral.com/collections/singlecellomics)
- [Single-Cell Genomics in the Journal Science](http://science.sciencemag.org/content/358/6359) - Special issue on Single-Cell Genomics
- [The emerging field of single-cell analysis](https://www.sciencedirect.com/journal/molecular-aspects-of-medicine/vol/59/suppl/C) - Special issue on  single cell analysis

### Big data approach overview
- [Single-cell Transcriptome Study as Big Data](http://www.sciencedirect.com/science/article/pii/S1672022916000437)

### Experimental design

- [Design and computational analysis of single-cell RNA-sequencing experiments](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y)
- [How to design a single-cell RNA-sequencing experiment: pitfalls, challenges and perspectives](https://doi.org/10.1093/bib/bby007)

### Methods comparisons

- [Comparative analysis of single-cell RNA sequencing methods](https://doi.org/10.1016/j.molcel.2017.01.023) - a comparison of wet lab protocols for scRNA sequencing.
- [Comparison of computational methods for imputing single-cell RNA-sequencing data](https://doi.org/10.1101/241190) - We compared eight imputation methods, evaluated their power in recovering original real data, and performed broad analyses to explore their effects on clustering cell types, detecting differentially expressed genes, and reconstructing lineage trajectories in the context of both simulated and real data. Simulated datasets and case studies highlight that there are no one method performs the best in all the situations.
- [Comparison of methods to detect differentially expressed genes between single-cell populations](https://doi.org/10.1093/bib/bbw057) - comparison of five statistical methods to detect differentially expressed genes between two distinct single-cell populations.
- [Bias, Robustness And Scalability In Differential Expression Analysis Of Single-Cell RNA-Seq Data](http://dx.doi.org/10.1038/nmeth.4612) - comparison of 36 statistical methods to detect differentially expressed genes between two annotated populations from the [conquer](http://imlspenticton.uzh.ch:3838/conquer/) database of consistently processed scRNA-seq datasets.
- [Single-Cell RNA-Sequencing: Assessment of Differential Expression Analysis Methods](https://www.frontiersin.org/articles/10.3389/fgene.2017.00062/full) - an assessment of main bulk and single-cell differential analysis methods used to analyze scRNA-seq data.
- [A comparison of single-cell trajectory inference methods: towards more accurate and robust tools](https://doi.org/10.1101/276907) - A comparison of 29 trajectory inference methods on real and synthetic data.  

## Similar lists and collections

- [CrazyHotTommy's RNA-seq analysis list](https://github.com/crazyhottommy/RNA-seq-analysis#single-cell-rna-seq) - Very broad list that includes some single cell RNA-seq packages and papers.
- [scRNA-tools.org](https://www.scrna-tools.org) - Database of scRNA-seq analysis tools and their functions. Managed through this [Github repository](https://github.com/Oshlack/scRNA-tools).
- [agitter's Pseudotime estimation list](https://github.com/agitter/single-cell-pseudotime) - An overview of algorithms for estimating pseudotime in single-cell RNA-seq data.

## People

Gender bias at conferences is a well known problem ([http://www.sciencemag.org/careers/2015/07/countering-gender-bias-conferences](http://www.sciencemag.org/careers/2015/07/countering-gender-bias-conferences)). Creating a list of potential speakers can help mitigate this bias and a community of people developing and maintaining helps to further diversify this list beyond smaller networks.

### Female

- [Rhonda Bacher (University of Wisconsin-Madison, USA)](https://twitter.com/rbacher)
- [Barbara Di Camillo (Information Engineering Department, University of Padova, Italy](https://www.dei.unipd.it/~dicamill/)
- [Jinmiao Chen (Singapore Immunology Network, A\*STAR, Singapore)](https://www.a-star.edu.sg/sign/PEOPLE/Principal-Investigators/Investigator-Details?givenName=Jinmiao&lastName=CHEN)
- [Sandrine Dudoit (UC Berkeley, USA)](http://www.stat.berkeley.edu/~sandrine/)
- [Lana X. Garmire, (University of Hawaii Cancer Center, USA)](http://garmiregroup.org/)
- [Laleh Haghverdi (EMBL, Germany)](https://www.embl.de/research/units/genome_biology/huber/members/index.php?s_personId=CP-60025160)
- [Stephanie Hicks (Dana Farber Cancer Institute, USA)](http://www.stephaniehicks.com/)
- [Christina Kendziorski (University of Wisconsin–Madison, USA)](https://www.biostat.wisc.edu/~kendzior/)
- [Keegan Korthauer (Dana Farber Cancer Institute, USA)](http://bcb.dfci.harvard.edu/~keegan/)
- [Ning Leng (Morgridge Institute for Research, USA)](https://www.biostat.wisc.edu/~ningleng/)
- [Elisabetta Mereu (Centre for Genomic Regulation, Barcelona)](http://www.cnag.cat/elisabetta-mereu)
- [Samantha Morris (Depts of Dev. Bio. and Genetics, Washington University, St. Louis)](http://morrislab.wustl.edu/)
- [Alicia Oshlack (Murdoch Children's Research Institute, Australia)](https://www.mcri.edu.au/users/dr-alicia-oshlack)
- [Dana Pe'er (Columbia University, USA)](http://www.c2b2.columbia.edu/danapeerlab/html/)
- [Emma Pierson (Stanford University, USA)](http://cs.stanford.edu/people/emmap1/)
- [Aviv Regev (Broad Institute, USA)](https://www.broadinstitute.org/scientific-community/science/core-faculty-labs/regev-lab/regev-lab-home)
- [Charlotte Soneson (Institute of Molecular Life Sciences, University of Zurich)](https://csoneson.github.io/)
- [Sarah Teichmann (Wellcome Trust Sanger Institute, UK)](http://www.teichlab.org/)
- [Barbara Treutlein (Max Planck Institute for Evolutionary Anthropology, Germany)](http://www.treutleinlab.org/)
- [Catalina Vallejos (The Alan Turing Institute & UCL, UK)](https://sites.google.com/view/catalinavallejos)

### Male

- [Stein Aerts (KU Leuven Center for Human Genetics, Belgium)](http://aertslab.org/)
- [Bart DePlancke (EPFL, School of Life sciences, Institute of Bioengineering, Switzerland)](https://deplanckelab.epfl.ch/)
- [Raphael Gottardo (Fred Hutchinson Cancer Research Center, USA)](https://www.fredhutch.org/en/labs/profiles/gottardo-raphael.html)
- [Chung Chau Hon (RIKEN Centre for Integrative Medical Sciences, Yokohama)](http://www.riken.jp/en/research/labs/ims/genom_inf/)
- [Martin Hemberg (Sanger Institute, UK)](http://www.sanger.ac.uk/science/groups/hemberg-group)
- [Holger Heyn (Centre for Genomic Regulation, Barcelona)](http://www.cnag.crg.eu/holger-heyn)
- [Peter Kharchenko (Department of Biomedical Informatics, Harvard Medical School, USA)](http://pklab.med.harvard.edu/)
- [Sten Linnarson (Karolinska Institutet, Sweden)](http://linnarssonlab.org/)
- [Aaron Lun (Cancer Research UK, UK)](http://www.cruk.cam.ac.uk/users/aaron-lun)
- [John Marioni (EBI, UK)](http://www.ebi.ac.uk/research/marioni)
- [Davis McCarthy (EBI, UK)](https://sites.google.com/site/davismcc/)
- [John Reid (MRC Biostatistics Unit, Cambridge University, UK)](http://johnreid.github.io/)
- [Mark Robinson (Institute of Molecular Life Sciences, University of Zurich)](http://www.imls.uzh.ch/en/research/robinson.html)
- [Yvan Saeys (Vlaams Instituut voor Biotechnologie, Ghent, Belgium)](http://www.vib.be/en/research/scientists/Pages/Yvan-Saeys-Lab.aspx)
- [Rahul Satija (New York Genome Center)](http://satijalab.org/)
- [Peter Sims (Columbia University, Department of Systems Biology)](http://www.columbia.edu/~pas2182/index.php/home-top.html)
- [Oliver Stegle (EBI, UK)](http://www.ebi.ac.uk/research/stegle)
- [Fabian Theis (Institute of Computational Biology, Helmholtz Zentrum München)](https://www.helmholtz-muenchen.de/icb/institute/staff/staff/ma/2494/index.html)
- [Cole Trapnell (University of Washington, Department of Genome Sciences)](http://cole-trapnell-lab.github.io/)
- [Itai Yanai (New York University, School of Medicine, Institute for Computational Medicine, USA)](https://yanailab.org/about/)

