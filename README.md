# awesome-single-cell

List of software packages for single-cell data analysis, including RNA-seq, ATAC-seq, etc. [Contributions welcome](https://github.com/seandavi/awesome-single-cell/blob/master/CONTRIBUTE.md)...

## Software packages

- [BPSC](https://github.com/nghiavtr/BPSC) - [R] - Beta-Poisson model for single-cell RNA-seq data analyses
- [clusterExperiment](https://github.com/epurdom/clusterExperiment) - [R] - Functions for running and comparing many different clusterings of single-cell sequencing data. Meant to work with SCONE and slingshot.
- [HocusPocus](https://github.com/joeburns06/hocuspocus) - [R] - Basic PCA-based workflow for analysis and plotting of single cell RNA-seq data.
- [MAST](https://github.com/RGLab/MAST) - [R] - Model-based Analysis of Single-cell Transcriptomics
(MAST) fits a two-part, generalized linear models that are specially adapted for bimodal and/or zero-inflated single cell gene expression data.
- [Monocle](http://cole-trapnell-lab.github.io/monocle-release/) - [R] - Differential expression and time-series analysis for single-cell RNA-Seq.
- [OEFinder](https://github.com/lengning/OEFinder) - [R] - Identify ordering effect genes in single cell RNA-seq data. OEFinder shiny impelemention depends on packages shiny, shinyFiles, gdata, and EBSeq.
- [SCDE](https://github.com/hms-dbmi/scde) - [R] - Differential expression using error models and overdispersion-based identification of important gene sets.
- [SCell](https://github.com/diazlab/SCell) - [matlab] - SCell is an integrated software tool for quality filtering, normalization, feature selection, iterative dimensionality reduction, clustering and the estimation of gene-expression gradients from large ensembles of single-cell RNA-seq datasets. SCell is open source, and implemented with an intuitive graphical interface.
- [scLVM](https://github.com/PMBio/scLVM) - [R] - scLVM is a modelling framework for single-cell RNA-seq data that can be used to dissect the observed heterogeneity into different sources, thereby allowing for the correction of confounding sources of variation.  scLVM was primarily designed to account for cell-cycle induced variations in single-cell RNA-seq data where cell cycle is the primary soure of variability. 
- [SCONE](https://github.com/YosefLab/scone) - [R] - SCONE (Single-Cell Overview of Normalized Expression), a package for single-cell RNA-seq data quality control (QC) and normalization. This data-driven framework uses summaries of expression data to assess the efficacy of normalization workflows.
- [SCOUP](https://github.com/hmatsu1226/SCOUP) - [C++] - Uses probabilistic model based on the Ornstein-Uhlenbeck process to analyze single-cell expression data during differentiation. 
- [sincell](http://bioconductor.org/packages/sincell) - [R] - Existing computational approaches for the assessment of cell-state hierarchies from single-cell data might be formalized under a general workflow composed of i) a metric to assess cell-to-cell similarities (combined or not with a dimensionality reduction step), and ii) a graph-building algorithm (optionally making use of a cells-clustering step). Sincell R package implements a methodological toolbox allowing flexible workflows under such framework.
- [sincera](https://research.cchmc.org/pbge/sincera.html) - [R] - R-based pipeline for single-cell analysis including clustering and visualization.
- [SinQC](http://www.morgridge.net/SinQC.html) - [R] - A Method and Tool to Control Single-cell RNA-seq Data Quality.
- [SLICER](https://github.com/jw156605/SLICER) - [R] - Selective Locally linear Inference of Cellular Expression Relationships (SLICER) algorithm for inferring cell trajectories.
- [slingshot](https://github.com/kstreet13/slingshot) - [R] - Functions for identifying and characterizing continuous developmental trajectories in single-cell sequencing data.
- [SPADE](http://www.nature.com/nprot/journal/v11/n7/full/nprot.2016.066.html) - [R] - Visualization and cellular hierarchy inference of single-cell data using SPADE.
- [TSCAN](https://github.com/zji90/TSCAN) - [R] - Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis.

## Tutorials and workflows
 
- [Aaron Lun's Single Cell workflow on Bioconductor](http://bioconductor.org/help/workflows/simpleSingleCell/) - [R] - This article describes a computational workflow for basic analysis of scRNA-seq data using software packages from the open-source Bioconductor project.
- [Bioconductor2016 Single-cell-RNA-sequencing workshop by Sandrine Dudoit lab](https://github.com/drisso/bioc2016singlecell) - [R] - SCONE, clusterExperiment, and slingshot tutorial.
- [BiomedCentral Single Cell Omics collectin](http://www.biomedcentral.com/collections/singlecellomics) - collection of papers describing techniques for single-cell analysis and protocols.
- [Hemberg Lab scRNA-seq course materials](http://hemberg-lab.github.io/scRNA.seq.course/index.html) 


Comparative analysis of methods: http://biorxiv.org/content/early/2016/01/05/035758

Review of experimental design and analysis: http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y

QC: http://www.morgridge.net/SinQC.html

# Normalization

http://michelebusby.tumblr.com/post/130202229486/the-ks-test-looks-pretty-good-for-single-cell

Accounting for technical variation: http://www.nature.com/ncomms/2015/151022/ncomms9687/full/ncomms9687.html#supplementary-information

ZIFA: Zero-inflated factor analysis https://github.com/epierson9/ZIFA

http://biorxiv.org/content/early/2016/04/22/049734.full.pdf+html

https://www.dropbox.com/s/pno78mmlj0exv7s/NODES_0.0.0.9010.tar.gz?dl=0

scran: http://bioconductor.org/packages/devel/bioc/html/scran.html

Correct for expression heterogeneity: https://github.com/PMBio/scLVM

# Transcript counting

Modified version of Kallisto: https://github.com/govinda-kamath/clustering_on_transcript_compatibility_counts

DISCO: https://pbtech-vc.med.cornell.edu/git/mason-lab/disco/tree/master

# Clustering

Comparative analysis: http://biorxiv.org/content/early/2016/04/07/047613

SC3: consensus clustering https://github.com/hemberg-lab/sc3

destiny: diffusion maps for single-cell data http://bioconductor.org/packages/release/bioc/html/destiny.html

https://github.com/govinda-kamath/clustering_on_transcript_compatibility_counts

GiniClust https://github.com/lanjiangboston/GiniClust

pcaReduce: https://github.com/JustinaZ/pcaReduce

https://github.com/BatzoglouLabSU/SIMLR

Vortex: http://web.stanford.edu/~samusik/vortex/

# Differential Expression

Monocle cole-trapnell-lab.github.io/monocle-release/

scDD: https://github.com/kdkorthauer/scDD

ISOP: comparison of isoform pairs in single cells https://github.com/nghiavtr/ISOP

D3E: http://hemberg-lab.github.io/D3E/

BASiCS: https://github.com/catavallejos/BASiCS

Beta Poisson: https://github.com/nghiavtr/BPSC

# Time-series/ordering/lineage prediction

Monocle

Analysis of pseudotime uncertainty: http://biorxiv.org/content/biorxiv/early/2016/04/05/047365.full.pdf

ECLAIR: cell lineage prediction https://github.com/GGiecold/ECLAIR

Identification of ordering effects: https://github.com/lengning/OEFinder

Slicer: non-linear trajectories https://github.com/jw156605/SLICER

Wishbone: identification of bifurcations in developmental trajectories http://www.c2b2.columbia.edu/danapeerlab/html/cyt-download.html

SCOUP: https://github.com/hmatsu1226/SCOUP

Ouija: https://github.com/kieranrcampbell/ouija

# Pipelines

Seurat http://www.satijalab.org/seurat.html

SINCERA https://research.cchmc.org/pbge/sincera.html

MAST: https://github.com/RGLab/MAST

scde (differential expression + gene set over-dispersion): https://github.com/hms-dbmi/scde

BaSiCs: Bayesian analysis of single cell data: https://github.com/catavallejos/BASiCS

FastProject: https://github.com/YosefLab/FastProject/wiki

Citrus: http://chenmengjie.github.io/Citrus/

Tools from Teichmann lab (cellity, celloline, scrnatb): https://github.com/Teichlab/

SCell: https://github.com/diazlab/SCell
# Other

Ginko: analysis of CNVs in single-cell data: http://qb.cshl.edu/ginkgo/?q=/XWxZEerqqY477b9i4V8F

CNV calling: http://genome.cshlp.org/content/early/2016/01/15/gr.198937.115.full.pdf

Gene co-expression: http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004892

DNA SNV calling: https://bitbucket.org/hamimzafar/monovar

# Methylation

Prediction of missing information: https://github.com/cangermueller/deepcpg
