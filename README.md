# cnvTree-scCNVcluster: scDNA-seq data cell clustering

## Overview

\

## Installation

```
    # In R version 4.2.1 is suitable
    install.packages('devtools')
    devtools::install_github('39652its/cnvTree_scDNAseqclustering')
```
\

## Inputs to cnvTree-scCNVcluster

### Copy number variation (CNV) data in single-cell DNA sequencing data

Two types of copy number variation (CNV) data for scCNVcluster is available.

\

1. CNV data output from R package Aneufinder

```
    # Installation 
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("AneuFinder")

    # Running Aneufinder
    # The following call produces plots and genome browser files for all BAM files in "my-data-folder"
    Aneufinder(inputfolder="my-data-folder", outputfolder="my-output-folder")
```

\

2.  Other CNV calling tools are available. A data frame of CNV data in either `.rds` or `.txt` file with the following required columns:\

    The first row of data frame is column name, please match the name with following required columns, with additional columns are acceptable. 
    Also, the data frame do not need the column indexes.

  -   cellID: Unique identifier for each cell.
  
  -   seqnames: Chromosome or sequence name.
  
  -   start: Start position of the segment.
  
  -   end: End position of the segment.
  
  -   copy.number: Copy number value for the segment.\

```
    # The parser for data from other CNV calling tools in dataframe with the following required columns

    Sample <- changeFormat(file = file_path)
```

\

### Cytoband information

A table for cytoband information seen on Giemsa-stained chromosomes. You could pick the suitable cytoband template on UCSC database at following download site.\
**hg38, GRCh38**: <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz>\
**hg19, GRCh37**: <https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz>\

\

## Running cnvTree-scCNVcluster

### Quick start

A single call allows the execution of the entire analysis of single-cell DNA cell clustering.

```
    UCSC_cytoband_file_path <- system.file("extdata", "hg38_cytoBand.txt.gz", package = "cnvTree")

    cnvTree_scDNAclustering(input = Sample,
                            pqArm_file = UCSC_cytoband_file_path)
```                            

### Step-by-step execution

Or you could exectue by steps. The analysis of scDNA cell clustering workflow, could be categorized into 4 steps:\

1.  **pqArm clustering step**: Grouping cells based on chromosomal arm-level (p/q arm) copy number variations, and dividing copy numbers into three states (Deletion, Neutral, and Amplification).\

2.  **Re-clustering step**: Re-clustering clusters based on differences in chromosomal arm-level copy number variations. Clusters were compared between clusters at the bin level, identifying and merging the most similar clusters to reduce cell exclusion.

    -   `difratio_chr`: A numeric value defining the threshold for acceptable difference ratios across different chromosomes. 
                        The higher value retains more cells for further analysis.\

3.  **Subclone clustering step**: Refining clusters based on copy number variations (CNVs) within subpopulations of cells. Within each cluster, significant breakpoints were identified based on absolute copy number values, delineating event regions and identifying segments with highly similar copy number ranges. These segments were used to further subdivide the clusters into subclones using segment CN sequences.

    -   `overlap_region`: An integer representing the genomic region where copy number frequently changes. 
                          The higher value masked out more region in calculation, for more rough clustering calculation. 
                          If the sample with more complicated CNV pattern, suggesting the the higher value to be setted.
    -   `dif_ratio`: A numeric value defining the tolerance threshold for copy number differences between cells within a cluster.

4.  **Output final results step**: Integrating clustering results and generating visualizations. The scDNA-seq data cell clustering approach used absolute copy number values to identified high-confidence CNV regions between clusters, and the presence of corresponding CNVs for each cluster in these regions was calculated for tumor evolutionary reconstruction.

    -   `consecutive_region`:
  
    -   `cellcutoff`:

```
    # Inputs

    UCSC_cytoband_file_path <- system.file("extdata", "hg38_cytoBand.txt.gz", package = "cnvTree")

    # Steps for scDNA cell clustering

    pqArm_result <- pqArmClustering(input = Sample, pqArm_file = UCSC_cytoband_file_path) 
    Reclustering_output <- Reclustering(input = Sample, pqArm_output = pqArm_result, pqArm_file = UCSC_cytoband_file_path) 
    Subclone_output <- SubcloneClustering(input = Example_data, Reclustering_output = Reclustering_output)

    # Outputs the cell clustering result

    scDNA_Output(input = Example_data, Summary = Subclone_output, pqArm_file = UCSC_cytoband_file_path)
```

## Outputs from cnvTree-scCNVcluster

The results of CNV-based clustering in scDNA-seq data will include the following output files:\
1. **cnvTree.scDNAseq_DefinedCNVregion.txt**: This file records all genomic regions where CNVs (Deletion/Amplification) have occurred in the sample and details the changes in copy number.\
2. **cnvTree.scDNAseq_DNAcluster.txt**: This file documents whether the clusters identified in the sample exhibit corresponding CNV changes within the defined CNV regions.\
3. **cnvTree.scDNAseq_grouping.txt**: This file records the clustering process at various stages, including the cell ID for each cluster, allowing researchers to study the composition of cells within each group.\
4. **cnvTree.scDNAseq_SubcloneRegionCN.txt**: This file details the final subclone clusters, dividing the genome into segments based on each event region and smoothing the copy number within these segments.\
5. **cnvTree.scDNAseq_Grouping_fig.pdf**: This file presents a heatmap of the final clustering results, showing the copy number profile for each cluster in the sample.\
6. **cnvTree.scDNAseq_heatmap.png**: This figure presents a clustering result between final clusters in CNV pattern, showing the copy number profile with dendrogram.

\

## Reference

\

## Contact

[twc@nycu.edu.tw](twc@nycu.edu.tw) and [it35101.bt11@nycu.edu.tw](it35101.bt11@nycu.edu.tw)
