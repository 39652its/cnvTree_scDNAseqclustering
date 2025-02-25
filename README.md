# cnvTree-scCNVcluster: scDNA-seq data cell clustering

## Overview

Single-cell DNA sequencing (scDNA-seq) detects copy number variations (CNVs) by analyzing the limited DNA content of individual cells. CNVs, involving gene or chromosome arm deletions or duplications, are common in tumors and can activate oncogenes or inactivate tumor suppressor genes, indirectly affecting gene expression. scDNAcluster is a method that utilizes CNV data from scDNA-seq to cluster cells based on their CNV patterns. This approach ultimately outputs cell clustering results and identifies high-confidence CNV regions, facilitating the reconstruction of tumor evolutionary history.

![image](https://github.com/39652its/cnvTree_scDNAseqclustering/blob/main/Github_scDNA_Overview.png)
 
## Installation

```
    # This package is designed for R version 4.2.1.
    install.packages('devtools')
    devtools::install_github('39652its/cnvTree_scDNAseqclustering')
```

## Inputs to cnvTree-scCNVcluster

### Copy number variation (CNV) data in single-cell DNA sequencing

Two types of copy number variation (CNV) data for scCNVcluster can be used with scCNVcluster.

1. CNV data output from the R package **Aneufinder**.

    The CNV data output from AneuFinder can be used directly. You can install and run AneuFinder as follows:

```
    # Installation 
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("AneuFinder")

    # Running Aneufinder
    # The following call produces plots and genome browser files for all BAM files in "my-data-folder"
    Aneufinder(inputfolder="my-data-folder", outputfolder="my-output-folder")
```

2.  CNV data from other CNV calling tools
    
    CNV data from other tools can also be used. The input should be a data frame saved as either an `.rds` or `.txt` file. The data frame **must** include the following required columns (additional columns are acceptable):

      -   cellID: Unique identifier for each cell.
      -   seqnames: Chromosome or sequence name.
      -   start: Start position of the segment.
      -   end: End position of the segment.
      -   copy.number: Copy number value for the segment.

    The first row of the data frame should contain the column names, matching the required names listed above. Row indices are not required.
    
```
    # Parsing data from other CNV calling tools into the required format
    Sample <- changeFormat(file = file_path)
```


### Cytoband information
  The parameter `pqArm_file` requires a table containing cytoband information, which corresponds to Giemsa-stained chromosome regions.

1.  Built-in cytoband templates

    There are four built-in cytoband templates available for direct selection. Below are their respective download links from the UCSC database for reference:

      -   **hg38 (GRCh38)**: <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz>
      -   **hg19 (GRCh37)**: <https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz>
      -   **mm39 (GRCm39)**: <https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/cytoBand.txt.gz>
      -   **mm10 (GRCm38)**: <https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/cytoBand.txt.gz>
    
2.  Custom cytoband template

    Alternatively, you may provide a custom cytoband file for analysis. The input file should be formatted according to the UCSC database cytoband specifications and must include the following columns:

      -   chrom: Reference sequence chromosome or scaffold.
      -   chromStart: Start position in genoSeq.
      -   chromEnd: End position in genoSeq.
      -   name: Name of cytogenetic band.
      -   gieStain: Giemsa stain results.

    If you need to use a custom cytoband file beyond the provided options, please refer to the [UCSC Genome Browser Table Schema](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=map&hgta_track=cytoBand&hgta_table=cytoBand&hgta_doSchema=describe+table+schema) for more details on the required format.


## Running cnvTree-scCNVcluster

### Quick start

A single function call executes the entire single-cell DNA clustering analysis:

```
    UCSC_cytoband_file_path <- system.file("extdata", "hg38_cytoBand.txt.gz", package = "cnvTree")

    cnvTree_scDNAclustering(input = Sample,
                            pqArm_file = UCSC_cytoband_file_path)
```                            

### Step-by-step execution

Alternatively, you can execute the analysis step by step. The scDNA cell clustering workflow consists of four main steps:

1.  **pqArm clustering step**: Cells are grouped based on chromosomal arm-level (p/q arm) copy number variations. Copy number values are categorized into three states: Deletion, Neutral, and Amplification.

2.  **Re-clustering step**: Clusters are re-analyzed based on differences in chromosomal arm-level copy number variations. Cells are compared at the bin level, and the most similar clusters are merged to minimize cell exclusion.

    -   `difratio_chr`: A numeric threshold for acceptable difference ratios across different chromosomes.
                        A higher value retains more cells for further analysis. If the sample has a complex CNV
                        pattern, a higher value is recommended.

3.  **Subclone clustering step**: Clusters are further refined based on subpopulation-specific copy number variations (CNVs). Significant breakpoints are identified within each cluster based on absolute copy number values, defining event regions and segmenting CN sequences for subclone classification.

    -   `overlap_region`: An integer representing genomic regions where CNVs frequently occur.
                          A higher value masks more regions, resulting in coarser clustering.
                          If the sample has a complex CNV pattern, a higher value is recommended.
    -   `dif_ratio`: A numeric threshold for tolerating copy number differences within a cluster.
                    A higher value clusters more cells, leading to a coarser clustering. 
                    If the sample has a complex CNV pattern, a higher value is recommended.

4.  **Output final results step**: IThe final clustering results are integrated, and visualizations are generated. The scDNA-seq cell clustering method identifies high-confidence CNV regions across clusters and reconstructs tumor evolution based on CNV presence.

    -   `consecutive_region`:  A numeric threshold defining the minimum CNV region length to be retained in the final output.
                              A higher value filters out shorter CNV regions.
    -   `cellcutoff`: A numeric threshold defining the minimum number of cells required for a cluster to be included in the final output.
                      A higher value increases clustering confidence. 
                      If the sample has a complex CNV pattern, a lower value is recommended to retain more cells.

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

The results of CNV-based clustering in scDNA-seq data will include the following output files:

1. **cnvTree.scDNAseq_DefinedCNVregion.txt**: This file records all genomic regions where CNVs (Deletion/Amplification) have occurred in the sample, along with details of copy number changes.

2. **cnvTree.scDNAseq_DNAcluster.txt**: This file documents whether the identified clusters exhibit corresponding CNV changes within the defined CNV regions (1/0).

3. **cnvTree.scDNAseq_grouping.txt**: This file records the clustering process at various stages, including the cell ID assigned to each cluster.

    - If you want to examine clustering results at each step, select the column name `"STEPname"+"_cluster"` along with `"cellID"`.
    - To visualize clustering results at any step, use the function `Totalcluster_pdf()`. In `Totalcluster_pdf()`, specify the corresponding `cluster_name = "STEPname"+"_cluster"` and `cellnum_name = "STEPname"+"_cellnum"`. Notice that same `"STEPname"` should be specified.


4. **cnvTree.scDNAseq_SubcloneRegionCN.txt**: This file provides details of the final subclone clusters by segmenting the genome into regions based on detected events and smoothing the copy number values within each segment.

5. **cnvTree.scDNAseq_Grouping_fig.pdf**: This file a heatmap visualization of the final clustering results, displaying the copy number profile for each cluster in the sample. 

    - If users wish to generate figures for specific clusters, they can select the relevant cells from **cnvTree.scDNAseq_grouping.txt**. Next, using the function `GenomeHeatmap()` to visualize CNV profile of the target cell group.

6. **cnvTree.scDNAseq_heatmap.png**: This figure presents a heatmap that illustrates clustering relationships between the final clusters based on CNV patterns. This figure includes a dendrogram to depict hierarchical relationships among clusters.


## Reference



## Contact

[twc@nycu.edu.tw](twc@nycu.edu.tw) and [it35101.bt11@nycu.edu.tw](it35101.bt11@nycu.edu.tw)
