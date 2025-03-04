# scDNAcluster

## Overview

Single-cell DNA sequencing (scDNA-seq) analyzes the limited DNA content of individual cells to detect copy number variations (CNVs). CNVs, which involve deletions or duplications of genes or chromosome arms, are common in tumors and can activate oncogenes or inactivate tumor suppressor genes, indirectly influencing gene expression. scDNAcluster is a method that leverages CNV data from scDNA-seq to cluster cells based on their CNV patterns. This approach provides cell clustering results and identifies high-confidence CNV regions, aiding in the reconstruction of tumor evolutionary history.

![image](https://github.com/39652its/cnvTree_scDNAseqclustering/blob/main/Github_scDNA_Overview.png)
 
## Installation

```
    # This package is designed for R version 4.2.1.
    install.packages('devtools')
    devtools::install_github('39652its/cnvTree_scDNAseqclustering')
```

## Inputs to scDNAcluster

### Copy number variation (CNV) data in single-cell DNA sequencing

scDNAcluster supports two types of CNV data inputs:

1. CNV data from **AneuFinder**

    CNV data generated by the R package **AneuFinder** can be used directly. To install and run AneuFinder, follow these steps:

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
    
    CNV data from other tools can also be used, provided the input is a data frame saved as either an .rds or .txt file. The data frame **must** contain the following required columns (additional columns are allowed):

      -   cellID: Unique identifier for each cell.
      -   seqnames: Chromosome or sequence name.
      -   start: Start position of the segment.
      -   end: End position of the segment.
      -   copy.number: Copy number value for the segment.

    The first row of the data frame should contain the column names exactly as listed above. Row indices are not required.
    
```
    # Converting CNV data from other tools into the required format 
    file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
    Sample <- changeFormat(file = file_path)
```


### Cytoband information
  The parameter `pqArm_file` requires a table containing cytoband information, which corresponds to Giemsa-stained chromosome regions.

1.  Built-in cytoband templates

    scDNAcluster provides four built-in cytoband templates, matching the required names. Below are their respective UCSC database download links for reference:

      -   **hg38**: [Download](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz)
      -   **hg19**: [Download](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz)
      -   **mm39**: [Download](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/cytoBand.txt.gz)
      -   **mm10**: [Download](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/cytoBand.txt.gz)
    
2.  Custom cytoband template

    Alternatively, you can provide a custom cytoband file for analysis. Specify the **file path** of your input file, which must follow the UCSC cytoband format and include the following columns:

      -   chrom: Reference sequence chromosome or scaffold.
      -   chromStart: Start position in genoSeq.
      -   chromEnd: End position in genoSeq.
      -   name: Name of cytogenetic band.
      -   gieStain: Giemsa stain results.

    For more details on the required format, refer to the [UCSC Genome Browser Table Schema](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=map&hgta_track=cytoBand&hgta_table=cytoBand&hgta_doSchema=describe+table+schema).


## Running scDNAcluster

### Quick start

You can execute the entire single-cell DNA clustering analysis with a single function call:

```
    # Input (Example using demo data)  
    file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
    Sample <- changeFormat(file = file_path)
    
    # Single-call execution
    cnvTree_scDNAclustering(input = Sample,
                            pqArm_file = "hg38")
```                            

### Step-by-step execution

Alternatively, you can run the analysis step by step. The scDNAcluster workflow consists of four main steps:

1.  **pqArm clustering step**: Cells are grouped based on chromosomal arm-level (p/q arm) copy number variations. Copy number values are categorized into three states: Deletion, Neutral, and Amplification.

2.  **Re-clustering step**: Clusters are re-analyzed by comparing different chromosomal arm-level CNVs at the bin level. The most similar clusters are merged to minimize cell exclusion.

    -   `difratio_chr`: A numeric threshold for acceptable CNV differences across chromosomes, where higher values
                        retain more cells for further analysis. If the sample has a complex CNV pattern, increasing 
                        this value is recommended. A reasonable range is between 0 and 0.5.

3.  **Subclone clustering step**: Clusters are further refined based on subpopulation-specific copy number variations (CNVs). Significant breakpoints are identified within each cluster based on absolute copy number values, defining event regions and segmenting CN sequences for subclone classification.

    -   `overlap_region`: A numeric value determines how much genomic space is masked when identifying frequent CNV
                          regions. A higher value results in coarser clustering, which is useful for complex CNV patterns.
    -   `dif_ratio`: A numeric threshold for tolerating CNV differences within a cluster.
                    Increasing this value allows more cells to be grouped together, making clustering coarser. 
                    If the sample has a complex CNV structure, a higher value is advisable. 
                    The recommended range is between 0 and 0.5.

4.  **Output final results step**: The final clustering results are integrated, and visualizations are generated. This step identifies high-confidence CNV regions and reconstructs tumor evolution based on CNV presence.

    -   `consecutive_region`:  A numeric threshold defining the minimum CNV region length to be retained in the final output.
                              A higher value filters out shorter CNV regions.
    -   `cellcutoff`: Minimum number of cells required for a cluster to be included in the final output.
                      A higher value increases clustering confidence. The reasonable value is > 2 to define a cluster.
                      If the sample has a complex CNV pattern, a lower value is recommended to retain more cells.

```
    # Input (Example using demo data)  
    file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
    Sample <- changeFormat(file = file_path)
    
    # Step-by-step execution  
    pqArm_result <- pqArmClustering(input = Sample, pqArm_file = "hg38") 
    Reclustering_output <- Reclustering(input = Sample, pqArm_output = pqArm_result, pqArm_file = "hg38") 
    Subclone_output <- SubcloneClustering(input = Example_data, Reclustering_output = Reclustering_output)

    # Output cell clustering results
    scDNA_Output(input = Example_data, Summary = Subclone_output, pqArm_file = "hg38")
```

## Outputs from scDNAcluster

The results of CNV-based clustering in scDNA-seq data will include the following output files:

1. **cnvTree.scDNAseq_DefinedCNVregion.txt**: This file contains a list of genomic regions where CNVs have been detected in the sample, along with details on copy number changes (deletions or amplifications).

2. **cnvTree.scDNAseq_DNAcluster.txt**: This file indicates whether each identified cluster exhibits CNV changes within the defined CNV regions. The presence of CNVs is recorded as 1, while the absence is recorded as 0.

3. **cnvTree.scDNAseq_grouping.txt**: This file tracks the clustering process across different stages, specifying which cluster each cell was assigned to.

    - To examine clustering results at any step, select the column `"STEPname"+"_cluster"` along with `"cellID"`.
    - To visualize clustering results, use the Totalcluster_pdf() function, specifying the relevant `step` name. The available `step` names are `"pqArm"`, `"Recluster"`, and `"Subclone"`, with `"Subclone"` set as the default.


4. **cnvTree.scDNAseq_SubcloneRegionCN.txt**: This file details the final subclone clusters, segmenting the genome into regions based on detected CNV events and smoothing copy number values within each segment.

5. **cnvTree.scDNAseq_Grouping_fig.pdf**: This PDF provides a heatmap visualization of the final clustering results, displaying the copy number profiles for each cluster. 

    -  If you wish to generate figures for specific clusters, first select the relevant cells from **cnvTree.scDNAseq_grouping.txt**, then use the function `GenomeHeatmap()` to visualize the CNV profile of the target cell group.

6. **cnvTree.scDNAseq_heatmap.png**: This figure presents a heatmap illustrating the clustering relationships between the final clusters based on CNV patterns. A dendrogram is included to depict hierarchical relationships among clusters.


## Reference



## Contact

[twc@nycu.edu.tw](twc@nycu.edu.tw) and [it35101.bt11@nycu.edu.tw](it35101.bt11@nycu.edu.tw)
