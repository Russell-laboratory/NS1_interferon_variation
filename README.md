# NS1_interferon_variation

The code contained within this repository describe the analysis of data, and generation of figures, in our attempts to understand how variation within a viral population lacking NS1 intersects with the induction of an interferon response.

Analysis are split across several different jupyter notebooks, with descriptions below, and figures inline.

All intermediate files less than 10mb were pushed to this repository to aid in understanding our analysis.
Larger files, such as sequencing, will need to be downloaded and reconstructed in the apppropriate folders in order to re-run this analysis.
Specifically, place NGS files within appropriate directories within Sequencing to rerun all sequencing analyses.

Data may be found in GEO (including some processed files) at accession GSE215914.

To rerun single-cell analyses, 10x genomics cell ranger must be first run to generate the appropriate files.

## Directories

The following directories, and their purpose, exist within this repository.

- <b>Database</b>       Repository for influenza genomic sequences. A/WSN/1933 BLAST database and STAR indicies provided. 
- <b>Results</b>        Final datafiles for analyses after processing. Most provided as supplemental files within manuscript.
- <b>Scripts</b>        Short scripts written for this analysis. Seperated from jupyter notebooks for readability and portability.
- <b>Sequencing</b>     Folder that would contain NGS samples and 10x genomics output to regenerate this pipeline. Folder achitecture essential to run off pipeline.
- <b>flowCytometry</b>  Flow cytometry raw data
- <b>Data</b>           qPCR, flow cytometry, and titer data
- <b>Kelly_et_al_2022_reanalysis</b> Reanalysis of Kelly et al. 2022 data. Typo where donor 1905 is actually 1904, so transpose for paper.
  

## Dependencies

Code within this repository was run with the following tools, and versions, installed and available from PATH. Specific websites and documentation are provided where available. 

- <b>Python</b>      run with version 3.7. Available from https://www.python.org/downloads/release/python-370/
- <b>Trimmomatic</b> run with version 0.39. Available from http://www.usadellab.org/cms/?page=trimmomatic
- <b>STAR</b>        run with version 2.7.1.a. Available from https://github.com/alexdobin/STAR
- <b>Samtools</b>    run with version 1.9. Available from http://www.htslib.org/
- <b>FastQC</b>      run with version 0.11.8. Available from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- <b>BLASTn</b>      run with version 2.9.0+. Available from https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
- <b>R</b>           run with version 4.0.2. Available from https://cran.r-project.org/bin/windows/base/
-<b>cellranger</b>   run with version 1.3.1. Available from https://www.10xgenomics.com/


The following python packages and versions were used. All were installed using conda. (https://docs.conda.io/en/latest/)
- <b>numpy</b>       run with version 1.19.5. (https://numpy.org/)
- <b>matplotlib</b>  run with version 3.3.3. (https://matplotlib.org/)
- <b>seaborn</b>     run with version 0.11.0. (https://seaborn.pydata.org/)
- <b>pandas</b>      run with version 1.2.0.(https://pandas.pydata.org/)
- <b>scipy</b>       run with version 1.6.0. (https://www.scipy.org/)
- <b>statsmodels</b> run with version 0.12.1. (https://www.statsmodels.org/stable/index.html)

The following R packages and versions were used. All were installed using CRAN (https://cran.r-project.org/)
- <b>dplyr</b>      run with version 1.0.9 (https://dplyr.tidyverse.org/)
-<b>patchwork</b>    run with version 1.1.1 (https://patchwork.data-imaginist.com/)
-<b>ggplot2</b>     run with version 3.3.6 (https://ggplot2.tidyverse.org/)
-<b>DropletUtils<\b> run with version 1.16.0 (https://bioconductor.org/packages/release/bioc/html/DropletUtils.html)
-<b>SoupX<\b> run with version SoupX 1.6.2 (https://cran.r-project.org/package=SoupX)

## Jupyter notebooks

General descriptions of pipelines within each notebook described below.

### SingleCellPython_thresholds.ipynb 

Notebook describing the empirical setting of influenza and interferon thresholds in our single-cell RNAseq datasets, as well as statistical analyses as presented in the paper. Uses a similar methodology to Russell et al. 2019 J. Virol., with the exception that rather than a flat percentage of contamination inferred directly from canine cells, co-emulsions were used to estimate the relative RNA content difference between human and canine cells in order to infer an adjusted contamination percentage. Additionally includes statistical analysis of deletions in single-cell RNAseq as identified in Single_cell_deletions.ipynb

### Single_cell_deletions.ipynb

Identification of deletions in single-cell RNAseq. Uses a similar method as Hamele et al. 2022 J. Virol. Remaps gapped reads from STAR output using BLAST and then calls deletions using a custom script. Additionally records non-gap spanning reads to infer mixed deletion/non-deletion infections.

### BulkSequencingMapping

Similar to Mendes et al. 2020 PLoS Pathogens.
Processing of natural diversity vRNA sequencing.Prior to any additional analysis, we first present our general quality metrics using fastQC and a custom script to parse results.

First, as all samples were tagmented using Nextera adapters, Trimmomatic is used to remove these sequences as well as generally process our fastQC files. Again, a custom script is used to parse output for readability in our notebook.

All reads are then mapped using STAR. Values were chosen to enforce ungapped mapping to the A/WSN/1933 genome. Unmapped reads are retained as seperate fastQ files.

For unmapped reads, files were converted to fasta files and run through BLASTn. Discontinuous junctions were identified wherein the discontinuity matched to the same segment, with the same polarity. These junctions were then used to initialize a new GTF, and generate a new STAR index, one for each biologically distinct viral population. 

Unmapped reads were then mapped using these new STAR indices, enforcing matching to annotated discontinuities only. The resulting four bamfiles (2 for read1, 2 for read2) were merged and flags were fixed to identify pairmates. To count the number of occurrences of a given deletion junction, bamfiles data were sorted on pairmates, and, using awk, only reads where one or the other read has a discontinuity were retained for counting. Thereafter, it was ensured that each read was consistent with the deletion; that is, if the two reads overlap in the region containing the deletion, they must both call the same deletion. It is also required that between the two reads, at least 3 bases of sequence are mapped to either side of the junction.

Deletions within replicates were maintained, however additional deletions, unmapped by STAR yet identified by BLAST, in NS, were retained owing to their relative rarity for the analyses described in the manuscript.

Remaining portions of this notebook describe analysis as depicted in the manuscript, and generation of those figures. 

### SeuratUMAP.Rmd

R-markdown file describing Seurat analysis and UMAP clustering of single-cell RNAseq data.

### FlowCytometry.ipynb

General analysis of flow cytometry data after initial debris gate and exporting to csv datafiles. Thresholds set on uninfected controls, expected to have no production of interferons or staining for influenza A virus proteins.

### qPCR.ipynb

Analysis of qCPR data exported to csv files. Also includes titer measurements for analysis of genome/infectivity ratios. 

## Note

The analysis of deletion junctions by BLASTn is very memory-intensive, and the intermediate, uncompressed, files take up a lot of space. I would recommend against running those components of the pipeline anywhere where you don't have access to ~100GB of space and ~32GB of RAM. Other components of this analysis can be run with very modest computational resources.

