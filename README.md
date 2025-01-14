
# Spatial profiling of chromatin accessibility in formalin-fixed paraffin-embedded tissues

## Introduction
This repository aims to share the raw data processing and visualization code used in the **"Spatial profiling of chromatin accessibility in formalin-fixed paraffin-embedded tissues"** paper.




![github_spatial_FFPE_ATAC](https://github.com/user-attachments/assets/25175ed9-11ca-4c47-b376-86a0cd320208)


## Data analysis
### 1. Preprocessing the sequencing data
 Next Generation Sequencing (NGS) was performed using the Illumina NovaSeq sequencer (paired-end 150 bp mode). 
 
In the Data_preprocessing folder, directories beginning with Snakemake_* contain the code for preprocessing different modalities. The preprocessing pipeline utilizes raw FASTQ data as input, where Read 1 comprises genomic sequences, and Read 2 contains the spatial barcodes, specifically Barcode A and Barcode B.

**The preprocessing pipeline we developed using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system is in the Data_preprocessing folder. After putting the input files in the correct directory, to run the pipeline, use the command:**

    sbatch Snakemake.sh


**Brief descriptions of the preprocessing pipeline in Snakefile:**

(1) Directory and File Setup
- Automates the creation of directories for storing raw and processed data per sample.
- List samples dynamically based on the provided raw data directory.

(2) Reads Filtering
- `filter_primer`: Utilize [bbduk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) to filter sequences with primers from the reads.
- `filter_L1` & `filter_L2`: Further filtering steps target specific linker sequences.

(3) Barcode Processing
- `bc_process`: Extract and reformat the data. (BC_process.py)
- `R1_rename`: Rename and reorganize the processed reads for consistency and further processing.

(4) Barcode Tagging and Matching
- `taggd`: Correct barcodes and prepare the data for demultiplexing.
- `add_BC`: Integrate barcode information into sequencing reads based on a pre-generated list of barcodes and their matches to specific reads. (add_BC.py)

(5) Adapter Trimming
- `adapter_trim`: Trims sequencing adapters from the reads using [Trim Galore](https://github.com/FelixKrueger/TrimGalore), preparing them for alignment.

(6) Sequence Alignment
- `bwa_align`: Map reads to a reference genome with [BWA](https://github.com/lh3/bwa), followed by sorting and indexing the alignments using [samtools](https://www.htslib.org/).

(7) Fragment File Generation
- `fragments`: Transforms BAM files into sorted, compressed, and indexed BED files using [sinto toolkit](https://timoast.github.io/sinto/), for downstream analysis. 




###  2. Preprocessing the microscope image

**Identify pixels on tissue from microscope image using Python**

See the files in Image_preprocess under Data_preprocessing folder.



### 3. Downstream analysis
Downstream analyses were completed with R and Python. 
The [SnapATAC2](https://github.com/kaizhang/SnapATAC2) was conducted for normalization and dimension reduction.
The R package used extensively the functions in [Seurat](https://github.com/satijalab/seurat) v.4.3.0.1, [ArchR](https://github.com/GreenleafLab/ArchR) v1.0.2. 

**Brief descriptions of analysis scripts:**
[SnapATAC2](https://github.com/kaizhang/SnapATAC2): Analysis of ATAC data.
[ArchR](https://github.com/GreenleafLab/ArchR): Analysis and visualization of ATAC data.


**Functions:**

spatial_data_visualization.R: Visualize spatially resolved data on tissue sections.


## References

Pengfei Guo, Yufan Chen, Liran Mao, Angelysia Cardilla, Chin Nien Lee, and Yanxiang Deng. "Spatial profiling of chromatin accessibility in formalin-fixed paraffin-embedded tissues."
