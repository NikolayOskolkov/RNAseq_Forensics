<p align="center">
  <img src="RNAseq_Forencics.png" alt="RNAseq Forensics Logo" height="320"/>
</p>

# RNAseq Forensics
This is a computational machine learning method for detecting unwanted tissue signals in RNAseq samples.

## Quick Start
To predict contamination, you need to provide a gene expression matrix `merged_tissues_organisms.txt` and run the following command line:

    Rscript RNAseq_Tissue_Predictor.R merged_tissues_organisms.txt

The GTEX v.7 gene expression matrix caan be downloaded from: https://drive.google.com/drive/folders/1AuHLbv_VqVgLmz9GptN80adjjbUaNj4I?usp=drive_link
