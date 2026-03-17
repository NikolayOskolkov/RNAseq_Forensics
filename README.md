<p align="center">
  <img src="images/RNAseq_Forencics.png" alt="RNAseq Forensics Logo" height="320"/>
</p>

# RNAseq Forensics
This is a computational machine learning method for detecting unwanted tissue signals in RNAseq samples.

## Quick start
To predict contamination, you need to provide a gene expression matrix `merged_tissues_organisms.txt` and run the following command line:

    Rscript RNAseq_Tissue_Predictor.R merged_tissues_organisms.txt

The output of the tool is the quantification of transcriptomic contribution from different tissues to each RNAseq sample:

<p align="center">
  <img src="images/DecisionTree_validation.png" alt="RNAseq Forensics Prediction" height="420"/>
</p>

The tool was validated on five independent RNAseq gene expression datasets originating from muscle, liver, heart and blood tissues from human, mouse and apes.

## Train data download
The GTEX v.7 gene expression matrix caan be downloaded from: https://drive.google.com/drive/folders/1AuHLbv_VqVgLmz9GptN80adjjbUaNj4I?usp=drive_link
