# DNA methylation analysis guides

This lab contains scripts and instructions for analyzing DNA methylation data using reduced representation bisulfite sequencing (RRBS) datasets. The analyses cover data preprocessing, quality control, alignment, methylation calling, and downstream analyses such as differential methylation analysis and visualization.

## Datasets

- Original paper: [Genome-wide DNA methylation profiling reveals cancer-associated changes within early colonic neoplasia](https://www.nature.com/articles/onc2017130) from **oncogene** journal published in Jan 23, 2018
- BioProject ID: [PRJNA377851](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA377851)
- GEO accession: [GSE95656](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95656)
- This subset consists of three colorectal cancer (CRC) samples and three matched normal colon tissue samples, restricted exclusively to reads from chromosome 16 to meet computational resource requirements in google colab. Consequently, these results are not intended to replicate the original study or represent the full biological complexity of the complete dataset. Please consider these constraints when performing the analysis and interpreting the final results.
- Subset samples include:
  - [SRR5311114](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR5311114): RRBS on adjacent normal colon, patient #15
  - [SRR5311115](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR5311115): RRBS on CRC, patient #15
  - [SRR5311116](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR5311116): RRBS on adjacent normal colon, patient #16
  - [SRR5311117](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR5311117): RRBS on CRC, patient #16
  - [SRR5311124](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR5311124): RRBS on CRC, patient #60
  - [SRR5311125](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR5311125): RRBS on adjacent normal colon, patient #60

## Analysis workflow

The analysis workflow is organized into the following Google Colab notebooks: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1MHwkoO7Jv0eNVaOkkU8Kd5gPxFWs4gs2?usp=sharing)
