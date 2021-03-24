# bioRFR: A MATLAB(R) Package for biological System Resilience Function Reconstruction
*by Yan Li (linay0124@outlook.com), Shaowu Zhang\* (zhangsw@nwpu.edu.cn)*

## Introduction
In this repository, we provide a MATLAB(R) package for reconstrusting resilience function of biological precesses from gene expression data and the code about biological analysis and figures plot associated with the paper: *Resilience function uncovers the critical transitions in cancer initiation*.

## Data Available
The data used in our paper is available at [Google Driver](https://drive.google.com/drive/folders/11VDCpGKDCT644WsrMZLJwQU2ChOEqzEc?usp=sharing). And the detail information can be found in [README](https://github.com/NWPU-903PR/bioRFR/blob/master/Data/README.md) file.

## How to use bioRFR

- Gene expression datasets for six cancer types: `expN_<cancer_type>.csv` and `expT_<cancer_type>.csv`
- Status and resilience paramaters for each samples: `x_effN_<cancer_type>.csv`, `beta_effN_<cancer_type>.csv`, `x_effT_<cancer_type>.csv` and `betaT_eff_<cancer_type>.csv`
- The maxima and minima value of landscape for each cancer type: `maxlocs_<cancer_type>.csv` and `minlocs_<cancer_type>.csv`
- Resilience Centrality of each gene for every samples: `resN_beta_<cancer_type>.csv` and `resT_beta_<cancer_type>.csv`
