# Data Available

The following data is available at [Google Drive](https://drive.google.com/drive/folders/11VDCpGKDCT644WsrMZLJwQU2ChOEqzEc?usp=sharing).
```
├─Clinical_info                    [containing clinical information for BRCA, LUAD and LUSC, source from TCGA]
├─Intermediate                     [containing the output of each step of bioRFR]
| ├─Landscape_Results                [the output of Estimate_landscape() for six type of cancer]
| ├─Parameters_Results               [the output of Parameter_calculate() for six type of cancer]
| └─Resilience_Centrality_Results    [the output of Resilience_centrality() for BRCA, LUAD and LUSC]
├─Others
| ├─cancer_gene_census.csv           [Known driver genes from CGC]
| └─fullNetwork.rda                  [the gene regulation network]
├─expN_<cancer_type>.csv           [the gene expression matrix of normal samples for specific <cancer_type>]
└─expT_<cancer_type>.csv           [the gene expression matrix of tumor samples for specific <cancer_type>]

```
