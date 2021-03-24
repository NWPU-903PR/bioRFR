# bioRFR: A MATLAB(R) Package for biological System Resilience Function Reconstruction
*by Yan Li (linay0124@outlook.com), Shaowu Zhang\* (zhangsw@nwpu.edu.cn)*

## Introduction
In this repository, we provide a MATLAB(R) package for reconstrusting resilience function of biological precesses from gene expression data and the code about biological analysis and figures plotting associated with the paper: *Resilience function uncovers the critical transitions in cancer initiation*.

## Data Available
The data used in our paper is available at [Google Drive](https://drive.google.com/drive/folders/11VDCpGKDCT644WsrMZLJwQU2ChOEqzEc?usp=sharing). And the detail information can be found in [README](https://github.com/NWPU-903PR/bioRFR/blob/master/Data/README.md) file.

## Directory Tree
```
├─Algorithm                   [the code of biRFR framework]
| ├─Parameter_calculate.m       [the function used to calculate status and resilience parameters]
| ├─Estimate_landscape.m        [the function used to estimate potential landscape]
| ├─Resilience_centrality.m     [the function used to calculate resilience centrality for each genes]
| ├─computeGrid.m               [a internal function called by Estimate_landscape()]
| └─example.m                   [a sample script illustrates the workflow of bioRFR]
├─BiologicalAnalysis          [the code of biological analysis and figures plotting]
| ├─PlotResilienceFunction.R    [plotting the resilience function shown in Fig.3 in main text of paper]
| ├─PredictSurvivalTime.R       [analyzing the relationship between status parameter and survival time (Fig. 4)]
| ├─SurvivalAnalysis_KGs.R      [survival Analysis based on KGs (Fig.5)]
| └─TopologicalDistance.R       [calculating topological distance between iKGs and CGC genes]
├─Data                        [the detail information about available data]
```

## How to use bioRFR
We provide a sample scrept [`Algorithm/example.m`](https://github.com/NWPU-903PR/bioRFR/blob/master/Algorithm/example.m) which contains the complete workflow of bioRFR. In brief, it consists with the following steps.

**1. Data preparation**

bioRFR needs two gene expression matrix whose rows and columns represent genes and samples respectively. One of them contains reference samples (`<exp_ref>`) and the other consists with the remaining samples (`<exp_other>`).

**2. Calculating status and resilience parameters**

One can call the following function to calculate parameters for every samples:
```matlab
[x_effN,beta_effN,x_effT,beta_effT] = Parameter_calculate(<exp_ref>,<exp_other>,th);
```

**3. Estimation potential landscape**

Based on the output of step 2, the potential landscape and corresponding minima and maxima can be estimate using `Estimate_landscape()` function.
```matlab
[minlocs,maxlocs] = Estimate_landscape(x_effN,beta_effN,x_effT,beta_effT,peakth);
```

**4. Calculation resilience centrality of each gene for every samples**

In order to assess the importance of each gene for state transition of individual patient, the resilience centrality of each gene can be calculated using `Resilience_centrality()` function.
```matlab
[resN_beta,resT_beta] = Resilience_centrality(expN,expT,th);
```

**5. Downstream biological analysis**
Based on resilience function and personalized key genes (KGs), one can perform downstream biological analysis (see `BiologicalAnalysis\` as some examples).

