# Flexible Algorithms and Latent Variables for Baseline Covariate Adjustment in Cross-Lagged Panel Models

## Introduction
This repository contains the R scripts and data to reproduce the results of the master's thesis *Flexible Algorithms and Latent Variables for Baseline Covariate Adjustment in Cross-Lagged Panel Models (CLPMs)*. 

In CLPMs, baseline confounders can be addressed explicitly as covariates or implicitly by latent variables. The former relies on observing and correctly modeling confounder effects, while the latter relies on the assumption of stationarity. Baseline Covariate Adjustment Structural Equation Modeling (BCA-SEM) combines both approaches, enabling the confounding functional form to be flexibly modeled. In this strategy, SEMs are fitted to residuals obtained from regressing the variables X and Y on observed covariates. This study examines whether a combined approach yields more accurate estimates of cross-lagged parameters in finite samples compared to models commonly used in psychology.  

Four simulation scenarios cumulatively introduce violations of the linearity assumption (s2), the stationarity assumption (s3) and the no unobserved confounding assumption (s4). Performance is compared for five models, namely:

1. The CLPM including covariates.
2. The CLPM fit to residuals based on Extreme Gradient Boosting (XGB) regression.
3. The Random Intercepts (RI) CLPM.
4. The RI-CLPM including covariates. 
5. The RI-CLPM fit to residuals based on XGB regression.

All scenarios are evaluated at sample sizes of 300, 1000, and 2000. 

## Contents 

This repository has the following structure. The `00_ignore` folder contains unordered project files, and may be ignored. 

```text
.
├── 01_scripts
│   ├── 01_data_generation
│   └── 02_plots_tables
│       ├── 01_thesis_results
│       └── 02_additional_results
├── 02_data
│   ├── 01_thesis_results
│   └── 02_additional_results
├── 03_output
├── renv
│   └── activate.R
├── .Rprofile
├── LICENSE
├── README.md
└── renv.lock
```

### `01_scripts`

This folder contains the scripts to generate the simulation results and to analyse these results. It is divided into two subfolders.  

#### `01_data_generation` 

This folder contains the scripts used to run the simulation study. All scripts are called from `12_simulation_runner.R`. The setup generates data under different simulation scenarios and efficiently applies a range of models to the generated data. The current version of `12_simulation_runner.R` is configured to reproduce the results presented in the thesis. For a more detailed explanation of the simulation scripts and their functionality, see the README included in this folder.

#### `02_plots_tables`

This folder is divided in two subfolders. `01_thesis_results`  contains the scripts that can be used to produce the output figures and tables for the thesis project. `02_additional_results` contains functions that can be used to explore any dataframe produced by the simulation runner. For further details, see the README included in this folder. 

### 02_data

This folder is again divided in two subfolders. `01_thesis_results` contains the ouput data from running the setup in `01_data_generation`, which is supplied here since the simulation study is computationally expensive. `02_additional_results` contains the output data from other specifications of the simulation runner. We refer the interested user to the included README included in this folder. 

### 03_output

This folder contains the four output plots and improper solutions table from the thesis. These are the results of requentially running the simulation study and plotting + tables scripts. 


## Reproducing Results

### Prerequisites

To reproduce the analyses, the following is needed.

* R (>4.5.2)
* The `renv` package (1.2.3)

> **Hardware note**
>
> The full simulation study was run on a high-performance machine with the following specifications:
>
> - CPU: 112 × Intel(R) Xeon(R) Platinum 8580
> - Threads: 224
> - Memory: 851 GB RAM
>
> These resources were used to generate the simulation results efficiently. Creating plots from simulation results can be done on any modern pc or laptop. 



## Ethics and Privacy

## Licence

## Permissions and Access 