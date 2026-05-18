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

### 01_scripts

### 02_data

### 03_output


## Reproducing Results


## Ethics and Privacy

## Licence

## Permissions and Access 