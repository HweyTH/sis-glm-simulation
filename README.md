# SIS-GLM Simulation Study

This repository contains code for reproducing the simulation study from the paper "Sure Independence Screening for Ultrahigh Dimensional Feature Space" by Fan and Lv (2008), specifically focusing on Setting 3 (S3) with logistic regression.

## Overview

The simulation study compares different variable screening methods for ultrahigh-dimensional logistic regression problems. The main methods evaluated are:

1. SIS with Marginal Likelihood Ratio (MLR) screening
2. SIS with Marginal Maximum Likelihood Estimation (MMLE) screening
3. MMLE ranking
4. LASSO with cross-validation
5. SCAD with cross-validation

## Simulation Setting (S3)

- Sample size (n): 600
- Number of predictors (p): 2000
- Number of simulations: 200
- True model sizes (s): 3, 6, 12, 24
- Screening size (nsis): floor(n/log(n))

The simulation setting S3 involves:
- First p-50 variables are independent standard normals
- Last 50 variables are generated based on the first s variables with alternating signs
- True coefficients are placed in positions p-50+1 to p-50+s
- Response is generated using logistic regression

## Dependencies

The code requires the following R packages:
- MASS
- glmnet
- ncvreg
- SIS

## Usage

1. Install the required R packages:
```R
install.packages(c("MASS", "glmnet", "ncvreg", "SIS"))
```

2. Run the simulation:
```R
source("simulation_code.R")
```

The code will:
- Generate data according to Setting S3
- Apply all screening methods
- Calculate the minimum model size (MMS) for each method
- Output results in a table format similar to Table 5 in the paper

## Output

The simulation results are presented in a table showing:
- The median MMS for each method
- The associated robust standard deviation (RSD)
- The maximum eigenvalue (M-lambda_max) and its RSD

## Notes

- The code includes error handling for various edge cases
- Results are reproducible with the set seed (123456789)
- The simulation may take some time to complete due to the number of iterations and methods being evaluated

## References

Fan, J., & Lv, J. (2008). Sure independence screening for ultrahigh dimensional feature space. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70(5), 849-911.