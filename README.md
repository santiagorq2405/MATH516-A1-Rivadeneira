# A1: Bi-Log-Normal Mixture for Snow-Flake Diameters

## MATH 516 - Applied Statistics, EPFL

### Author
Santiago Rivadeneira Quintero

### Description
Fitting a mixture of two log-normal distributions to binned snow-flake diameter data using maximum likelihood (EM algorithm and direct optimization) and Bayesian inference, with parametric bootstrap goodness-of-fit testing.

### Requirements
- R >= 4.4
- Quarto >= 1.4

### R Packages
```r
install.packages(c("ggplot2", "dplyr", "tidyr",
                   "gridExtra", "scales", "knitr", "kableExtra"))
```

### How to Reproduce
1. Clone the repository
2. Install the required R packages listed above
3. Render the report:
```bash
cd A1
quarto render report.qmd
```
4. The output `report.pdf` will be generated in the same directory

### Project Structure
```
A1/
  report.qmd              Quarto report with all analysis code
  report.pdf              Final PDF report (max 15 pages)
  README.md               This file
  src/
    utils.R               Helper functions (densities, likelihoods, EM, bootstrap)
  plots/                  Generated figures
  1_snow_particles.csv    Raw data
```

### Data
Snow-flake diameter measurements from the Laboratory of Cryospheric Sciences at EPFL. Data consists of 52 bins with particle counts (N = 705,044 total particles). Shared with permission of the authors of Melo et al. (2021).
