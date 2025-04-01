
# standart <a href='https://thierrygosselin.github.io/standart/'><img src='man/figures/logo.png' align="right" height="200" /></a>

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Last-changedate](https://img.shields.io/badge/last%20change-2025--04--01-brightgreen.svg)](/commits/master)
<!-- badges: end -->

standart provides R functions that help you solve the most common
quality control challenges of genome reduction sequences (DArT, RADseq,
GBS, etc).

standart really shines in projects that requires consistency among
sequencing runs and when samples are from animals with large and
complicated genomes.

### Who is it for ?

The package is currently developed for projects at CSIRO.

## Installation

You can install the development version of standart like so:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("thierrygosselin/standart")
library(standart)
```

## Overview

- Reads/sequences counter and length distribution.
- Read depth plot that highlights read coverage groups. Distinct and
  unique reads with high coverage are repetitive elements that when
  assembled in locus are usually paralogs, retrotransposons,
  transposable elements, etc.
- Noise reduction and fastq cleaning.
- Data normalization/rarefaction to remove or reduce ascertainment bias
  driving populations polymorphism discovery bias, missing data
  patterns, individual heterozygosity problems or patterns.
- Noise reduction and normalization will help to answer the question:
  are the differences between your samples: artifactual or biological.

## Getting help

- [Vignettes](https://thierrygosselin.github.io/standart/articles/index.html)
- [Computer setup and
  troubleshooting](https://thierrygosselin.github.io/radiator/articles/rad_genomics_computer_setup.html)
- My [github](https://github.com/thierrygosselin) as other packages to
  help in RAD/DArTseq analysis

[<img src="man/figures/stackr_logo.png" width="100" alt="stackr" />](https://thierrygosselin.github.io/stackr/)
[<img src="man/figures/radr_logo.png" width="100" alt="radr" />](https://thierrygosselin.github.io/radiator/)
[<img src="man/figures/assigner_logo.png" width="100" alt="assigner" />](https://thierrygosselin.github.io/assigner/)
[<img src="man/figures/grur_logo.png" width="100" alt="grur" />](https://thierrygosselin.github.io/grur/)
