
# standart <a href='https://thierrygosselin.github.io/standart/'><img src='man/figures/logo.png' align="right" height="200" /></a>

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Last-changedate](https://img.shields.io/badge/last%20change-2025--04--02-brightgreen.svg)](/commits/master)
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

- **Detect**
  - `read_counter`: Count reads/sequences and generate distributions by
    stratifications or not (think populations, sampling sites,
    sequencing runs/technologies)
  - `read_depth_plot`: Generate a figure of read coverage group to
    highlight the different categories of read depth. Reads and
    sequences are not created equals, some are distinct, unique and with
    high coverage. These repetitive elements can potentially be
    problematic in downstream analysis (think paralogs,
    retrotransposons, transposable elements).
- **Selection & Correction**
  - `clean_fq`: Noisy fastq files identified with the functions above
    can be cleaned. Different options available to tidy those files and
    have good bookkeeping habits and reporting.
  - `normalize_reads`: Normalize/Rarefy the fastq sequences. This
    technique can be used for a number of reasons, e.g. to remove or
    reduce ascertainment bias, discovery bias driving populations
    polymorphism, missing data patterns, individual heterozygosity
    problems or patterns, etc.
- **Artifactual or Biological signals:** Not certain if the differences
  between your samples or your sequencing runs is artifactual or
  biological ? Both noise reduction and normalization can help answer
  the question.

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
