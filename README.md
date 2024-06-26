
# GOeval

<!-- badges: start -->
[![R-CMD-check](https://github.com/BrentLab/GOeval/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BrentLab/GOeval/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/BrentLab/GOeval/branch/main/graph/badge.svg)](https://app.codecov.io/gh/BrentLab/GOeval?branch=main)
<!-- badges: end -->

Welcome to GOeval - a package for gene regulatory network evaluation using Gene Ontology and other gene set databases.

This package provides functions to systematically run over-representation analysis (ORA) on multiple subsets of a gene regulatory network and its permutations using the 'WebGestaltR' package. It also provides functions to assess the quality of biological information captured by the network via metrics calculated from the ORA results.

## Installation

You can install GOeval with:

``` r
devtools::install_github("BrentLab/GOeval")
```

You can then load the package with:

``` r
library("GOeval")
```

## Documentation

https://BrentLab.github.io/GOeval/
