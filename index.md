<!-- badges: start -->
[![R-CMD-check](https://github.com/westbrooktm/GOeval/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/westbrooktm/GOeval/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/westbrooktm/GOeval/branch/main/graph/badge.svg)](https://app.codecov.io/gh/westbrooktm/GOeval?branch=main)
<!-- badges: end -->

Welcome to GOeval - a package for gene regulatory network evaluation using Gene Ontology and other gene set databases.

This package provides functions to systematically run over-representation analysis (ORA) on multiple subsets of a gene regulatory network and its permutations using the 'WebGestaltR' package. It also provides functions to assess the quality of biological information captured by the network via metrics calculated from the ORA results.

Here, you will find documentation on the purposes of the package's functions and their arguments, as well as a "Get started" vignette showing example usage of the package's functions. Link to vignette: https://westbrooktm.github.io/GOeval/articles/GOeval.html.

Input format: You will need a network file in which each row represents a directed edge. The tab-separated entries in each row are the the source node (e.g. transcription factor), the target node (e.g. the regulated gene), and optionally an edge score. The subset_network function requires edge scores. The webgestalt_network function requires a “gene universe” file that has a single column containing all the genes that could possibly be present in the network.
