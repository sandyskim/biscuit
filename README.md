<!-- README.md is generated from README.Rmd. Please edit that file -->

# biscuit

<!-- badges: start -->

<!-- badges: end -->

<p align="left">

<img src="man/figures/biscuit_logo.png" width="150"/>

</p>

`biscuit` is an R package for analyzing CRISPR screens utilizing a Bayesian hierarchical inference model.

## Installation

`biscuit` requires `cmdstanr` to perform inference. To install `cmdstanr` and `cmdstan`,

``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

library(cmdstanr)

install_cmdstan(cores = 2)
```

To install `biscuit`,

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("sandyskim/biscuit")
```

## Usage

We provide a full tutorial that demonstrates how to use `biscuit` with an example data set in the introductory [vignette](https://github.com/sandyskim/biscuit/blob/main/vignettes/introduction.Rmd/).
