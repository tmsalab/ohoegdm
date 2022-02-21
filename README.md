
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ohoegdm

<!-- badges: start -->

[![R-CMD-check](https://github.com/tmsalab/ohoegdm/workflows/R-CMD-check/badge.svg)](https://github.com/tmsalab/ohoegdm/actions)
[![Package-License](https://img.shields.io/badge/license-GPL%20(%3E=2)-brightgreen.svg?style=flat)](https://www.gnu.org/licenses/gpl-2.0.html)
<!-- badges: end -->

The goal of `ohoegdm` is to provide an implementation of the Ordinal
Higher-order Exploratory General Diagnostic Model for Polytomous Data as
described by Culpepper and Balamuta (In Press).

## Installation

You can install the released version of ohoegdm from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ohoegdm")
```

Or, you can be on the cutting-edge development version on
[GitHub](https://github.com/) using:

``` r
# install.packages("devtools")
devtools::install_github("tmsalab/ohoegdm")
```

## Usage

To use `ohoegdm`, load the package using:

``` r
library("ohoegdm")
```

From here, the OHO-EGDM model can be estimated using:

``` r
my_model = ohoegdm::ohoegdm(
  y = <data>,
  k = <k>,
  m = <item-responses-categories>,
  order = <model-interaction-order>
)
```

## Authors

Steven Andrew Culpepper and James Joseph Balamuta

## Citing the `ohoegdm` package

To ensure future development of the package, please cite `ohoegdm`
package if used during an analysis or simulation study. Citation
information for the package may be acquired by using in *R*:

``` r
citation("ohoegdm")
```

## License

GPL (>= 2)
