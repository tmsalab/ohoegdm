# ohoegdm 0.1.1

## Changes

- Added explicit dependencies on R (>= 4.3.0), Rcpp (>= 1.1.0), and RcppArmadillo (>= 15.0.2-2)
- Removed CXX11 from `src/Makevars` and `src/Makevars.win` to avoid potential compilation issues
  with newer versions of Armadillo through RcppArmadillo.
- Switched README.Rmd to README.qmd to use Quarto for rendering.
- Fixed CITATION file to use `bibentry()` instead of `citEntry()` to 
  avoid CRAN check notes.
- Updated GitHub Action workflows.


# ohoegdm 0.1.0

## New Features

- Enabled the Ordinal Higher-order Exploratory General Diagnostic Model.

## Deployment

- Added a GitHub Action's deployment.
