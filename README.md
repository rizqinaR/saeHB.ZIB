
<!-- README.md is generated from README.Rmd. Please edit that file -->

# saeHB.ZIB

<!-- badges: start -->
<!-- badges: end -->

Provides function for area level of small area estimation using
hierarchical Bayesian (HB) method with Zero-Inflated Binomial
distribution for variables of interest. Some dataset produced by a data
generation are also provided. The ‘rjags’ package is employed to obtain
parameter estimates. Model-based estimators involves the HB estimators
which include the mean, the variation of mean, and the quantile.

## Author

Rizqina Rahmati , Azka Ubaidillah

## Maintaner

Rizqina Rahmati <221810583@stis.ac.id>

## Installation

You can install the development version of saeHB.ZIB from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rizqinaR/saeHB.ZIB")
```

## Function

-   `ziBinomial()` This function gives small area estimator under Zero
    Inflated Binomial Model and is implemented to variable of
    interest (y) that assumed to be a Zero Inflated Binomial
    Distribution. The range of data is (y >= 0)

## Example

This is a basic example of using `ziBinomial()` function to make an
estimate based on synthetic data in this package

``` r
library(saeHB.ZIB)
## For data without any non-sampled area
data(dataZIB)       # Load dataset
## For data with non-sampled area use dataZIB.ns
## Fitting model
result <- ziBinomial(y ~ x1 + x2, n.samp = "n.samp", data=dataZIB)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 98
#>    Total graph size: 632
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 98
#>    Total graph size: 632
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 98
#>    Total graph size: 632
#> 
#> Initializing model
```

Small Area mean Estimates

``` r
result$Est
```

Estimated model coefficient

``` r
result$coefficient
```

Estimated random effect variances

``` r
result$refVar
```

## References

*W. Bodromurti, “ZERO INFLATED BINOMIAL MODELS IN SMALL AREA ESTIMATION WITH APPLICATION TO INFANT MORTALITY DATA IN INDONESIA”, in Thesis IPB University Scientific Repository, 2017. 
*Matranga, D, Firenze, A, Vullo, A. Can bayesian models play a role in dental caries epidemiology? Evidence from an application to the BELCAP data set. Community Dent Oral Epidemiol 2013; 41: 473– 480. © 2013 John Wiley & Sons A/S. Published by John Wiley & Sons Ltd. 
*B. Hartono, “Kajian Pendugaan Area Kecil dengan Model Zero Inflated Binomial (Studi Kasus Pengangguran di Provinsi Jambi)”, in Thesis IPB University Scientific Repository, 2018. 
*Hall DB. Zero-inflated Poisson and binomial regression with random effects: a case study. Biometrics 2000;56(4):1030–9.
