#' @title Sample Data for Small Area Estimation using Hierarchical Bayesian Method under Zero-Inflated Binomial Distribution
#' @description Dataset to simulate Small Area Estimation using Hierarchical Bayesian Method under Zero-Inflated Binomial distribution with non-sampled areas
#'
#' This data contains NA values that indicates no sampled at one or more small areas. It uses the dataZIB.ns with the direct estimates and the related variances in 3 small areas are missing.
#' @usage data(dataZIB.ns)
#'
#' @format A data frame with 30 rows and 4 variables :
#' \describe{
#' \item{y}{Direct Estimation of y}
#' \item{x1}{Auxiliary variable of x1}
#' \item{x2}{Auxiliary variable of x2}
#' \item{vardir}{sampling variance of y}
#' \item{n.samp}{number of sample}
#' }
#'
#'
"dataZIB.ns"
