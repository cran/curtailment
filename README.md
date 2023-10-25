
<!-- README.md is generated from README.Rmd. Please edit that file -->

# curtailment

<!-- badges: start -->
<!-- badges: end -->

The goal of curtailment is to allow users to create single- or two-arm
binary outcome clinical trial designs that use a form of early stopping
known as stochastic curtailment. Using stochastic curtailment can result
in trial designs with a lower sample size on average, compared to
typical trial designs or even multi-stage designs. The package contains
functions that will search for designs that are suitable for your trial,
and draw a diagram showing all points in the trial that a decision may
be made to stop early. The package can also find two-stage designs that
stop for futility only (Simon designs) or for benefit or futility
(Mander and Thompson designs). For all designs, a maximum conditional
power for futility can be specified, above which a trial would not be
permitted to stop early for futility.

## Installation

You can install the released version of curtailment from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("curtailment")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("martinlaw/curtailment")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(curtailment)
## Find a single-arm design with a maximum sample size in the range 10-25, allowing early stopping after every 5 results, with a type-I error-rate of 0.05 when the response rate equals 0.1 and a power of 0.8 when the response rate equals 0.4:
output <- singlearmDesign(nmin = 10,
                          nmax = 25,
                          C = 5,
                          p0 = 0.1,
                          p1 = 0.4,
                          power = 0.8,
                          alpha = 0.05)
# Obtain the stopping boundaries and a diagram of the above trial:
fig <- drawDiagram(output)
```

## Future plans

Speed up Mander and Thompson and Simon design searches by reducing the
number of combinations of n1/n2/r1/r/e1 search over, using Wald’s
sequential probability ratio test as in the multi-stage/continuous
functions:

<!-- denom <- log(p1/p0) - log((1-p1)/(1-p0)) -->
<!-- accept.null <- log((1-power)/(1-alpha)) / denom  + nposs * log((1-p0)/(1-p1))/denom -->
<!-- accept.null <- floor(accept.null) -->
<!-- reject.null <- log((power)/alpha) / denom  + nposs * log((1-p0)/(1-p1))/denom -->
<!-- reject.null <- ceiling(reject.null) -->
<!-- r.wald <- NULL -->
<!-- ns.wald <- NULL -->
<!-- for(i in 1:length(nposs)){ -->
<!-- r.wald <- c(r.wald, accept.null[i]:reject.null[i]) -->
<!-- ns.wald <- c(ns.wald, rep(nposs[i], length(accept.null[i]:reject.null[i]))) -->
<!-- } -->
<!-- sc.subset <- data.frame(n=ns.wald, r=r.wald) -->
