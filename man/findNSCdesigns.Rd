% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NSC_functions.R
\name{findNSCdesigns}
\alias{findNSCdesigns}
\title{findNSCdesigns}
\usage{
findNSCdesigns(nmin, nmax, p0, p1, alpha, power, progressBar = FALSE)
}
\arguments{
\item{nmin}{Minimum permitted sample size.}

\item{nmax}{Maximum permitted sample size.}

\item{p0}{Probability for which to control the type-I error-rate}

\item{p1}{Probability for which to control the power}

\item{alpha}{Significance level}

\item{power}{Required power (1-beta).}

\item{progressBar}{Logical. If TRUE, shows progress bar. Defaults to FALSE.}
}
\value{
Output is a list of two dataframes. The first, $input, is a one-row
data frame that contains important arguments used in the call. The second,
$all.des,contains the operating characteristics of all admissible designs found.
}
\description{
This function finds admissible design realisations for single-arm binary outcome trials, using non-stochastic curtailment.
The output is a data frame of admissible design realisations.
}
\examples{
findNSCdesigns(nmin=20, nmax=21, p0=0.1, p1=0.4, alpha=0.1, power=0.8)
}
