% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findSCdesigns.R
\name{findSCdesigns}
\alias{findSCdesigns}
\title{findSCDesigns}
\usage{
findSCdesigns(
  nmin,
  nmax,
  p0,
  p1,
  alpha,
  power,
  minthetaE = p1,
  maxthetaF = p1,
  bounds = "wald",
  fixed.r1 = NA,
  fixed.r = NA,
  fixed.n1 = NA,
  max.combns = 1e+06,
  maxthetas = NA,
  exact.thetaF = NA,
  exact.thetaE = NA,
  progressBar = FALSE
)
}
\arguments{
\item{nmin}{Minimum permitted sample size.}

\item{nmax}{Maximum permitted sample size.}

\item{p0}{Probability for which to control the type-I error-rate}

\item{p1}{Probability for which to control the power}

\item{alpha}{Significance level}

\item{power}{Required power (1-beta).}

\item{minthetaE}{Minimum value of upper threshold theta_E_min. Defaults to p.}

\item{maxthetaF}{Maximum value of lower CP threshold theta_F_max. Defaults to p.}

\item{bounds}{choose what final rejection boundaries should be searched over: Those of A'Hern ("ahern"), Wald ("wald") or no constraints (NA). Defaults to "wald".}

\item{fixed.r1}{Choose what interim rejection boundaries should be searched over. Useful for reproducing a particular design realisation. Defaults to NA.}

\item{fixed.r}{Choose what final rejection boundaries should be searched over. Useful for reproducing a particular design realisation. Defaults to NA.}

\item{fixed.n1}{Choose what interim sample size values n1 should be searched over. Useful for reproducing a particular design realisation. Defaults to NA.}

\item{max.combns}{Provide a maximum number of ordered pairs (theta_F, theta_E). Defaults to 1e6.}

\item{maxthetas}{Provide a maximum number of CP values used to create ordered pairs (theta_F, theta_E). Can be used instead of max.combns. Defaults to NA.}

\item{exact.thetaF}{Provide an exact value for lower threshold theta_F. Useful for reproducing a particular design realisation. Defaults to NA.}

\item{exact.thetaE}{Provide an exact value for upper threshold theta_E. Useful for reproducing a particular design realisation. Defaults to NA.}

\item{progressBar}{Logical. If TRUE, shows progress bar. Defaults to FALSE.}
}
\value{
Output is a list of two dataframes. The first, $input, is a one-row
data frame that contains important arguments used in the call. The second,
$all.des,contains the operating characteristics of all admissible designs found.
}
\description{
This function finds admissible design realisations for single-arm binary outcome trials, using stochastic curtailment.
This function differs from singlearmDesign in that it includes a Simon-style interim analysis after some n1 participants.
The output is a data frame of admissible design realisations.
}
\examples{
\donttest{findSCdesigns(nmin = 20, nmax = 20, p0 = 0.1, p1 = 0.4, power = 0.8, alpha = 0.1, max.combns=1e2)}
}
