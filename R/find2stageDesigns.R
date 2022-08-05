#' Find two-stage designs
#'
#' This function finds two-stage designs for a given set of design parameters, allowing
#' stopping for benefit at the interim (Mander and Thompson's design) or no stopping
#' for benefit at the interim (Simon's design). It returns not
#' only the optimal and minimax design realisations, but all design realisations that could
#' be considered "best" in terms of expected sample size under p=p0 (EssH0), expected
#' sample size under p=p1 (Ess), maximum sample size (n) or any weighted combination of these
#' three optimality criteria.
#'
#' @param nmin Minimum permitted sample size. Should be a multiple of block size or number of stages.
#' @param nmax Maximum permitted sample size. Should be a multiple of block size or number of stages.
#' @param p0 Probability for which to control the type-I error-rate
#' @param p1 Probability for which to control the power
#' @param alpha Significance level
#' @param power Required power (1-beta)
#' @param maxthetaF Maximum permitted conditional power for futility stopping. Optional.
#' @param benefit Allow the trial to end for a go decision and reject the null hypothesis at the interim analysis (i.e., the design of Mander and Thompson)
#' @return A list of class "curtailment_simon" containing two data frames. The first data frame, $input,
#' has a single row and contains all the inputted values. The second data frame, $all.des, contains one
#' row for each design realisation, and contains the details of each design, including sample size,
#' stopping boundaries and operating characteristics. To see a diagram of any obtained design realisation
#' and its corresponding stopping boundaries, simply call the function drawDiagram with this output as the only argument.
#' @author Martin Law, \email{martin.law@@mrc-bsu.cam.ac.uk}
#' @examples
#' \donttest{
#' find2stageDesigns(nmin=23,
#'  nmax=27,
#'  p0=0.75,
#'  p1=0.92,
#'  alpha=0.22,
#'  power=0.95,
#'  benefit=TRUE)
#'  }
#' @references
#' \doi{10.1016/j.cct.2010.07.008}{A.P. Mander, S.G. Thompson,
#' Two-stage designs optimal under the alternative hypothesis for phase II cancer clinical trials,
#' Contemporary Clinical Trials,
#' Volume 31, Issue 6,
#' 2010,
#' Pages 572-578}
#'
#' \doi{10.1016/0197-2456(89)90015-9}{Richard Simon,
#' Optimal two-stage designs for phase II clinical trials,
#' Controlled Clinical Trials,
#' Volume 10, Issue 1,
#' 1989,
#' Pages 1-10}
#' @export
find2stageDesigns <- function(nmin, nmax, p0, p1, alpha, power, maxthetaF=NA, benefit=FALSE)
{

  if(benefit==FALSE){
    nr.lists <- findSimonN1N2R1R2(nmin=nmin, nmax=nmax, e1=FALSE)
    simon.df <- apply(nr.lists, 1, function(x) {findSingleSimonDesign(n1=x[1], n2=x[2], r1=x[4], r=x[5], p0=p0, p1=p1)})
    simon.df <- t(simon.df)
    simon.df <- as.data.frame(simon.df)
    names(simon.df) <- c("n1", "n2", "n", "r1", "r", "alpha", "power", "EssH0", "Ess")
  } else{
    nr.lists <- findSimonN1N2R1R2(nmin=nmin, nmax=nmax, e1=TRUE)
    simon.df <- apply(nr.lists, 1, function(x) {simonEfficacy(n1=x[1], n2=x[2], r1=x[4], r=x[5], p0=p0, p1=p1, e1=x[6])})
    simon.df <- t(simon.df)
    simon.df <- as.data.frame(simon.df)
    names(simon.df) <- c("n1", "n2", "n", "r1", "r", "alpha", "power", "EssH0", "Ess", "e1")
  }
  correct.alpha.power <- simon.df$alpha < alpha & simon.df$power > power
  simon.df <- simon.df[correct.alpha.power, ]

# Find max CP among all terminal points for futility:
thetaF <- vector("numeric", nrow(simon.df))
for(i in 1:nrow(simon.df)){
  thetaF[i] <- find2stageThetaF(r=simon.df$r[i], r1=simon.df$r1[i], n2=simon.df$n2[i], p=p1)
}
simon.df$thetaF <- thetaF

if(!is.na(maxthetaF)){
  simon.df <- simon.df[simon.df$thetaF <= maxthetaF, ]
}




  if(nrow(simon.df)==0){
    stop("No suitable designs exist for these design parameters.")
  }
  # Discard all dominated designs. Note strict inequalities as EssH0 and Ess will (almost?) never be equal for two designs:
  discard <- rep(NA, nrow(simon.df))
  for(i in 1:nrow(simon.df))
  {
    discard[i] <- sum(simon.df$EssH0[i] > simon.df$EssH0 & simon.df$Ess[i] > simon.df$Ess & simon.df$n[i] >= simon.df$n)
  }

  simon.df <- simon.df[discard==0,]
  simon.input <- data.frame(nmin=nmin, nmax=nmax, p0=p0, p1=p1, alpha=alpha, power=power, maxthetaF=maxthetaF)
  simon.output <- list(input=simon.input,
                       all.des=simon.df)
  class(simon.output) <- append(class(simon.output), "curtailment_simon")
  return(simon.output)
}
