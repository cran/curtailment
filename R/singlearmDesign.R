#' Find single-arm trials using stochastic curtailment
#'
#' This function finds admissible design realisations for single-arm binary outcome trials, using stochastic curtailment.
#' The output can be used as the sole argument in the function 'drawDiagram', which will return the stopping boundaries for the
#' admissible design of your choice. Monitoring frequency can set in terms of block(/cohort) size ("C") or number of stages ("stages").
#' @param nmin Smallest maximum sample size. Should be a multiple of block size or number of stages.
#' @param nmax Largest maximum sample size. Should be a multiple of block size or number of stages.
#' @param C Block size. Vectors, i.e., multiple values, are permitted.
#' @param stages Number of interim analyses or "stages". Only required if not setting block size C. Vectors, i.e., multiple values, are permitted.
#' @param p0 Probability for which to control the type-I error-rate
#' @param p1 Probability for which to control the power
#' @param alpha Significance level
#' @param power Required power (1-beta).
#' @param minstop Minimum permitted sample size at the first interim analysis
#' @param maxthetaF Maximum value of lower CP threshold theta_F_max. Defaults to p1.
#' @param minthetaE Minimum value of upper threshold theta_E_min. Defaults to p1.
#' @param bounds choose what final rejection boundaries should be searched over: Those of A'Hern ("ahern"), Wald ("wald") or no constraints (NA). Defaults to "wald".
#' @param return.only.admissible Logical. Returns only admissible design realisations if TRUE, otherwise returns all feasible designs. Defaults to TRUE.
#' @param max.combns Provide a maximum number of ordered pairs (theta_F, theta_E). Defaults to 1e6.
#' @param maxthetas Provide a maximum number of CP values used to create ordered pairs (theta_F, theta_E). Can be used instead of max.combns. Defaults to NA.
#' @param fixed.r Choose what final rejection boundaries should be searched over. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param exact.thetaF Provide an exact value for lower threshold theta_F. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param exact.thetaE Provide an exact value for upper threshold theta_E. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param progressBar Logical. If TRUE, shows progress bar. Defaults to FALSE.
#' @export
#' @author Martin Law, \email{martin.law@@mrc-bsu.cam.ac.uk}
#' @return Output is a list of two dataframes. The first, $input, is a one-row data frame that contains all the arguments used in the call. The second, $all.des, contains the operating characteristics of all admissible designs found.
#' @examples output <- singlearmDesign(nmin = 30,
#'  nmax = 30,
#'  C = 5,
#'  p0 = 0.1,
#'  p1 = 0.4,
#'  power = 0.8,
#'  alpha = 0.05)
singlearmDesign <- function(nmin,
                        nmax,
                        C=NA,
                        stages=NA,
                        p0,
                        p1,
                        alpha,
                        power,
                        minstop=1,
                        maxthetaF=p1,
                        minthetaE=p1,
                        bounds="wald",
                        return.only.admissible=TRUE,
                        max.combns=1e6,
                        maxthetas=NA,
                        fixed.r=NA,
                        exact.thetaF=NA,
                        exact.thetaE=NA,
                        progressBar=FALSE){
  use.stages <- any(is.na(C))
  if(any(!is.na(C)) & any(!is.na(stages))) stop("Values given for both cohort/block size C and number of stages. Please choose one only.")
  if(minstop>=nmin) stop("earliest stopping point (minstop) must be smaller than all potential maximum sample sizes, i.e. smaller than nmin.")

  if(use.stages==TRUE){
    intermediate.output <- lapply(stages,
                                  FUN = function(x) findDesignsGivenCohortStage(stages=x,
                                                         C=NA,
                                                         nmin=nmin,
                                                         nmax=nmax,
                                                         p0=p0,
                                                         p1=p1,
                                                         alpha=alpha,
                                                         power=power,
                                                         minstop=minstop,
                                                         maxthetaF=maxthetaF,
                                                         minthetaE=minthetaE,
                                                         bounds=bounds,
                                                         return.only.admissible=return.only.admissible,
                                                         max.combns=max.combns,
                                                         maxthetas=maxthetas,
                                                         fixed.r=fixed.r,
                                                         exact.thetaF=exact.thetaF,
                                                         exact.thetaE=exact.thetaE,
                                                         progressBar=progressBar)
           )
  }else{
    intermediate.output <- lapply(C,
                                  FUN = function(x) findDesignsGivenCohortStage(stages=NA,
                                                         C=x,
                                                         nmin=nmin,
                                                         nmax=nmax,
                                                         p0=p0,
                                                         p1=p1,
                                                         alpha=alpha,
                                                         power=power,
                                                         minstop=minstop,
                                                         maxthetaF=maxthetaF,
                                                         minthetaE=minthetaE,
                                                         bounds=bounds,
                                                         return.only.admissible=return.only.admissible,
                                                         max.combns=max.combns,
                                                         maxthetas=maxthetas,
                                                         fixed.r=fixed.r,
                                                         exact.thetaF=exact.thetaF,
                                                         exact.thetaE=exact.thetaE,
                                                         progressBar=progressBar)
    )
  }
  input <- data.frame(nmin=nmin, nmax=nmax,
                      Cmin=C[1], Cmax=C[length(C)],
                      stagesmin=stages[1], stagesmax=stages[length(stages)],
                      minstop=minstop,
                      p0=p0, p1=p1, alpha=alpha, power=power,
                      maxthetaF=maxthetaF, minthetaE=minthetaE, bounds=bounds, fixed.r=fixed.r,
                      return.only.admissible=return.only.admissible, max.combns=max.combns,
                      maxthetas=maxthetas, exact.thetaF=exact.thetaF, exact.thetaE=exact.thetaE)
  intermediate.output.combined <- do.call(rbind, intermediate.output)
  final.output <- list(input=input,
                       all.des=intermediate.output.combined)
  class(final.output) <- append(class(final.output), "curtailment_single")
  return(final.output)
}
