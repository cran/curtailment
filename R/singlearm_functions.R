#' @import data.table
#' @import ggplot2
#' @import ggthemes
#' @import pkgcond
#' @import stats
#' @importFrom grDevices dev.off pdf
#' @importFrom utils combn head setTxtProgressBar txtProgressBar
#requireNamespace("data.table")

#### EXPLANATION OF CODE
# singlearmDesign: used to search for m-stage designs
# findDesignOCs: called within singlearmDesign
# findCPmatrix: called within singlearmDesign
# findSingleSimonDesign: finds operating characteristics of a single Simon design
# findWald: finds ESS under H0 and H1 for Wald's design



##################### FUNCTION findCPmatrix
# n = sample size
# r = stopping boundary
# Csize = Block size

findCPmatrix <- function(n, r, Csize, p0, p1){
  q1 <- 1-p1
  mat <- matrix(3, ncol=n, nrow=max(3, min(r+Csize, n)+1)) #  nrow=r+Csize+1 unless number of stages equals 1, minimum of 3.
  rownames(mat) <- 0:(nrow(mat)-1)
  mat[(r+2):nrow(mat),] <- 1
  mat[1:(r+1),n] <- 0

  pat.cols <- seq(n, 1, by=-Csize)[-1]

  for(i in (r+1):1){
    for(j in pat.cols){  # Only look every C patients (no need to look at final col)
      if((r+1-(i-1) > n-j)) { # IF success is not possible [Total responses needed (r+1) - responses so far (i-1)] > [no. of patients remaining (n-j)], THEN
        mat[i,j] <- 0
      }else{
        if(i-1<=j){ # Condition: Sm<=m
          newcp <- sum(choose(Csize,0:Csize)*p1^(0:Csize)*q1^(Csize:0)*mat[i:(i+Csize),j+Csize])
          mat[i,j] <- newcp
        }
      }
    }
  }


for(i in 3:nrow(mat)) {
    mat[i, 1:(i-2)] <- NA
  }
mat
}


##################### FUNCTION findWald

findWald <- function(alpha, power, p0, p1){
  beta <- 1-power
  h0.top <- (1-alpha)*log(beta/(1-alpha))+alpha*log((1-beta)/alpha)
  h0.bottom <- p0*log(p1/p0)+(1-p0)*log((1-p1)/(1-p0))
  Ep0n <- h0.top/h0.bottom
  h1.top <- beta*log(beta/(1-alpha))+(1-beta)*log((1-beta)/alpha)
  h1.bottom <- p1*log(p1/p0)+(1-p1)*log((1-p1)/(1-p0))
  Ep1n <- h1.top/h1.bottom
  round(c(Ep0n, Ep1n),1)
}

##################### FUNCTION findSingleSimonDesign
# Simon's design: No curtailment -- only stopping is at end of S1:

#' Find a single Simon design
#'
#' This function finds the operating characteristics of a single realisation of a Simon design. It returns
#' the inputted values along with the type-I error-rate (alpha), power,  expected sample size under p=p0 (EssH0)
#' and expected sample size under p=p1 (Ess).
#'
#' @param n1 Number of participants in stage 1
#' @param n2 Number of participants in stage 2
#' @param r1 Interim stopping boundary that must be exceeded to avoid no go stopping
#' @param r Final rejection boundary that must be exceeded to reject the null hypothesis.
#' @param p0 Anticipated response probability for inefficacious treatment
#' @param p1 Anticipated response probability for efficacious treatment
#' @return A vector containing all the inputted values and corresponding operating characteristics.
#' @author Martin Law, \email{martin.law@@mrc-bsu.cam.ac.uk}
#' @references
#' \doi{10.1016/0197-2456(89)90015-9}{Richard Simon,
#' Optimal two-stage designs for phase II clinical trials,
#' Controlled Clinical Trials,
#' Volume 10, Issue 1,
#' 1989,
#' Pages 1-10}
#' @examples findSingleSimonDesign(n1 = 15, n2 = 11, r1 = 1, r = 4, p0 = 0.1, p1 = 0.3)
#' @export
findSingleSimonDesign <- function(n1, n2, r1, r, p0, p1)
  {
  n <- n1+n2
  # Create Pascal's triangle for S1: these are the coefficients (before curtailment) A, where A * p^b * q*c
  pascal.list.s1 <- list(1)
  for (i in 2:(n1+1)) pascal.list.s1[[i]] <- c(0, pascal.list.s1[[i-1]]) + c(pascal.list.s1[[i-1]], 0)
  pascal.list.s1[[1]] <- NULL
  # For Simon's design, only need the final line:
  pascal.list.s1 <- pascal.list.s1[n1]

  # Curtail at n1 only:
  curtail.index <- 1:(r1+1)
  curtail.coefs.s1 <- pascal.list.s1[[1]][curtail.index] # from k=0 to k=r1

  # Use final column from S1:
  pascal.list.s2 <- pascal.list.s1
  pascal.list.s2[[1]][curtail.index] <- 0

  for (i in 2:(n2+1)) pascal.list.s2[[i]] <- c(0, pascal.list.s2[[i-1]]) + c(pascal.list.s2[[i-1]], 0)
  pascal.list.s2[[1]] <- NULL

  # Now obtain the rest of the probability -- the p^b * q^c :
  # S1
  q1 <- 1-p1
  coeffs <- p1^(0:n1)*q1^(n1:0)
  coeffs <- coeffs[curtail.index]

  q0 <- 1-p0
  coeffs.p0 <- p0^(0:n1)*q0^(n1:0)
  coeffs.p0 <- coeffs.p0[curtail.index]

  # Multiply the two vectors (A and p^b * q^c):
  prob.curt.s1 <- curtail.coefs.s1*coeffs

  # for finding type I error prob:
  prob.curt.s1.p0 <- curtail.coefs.s1*coeffs.p0

  # The (S1) curtailed paths:
  k.curt.s1 <- 0:r1
  n.curt.s1 <- rep(n1, length(k.curt.s1))
  curtail.s1 <- cbind(k.curt.s1, n.curt.s1, prob.curt.s1, prob.curt.s1.p0)


  ############## S2

  # Pick out the coefficients for the S2 paths (A, say):
  s2.index <- (r1+2):(n+1)
  curtail.coefs.s2 <- pascal.list.s2[[n2]][s2.index]

  # Now obtain the rest of the probability -- the p^b * q^c :
  coeffs.s2 <- p1^(0:n)*q1^(n:0)
  coeffs.s2 <- coeffs.s2[s2.index]

  coeffs.s2.p0 <- p0^(0:n)*q0^(n:0)
  coeffs.s2.p0 <- coeffs.s2.p0[s2.index]

  # Multiply the two vectors (A and p^b * q^c):
  prob.go <- curtail.coefs.s2*coeffs.s2

  # for finding type I error prob:
  prob.go.p0 <- curtail.coefs.s2*coeffs.s2.p0

  # Paths that reach the end:
  k.go <- (r1+1):n
  n.go <- rep(n, length(k.go))

  go <- cbind(k.go, n.go, prob.go, prob.go.p0)

  final <- rbind(curtail.s1, go)

  ############## WRAPPING UP THE RESULTS

  output <- data.frame(k=final[,1], n=final[,2], prob=final[,3], prob.p0=final[,4])
  output$success <- "Fail"
  output$success[output$k > r] <- "Success"
  power <- sum(output$prob[output$success=="Success"])
  sample.size.expd <- sum(output$n*output$prob)
  sample.size.expd.p0 <- sum(output$n*output$prob.p0)
  alpha <- sum(output$prob.p0[output$success=="Success"])
  to.return <- c(n1=n1, n2=n2, n=n, r1=r1, r=r, alpha=alpha, power=power, EssH0=sample.size.expd.p0, Ess=sample.size.expd)
  to.return
}



##################### FUNCTION findDesignOCs
# n = sample size
# r = stopping boundary
# C = Block size
# thetaF = theta_F (lower threshold)
# thetaE = theta_E (upper threshold)
# mat = Conditional Power matrix (before application of stochastic curtailment)
# alpha = significance level
# power = required power (1-beta)
# coeffs.p0 = coefficients under H0
# coeffs = coefficients under H1

findDesignOCs <- function(n, r, C, thetaF, thetaE, mat, power, alpha, coeffs, coeffs.p0, p0, p1, return.tps=FALSE){

  q0 <- 1-p0
  q1 <- 1-p1

  interims <- seq(C, n, by=C)
  pat.cols <- rev(interims)[-1]

  # Amend CP matrix, rounding up to 1 when CP>theta_1 and rounding down to 0 when CP<thetaF:

  const <- choose(C,0:C)*p1^(0:C)*q1^(C:0)
  for(i in (r+1):1){
    for(j in pat.cols){  # Only look every C patients (no need to look at final col)
      if(i-1<=j){ # Condition: Sm<=m
        newcp <- sum(const*mat[i:(i+C), j+C])
        if(newcp > thetaE){
          mat[i,j] <- 1
          }else{
            if(newcp < thetaF){
              mat[i,j] <- 0
            }else{
              mat[i,j] <- newcp
            }
          }
      }
    }
  }

  # for(i in (r+1):1){
  #   for(j in pat.cols){  # Only look every C patients (no need to look at final col)
  #     if(i-1<=j){ # Condition: Sm<=m
  #       newcp <- sum(choose(C,0:C)*p1^(0:C)*q1^(C:0)*mat[i:(i+C),j+C])
  #       if(newcp > thetaE){
  #         mat[i,j] <- 1
  #       }else{
  #         if(newcp < thetaF){
  #           mat[i,j] <- 0
  #         }else{
  #         mat[i,j] <- newcp
  #       }
  #     }
  #   }
  #   }
  # }

 #if(abs(thetaF-0.352)<0.0001 & abs(thetaE-0.98823057)<0.0001) browser()
  ###### STOP if design is pointless, i.e either failure or success is not possible:

  # IF DESIGN GUARANTEES FAILURE (==0) or SUCCESS (==2) at n=C:
  first.cohort <- sum(mat[,C], na.rm = T)

  if(first.cohort==C+1){
    return(c(n, r, C, 1, 1,  NA, NA, thetaF, thetaE, NA))
  }

  if(first.cohort==0){
    return(c(n, r, C, 0, 0,  NA, NA, thetaF, thetaE, NA))
  }

  ############# Number of paths to each point:
  pascal.list <- list(c(1,1))
  if(C>1){
    for(i in 2:C){
      pascal.list[[i]] <- c(0, pascal.list[[i-1]]) + c(pascal.list[[i-1]], 0)
    }
  }

  for(i in min(C+1,n):n){ # Adding min() in case of stages=1, i.e. cohort size C=n.
    if(i %% C == 1 || C==1){
      column <- as.numeric(mat[!is.na(mat[,i-1]), i-1])
      newnew <- pascal.list[[i-1]]
      CPzero.or.one <- which(column==0 | column==1)
      newnew[CPzero.or.one] <- 0
      pascal.list[[i]] <- c(0, newnew) + c(newnew, 0)
    } else {
      pascal.list[[i]] <- c(0, pascal.list[[i-1]]) + c(pascal.list[[i-1]], 0)
    }
  }

  for(i in 1:length(pascal.list)){
    pascal.list[[i]] <- c(pascal.list[[i]], rep(0, length(pascal.list)+1-length(pascal.list[[i]])))
  }


  pascal.t <- t(do.call(rbind, args = pascal.list))
  pascal.t <- pascal.t[1:min(r+C+1, n+1), ] # Adding min() in case of stages=1, i.e. cohort size C=n.

  # Multiply the two matrices (A and p^b * q^c). This gives the probability of reaching each point:
  # Note: Only need the first r+C+1 rows -- no other rows can be reached.
  coeffs <- coeffs[1:min(r+C+1, n+1), 1:n]
  coeffs.p0 <- coeffs.p0[1:min(r+C+1, n+1), 1:n]
  pascal.t <- pascal.t[1:min(r+C+1, n+1), ]
  final.probs.mat <- pascal.t*coeffs
  final.probs.mat.p0 <- pascal.t*coeffs.p0


  ##### Now find the terminal points: m must be a multiple of C, and CP must be 1 or 0:
  ##### SUCCESS:
  ## Only loop over rows that contain CP=1:
  #rows.with.cp1 <- which(apply(mat, 1, function(x) {any(x==1, na.rm = T)}))
  rows.with.cp1 <- which(rowSums(mat==1, na.rm = T)>0)
  m.success <- NULL
  Sm.success <- NULL
  prob.success <- NULL
  prob.success.p0 <- NULL
  #

  for(i in rows.with.cp1){
    for(j in interims[interims >= i-1]){ # condition ensures Sm <= m
      if(mat[i,j]==1){
        m.success <- c(m.success, j)
        Sm.success <- c(Sm.success, i-1)
        prob.success <- c(prob.success, final.probs.mat[i, j])
        prob.success.p0 <- c(prob.success.p0, final.probs.mat.p0[i, j])
      }
    }
  }

  # mat.equals.1 <- mat==1
  # row.index.list <- alply(mat.equals.1[, interims], 2, which)
  # # alply always returns a list, whereas apply returns a matrix when all cols have same number of results.
  # row.index <- unlist(row.index.list)
  # total.length <- unlist(lapply(row.index.list, length))
  # m.success <- rep(interims, c(total.length))
  # index <- cbind(row.index, m.success)
   # for(i in 1:nrow(index)){
   #   prob.success[i] <- final.probs.mat[index[i, 1], index[i, 2]]
   #   prob.success.p0[i] <- final.probs.mat.p0[index[i, 1], index[i, 2]]
   #   }
   # prob.success <- apply(index, 1, function(x) final.probs.mat[x[1], x[2]])
  # prob.success.p0 <- apply(index, 1, function(x) final.probs.mat.p0[x[1], x[2]])

  success <- data.frame(Sm=Sm.success, m=m.success, prob=prob.success, success=rep("Success", length(m.success)))


  #
  # CAN FIND POWER AND TYPE I ERROR AT THIS POINT.
  # IF TOO LOW OR HIGH, STOP AND MOVE ON:
  #

  # pwr <- sum(success[,"prob"])
  # typeI <- sum(success[,"prob.p0"])
  pwr <- sum(prob.success)
  typeI <- sum(prob.success.p0)


  if(pwr < power | typeI > alpha) # | pwr > power+tol | alpha < alpha-tol)
  {
    return(c(n, r, C, typeI, pwr,  NA, NA, thetaF, thetaE, NA))
  }


  ##### FAILURE:

  #first.zero.in.row <- apply(mat[1:(r+1),], 1, function(x) {which.max(x[interims]==0)})
  submat <- mat[1:(r+1), interims]
  first.zero.in.row <- apply(submat, 1, function(x) match(0, x))
  m.fail <- C * as.numeric(first.zero.in.row)
  Sm.fail <- as.numeric(names(first.zero.in.row))

  #fail.deets <- data.frame(Sm=Sm.fail, m=m.fail)
  #fail.deets$prob <- apply(fail.deets, 1, function(x) {final.probs.mat[x["Sm"]+1, x["m"]]})
  #fail.deets$prob.p0 <- apply(fail.deets, 1, function(x) {final.probs.mat.p0[x["Sm"]+1, x["m"]]})
  fail.deets.prob <- numeric(length(m.fail))
  fail.deets.prob.p0 <- numeric(length(m.fail))
  for(i in 1:length(m.fail)){
    fail.deets.prob[i] <- final.probs.mat[Sm.fail[i]+1, m.fail[i]]
    fail.deets.prob.p0[i] <- final.probs.mat.p0[Sm.fail[i]+1, m.fail[i]]
  }
  # Make df of terminal points:
  fail.deets <- data.frame(Sm=Sm.fail, m=m.fail, prob=fail.deets.prob, success=rep("Fail", length(Sm.fail)))
  tp <- rbind(fail.deets, success)
  tp <- tp[tp$prob>0, ]

  #sample.size.expd <- sum(output$m*output$prob)
  #sample.size.expd.p0 <- sum(output$m*output$prob.p0)

  #sample.size.expd.fail <- sum(fail.deets$m*fail.deets$prob)
  #sample.size.expd.p0.fail <- sum(fail.deets$m*fail.deets$prob.p0)
  sample.size.expd.fail <- sum(m.fail*fail.deets.prob)
  sample.size.expd.p0.fail <- sum(m.fail*fail.deets.prob.p0)

  sample.size.expd.success <- sum(m.success*prob.success)
  sample.size.expd.p0.success <- sum(m.success*prob.success.p0)

  sample.size.expd <- sample.size.expd.fail + sample.size.expd.success
  sample.size.expd.p0 <- sample.size.expd.p0.fail + sample.size.expd.p0.success



  ############ EFFECTIVE N. THE "EFFECTIVE N" OF A STUDY IS THE "REAL" MAXIMUM SAMPLE SIZE
  ######## The point is where every Sm for a given m equals zero or one is necessarily where a trial stops

  #cp.colsums <- apply(mat, 2, function(x) { sum(x==0, na.rm=TRUE)+sum(x==1, na.rm=TRUE)} ) # Sum the CP values that equal zero or one in each column
  cp.colsums <- colSums(mat==0, na.rm = T) + colSums(mat==1, na.rm = T)
  # possible.cps <- apply(mat, 2, function(x) {sum(!is.na(x))})
  possible.cps <- colSums(!is.na(mat))
  effective.n <- min(which(cp.colsums==possible.cps))
  return.vec <- c(n, r, C, typeI, pwr, sample.size.expd.p0, sample.size.expd, thetaF, thetaE, effective.n)
  if(return.tps==TRUE){
    return.list <- list(ocs=return.vec, tp=tp)
    return(return.list)
    }
  return.vec
}

##################### FUNCTION findDesignsGivenCohortStage
findDesignsGivenCohortStage <- function(nmin,
                        nmax,
                        C=NA,
                        stages=NA,
                        p0,
                        p1,
                        alpha,
                        power,
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
  grp <- V2 <- NULL
  q0 <- 1-p0
  q1 <- 1-p1

#use.stages <- ifelse(test=is.na(C), yes = TRUE, no = FALSE)
use.stages <- is.na(C)
if(!is.na(C) & !is.na(stages)) stop("Values given for both cohort/block size C and number of stages. Please choose one only.")

# if(use.stages==TRUE){
#   if(nmin%%stages!=0) stop("nmin must be a multiple of number of stages")
#   if(nmax%%stages!=0) stop("nmax must be a multiple of number of stages")
#   nposs <- seq(from=nmin, to=nmax, by=stages)
# } else {
#   if(nmin%%C!=0) stop("nmin must be a multiple of cohort size")
#   if(nmax%%C!=0) stop("nmax must be a multiple of cohort size")
#   nposs <- seq(from=nmin, to=nmax, by=C)
#   }
n.initial.range <- nmin:nmax
if(use.stages==TRUE){
  nposs.min.index <- min(which((nmin:nmax)%%stages==0)) # nmin must be a multiple of number of stages
  nposs.max.index <- max(which((nmin:nmax)%%stages==0)) # nmax must be a multiple of number of stages
  if(any(is.infinite(nposs.min.index), is.infinite(nposs.max.index))) stop("There must be at least one value within [nmin, nmax] that is a multiple of the number of stages")
  nposs.min <- n.initial.range[nposs.min.index]
  nposs.max <- n.initial.range[nposs.max.index]
  nposs <- seq(from=nposs.min, to=nposs.max, by=stages)
}else{
  nposs.min.index <- min(which((nmin:nmax)%%C==0)) # nmin must be a multiple of cohort size
  nposs.max.index <- max(which((nmin:nmax)%%C==0)) # nmax must be a multiple of cohort size
  if(any(is.infinite(nposs.min.index), is.infinite(nposs.max.index))) stop("There must be at least one value within [nmin, nmax] that is a multiple of the cohort size")
  nposs.min <- n.initial.range[nposs.min.index]
  nposs.max <- n.initial.range[nposs.max.index]
  nposs <- seq(from=nposs.min, to=nposs.max, by=C)
}



  r <- list()
  for(i in 1:length(nposs))
  {
    r[[i]] <- 0:(nposs[i]-1) # r values: 0 to nposs[i]-1
  }

  ns <- NULL

  for(i in 1:length(nposs)){
    ns <- c(ns, rep(nposs[i], length(r[[i]])))
  }

  sc.subset <- data.frame(n=ns, r=unlist(r))

  if(!is.na(bounds)){
    # Incorporate A'Hern's bounds:
    if(bounds=="ahern")  {
      sc.subset <- sc.subset[sc.subset$r >= p0*sc.subset$n & sc.subset$r <= p1*sc.subset$n, ]
    }

    if(bounds=="wald"){
    # Faster to incorporate Wald's bounds:
      denom <- log(p1/p0) - log((1-p1)/(1-p0))
      accept.null <- log((1-power)/(1-alpha)) / denom  + nposs * log((1-p0)/(1-p1))/denom
      accept.null <- floor(accept.null)

      reject.null <- log((power)/alpha) / denom  + nposs * log((1-p0)/(1-p1))/denom
      reject.null <- ceiling(reject.null)

      r.wald <- NULL
      ns.wald <- NULL
      for(i in 1:length(nposs)){
        r.wald <- c(r.wald, accept.null[i]:reject.null[i])
        ns.wald <- c(ns.wald, rep(nposs[i], length(accept.null[i]:reject.null[i])))
      }
      sc.subset <- data.frame(n=ns.wald, r=r.wald)
    }
  }

  # In case user wishes to specify values for r:
  if(!is.na(fixed.r)){
    sc.subset <- sc.subset[sc.subset$r %in% fixed.r,]
  }

  sc.subset <- sc.subset[sc.subset$r>=0, ]
  sc.subset <- sc.subset[sc.subset$r<sc.subset$n, ]

  if(use.stages==TRUE){
    sc.subset$C <- sc.subset$n/stages
  } else {
    sc.subset$C <- C
  }

  ###### Find thetas for each possible {r, N} combn:
  mat.list <- vector("list", nrow(sc.subset))
  for(i in 1:nrow(sc.subset)){
    mat.list[[i]] <- findCPmatrix(n=sc.subset[i,"n"], r=sc.subset[i,"r"], Csize=sc.subset[i,"C"], p0=p0, p1=p1)
  }

  store.all.thetas <- lapply(mat.list, function(x) {sort(unique(c(x))[unique(c(x)) <= 1])})

if(!is.na(max.combns) & is.na(maxthetas)){ # Only use max.combns if maxthetas is not specified.
  maxthetas <- sqrt(2*max.combns)
}

  if(is.na(exact.thetaF)){
    for(i in 1:nrow(sc.subset)){
      current.theta.vec <- store.all.thetas[[i]]
      store.all.thetas[[i]] <- current.theta.vec[current.theta.vec >= minthetaE | current.theta.vec <= maxthetaF] # Subset thetas BEFORE thinning
      while(length(store.all.thetas[[i]]) > maxthetas){
        every.other.element <- rep(c(FALSE, TRUE), 0.5*(length(store.all.thetas[[i]])-2))
        store.all.thetas[[i]] <-  store.all.thetas[[i]][c(TRUE, every.other.element, TRUE)]
      }
    }
    all.theta.combns <- lapply(store.all.thetas, function(x) {t(combn(x, 2)) })
    for(i in 1:nrow(sc.subset)){
      all.theta.combns[[i]] <- all.theta.combns[[i]][all.theta.combns[[i]][,2] >= minthetaE & all.theta.combns[[i]][,1] <= maxthetaF, ] # Remove ordered pairs where both are low or high
      if(length(all.theta.combns[[i]])==2) all.theta.combns[[i]] <- matrix(c(0, 0.999, 0, 1), nrow=2, byrow=T) # To avoid a crash. See (*) below
      all.theta.combns[[i]] <- all.theta.combns[[i]][order(all.theta.combns[[i]][,2], decreasing=T),]
    }
    # (*) If the only thetas for a design are 0 and 1, then some code below will crash. This could be resolved by adding an if or ifelse statement,
    # but probably at unnecessary computational expense. We deal with these exceptions here, by adding a "false" value for theta.
    }else{ # if exact thetas are given (to speed up a result check):
      for(i in 1:length(store.all.thetas)){
        keep <- abs(store.all.thetas[[i]]-exact.thetaF)<1e-3 | abs(store.all.thetas[[i]]-exact.thetaE)<1e-3
        store.all.thetas[[i]] <- store.all.thetas[[i]][keep]
      }
      all.theta.combns <- lapply(store.all.thetas, function(x) {t(combn(x, 2)) })
      }

 # find


  ##### Matrix of coefficients: Probability of a SINGLE path leading to each point:
  coeffs <- matrix(NA, ncol=nmax, nrow=nmax+1)
  coeffs.p0 <- coeffs

  for(i in 1:(nmax+1)){
    coeffs[i, ] <- p1^(i-1) * q1^((2-i):(2-i+nmax-1))
    coeffs.p0[i, ] <- p0^(i-1) * q0^((2-i):(2-i+nmax-1))
  }


  h.results.list <- vector("list", nrow(sc.subset)) #

  if(progressBar==TRUE)  pb <- txtProgressBar(min = 0, max = nrow(sc.subset), style = 3)

  ######### START HERE

  # Now, find the designs, looping over each possible {r, N} combination, and within each {r, N} combination, loop over all combns of {thetaF, thetaE}:
  for(h in 1:nrow(sc.subset)){
    k <- 1
    h.results <- vector("list", nrow(all.theta.combns[[h]]))

    current.theta.combns <- all.theta.combns[[h]] # Take just the thetaF/1 combns for that design.
    # Don't need a matrix of all thetaF and thetaE combns -- potentially quicker to have thetaF as a vector (already defined), and the list can be a list of thetaE vectors:
    current.theta.combns <- data.table::data.table(current.theta.combns)
    current.theta.combns[, grp := .GRP, by = V2]
    data.table::setkey(current.theta.combns, grp)
    split.thetas <- current.theta.combns[, list(list(.SD)), by = grp]$V1
    thetaF.list <- lapply(split.thetas, function(x) x[, 1])
    all.thetas <- rev(store.all.thetas[[h]])[-length(store.all.thetas[[h]])] # All thetaE values, decreasing, not including the final value, thetaE=0.
    all.thetas <- all.thetas[all.thetas>=minthetaE]

    for(i in 1:length(all.thetas)){ # For each thetaE,
      thetaFs.current <- thetaF.list[[i]] # thetaF.list is a list of i vectors. The i'th vector in the list contains all possible values of thetaF for the i'th thetaE (stored as all.thetas[i])
      rows <- nrow(thetaFs.current)
      # Begin bisection method:
      a <- 1
      b <- nrow(thetaFs.current)
      d <- ceiling((b-a)/2)

      while((b-a)>1){
        output <- findDesignOCs(n=sc.subset$n[h], r=sc.subset$r[h], C=sc.subset$C[h], thetaF=as.numeric(thetaFs.current[d]), thetaE=all.thetas[i], mat=mat.list[[h]],
                                   power=power, alpha=alpha, coeffs=coeffs, coeffs.p0=coeffs.p0, p0=p0, p1=p1)

        if(output[4] <= alpha) { # type I error decreases as index (and thetaF) increases. Jump backwards if type I error is smaller than alpha, o/w forwards.
          b <- d
        } else {
          a <- d
        }
        d <- a + floor((b-a)/2)
      }

      # Take care of "edge case" where feasible design exists in the first row:
      if(a==1) { # (and also by necessity, b==2)
        output <-  findDesignOCs(n=sc.subset$n[h], r=sc.subset$r[h], C=sc.subset$C[h], thetaF=as.numeric(thetaFs.current[1]), thetaE=all.thetas[i], mat=mat.list[[h]],
                                    power=power, alpha=alpha, coeffs=coeffs, coeffs.p0=coeffs.p0, p0=p0, p1=p1)
        if(output[4] <= alpha) b <- 1 # Should we start at the first row or the second?
      }

      # We can now proceed moving sequentially from index==b (or cause a break if we wish).
      first.result <-  findDesignOCs(n=sc.subset$n[h], r=sc.subset$r[h], C=sc.subset$C[h], thetaF=as.numeric(thetaFs.current[b]), thetaE=all.thetas[i],
                                        mat=mat.list[[h]], coeffs.p0=coeffs.p0, coeffs=coeffs, power=power, alpha=alpha, p0=p0, p1=p1)
          if((first.result[4]!=0 | first.result[4]!=2) ) {
        pwr <- first.result[5]
        while(pwr >= power & b <= rows) # Keep going until power drops below 1-beta, i.e. no more feasible designs, or we reach the end of the data frame.
        {
          h.results[[k]] <-  findDesignOCs(n=sc.subset$n[h], r=sc.subset$r[h], C=sc.subset$C[h], thetaF=as.numeric(thetaFs.current[b]), thetaE=all.thetas[i],
                                              mat=mat.list[[h]], coeffs.p0=coeffs.p0, coeffs=coeffs, power=power, alpha=alpha, p0=p0, p1=p1)
          pwr <- h.results[[k]][5] ####
          k <- k+1
          b <- b+1
        }
      } else { # if first.result[4]==0 or 2, i.e., there are no feasible designs for this thetaE (and hence no remaining feasible designs for this n/r combn), break and move on to next n/r combn:
        break
      }

      # If there is only one possible combination of thetaF and thetaE, this is the result:
      if(rows==1){
        h.results[[k]] <- output
        k <- k+1
      }
    } # end of "i" loop

    if(progressBar==TRUE) setTxtProgressBar(pb, h)

    ####### END HERE

    h.results.df <- do.call(rbind, h.results)

    if(!is.null(h.results.df)){
      # Remove all "skipped" results:
      colnames(h.results.df) <- c("n", "r", "C", "alpha", "power", "EssH0", "Ess", "thetaF", "thetaE", "eff.n")
      h.results.df <- h.results.df[!is.na(h.results.df[, "Ess"]),]
      if(nrow(h.results.df)>0){
        # Remove dominated and duplicated designs:
        discard <- rep(NA, nrow(h.results.df))
        for(i in 1:nrow(h.results.df)){
          discard[i] <- sum(h.results.df[i, "EssH0"] > h.results.df[, "EssH0"] & h.results.df[i, "Ess"] > h.results.df[, "Ess"] & h.results.df[i, "n"] >= h.results.df[, "n"])
        }
        h.results.df <- h.results.df[discard==0,, drop=FALSE]
        h.results.list[[h]] <- h.results.df
      }
    }
  } # End of "h" loop


  full.results <- do.call(rbind, h.results.list)
  if(length(full.results)>0) {
    if(return.only.admissible==TRUE){
      # Discard all "inferior" designs:
      discard <- rep(NA, nrow(full.results))
      for(i in 1:nrow(full.results)){
        discard[i] <- sum(full.results[i, "EssH0"] > full.results[, "EssH0"] & full.results[i, "Ess"] > full.results[, "Ess"] & full.results[i, "n"] >= full.results[, "n"])
        }
      full.results <- full.results[discard==0,,drop=FALSE]
    }
    # Remove duplicates:
    duplicates <- duplicated(full.results[, c("n", "Ess", "EssH0"), drop=FALSE])
    final.results <- full.results[!duplicates,,drop=FALSE]
    stage <- final.results[,"n"]/final.results[,"C"]
    final.results <- cbind(final.results, stage)
  } else {
    final.results <- NULL
    if(nmin!=nmax){
      print("There are no feasible designs for this combination of design parameters")
      return(final.results)
    }
  }

  # input <- data.frame(nmin=nmin, nmax=nmax,
  #                     #nmin.used=nposs.min, nmax.used=nposs.max,
  #                     C=C, stages=stages, p0=p0, p1=p1, alpha=alpha, power=power,
  #                     maxthetaF=maxthetaF, minthetaE=minthetaE, bounds=bounds, fixed.r=fixed.r,
  #                     return.only.admissible=return.only.admissible, max.combns=max.combns,
  #                     maxthetas=maxthetas, exact.thetaF=exact.thetaF, exact.thetaE=exact.thetaE)
  # return(list(input=input,
  #             all.des=final.results,
  #             bounds.mat=single.des$bounds.mat,
  #             diagr=single.des$diagram
  #             ))
  return(final.results)
}


# nmin = minimum permitted sample size.
# nmax = maximum permitted sample size
# C = Block size. Default=1, i.e. continuous monitoring
# stages = number of interim analyses or "stages". Only required if not setting block size C. If using this argument, set C=NULL
# alpha = significance level
# power = required power (1-beta)
# maxthetaF = theta_F_max (max. value of lower threshold)
# minthetaE = theta_E_min (min. value of upper threshold)
# bounds = choose what set of stopping boundaries should be searched over: Those of A'Hern ("ahern"), Wald ("wald") or any (NULL).
# fixed.r = choose what stopping boundaries should be searched over (generally used only to check results for a single scalar value of r)

#' Find single-arm trials using stochastic curtailment
#'
#' This function finds admissible design realisations for single-arm binary outcome trials, using stochastic curtailment.
#' The output can be used as the sole argument in the function 'drawDiagram', which will return the stopping boundaries for the
#' admissible design of your choice. Monitoring frequency can set in terms of block(/cohort) size ("C") or number of stages ("stages").
#' @param nmin Minimum permitted sample size. Should be a multiple of block size or number of stages.
#' @param nmax Maximum permitted sample size. Should be a multiple of block size or number of stages.
#' @param C Block size. Vectors, i.e., multiple values, are permitted.
#' @param stages Number of interim analyses or "stages". Only required if not setting block size C. Vectors, i.e., multiple values, are permitted.
#' @param p0 Probability for which to control the type-I error-rate
#' @param p1 Probability for which to control the power
#' @param alpha Significance level
#' @param power Required power (1-beta).
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

#out <- bounds(all.results[6,], type="mstage", p0=0.1, p=0.3)
#test2<- build_gs(J = out$N, n = 1:out$N, a = out$a, r = out$r, pi0 = 0.1, pi1 = 0.3, alpha = out$alpha, beta = out$beta)
#write.csv(all.results, file="scen1.csv")

sc_bounds <- function(n1, n2, r1, r2, p0, p, thetaF, thetaE, alpha, beta)
{

  n <- n1+n2
  q=1-p
  q0=1-p0

  # Start with all zeros (why not):
  mat <- matrix(0, nrow = n+1, ncol = n)
  rownames(mat) <- 0:n

  # Add the 1's for CP=1:
  mat[(r2+2):(n+1),] <- 1  # r+2 === Sm=r+1, n+1 === Sm=n


  # Now input CP for points where r1 < k < r+1 (ie progression to Stage 2 is guaranteed):
  Sm <- 0:n # possible k~ values.
  cp.subset1.Sm <- which(r1 < Sm & Sm < (r2+1)) - 1  # Values for Sm where r1 < Sm < r2+1, i.e. progression to S2 certain, but success not guaranteed.
  ### ^^^ XXX CHECK THIS XXX -- the inequality!!!


  # Which columns? Trial success must still be possible, ie  n-n~ > r-Sm+1
  m <- 1:n # possible numbers of patients

  cp.subset1.m <- list()

  i <- 1

  for(k in cp.subset1.Sm)
  {
    cp.subset1.m[[i]] <- which(n-m >= r2-k+1 & m>=k)
    i <- i+1
  }
  # These are the m's for which trial success is still possible when Sm is cp.subset1.Sm
  # Note: This code seems faster than using "cp.subset1.m <- lapply(cp.subset1.Sm, function(x) {which(n-m >= r2-x+1 & m>=x)})"


  ##### Now calculate CP for all these points:

  # For each Sm, begin with the earliest point, ie lowest m, and work row-wise to the latest point, ie highest m:
  # Remember that row number is Sm+1.
  # How many rows to do? length(cp.subset1.m):

  summation <- list()

  for(j in 1:length(cp.subset1.m))
  {
    current.Sm <- cp.subset1.Sm[j]

    l <- seq(from=r2-current.Sm, length=length(cp.subset1.m[[j]]))

    summation[[j]] <- choose(l, r2-current.Sm) * p^(r2-current.Sm+1) * q^(l-(r2-current.Sm))
  }

  cp1 <- lapply(summation, cumsum)
  cp1 <- lapply(cp1, rev)


  # Now, insert CPs into matrix of CPs:

  for(i in 1:length(cp.subset1.Sm))  {  mat[cp.subset1.Sm[i]+1, cp.subset1.m[[i]]] <- cp1[[i]] }


  # Final region to be addressed: Points where k <= r1 but progression is still possible (S1 only):
  cp.subset2.Sm <- 0:r1  # Values of Sm <= r1

  # Which columns? Progression to S2 must still be possible, ie  n1-m > r1-Sm+1
  n.to.n1 <- 1:n1 # possible S1 values

  cp.subset2.m <- list()
  i <- 1
  for(k in cp.subset2.Sm)
  {
    cp.subset2.m[[i]] <- which(n1-n.to.n1 >= r1-k+1 & n.to.n1 >= k)
    i <- i+1
  }

  # Now we have the rows and columns of this region:

  # cp.subset2.Sm # Values of Sm. Rows are this plus 1.
  # cp.subset2.m # Values of m (columns)

  # Have to propagate "backwards", ending at Sm=0, n=1.
  # As with previous region, proceed row-wise -- start with greatest Sm.

  for(k in length(cp.subset2.Sm):1) # Note that we END at the earliest row.
  {
    for(j in rev(cp.subset2.m[[k]]))
    {
      mat[k, j] <- q*mat[k, j+1] + p*mat[k+1, j+1]
    }
  }
  # This *must* be done in this way -- cannot be vectorised.

  # Create the "blanks":
  for(i in 1:(ncol(mat)-1))
  {
    mat[(i+2):nrow(mat), i] <- NA
  }

  # all.thetas <- unique(c(mat))
  # subset.thetas <- all.thetas[all.thetas < 0.2]

  # We now have a matrix m containing the CPs in (more or less) the upper triangle.



  ######IMPORTANT


  # If a point in a path has CP < theta, we stop the trial due to stochastic curtailment.
  # If an earlier point in the path leads only to curtailment of some kind -- non-stochastic, stochastic, or both --
  # then the CP of that earlier point is essentially zero. There is no point in "reaching" it.
  # Similarly, if an earlier point in the path *may* lead to curtailment of some kind, the CP of that earlier point must change.

  # With the above in mind, it seems that it is not possible to separate the calculation of CP from the comparison of CP against theta;
  # They must be undertaken together.


  ###### NEXT SECTION: CHECKING IF EACH CP<thetaF OR CP>thetaE

  # Begin at row r+1, i.e. where Sm=r, and at column n-1.
  # Proceed right to left then bottom to top, ignoring cases where CP=0 or CP=1.
  # As CP increases from right to left, if CP>=theta then can move on to next row (because all CPs will be >=theta).

  # The CPs in regions A and B have already been defined above; these are the points
  # for which 0 < CP < 1. Combine these regions and examine:

  cp.sm <- c(cp.subset2.Sm, cp.subset1.Sm)
  cp.m <- c(cp.subset2.m, cp.subset1.m)
  # ^ These are all the points that it is necessary to cycle through.


  coeffs <- list()
  coeffs.p0 <- list()

  for(i in 1:n){
    j <- 1:(i+1)
    coeffs[[i]] <- p^(j-1)*q^(i+1-j)
    coeffs.p0[[i]] <- p0^(j-1)*q0^(i+1-j)
  }


  ########### END OF OUTSIDE THETA FN


  # To reduce looping, identify the rows which contain a CP < thetaF OR CP > thetaE:

  theta.test <- function(SM, M)
  {
    sum(mat[SM+1, M] < thetaF | mat[SM+1, M] > thetaE)
  }


  low.cp.rows <- mapply(theta.test, SM=cp.sm, M=cp.m)
  # Because of how CP changes, propagating right to left and upwards, we only
  # need to identify the "lowest" row with CP < theta or CP > thetaE, ie the row with greatest Sm
  # that contains a CP < theta or CP > thetaE. This is the row where we begin. Note: This may be all rows, or even none!
  # Recall that row number = Sm-1

  # mat.old <- mat

  # We only need to act if there are CPs < thetaF or > thetaE -- otherwise, the matrix of CPs does not change.
  if(sum(low.cp.rows)>0)    # There may be some cases where there are no CPs to be changed; i.e. no stopping due to stochastic curtailment
  {
    begin.at.row <- max(which(low.cp.rows>0))

    cp.changing.sm <- cp.sm[1:begin.at.row]
    cp.changing.m <- cp.m[1:begin.at.row]


    # We calculate the CP, truncate to zero if CP<thetaF (or to 1 if CP > thetaE), then move on:

    for(rown in rev(cp.changing.sm+1))
    {
      for(coln in rev(cp.changing.m[[rown]]))
      {
        currentCP <- q*mat[rown, coln+1] + p*mat[rown+1, coln+1]
        if(currentCP > thetaE)   mat[rown, coln] <- 1   # If CP > thetaE, amend to equal 1
        else  mat[rown, coln] <- ifelse(test = currentCP < thetaF, yes=0, no=currentCP)  # Otherwise, test if CP < thetaF. If so, amend to 0, otherwise calculate CP as normal
      }
    } # Again, this *must* be done one entry at a time -- cannot vectorise.
  } # End of IF statement


  ###### At this point, the matrix "mat" contains the CPs, adjusted for stochastic curtailment.
  ###### The points in the path satisfying 0 < CP < 1 are the possible points for this design;
  ###### All other points are either terminal (i.e. points of curtailment) or impossible.
  ###### We now need the characteristics of this design: Type I error, expected sample size, and so on.

  pascal.list <- list(1, c(1,1))


  for(i in 3:(n+2))
  {
    column <- as.numeric(mat[!is.na(mat[,i-2]), i-2])
    CPzero.or.one <- which(column==0 | column==1)
    newnew <- pascal.list[[i-1]]
    newnew[CPzero.or.one] <- 0
    pascal.list[[i]] <- c(0, newnew) + c(newnew, 0)
  }

  pascal.list <- pascal.list[c(-1, -length(pascal.list))]


  # Multiply the two triangles (A and p^b * q^c):
  final.probs <- Map("*", pascal.list, coeffs)

  # for finding type I error prob:
  final.probs.p0 <- Map("*", pascal.list, coeffs.p0)


  ###### We have the probability of each path, taking into account stochastic and NS curtailment.
  ###### We now must tabulate these paths.
  # Might be best to sort into Failure (Sm <= r) and Success (Sm = r+1):

  final.probs.mat <- matrix(unlist(lapply(final.probs, '[', 1:max(sapply(final.probs, length)))), ncol = n, byrow = F)
  #rownames(final.probs.mat) <- 0:n

  final.probs.mat.p0 <- matrix(unlist(lapply(final.probs.p0, '[', 1:max(sapply(final.probs.p0, length)))), ncol = n, byrow = F)
  #rownames(final.probs.mat.p0) <- 0:n


  # FIRST: Search for terminal points of success. These can only exist in rows where (updated) CP=1, and where Sm<=r+1:
  potential.success.rows <- rowSums(mat[1:(r2+2), ]==1, na.rm = TRUE)
  rows.with.cp1 <- which(as.numeric(potential.success.rows)>0)
  # ^ These are the rows containing possible terminal points of success. They must have CP=1:

  columns.of.rows.w.cp1 <- list()
  j <- 1

  for(i in rows.with.cp1)
  {
    columns.of.rows.w.cp1[[j]] <- which(mat[i, ]==1 & !is.na(mat[i, ]))
    j <- j+1
  }

  # These rows and columns contain all possible terminal points of success.
  # The point CP(Sm, m) is terminal if CP(Sm-1, m-1) < 1 .
  # Strictly speaking, CP(Sm, m) is also terminal if CP(Sm, m-1) < 1 .
  # However, CP(Sm, m-1) >= CP(Sm-1, m-1)  [I think], so the case of
  # CP(Sm, m) == 1  AND  CP(Sm, m-1) < 1 is not possible.


  index1 <- 1
  success <- NULL


  for(i in rows.with.cp1)
  {
    for(j in columns.of.rows.w.cp1[[index1]])
    {
      #print(i);print(j)
      if(j==1) success <- rbind(success, c(i-1, j, final.probs.mat[i, j], final.probs.mat.p0[i, j]))
      else
        if(mat[i-1, j-1] < 1) success <- rbind(success, c(i-1, j, final.probs.mat[i, j], final.probs.mat.p0[i, j]))
    }
    index1 <- index1 + 1
  }

  colnames(success) <- c("Sm", "m", "prob", "prob.p0")



  # Now failure probabilities. Note that there AT MOST one failure probability in each row, and in that
  # row the failure probability is the one that has the greatest m (i.e. the "furthest right" non-zero entry):

  # First, in the matrix of CPs, find the earliest column containing only zeros and ones (and NAs): This is the latest point ("m") at which the trial stops:
  first.zero.one.col <- which(apply(mat, 2, function(x) sum(!(x %in% c(0,1, NA))))==0)[1]

  # Second, find the the row at which the final zero exists in this column:
  final.zero.in.col <- max(which(mat[, first.zero.one.col]==0))

  # Finally, find the rightmost non-zero probability in each row up and including the one above:

  m.fail <- NULL
  prob.fail <- NULL
  prob.fail.p0 <- NULL

  for(i in 1:final.zero.in.col)
  {
    m.fail[i] <- max(which(final.probs.mat[i ,]!=0))
    prob.fail[i] <- final.probs.mat[i, m.fail[i]]
    prob.fail.p0[i] <- final.probs.mat.p0[i, m.fail[i]]
  }
  Sm.fail <- 0:(final.zero.in.col-1)


  eff.bound.Sm <- rep(Inf, n) # Inf instead of n+1 -- either way, no stopping
  eff.bound.Sm[success[,"m"]] <- success[,"Sm"]

  fut.bound.Sm <- rep(-Inf, n) # -Inf instead of 0 -- either way, no stopping
  fut.bound.Sm[m.fail] <- Sm.fail

  #to.return <- c(n1, n2, n, r1, r2, typeI, pwr, sample.size.expd.p0, sample.size.expd, thetaF, thetaE, effective.n)
  to.return <- list(J=2, n=c(n1, n2), N=n, a=c(r1, r2), r=c(Inf, r2+1), a_curt=fut.bound.Sm, r_curt=eff.bound.Sm, alpha=alpha, beta=beta)
  to.return
}


mstage_bounds <- function(n, r,thetaF, thetaE, p, p0, alpha, power, returnCPmat=FALSE)
{

  q <- 1-p
  q0<- 1-p0
  # Start with all zeros (why not):
  mat <- matrix(c(rep(0, n*(r+1)), rep(1, n*((n+1)-(r+1)))), nrow = n+1, byrow=T)
  rownames(mat) <- 0:n

  # Create the "blanks":
  NAs <- rbind(FALSE, lower.tri(mat)[-nrow(mat),])
  mat[NAs] <- NA

  # Now input CP for all points
  Sm <- 0:n # possible k~ values. # XXX not 0:(r+1) ???
  cp.subset1.Sm <- 0:r

  # Which columns? Trial success must still be possible, ie  n-m > r-Sm+1
  m <- 1:n # possible numbers of patients


  cp.subset1.m <- list()

  i <- 1

  for(k in cp.subset1.Sm)
  {
    cp.subset1.m[[i]] <- which(n-m >= r-k+1 & m>=k)
    i <- i+1
  }
  # These are the m's for which trial success is still possible when Sm is cp.subset1.Sm
  # Note: This code seems faster than using "cp.subset1.m <- lapply(cp.subset1.Sm, function(x) {which(n-m >= r2-x+1 & m>=x)})"


  ##### Now calculate CP for all these points:

  # For each Sm, begin with the earliest point, ie lowest m, and work row-wise to the latest point, ie highest m:
  # Remember that row number is Sm+1.
  # How many rows to do? length(cp.subset1.m):

  summation <- list()

  for(j in 1:length(cp.subset1.m))
  {
    current.Sm <- cp.subset1.Sm[j] # Sm of current row

    l <- seq(from=r-current.Sm, length=length(cp.subset1.m[[j]]))

    summation[[j]] <- choose(l, r-current.Sm) * p^(r-current.Sm+1) * q^(l-(r-current.Sm))
  }

  cp1 <- lapply(summation, cumsum)
  cp1 <- lapply(cp1, rev)


  # Now, insert CPs into matrix of CPs:
  for(i in 1:length(cp.subset1.Sm))  {  mat[cp.subset1.Sm[i]+1, cp.subset1.m[[i]]] <- cp1[[i]] }

  cp.sm <- cp.subset1.Sm
  cp.m  <- cp.subset1.m


  coeffs <- list()
  coeffs.p0 <- list()

  for(i in 1:n){
    j <- 1:(i+1)
    coeffs[[i]] <- p^(j-1)*q^(i+1-j)
    coeffs.p0[[i]] <- p0^(j-1)*q0^(i+1-j)
  }


  ########### insideTheta


  theta.test <- function(SM, M)
  {
    sum(mat[SM+1, M] < thetaF-1e-6 | mat[SM+1, M] > thetaE + 1e-6) # The 1e-6 terms are added to deal with floating point errors.
  }


  low.cp.rows <- mapply(theta.test, SM=cp.sm, M=cp.m)
  # Because of how CP changes, propagating right to left and upwards, we only
  # need to identify the "lowest" row with CP < thetaF or CP > thetaE, ie the row with greatest Sm
  # that contains a CP < thetaF or CP > thetaE. This is the row where we begin. Note: This may be all rows, or even none!
  # Recall that row number = Sm-1



  # We only need to act if there are CPs < thetaF or > thetaE -- otherwise, the matrix of CPs does not change.
  if(sum(low.cp.rows)>0)    # There may be some cases where there are no CPs to be changed; i.e. no stopping due to stochastic curtailment
  {
    begin.at.row <- max(which(low.cp.rows>0))

    cp.changing.sm <- cp.sm[1:begin.at.row]
    cp.changing.m <- cp.m[1:begin.at.row]


    # We calculate the CP, truncate to zero if CP<thetaF (or to 1 if CP > thetaE), then move on:

    for(rown in rev(cp.changing.sm+1))
    {
      for(coln in rev(cp.changing.m[[rown]]))
      {
        currentCP <- q*mat[rown, coln+1] + p*mat[rown+1, coln+1]
        if(currentCP > thetaE + 1e-6)   mat[rown, coln] <- 1   # If CP > thetaE, amend to equal 1. The term 1e-6 is added to deal with floating point errors.
        #if(currentCP > thetaE) {  mat[rown, coln] <- 1 ; print(currentCP) } # If CP > thetaE, amend to equal 1
        else  mat[rown, coln] <- ifelse(test = currentCP < thetaF - 1e-6, yes=0, no=currentCP)  # Otherwise, test if CP < thetaF. If so, amend to 0, otherwise calculate CP as normal. Again, the term 1e-6 is added to deal with floating point errors.
      }
    } # Again, this *must* be done one entry at a time -- cannot vectorise.
  } # End of IF statement



  pascal.list <- list(1, c(1,1))


  for(i in 3:(n+2))
  {
    column <- as.numeric(mat[!is.na(mat[,i-2]), i-2])
    CPzero.or.one <- which(column==0 | column==1)
    newnew <- pascal.list[[i-1]]
    newnew[CPzero.or.one] <- 0
    pascal.list[[i]] <- c(0, newnew) + c(newnew, 0)
  }

  pascal.list <- pascal.list[c(-1, -length(pascal.list))]


  # Multiply the two triangles (A and p^b * q^c):
  final.probs <- Map("*", pascal.list, coeffs)

  # for finding type I error prob:
  final.probs.p0 <- Map("*", pascal.list, coeffs.p0)


  ###### We have the probability of each path, taking into account stochastic and NS curtailment.
  ###### We now must tabulate these paths.
  final.probs.mat <- matrix(unlist(lapply(final.probs, '[', 1:max(sapply(final.probs, length)))), ncol = n, byrow = F)

  final.probs.mat.p0 <- matrix(unlist(lapply(final.probs.p0, '[', 1:max(sapply(final.probs.p0, length)))), ncol = n, byrow = F)


  # FIRST: Search for terminal points of success. These can only exist in rows where (updated) CP=1, and where Sm<=r+1:
  potential.success.rows <- rowSums(mat[1:(r+2), ]==1, na.rm = TRUE)
  rows.with.cp1 <- which(as.numeric(potential.success.rows)>0)
  # ^ These are the rows containing possible terminal points of success. They must have CP=1:

  columns.of.rows.w.cp1 <- list()
  j <- 1

  for(i in rows.with.cp1)
  {
    columns.of.rows.w.cp1[[j]] <- which(mat[i, ]==1 & !is.na(mat[i, ]))
    j <- j+1
  }


  index1 <- 1
  success <- NULL

  for(i in rows.with.cp1)
  {
    for(j in columns.of.rows.w.cp1[[index1]])
    {
      if(j==1) success <- rbind(success, c(i-1, j, final.probs.mat[i, j], final.probs.mat.p0[i, j]))
      else
        if(mat[i-1, j-1] < 1) success <- rbind(success, c(i-1, j, final.probs.mat[i, j], final.probs.mat.p0[i, j]))
    }
    index1 <- index1 + 1
  }

  colnames(success) <- c("Sm", "m", "prob", "prob.p0")


  # Now failure probabilities. Note that there AT MOST one failure probability in each row, and in that
  # row the failure probability is the one that has the greatest m (i.e. the "furthest right" non-zero entry):

  # First, in the matrix of CPs, find the earliest column containing only zeros and ones (and NAs): This is the latest point ("m") at which the trial stops:
  first.zero.one.col <- which(apply(mat, 2, function(x) sum(!(x %in% c(0,1, NA))))==0)[1]

  # Second, find the the row at which the final zero exists in this column:
  final.zero.in.col <- max(which(mat[, first.zero.one.col]==0))

  # Finally, find the rightmost non-zero probability in each row up and including the one above:

  m.fail <- NULL
  prob.fail <- NULL
  prob.fail.p0 <- NULL

  for(i in 1:final.zero.in.col)
  {
    m.fail[i] <- max(which(final.probs.mat[i ,]!=0))
    prob.fail[i] <- final.probs.mat[i, m.fail[i]]
    prob.fail.p0[i] <- final.probs.mat.p0[i, m.fail[i]]
  }
  Sm.fail <- 0:(final.zero.in.col-1)

  fail.deets <- cbind(Sm.fail, m.fail, prob.fail, prob.fail.p0)

  eff.bound.Sm <- rep(Inf, n) # Inf instead of n+1 -- either way, no stopping
  eff.bound.Sm[success[,"m"]] <- success[,"Sm"]

  fut.bound.Sm <- rep(-Inf, n) # -Inf instead of n+1 -- either way, no stopping
  fut.bound.Sm[m.fail] <- Sm.fail

  to.return <- list(J=1, n=n, N=n, a=r, r=r+1, a_curt=fut.bound.Sm, r_curt=eff.bound.Sm, alpha=alpha, beta=1-power)
  if(returnCPmat==TRUE) to.return$mat <- mat
  to.return
}



bounds <- function(des, type=c("simon", "simon_e1", "nsc", "sc", "mstage"), p0, p)
{
  n1 <- des[["n1"]]
  n2 <- des[["n2"]]
  r1 <- des[["r1"]]
  r <- des[["r"]]
  N <- des[["n"]]
  alpha <- des[["alpha"]]
  beta <- 1-des[["power"]]

  if(type=="simon"){
    result <- list(J=2, n=c(n1, n2), N=N, a=c(r1, r), r=c(Inf, r+1), alpha=alpha, beta=beta, n1=n1) # Or N+1 instead of Inf -- either way, no stopping
  }
  if(type=="simon_e1"){
    result <- list(J=2, n=c(n1, n2), N=N, a=c(r1, r), r=c(des$e1, r+1), alpha=alpha, beta=beta, n1=n1)
  }
  if(type=="nsc"){
    eff.bound.Sm <- c(rep(Inf, r), rep(r+1, N-r)) # Stop for efficacy at r+1, at every point. Also, or N+1 instead of Inf -- either way, no stopping
    # Stage 1
    fut.bound.s1.m.partial <- (n1-r1):n1 # Values of m where S1 stopping for futility is possible.
    fut.bound.s1.m <- 1:n1
    fut.bound.s1.Sm <- c(rep(-Inf, fut.bound.s1.m.partial[1]-1), seq(from=0, length=length(fut.bound.s1.m.partial))) # or 0 instead of -Inf -- either way, no stopping
    # Stage 2:
    fut.bound.s2.m.partial <- max(N+r1-r+1, n1+1):N
    fut.bound.s2.m <- (n1+1):N
    fut.bound.s2.Sm <- c(rep(-Inf, length(fut.bound.s2.m)-length(fut.bound.s2.m.partial)), seq(from=r1+1, length=length(fut.bound.s2.m.partial))) # or 0 instead of -Inf -- either way, no stopping
    # Output
    result <- list(J=2, n=c(n1, n2), N=N, a=c(r1, r), r=c(Inf, r+1), a_curt=c(fut.bound.s1.Sm, fut.bound.s2.Sm), r_curt=eff.bound.Sm, alpha=alpha, beta=beta)
  }

  if(type=="sc"){
    result <- sc_bounds(n1=n1, n2=des$n2, r1=r1, r2=r, p0=p0, p=p, thetaF=des$thetaF, thetaE=des$thetaE, alpha=des$alpha, beta=1-des$power)
  }
  if(type=="mstage"){
    result <- mstage_bounds(n=des$n, r=des$r, thetaF=des$thetaF, thetaE=des$thetaE, p0=p0, p=p, power=des$power, alpha=des$alpha)
  }

  result
}

#### Below: Fns for finding SC designs ####

#### Part 1: from "find_N1_N2_R1_R_df.R" ####

################################# THIS CODE CREATES ALL COMBNS OF N1, N2, N, R, R


###### OBJECTS
# ns: data frame containing all combinations of n1/n2/n
# n1.list2: Vector of n1's from this data frame
# n2.list2: Vector of n2's from this data frame
# n.list2:  Vector of n's from this data frame
#
# Note that we start from the greatest n and decrease. This is from a "deprecated" idea of mine re. searching,
#   but doesn't affect anything that follows.
#
# r1: List of vectors. Each vector contains all possibilities of r1 for the same row in ns. For example,
#     r1[[i]] is a vector containing all possibilities of r1 for the combination of n1/n2/n in the row ns[i, ].
#
# r: List of list of vectors. Each vector contains all possibilities of r for the same row in ns and single element in r1.
#    For example, r[[1]][[1]] contains all possibilities for r for ns[1, ] and r1[[1]][1], while
#    r[[1]] will contain all possibilities for r for n1[1] and n2[1] for all r1[[1]].
#
#


findSimonN1N2R1R2 <- function(nmin, nmax, e1=TRUE)
{

  nposs <- nmin:nmax

  n1.list <- list()
  n2.list <- list()

  for(i in 1:length(nposs))
  {
    n1.list[[i]] <- 1:(nposs[i]-1)
    n2.list[[i]] <- nposs[i]-n1.list[[i]]
  }

  # All possibilities together:

  n1 <- rev(unlist(n1.list))
  n2 <- rev(unlist(n2.list))

  n <- n1 + n2

  ns <- cbind(n1, n2, n)

  #colnames(ns) <- c("n1", "n2", "n")



  ################################ FIND COMBNS OF R1 AND R
  #
  # For each possible total N, there is a combination of possible n1 and n2's (above).
  #
  # For each possible combination of n1 and n2, there is a range of possible r1's.
  #
  # Note constraints for r1 and r:
  #
  # r1 < n1
  # r1 < r < n
  # r2-r1 <= n2
  #
  #
  # Note: Increasing r1 (or r) increases P(failure), and so decreases power and type I error;
  #       Decreasing r1 (or r) decreases P(failure), and so increases power and type I error.
  #
  # But how does power and type I error change WRT to each individually? Must find out to be able to search more intelligently.
  #



  r1 <- NULL

  for(i in 1:nrow(ns))
  {
    r1 <- c(r1, 0:(n1[i]-1)) # r1 values: 0 to n1-1
  }

  rownames(ns) <- 1:nrow(ns)

  ns <- ns[rep(row.names(ns), n1), ] # duplicate each row n1 times (0 to (n1-1))

  ns <- cbind(ns, r1)



  ######### Add possible r values


  r <- apply(ns, 1, function(x) {(x["r1"]+1):min(x["n"]-1, x["n"]-x["n1"]+x["r1"])})
  # r must satisfy r > r1 and r < n. Also, number of responses required in stage 2 (r2-r1) must be at most n2

  how.many.rs <- sapply(r, length)

  row.names(ns) <- 1:nrow(ns)

  ns <- ns[rep(row.names(ns), how.many.rs), ] # duplicate each row a certain number of times

  ns <- cbind(ns, unlist(r))

  colnames(ns)[5] <- "r"

  ### Finally, add e1 for stopping for benefit:

  if(e1==TRUE)
  {
    # If stopping for benefit, r1 values must run only from 0 to (n1-2), rather than 0 to (n1-1) if not stopping for benefit, otherwise there is no possible value for e1:
    ns <- ns[ns[, "n1"]-ns[, "r1"] >= 2, ]
    e1 <- apply(ns, 1, function(x) {(x["r1"]+1):(x["n1"]-1)})  # e1 must be constrained to the interval [r1+1, n1-1]
    how.many.e1s <- sapply(e1, length)

    row.names(ns) <- 1:nrow(ns)

    ns <- ns[rep(row.names(ns), how.many.e1s), ] # duplicate each row a certain number of times

    row.names(ns) <- 1:nrow(ns)

    ns <- data.frame(ns, e1=unlist(e1))
  } else {
    rownames(ns) <- 1:nrow(ns)
    ns <- data.frame(ns)
  }

  return(ns)
}

# Could possibly speed things up

#### Part 2: from "outside_theta.R" ####
# Note: Can run much of the code outside the loops for theta -- independent of theta:


outsideTheta <- function(n1=4, n2=4, r1=1, r2=4, p0, p)
{

  n <- n1+n2
  q=1-p
  q0=1-p0

  # Start with all zeros (why not):
  mat <- matrix(0, nrow = n+1, ncol = n)
  rownames(mat) <- 0:n

  # Add the 1's for CP=1:
  mat[(r2+2):(n+1),] <- 1  # r+2 === Sm=r+1, n+1 === Sm=n




  # Now input CP for points where r1 < k < r+1 (ie progression to Stage 2 is guaranteed):
  Sm <- 0:n # possible k~ values.
  cp.subset1.Sm <- which(r1 < Sm & Sm < (r2+1)) - 1  # Values for Sm where r1 < Sm < r2+1, i.e. progression to S2 certain, but success not guaranteed.
  ### ^^^ XXX CHECK THIS XXX -- the inequality!!!


  # Which columns? Trial success must still be possible, ie  n-n~ > r-Sm+1
  m <- 1:n # possible numbers of patients

  cp.subset1.m <- list()

  i <- 1

  for(k in cp.subset1.Sm)
  {
    cp.subset1.m[[i]] <- which(n-m >= r2-k+1 & m>=k)
    i <- i+1
  }
  # These are the m's for which trial success is still possible when Sm is cp.subset1.Sm
  # Note: This code seems faster than using "cp.subset1.m <- lapply(cp.subset1.Sm, function(x) {which(n-m >= r2-x+1 & m>=x)})"


  ##### Now calculate CP for all these points:

  # For each Sm, begin with the earliest point, ie lowest m, and work row-wise to the latest point, ie highest m:
  # Remember that row number is Sm+1.
  # How many rows to do? length(cp.subset1.m):

  summation <- list()

  for(j in 1:length(cp.subset1.m))
  {
    current.Sm <- cp.subset1.Sm[j]

    l <- seq(from=r2-current.Sm, length=length(cp.subset1.m[[j]]))

    summation[[j]] <- choose(l, r2-current.Sm) * p^(r2-current.Sm+1) * q^(l-(r2-current.Sm))
  }

  cp1 <- lapply(summation, cumsum)
  cp1 <- lapply(cp1, rev)


  # Now, insert CPs into matrix of CPs:

  for(i in 1:length(cp.subset1.Sm))  {  mat[cp.subset1.Sm[i]+1, cp.subset1.m[[i]]] <- cp1[[i]] }


  # Final region to be addressed: Points where k <= r1 but progression is still possible (S1 only):
  cp.subset2.Sm <- 0:r1  # Values of Sm <= r1

  # Which columns? Progression to S2 must still be possible, ie  n1-m > r1-Sm+1
  n.to.n1 <- 1:n1 # possible S1 values

  cp.subset2.m <- list()
  i <- 1
  for(k in cp.subset2.Sm)
  {
    cp.subset2.m[[i]] <- which(n1-n.to.n1 >= r1-k+1 & n.to.n1 >= k)
    i <- i+1
  }

  # Now we have the rows and columns of this region:

  # cp.subset2.Sm # Values of Sm. Rows are this plus 1.
  # cp.subset2.m # Values of m (columns)

  # Have to propagate "backwards", ending at Sm=0, n=1.
  # As with previous region, proceed row-wise -- start with greatest Sm.

  for(k in length(cp.subset2.Sm):1) # Note that we END at the earliest row.
  {
    for(j in rev(cp.subset2.m[[k]]))
    {
      mat[k, j] <- q*mat[k, j+1] + p*mat[k+1, j+1]
    }
  }
  # This *must* be done in this way -- cannot be vectorised.

  # Create the "blanks":
  for(i in 1:(ncol(mat)-1))
  {
    mat[(i+2):nrow(mat), i] <- NA
  }

  # all.thetas <- unique(c(mat))
  # subset.thetas <- all.thetas[all.thetas < 0.2]

  # We now have a matrix m containing the CPs in (more or less) the upper triangle.



  ######IMPORTANT


  # If a point in a path has CP < theta, we stop the trial due to stochastic curtailment.
  # If an earlier point in the path leads only to curtailment of some kind -- non-stochastic, stochastic, or both --
  # then the CP of that earlier point is essentially zero. There is no point in "reaching" it.
  # Similarly, if an earlier point in the path *may* lead to curtailment of some kind, the CP of that earlier point must change.

  # With the above in mind, it seems that it is not possible to separate the calculation of CP from the comparison of CP against theta;
  # They must be undertaken together.


  ###### NEXT SECTION: CHECKING IF EACH CP<thetaF OR CP>thetaE

  # Begin at row r+1, i.e. where Sm=r, and at column n-1.
  # Proceed right to left then bottom to top, ignoring cases where CP=0 or CP=1.
  # As CP increases from right to left, if CP>=theta then can move on to next row (because all CPs will be >=theta).

  # The CPs in regions A and B have already been defined above; these are the points
  # for which 0 < CP < 1. Combine these regions and examine:

  cp.sm <- c(cp.subset2.Sm, cp.subset1.Sm)
  cp.m <- c(cp.subset2.m, cp.subset1.m)
  # ^ These are all the points that it is necessary to cycle through.


  coeffs <- list()
  coeffs.p0 <- list()

  for(i in 1:n){
    j <- 1:(i+1)
    coeffs[[i]] <- p^(j-1)*q^(i+1-j)
    coeffs.p0[[i]] <- p0^(j-1)*q0^(i+1-j)
  }

  return(list(n1=n1, n2=n2, n=n, r1=r1, r2=r2, cp.sm=cp.sm, cp.m=cp.m, mat=mat, coeffs.p0=coeffs.p0, coeffs=coeffs))
}


#### Part 3: from findThetas_apr.R ####
############### Varying theta

# For each design, find the potential values of theta
# Note: Function is merely a subset of sc.fast, the
# stochastic curtailment function


findThetas <- function(n1=4, n2=4, n, r1=1, r2=4, p0, p, unique.only=TRUE)
{

  n <- n1+n2
  q <- 1-p
  q0 <- 1-p0

  # Start with all zeros (why not):
  mat <- matrix(0, nrow = n+1, ncol = n)
  rownames(mat) <- 0:n

  # Add the 1's for CP=1:
  mat[(r2+2):(n+1),] <- 1  # r+2 === Sm=r+1, n+1 === Sm=n

  # Now input CP for points where r1 < k < r+1 (ie progression to Stage 2 is guaranteed):
  Sm <- 0:n # possible k~ values.
  cp.subset1.Sm <- which(r1 < Sm & Sm < (r2+1)) - 1  # Values for Sm where r1 < Sm < r2+1, i.e. progression to S2 certain, but success not guaranteed.

  # Which columns? Trial success must still be possible, ie  n-n~ > r-Sm+1
  m <- 1:n # possible numbers of patients

  cp.subset1.m <- list()

  i <- 1

  for(k in cp.subset1.Sm)
  {
    cp.subset1.m[[i]] <- which(n-m >= r2-k+1 & m>=k)
    i <- i+1
  }
  # These are the m's for which trial success is still possible when Sm is cp.subset1.Sm
  # Note: This code seems faster than using "cp.subset1.m <- lapply(cp.subset1.Sm, function(x) {which(n-m >= r2-x+1 & m>=x)})"


  ##### Now calculate CP for all these points:

  # For each Sm, begin with the earliest point, ie lowest m, and work row-wise to the latest point, ie highest m:
  # Remember that row number is Sm+1.
  # How many rows to do? length(cp.subset1.m):

  summation <- list()

  for(j in 1:length(cp.subset1.m))
  {
    current.Sm <- cp.subset1.Sm[j]

    l <- seq(from=r2-current.Sm, length=length(cp.subset1.m[[j]]))

    summation[[j]] <- choose(l, r2-current.Sm) * p^(r2-current.Sm+1) * q^(l-(r2-current.Sm))
  }

  cp1 <- lapply(summation, cumsum)
  cp1 <- lapply(cp1, rev)


  # Now, insert CPs into matrix of CPs:

  for(i in 1:length(cp.subset1.Sm))  {  mat[cp.subset1.Sm[i]+1, cp.subset1.m[[i]]] <- cp1[[i]] }


  # Final region to be addressed: Points where k <= r1 but progression is still possible (S1 only):
  cp.subset2.Sm <- 0:r1  # Values of Sm <= r1

  # Which columns? Progression to S2 must still be possible, ie  n1-m > r1-Sm+1
  n.to.n1 <- 1:n1 # possible S1 values

  cp.subset2.m <- list()
  i <- 1
  for(k in cp.subset2.Sm)
  {
    cp.subset2.m[[i]] <- which(n1-n.to.n1 >= r1-k+1 & n.to.n1 >= k)
    i <- i+1
  }

  # Now we have the rows and columns of this region:

  # cp.subset2.Sm # Values of Sm. Rows are this plus 1.
  # cp.subset2.m # Values of m (columns)

  # Have to propagate "backwards", ending at Sm=0, n=1.
  # As with previous region, proceed row-wise -- start with greatest Sm.

  for(k in length(cp.subset2.Sm):1) # Note that we END at the earliest row.
  {
    for(j in rev(cp.subset2.m[[k]]))
    {
      mat[k, j] <- q*mat[k, j+1] + p*mat[k+1, j+1]
    }
  }
  # This *must* be done in this way -- cannot be vectorised.

  # Create the "blanks":
  for(i in 1:(ncol(mat)-1))
  {
    mat[(i+2):nrow(mat), i] <- NA
  }
  all.thetas <- c(mat)
  if(unique.only==TRUE)  all.thetas <- unique(all.thetas)
  # We now have a matrix m containing the CPs in (more or less) the upper triangle.
  all.thetas <- sort(all.thetas[!is.na(all.thetas)])
  all.thetas
}


#### Part 4: from "inside_theta_oct2020.R" ####


#rm(theta.test)

insideTheta <- function(n1, n2, n, r1, r2, thetaF, thetaE, cp.sm, cp.m, p0, p, q0=1-p0, q=1-p, mat, coeffs.p0, coeffs, power, alpha){
  # To reduce looping, identify the rows which contain a CP < thetaF OR CP > thetaE:

  theta.test <- function(SM, M)
  {
    current.vals <- mat[SM+1, M]
    any(current.vals < thetaF | current.vals > thetaE)
  }


  low.cp.rows <- mapply(theta.test, SM=cp.sm, M=cp.m)
  # low.cp.rows <- numeric(length(cp.sm))
  # for(i in 1:length(low.cp.rows)){
  #   current.vals <- mat[cp.sm[i]+1, cp.m[[i]]]
  #   low.cp.rows[i] <- any(current.vals < thetaF | current.vals > thetaE)
  # }

  # Because of how CP changes, propagating right to left and upwards, we only
  # need to identify the "lowest" row with CP < theta or CP > thetaE, ie the row with greatest Sm
  # that contains a CP < theta or CP > thetaE. This is the row where we begin. Note: This may be all rows, or even none!
  # Recall that row number = Sm-1

  # mat.old <- mat

  # We only need to act if there are CPs < thetaF or > thetaE -- otherwise, the matrix of CPs does not change.
  if(sum(low.cp.rows)>0)    # There may be some cases where there are no CPs to be changed; i.e. no stopping due to stochastic curtailment
  {
    begin.at.row <- max(which(low.cp.rows>0))

    cp.changing.sm <- cp.sm[1:begin.at.row]
    cp.changing.m <- cp.m[1:begin.at.row]


    # We calculate the CP, round to zero if CP<thetaF (or to 1 if CP > thetaE), then move on:

    ###### OCT 2020: IMPORTANT: This clearly only works for continuous monitoring, i.e. C=1.
    # If you want to use less frequent monitoring, this code (and probably more) will have to change.
    for(rown in rev(cp.changing.sm+1))
    {
      mat.rown <- mat[rown, ]
      mat.rown.plus1 <- mat[rown+1, ]
      for(coln in rev(cp.changing.m[[rown]]))
      {
        currentCP <- q*mat.rown[coln+1] + p*mat.rown.plus1[coln+1]
        if(currentCP > thetaE){
          mat[rown, coln] <- 1   # If CP > thetaE, amend to equal 1
        }else{
          if(currentCP < thetaF){
            mat[rown, coln] <- 0
          }else{
            mat[rown, coln] <- currentCP
          }
        }
        mat.rown[coln] <- mat[rown, coln]
      }
    } # Again, this *must* be done one entry at a time -- cannot vectorise.
  } # End of IF statement


  ###### STOP if design is pointless, i.e either failure or success is not possible:

  # IF DESIGN GUARANTEES FAILURE (==0) or SUCCESS (==2) at n=1:
  first.patient <- sum(mat[,1], na.rm = T)

  if(first.patient==2)
  {
    return(c(n1, n2, n, r1, r2, 1, 1,  NA, NA, thetaF, thetaE, NA))
  }

  if(first.patient==0)
  {
    return(c(n1, n2, n, r1, r2, 0, 0,  NA, NA, thetaF, thetaE, NA))
  }





  ###### At this point, the matrix "mat" contains the CPs, adjusted for stochastic curtailment.
  ###### The points in the path satisfying 0 < CP < 1 are the possible points for this design;
  ###### All other points are either terminal (i.e. points of curtailment) or impossible.
  ###### We now need the characteristics of this design: Type I error, expected sample size, and so on.

  pascal.list <- list(1, c(1,1))

  for(i in 3:(n+2)){
    #column <- as.numeric(mat[!is.na(mat[,i-2]), i-2])
    current.col <- mat[, i-2]
    column <- as.numeric(current.col[!is.na(current.col)])
    #CPzero.or.one <- which(column==0 | column==1)
    CPzero.or.one <- column==0 | column==1
    newnew <- pascal.list[[i-1]]
    newnew[CPzero.or.one] <- 0
    pascal.list[[i]] <- c(0, newnew) + c(newnew, 0)
  }

  pascal.list <- pascal.list[c(-1, -length(pascal.list))]


  # Multiply the two triangles (A and p^b * q^c):
  final.probs <- Map("*", pascal.list, coeffs)

  # for finding type I error prob:
  final.probs.p0 <- Map("*", pascal.list, coeffs.p0)


  ###### We have the probability of each path, taking into account stochastic and NS curtailment.
  ###### We now must tabulate these paths.
  # Might be best to sort into Failure (Sm <= r) and Success (Sm = r+1):

  final.probs.mat <- matrix(unlist(lapply(final.probs, '[', 1:max(sapply(final.probs, length)))), ncol = n, byrow = F)
  #rownames(final.probs.mat) <- 0:n

  final.probs.mat.p0 <- matrix(unlist(lapply(final.probs.p0, '[', 1:max(sapply(final.probs.p0, length)))), ncol = n, byrow = F)
  #rownames(final.probs.mat.p0) <- 0:n

  # Find successful probabilities first:
  # m.success <- (r2+1):n
  # Sm.success <- rep(r2+1, length(m.success))
  # prob.success <- final.probs.mat[r2+2, m.success]
  # prob.success.p0 <- final.probs.mat.p0[r2+2, m.success]
  # success.deets <- cbind(Sm.success, m.success, prob.success, prob.success.p0)

  # FIRST: Search for terminal points of success. These can only exist in rows where (updated) CP=1, and where Sm<=r+1:
  potential.success.rows <- rowSums(mat[1:(r2+2), ]==1, na.rm = TRUE)
  rows.with.cp1 <- which(as.numeric(potential.success.rows)>0)
  # ^ These are the rows containing possible terminal points of success. They must have CP=1:

  columns.of.rows.w.cp1 <- list()
  j <- 1

  for(i in rows.with.cp1)
  {
    columns.of.rows.w.cp1[[j]] <- which(mat[i, ]==1 & !is.na(mat[i, ]))
    j <- j+1
  }

  # These rows and columns contain all possible terminal points of success.
  # The point CP(Sm, m) is terminal if CP(Sm-1, m-1) < 1 .
  # Strictly speaking, CP(Sm, m) is also terminal if CP(Sm, m-1) < 1 .
  # However, CP(Sm, m-1) >= CP(Sm-1, m-1)  [I think], so the case of
  # CP(Sm, m) == 1  AND  CP(Sm, m-1) < 1 is not possible.

  # DEPRECATED -- TOO SLOW. FASTER CODE BELOW
  # success <- NULL
  # for(i in 1:length(rows.with.cp1))
  # {
  #   for(j in 1:length(columns.of.rows.w.cp1[[i]]))
  #   {
  #     if(mat[rows.with.cp1[i] - 1, columns.of.rows.w.cp1[[i]][j] - 1] < 1) success <- rbind(success, c(rows.with.cp1[i]-1, columns.of.rows.w.cp1[[i]][j]))
  #   }
  # }


  index1 <- 1
  #success <- NULL

  # for(i in rows.with.cp1)
  # {
  #   for(j in columns.of.rows.w.cp1[[index1]])
  #   {
  #     if(j==1) success <- rbind(success, c(i-1, j, final.probs.mat[i, j], final.probs.mat.p0[i, j]))
  #     else
  #     if(mat[i-1, j-1] < 1) success <- rbind(success, c(i-1, j, final.probs.mat[i, j], final.probs.mat.p0[i, j]))
  #   }
  #   index1 <- index1 + 1
  # }

  success <- vector("list", length(rows.with.cp1)*length(columns.of.rows.w.cp1))
  k <- 1
  for(i in rows.with.cp1)
  {
    current.row.probs.mat <- final.probs.mat[i, ]
    current.row.probs.mat.p0 <- final.probs.mat.p0[i, ]
    for(j in columns.of.rows.w.cp1[[index1]])
    {
      if(mat[i-1, j-1] < 1 || j==1){
        success[[k]] <- c(i-1, j, current.row.probs.mat[j], current.row.probs.mat.p0[j])
        k <- k+1
      }
    }
    index1 <- index1 + 1
  }
  success <- do.call(rbind, success)

  colnames(success) <- c("Sm", "m", "prob", "prob.p0")



  #
  # CAN FIND POWER AND TYPE I ERROR AT THIS POINT.
  #
  # IF TOO LOW OR HIGH, STOP AND MOVE ON:
  #


  pwr <- sum(success[,"prob"])
  typeI <- sum(success[,"prob.p0"])

  if(pwr < power | typeI > alpha) # | pwr > power+tol | alpha < alpha-tol)
  {
    return(c(n1, n2, n, r1, r2, typeI, pwr,  NA, NA, thetaF, thetaE, NA))
  }



  # Now failure probabilities. Note that there AT MOST one failure probability in each row, and in that
  # row the failure probability is the one that has the greatest m (i.e. the "furthest right" non-zero entry):

  # First, in the matrix of CPs, find the earliest column containing only zeros and ones (and NAs): This is the latest point ("m") at which the trial stops:
  first.zero.one.col <- which(apply(mat, 2, function(x) sum(!(x %in% c(0,1, NA))))==0)[1]

  # Second, find the the row at which the final zero exists in this column:
  final.zero.in.col <- max(which(mat[, first.zero.one.col]==0))

  # Finally, find the rightmost non-zero probability in each row up and including the one above:

  m.fail <- NULL
  prob.fail <- NULL
  prob.fail.p0 <- NULL

  for(i in 1:final.zero.in.col)
  {
    m.fail[i] <- max(which(final.probs.mat[i ,]!=0))
    prob.fail[i] <- final.probs.mat[i, m.fail[i]]
    prob.fail.p0[i] <- final.probs.mat.p0[i, m.fail[i]]
  }
  Sm.fail <- 0:(final.zero.in.col-1)



  # # Identify non-zero terms in each row:
  # m.fail <- rep(NA, r2+1)
  # prob.fail <- rep(NA, r2+1)
  # prob.fail.p0 <- rep(NA, r2+1)
  #
  #
  #
  # for(i in 1:(r2+1))
  # {
  # m.fail[i] <- max(which(final.probs.mat[i ,]!=0))
  # prob.fail[i] <- final.probs.mat[i, m.fail[i]]
  # prob.fail.p0[i] <- final.probs.mat.p0[i, m.fail[i]]
  # }
  # Sm.fail <- 0:r2


  fail.deets <- cbind(Sm.fail, m.fail, prob.fail, prob.fail.p0)

  all.deets <- rbind(fail.deets, success)
  all.deets <- as.data.frame(all.deets)
  all.deets$success <- c(rep("Fail", length(m.fail)), rep("Success", nrow(success)))
  names(all.deets) <- c("Sm", "m", "prob", "prob.p0", "success")

  output <- all.deets


  ##################### Now find characteristics of design


  #PET.failure <- sum(output$prob[output$m < n & output$Sm <= r2]) # Pr(early termination due to failure)
  #PET.failure.p0 <- sum(output$prob.p0[output$m < n & output$Sm <= r2])

  #PET.success <- sum(output$prob[output$m < n & output$Sm > r2]) # Pr(early termination due to success)
  #PET.success.p0 <- sum(output$prob.p0[output$m < n & output$Sm > r2])


  #PET <- PET.failure + PET.success # Pr(Early termination due to failure or success)
  #PET.p0 <- PET.failure.p0 + PET.success.p0

  #PET.df <- data.frame(PET.failure, PET.success, PET, PET.failure.p0, PET.success.p0, PET.p0)

  #power <- sum(output$prob[output$success=="Success"])

  #output$obsd.p <- output$Sm/output$m

  #output$bias <- output$obsd.p - p


  #bias.mean <- sum(output$bias*output$prob)
  # wtd.mean(output$bias, weights=output$prob, normwt=TRUE) # Check
  #bias.var <- wtd.var(output$bias, weights=output$prob, normwt=TRUE)

  sample.size.expd <- sum(output$m*output$prob)
  # wtd.mean(output$m, weights=output$prob, normwt=TRUE) # Check
  #sample.size <- wtd.quantile(output$m, weights=output$prob, normwt=TRUE, probs=c(0.25, 0.5, 0.75))

  sample.size.expd.p0 <- sum(output$m*output$prob.p0)
  # wtd.mean(output$m, weights=output$prob.p0, normwt=TRUE) # Check
  #sample.size.p0 <- wtd.quantile(output$m, weights=output$prob.p0, normwt=TRUE, probs=c(0.25, 0.5, 0.75))


  #alpha <- sum(output$prob.p0[output$success=="Success"])

  #print(output)


  ############ MARCH 2018: EFFECTIVE N. THE "EFFECTIVE N" OF A STUDY IS THE "REAL" MAXIMUM SAMPLE SIZE
  ######## The point is where every Sm for a given m equals zero or one is necessarily where a trial stops

  cp.colsums <- apply(mat, 2, function(x) { sum(x==0, na.rm=TRUE)+sum(x==1, na.rm=TRUE)} ) # Sum the CP values that equal zero or one in each column
  # Test these against m+1 -- i.e. for i=m, need sum to be m+1:
  effective.n <- min(which(cp.colsums==2:(n+1)))



  #output <- list(output, PET.df, mean.bias=bias.mean, var.bias=bias.var, sample.size=sample.size, expd.sample.size=sample.size.expd,
  #               sample.size.p0=sample.size.p0, expd.sample.size.p0=sample.size.expd.p0, alpha=alpha, power=power)
  to.return <- c(n1, n2, n, r1, r2, typeI, pwr, sample.size.expd.p0, sample.size.expd, thetaF, thetaE, effective.n)
  #to.return <- mat
  to.return
}

#### Part 5: from "SC_bisection_oct2020.R" ####

############## THIS CODE FINDS ALPHA AND POWER FOR ALL COMBNS OF n1/n2/n/r1/r, given p/p0, UNDER NSC



#' findSCDesigns
#'
#' This function finds admissible design realisations for single-arm binary outcome trials, using stochastic curtailment.
#' This function differs from singlearmDesign in that it includes a Simon-style interim analysis after some n1 participants.
#' The output is a data frame of admissible design realisations.
#' @param nmin Minimum permitted sample size.
#' @param nmax Maximum permitted sample size.
#' @param p0 Probability for which to control the type-I error-rate
#' @param p1 Probability for which to control the power
#' @param alpha Significance level
#' @param power Required power (1-beta).
#' @param maxthetaF Maximum value of lower CP threshold theta_F_max. Defaults to p.
#' @param minthetaE Minimum value of upper threshold theta_E_min. Defaults to p.
#' @param bounds choose what final rejection boundaries should be searched over: Those of A'Hern ("ahern"), Wald ("wald") or no constraints (NA). Defaults to "wald".
#' @param max.combns Provide a maximum number of ordered pairs (theta_F, theta_E). Defaults to 1e6.
#' @param maxthetas Provide a maximum number of CP values used to create ordered pairs (theta_F, theta_E). Can be used instead of max.combns. Defaults to NA.
#' @param fixed.r1 Choose what interim rejection boundaries should be searched over. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param fixed.r Choose what final rejection boundaries should be searched over. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param fixed.n1 Choose what interim sample size values n1 should be searched over. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param exact.thetaF Provide an exact value for lower threshold theta_F. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param exact.thetaE Provide an exact value for upper threshold theta_E. Useful for reproducing a particular design realisation. Defaults to NA.
#' @param progressBar Logical. If TRUE, shows progress bar. Defaults to FALSE.
#' @return Output is a list of two dataframes. The first, $input, is a one-row
#' data frame that contains important arguments used in the call. The second,
#' $all.des,contains the operating characteristics of all admissible designs found.
#' @export
#' @examples \donttest{findSCdesigns(nmin = 20, nmax = 20, p0 = 0.1, p1 = 0.4, power = 0.8, alpha = 0.1, max.combns=1e2)}
findSCdesigns <- function(nmin,
                 nmax,
                 p0,
                 p1,
                 alpha,
                 power,
                 minthetaE=p1,
                 maxthetaF=p1,
                 bounds="wald",
                 fixed.r1=NA,
                 fixed.r=NA,
                 fixed.n1=NA,
                 max.combns=1e6,
                 maxthetas=NA,
                 exact.thetaF=NA,
                 exact.thetaE=NA,
                 progressBar=FALSE){
  grp <- V2 <- NULL
  #sc.all <- findN1N2R1R2_df(nmin=nmin, nmax=nmax, e1=FALSE)
  sc.all <- findSimonN1N2R1R2(nmin=nmin, nmax=nmax, e1=FALSE)

  if(is.numeric(bounds)){
    sc.subset <- sc.all[sc.all$r %in% bounds,]
  }else{
    if(bounds=="ahern"){
      ### AHERN'S BOUNDS: ###
      sc.subset <- sc.all[sc.all$r>=p0*sc.all$n & sc.all$r<=p1*sc.all$n, ]
    }else{
      if(bounds=="wald"){
        # Even better to incorporate Wald's bounds:
        nposs <- nmin:nmax
        denom <- log(p1/p0) - log((1-p1)/(1-p0))
        accept.null <- log((1-power)/(1-alpha)) / denom  + nposs * log((1-p0)/(1-p1))/denom
        accept.null <- floor(accept.null)
        reject.null <- log((power)/alpha) / denom  + nposs * log((1-p0)/(1-p1))/denom
        reject.null <- ceiling(reject.null)
        r.wald <- NULL
        ns.wald <- NULL
        for(i in 1:length(nposs)){
          r.wald <- c(r.wald, accept.null[i]:reject.null[i])
          ns.wald <- c(ns.wald, rep(nposs[i], length(accept.null[i]:reject.null[i])))
        }
        #sc.subset <- data.frame(n=ns.wald, r=r.wald)
        keep.index <- NULL
        for(i in 1:nrow(sc.all)){
          keep.index[i] <- sum(sc.all$n[i]==ns.wald & sc.all$r[i]==r.wald) > 0
        }
        sc.subset <- sc.all[keep.index,]
      }
    }
  }


  if(!is.na(fixed.r1)){
    sc.subset <- sc.subset[sc.subset$r1 %in% fixed.r1,]
  }
  if(!is.na(fixed.r)){
    sc.subset <- sc.subset[sc.subset$r %in% fixed.r,]
  }
  if(!is.na(fixed.n1)){
    sc.subset <- sc.subset[sc.subset$n1 %in% fixed.n1,]
  }

  # For each design, find the potential values of theta:
  n.design.rows <- nrow(sc.subset)
  store.all.thetas <- vector("list", n.design.rows)

  for(i in 1:n.design.rows)
  {
    store.all.thetas[[i]] <- findThetas(n1=sc.subset[i,1], n2=sc.subset[i,2], n=sc.subset[i,3], r1=sc.subset[i,4], r2=sc.subset[i,5], p0=p0, p=p1)
  }
  # Takes seconds




  # Much of the required computing is independent of thetaF/1, and therefore only needs to be computed once per design.
  # The function
  #
  #       outsideTheta
  #
  # undertakes this computation, and stores all the necessary objects for obtaining power, type I error, ESS, etc. for later
  # use in the function
  #
  #       insideTheta
  #


  outside.theta.data <- vector("list", nrow(sc.subset))
  all.theta.combns <- vector("list", nrow(sc.subset))

if(!is.na(max.combns) & is.na(maxthetas)){ # Only use max.combns if maxthetas is not specified.
    maxthetas <- sqrt(2*max.combns)
}

  if(is.na(exact.thetaF)){ # if exact theta values have NOT been supplied:
  ##### To cut down on computation, try cutting down the number of thetas used.
  ##### If a design has > maxthetas theta values (CP values), cut this down:
    for(i in 1:nrow(sc.subset)){
      current.theta.vec <- store.all.thetas[[i]]
      store.all.thetas[[i]] <- current.theta.vec[current.theta.vec >= minthetaE | current.theta.vec <= maxthetaF] # Subset thetas BEFORE thinning
      while(length(store.all.thetas[[i]]) > maxthetas){
        every.other.element <- rep(c(FALSE, TRUE), 0.5*(length(store.all.thetas[[i]])-2))
        store.all.thetas[[i]] <-  store.all.thetas[[i]][c(TRUE, every.other.element, TRUE)]
      }
    }
    all.theta.combns <- lapply(store.all.thetas, function(x) {t(combn(x, 2)) })
  for(i in 1:nrow(sc.subset)){
    outside.theta.data[[i]] <- outsideTheta(n1=sc.subset[i,1], n2=sc.subset[i,2], r1=sc.subset[i,4], r2=sc.subset[i,5], p0=p0, p=p1)
    all.theta.combns[[i]] <- t(combn(store.all.thetas[[i]], 2))
    all.theta.combns[[i]] <- all.theta.combns[[i]][all.theta.combns[[i]][,2] >= minthetaE & all.theta.combns[[i]][,1] <= maxthetaF, ] # Reduce number of combinations
    if(length(all.theta.combns[[i]])==2) all.theta.combns[[i]] <- matrix(c(0, 0.999, 0, 1), nrow=2, byrow=T) # To avoid a crash. See (*) below
    all.theta.combns[[i]] <- all.theta.combns[[i]][order(all.theta.combns[[i]][,2], decreasing=T),]  ### XXX
  }
  # (*) If the only thetas for a design are 0 and 1, then some code below will crash. This
  # could be elegantly resolved by adding an if or ifelse statement, but probably at some
  # computational expense. Will deal with these exceptions here, by adding a "false" value
  # for theta.
  }else{ # if exact theta values have been supplied:
    for(i in 1:length(store.all.thetas)){
      keep <- abs(store.all.thetas[[i]]-exact.thetaF)<1e-2 | abs(store.all.thetas[[i]]-exact.thetaE)<1e-2
      store.all.thetas[[i]] <- store.all.thetas[[i]][keep]
    }
    for(i in 1:nrow(sc.subset)){
      outside.theta.data[[i]] <- outsideTheta(n1=sc.subset[i,1], n2=sc.subset[i,2], r1=sc.subset[i,4], r2=sc.subset[i,5], p0=p0, p=p1)
      all.theta.combns[[i]] <- t(combn(store.all.thetas[[i]], 2))
      all.theta.combns[[i]] <- all.theta.combns[[i]][order(all.theta.combns[[i]][,2], decreasing=T),]  ### XXX
    }
  }







  # NOTE: thetaE is fixed while thetaF increases. As thetaF increases, power decreases.

  if(progressBar==TRUE) pb <- txtProgressBar(min = 0, max = length(outside.theta.data), style = 3)
  h.results.list <- vector("list", length(outside.theta.data)) #
  y <- NULL

  for(h in 1:length(outside.theta.data)) # For every possible combn of n1/n2/n/r1/r:
  {
    k <- 1
    h.results <- vector("list", nrow(all.theta.combns[[h]]))
    current.theta.combns <- all.theta.combns[[h]] # Take just the thetaF/1 combns for that design.
    # current.theta.combns <- current.theta.combns[current.theta.combns[,2]>0.7, ] # remove combns where thetaE < 0.7 (as not justifiable)

    # Don't need a matrix of all thetaF and thetaE combns -- potentially quicker to have thetaF as a vector (already defined), and the list can be a list of thetaE vectors:
    current.theta.combns <- data.table::data.table(current.theta.combns)
    current.theta.combns[, grp := .GRP, by = V2]
    data.table::setkey(current.theta.combns, grp)
    split.thetas <- current.theta.combns[, list(list(.SD)), by = grp]$V1
    thetaF.list <- lapply(split.thetas, function(x) x[, 1])
    all.thetas <- rev(store.all.thetas[[h]])[-length(store.all.thetas[[h]])] # All thetaE values, decreasing, not including the final value, thetaE=0.
    all.thetas <- all.thetas[all.thetas>=minthetaE]

    # Of course we need the other characteristics and so on of the design, which are stored in the list "outside.theta.data":
    y <- outside.theta.data[[h]]
    eff.n <- y$n

    for(i in 1:length(all.thetas)) # For each thetaE,
    {
      thetaFs.current <- thetaF.list[[i]] # thetaF.list is a list of i vectors. The i'th vector in the list contains all possible values of thetaF for the i'th thetaE (stored as all.thetas[i])
      rows <- nrow(thetaFs.current)

      # Begin bisection method: Use the bisection method to find where to "start":
      a <- 1
      b <- nrow(thetaFs.current)
      d <- ceiling((b-a)/2)

      while((b-a)>1)
      {
        output <- insideTheta(n1=y$n1, n2=y$n2, n=y$n, r1=y$r1, r2=y$r2, thetaF=as.numeric(thetaFs.current[d]), thetaE=all.thetas[i], cp.sm=y$cp.sm, cp.m=y$cp.m, p0=p0, p=p1,
                              mat=y$mat, coeffs.p0=y$coeffs.p0, coeffs=y$coeffs, power=power, alpha=alpha)
        if(output[6] <= alpha) { # type I error decreases as index (and thetaF) increases. Jump backwards if type I error is smaller than alpha, o/w forwards.
          b <- d } else {
            a <- d
          }
        #  print(paste("a is ", a, "... b is ", b, "... d is ", d, sep=""))
        d <- a + floor((b-a)/2)
      }

      # Take care of "edge case" where feasible design exists in the first row:
      if(a==1) # (and also by necessity, b==2)
      {
        output <- insideTheta(n1=y$n1, n2=y$n2, n=y$n, r1=y$r1, r2=y$r2, thetaF=as.numeric(thetaFs.current[1]), thetaE=all.thetas[i], cp.sm=y$cp.sm, cp.m=y$cp.m, p0=p0, p=p1,
                              mat=y$mat, coeffs.p0=y$coeffs.p0, coeffs=y$coeffs, power=power, alpha=alpha)
        if(output[6] <= alpha) b <- 1 # Should we start at the first row or the second?
      }

      # We can now proceed, moving sequentially from index==b (or cause a break if we wish):

      first.result <-  insideTheta(n1=y$n1, n2=y$n2, n=y$n, r1=y$r1, r2=y$r2, thetaF=as.numeric(thetaFs.current[b]), thetaE=all.thetas[i], cp.sm=y$cp.sm, cp.m=y$cp.m, p0=p0, p=p1,
                                   mat=y$mat, coeffs.p0=y$coeffs.p0, coeffs=y$coeffs, power=power, alpha=alpha)
      # If the type I error equals zero, then the trial is guaranteed to fail at m=1. Further, this will also be the case for all subsequent thetaF values
      # (as thetaF is increasing). Further still, there are no feasible designs for subsequent thetaE values, and so we can skip to the next n/r combination
      # If type I error doesn't equal zero, then feasible designs exist and they should be recorded. Note: type I error equals 2 means success guaranteed
      # at m=1, and the conclusion is the same.
      if(first.result[6]!=0 | first.result[6]!=2)
      {
        pwr <- first.result[7]
        while(pwr >= power & b <= rows) # Keep going until power drops below 1-beta, i.e. no more feasible designs, or we reach the end of the data frame.
        {
          h.results[[k]] <-  insideTheta(n1=y$n1, n2=y$n2, n=y$n, r1=y$r1, r2=y$r2, thetaF=as.numeric(thetaFs.current[b]), thetaE=all.thetas[i], cp.sm=y$cp.sm, cp.m=y$cp.m, p0=p0, p=p1,
                                         mat=y$mat, coeffs.p0=y$coeffs.p0, coeffs=y$coeffs, power=power, alpha=alpha)
          pwr <- h.results[[k]][7] #### Added 3rd July -- should have been there all the time ¬_¬
          k <- k+1
          b <- b+1
        }
      } else { # if first.result[3]==0 or 2, i.e., there are no feasible designs for this thetaE (and hence no remaining feasible designs for this n/r combn), break:
        break
      }

    } # end of "i" loop

    if(progressBar==TRUE) setTxtProgressBar(pb, h)

    h.results.df <- do.call(rbind.data.frame, h.results)

    # Before moving on to the next set of values, avoid creating a huge object by cutting out all dominated and duplicated designs for this combn of n1/n2/n/r1/r:
    if(nrow(h.results.df)>0)
    {
      names(h.results.df) <- c("n1", "n2", "n", "r1", "r2", "alpha", "power", "EssH0", "Ess", "thetaF", "thetaE", "eff.n")

      # Remove all "skipped" results:
      h.results.df <- h.results.df[!is.na(h.results.df$Ess),]
      # Remove dominated and duplicated designs:
      discard <- rep(NA, nrow(h.results.df))
      for(i in 1:nrow(h.results.df))
      {
        discard[i] <- sum(h.results.df$EssH0[i] > h.results.df$EssH0 & h.results.df$Ess[i] > h.results.df$Ess & h.results.df$n[i] >= h.results.df$n)
      }
      h.results.df <- h.results.df[discard==0,]

      duplicates <- duplicated(h.results.df[, c("n", "Ess", "EssH0")])
      h.results.df <- h.results.df[!duplicates,]

      h.results.list[[h]] <- h.results.df
    }

  }

  # All designs are combined in the dataframe "results".
  results <- do.call(rbind.data.frame, h.results.list)

  # Remove all "skipped" results:
  results <- results[!is.na(results$Ess),]

  # Remove dominated and duplicated designs:
  discard <- rep(NA, nrow(results))
  for(i in 1:nrow(results))
  {
    discard[i] <- sum(results$EssH0[i] > results$EssH0 & results$Ess[i] > results$Ess & results$n[i] >= results$n)
  }
  subset.results <- results[discard==0,]

  # Remove duplicates:
  duplicates <- duplicated(subset.results[, c("n", "Ess", "EssH0")])
  subset.results <- subset.results[!duplicates,]
  sc.input <- data.frame(nmin=nmin, nmax=nmax, p0=p0, p1=p1, alpha=alpha, power=power)
  sc.output <- list(input=sc.input,
                       all.des=subset.results)
  return(sc.output)
}

#### Plotting ####
#### Plotting omni-admissible designs (Fig 6 in manuscript):

design.plot <- function(results.df, loss.df, scenario, design)
{
  qs <- q0 <- q1 <- NULL
  index <- which(results.df$design==design)
  all.loss.subset <- t(loss.df[index,])
  designs <- results.df[index,]
  designs.vec <- as.character(as.vector(apply(all.loss.subset, 1, which.min)))
  #code <- c("1"=design.details[1], "2"=design.details[2], "3"=design.details[3])
  #simon.designs.scen1.vec <- code[simon.designs.scen1.vec]
  designs.used <- sort(as.numeric(unique(designs.vec)))

  if(design=="simon" | design=="nsc")  design.details <- apply(results.df[index[designs.used],], 1, function(x) paste("{", x["r1"], "/", x["n1"], ", ", x["r"], "/", x["n"], "}", sep=""))
  if(design=="simon_e1") design.details <- apply(results.df[index[designs.used],], 1, function(x) paste("{(", x["r1"], " ", x["e1"], ")/", x["n1"], ", ", x["r"], "/", x["n"], "}", sep=""))
  if(design=="sc") design.details <- apply(results.df[index[designs.used],], 1, function(x) paste("{", x["r1"], "/", x["n1"], ", ", x["r"], "/", x["n"], ", ", format(round(as.numeric(x["thetaF"]),2), nsmall=2), "/", format(round(as.numeric(x["thetaE"]),2), nsmall=2), "}", sep=""))
  if(design=="mstage")  design.details <- apply(results.df[index[designs.used],], 1, function(x) paste("{", x["r"], "/", x["n"], ", ", format(round(as.numeric(x["thetaF"]),2), nsmall=2), "/", format(round(as.numeric(x["thetaE"]),2), nsmall=2), "}", sep=""))

  designs.df <- data.frame(design=designs.vec, q0=qs[,"w0"], q1=qs[, "w1"])
  designs.df$design <- factor(designs.df$design, levels = as.character(sort(as.numeric(levels(designs.df$design)))))
  design.type <- switch(design, "simon"="Simon", "simon_e1"="Mander Thompson", "nsc"="NSC", "sc"="SC", "mstage"="m-stage")

  output <-  ggplot2::ggplot(designs.df,  ggplot2::aes(x=q1, y=q0)) +
     ggplot2::geom_raster(ggplot2::aes(fill = design)) +
     ggplot2::scale_fill_discrete(name="Design",labels=as.vector(design.details)) +
     ggplot2::ggtitle(paste(design.type, ", scenario ", scenario, sep="")) +
     ggplot2::xlab(expression(w[1])) +
     ggplot2::ylab(expression(w[0])) +
     ggplot2::theme_bw()+
     ggplot2::theme(#axis.text.x=element_text(),
      legend.text= ggplot2::element_text(size=7),
      plot.title =  ggplot2::element_text(size =  ggplot2::rel(1.4)),
      axis.ticks= ggplot2::element_blank(),
      axis.line= ggplot2::element_blank(),
      axis.title= ggplot2::element_text(size=rel(1.5)),
      axis.title.y= ggplot2::element_text(angle = 0, vjust=0.5),
      axis.text= ggplot2::element_text(size=rel(1.5)),
      legend.key.size =  ggplot2::unit(0.4, "cm"),
      legend.position = c(0.8,0.75),
      panel.border= ggplot2::element_blank(),
      panel.grid.major= ggplot2::element_line(color='#eeeeee'))
  return(output)
}

plot.all <- function(results, loss, scen){
  design.types <- c("simon", "simon_e1", "nsc", "sc", "mstage")
  plots.list <- vector("list", 5)
  for(i in 1:5){
    plots.list[[i]] <- design.plot(results.df=results, loss.df=loss, scenario=scen, design=design.types[i])
  }
  plots.list
}




# Expected loss:

plotExpdLoss <- function(loss.df, design.type, scenario, design.loss.range){
  title <- paste("Expected loss: ", design.type, ", Scenario ", scenario, sep="")
  ggplot2::ggplot(loss.df, ggplot2::aes_string("w1", "w0")) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw() +
    ggplot2::xlab(expression(w[1])) +
    ggplot2::ylab(expression(w[0])) +
    ggplot2::geom_raster(ggplot2::aes_string(fill="loss")) +
    ggplot2::scale_fill_gradient(low="red", high="darkblue", limits=design.loss.range) +
    ggplot2::theme(#axis.text.x=element_text(),
      axis.text=ggplot2::element_text(size=15),
      axis.title.x = ggplot2::element_text(size=15),
      axis.title.y = ggplot2::element_text(size=15, angle = 0, vjust = 0.5),
      legend.text=ggplot2::element_text(size=12),
      plot.title = ggplot2::element_text(size = rel(1.5)),
      axis.ticks=ggplot2::element_blank(),
      axis.line=ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(0.75, "cm"),
      legend.position = c(0.9,0.65),
      legend.title = ggplot2::element_text(size=15),
      panel.border=ggplot2::element_blank(),
      panel.grid.major=ggplot2::element_line(color='#eeeeee'))
}

plotAllExpdLoss <- function(loss.list, scen, loss.range){
  design.types <- c("Simon", "MT", "NSC", "SC", "m-stage")
  plots.list <- vector("list", 5)
  for(i in 1:5){
    plots.list[[i]] <- plotExpdLoss(loss.df=loss.list[[i]],
                                    design.type = design.types[i],
                                    scenario=scen,
                                    design.loss.range=loss.range)
  }
  plots.list
}

# Discard dominated designs:
rmDominatedDesigns <- function(data, ess0="EssH0", ess1="Ess", n="n"){
  discard <- rep(NA, nrow(data))
  for(i in 1:nrow(data)){
    discard[i] <- sum(data[i, ess0] > data[, ess0] & data[i, ess1] > data[, ess1] & data[i, n] > data[, n])
    }
  admissible.designs <- data[discard==0, ]
  admissible.designs
}


##### Plotting quantiles of ESS #####
cohortFindMatrix <- function(n, r, Csize, p0, p1){
  q1 <- 1-p1
  mat <- matrix(3, ncol=n, nrow=min(r+Csize, n)+1) #  nrow=r+Csize+1 unless number of stages equals 1.
  rownames(mat) <- 0:(nrow(mat)-1)
  mat[(r+2):nrow(mat),] <- 1
  mat[1:(r+1),n] <- 0

  pat.cols <- seq(n, 1, by=-Csize)[-1]

  for(i in (r+1):1){
    for(j in pat.cols){  # Only look every C patients (no need to look at final col)
      if((r+1-(i-1) > n-j)) { # IF success is not possible [Total responses needed (r+1) - responses so far (i-1)] > [no. of patients remaining (n-j)], THEN
        mat[i,j] <- 0
      }else{
        if(i-1<=j){ # Condition: Sm<=m
          newcp <- sum(choose(Csize,0:Csize)*p1^(0:Csize)*q1^(Csize:0)*mat[i:(i+Csize),j+Csize])
          # NOTE: No rounding of CP up to 1 or down to 0: This is the basic matrix to which many pairs of (thetaF, thetaE) will be applied in the function cohortFindDesign(). Lines below can be deleted.
          # if(newcp > thetaE) mat[i,j] <- 1
          # if(newcp < thetaF) mat[i,j] <- 0
          # if(newcp <= thetaE & newcp >= thetaF) mat[i,j] <- newcp
          mat[i,j] <- newcp
        }
      }
    }
  }


  for(i in 3:nrow(mat)) {
    mat[i, 1:(i-2)] <- NA
  }
  return(mat)
}


cohortFindQuantileSS <- function(n, r, Csize, thetaF, thetaE, p0, p1, coarseness=0.01, lower=0.1, upper=0.9, return.only.quantiles=TRUE){

  mat <- cohortFindMatrix(n = n, r=r, Csize=Csize, p0=p0, p1=p1)

  q1 <- 1-p1
  q0 <- 1-p0

  # Note: Only need the first r+Csize+1 rows -- no other rows can be reached.
  coeffs <- matrix(NA, ncol=n, nrow=r+Csize+1)
  coeffs.p0 <- coeffs

  for(i in 1:(r+Csize+1)){
    coeffs[i, ] <- p1^(i-1) * q1^((2-i):(2-i+n-1))
    coeffs.p0[i, ] <- p0^(i-1) * q0^((2-i):(2-i+n-1))
  }


  ##### LIST of matrix of coefficients: Probability of a path leading to each point:
  p.vec <- seq(0, 1, by=coarseness)
  coeffs.list <- vector("list", length=length(p.vec))

  for(k in 1:length(coeffs.list)){
    coeffs.list[[k]] <- matrix(NA, ncol=n, nrow=r+Csize+1)
    for(i in 1:(r+Csize+1)){
      coeffs.list[[k]][i, ] <- p.vec[k]^(i-1) * (1-p.vec[k])^((2-i):(2-i+n-1))
    }
  }


  pat.cols <- seq(n, 1, by=-Csize)[-1]

  # Amend CP matrix, rounding up to 1 when CP>theta_1 and rounding down to 0 when CP<thetaF:

  for(i in (r+1):1){
    for(j in pat.cols){  # Only look every Csize patients (no need to look at final col)
      if(i-1<=j){ # Condition: Sm<=m
        newcp <- sum(choose(Csize,0:Csize)*p1^(0:Csize)*q1^(Csize:0)*mat[i:(i+Csize),j+Csize])
        if(newcp > thetaE) mat[i,j] <- 1
        if(newcp < thetaF) mat[i,j] <- 0
        if(newcp <= thetaE & newcp >= thetaF) mat[i,j] <- newcp
      }
    }
  }



  ############# Number of paths to each point:

  pascal.list <- list(c(1,1))
  if(Csize>1){
    for(i in 2:Csize){
      pascal.list[[i]] <- c(0, pascal.list[[i-1]]) + c(pascal.list[[i-1]], 0)
    }
  }

  for(i in (Csize+1):n){
    if(i %% Csize == 1 | Csize==1){
      column <- as.numeric(mat[!is.na(mat[,i-1]), i-1])
      newnew <- pascal.list[[i-1]]
      CPzero.or.one <- which(column==0 | column==1)
      newnew[CPzero.or.one] <- 0
      pascal.list[[i]] <- c(0, newnew) + c(newnew, 0)
    } else {
      pascal.list[[i]] <- c(0, pascal.list[[i-1]]) + c(pascal.list[[i-1]], 0)
    }
  }

  for(i in 1:length(pascal.list)){
    pascal.list[[i]] <- c(pascal.list[[i]], rep(0, length(pascal.list)+1-length(pascal.list[[i]])))
  }


  pascal.t <- t(do.call(rbind, args = pascal.list))
  pascal.t <- pascal.t[1:(r+Csize+1), ]

  # Multiply the two matrices (A and p^b * q^c). This gives the probability of reaching each point:
  # Note: Only need the first r+Csize+1 rows -- no other rows can be reached.
  pascal.t <- pascal.t[1:(r+Csize+1), ]

  final.probs.mat <- pascal.t*coeffs
  final.probs.mat.p0 <- pascal.t*coeffs.p0

  final.probs.mat.list <- vector("list", length = length(p.vec))
  for(i in 1:length(p.vec)){
    final.probs.mat.list[[i]] <- pascal.t*coeffs.list[[i]]
  }


  ##### Now find the terminal points: m must be a multiple of Csize, and CP must be 1 or 0:
  ##### SUCCESS:
  ## Only loop over rows that contain CP=1:
  rows.with.cp1 <- which(apply(mat, 1, function(x) {any(x==1, na.rm = T)}))

  m.success <- NULL
  Sm.success <- NULL
  prob.success <- NULL
  prob.success.p0 <- NULL
  prob.success.list <- vector("list", length(p.vec))

  interims <- seq(Csize, n, by=Csize)

  for(i in rows.with.cp1){
    for(j in interims[interims >= i-1]){ # condition ensures Sm >= m
      if(mat[i,j]==1){
        m.success <- c(m.success, j)
        Sm.success <- c(Sm.success, i-1)
        prob.success <- c(prob.success, final.probs.mat[i, j])
        prob.success.p0 <- c(prob.success.p0, final.probs.mat.p0[i, j])
        for(k in 1:length(p.vec)){
          prob.success.list[[k]] <- c(prob.success.list[[k]], final.probs.mat.list[[k]][i,j])
        }
      }
    }
  }

  prob.success.mat <- do.call(cbind, prob.success.list)
  success <- data.frame(Sm=Sm.success, m=m.success, prob=prob.success, prob.p0=prob.success.p0, success=rep("Success", length(m.success)), prob.success.mat)
  names(success)[6:ncol(success)]


  pwr <- sum(success[,"prob"])
  typeI <- sum(success[,"prob.p0"])


  ##### FAILURE:
  first.zero.in.row <- apply(mat[1:(r+1),], 1, function(x) {which.max(x[interims]==0)})
  m.fail <- Csize * as.numeric(first.zero.in.row)
  Sm.fail <- as.numeric(names(first.zero.in.row))

  fail.deets <- data.frame(Sm=Sm.fail, m=m.fail)
  fail.deets$prob <- apply(fail.deets, 1, function(x) {final.probs.mat[x["Sm"]+1, x["m"]]})
  fail.deets$prob.p0 <- apply(fail.deets, 1, function(x) {final.probs.mat.p0[x["Sm"]+1, x["m"]]})

  prob.fail.list <- vector("list", length(p.vec))
  for(i in 1:length(p.vec)){
    prob.fail.list[[i]] <- apply(fail.deets, 1, function(x) {final.probs.mat.list[[i]][x["Sm"]+1, x["m"]]})
  }
  prob.fail.mat <- do.call(cbind, prob.fail.list)

  fail.deets$success <- rep("Fail", nrow(fail.deets))

  fail.deets <- cbind(fail.deets, prob.fail.mat)

  names(fail.deets) <- names(success)
  output <- rbind(fail.deets, success)

  ###### ADDED

  ordered.outp <- output[order(output$m),]
  cum.probs.mat <- apply(ordered.outp[,c(6:ncol(ordered.outp))], 2, cumsum)

  index.median <- apply(cum.probs.mat, 2, function(x) which.max(x>=0.5))
  actual.median <- ordered.outp$m[index.median]

  index.lower.percentile <- apply(cum.probs.mat, 2, function(x) which.max(x>=lower))
  actual.lower.percentile <- ordered.outp$m[index.lower.percentile]

  index.upper.percentile <- apply(cum.probs.mat, 2, function(x) which.max(x>=upper))
  actual.upper.percentile <- ordered.outp$m[index.upper.percentile]

  ###### END OF ADDED


  sample.size.expd <- sum(output$m*output$prob)
  sample.size.expd.p0 <- sum(output$m*output$prob.p0)


  ############ EFFECTIVE N. THE "EFFECTIVE N" OF A STUDY IS THE "REAL" MAXIMUM SAMPLE SIZE
  ######## The point is where every Sm for a given m equals zero or one is necessarily where a trial stops

  cp.colsums <- apply(mat, 2, function(x) { sum(x==0, na.rm=TRUE)+sum(x==1, na.rm=TRUE)} ) # Sum the CP values that equal zero or one in each column
  possible.cps <- apply(mat, 2, function(x) {sum(!is.na(x))})

  effective.n <- min(which(cp.colsums==possible.cps))
  if(return.only.quantiles==FALSE){
    to.return <- list(row=c(n, r, Csize, typeI, pwr, sample.size.expd.p0, sample.size.expd, thetaF, thetaE, effective.n),
                      quantiles=data.frame(p.vec=p.vec, median=actual.median, lower=actual.lower.percentile, upper=actual.upper.percentile, design="m-stage", cohort.size=Csize, stages=n/Csize))
  } else {
    to.return <- data.frame(p.vec=p.vec, median=actual.median, lower=actual.lower.percentile, upper=actual.upper.percentile, design="m-stage")
    to.return$cohort.size <- Csize
    to.return$stages <- n/to.return$cohort.size
  }
  return(to.return)
}


simonQuantiles <- function(r1, n1, n2, coarseness=0.1, lower=0.1, upper=0.9){
  p.vec <- seq(from=0, to=1, by=coarseness)
  pet <- pbinom(r1, n1, p.vec)
  actual.median <- rep(NA, length(p.vec))
  actual.median[pet > 0.5] <- n1
  actual.median[pet == 0.5] <- n1 + 0.5*n2
  actual.median[pet < 0.5] <- n1 + n2
  actual.10.percentile <- rep(NA, length(p.vec))
  actual.10.percentile[pet > lower] <- n1
  actual.10.percentile[pet == lower] <- n1 + 0.5*n2
  actual.10.percentile[pet < lower] <- n1 + n2
  actual.90.percentile <- rep(NA, length(p.vec))
  actual.90.percentile[pet > upper] <- n1
  actual.90.percentile[pet == upper] <- n1 + 0.5*n2
  actual.90.percentile[pet < upper] <- n1 + n2
  output <- data.frame(p.vec=p.vec, median=actual.median, lower=actual.10.percentile, upper=actual.90.percentile,
                       design="simon", cohort.size=0.5*(n1+n2), stages=2)
  return(output)
}


plotMedianSSDesign <- function(dataf, p0, p1, factor, title){
  required.columns <- c("p.vec", "median", "lower", "upper")
  if(!all(required.columns %in% names(dataf))){
    stop(cat("Data frame must have the following columns: ", required.columns))
  }
  dataf[, factor] <- factor(dataf[, factor])

  ggplot2::ggplot(dataf, ggplot2::aes_string(x="p.vec", y="median", ymin="lower", ymax="upper"), ggplot2::aes_string(group="factor")) +
    ggplot2::geom_line(ggplot2::aes_string(color="factor"), size=1) +
    ggplot2::geom_ribbon(ggplot2::aes_string(fill="factor"), alpha = 0.3) +
    ggplot2::labs(title=title,
         x ="p", y = "Sample size", fill="Design\n(block size B)", color="Design\n(block size B)")+
    ggplot2::geom_vline(xintercept = p0, size=0.5, color="grey60", linetype="longdash") +
    ggplot2::geom_vline(xintercept = p1, size=0.5, color="grey60", linetype="longdash") +
    ggplot2::theme_bw()+
    ggplot2::theme(legend.text=ggplot2::element_text(size=rel(1.5)),
          legend.title=ggplot2::element_text(size=rel(1.5)),
          plot.title = ggplot2::element_text(size = rel(1.5)),
          axis.ticks=ggplot2::element_blank(),
          axis.line=ggplot2::element_blank(),
          axis.text.x=ggplot2::element_text(angle=0),
          axis.title=ggplot2::element_text(size=rel(1.5)),
          axis.text=ggplot2::element_text(size=rel(1.5))
    )
}

# Compare boundaries to Wald/A'Hern ####
# Find the boundaries:
findWaldAhernBounds <- function(N.range, beta, alpha, p0, p1, round=F){
  findWaldBounds <- function(N.vec, bet, alph, p0., p1., rounding){
    denom <- log(p1./p0.) - log((1-p1.)/(1-p0.))
    accept.null <- log((1-bet)/alph) / denom  +  N.vec * log((1-p0.)/(1-p1.))/denom
    reject.null <- log(bet/(1-alph)) / denom  +  N.vec * log((1-p0.)/(1-p1.))/denom
    if(rounding==TRUE){
      accept.null <- ceiling(accept.null)
      reject.null <- floor(reject.null)
    }
    output <- data.frame(N=N.vec, lower=reject.null, upper=accept.null, design="Wald")
    output
  }
  findAhernBounds <- function(N.vec, p0., p1., rounding){
    ahern.lower <- N.vec*p0.
    ahern.upper <- N.vec*p1.
    if(rounding==TRUE){
      ahern.lower <- floor(ahern.lower)
      ahern.upper <- ceiling(ahern.upper)
    }
    output <- data.frame(N=N.vec, lower=ahern.lower, upper=ahern.upper, design="A'Hern")
    output
  }

  wald.output <- findWaldBounds(N.vec=N.range, bet=beta, alph=alpha, p0.=p0, p1.=p1, rounding=round)
  ahern.output <- findAhernBounds(N.vec=N.range, p0.=p0, p1.=p1, rounding=round)
  output <- rbind(ahern.output, wald.output)
  output
}

# Plot the boundaries along with the m-stage boundaries (two different data frames):
plotBounds <- function(wald.ahern.bounds, mstage.bounds, title="XXX"){
  print(range(wald.ahern.bounds$N))
  plot.output <- ggplot2::ggplot()+
    ggplot2::geom_ribbon(data=wald.ahern.bounds,
                         ggplot2::aes(x=wald.ahern.bounds$N,
                                      y=NULL,
                                      ymin=wald.ahern.bounds$lower,
                                      ymax=wald.ahern.bounds$upper,
                                      fill=wald.ahern.bounds$design),
                         alpha = 0.3)+
    ggplot2::geom_point(data=mstage.bounds,
                        ggplot2::aes(x=mstage.bounds$n,
                                     y=mstage.bounds$r))+
    ggplot2::labs(title=title,
                  x ="Maximum sample size",
                  y = "Final rejection boundary",
                  fill="Design")+
    ggplot2::theme_bw()+
    ggplot2::scale_x_continuous(limits=range(wald.ahern.bounds$N), expand=c(0,1))
  plot.output
}



# Constraining thetaF and thetaE: ####
countCPsSingleStageNSC <- function(r,N){
  total.CPs <- (r+1)*(N-r)-1
  total.CPs
}

countCPsTwoStageNSC <- function(r1,n1,r,N){
  total.CPs <- (r1+1)*(n1-r1) + (r-r1)*(N-r) - 1
  total.CPs
}

countOrderedPairsSimple <- function(total.CPs, SI.units=FALSE){
  total.ordered.pairs <- t(combn(c(1:total.CPs), 2))
  if(SI.units)  total.ordered.pairs <- format(total.ordered.pairs, scientific=TRUE)
  total.ordered.pairs
}

countOrderedPairsComplex <- function(r, n, p0, p1, thetaFmax=1, thetaEmin=0){
  CPmat <- findCPmatrix(r=r, n=n, Csize=1, p0=p0, p1=p1)
  CP.vec <- sort(unique(c(CPmat)))
  CP.vec <- CP.vec[CP.vec<=thetaFmax | CP.vec>=thetaEmin]
  ordered.pairs <- t(combn(CP.vec, 2))
  ordered.pairs <- ordered.pairs[ordered.pairs[,1]<=thetaFmax & ordered.pairs[,2]>=thetaEmin, ]
  output <- list(CPs=length(CP.vec),
                 OPs=nrow(ordered.pairs)
                 )
  return(output)
  }

# findTotalNoPairs <- function(nmin, nmax, p0, p1, ahern.plus=FALSE, SI.units=FALSE){
#   n.vec <- nmin:nmax
#   OP.list <- vector("list", length=length(n.vec))
#   for(i in 1:length(n.vec)){
#     rmin <- floor(N*p0)
#     rmax <- ceiling(N*p1)
#     if(ahern.plus==TRUE) rmax <- rmax + 1
#     CPs <- countCPsSingleStageNSC(r=rmin:rmax, N=n.vec[i])
#     OP.list[[i]] <- countOrderedPairs(CPs)
#   }
#   total.OPs <- sum(unlist(OP.list))
#   if(SI.units)  total.OPs <- format(total.OPs, scientific=TRUE)
#   total.OPs
# }



# Loss function: find individual components ####
# Data frames required: all.results for scen 1 only, all.scen2.results for scen 2, all.scen3.results for scen 3.
# These data frames can be found in the file "scen123_allresults.RData".
# Input: w0, w1, full dataframe (all.results).
# Output: the admissible design for each design type, with expected loss and each loss component.
# findLossComponents <- function(df, w0=1, w1=0){
#   df$loss.ESS0 <- w0*df$EssH0
#   df$loss.ESS1 <- w1*df$Ess
#   df$loss.N <- (1-w0-w1)*df$n
#   df$loss.total <- df$loss.ESS0 + df$loss.ESS1 + df$loss.N
#   df$loss.diff <- df$loss.total - min(df$loss.total)
#   admiss.des <- by(data=df, INDICES = df$design, FUN = function(x) x[x$loss.total==min(x$loss.total), ])
#   admiss.des <- do.call(rbind, admiss.des)
#   admiss.des <- admiss.des[order(match(admiss.des$design, des.order)), ]
#   rownames(admiss.des) <- c("Simon", "MT", "NSC", "SC", "mstage")
#   output <- admiss.des[, c("EssH0", "loss.ESS0", "Ess", "loss.ESS1", "n", "loss.N", "loss.total", "loss.diff")]
#   output
# }



# Simon's design: No curtailment -- only stopping is at end of S1:
simonEfficacy <- function(n1, n2, r1, r, e1, p0, p1)
{

  n <- n1+n2

  # Create Pascal's triangle for S1: these are the coefficients (before curtailment) A, where A * p^b * q*c
  pascal.list.s1 <- list(1)
  for (i in 2:(n1+1)) pascal.list.s1[[i]] <- c(0, pascal.list.s1[[i-1]]) + c(pascal.list.s1[[i-1]], 0)
  pascal.list.s1[[1]] <- NULL
  # For Simon's design, only need the final line:
  pascal.list.s1 <- pascal.list.s1[n1]

  # Curtail at n1 only:
  curtail.index <- c(1:(r1+1), (e1+2):(n1+1)) # We curtail at these indices -- Sm=[0, r1] and Sm=[e1+1, n1] (remembering row 1 is Sm=0)
  curtail.coefs.s1 <- pascal.list.s1[[1]][curtail.index] # from k=0 to k=r1

  # Use final column from S1:
  pascal.list.s2 <- pascal.list.s1
  pascal.list.s2[[1]][curtail.index] <- 0

  for (i in 2:(n2+1)) pascal.list.s2[[i]] <- c(0, pascal.list.s2[[i-1]]) + c(pascal.list.s2[[i-1]], 0)
  pascal.list.s2[[1]] <- NULL



  # Now obtain the rest of the probability -- the p^b * q^c :
  # S1
  q1 <- 1-p1
  coeffs <- p1^(0:n1)*q1^(n1:0)
  coeffs <- coeffs[curtail.index]

  q0 <- 1-p0
  coeffs.p0 <- p0^(0:n1)*q0^(n1:0)
  coeffs.p0 <- coeffs.p0[curtail.index]


  # Multiply the two vectors (A and p^b * q^c):
  prob.curt.s1 <- curtail.coefs.s1*coeffs

  # for finding type I error prob:
  prob.curt.s1.p0 <- curtail.coefs.s1*coeffs.p0

  # The (S1) curtailed paths:
  k.curt.s1 <- c(0:r1, (e1+1):n1)
  n.curt.s1 <- rep(n1, length(k.curt.s1))
  curtail.s1 <- cbind(k.curt.s1, n.curt.s1, prob.curt.s1, prob.curt.s1.p0)


  ############## S2 ###############

  # Pick out the coefficients for the S2 paths (A, say):
  s2.index <- (r1+2):(n+1)
  curtail.coefs.s2 <- pascal.list.s2[[n2]][s2.index]

  # Now obtain the rest of the probability -- the p^b * q^c :
  coeffs.s2 <- p1^(0:n)*q1^(n:0)
  coeffs.s2 <- coeffs.s2[s2.index]

  coeffs.s2.p0 <- p0^(0:n)*q0^(n:0)
  coeffs.s2.p0 <- coeffs.s2.p0[s2.index]

  # Multiply the two vectors (A and p^b * q^c):
  prob.go <- curtail.coefs.s2*coeffs.s2

  # for finding type I error prob:
  prob.go.p0 <- curtail.coefs.s2*coeffs.s2.p0


  # Paths that reach the end:
  k.go <- (r1+1):n
  n.go <- rep(n, length(k.go))

  go <- cbind(k.go, n.go, prob.go, prob.go.p0)

  final <- rbind(curtail.s1, go)


  ############## WRAPPING UP THE RESULTS ##############

  output <- data.frame(k=final[,1], n=final[,2], prob=final[,3], prob.p0=final[,4])

  output$success <- "Fail"
  output$success[output$k > r] <- "Success"
  output$success[output$n==n1 & output$k > e1] <- "Success"


  # Pr(early termination):
  #PET <- sum(output$prob[output$n < n])
  #PET.p0 <- sum(output$prob.p0[output$n < n])

  power <- sum(output$prob[output$success=="Success"])

  #output$obsd.p <- output$k/output$n

  #output$bias <- output$obsd.p - p

  #bias.mean <- wtd.mean(output$bias, weights=output$prob, normwt=TRUE)
  #bias.var <- wtd.var(output$bias, weights=output$prob, normwt=TRUE)

  #sample.size <- wtd.quantile(output$n, weights=output$prob, normwt=TRUE, probs=c(0.25, 0.5, 0.75))
  #sample.size.expd <- wtd.mean(output$n, weights=output$prob, normwt=TRUE)
  sample.size.expd <- sum(output$n*output$prob)

  #sample.size.p0 <- wtd.quantile(output$n, weights=output$prob.p0, normwt=TRUE, probs=c(0.25, 0.5, 0.75))
  #sample.size.expd.p0 <- wtd.mean(output$n, weights=output$prob.p0, normwt=TRUE)
  sample.size.expd.p0 <- sum(output$n*output$prob.p0)


  alpha <- sum(output$prob.p0[output$success=="Success"])

  #output <- list(output, mean.bias=bias.mean, var.bias=bias.var, sample.size=sample.size, expd.sample.size=sample.size.expd, PET=PET,
  #               sample.size.p0=sample.size.p0, expd.sample.size.p0=sample.size.expd.p0, PET.p0=PET.p0, alpha=alpha, power=power)
  to.return <- c(n1=n1, n2=n2, n=n, r1=r1, r=r, alpha=alpha, power=power, EssH0=sample.size.expd.p0, Ess=sample.size.expd, e1=e1)
  to.return
}

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
find2stageDesigns <- function(nmin, nmax, p0, p1, alpha, power, benefit=FALSE)
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
  simon.input <- data.frame(nmin=nmin, nmax=nmax, p0=p0, p1=p1, alpha=alpha, power=power)
  simon.output <- list(input=simon.input,
                       all.des=simon.df)
  class(simon.output) <- append(class(simon.output), "curtailment_simon")
  return(simon.output)
}


findCoeffs <- function(n, p0, p1){
  ##### Matrix of coefficients: Probability of a SINGLE path leading to each point:
  coeffs <- matrix(NA, ncol=n, nrow=n+1)
  coeffs.p0 <- coeffs
  q0 <- 1-p0
  q1 <- 1-p1
  for(i in 1:(n+1)){
    coeffs[i, ] <- p1^(i-1) * q1^((2-i):(2-i+n-1))
    coeffs.p0[i, ] <- p0^(i-1) * q0^((2-i):(2-i+n-1))
  }
  coeffs.list <- list(coeffs.p0=coeffs.p0,
                      coeffs=coeffs)
}







#' drawDiagram
#'
#' This function produces both a data frame and a diagram of stopping boundaries.
#' The function takes a single argument: the output from the function singlearmDesign.
#' If the supplied argument contains more than one admissible designs, the user is offered a choice of which design to use.
#' @param  findDesign.output Output from either the function singlearmDesign or find2stageDesigns
#' @param  print.row Choose a row number to directly obtain a plot and stopping boundaries for a particular design realisation. Default is NULL.
#' @param  xmax,ymax Choose values for the upper limits of the x- and y-axes respectively. Helpful for comparing two design realisations. Default is NULL.
#' @export
#' @author Martin Law, \email{martin.law@@mrc-bsu.cam.ac.uk}
#' @return The output is a list of two elements. The first, $diagram, is a ggplot2 object showing how the trial should proceed: when to to undertake an interim analysis, that is, when to check if a stopping boundary has been reached (solid colours) and what decision to make at each possible point (continue / go decision / no go decision). The second list element, $bounds.mat, is a data frame containing three columns: the number of participants at which to undertake an interim analysis (m), and the number of responses at which the trial should stop for a go decision (success) or a no go decision (fail).
#' @examples output <- singlearmDesign(nmin = 30,
#'  nmax = 30,
#'  C = 5,
#'  p0 = 0.1,
#'  p1 = 0.4,
#'  power = 0.8,
#'  alpha = 0.05)
#'  dig <- drawDiagram(output, print.row = 2)
drawDiagram <- function(findDesign.output, print.row=NULL, xmax=NULL, ymax=NULL){
  UseMethod("drawDiagram")
}

#' @export
drawDiagram.curtailment_single <- function(findDesign.output, print.row=NULL, xmax=NULL, ymax=NULL){
    des <- findDesign.output$all.des
  row.names(des) <- 1:nrow(des)
  if(!is.null(print.row)){
    des <- des[print.row, , drop=FALSE]
  }
  print(des)
  des.input <- findDesign.output$input
  if(nrow(des)>1){
    rownum <- 1
    while(is.numeric(rownum)){
      rownum <- readline("Input a row number to choose a design and see the trial design diagram. Press 'q' to quit: ")
      if(rownum=="q"){
        if(exists("plot.and.bounds")){
          return(plot.and.bounds)
        }else{
          print("No designs selected, nothing to return", q=F)
          return()
        }
      }else{
        rownum <- as.numeric(rownum)
        plot.and.bounds <- createPlotAndBounds(des=des, des.input=des.input, rownum=rownum, xmax=xmax, ymax=ymax)
      }
    } # end of while
  }else{
    print("Returning diagram and bounds for single design.", quote = F)
    plot.and.bounds <- createPlotAndBounds(des=des, des.input=des.input, rownum=1, xmax=xmax, ymax=ymax)
    return(plot.and.bounds)
  }
} # end of drawDiagram()

#' @export
drawDiagram.curtailment_simon <- function(findDesign.output, print.row=NULL, xmax=NULL, ymax=NULL){
  des <- findDesign.output$all.des
  row.names(des) <- 1:nrow(des)
  if(!is.null(print.row)){
    des <- des[print.row, , drop=FALSE]
  }
  print(des)
  des.input <- findDesign.output$input
  if(nrow(des)>1){
    rownum <- 1
    while(is.numeric(rownum)){
      rownum <- readline("Input a row number to choose a design and see the trial design diagram. Press 'q' to quit: ")
      if(rownum=="q"){
        if(exists("plot.and.bounds")){
          return(plot.and.bounds)
        }else{
          print("No designs selected, nothing to return", q=F)
          return()
        }
      }else{
        rownum <- as.numeric(rownum)
        plot.and.bounds <- createPlotAndBoundsSimon(des=des, des.input=des.input, rownum=rownum, xmax=xmax, ymax=ymax)
      }
    } # end of while
  }else{
    print("Returning diagram and bounds for single design.", quote = F)
    plot.and.bounds <- createPlotAndBoundsSimon(des=des, des.input=des.input, rownum=1, xmax=xmax, ymax=ymax)
    return(plot.and.bounds)
  }
}


#function generator
defunct <- function(msg = "This function is depreciated") function(...) return(stop(msg))
# @export
#' SCfn = defunct("SCfn changed name to findSCdesigns")
