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

findDesignOCs <- function(n, r, C, thetaF, thetaE, mat, power, alpha, coeffs, coeffs.p0, p0, p1, return.tps=FALSE, minstop){

  q0 <- 1-p0
  q1 <- 1-p1

  interims <- seq(minstop, n, by=C)
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

  ###### STOP if design is pointless, i.e either failure or success is not possible:
  # IF DESIGN GUARANTEES FAILURE (==0) or SUCCESS (==2) at n=minstop:
  #first.cohort <- sum(mat[,C], na.rm = T)
  first.cohort <- sum(mat[,minstop], na.rm = T)
  if(first.cohort==minstop+1){
    return(c(n, r, C, 1, 1,  NA, NA, thetaF, thetaE, NA, minstop))
  }
  if(first.cohort==0){
    return(c(n, r, C, 0, 0,  NA, NA, thetaF, thetaE, NA, minstop))
  }

  ############# Number of paths to each point:
  pascal.list <- list(c(1,1))
  if(C>1 | minstop>1){
    for(i in 2:minstop){
      pascal.list[[i]] <- c(0, pascal.list[[i-1]]) + c(pascal.list[[i-1]], 0)
    }
  }

  for(i in min(minstop+1,n):n){ # Adding min() in case of stages=1, i.e. cohort size C=n.
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
  pascal.t <- pascal.t[1:max(min(r+C+1, n+1), minstop+1), ] # Adding min() in case of stages=1, i.e. cohort size C=n.

  # Multiply the two matrices (A and p^b * q^c). This gives the probability of reaching each point:
  # Note: Only need the first r+C+1 (or first minstop+1) rows -- no other rows can be reached.
  coeffs <- coeffs[1:max(min(r+C+1, n+1), minstop+1), 1:n]
  coeffs.p0 <- coeffs.p0[1:max(min(r+C+1, n+1), minstop+1), 1:n]
  pascal.t <- pascal.t[1:max(min(r+C+1, n+1), minstop+1), ]
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

  success <- data.frame(Sm=Sm.success, m=m.success, prob=prob.success, success=rep("Success", length(m.success)))
  # CAN FIND POWER AND TYPE I ERROR AT THIS POINT.
  # IF TOO LOW OR HIGH, STOP AND MOVE ON:
  pwr <- sum(prob.success)
  typeI <- sum(prob.success.p0)
  if(pwr < power | typeI > alpha)
  {
    return(c(n, r, C, typeI, pwr,  NA, NA, thetaF, thetaE, NA, minstop))
  }


  ##### FAILURE:
  submat <- mat[1:(r+1), ]
  first.zero.in.row <- apply(submat, 1, function(x) match(0, x))
  m.fail <- as.numeric(first.zero.in.row)
  Sm.fail <- as.numeric(names(first.zero.in.row))

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

  sample.size.expd.fail <- sum(m.fail*fail.deets.prob)
  sample.size.expd.p0.fail <- sum(m.fail*fail.deets.prob.p0)

  sample.size.expd.success <- sum(m.success*prob.success)
  sample.size.expd.p0.success <- sum(m.success*prob.success.p0)

  sample.size.expd <- sample.size.expd.fail + sample.size.expd.success
  sample.size.expd.p0 <- sample.size.expd.p0.fail + sample.size.expd.p0.success


  ############ EFFECTIVE N. THE "EFFECTIVE N" OF A STUDY IS THE "REAL" MAXIMUM SAMPLE SIZE
  ######## The point is where every Sm for a given m equals zero or one is necessarily where a trial stops
  cp.colsums <- colSums(mat==0, na.rm = T) + colSums(mat==1, na.rm = T)
  possible.cps <- colSums(!is.na(mat))
  effective.n <- min(which(cp.colsums==possible.cps))
  return.vec <- c(n, r, C, typeI, pwr, sample.size.expd.p0, sample.size.expd, thetaF, thetaE, effective.n, minstop)
  if(return.tps==TRUE){
    return.list <- list(ocs=return.vec, tp=tp)
    return(return.list)
    }
  return.vec
}

