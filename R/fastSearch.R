fastSearch <- function(thetas,
                       maxthetaF,
                       minthetaE,
                       sc.h,
                       Bsize,
                       mat.h,
                       blank.mat,
                       zero.mat,
                       power,
                       alpha,
                       pat.cols.single,
                       prob.vec,
                       prob.vec.p0,
                       prob.vec.minstop,
                       prob.vec.p0.minstop,
                       minstop){
  k <- 1
  ########### START 2D BISECTION
  thetaF.vec <- thetas[thetas<=maxthetaF]
  thetaE.vec <- thetas[thetas>=minthetaE]
  h.results <- vector("list", length(thetaF.vec)*length(thetaE.vec))
  # Bounds for thetaF:
  a0 <- 1
  b0 <- length(thetaF.vec)
  d0 <- ceiling((b0-a0)/2)
  # Bounds for thetaE:
  a1 <- 1
  b1 <- length(thetaE.vec)
  d1 <- ceiling((b1-a1)/2)
  minthetaF <- NA
  maxthetaE <- NA
  while(min((b0-a0),(b1-a1))>1 & is.na(minthetaF)){ # Break/move on when bisection method fails to find anything OR when final feasible design is found.
    output <- find2armBlockOCs(n=sc.h$n, r=sc.h$r, Bsize=Bsize, thetaF=thetaF.vec[d0], thetaE=thetaE.vec[d1], mat=mat.h,  power=power, alpha=alpha,
                               pat.cols=pat.cols.single, prob.vec=prob.vec, prob.vec.p0=prob.vec.p0, prob.vec.minstop=prob.vec.minstop, prob.vec.p0.minstop=prob.vec.p0.minstop, blank.mat=blank.mat, zero.mat=zero.mat, minstop=minstop)
    if(!is.na(output[6])){ # If ESS is not NA, then design IS feasible, and do:
      feasible <- TRUE
      maxthetaE <- thetaE.vec[d1]
      while((feasible==TRUE) & d0<length(thetaF.vec)){
        d0 <- d0+1
        h.results[[k]] <- find2armBlockOCs(n=sc.h$n, r=sc.h$r, Bsize=Bsize, thetaF=thetaF.vec[d0], thetaE=maxthetaE, mat=mat.h,  power=power, alpha=alpha,
                                           pat.cols=pat.cols.single, prob.vec=prob.vec, prob.vec.p0=prob.vec.p0, prob.vec.minstop=prob.vec.minstop, prob.vec.p0.minstop=prob.vec.p0.minstop,
                                           blank.mat=blank.mat, zero.mat=zero.mat, minstop=minstop)
        feasible <- !is.na(h.results[[k]][6])
        k <- k+1
      } # Once the final feasible design for the given thetaF/1 is found (or we reach the largest thetaF), record thetaF and make it a limit:
      minthetaF <- thetaF.vec[d0-1]
    } else { # If design isn't feasible, decrease thetaF, increase thetaE and test again:
      b0 <- d0
      a1 <- d1
      d0 <- a0 + floor((b0-a0)/2)
      d1 <- a1 + floor((b1-a1)/2)
    }
  }
  if(!is.na(minthetaF)){ # If at least one feasible design was found, then minthetaF exists, and we search over all thetaF/1 combinations subject to the new limits we have just created:
    thetaF.vec <- thetaF.vec[thetaF.vec>=minthetaF]
    thetaE.vec <- thetaE.vec[thetaE.vec<=maxthetaE]
    for(i in 1:length(thetaE.vec)){
      for(j in 1:length(thetaF.vec)){
        #  print(paste(thetaF.vec[i], thetaE.vec[j]))
        h.results[[k]] <- find2armBlockOCs(n=sc.h$n, r=sc.h$r, Bsize=Bsize, thetaF=thetaF.vec[j], thetaE=thetaE.vec[i], mat=mat.h,
                                           power=power, alpha=alpha, pat.cols=pat.cols.single, prob.vec=prob.vec, prob.vec.p0=prob.vec.p0, prob.vec.minstop=prob.vec.minstop, prob.vec.p0.minstop=prob.vec.p0.minstop,
                                           blank.mat=blank.mat, zero.mat=zero.mat, minstop=minstop)
        k <- k+1
      }
    }
  } # if no feasible designs found, do nothing and let loop end.
  return(h.results)
}
