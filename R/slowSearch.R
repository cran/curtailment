slowSearch <- function(thetas,
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
                       minstop
                       ){
  k <- 1
  h.results <- vector("list", ceiling(0.5*length(thetas)^2))
  all.thetas <- rev(thetas)[-length(thetas)] # All thetaE values, decreasing, not including the final value, thetaE=0.
  all.thetas <- all.thetas[all.thetas>=minthetaE]
  for(i in 1:length(all.thetas)){ # For each thetaE,
    thetaFs.current <- thetas[thetas<all.thetas[i] & thetas<=maxthetaF]
    for(m in 1:length(thetaFs.current)){
      h.results[[k]] <- find2armBlockOCs(n=sc.h$n, r=sc.h$r, Bsize=Bsize, thetaF=as.numeric(thetaFs.current[m]), thetaE=all.thetas[i], mat=mat.h,
                                         power=power, alpha=alpha, pat.cols=pat.cols.single,
                                         prob.vec=prob.vec, prob.vec.p0=prob.vec.p0, prob.vec.minstop=prob.vec.minstop, prob.vec.p0.minstop=prob.vec.p0.minstop,
                                         blank.mat=blank.mat, zero.mat=zero.mat, minstop=minstop)
      k <- k+1
      # Add lines here: if power decreases below desired value, break:
      if(h.results[[k-1]][5] < power){
        break
      }
    }
  } # end of "i" loop
  return(h.results)
}
