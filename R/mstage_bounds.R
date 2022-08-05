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
  Sm <- 0:n # possible k~ values. # XXX not 0:(r+1) ?
  cp.subset1.Sm <- 0:r

  # Which columns? Trial success must still be possible, ie  n-m > r-Sm+1
  m <- 1:n # possible numbers of patients
  cp.subset1.m <- list()
  i <- 1
  for(k in cp.subset1.Sm){
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
  for(j in 1:length(cp.subset1.m)){
    current.Sm <- cp.subset1.Sm[j] # Sm of current row
    l <- seq(from=r-current.Sm, length=length(cp.subset1.m[[j]]))
    summation[[j]] <- choose(l, r-current.Sm) * p^(r-current.Sm+1) * q^(l-(r-current.Sm))
    }
  cp1 <- lapply(summation, cumsum)
  cp1 <- lapply(cp1, rev)
  # Now, insert CPs into matrix of CPs:
  for(i in 1:length(cp.subset1.Sm)){
    mat[cp.subset1.Sm[i]+1, cp.subset1.m[[i]]] <- cp1[[i]]
    }
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
  theta.test <- function(SM, M){
    sum(mat[SM+1, M] < thetaF-1e-6 | mat[SM+1, M] > thetaE + 1e-6) # The 1e-6 terms are added to deal with floating point errors.
    }

  low.cp.rows <- mapply(theta.test, SM=cp.sm, M=cp.m)
  # Because of how CP changes, propagating right to left and upwards, we only
  # need to identify the "lowest" row with CP < thetaF or CP > thetaE, ie the row with greatest Sm
  # that contains a CP < thetaF or CP > thetaE. This is the row where we begin. Note: This may be all rows, or even none!
  # Recall that row number = Sm-1

  # We only need to act if there are CPs < thetaF or > thetaE -- otherwise, the matrix of CPs does not change.
  if(sum(low.cp.rows)>0){    # There may be some cases where there are no CPs to be changed; i.e. no stopping due to stochastic curtailment
    begin.at.row <- max(which(low.cp.rows>0))
    cp.changing.sm <- cp.sm[1:begin.at.row]
    cp.changing.m <- cp.m[1:begin.at.row]
    # We calculate the CP, truncate to zero if CP<thetaF (or to 1 if CP > thetaE), then move on:
    for(rown in rev(cp.changing.sm+1)){
      for(coln in rev(cp.changing.m[[rown]])){
        currentCP <- q*mat[rown, coln+1] + p*mat[rown+1, coln+1]
        if(currentCP > thetaE + 1e-6)   mat[rown, coln] <- 1   # If CP > thetaE, amend to equal 1. The term 1e-6 is added to deal with floating point errors.
        else  mat[rown, coln] <- ifelse(test = currentCP < thetaF - 1e-6, yes=0, no=currentCP)  # Otherwise, test if CP < thetaF. If so, amend to 0, otherwise calculate CP as normal. Again, the term 1e-6 is added to deal with floating point errors.
      }
    } # Again, this *must* be done one entry at a time -- cannot vectorise.
  } # End of IF statement

  pascal.list <- list(1, c(1,1))
  for(i in 3:(n+2)){
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
  for(i in rows.with.cp1){
    columns.of.rows.w.cp1[[j]] <- which(mat[i, ]==1 & !is.na(mat[i, ]))
    j <- j+1
    }

    index1 <- 1
  success <- NULL

  for(i in rows.with.cp1){
    for(j in columns.of.rows.w.cp1[[index1]]){
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

  for(i in 1:final.zero.in.col){
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
