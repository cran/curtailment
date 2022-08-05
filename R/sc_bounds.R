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
  ### ^^^ check this inequality

  # Which columns? Trial success must still be possible, ie  n-n~ > r-Sm+1
  m <- 1:n # possible numbers of patients
  cp.subset1.m <- list()
  i <- 1

  for(k in cp.subset1.Sm){
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
  for(j in 1:length(cp.subset1.m)){
    current.Sm <- cp.subset1.Sm[j]
    l <- seq(from=r2-current.Sm, length=length(cp.subset1.m[[j]]))
    summation[[j]] <- choose(l, r2-current.Sm) * p^(r2-current.Sm+1) * q^(l-(r2-current.Sm))
    }
  cp1 <- lapply(summation, cumsum)
  cp1 <- lapply(cp1, rev)

  # Now, insert CPs into matrix of CPs:
  for(i in 1:length(cp.subset1.Sm)){
    mat[cp.subset1.Sm[i]+1, cp.subset1.m[[i]]] <- cp1[[i]]
    }

  # Final region to be addressed: Points where k <= r1 but progression is still possible (S1 only):
  cp.subset2.Sm <- 0:r1  # Values of Sm <= r1
  # Which columns? Progression to S2 must still be possible, ie  n1-m > r1-Sm+1
  n.to.n1 <- 1:n1 # possible S1 values
  cp.subset2.m <- list()
  i <- 1
  for(k in cp.subset2.Sm){
    cp.subset2.m[[i]] <- which(n1-n.to.n1 >= r1-k+1 & n.to.n1 >= k)
    i <- i+1
    }

  # Now we have the rows and columns of this region:
  # cp.subset2.Sm # Values of Sm. Rows are this plus 1.
  # cp.subset2.m # Values of m (columns)
  # Have to propagate "backwards", ending at Sm=0, n=1.
  # As with previous region, proceed row-wise -- start with greatest Sm.
  for(k in length(cp.subset2.Sm):1){ # Note that we END at the earliest row.
    for(j in rev(cp.subset2.m[[k]])){
      mat[k, j] <- q*mat[k, j+1] + p*mat[k+1, j+1]
      }
    }
  # This *must* be done in this way -- cannot be vectorised.

  # Create the "blanks":
  for(i in 1:(ncol(mat)-1)){
    mat[(i+2):nrow(mat), i] <- NA
    }
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

  theta.test <- function(SM, M){
    sum(mat[SM+1, M] < thetaF | mat[SM+1, M] > thetaE)
    }

  low.cp.rows <- mapply(theta.test, SM=cp.sm, M=cp.m)
  # Because of how CP changes, propagating right to left and upwards, we only
  # need to identify the "lowest" row with CP < theta or CP > thetaE, ie the row with greatest Sm
  # that contains a CP < theta or CP > thetaE. This is the row where we begin. Note: This may be all rows, or even none!
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
  # Might be best to sort into Failure (Sm <= r) and Success (Sm = r+1):
  final.probs.mat <- matrix(unlist(lapply(final.probs, '[', 1:max(sapply(final.probs, length)))), ncol = n, byrow = F)
  final.probs.mat.p0 <- matrix(unlist(lapply(final.probs.p0, '[', 1:max(sapply(final.probs.p0, length)))), ncol = n, byrow = F)

  # FIRST: Search for terminal points of success. These can only exist in rows where (updated) CP=1, and where Sm<=r+1:
  potential.success.rows <- rowSums(mat[1:(r2+2), ]==1, na.rm = TRUE)
  rows.with.cp1 <- which(as.numeric(potential.success.rows)>0)
  # ^ These are the rows containing possible terminal points of success. They must have CP=1:

  columns.of.rows.w.cp1 <- list()
  j <- 1
  for(i in rows.with.cp1){
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

  eff.bound.Sm <- rep(Inf, n) # Inf instead of n+1 -- either way, no stopping
  eff.bound.Sm[success[,"m"]] <- success[,"Sm"]

  fut.bound.Sm <- rep(-Inf, n) # -Inf instead of 0 -- either way, no stopping
  fut.bound.Sm[m.fail] <- Sm.fail

  to.return <- list(J=2, n=c(n1, n2), N=n, a=c(r1, r2), r=c(Inf, r2+1), a_curt=fut.bound.Sm, r_curt=eff.bound.Sm, alpha=alpha, beta=beta)
  to.return
}
