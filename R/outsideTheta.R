#### Part 2: from "outside_theta.R" ####
# Note: Can run much of the code outside the loops for theta -- independent of theta:
outsideTheta <- function(n1=4, n2=4, r1=1, r2=4, p0, p){
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
  ### ^^^ check this inequality!
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
