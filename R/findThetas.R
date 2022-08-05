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

