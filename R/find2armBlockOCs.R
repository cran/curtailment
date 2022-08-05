find2armBlockOCs <- function(n,r, Bsize, mat, thetaF, thetaE, power, alpha, pat.cols, prob.vec, prob.vec.p0, prob.vec.minstop, prob.vec.p0.minstop, blank.mat, zero.mat, minstop){
######################## UPDATE CP MATRIX USING thetaF/1 VALUES:
for(i in (n+r+1):1){
  for(j in pat.cols){  # Only look every Bsize patients (no need to look at final col)
    if(i-1<=j){ # Condition: Sm<=m
      newcp <- sum(prob.vec*mat[i:(i+Bsize), j+Bsize])
      if(newcp > thetaE) mat[i,j] <- 1
      if(newcp < thetaF) mat[i,j] <- 0
      if(newcp <= thetaE & newcp >= thetaF) mat[i,j] <- newcp
    }
  }
}
###### STOP if design is pointless, i.e either failure or success is not possible:
# IF DESIGN GUARANTEES FAILURE (==0) or SUCCESS (==2) at n=C:
first.cohort <- sum(mat[, minstop], na.rm = T)
if(first.cohort==minstop+1){
  return(c(n, r, Bsize, 1, 1,  NA, NA, thetaF, thetaE, NA))
}
if(first.cohort==0){
  return(c(n, r, Bsize, 0, 0,  NA, NA, thetaF, thetaE, NA))
}

########################### FIND PROB. OF REACHING EACH POINT:
################# START WITH AN INDICATOR MATRIX OF NON-TERMINAL POINTS:
tp.mat <- blank.mat
tp.mat[which(mat==0 | mat==1)] <- 0
tp.mat[which(mat>0 & mat<1)] <- 1

############ CREATE MATRIX OF "POSSIBLE POINTS" -- IE, LIKE MATRIX OF NON-TERMINAL POINTS ***PLUS*** THE TERMINAL POINTS:
##### THIS WILL BE USED AS AN INDICATOR MATRIX FOR WHICH POINTS TO CALCULATE THE PROB'Y OF REACHING.
# Start with non-terminal points and add the terminal points
poss.mat <- tp.mat
poss.mat[1:(minstop+1), minstop] <- 1 # It is of course possible to reach 0,1,..., minstop responses after minstop patients. This line is included in case CP=0 or CP=1 at first check (i.e. m=minstop)
# Failures first:
fail.mat <- zero.mat
rows.with.cp0 <- which(apply(mat, 1, function(x) {any(x==0, na.rm = T)}))
fail.n <- apply(mat[rows.with.cp0,], 1, which.min)
for(i in 1:length(fail.n)){
  poss.mat[names(fail.n)[i],fail.n[i]] <- 1
  fail.mat[names(fail.n)[i],fail.n[i]] <- 1
}
# Now successes: what are the successful terminal points?
# Points with mat[i,j]==1 AND (0<mat[i,j-2]<1 OR 0<mat[i-1,j-2]<1 OR 0<mat[i-2,j-2]<1)
success.mat <- zero.mat
rows.with.cp1 <- which(apply(mat, 1, function(x) {any(x==1, na.rm = T)}))
for(i in rows.with.cp1){
  for(j in seq(minstop, 2*n, by=Bsize)){
    if(i-1<=j & mat[i,j]==1 & (j==minstop | any(tp.mat[i:max(1,(i-Bsize)), j-Bsize]==1, na.rm=TRUE))){ # max() condition to take care of cases where
        # Conditions ensure CP=1 and that it is possible to actually reach the point (tp.mat==1 indicates non terminal point)
        poss.mat[i, j] <- 1
        success.mat[i,j] <- 1
    }
  }
}
######################## PROBABILITY OF REACHING EACH POINT
############################## FIRSTLY, UNDER PT=PT
final.probs.mat <- poss.mat
# First n=Bsize rows (0,1,... Bsize-1) are special cases:
# fill in first column:
#final.probs.mat[1:(Bsize+1), Bsize] <- prob.vec
final.probs.mat[1:(minstop+1), minstop] <- prob.vec.minstop
for(i in 1:Bsize){
  poss.cols.index <- which(poss.mat[i,]==1)[-1] # Columns in row i that are possible to reach.
  # Note: first entry (ie first column) has been inputted already, directly above, hence [-1]
  for(j in poss.cols.index){
    final.probs.mat[i, j] <- sum(prob.vec[1:i]*final.probs.mat[i:1, j-Bsize]*tp.mat[i:1, j-Bsize])
    # j-Bsize : column of previous interim
  }
}
# For the remaining rows:
for(i in (Bsize+1):nrow(final.probs.mat)){ # Skip first Bsize rows; they have been taken care of above.
  for(j in seq(minstop+Bsize, 2*n, by=Bsize)){ # skipping first interim (at minstop), aka first column of patients -- again, they have been taken care of above.
    if(i-1<=j & poss.mat[i,j]==1){
      final.probs.mat[i,j] <- sum(prob.vec*final.probs.mat[i:(i-Bsize), j-Bsize]*tp.mat[i:(i-Bsize), j-Bsize], na.rm = TRUE)
    }
  }
}
# IMPORTANT: end early if pwr < power.
# Note 2: Could comment this out to ensure that the final test in undertaken for all runs, not just feasible ones.
prob.success <- final.probs.mat[success.mat==1]
pwr <- sum(prob.success)
if(is.na(pwr)) browser()
if(pwr < power) # | pwr < power+tol )
  {
  return(c(n, r, Bsize, NA, pwr,  NA, NA, thetaF, thetaE, NA))
}

############################## SECONDLY, UNDER PT=PC
final.probs.mat.p0 <- poss.mat
# First n=Bsize rows (0,1,... Bsize-1) are special cases:
# fill in first column:
#final.probs.mat.p0[1:(Bsize+1), Bsize] <- prob.vec.p0
final.probs.mat.p0[1:(minstop+1), minstop] <- prob.vec.p0.minstop
for(i in 1:Bsize){
  poss.cols.index <- which(poss.mat[i,]==1)[-1] # Columns in row i that are possible to reach.
  # Note: first entry (ie first column) has been inputted already, directly above, hence [-1]
  for(j in poss.cols.index){
    final.probs.mat.p0[i, j] <- sum(prob.vec.p0[1:i]*final.probs.mat.p0[i:1, j-Bsize]*tp.mat[i:1, j-Bsize])
  }
}

# For the remaining rows:
for(i in (Bsize+1):nrow(final.probs.mat.p0)){ # Skip first Bsize rows; they have been taken care of above.
  for(j in seq(minstop+Bsize, 2*n, by=Bsize)){ # skipping first interim (at minstop), aka first column of patients -- again, they have been taken care of above.
    if(i-1<=j & poss.mat[i,j]==1){
      final.probs.mat.p0[i,j] <- sum(prob.vec.p0*final.probs.mat.p0[i:(i-Bsize), j-Bsize]*tp.mat[i:(i-Bsize), j-Bsize], na.rm = TRUE)
    }
  }
}
prob.success.p0 <- final.probs.mat.p0[success.mat==1]
typeIerr <- sum(prob.success.p0)
# IMPORTANT: end early if type I error > alpha
# Note 2: Could comment this out to ensure that the final test in undertaken for all runs, not just feasible ones.
if( typeIerr > alpha) # | alpha < alpha-tol)
{
  return(c(n, r, Bsize, typeIerr, pwr,  NA, NA, thetaF, thetaE, NA))
}
########################## ESS FOR SUCCESS POINTS
success.n <- which(success.mat==1, arr.ind = T)[,"col"] # Note: this is potentially dangerous if this and final.probs.mat[success.mat==1] are found in a different order, though this shouldn't be the case.
success.df <- data.frame(prob=prob.success, prob.p0=prob.success.p0, ess=prob.success*success.n, essH0=prob.success.p0*success.n)

################## ESS FOR FAILURE POINTS
fail.n <- which(fail.mat==1, arr.ind = T)[,"col"] # Note: this is potentially dangerous if this and final.probs.mat[success.mat==1] are found in a different order, though this shouldn't be the case.
prob.fail <- final.probs.mat[fail.mat==1]
prob.fail.p0 <- final.probs.mat.p0[fail.mat==1]
fail.df <- data.frame(prob=prob.fail, prob.p0=prob.fail.p0, ess=prob.fail*fail.n, essH0=prob.fail.p0*fail.n)
all.df <- rbind(fail.df, success.df)
###### CHECK PROBS ALL SUM TO 1
if(sum(all.df$prob)+sum(all.df$prob.p0)-2 > 1e-8)  stop("Total probability of failure + success =/= 1. Something has gone wrong." , call. = FALSE)
ess <- sum(all.df$ess)
essH0 <- sum(all.df$essH0)

############ EFFECTIVE N. THE "EFFECTIVE N" OF A STUDY IS THE "REAL" MAXIMUM SAMPLE SIZE
######## The point is where every Sm for a given m equals zero or one is necessarily where a trial stops
cp.colsums <- apply(mat, 2, function(x) { sum(x==0, na.rm=TRUE)+sum(x==1, na.rm=TRUE)} ) # Sum the CP values that equal zero or one in each column
possible.cps <- apply(mat, 2, function(x) {sum(!is.na(x))})
effective.n <- min(which(cp.colsums==possible.cps))
return(data.frame(n=n, r=r, Bsize=Bsize, typeIerr=typeIerr, pwr=pwr, EssH0=essH0, Ess=ess, thetaF=thetaF, thetaE=thetaE, eff.n=effective.n))
}
