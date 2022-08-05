
#### Part 4: from "inside_theta_oct2020.R" ####
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
