
#' findNSCdesigns
#'
#' This function finds admissible design realisations for single-arm binary outcome trials, using non-stochastic curtailment.
#' The output is a data frame of admissible design realisations.
#' @param nmin Minimum permitted sample size.
#' @param nmax Maximum permitted sample size.
#' @param p0 Probability for which to control the type-I error-rate
#' @param p1 Probability for which to control the power
#' @param alpha Significance level
#' @param power Required power (1-beta).
#' @param progressBar Logical. If TRUE, shows progress bar. Defaults to FALSE.
#' @return Output is a list of two dataframes. The first, $input, is a one-row
#' data frame that contains important arguments used in the call. The second,
#' $all.des,contains the operating characteristics of all admissible designs found.
#' @export
#' @examples findNSCdesigns(nmin=20, nmax=21, p0=0.1, p1=0.4, alpha=0.1, power=0.8)
findNSCdesigns <- function(nmin, nmax, p0, p1, alpha, power, progressBar=FALSE)
{
  #system.time({
  nr.lists <- findN1N2R1R2NSC(nmin, nmax)

  n1.vec <- nr.lists$n1
  n2.vec <- nr.lists$n2
  n.vec <- nr.lists$n

  r1 <- nr.lists$r1
  r <-  nr.lists$r

  alpha.power.nsc <-  vector("list", length(n.vec))
  l <- 1

  ns <- nr.lists$ns

  if(progressBar==TRUE) pb <- txtProgressBar(min = 0, max = nrow(ns), style = 3)

  for(i in 1:nrow(ns)) # For every combination of n1/n2/n,
  {
    # print(i)
    n1 <- n1.vec[i]
    n2 <- n2.vec[i]
    n <- n.vec[i]

    n.to.n1 <- 1:n1 # possible S1 values
    Sm <- 0:n
    m <- 1:n

    for(j in 1:length(r1[[i]])) # For every possible r1,
    {
      r1.j <- r1[[i]][[j]]
      r2.j <- r[[i]][[j]]

      cp.subset1.Sm.list <- lapply(r2.j, function(x) {which(r1.j < Sm & Sm < (x+1)) - 1})

      cp.subset2.Sm <- 0:r1.j
      cp.subset2.m <- lapply(cp.subset2.Sm, function(x) {which(n1-n.to.n1 >= r1.j-x+1 & n.to.n1 >= x)})


      for(k in 1:length(r[[i]][[j]])) # For every possible r,
      {
        alpha.power.nsc[[l]] <- findNSCerrorRates(n1=n1, n2=n2, r1=r1.j, r2=r[[i]][[j]][k], p0=p0, p1=p1, Sm=Sm, m=m, n.to.n1=n.to.n1,
                                                  cp.subset2.Sm=cp.subset2.Sm, cp.subset2.m=cp.subset2.m, cp.subset1.Sm=cp.subset1.Sm.list[[k]])
        l <- l+1
      }
    }
     if(progressBar==TRUE) setTxtProgressBar(pb, i)
  }
  #})

  alpha.power.nsc <- do.call(rbind.data.frame, alpha.power.nsc)
  names(alpha.power.nsc) <- c("n1", "n2", "n", "r1", "r", "alpha", "power")
  #nrow(alpha.power.nsc)
  # nmax = 80 : 1.7 million rows, takes 20 mins
  # nmax = 60 : 600,000 rows, takes 4 mins
  # nmax = 50 : 300,000 rows, takes 2 mins
  # nmax = 40 : 100,000 rows, takes 1 min
  # n=[61,65] :
  #save.image("march09.RData")


  # NOW THAT ALPHA AND POWER IS KNOWN FOR ALL COMBNS OF n1/n2/n/r1/r (given p/p0),
  # FIND THE OPTIMAL DESIGN (UNDER NSC)

  nsc.search <- alpha.power.nsc$power > power & alpha.power.nsc$alpha < alpha
  results.nsc.search <- alpha.power.nsc[nsc.search, ]
  if(sum(nsc.search)>0)##
  {##
    ess.nsc.search <- apply(results.nsc.search, 1, function(x) {findNSCdesignOCs(n1=x[1], n2=x[2], r1=x[4], r2=x[5], p0=p0, p1=p1)})
    # <2mins to do 60k designs; 1 second to do 700 designs

    ess.n.search <- as.data.frame(t(ess.nsc.search))

    # For NSC df above, remove dominated and duplicated designs:
    discard <- rep(NA, nrow(ess.n.search))
    for(i in 1:nrow(ess.n.search))
    {
      discard[i] <- sum(ess.n.search$EssH0[i] > ess.n.search$EssH0 & ess.n.search$Ess[i] > ess.n.search$Ess & ess.n.search$n[i] >= ess.n.search$n)
    }
    subset.nsc <- ess.n.search[discard==0,]

    # Remove duplicates:
    duplicates <- duplicated(subset.nsc[, c("n", "Ess", "EssH0")])
    subset.nsc <- subset.nsc[!duplicates,]

    rm(results.nsc.search)

  } else {subset.nsc <- rep(NA, 9)}
  names(subset.nsc) <- c("n1", "n2", "n", "r1", "r", "alpha", "power", "EssH0", "Ess")
  nsc.input <- data.frame(nmin=nmin, nmax=nmax, p0=p0, p1=p1, alpha=alpha, power=power)
  nsc.output <- list(input=nsc.input,
                     all.des=subset.nsc)
  return(nsc.output)

}


findNSCerrorRates <- function(n1, n2, r1, r2, p1, p0, theta0=0, theta1=1,
                              cp.subset2.Sm=cp.subset2.Sm, cp.subset2.m=cp.subset2.m, Sm=Sm, m=m, n.to.n1=n.to.n1, cp.subset1.Sm=cp.subset1.Sm)
  #nsc.power.alpha.only2 <- function(n1=4, n2=4, r1=1, r2=4, p=0.25, p0=0.1, theta0=0, theta1=1)
{

  # if(theta0 >= theta1) stop("theta1 must be greater than theta0")
  # if(p0 >= p) stop("p1 must be greater than p")
  # if(r1 >= r2) stop("r2 must be greater than r1")
  #


  n <- n1+n2
  q1 <- 1-p1
  q0 <- 1-p0

  # Start with all zeros (why not):
  #mat <- matrix(0, nrow = n+1, ncol = n)
  #rownames(mat) <- 0:n
  # Add the 1's for CP=1:
  #mat[(r2+2):(n+1),] <- 1  # r+2 === Sm=r+1, n+1 === Sm=n

  # Quicker than the above:
  mat <- matrix(c(rep(0, n*(r2+1)), rep(1, n*((n+1)-(r2+1)))), nrow = n+1, byrow=T)


  # Now input CP for points where r1 < k < r+1 (ie progression to Stage 2 is guaranteed):
  # Sm <- 0:n # possible k~ values.
  # cp.subset1.Sm <- which(r1 < Sm & Sm < (r2+1)) - 1  # Values for Sm where r1 < Sm < r2+1, i.e. progression to S2 certain, but overall trial success not guaranteed.

  # Which columns? Trial success must still be possible, ie  n-m > r-Sm+1
  # m <- 1:n # possible numbers of patients

  cp.subset1.m <- list()

  i <- 1

  for(k in cp.subset1.Sm)
  {
    cp.subset1.m[[i]] <- which(n-m >= r2-k+1 & m>=k)
    i <- i+1
  }
  # These are the m's for which trial success is still possible when Sm is cp.subset1.Sm
  # Note: This code seems faster than using "cp.subset1.m <- lapply(cp.subset1.Sm, function(x) {which(n-m >= r2-x+1 & m>=x)})"


  ##### Now calculate CP for all these points: #####

  # For each Sm, begin with the earliest point, ie lowest m, and work row-wise to the latest point, ie highest m:
  # Remember that row number is Sm+1.
  # How many rows to do? length(cp.subset1.m):

  # summation <- list()
  #
  # for(j in 1:length(cp.subset1.m))
  # {
  #   current.Sm <- cp.subset1.Sm[j]
  #
  #   l <- seq(from=r2-current.Sm, length=length(cp.subset1.m[[j]]))
  #
  #   summation[[j]] <- choose(l, r2-current.Sm) * p^(r2-current.Sm+1) * q^(l-(r2-current.Sm))
  # }
  #
  # cp1 <- lapply(summation, cumsum)
  # cp1 <- lapply(cp1, rev)
  #
  #
  # # Now, insert CPs into matrix of CPs:
  #
  # for(i in 1:length(cp.subset1.Sm))  {  mat[cp.subset1.Sm[i]+1, cp.subset1.m[[i]]] <- cp1[[i]] }
  #
  #
  # # Final region to be addressed: Points where k <= r1 but progression is still possible (S1 only):
  # cp.subset2.Sm <- 0:r1  # Values of Sm <= r1
  #
  # # Which columns? Progression to S2 must still be possible, ie  n1-m > r1-Sm+1
  # n.to.n1 <- 1:n1 # possible S1 values
  #
  # cp.subset2.m <- list()
  # i <- 1
  # for(k in cp.subset2.Sm)
  # {
  #   cp.subset2.m[[i]] <- which(n1-n.to.n1 >= r1-k+1 & n.to.n1 >= k)
  #   i <- i+1
  # }

  # Combine all the non-zero, non-one CP points:
  cp.subset21.Sm <- c(cp.subset2.Sm, cp.subset1.Sm)
  cp.subset21.m <- c(cp.subset2.m, cp.subset1.m)

  for(i in 1:length(cp.subset21.Sm))
  {
    mat[cp.subset21.Sm[i]+1, cp.subset21.m[[i]]] <- 0.5
  }

  # Now we have the rows and columns of this region:

  # cp.subset2.Sm # Values of Sm. Rows are this plus 1.
  # cp.subset2.m # Values of m (columns)

  # Have to propagate "backwards", ending at Sm=0, n=1.
  # As with previous region, proceed row-wise -- start with greatest Sm.

  # for(k in length(cp.subset2.Sm):1) # Note that we END at the earliest row.
  # {
  #   for(j in rev(cp.subset2.m[[k]]))
  #     {
  #     mat[k, j] <- q*mat[k, j+1] + p*mat[k+1, j+1]
  #     }
  # }
  # This *must* be done in this way -- cannot be vectorised.

  # Create the "blanks":
  # for(i in 1:(ncol(mat)-1))
  # {
  #   mat[(i+2):nrow(mat), i] <- NA
  # }
  NAs <- rbind(FALSE, lower.tri(mat)[-nrow(mat),])
  mat[NAs] <- NA

  # all.thetas <- unique(c(mat))
  # subset.thetas <- all.thetas[all.thetas < 0.2]

  # We now have a matrix m containing the CPs in (more or less) the upper triangle.


  ######################
  ######IMPORTANT ######
  ######################

  # If a point in a path has CP < theta, we stop the trial due to stochastic curtailment.
  # If an earlier point in the path leads only to curtailment of some kind -- non-stochastic, stochastic, or both --
  # then the CP of that earlier point is essentially zero. There is no point in "reaching" it.
  # Similarly, if an earlier point in the path *may* lead to curtailment of some kind, the CP of that earlier point must change.

  # With the above in mind, it seems that it is not possible to separate the calculation of CP from the comparison of CP against theta;
  # They must be undertaken together.


  ###### NEXT SECTION: CHECKING IF EACH CP<THETA0 OR CP>THETA1 ######

  # Begin at row r+1, i.e. where Sm=r, and at column n-1.
  # Proceed right to left then bottom to top, ignoring cases where CP=0 or CP=1.
  # As CP increases from right to left, if CP>=theta then can move on to next row (because all CPs will be >=theta).

  # The CPs in regions A and B have already been defined above; these are the points
  # for which 0 < CP < 1. Combine these regions and examine:

  #cp.sm <- c(cp.subset2.Sm, cp.subset1.Sm)
  #cp.m <- c(cp.subset2.m, cp.subset1.m)
  # ^ These are all the points that it is necessary to cycle through.

  # To reduce looping, identify the rows which contain a CP < theta0 OR CP > theta1:

  # theta.test <- function(Sm, m)
  #   {
  #   sum(mat[Sm+1, m] < theta0 | mat[Sm+1, m] > theta1)
  # }
  #
  # low.cp.rows <- mapply(theta.test, Sm=cp.sm, m=cp.m)
  # Because of how CP changes, propagating right to left and upwards, we only
  # need to identify the "lowest" row with CP < theta or CP > theta1, ie the row with greatest Sm
  # that contains a CP < theta or CP > theta1. This is the row where we begin. Note: This may be all rows, or even none!
  # Recall that row number = Sm-1

  # mat.old <- mat

  # We only need to act if there are CPs < theta0 or > theta1 -- otherwise, the matrix of CPs does not change.
  # if(sum(low.cp.rows)>0)    # There may be some cases where there are no CPs to be changed; i.e. no stopping due to stochastic curtailment
  # {
  #   begin.at.row <- max(which(low.cp.rows>0))
  #
  #   cp.changing.sm <- cp.sm[1:begin.at.row]
  #   cp.changing.m <- cp.m[1:begin.at.row]
  #
  #
  #   # We calculate the CP, truncate to zero if CP<theta0 (or to 1 if CP > theta1), then move on:
  #
  # for(rown in rev(cp.changing.sm+1))
  #   {
  #   for(coln in rev(cp.changing.m[[rown]]))
  #     {
  #     currentCP <- q*mat[rown, coln+1] + p*mat[rown+1, coln+1]
  #     if(currentCP > theta1)   mat[rown, coln] <- 1   # If CP > theta1, amend to equal 1
  #     else  mat[rown, coln] <- ifelse(test = currentCP < theta0, yes=0, no=currentCP)  # Otherwise, test if CP < theta0. If so, amend to 0, otherwise calculate CP as normal
  #     }
  #   } # Again, this *must* be done one entry at a time -- cannot vectorise.
  #
  #
  # } # End of IF statement
  #
  # mat

  ###### STOP if design is pointless, i.e either failure or success is not possible from the beginning:
  # if(sum(mat[,1], na.rm = T)==2) stop("Design guarantees success")
  # if(sum(mat[,1], na.rm = T)==0) stop("Design guarantees failure")



  ###### At this point, the matrix "mat" contains the CPs, adjusted for stochastic curtailment.
  ###### The points in the path satisfying 0 < CP < 1 are the possible points for this design;
  ###### All other points are either terminal (i.e. points of curtailment) or impossible.
  ###### We now need the characteristics of this design: Type I error, expected sample size, and so on.



  pascal.list <- list(1, c(1,1))

  # for(i in 3:(n+2))
  # {
  # column <- as.numeric(mat[!is.na(mat[,i-2]), i-2])
  # CPzero.or.one <- which(column==0 | column==1)
  # newnew <- pascal.list[[i-1]]
  # newnew[CPzero.or.one] <- 0
  # pascal.list[[i]] <- c(0, newnew) + c(newnew, 0)
  # }

  for(i in 3:(n+2))
  {
    column <- mat[!is.na(mat[,i-2]), i-2]
    CPzero.or.one <- which(column!=0.5)
    newnew <- pascal.list[[i-1]]
    newnew[CPzero.or.one] <- 0
    pascal.list[[i]] <- c(0, newnew) + c(newnew, 0)
  }


  pascal.list <- pascal.list[c(-1, -length(pascal.list))]


  # Now obtain the rest of the probability -- the p^b * q^c :

  # coeffs <- list()
  # coeffs.p0 <- list()
  #
  # for(i in 1:n){
  #     j <- 1:(i+1)
  #     coeffs[[i]] <- p^(j-1)*q^(i+1-j)
  #     coeffs.p0[[i]] <- p0^(j-1)*q0^(i+1-j)
  # }

  needed <- (r2+1):n
  coeffs2 <- p1^(r2+1)*q1^(needed-(r2+1))
  coeffs2.p0 <- p0^(r2+1)*q0^(needed-(r2+1))

  # We only want the (r2+2)th element of each list (equivalent to Sm=r2+1, as element 1 is Sm=0), from m=r2+1:
  pascal.element.r2plus1 <- sapply(pascal.list, "[", (r2+2))[(r2+1):n]
  #pascal.element.r2plus1*coeffs2
  # single <- lapply(pascal.list, function(x) {x[r2+2]}) # Another way of getting the (r2+1)th element of each list.


  # Multiply the two triangles (A and p^b * q^c):
  #final.probs <- Map("*", pascal.list, coeffs)


  # for finding type I error prob:
  #final.probs.p0 <- Map("*", pascal.list, coeffs.p0)



  ###### We have the probability of each path, taking into account stochastic and NS curtailment. ######
  ###### We now must tabulate these paths.

  #final.probs.mat <- matrix(unlist(lapply(final.probs, '[', 1:max(sapply(final.probs, length)))), ncol = n, byrow = F)
  #rownames(final.probs.mat) <- 0:n

  #final.probs.mat.p0 <- matrix(unlist(lapply(final.probs.p0, '[', 1:max(sapply(final.probs.p0, length)))), ncol = n, byrow = F)
  #rownames(final.probs.mat.p0) <- 0:n

  # Find successful probabilities first:
  # m.success <- (r2+1):n
  # Sm.success <- rep(r2+1, length(m.success))
  # prob.success <- final.probs.mat[r2+2, m.success]
  # prob.success.p0 <- final.probs.mat.p0[r2+2, m.success]
  # success.deets <- cbind(Sm.success, m.success, prob.success, prob.success.p0)

  # FIRST: Search for terminal points of success. These can only exist in rows where (updated) CP=1, and where Sm<=r+1:
  #potential.success.rows <- rowSums(mat[1:(r2+2), ]==1, na.rm = TRUE)
  #rows.with.cp1 <- which(as.numeric(potential.success.rows)>0)
  #rows.with.cp1 <- r2+1+1
  # ^ These are the rows containing possible terminal points of success. They must have CP=1:

  #columns.of.rows.w.cp1 <- list()
  #columns.of.rows.w.cp1 <- which(mat[rows.with.cp1, ]==1 & !is.na(mat[rows.with.cp1, ]))


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


  # success <- NULL
  #
  # for(j in columns.of.rows.w.cp1)
  #   {
  #     if(mat[rows.with.cp1-1, j-1] < 1) success <- rbind(success, c(rows.with.cp1-1, j, final.probs.mat[rows.with.cp1, j], final.probs.mat.p0[rows.with.cp1, j]))
  #   }

  #success.n <- (r2+1):n
  #success.Sm <- rep(r2+1, length(columns.of.rows.w.cp1))
  #success.prob <- final.probs.mat[rows.with.cp1, columns.of.rows.w.cp1]
  #success.prob.p0 <- final.probs.mat.p0[rows.with.cp1, columns.of.rows.w.cp1]

  #success <- cbind(success.Sm, success.n, success.prob, success.prob.p0)

  #colnames(success) <- c("Sm", "m", "prob", "prob.p0")

  # Now failure probabilities. Note that there is one failure probability in each row, and in that
  # row the failure probability is the one that has the greatest m (i.e. the "furthest right" non-zero entry):

  # Identify non-zero terms in each row:
  # m.fail <- rep(NA, r2+1)
  # prob.fail <- rep(NA, r2+1)
  # prob.fail.p0 <- rep(NA, r2+1)


  # for(i in 1:(r2+1))
  # {
  # m.fail[i] <- max(which(final.probs.mat[i ,]!=0))
  # #prob.fail[i] <- final.probs.mat[i, m.fail[i]]
  # prob.fail.p0[i] <- final.probs.mat.p0[i, m.fail[i]]
  # }
  # Sm.fail <- 0:r2


  #fail.deets <- cbind(Sm.fail, m.fail, prob.fail, prob.fail.p0)

  #output <- rbind(fail.deets, success)
  #output <- as.data.frame(output)
  #output$success <- c(rep("Fail", length(m.fail)), rep("Success", nrow(success)))
  #names(output) <- c("Sm", "m", "prob", "prob.p0", "success")



  ##################### Now find characteristics of design #####################


  #sample.size.expd <- sum(output$m*output$prob)

  #sample.size.expd.p0 <- sum(output$m*output$prob.p0)

  #alpha <- sum(output$prob.p0[output$success=="Success"])

  #power <- sum(output$prob[output$success=="Success"])

  #output <- c(n1=n1, n2=n2, n=n, r1=r1, r=r2, alpha=sum(success.prob.p0), power=sum(success.prob))
  output <- c(n1=n1, n2=n2, n=n, r1=r1, r=r2, alpha=sum(pascal.element.r2plus1*coeffs2.p0), power=sum(pascal.element.r2plus1*coeffs2))
  output
}

findN1N2R1R2NSC <- function(nmin, nmax)
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

  n1.list2 <- rev(unlist(n1.list))
  n2.list2 <- rev(unlist(n2.list))
  length(n1.list2)
  # There are approx. 3000 combinations of n1 and n2.

  ns <- cbind(n1.list2, n2.list2)

  n.list2 <- apply(ns, 1, sum)

  ns <- cbind(ns, n.list2)

  colnames(ns) <- c("n1", "n2", "n")



  ########################################################################################
  ################################ FIND COMBNS OF R1 AND R ###############################
  ########################################################################################
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


  r1 <- list()

  for(i in 1:nrow(ns))
  {
    r1[[i]] <- 0:(ns[i,"n1"]-1) # r1 values: 0 to n1-1
  }

  # How many combinations so far?
  length(unlist(r1)) # Approx. 85,000 designs, before r is varied.

  # For each possible combination of n1, n2 and r1, there is a combination of possible r's:

  r <- vector("list", (nrow(ns)))
  # ^ r will be a list of lists: e.g., r[[1]] will contain all possibilities for r for n1[1] and n2[1],
  # while r[[1]][[1]] will contain all possibilities for r for n1[1] and n2[1] AND r1[[1]][1]

  for(i in 1:nrow(ns))
  {
    #print(i)
    for(j in 1:length(r1[[i]]))
    {
      r[[i]][[j]] <- (r1[[i]][j]+1):(ns[i,"n"]-1)  # r must satisfy r > r1 and r < n

      keep.these <- r[[i]][[j]]-r1[[i]][j] <= ns[i,"n2"] # <- Note constraint: The number of responses required in stage 2 (r2-r1) must be at most n2:
      r[[i]][[j]] <- r[[i]][[j]][keep.these]
    }
  }


  #length(unlist(r)) # Approx. 300,000 designs for n [2, 50]
  # Approx. 200,000 designs for n [61, 65] alone

  ######### REMOVE UNNEEDED OBJECTS:
  rm(nposs, keep.these)
  return(list(ns=ns, n1=n1.list2, n2=n2.list2, n=n.list2, r1=r1, r=r))
}

############ Incorporating non-stochastic curtailment for both futility and benefit ############

# There are 4 regions:
#   - cp=0 (due to failure in S1 or failure in S2 or because path is impossible due to design (S2 only))
#   - cp=1 (in S1 or S2)
#   - cp where k > r1 (in S1 or S2)
#   - cp where k<= r1 but progression is still possible (S1 only)
#
# Obtain conditional power (CP) for these regions in this order: CP for the final region is obtained from the CP in the other regions.


findNSCdesignOCs <- function(n1, n2, r1, r2, p0=p0, p1=p1, theta0=0, theta1=1)
{

  # if(theta0 >= theta1) stop("theta1 must be greater than theta0")
  # if(p0 >= p) stop("p1 must be greater than p")
  # if(r1 >= r2) stop("r2 must be greater than r1")


  n <- as.numeric(n1+n2)
  q1 <- 1-p1
  q0 <- 1-p0

  # Start with all zeros (why not):
  mat <- matrix(c(rep(0, n*(r2+1)), rep(1, n*((n+1)-(r2+1)))), nrow = n+1, byrow=T)


  # Now input CP for points where r1 < k < r+1 (ie progression to Stage 2 is guaranteed):
  Sm <- 0:n # possible k~ values.
  cp.subset1.Sm <- which(r1 < Sm & Sm < (r2+1)) - 1  # Values for Sm where r1 < Sm < r2+1, i.e. progression to S2 certain, but overall trial success not guaranteed.




  # Which columns? Trial success must still be possible, ie  n-m > r-Sm+1
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


  ##### Now calculate CP for all these points: #####

  # For each Sm, begin with the earliest point, ie lowest m, and work row-wise to the latest point, ie highest m:
  # Remember that row number is Sm+1.
  # How many rows to do? length(cp.subset1.m):

  # summation <- list()
  #
  # for(j in 1:length(cp.subset1.m))
  # {
  #   current.Sm <- cp.subset1.Sm[j]
  #
  #   l <- seq(from=r2-current.Sm, length=length(cp.subset1.m[[j]]))
  #
  #   summation[[j]] <- choose(l, r2-current.Sm) * p^(r2-current.Sm+1) * q^(l-(r2-current.Sm))
  # }
  #
  # cp1 <- lapply(summation, cumsum)
  # cp1 <- lapply(cp1, rev)
  #
  #
  # # Now, insert CPs into matrix of CPs:
  #
  # for(i in 1:length(cp.subset1.Sm))  {  mat[cp.subset1.Sm[i]+1, cp.subset1.m[[i]]] <- cp1[[i]] }
  #
  #
  # # Final region to be addressed: Points where k <= r1 but progression is still possible (S1 only):
  cp.subset2.Sm <- 0:r1  # Values of Sm <= r1
  #
  # # Which columns? Progression to S2 must still be possible, ie  n1-m > r1-Sm+1
  n.to.n1 <- 1:n1 # possible S1 values

  cp.subset2.m <- list()
  i <- 1
  for(k in cp.subset2.Sm)
  {
    cp.subset2.m[[i]] <- which(n1-n.to.n1 >= r1-k+1 & n.to.n1 >= k)
    i <- i+1
  }

  # Combine all the non-zero, non-one CP points:
  cp.subset21.Sm <- c(cp.subset2.Sm, cp.subset1.Sm)
  cp.subset21.m <- c(cp.subset2.m, cp.subset1.m)

  for(i in 1:length(cp.subset21.Sm))
  {
    mat[cp.subset21.Sm[i]+1, cp.subset21.m[[i]]] <- 0.5
  }

  # Now we have the rows and columns of this region:

  # cp.subset2.Sm # Values of Sm. Rows are this plus 1.
  # cp.subset2.m # Values of m (columns)

  # Have to propagate "backwards", ending at Sm=0, n=1.
  # As with previous region, proceed row-wise -- start with greatest Sm.

  # for(k in length(cp.subset2.Sm):1) # Note that we END at the earliest row.
  # {
  #   for(j in rev(cp.subset2.m[[k]]))
  #     {
  #     mat[k, j] <- q*mat[k, j+1] + p*mat[k+1, j+1]
  #     }
  # }
  # This *must* be done in this way -- cannot be vectorised.

  # Create the "blanks":
  # for(i in 1:(ncol(mat)-1))
  # {
  #   mat[(i+2):nrow(mat), i] <- NA
  # }
  NAs <- rbind(FALSE, lower.tri(mat)[-nrow(mat),])
  mat[NAs] <- NA

  # all.thetas <- unique(c(mat))
  # subset.thetas <- all.thetas[all.thetas < 0.2]

  # We now have a matrix m containing the CPs in (more or less) the upper triangle.


  ######################
  ######IMPORTANT ######
  ######################

  # If a point in a path has CP < theta, we stop the trial due to stochastic curtailment.
  # If an earlier point in the path leads only to curtailment of some kind -- non-stochastic, stochastic, or both --
  # then the CP of that earlier point is essentially zero. There is no point in "reaching" it.
  # Similarly, if an earlier point in the path *may* lead to curtailment of some kind, the CP of that earlier point must change.

  # With the above in mind, it seems that it is not possible to separate the calculation of CP from the comparison of CP against theta;
  # They must be undertaken together.


  ###### NEXT SECTION: CHECKING IF EACH CP<THETA0 OR CP>THETA1 ######

  # Begin at row r+1, i.e. where Sm=r, and at column n-1.
  # Proceed right to left then bottom to top, ignoring cases where CP=0 or CP=1.
  # As CP increases from right to left, if CP>=theta then can move on to next row (because all CPs will be >=theta).

  # The CPs in regions A and B have already been defined above; these are the points
  # for which 0 < CP < 1. Combine these regions and examine:

  cp.sm <- c(cp.subset2.Sm, cp.subset1.Sm)
  cp.m <- c(cp.subset2.m, cp.subset1.m)
  # ^ These are all the points that it is necessary to cycle through.

  # To reduce looping, identify the rows which contain a CP < theta0 OR CP > theta1:

  # theta.test <- function(Sm, m)
  #   {
  #   sum(mat[Sm+1, m] < theta0 | mat[Sm+1, m] > theta1)
  # }
  #
  # low.cp.rows <- mapply(theta.test, Sm=cp.sm, m=cp.m)
  # Because of how CP changes, propagating right to left and upwards, we only
  # need to identify the "lowest" row with CP < theta or CP > theta1, ie the row with greatest Sm
  # that contains a CP < theta or CP > theta1. This is the row where we begin. Note: This may be all rows, or even none!
  # Recall that row number = Sm-1

  # mat.old <- mat

  # We only need to act if there are CPs < theta0 or > theta1 -- otherwise, the matrix of CPs does not change.
  # if(sum(low.cp.rows)>0)    # There may be some cases where there are no CPs to be changed; i.e. no stopping due to stochastic curtailment
  # {
  #   begin.at.row <- max(which(low.cp.rows>0))
  #
  #   cp.changing.sm <- cp.sm[1:begin.at.row]
  #   cp.changing.m <- cp.m[1:begin.at.row]
  #
  #
  #   # We calculate the CP, truncate to zero if CP<theta0 (or to 1 if CP > theta1), then move on:
  #
  # for(rown in rev(cp.changing.sm+1))
  #   {
  #   for(coln in rev(cp.changing.m[[rown]]))
  #     {
  #     currentCP <- q*mat[rown, coln+1] + p*mat[rown+1, coln+1]
  #     if(currentCP > theta1)   mat[rown, coln] <- 1   # If CP > theta1, amend to equal 1
  #     else  mat[rown, coln] <- ifelse(test = currentCP < theta0, yes=0, no=currentCP)  # Otherwise, test if CP < theta0. If so, amend to 0, otherwise calculate CP as normal
  #     }
  #   } # Again, this *must* be done one entry at a time -- cannot vectorise.
  #
  #
  # } # End of IF statement
  #
  # mat

  ###### STOP if design is pointless, i.e either failure or success is not possible from the beginning:
  # if(sum(mat[,1], na.rm = T)==2) stop("Design guarantees success")
  # if(sum(mat[,1], na.rm = T)==0) stop("Design guarantees failure")



  ###### At this point, the matrix "mat" contains the CPs, adjusted for stochastic curtailment.
  ###### The points in the path satisfying 0 < CP < 1 are the possible points for this design;
  ###### All other points are either terminal (i.e. points of curtailment) or impossible.
  ###### We now need the characteristics of this design: Type I error, expected sample size, and so on.



  pascal.list <- list(1, c(1,1))

  # for(i in 3:(n+2))
  # {
  # column <- as.numeric(mat[!is.na(mat[,i-2]), i-2])
  # CPzero.or.one <- which(column==0 | column==1)
  # newnew <- pascal.list[[i-1]]
  # newnew[CPzero.or.one] <- 0
  # pascal.list[[i]] <- c(0, newnew) + c(newnew, 0)
  # }

  for(i in 3:(n+2))
  {
    column <- mat[!is.na(mat[,i-2]), i-2]
    CPzero.or.one <- which(column!=0.5)
    newnew <- pascal.list[[i-1]]
    newnew[CPzero.or.one] <- 0
    pascal.list[[i]] <- c(0, newnew) + c(newnew, 0)
  }


  pascal.list <- pascal.list[c(-1, -length(pascal.list))]


  # Now obtain the rest of the probability -- the p^b * q^c :

  coeffs <- list()
  coeffs.p0 <- list()

  for(i in 1:n){
    j <- 1:(i+1)
    coeffs[[i]] <- p1^(j-1)*q1^(i+1-j)
    coeffs.p0[[i]] <- p0^(j-1)*q0^(i+1-j)
  }

  needed <- (r2+1):n
  coeffs2 <- p1^(r2+1)*q1^(needed-(r2+1))
  coeffs2.p0 <- p0^(r2+1)*q0^(needed-(r2+1))

  # We only want the (r2+2)th element of each list (equivalent to Sm=r2+1, as element 1 is Sm=0), from m=r2+1:
  pascal.element.r2plus1 <- sapply(pascal.list, "[", (r2+2))[(r2+1):n]
  # single <- lapply(pascal.list, function(x) {x[r2+2]}) # Another way of getting the (r2+1)th element of each list.


  # Multiply the two triangles (A and p^b * q^c), needed for Ess under H1:
  final.probs <- Map("*", pascal.list, coeffs)

  # for finding Ess under H0:
  final.probs.p0 <- Map("*", pascal.list, coeffs.p0)



  ###### We have the probability of each path, taking into account stochastic and NS curtailment. ######
  ###### We now must tabulate these paths.

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
  #potential.success.rows <- rowSums(mat[1:(r2+2), ]==1, na.rm = TRUE)
  #rows.with.cp1 <- which(as.numeric(potential.success.rows)>0)
  rows.with.cp1 <- r2+1+1
  # ^ These are the rows containing possible terminal points of success. They must have CP=1:
  columns.of.rows.w.cp1 <- (r2+1):n

  #columns.of.rows.w.cp1 <- list()
  #columns.of.rows.w.cp1 <- which(mat[rows.with.cp1, ]==1 & !is.na(mat[rows.with.cp1, ]))


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


  # success <- NULL
  #
  # for(j in columns.of.rows.w.cp1)
  #   {
  #     if(mat[rows.with.cp1-1, j-1] < 1) success <- rbind(success, c(rows.with.cp1-1, j, final.probs.mat[rows.with.cp1, j], final.probs.mat.p0[rows.with.cp1, j]))
  #   }


  success.n <- (r2+1):n
  success.Sm <- rep(r2+1, length(columns.of.rows.w.cp1))
  #success.prob <- final.probs.mat[rows.with.cp1, columns.of.rows.w.cp1]
  success.prob <- pascal.element.r2plus1*coeffs2
  #success.prob.p0 <- final.probs.mat.p0[rows.with.cp1, columns.of.rows.w.cp1]
  success.prob.p0 <-pascal.element.r2plus1*coeffs2.p0

  success <- cbind(success.Sm, success.n, success.prob, success.prob.p0)


  colnames(success) <- c("Sm", "m", "prob", "prob.p0")

  # Now failure probabilities. Note that there is one failure probability in each row, and in that
  # row the failure probability is the one that has the greatest m (i.e. the "furthest right" non-zero entry):

  # Identify non-zero terms in each row:
  m.fail <- rep(NA, r2+1)
  prob.fail <- rep(NA, r2+1)
  prob.fail.p0 <- rep(NA, r2+1)


  for(i in 1:(r2+1))
  {
    m.fail[i] <- max(which(final.probs.mat[i ,]!=0))
    prob.fail[i] <- final.probs.mat[i, m.fail[i]]
    prob.fail.p0[i] <- final.probs.mat.p0[i, m.fail[i]]
  }
  Sm.fail <- 0:r2


  fail.deets <- cbind(Sm.fail, m.fail, prob.fail, prob.fail.p0)

  output <- rbind(fail.deets, success)
  rownames(output) <- NULL
  output <- as.data.frame(output)
  output$success <- c(rep("Fail", length(m.fail)), rep("Success", nrow(success)))
  names(output) <- c("Sm", "m", "prob", "prob.p0", "success")



  #print(output)


  ##################### Now find characteristics of design #####################


  sample.size.expd <- sum(output$m*output$prob)

  sample.size.expd.p0 <- sum(output$m*output$prob.p0)

  #alpha <- sum(output$prob.p0[output$success=="Success"])

  #power <- sum(output$prob[output$success=="Success"])

  #output <- c(n1=n1, n2=n2, n=n, r1=r1, r=r2, alpha=sum(success.prob.p0), power=sum(success.prob))
  output <- c(n1=n1, n2=n2, n=n, r1=r1, r=r2, alpha=sum(success.prob.p0), power=sum(success.prob), EssH0=sample.size.expd.p0, Ess=sample.size.expd)
  output
}
#
# system.time({
# for(i in 1:1000) nsc.power.alpha.ess(n1=20, n2=20, r1=5, r2=10, p=0.3, p0=0.1)
# })
#
#   n1=20; n2=20; r1=5; r2=10; p=0.3; p0=0.1
# warnings()


