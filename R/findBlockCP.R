 findBlockCP <- function(n, r, Bsize, pc, pt, thetaF, thetaE, minstop){
   pat.cols <- seq(from=2*n, to=minstop, by=-Bsize)[-1]
   qc <- 1-pc
   qt <- 1-pt
   prob.vec <- findProbVec(Bsize=Bsize,
                           pt=pt,
                           qt=qt,
                           pc=pc,
                           qc=qc)
   # CREATE UNCURTAILED MATRIX
   mat <- matrix(3, ncol=2*n, nrow=min(n+r+Bsize+2, 2*n+1))
   rownames(mat) <- 0:(nrow(mat)-1)
   mat[(n+r+2):nrow(mat),] <- 1
   mat[1:(n+r+1),2*n] <- 0 # Fail at end
   for(i in (n+r+1):1){
     for(j in pat.cols){  # Only look every Bsize patients (no need to look at final col)
       if(i-1<=j){ # Condition: Sm<=m
         #    browser()
         #    print(paste("Rows:", i:(i+Bsize), ", Columns: ", j+Bsize, sep=""))
         #    print(mat[i:(i+Bsize), j+Bsize])
         mat[i,j] <- ifelse(test=j-(i-1) > n-r+1, yes=0, no=sum(prob.vec*mat[i:(i+Bsize), j+Bsize]))
         # IF success is not possible (i.e. [total no. of pats-Xa+Ya-Xb] > n-r+1), THEN set CP to zero. Otherwise, calculate it based on "future" CPs.
       }
     }
   }
   for(i in 3:nrow(mat)){
     mat[i, 1:(i-2)] <- NA
   }
   #uncurt <- mat
   ### CREATE CURTAILED MATRIX
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
   return(mat)
 }

