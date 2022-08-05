findBounds <- function(des, des.input){
   Bsize <- des.input$block
   mat <- findBlockCP(n=des$n.arm,
                      r=des$r,
                      Bsize=Bsize,
                      pc=des.input$pc,
                      pt=des.input$pt,
                      thetaF=des$thetaF,
                      thetaE=des$thetaE,
                      minstop=des$eff.minstop)
   interims <- seq(from=des$eff.minstop, to=2*des$n.arm, by=Bsize)
   #boundaries <- matrix(NA, nrow=2, ncol=ncol(mat)/Bsize)
   boundaries <- matrix(NA, nrow=2, ncol=length(interims))
   rownames(boundaries) <- c("lower", "upper")
   #interims <- seq(from=Bsize, to=ncol(mat), by=Bsize)
   colnames(boundaries) <- paste(interims)
   for(i in 1:length(interims)){
     j <- interims[i]
     lower <- if (any(mat[,j]==0, na.rm=TRUE) ) max(which(mat[,j]==0))-1 else -Inf
     upper <- if (any(mat[,j]==1, na.rm=TRUE) ) which.max(mat[,j])-1 else Inf
     # -1 terms to account for the fact that row 1 is equivalent to zero successes.
     boundaries[, i] <- c(lower, upper)
   }
   return(boundaries)
 }
