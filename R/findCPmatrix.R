findCPmatrix <- function(n, r, Csize, p0, p1, minstop){
  q1 <- 1-p1
  mat <- matrix(3, ncol=n, nrow=max(3, min(r+Csize, n)+1, minstop+1)) #  nrow=r+Csize+1 unless number of stages equals 1, minimum of 3.
  rownames(mat) <- 0:(nrow(mat)-1)
  mat[(r+2):nrow(mat),] <- 1
  mat[1:(r+1),n] <- 0

  #pat.cols <- seq(n, 1, by=-Csize)[-1]
  pat.cols <- seq(n, minstop, by=-Csize)[-1]

  for(i in (r+1):1){
    for(j in pat.cols){  # Only look every C patients (no need to look at final col)
      if((r+1-(i-1) > n-j)) { # IF success is not possible [Total responses needed (r+1) - responses so far (i-1)] > [no. of patients remaining (n-j)], THEN
        mat[i,j] <- 0
      }else{
        if(i-1<=j){ # Condition: Sm<=m
          newcp <- sum(choose(Csize,0:Csize)*p1^(0:Csize)*q1^(Csize:0)*mat[i:(i+Csize),j+Csize])
          mat[i,j] <- newcp
        }
      }
    }
  }
for(i in 3:nrow(mat)) {
    mat[i, 1:(i-2)] <- NA
  }
mat
}
