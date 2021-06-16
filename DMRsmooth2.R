
#In this version, h is a vector
LQfit <- function(data, h) {
   x <- data[ , 1]
   y <- data[ , 2]
   l <- length(x)
   beta0 <- rep( 0, l)
   beta1 <- rep( 0, l)
   for( i in 1: l) {
       x.reg <- x - x[i]
       w    <- dnorm(x - x[i], 0, h[i])
       fit  <- lm(y~x.reg + I( x.reg^ 2), weights = w)
       beta0[i] <- fit$coe[1]
       beta1[i] <- fit$coe[2]
   }
   beta <- cbind( beta0, beta1)
   return(beta)
}

AddCI <- function(data, h, beta) {
  x <- data[ , 1]
  y <- data[ , 2]
  l <- length(x)
  beta0 <- beta[ , 1]
  beta1 <- beta[ , 2]
 
  ##Test for equilibruim points
  w    <- rep( 0, l)
  diag <- rep( 0, l)
  upperCI <- rep( 0, l)
  lowerCI <- rep( 0, l)
  se  <- rep( 0, l)
  z   <- rep( 0, l)
  p   <- rep( 0 ,l)

  ##Estimate sigma^2
  for(i in 1:l){
    Xi  <- cbind(rep( 1, l), x - x[i], (x - x[i]) ^ 2)
    ker <- dnorm(Xi[,2], 0, h[i])
    A   <- matrix(0, ncol = 3, nrow = 3)
    A[1,1] <- sum(ker)
    A[1,2] <- ker %*%Xi[,2]
    A[2,1] <- A[1,2]
    A[1,3] <- ker %*%Xi[,3]
    A[2,2] <- A[1,3]
    A[3,1] <- A[1,3]
    A[2,3] <- ker %*%Xi[,2]^3
    A[3,2] <- A[2,3]
    A[3,3] <- ker %*%Xi[,3]^2
    B <- solve(A)[1,]
    C <- rbind(ker,ker*Xi[,2],ker*Xi[,3])
    wi<- B %*% C
    diag[i] <- wi[i]
    w <- rbind(w, wi)
  }
  w       <- w[2:(l + 1),]
  second  <- sum(w^2)
  first   <- 2*sum(diag)
  v       <- first-second
  vari    <- 1/(l-v)*sum((y-beta0)^2)

  for(i in 1:l) {
       X <- cbind( rep( 1, l), x - x[i], ( x - x[i]) ^ 2)
       kernel <- dnorm( X[, 2], 0, h[i])
       An <- matrix( 0, ncol = 3, nrow = 3)
       Bn <- matrix( 0, ncol = 3, nrow = 3)
       An[1,1] <- sum(kernel)/l
       An[1,2] <- kernel %*% X[,2]/l
       An[2,1] <- An[1,2]
       An[1,3] <- kernel %*% X[,3]/l
       An[2,2] <- An[1,3]
       An[3,1] <- An[1,3]
       An[2,3] <- kernel %*% X[,2] ^ 3 / l
       An[3,2] <- An[2,3]
       An[3,3] <- kernel %*% X[,3] ^ 2 / l
       kernel2 <- kernel ^ 2
       Bn[1,1] <- sum(kernel2)/l/l
       Bn[1,2] <- kernel2 %*% X[,2] / l / l
       Bn[2,1] <- Bn[1,2]
       Bn[1,3] <- kernel2 %*% X[,3]/l/l
       Bn[2,2] <- Bn[1,3]
       Bn[3,1] <- Bn[1,3]
       Bn[2,3] <- kernel2 %*% X[,2]^3/l/l
       Bn[3,2] <- Bn[2,3]
       Bn[3,3] <- kernel2 %*% X[,3] ^ 2 / l / l
       sol <- solve( An)
       temp <- sol%*%Bn%*%sol
       temp2 <- temp[2,2]
       se[i] <- sqrt(vari*temp2)
       z[i] <- abs(beta1[i]/se[i])
       p[i] <-(1-pnorm(z[i]))*2
       upperCI[i] <- beta1[i] + 1.96 * se[i]
       lowerCI[i] <- beta1[i] - 1.96 * se[i]
  }
  upperCI <- round(upperCI, 5)
  lowerCI <- round(lowerCI, 5)
  p   <- round(p, 5)
  CIp <- cbind(upperCI,lowerCI,p)
  return(CIp)
}


## window size calculation
WinCal <- function(idx,h0=0.5,loc){
  if(idx <= 70){
    j   <- 1:140 
    d   <- abs(loc[j]-loc[idx])
    h   <- (loc[max(j[order(d)][1:70])]-loc[min(j[order(d)][1:70])])/4 
  }else if(idx >= (length(loc)-70)){
    j   <- (length(loc)-140):(length(loc)) 
    d   <- abs(loc[j]-loc[idx]) 
    h   <- (loc[max(j[order(d)][1:70])]-loc[min(j[order(d)][1:70])])/4 
  }else{
    j   <- (idx-70):(idx+70) # candidate sites to be included in window
    d   <- abs(loc[j]-loc[idx]) # calculate the distance for the upper 70 and lower 70
    h   <- (loc[max(j[order(d)][1:70])]-loc[min(j[order(d)][1:70])])/4 # calculate the h get from the 
  }
  h <- max(h,h0)
  return(h)
}


