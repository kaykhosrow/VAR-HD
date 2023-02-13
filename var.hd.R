var.historical.decomposition <- function(Estimation){
  
  ## Retrieve and initialize variables 
  
  invA    <- t(chol(as.matrix(summary(Estimation)$covres)))   # inverse of the A matrix
  Fcomp   <- companionmatrix(Estimation)                      # Companion matrix
  eps     <- ginv(invA) %*% t(residuals(Estimation))          # structural errors 
  nvar    <- Estimation$K                                     # number of endogenous variables
  nvarXeq <- nvar * Estimation$p                              # number of lagged endogenous per equation
  nvar_ex <- 0                                                # number of exogenous (excluding constant and trend)
  Y       <- VARmakexy(Estimation$y,Estimation$p,1)$Y         # Left-hand side
  nobs    <- nrow(Y)                                          # number of observations
  
  ## Compute historical decompositions
  
  invA_big <- matrix(0,nvarXeq,nvar)
  invA_big[1:nvar,] <- invA
  Icomp <- cbind(diag(nvar), matrix(0,nvar,(Estimation$p-1)*nvar))
  HDshock_big <- array(0, dim=c(Estimation$p*nvar,nobs+1,nvar))
  HDshock <- array(0, dim=c(nvar,(nobs+1),nvar))
  
  for (j in 1:nvar){  
    eps_big <- matrix(0,nvar,(nobs+1)) 
    eps_big[j,2:ncol(eps_big)] <- eps[j,]
    for (i in 2:(nobs+1)){
      HDshock_big[,i,j] <- invA_big %*% eps_big[,i] + Fcomp %*% HDshock_big[,(i-1),j]
      HDshock[,i,j] <-  Icomp %*% HDshock_big[,i,j]
    } 
  } 
  
  HD.shock <- array(0, dim=c((nobs+Estimation$p),nvar,nvar)) 
  
  for (i in 1:nvar){
    for (j in 1:nvar){
      HD.shock[,j,i] <- c(rep(NA,Estimation$p), HDshock[i,(2:dim(HDshock)[2]),j])
    }
  }
  
  return(HD.shock)
}
