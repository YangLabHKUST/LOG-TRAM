estimate_Oemga <- function(Beta,Sigma,verbose=T){
  p <- nrow(Beta)
  OmegaHat0 <- t(Beta) %*% Beta / p - Sigma
  omega1 <- OmegaHat0[1,1]
  omega2 <- OmegaHat0[2,2]
  if(omega1<1e-30|omega2<1e-30){
    rhog <- 0
  } else {
    rhog <- OmegaHat0[1,2]/sqrt(omega1*omega2)
  }
  rhog <- sign(rhog)*min(0.99,abs(rhog))
  
  omega1 <- max(omega1,1e-30)
  omega2 <- max(omega2,1e-30)
  
  OmegaHat <- matrix(c(omega1,rhog*sqrt(omega1*omega2),rhog*sqrt(omega1*omega2),omega2),2,2)
  
  return(list(OmeOmegaHatga=OmegaHat,OmegaHat0=OmegaHat0))
}

meta_blup <- function(Omega,Sigma,beta){
  # K: number of pop; p: number of SNPs
  # Omega: K by K by p per-SNP ld weighted genetic covariance
  # Sigms: K by K by p SNP-specific error covariance
  # beta: p by K vector of marginal GWAS effect size
  
  K <- ncol(beta)
  p <- nrow(beta)
  
  betahat <- matrix(0,p,K)
  varMat_beta <- array(0,dim = c(K,K,p))
  for(j in 1:p){
    varMat_beta[,,j] <- solve(solve(Omega[,,j])+solve(Sigma[,,j]))
    betahat[j,] <- varMat_beta[,,j]%*%solve(Sigma[,,j],beta[j,])
  }
  var_beta <- t(apply(varMat_beta,3,diag))
  z_score <- betahat/sqrt(var_beta)
  pval <- 2*pnorm(abs(z_score),lower.tail = F)
  return(list(beta=betahat,var_beta=var_beta,varMat_beta=varMat_beta,z_score=z_score,pval=pval))
}
