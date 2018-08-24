#NEED TO HAVE THIS FUNCTION DEFINED
# CREATING THE LOG-POSTERIOR OBJECTIVE FUNCTION
krige.posterior <- function(tau2,phi,sigma2,beta,y,X,east,north){ 
  distance <- as.matrix(dist(cbind(east,north), method = "euclidean", diag = TRUE, upper = TRUE) )
  Sigma <- ifelse(distance>0, sigma2*exp(-(phi^2)*(distance^2)), tau2+sigma2) # GAUSSIAN VAR/COV
  mu<-X%*%beta
  log.det.Sigma <- sum(log(eigen(Sigma)$values))
  chol.Sigma <- chol(Sigma)
  loglik <- -(1/2)*log.det.Sigma - 0.5*t(y-mu)%*%chol2inv(chol.Sigma)%*%(y-mu) 
  # PRIOR LINE: by proportionality we leave out -(N/2)*log(2) -(N/2)*log(pi) 
  # define prior on linear coefficients
  log.mu.prior <- sum(log(dmvnorm(beta, mean=rep(0,length(beta)), sigma = diag(length(beta))*10)))
  logpost <- loglik-((2*(sigma2^2)+1)/sigma2)-((2*(tau2^2)+1)/tau2)+log.mu.prior
  # PRIOR LINE: Adds Inverse Gamma(1,1) priors for tau2 and sigma2. log.mu.prior captures the priors on the regrssion coefficients.
  # PRIOR LINE: Note that we are using a Uniform(0.01,0.99) prior on phi. PRIOR ON PHI. We left that out here as a proportionality 
  # issue. If you don't use the probit transform, then the space of sampling has to be restricted.
  return(logpost)
}