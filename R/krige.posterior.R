# CREATING THE LOG-POSTERIOR OBJECTIVE FUNCTION
krige.posterior <- function(tau2,phi,sigma2,beta,y,X,east,north,p.spatial.share=0.5,p.range.share=0.5,p.beta.var=10,tot.var=var(y)){ 
  distance <- as.matrix(dist(cbind(east,north), method = "euclidean", diag = TRUE, upper = TRUE) )
  tau2.shape<-1+1/(1-p.spatial.share)
  sigma2.shape<-1+1/p.spatial.share
  phi.shape<-1+(1/p.range.share)
  phi.rate<-max(distance)
  Sigma <- ifelse(distance>0, sigma2*exp(-(phi^2)*(distance^2)), tau2+sigma2) # GAUSSIAN VAR/COV
  mu<-X%*%beta
  log.det.Sigma <- sum(log(eigen(Sigma)$values))
  chol.Sigma <- chol(Sigma)
  loglik <- -(1/2)*log.det.Sigma - 0.5*t(y-mu)%*%chol2inv(chol.Sigma)%*%(y-mu) 
  log.mu.prior <- sum(log(dmvnorm(beta, mean=rep(0,length(beta)), sigma = diag(length(beta))*p.beta.var)))
  log.tau2.prior <-dinvgamma(tau2,shape=tau2.shape,rate=tot.var,log=T) 
  log.sigma2.prior <-dinvgamma(sigma2,shape=sigma2.shape,rate=tot.var,log=T)
  log.phi.prior <- dinvgamma(phi,shape=phi.shape,rate=phi.rate,log=T)
  logpost <- loglik+log.tau2.prior+log.sigma2.prior+log.phi.prior+log.mu.prior
  return(logpost)
}
