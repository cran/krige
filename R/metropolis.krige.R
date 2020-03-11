# SAMPLER
metropolis.krige <- function(y,X,east,north,mcmc.samples=100,spatial.share=0.5,range.share=0.5,beta.var=10,a.tune=10000){
	
  #DEFINE STANDING PARAMETERS
  dist.mat <- as.matrix( dist(cbind(east,north), method = "euclidean", diag = TRUE, upper = TRUE) )
  ideo.out <- lm(y ~ X[,-1])
  err.var<-var(resid(ideo.out))
  s.tau2.shape<-1+1/(1-spatial.share)
  s.sigma2.shape<-1+1/spatial.share
  s.phi.shape<-1+1/range.share
  s.phi.rate<-max(dist.mat)
  s2<-0.05
  
  # STARTING VALUES
  mcmc.mat <- matrix(NA,nrow=mcmc.samples,ncol=3+ncol(X))
  dimnames(mcmc.mat)[[2]] <- c("tau2","phi","sigma2",dimnames(X)[[2]])
  n.vars <- ncol(X)
  accepted <- 0
  set.seed(runif(1))
  mcmc.mat[1,] <- c((1-spatial.share)*err.var,range.share*s.phi.rate,spatial.share*err.var,ideo.out$coef)
  
  # START OF MARKOV CHAIN
  for (i in 1:(nrow(mcmc.mat)-1))  {
    if(i%%100==0){cat("Iteration ",i,": Nugget=",mcmc.mat[i,1],". Partial sill=",mcmc.mat[i,3],".\n",sep="")}
    current <- krige.posterior(mcmc.mat[i,1],mcmc.mat[i,2],mcmc.mat[i,3],mcmc.mat[i,4:ncol(mcmc.mat)],
    		y,X,east,north,p.spatial.share=spatial.share,p.range.share=range.share,p.beta.var=beta.var,tot.var=err.var)
    tau2 <- rinvgamma(1,shape=s.tau2.shape,rate=err.var)
    sigma2 <- rinvgamma(1,shape=s.sigma2.shape,rate=err.var)
    phi <-  rinvgamma(1,shape=s.phi.shape,rate=s.phi.rate)
    #beta <- mvrnorm(n = 1, mu=mcmc.mat[i,4:ncol(mcmc.mat)], diag(diag(vcov(ideo.out)))^s2)
    #beta <- mvrnorm(n = 1, mu=mcmc.mat[i,4:ncol(mcmc.mat)], s2*diag(diag(vcov(ideo.out))))
    beta <- mvrnorm(n = 1, mu=mcmc.mat[i,4:ncol(mcmc.mat)], s2*vcov(ideo.out))
    candidate <- krige.posterior(tau2,phi,sigma2,beta,
    		y,X,east,north,p.spatial.share=spatial.share,p.range.share=range.share,p.beta.var=beta.var,tot.var=err.var)
    a.ratio <- ifelse(candidate>current,1.1,(exp(candidate/a.tune)/exp(current/a.tune))^a.tune) 
    if (is.na(a.ratio) == TRUE) { print(paste("is.na:",a.ratio)); a.ratio <- 0 }
    if (a.ratio > runif(1)) {
      accepted <- accepted + 1
      mcmc.mat[(i+1),] <- c(tau2,phi,sigma2,beta)
    }
    else (mcmc.mat[(i+1),] <- mcmc.mat[i,])
  }
  print(paste(accepted,"accepted out of",nrow(mcmc.mat)-1,"generated", round(100*accepted/(nrow(mcmc.mat)-1)), "percent"))
  return(mcmc.mat)
}
