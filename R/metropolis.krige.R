# SAMPLER
metropolis.krige <- function(y,X,east,north,mcmc.samples=100,tau2.shape=5000,tau2.rate=10,phi.shape=5,phi.rate=2000,sigma2.shape=1600,sigma2.rate=10,a.tune=10000)  {
  dist.mat <- as.matrix( dist(cbind(east,north), method = "euclidean", diag = TRUE, upper = TRUE) )
  ideo.out <- lm(y ~ X[,-1])
  s2 <- summary(ideo.out)$coef[,2]^2/3
  # STARTING VALUES
  mcmc.mat <- matrix(NA,nrow=mcmc.samples,ncol=3+ncol(X))
  dimnames(mcmc.mat)[[2]] <- c("tau2","phi","sigma2",dimnames(X)[[2]])
  n.vars <- ncol(X)
  accepted <- 0
  set.seed(runif(1))
  mcmc.mat[1,] <- c(1,1,1,rep(0,n.vars))
  # START OF MARKOV CHAIN
  for (i in 1:(nrow(mcmc.mat)-1))  {
    if(i%%100==0){cat("Iteration ",i,": Nugget=",mcmc.mat[i,1],". Partial sill=",mcmc.mat[i,3],".\n",sep="")}
    current <- krige.posterior(mcmc.mat[i,1],mcmc.mat[i,2],mcmc.mat[i,3],mcmc.mat[i,4:ncol(mcmc.mat)],y,X,east,north)
    tau2 <- rgamma(1,shape=tau2.shape,rate=tau2.rate)
    phi <-  rgamma(1,shape=phi.shape,rate=phi.rate)
    sigma2 <- rgamma(1,shape=sigma2.shape,rate=sigma2.rate)
    beta <- mvrnorm(n = 1, mu=mcmc.mat[i,4:ncol(mcmc.mat)], s2*diag(n.vars))
    candidate <- krige.posterior(tau2,phi,sigma2,beta,y,X,east,north)
    a.ratio <- (exp(candidate/a.tune)/exp(current/a.tune))
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


