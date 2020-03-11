% File src/library/base/man/load.Rd
     \name{metropolis.krige}
     \alias{metropolis.krige}
     \title{Sampling Technique Using Metropolis-Hastings}
     \description{
       This function performs Metropolis-Hastings sampling for a linear model specified over point-referenced geospatial data. It returns MCMC iterations, with which results of the geospatial linear model can be summarized.       
     }
     \usage{
metropolis.krige(y,X,east,north,mcmc.samples=100,
     spatial.share=0.5,range.share=0.5,beta.var=10,a.tune=10000)
     }
     \arguments{
            \item{y}{The dependent variable that is used in the kriging model.}
             \item{X}{The matrix of independent variables used in the kriging model.}
              \item{east}{Vector of eastings for all observations.}
               \item{north}{Vector of northings for all observations.}
               \item{mcmc.samples}{Number of MCMC iterations.}
               \item{spatial.share}{Prior for proportion of unexplained variance that is spatial in nature. Must be greater than 0 and less than 1. Defaults to an even split.}
               \item{range.share}{Prior for the range term, as a proportion of the maximum distance in the data. Defaults to half the maximum distance.}
               \item{beta.var}{Prior for the variance on zero-meaned normal priors on the regression coefficients. Defaults to 10.}
               \item{a.tune}{Tuning parameter for the acceptance ratio of the Metropolis algorithm. Defaults to 10000.}
     }
     \details{
      Analysts should use this function if they want to estimate a linear regression model in which each observation can be located at points in geographic space. That is, each observation is observed for a set of coordinates in eastings & northings or longitude & latitude. 
      
      Researchers must specify their model in the following manner: \code{y} should be a column vector of the dependent variable. \code{X} should be a matrix that includes all independent variables in the model, including a constant vector to estimate an intercept term. \code{east} should be a vector of all west-east coordinates for observations (ideally eastings but possibly longitude). \code{north} should be a vector of all north-south coordinates for observations (ideally northings but possibly latitude). \code{mcmc.samples} is the number of iterations to sample from the posterior distribution using the Metropolis-Hastings algorithm. This defaults to 100 iterations, but many more iterations would normally be preferred.
      
      The returned value is a matrix of sampled values from the posterior distribution. Rows represent the iteration number, and ideally the user will discard the first several rows as burn-in. Columns represent the parameter, so summarizing the matrix by column offers summaries of the model's results.
     }
     \value{
     Returns an object of class "matrix" that includes all iterations of the Metropolis-Hastings sampler. Each column of the matrix represents a different parameter--starting with \code{tau2}, \code{phi}, and \code{sigma2} before listing regression coefficients. Each row represents another iteration of the MCMC sampler. Summarizing the matrix by column offers summaries of the marginal posterior distribution by parameter.
     }
     
     \references{
      James E. Monogan III & Jeff Gill. 2016. "Measuring State and District Ideology with Spatial Realignment." \emph{Political Science Research and Methods} 4(1):97-121.
     }

     \examples{
#Data
ny <- NYcity_subset

# Explanatory Variable Matrix
psrm.data <-cbind(1, ny$age, ny$educ, I(ny$age*ny$educ), as.numeric(ny$race==2), 
     as.numeric(ny$race==3), ny$female, I(as.numeric(ny$race==2)*ny$female), 
     I(as.numeric(ny$race==3)*ny$female), ny$catholic, ny$mormon, ny$orthodox, 
     ny$jewish, ny$islam, ny$mainline, ny$evangelical, ny$rural, 
     ny$ownership, as.numeric(ny$empstat==2), as.numeric(ny$empstat==3),
     as.numeric(ny$inc14==2), as.numeric(ny$inc14==3), as.numeric(ny$inc14==4),
     as.numeric(ny$inc14==5), as.numeric(ny$inc14==6), as.numeric(ny$inc14==7),
     as.numeric(ny$inc14==8), as.numeric(ny$inc14==9), as.numeric(ny$inc14==10),
     as.numeric(ny$inc14==11), as.numeric(ny$inc14==12), as.numeric(ny$inc14==13),
     as.numeric(ny$inc14==14))
dimnames(psrm.data)[[2]] <- c("Intercept", "Age", "Education", "Age.education", 
     "African.American", "Nonwhite.nonblack","Female", 
     "African.American.female", "Nonwhite.nonblack.female", 
     "Income1", "Income2", "Income3", "Income4","Income5", 
     "Income6", "Income7", "Income8", "Income9", "Income10", 
     "Income11", "Income12", "Income13", "Catholic", "Mormon", 
     "Orthodox", "Jewish", "Islam", "Mainline", "Evangelical", 
     "Ruralism", "Homeowner", "Unemployed","Not.in.workforce")

# Outcome Variable
ideo <- matrix(ny$ideology,ncol=1)

# Set Number of Iterations (5 for Illustration, But Want Many More)
M<-5

# Trial Run, Linear Model of Ideology with Geospatial Errors Using Metropolis-Hastings:
test.mat <- metropolis.krige(y=ideo,X=psrm.data,east=ny$eastings,north=ny$northings,
     mcmc.samples=M,spatial.share=.25,range.share=.9) 

# Summarize Results
summary(test.mat)

\donttest{
# This second example uses 100 iterations instead. 
# The command takes much longer, but many more than 100 are desired.
# Set Number of Iterations: 
M<-100

# Estimate Linear Model of Ideology with Geospatial Errors Using Metropolis-Hastings
# Again, even 100 iterations is short:
out.mat <- metropolis.krige(y=ideo,X=psrm.data,east=ny$eastings,north=ny$northings,
     mcmc.samples=M,spatial.share=.25,range.share=.9) 

# Discard first 10\% of Iterations as Burn-In (User Discretion Advised).
out.mat <- out.mat[(.1*M+1):M,]

# Summarize Results
summary(out.mat)
}
}
\keyword{estimation}
\keyword{posterior}

