% File src/library/base/man/load.Rd
     \name{krige.posterior}
     \alias{krige.posterior}
     \title{Posterior Distribution for the Kriging Process}
     \description{
       This function finds the posterior density of a geospatial linear regression model given a point-referenced geospatial dataset and a set of parameter values. The function is useful for finding the optimum of or for sampling from the posterior distribution.
     }
     \usage{
krige.posterior(tau2,phi,sigma2,beta,y,X,east,north)
     }
     \arguments{
        \item{tau2}{Value of the nugget, or non-spatial error variance.}
         \item{phi}{Value of the decay term, driving the level of spatial correlation.}
          \item{sigma2}{Value of the partial sill, or maximum spatial error variance.}
           \item{beta}{Coefficients from linear model.}
            \item{y}{The dependent variable that is used in the kriging model.}
             \item{X}{The matrix of independent variables used in the kriging model.}
              \item{east}{Vector of eastings for all observations.}
               \item{north}{Vector of northings for all observations.}
          
     }
     \details{
      This function finds the posterior densisty for a kriging model. The function utilizes information provided about the parameters \code{tau2}, \code{phi}, \code{sigma2}, and \code{beta}. It also utilizes the observed data \code{y}, \code{X}, \code{east}, and \code{north}. Given a set of parameter values as well as the observed data, the function returns the posterior density for the specified model.
     }
     \value{
     Returns a single number that is the posterior density of the function, which is stored in object of class "matrix."
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

# Find the posterior density for New York City data if all parameters were 1:
s.test <- krige.posterior(tau2=1, phi=1, sigma2=1, beta=rep(1,33), y=ideo,
     X=psrm.data,east=ny$eastings, north=ny$northings)

# Print posterior density
s.test
}
\keyword{posterior}

