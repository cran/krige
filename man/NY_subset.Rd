% File src/library/base/man/load.Rd
     \docType{data}
     \name{NY_subset}
     \alias{NY_subset}
     \title{New York State CCES Respondents in 2008}

     \description{
       These data are a subset of the 2008 Cooperative Congressional Election Survey (CCES) Common Content. Only 1108 respondents from the state of New York are included, with predictors drawn from Monogan & Gill's (2016) model of self-reported ideology. The CCES data are merged with predictors on geographic location based on ZIP codes (from ArcGIS & TomTom) and county ruralism (from the USDA).
     }
     \usage{
data(NY_subset)
     }

     \format{
     	       The NY_subset dataset has 1108 observations and 26 variables.

     	\describe{
      \item{\code{state}}{The state abbreviation of the respondent's residence.}
       \item{\code{zip}}{The respondent's ZIP code.}
        \item{\code{age}}{The age of the respondent in years.}
         \item{\code{female}}{An indicator of whether the respondent is female.}
          \item{\code{ideology}}{The respondent's self-reported ideology on a scale of 0 (liberal) to 100 (conservative).}
           \item{\code{educ}}{The respondent's level of education. 
           0=No Highschool,
           1=High School Graduate,
           2=Some College,
           3=2-year Degree,
           4=4-year degree,
           5=Post-Graduate.}
             \item{\code{race}}{The respondent's race. 1=White, 2=African American, 3=Nonwhite & nonblack.}
               \item{\code{empstat}}{The respondent's employment status. 1=employed, 2=unemployed, 3=not in workforce.}
                 \item{\code{ownership}}{Indicator for whether the respondent owns his or her own home.}
                   \item{\code{inc14}}{The respondent's self reported income. 1=Less than $10,000, 2=$10,000-$14,999, 3=$15,000-$19,000, 4=$20,000-$24,999, 5=$25,000-$29,999, 6=$30,000-$39,999, 7=$40,000-$49,999, 8=$50,000-$59,999, 9=$60,000-$69,999, 10=$70,000-$79,999, 11=$80,000-$89,999, 12=$100,000-$119,999, 13=$120,000-$149,999, 14=$150,000 or more.}
                    \item{\code{catholic}}{Indicator for whether the respondent is Catholic.}
                     \item{\code{mormon}}{Indicator for whether the respondent is Mormon.}
                      \item{\code{orthodox}}{Indicator for whether the respondent is Orthodox Christian.}
                       \item{\code{jewish}}{Indicator for whether the respondent is Jewish.}
                        \item{\code{islam}}{Indicator for whether the respondent is Muslim.}
                         \item{\code{mainline}}{Indicator for whether the respondent is Mainline Christian.}
                          \item{\code{evangelical}}{Indicator for whether the respondent is Evangelical Christian.}
                           \item{\code{FIPS_Code}}{FIPS code of the repondent's state.}
                            \item{\code{rural}}{Nine-point USDA scale of the ruralism of each county, with 0 meaning the most urban and 8 meaning the most rural.}
                             \item{\code{zipPop}}{Indicates the population of the repondent's ZIP code.}
                              \item{\code{zipLandKM}}{Indicates the land area in square kilometers of the repondent's ZIP code.}
                               \item{\code{weight}}{Survey weights created by the CCES.}
                                \item{\code{cd}}{The congressional district the respondent resides in.}
                                 \item{\code{fipsCD}}{Index that fuses the state FIPS code in the first two digits and the congressional district number in the last two digits.}
                                  \item{\code{northings}}{Indicates the geographical location of the respondent in kilometer-based northings.}
                                   \item{\code{eastings}}{Indicates the geographical location of the respondent in kilometer-based eastings.}
	    }
     }
\examples{
# Descriptive Statistics
summary(NY_subset)

# Simple Linear Model of Ideology as a Function of Race and Sex
simple.mod <- lm(ideology~as.factor(race)+female, data=NY_subset)
summary(simple.mod)
}
      \source{
     Ansolabehere, Stephen. 2011. "CCES, Common Content, 2008." Ver. 4. \cr
     \url{http://hdl.handle.net/1902.1/14003}
     \cr\cr
     ArcGIS. 2012. "USA ZIP Code Areas." \cr
     \url{https://www.arcgis.com/home/item.html?id=8d2012a2016e484dafaac0451f9aea24}
     \cr\cr
     United States Department of Agriculture. 2013. "2013 Rural-Urban Continuum Codes." \cr
     \url{http://www.ers.usda.gov/data-products/rural-urban-continuum-codes.aspx}
     }
     \references{
      James E. Monogan III & Jeff Gill. 2016. "Measuring State and District Ideology with Spatial Realignment." \emph{Political Science Research and Methods} 4(1):97-121.
     }
\keyword{datasets}

