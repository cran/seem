\name{Weather}

\alias{rain.day}
\alias{markov.rain.seq}
\alias{markov.rain}
\alias{pet}

\title{Weather simulation}
\description{Rain simulation and PET calculations}

\usage{
rain.day(param)
markov.rain.seq(P, ndays)
markov.rain(P, ndays, amount.param, plot.out = T)
pet(eleva, Rad, Temp, RH, wind, plot.out = F)
}

\arguments{
  \item{param}{parameters: a list that contains the mean, standard deviation, skew coefficient, shape, and the distribution to be sampled}
  \item{P}{Markov transition probability matrix}
  \item{ndays}{number of days to simulate}
  \item{amount.param}{same as param above}
  \item{plot.out}{logical variable to decide whether to plot}
  \item{eleva}{elevation amsl}
  \item{Rad}{solar radiation}
  \item{Temp}{air temperature}
  \item{RH}{relative humidity}
  \item{wind}{wind speed}
}

\details{


rain.day: This function uses one argument, which is composed of parameters, to generate the amount of rain in a day and it is given as a list. We call the function rain.day once we determine that the current day is rainy. This function
allows to select from a set of distributions: exponential, Weibull, gamma, and skewed.

markov.rain.seq: This function implements a Markov chain for rainfall generation. Its arguments include the transition probability of the Markov chain and the number of days. This function is the first step of markov.rain.

markov.rain: This function implements two-step rainfall generation. Its arguments include the transition probability of the Markov chain and the number of days, which are used for the first step. In addition, one argument is the list of amount parameters, which is passed internally to the rain.day function described earlier.

The function pet calculates the Priestley-Taylor model as a function of radiation values under various temperature values. It also calculates the Penman model (radiation term, aerodynamic term, and total) as a function of radiation for various temperature values and for fixed RH and wind speed conditions.
}

\value{
rain.day: amount of rain in a day
markov.rain.seq:
\item{x}{amount of rain in a day}
\item{wet.days}{number of days with rain}
\item{expected.wet.days}{expected number of days with rain}
\item{dry.days}{number of days without rain}
\item{expected.dry.days}{expected number of days without rain}
markov.rain:
\item{x}{amount of rain in a day}
\item{wet.days}{number of days with rain}
\item{expected.wet.days}{expected number of days with rain}
\item{dry.days}{number of days without rain}
\item{expected.dry.days}{expected number of days without rain}
\item{rain.tot}{total amount of rain}
\item{expec.rain.tot}{expected total amount of rain}
\item{rain.avg}{average amount of rain}
\item{expected.rain.avg}{expected average amount of rain}
\item{rain.wet.avg}{average amount of rain for rainy days only}
\item{expected.wet.avg}{expected average amount of rain for rainy days only}
pet:
data.frame(Rad,PETR.taylor,PET.penman,PETR.penman,PETA.penman)

}

\references{
Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.

Neitsch, S. L., J. G. Arnold, J. R. Kiniry, J. R. Williams, and K. W. King. 2002. Soil and Water Assessment Tool. Theoretical Documentation Version 2000. Temple, TX: Grassland, Soil And Water Research Laboratory, Agricultural Research Service.

Richardson, C. W., and A. D. Nicks. 1990. Weather generator description. In EPIC-Erosion/Productivity Impact Calculator: 1 Model Documentation, eds. A. N. Sharpley, and J. R. Williams, 93-104. Durant, OK and Temple, TX: United States Department of Agriculture, Agricultural Research Center. Technical Bulletin No: 1768.

Williams, J. R., C. A. Jones, and P. T. Dyke. 1990. The EPIC model. In EPIC-Erosion/Productivity Impact Calculator: 1 Model Documentation, eds. A. N. Sharpley, and J. R. Williams, 3-92. Durant, OK and Temple, TX: United States Department of Agriculture, Agricultural Research Center, USDA Technical Bulletin No: 1768.


}
\author{
Miguel F. Acevedo \email{Acevedo@unt.edu}
}

\note{
Model functions are employed mainly to define the ODE to be simulated by sim, simd, simt and other simulation functions.

Nominal parameter values are defined in input files. Variation of param values are defined in lists.

Input files are in 'datafiles.zip' in directory 'datafiles' and organized by chapters of Acevedo (2012). Input files are required to run the examples below.    
}


\seealso{
Soil water functions e.g., \code{\link{infilt.rate}}, Solar functions e.g., \code{\link{sun.rad.hr}}
}
\examples{
\dontrun{
amount.param=list(mu=5,std=8,skew=2, shape=1.3, model.pdf ="w")
ndays=30
# rainy followed by rainy
P <- matrix(c(0.4,0.2,0.6,0.8), ncol=2, byrow=T)
rainy1 <- markov.rain(P, ndays, amount.param)
mtext(side=3,line=-1,paste("Rainy after rainy Pr=0.80","Prop rainy=",round(rainy1$wet.days/ndays,2)),cex=0.8)

pet.test <- pet(0, Rad =seq(10,30), Temp=c(10,20,30), RH=70, wind=2, plot=T)

}
}

\keyword{ math }
\keyword{ iteration }
