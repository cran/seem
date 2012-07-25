\name{Soil water}

\alias{soil.water}
\alias{infilt.rate}
\alias{soilwat.1}
\alias{soilwat.n}
\alias{green.ampt.ramos}
\alias{green.ampt.maila}
\alias{sim.ga}


\title{Soil water dynamics}
\description{Soil hydraulic properties and model functions to simulate soil water dynamics}

\usage{
soil.water(theta, param, title = "", plot = F, pdfout = F)
infilt.rate(param, rain, v.soil)
soilwat.1(t, p, x)
soilwat.n(t, p, x)
green.ampt.ramos(t, p, x)
green.ampt.maila(t, p, x)
sim.ga(t, p)
}

\arguments{
  \item{theta}{seq of soil water content as wetness (vol fraction)}
  \item{param}{For soil water: Pb bubbling suction in cm, b pore size distribution index (adim), theta.r residual water, theta.s total porosity. For Infiltration rate: conductivities Kd, Ks, saturation capacity Z, field capacity Fc.coeff and porosity}
  \item{title}{string for title}
  \item{plot}{logical variable to decide whether to plot}
  \item{pdfout}{logical avriable to decide whther to generate PDF file}
  \item{rain}{For infilt.rate: water input at the surface}
  \item{v.soil}{for infilt.rate: soil water content}

For model functions
  \item{t}{time}
  \item{p}{parameters, an array}
  \item{x}{state variable}
}

\details{


soil.water: soil hydraulic properties based on Brooks Corey. Effective saturation is calculated as a function of theta. Then, matric suction P, capacity C, and conductivity K are calculated as a function of effective saturation. Calculations are plotted in both regular and semilogarithmic axes (the y-axis).

infilt.rate: infiltration rate function. Infiltration rate capacity based on deficit and limited by dry and sat infiltration, infiltration capacity Dinkin & Nazimov (1995).

Model functions are used as part of argument model to call simulation functions.

soilwat.1: This function implements the combination of infiltration and percolation as in Equation 15.77 of Acevedo (2012) and as a function or input intensity and duration.

soilwat.n: This model function is an extension of soilwat.1 adding the ability to calculate percolation for multiple layers. Similarly to soilwat.1, this function is called from RK4 for simulation using sim.

green.ampt.ramos: is utilized as a model function to define the Green-Ampt model and to respond to a call from sim using ramos that requires f and df/dx.

green.ampt.maila: This function is an alternative to calculate the rate dF for the numerical integration of F in the Green-Ampt model using the Mailapalli et al. (2009) method given in Equation 15.99 and based on the Ramos algorithm. The function is called by the simulation function sim.ga.

sim.ga: this function implements a numerical simulation using green.ampt.maila, described above. It uses a simple loop for time. Arguments are time sequence and a parameter set p. This function is an alternative to utilize sim with a call to the ramos function instead of RK4.

}

\value{
soil water:
\item{Ks}{saturated hydraulic conductivity}
\item{Pf}{Potential wet front}
\item{var.theta}{data.frame(theta,theta.e,P,C,K,LP,LC,LK), where theta.e is theta effective, P is matric potential, C is capacity, K is conductivity. LP,LC, ad LK are the log of P,C, and K}

Infiltration rate:
\item{infilt.cap}{infiltration rate capacity f}
\item{infilt.rate}{actual infiltration rate q} 
\item{percol}{percolation rate}

Model functions: Rate of change or derivative of model.

Green.ampt.ramos also returns the rate of change of flow

sim.ga returns a data.frame(t,x,f,q,Q) where t=time,x=depth, f=infiltration capacity, q=infiltration rate, Q=infiltrated depth.

}

\references{
Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.

Diskin, M. H., and N. Nazimov. 1995. Linear reservoir with feedback regulated inlet as a model for the infiltration process. Journal of Hydrology 172:313-330.

Mailapalli, D. R., W. W. Wallender, R. Singh, and N. S. Raghuwanshi. 2009. Application of a nonstandard explicit integration to solve Green and Ampt infiltration equation. Journal of Hydrologic Engineering 14:203-206.

van Genutchen, M. T. 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils. Soil Science Social of America Journal 44:892-898.

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
\code{\link{sim}}, \code{\link{simd}}, \code{\link{simr}}
}
\examples{
# sandy clay
param=c(Pb=29.17,b=10.4, theta.r=0.109,theta.s=0.43)
theta <- seq(0.2,0.43,0.01)
sw <- soil.water(theta, param, title="Sandy Clay")
\dontrun{
sw <- list(f=soilwat.1);param <- list(plab="Ic",pval=c(10,20,30,40))
t.X <- sim(sw,"chp15/sw-inp.csv", param, pdfout=TRUE)

sw2 <- list(f=soilwat.n);param <- list(plab="Ic",pval=c(10,20,30,40))
t.X <- sim(sw2,"chp15/sw-2layer-inp.csv", param, pdfout=TRUE)

# sandy loam
param=c(Pb=14.66,b=4.9, theta.r=0.041,theta.s=0.453)
init = 0.1 # volumetric water content m3/m3
sw <- soil.water(theta=init,param)
Ks<- sw$Ks*3600*10 # cm/s to mm/h
Pf <- sw$Pf*10 # cm to mm
porosity<- param[4]
dt=0.00001
# simulate using simr
ga <- list(f=green.ampt.ramos);param <- list(plab="Ks",pva
l=c(10,20,30,40))
t.X <- simr(ga,"chp15/ga-inp.csv", param, pdfout=TRUE)
# simulate using sim.ga
rain <- 200
p <- c(Pf,Ks,porosity,init,dt,rain)
t <- seq(0,0.2,dt)
sga <- sim.ga(t,p)

}}

\keyword{ math }
\keyword{ iteration }
