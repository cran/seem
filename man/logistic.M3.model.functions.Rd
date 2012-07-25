\name{Logistic and M3 Model Functions}

\alias{logistic}
\alias{logistic.forced}
\alias{logistic.var}
\alias{logistic.z}
\alias{logistic.g}
\alias{M3}
\alias{monod.batch}

\title{Functions to define logistic and M3 (Michaelis-Menten-Monod) models for simulation}

\description{
Functions used as part of argument model to call simulation functions. 
}
\usage{
logistic(t,p,x)
logistic.forced(t,p,x)
logistic.var(t,p,x)
logistic.z(t,p,x)
logistic.g(t,p,x,tz)
M3(t,p,x)
monod.batch(t,p,x)

}
\arguments{
  \item{t}{ time}
  \item{p}{ parameters, an array}
  \item{x}{ state variable}
  \item{tz}{ times of sudden disturbance}
}
\details{
State variable can be population density (logistic family), or chemical concentration (M3). Two states for monod.batch: chemical concentration and population density. 

For logistic, p is dim 2: the intrinsic growth coeff, and the carrying capacity. For logistic.forced, p is dim 4:  the intrinsic growth coeff, the carrying capacity, and disturbance parameters. For logistic.var, p is dim 5: mean of the intrinsic growth coeff, the carrying capacity, and, drift coeff, amplitude, and period of sinusoidal variation.

Functions logistic.z and logistic.g define discontinuous disturbance regime. Function logistic.z defines times of discontinuities, and logistic.g applies a linear disturbance at tz times defined by logistic.z.   

For function M3, p is dim 2: half-rate concentration and maximum rate. For monod-batch p is dim 4: half-rate concentration, maximum rate, yield coefficient, and death rate.  
}


\value{
Rate of change or derivative of model. Except logistic.z and logistic.g that return disturbance regime. 
}
\references{
Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.

Hallam, T.G., 1986b. Population dynamics in a homogeneous environment. In Mathetmatical Ecology, eds. T. G. Hallam, and S. A. Levin, 61-94. New York: Springer-Verlag.

Hanson, F.B., and H.C. Tuckwell. 1981. Logistic growth with random density independent disasters. Theoretical Population Biology 19:1-18.

}
\author{
Miguel F. Acevedo \email{Acevedo@unt.edu}
}
\note{
Model functions are employed mainly to define the ODE to be simulated by sim.comp, sim.mruns, sim.rnum, sim, simd and other simulation functions.

Nominal parameter values are defined in input files. Variation of param values are defined in lists.

Input files are in 'datafiles.zip' in directory 'datafiles' and organized by chapters of Acevedo (2012). Input files are required to run the examples below.    
}

\seealso{Simulation functions \code{\link{sim.comp}}, \code{\link{sim.rnum}}, \code{\link{sim.mruns}}, \code{\link{sim}}, \code{\link{simd}} Methods \code{\link{euler}}, \code{\link{RK4}}
}
\examples{
\dontrun{
# logistic
logis <- list(f=logistic)
# single run
t.X <- sim(logis,file="chp6/logistic-inp.csv")
# multiple runs
param <- list(plab="r", pval = seq(0.2,0.6,0.2))
t.X <- sim(logis,"chp6/logistic-inp.csv",param)
# Harvest, multiple runs 
logis.f <- list(f=logistic.forced)
param <- list(plab="Ha",pval=c(0,-0.1,-0.2,-0.3))
t.X <- sim(logis.f,"chp7/logis-harvest-inp.csv",param)
# Seasonality 
logis.v <- list(f=logistic.var)
param <- list(plab="Kd",pval=c(0,-0.5,-0.6))
t.X <- sim(logis.v,"chp7/logis-var-inp.csv",param)
# sudden disturbance
logis.sud <- list(f=logistic,z=logistic.z,g=logistic.g)
t.X <- simd(logis.sud,file="chp7/logis-sud-inp.csv")
# M3 decay multiple runs
m3decay<- list(f=M3)
param <- list(plab="Kmax", pval = seq(-20,-60,-10))
t.X <- sim(m3decay,"chp6/m3decay-inp.csv", param)
# monod batch multiple runs
monod <- list(f=monod.batch)
param <- list(plab="Kmax", pval = seq(0.5,1.5,0.5))
t.X <- sim(monod,"chp6/monod-batch-inp.csv", param)

}
}
\keyword{ math }
\keyword{ iteration }
