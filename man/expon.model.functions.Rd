\name{Exponential Model Functions}

\alias{expon}
\alias{expon.rand}
\alias{expon.forced}
\alias{expon.var}
\alias{expon.z}
\alias{expon.g}

\title{Functions to define exponential models for simulation}

\description{
Functions used as part of argument model to call simulation functions. 
}
\usage{
expon(t,p,x)
expon.rand(t,p,x)
expon.forced(t,p,x)
expon.var(t,p,x)
expon.z(t,p,x)
expon.g(t,p,x,tz)
}
\arguments{
  \item{t}{ time}
  \item{p}{ parameters, an array}
  \item{x}{ state variable}
  \item{tz}{ times of sudden disturbance}
}
\details{
State variable can be population density, or chemical concentration. For expon, p is dim 1: the rate coeff. For expon.rand, p is dim 2: mean and sdev of rate coeff. For expon.forced, p is dim 3: rate coeff and disturbance parameters. For expon.var, p is dim 4: mean, drift coeff, amplitude, and period of sinusoidal variation.

Functions expon.z and expon.g define discontinuous disturbance regime. Function expon.z defines times of discontinuities, and expon.g applies a linear disturbance at tz times defined by expon.z.   
}
\value{
Rate of change or derivative of model. Except expon.z and expon.g that return disturbance regime. 
}
\references{
Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.
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
# exponential
model <- list(f=expon)
# single run
file<-"chp4/exp-pop-inp.csv"
out <- sim.comp(model,file)
# multiple runs
param <- list(plab="r", pval = seq(-0.02,0.02,0.01))
out.r <- sim.mruns(model, file, param)
# Monte Carlo run
model <- list(f=expon.rand)
file<-"chp4/exp-rnum-inp.csv"
out.rk <- sim.rnum(model,file)
# Harvest, multiple runs 
exp.f <- list(f=expon.forced)
param <- list(plab="Ha",pval=c(0,-0.1,-0.2))
t.X <- sim(exp.f,"chp7/exp-harvest-inp.csv",param)
# Seasonality 
exp.v <- list(f=expon.var)
param <- list(plab="rd",pval=c(0,-0.005,-0.01))
t.X <- sim(exp.v,"chp7/exp-var-inp.csv",param)
# sudden disturbance
exp.sud <- list(f=expon,z=expon.z,g=expon.g)
param <- list(plab="Ha",pval=c(-10,-20,-30))
t.X <- simd(exp.sud,file="chp7/exp-sud-inp.csv",param)
}
}
\keyword{ math }
\keyword{ iteration }
