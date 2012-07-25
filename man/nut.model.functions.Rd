\name{Nutrient Model Functions}

\alias{nut.cycle}
\alias{nut.cycle.z}
\alias{nut.cycle.g}
\alias{nut.river}

\title{Functions to define ecosystem models for simulation.}

\description{
Model functions are used as part of argument model to call simulation functions.
}
\usage{
nut.cycle(t,p,x)
nut.cycle.z(t,p,x)
nut.cycle.g(t,p,x,tz)
nut.river(t,p,x)
}
\arguments{
  \item{t}{ time}
  \item{p}{ parameters, an array}
  \item{x}{ state variable}
  \item{tz}{ times of sudden disturbance}
}

\details{
Function nut.cycle has a function f with the ODE resulting from the multiplication of the matrix k of coefficients
by the X vector and the addition of a proportional u input/output (e.g., donor) component and
the independent input/output forcing terms U. State varibale is concentration.

Functions nut.cycle.z and nut.cycle.g define discontinuous disturbance regime. Function nut.cycle.z defines times of discontinuities, and nut.cycle.g applies a linear disturbance at tz times defined by nut.cycle.z.

Function nut.river: Nitrogen is controlled by load, which is calculated as discharge over volume, inflow from upstream conditions (when multiplied by flow/volume), and depuration or washout (which is the concentration multiplied by flow/volume).

}

\value{
Model functions: Rate of change or derivative of model.
}
\references{
Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.
}
\author{
Miguel F. Acevedo \email{Acevedo@unt.edu}
}
\note{
Model functions are employed mainly to define the ODE to be simulated by sim, simd, simt and other simulation functions.

Nominal parameter values are defined in input files. Variation of param values are defined in lists.

Input files are in 'datafiles.zip' in directory 'datafiles' and organized by chapters of Acevedo (2012). Input files are required to run the examples below.    
}

\seealso{Simulation functions \code{\link{sim}}, \code{\link{simd}} 
}
\examples{
\dontrun{

# Nitrogen cycle
cycle <- list(f=nut.cycle, z=nut.cycle.z, g=nut.cycle.g)
param <- list(plab="uh",pval=c(-0.5,-0.8))
t.X <- simd(cycle,"chp13/Ncycle-inp.csv",param)

# river no algae no flow and low flow,
river <- list(f=nut.river); param <- list(plab="Q",pval=c(0,0.060,0.600))
t.X <- sim(river,file="chp14/river-no-algae-inp.csv", param, pdfout=TRUE)

}
}
\keyword{ math }
\keyword{ iteration }
