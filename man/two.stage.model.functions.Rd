\name{Structured Population Model Functions}

\alias{two.stage.cont}
\alias{two.stage.cont.delay}
\alias{two.stage.x0}


\title{Functions to define structured population models for simulation}

\description{
Functions used as part of argument model to call simulation functions. 
}
\usage{
two.stage.cont(t,p,x)
two.stage.cont.delay(t,p,x) 
two.stage.x0(p,X0)
}
\arguments{
  \item{t}{ time}
  \item{p}{ parameters, an array}
  \item{x}{ state variable}
  \item{X0}{ initial condition for class with delay}
}
\details{
State variable is population density.  
}


\value{
Rate of change or derivative of model. Except two.stage.x0 that returns distributed initial condition by allocation.
}
\references{
Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.

Lewis, E.R., 1977. Network Models in Population Biology. New York: Springer Verlag.

McDonald, N. 1978. Time Lags in Biological Models. New York: Springer Verlag.
}
\author{
Miguel F. Acevedo \email{Acevedo@unt.edu}
}
\note{
Model functions are employed mainly to define the ODE to be simulated by sim, simd, simt and other simulation functions.

Nominal parameter values are defined in input files. Variation of param values are defined in lists.

Input files are in 'datafiles.zip' in directory 'datafiles' and organized by chapters of Acevedo (2012). Input files are required to run the examples below.    
}

\seealso{Simulation functions \code{\link{sim.comp}}, \code{\link{sim.rnum}}, \code{\link{sim.mruns}}, \code{\link{sim}}, \code{\link{simd}} Methods \code{\link{euler}}, \code{\link{RK4}}
}
\examples{
\dontrun{
two.stage <-list(f=two.stage.cont)
param <- list(plab="b", pval = c(0.15,0.25,0.3,0.31))
t.X <- simt(two.stage,"chp10/two-stage-inp.csv", param,pdfout=TRUE)

two.stage.delay <-list(f=two.stage.cont.delay,x0=two.stage.x0)
param <- list(plab="nd", pval = c(2,5,6,7))
t.X <- simt(two.stage.delay,"chp10/two-stage-delay-inp.csv",
param,pdfout=TRUE)

}
}
\keyword{ math }
\keyword{ iteration }
