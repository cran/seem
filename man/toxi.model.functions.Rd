\name{Ecotox Model Functions}

\alias{toxi1.bioacc}
\alias{toxi2.bioacc}
\alias{toxi.multi.bioacc}

\title{Functions to define ecotoxicology models for simulation}

\description{
Functions used as part of argument model to call simulation functions. 

}
\usage{
toxi1.bioacc(t,p,x)
toxi2.bioacc(t,p,x) 
toxi.multi.bioacc(t,p,x)
}
\arguments{
  \item{t}{time}
  \item{p}{parameters, an array}
  \item{x}{state variable}
}
\details{
State variable is concentration in a compartment.

Function toxi1.bioacc is for one compartment. It is essentially the same as the expon.forced function, but adapted to more parameters to calculate the uptake rate.

Function toxi2.bioacc is for 2 compartments. The diagonal elements of the inter-compartment transfer matrix, are initialized to zero and then completed assuming that they should be equal to the negative of the sum of the column elements. Uptake and depuration are calculated separately.

Function Toxi.multi.bioacc is . Variation of the tox2.bioacc function but includes the pulse exposure, one more compartment, and M3 breakdown kinetics.
}

\value{
Rate of change or derivative of model.
}
\references{
Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.

Mancini, J.L. 1983. A method for calculating effects, on aquatic organisms, of time varying concentrations. Water Research 17:1355-1362.

Thomann, R.V. 1989. Deterministic and statistical models of chemical fate in aquatic systems. In Ecotoxicology: Problems and Approaches, eds. S.A. Levin, M.A. Harwell, J.R. Kelly, and K.D. Kimball, 245-277. New York: Springer-Verlag.

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
tox1 <- list(f=toxi1.bioacc)
param <- list(plab="Cw",pval=c(0,0.2,0.5,1))
t.X <- sim(tox1,"chp11/simple-bioaccu-inp.csv",param,pdfout=TRUE)

tox2 <- list(f=toxi2.bioacc)
param <- list(plab="k12",pval=c(0.1,0.2))
t.X <- sim(tox2,"chp11/two-compart-closed-inp.csv",param,pdfout=TRUE)

toxm <- list(f=toxi.multi.bioacc)
param <- list(plab="Cw",pval=c(1,5,10))
t.X <- sim(toxm,"chp11/multi-compart-inp.csv",param,pdfout=TRUE)
}
}
\keyword{ math }
\keyword{ iteration }
