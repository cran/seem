\name{Temp response functions}

\alias{parab}
\alias{kT.rate}
\alias{Q10}
\alias{Q10.mruns}

\title{Functions to define temperature response.}

\description{
Functions parab, kt.rate, Q10 are for single temperature response. Function Q10.mruns explores sensitivity to parameters of Q10. 
}
\usage{
parab(Temp,a,b)
kT.rate(Temp,k0,q10)
Q10(Temp, param)
Q10.mruns(Temp, param.nom, Topt.sens, Tmax.sens, q10.sens)
}
\arguments{
  \item{Temp}{temperature}
  \item{a}{parameter}
  \item{b}{parameter}
  \item{k0}{parameter}
  \item{q10}{parameter}
  \item{param}{parameter array c(Topt,Tmax,q10) where Topt is temp for which optimum occurs, Tmax is maximum temperature, and q10 is coefficient}
  \item{param.nom}{nominal values of parameters c(Topt,Tmax,q10)}
  \item{Topt.sens}{values to vary parameter Topt}
  \item{Tmax.sens}{values to vary parameter Tmax}
  \item{q10.sens}{values to vary parameter q10}
}

\details{
The three functions parab, kT.rate, and Q10, calculate the response to temperature using the parabolic, kT by q10, and the gamma-like q10.

Function Q10.mruns explores sensitivity to parameters of Q10.
}

\value{
Temperature response
}
\references{
Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.
}
\author{
Miguel F. Acevedo \email{Acevedo@unt.edu}
}
\note{
Nominal parameter values are defined in input files. Variation of param values are defined in lists.

Input files are in 'datafiles.zip' in directory 'datafiles' and organized by chapters of Acevedo (2012). Input files are required to run the examples below.    
}

\seealso{Simulation functions \code{\link{sim}}, \code{\link{simd}} 
}
\examples{
\dontrun{

# parabolic temp response
Temp <- seq(0,50); nT <- length(Temp)
Tu <- c(40,45,50);np <- length(Tu)
Tl <- 0
# dimension array
y <- matrix(nrow=nT,ncol=np)
# loop for the various optimal temperatures and then for temp range
for(j in 1:np){
 for(i in 1:nT) y[i,j] <- parab(Temp[i],Tl,Tu[j])
}

# rate kT
Temp <- seq(20,30); nT <- length(Temp)
k0 <- 1
q10 <- c(1,1.5,2); np <- length(q10)
kT <- matrix(nrow=nT,ncol=np)
for(i in 1:np){
kT[,i] <- kT.rate(Temp,k0,q10[i])/kT.rate(Temp[1],k0,q10[i])
}

# Q10
Temp <- seq(0,50); nT <- length(Temp)
Topt= 25; Tmax <- 40; q10 <- 1.5
fT <- Q10(Temp, param=c(Topt,Tmax,q10))

Temp <- seq(0,50); nT <- length(Temp)
# param Topt, Tmax,q10
param.nom <- c(25,40,1.5)
Topt.sens <- c(20,25,30)
Tmax.sens <- c(40,45,50)
q10.sens <- c(1,1.5,2)
np <- length(q10)
test <- Q10.mruns(Temp,param.nom,Topt.sens,Tmax.sens,q10.sens)

}
}
\keyword{ math }
\keyword{ iteration }
