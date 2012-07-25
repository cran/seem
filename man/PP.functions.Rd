\name{Primary productivity functions}

\alias{PP.Smith}
\alias{PP.Thornley}
\alias{PP.Steele}
\alias{PP.Eilers}
\alias{PPT.Smith}
\alias{PPT.Steele}
\alias{PPrates.depth}
\alias{DO.PP.pond}

\title{Functions to determine primary productivity}

\description{
Functions to calculate PP rates in aquatic ecosystems. DO.PP.pond is a model functions are used as part of argument model to call simulation functions.
}

\usage{
PP.Smith(L, Pmax, alpha)
PP.Thornley(L, Pmax, alpha, ksi)
PP.Steele(L, Pmax, Lopt)
PP.Eilers(L, Pmax, alpha, Lopt)
PPT.Smith(Ls,k,z,Pmax,alpha)
PPT.Steele(Ls,k,z,Pmax,Lopt)
PPrates.depth(z, param,sw.plot=T)
DO.PP.pond(t,p,x)
}
\arguments{
  \item{L}{sunlight intensity}
  \item{Pmax}{maximum rate}
  \item{alpha}{slope at L=0}
  \item{ksi}{shape parameter}
  \item{Lopt}{sunlight at which PP is maximum}
  \item{Ls}{sunlight at the surface (actually just below surface)}
  \item{k}{sunlight attenuation coefficient}
  \item{z}{depth}
  \item{param}{parameter array c(k,Ls, Pmax,alpha,Lopt)}
  \item{sw.plot}{logical variable to decide whether to plot}
For DO.PP.pond
  \item{t}{time}
  \item{p}{parameters, an array}
  \item{x}{state variable}
}

\details{
The four functions PP.Smith, PP.Thornley, PP.Steele, and PP.Eilers calculate the response to sunlight by the Smith, Thornley, Steele, and Eilers-Peeters models.

PPT.Smith and PPT.Steele functions implement depth-averaged Smith and Steele models. In addition to the model parameters, they use sunlight at the subsurface, attenuation coefficient, and depth as arguments.

PPrates.depth calculates productivity as a function of depth using PPT.Smith and PPT.Steele.
}

\value{
PP and PPT functions return primary productivity rate.
PPrates.depth returns a data.frame with z,L,PP,PPT
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

L <- seq(0,300,10);nL <- length(L)

# Smith model
Pmax=1; alpha <- c(0.01,0.02,0.03)
np <- length(alpha)
PP <- matrix(nrow=nL,ncol=np)
for(i in 1:np)
PP[,i] <- PP.Smith(L,Pmax,alpha[i])

# Thornley NRH model
L <- seq(0,300,1);nL <- length(L)
Pmax=1; alpha <-c(0.01,0.02,0.03); ksi <- 0.1
np <- length(alpha)
PP <- matrix(nrow=nL,ncol=np)
for(i in 1:np)
 PP[,i] <- PP.Thornley(L,Pmax,alpha[i],ksi)

# Steele model
Pmax=1; Lopt <- c(150,200,250)
np <- length(Lopt)
PP <- matrix(nrow=nL,ncol=np)
for(i in 1:np)
PP[,i] <- PP.Steele(L,Pmax,Lopt[i])

#  Eiler-Peeters
Pmax=1; Lopt <- c(150,200,250); alpha=0.02
np <- length(Lopt)
PP <- matrix(nrow=nL,ncol=np)
for(i in 1:np)
PP[,i] <- PP.Eilers(L,Pmax,alpha,Lopt[i])

z <- seq(0,3,0.1)
k = 0.5; Ls=300
Pmax =1; alpha =0.02
Pmax =1; Lopt = 200
PPrates.depth(z, param=c(k,Ls,Pmax,alpha,Lopt))

# cool water
pond <- list(f=DO.PP.pond);param <- list(plab="Rsp",pval=c(0.18,0.2,0.22))
t.X <- sim(pond,"chp14/DO-PP-inp.csv", param)
}
}

\keyword{ math }
\keyword{ iteration }
