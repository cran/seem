\name{Mat projection}
\alias{sim.proj}
\title{
Matrix projection
}
\description{
Simulates structured population by projection matrix
}
\usage{
sim.proj(t, x0, proj.mat, vars.plot, tunit, xunit, pdfout = F)
}
\arguments{
  \item{t}{projection times}
  \item{x0}{initial condition}
  \item{proj.mat}{projection matrix}
  \item{vars.plot}{variables to plot}
  \item{tunit}{units for time}
  \item{xunit}{units for population density}
  \item{pdfout}{switch to produce pdf files, (optional, default F)}
}

\details{
To perform the simulation, we write a loop for time to iterate the multiplication of the matrix proj.mat by the state vector x. As time progresses in the loop, it calculates the total population xtot, the proportions xp in each class with respect to the total, and an approximation lambda to the dominant eigenvalue as the ratio of the total population at adjacent times, and the population is scaled with respect to the last class. Note that the calculation of lambda is the ratio of current year to previous. 
}
\value{
\item{x}{results for all age classes}
\item{xp}{results for all age classes, as proportion of total}
\item{xs}{results for all age classes, as proportion of last class}
\item{lambda}{ratio of two consecutive values of xtot}
\item{xtot}{total population}
}
\references{
Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.
}
\author{
Miguel F. Acevedo \email{Acevedo@unt.edu}
}
\note{
Input files are in 'datafiles.zip' in directory 'datafiles' and organized by chapters of Acevedo (2012). Input files are required to run the examples below.
}

\seealso{
Simulation functions \code{\link{sim.comp}}, \code{\link{sim.mruns}}, \code{\link{sim.rnum}}, \code{\link{sim}}
}
\examples{
P <- matrix(c(0,2,3, 0.5,0,0, 0,0.5,0), byrow=TRUE, ncol=3)
x0 <- c(100,50,25); xunit="(Indiv)"
t=seq(0,10,1); tunit="(Year)"
x <- sim.proj(t, x0, P, vars.plot=c(1:3),tunit, xunit)
}

\keyword{ math }
\keyword{ iteration }
