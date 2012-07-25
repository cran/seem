\name{Stability}
\alias{sta.plot}
\alias{phase.plot}
\alias{root.plot}
\alias{root.locus.quad}
\alias{ReIm.plot}
\alias{D.plot}

\title{
Stability analysis and plots (didactic purpose)
}
\description{
For a 2D linear dynamical system, function sta.plot is used to plot state trajectories vs. time. Function sta.plot is used to plot phase-plane and illustrate orbit by arrow. Function root.plot just plots the roots of A. Function root.locus.quad plots is used to plot the locus of the root by sweeping a parameter. It calculates eigenvalues in four quadrants, to be used by functions ReIm.plot and D.plot. 
}
\usage{
sta.plot(A, x0, t)
phase.plot(A, x0, t, long = 20)
root.plot(A)
root.locus.quad(det, Tr)
ReIm.plot(x, j, k)
D.plot(x, i)
}
\arguments{
  \item{A}{matrix}
  \item{x0}{initial condition}
  \item{t}{time sequence}
  \item{long}{determine size of arrows}
  \item{det}{determinant of matrix}
  \item{Tr}{trace of matrix}
  \item{x}{list containing sequences of det and Tr and their discriminant, and eigenvalues (see output from root.locus.quad)}
  \item{j}{quadrant number}
  \item{i}{quadrant number}
  \item{k}{position in array}
}
\details{
Given matrix A of 2D linear dynamical system, function sta.plot is used to plot X trajectories vs. time. Function sta.plot is used to plot phase-plane and illustrate orbit by arrow. Determinant, trace, discriminant, and eigenvalues are printed on the plot.

Function root.plot just plots the roots of A. Function root.locus.quad calculates eigenvalues in four quadrants, to be used by plotting functions ReIm.plot and D.plot.
}

\value{
Function root.locus.quad:
\item{det}{sequence of det values}
\item{Tr}{sequence of trace values}
\item{D}{ discriminant for each pair of det and Tr}
\item{ReL}{real part of eigenvalues for each pair of det and Tr} 
\item{ImL}{imaginary part of eigenvalues for each pair of det and Tr} 
Other functions do not return values but produce plots.
}
\references{
Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.
}
\author{
Miguel F. Acevedo \email{Acevedo@unt.edu}
}
\note{
For didactic purposes.
}

\seealso{
Matrix operations \code{\link{eigen}}, \code{\link{solve}}

}
\examples{
\dontrun{
t <- seq(0,100,0.1)
# pos det, neg trace, pos D
A <- matrix(c(-0.1,-0.1,-0.2,-0.3),2,byrow=T);x0 <-c(10,5)
sta.plot(A,x0,t)

# pos det, neg trace, pos D
A <- matrix(c(-0.1,-0.1,-0.2,-0.3),2,byrow=T);x0 <-c(10,5)
t <- seq(0,10,0.1)
phase.plot(A,x0,t,long=12)

x <- list()
#pos det, neg trace
det <- seq(0,1,0.5);Tr <- seq(-3,0,0.001) 
x[[1]] <- root.locus.quad(det,Tr)
# pos det, pos trace
det <- seq(0,1,0.5);Tr <- seq(0,3,0.001) 
x[[3]] <- root.locus.quad(det,Tr)
#neg det, neg trace
det <- seq(-1,0,0.5);Tr <- seq(-3,0,0.001) 
x[[2]] <- root.locus.quad(det,Tr)
#neg det, pos trace
det <- seq(-1,0,0.5);Tr <- seq(0,3,0.001) 
x[[4]] <- root.locus.quad(det,Tr)

 mat <- matrix(1:4,2,2,byrow=T)
 nf <- layout(mat, widths=rep(7/2,2), heights=rep(7/2,2), TRUE)
 par(mar=c(4,4,1,.5),xaxs="i",yaxs="i")
for(i in 1:4) D.plot(x[[i]],i)

 mat <- matrix(1:4,2,2,byrow=T)
 nf <- layout(mat, widths=rep(7/2,2), heights=rep(7/2,2), TRUE)
 par(mar=c(4,4,1,.5),xaxs="r",yaxs="r")

  j=1
  ReIm.plot(x[[j]],j,1)
  ReIm.plot(x[[j]],j,2)
  j=3
  ReIm.plot(x[[j]],j,1)
  ReIm.plot(x[[j]],j,2)
  j=2
  ReIm.plot(x[[j]],j,1)
  ReIm.plot(x[[j]],j,2)
  j=4
  ReIm.plot(x[[j]],j,1)
  ReIm.plot(x[[j]],j,2)
}
}

\keyword{ math }
\keyword{ dplot }

