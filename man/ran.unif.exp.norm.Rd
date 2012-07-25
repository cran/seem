\name{Random number gen}

\alias{ran.unif}
\alias{ran.exp}
\alias{ran.norm}

\title{
Custom random number generators for didactic purposes
}
\description{
Random generator of uniform RV between a and b, exponential, and standard normal RV
}
\usage{
ran.unif(a,b)
ran.exp(rate)
ran.norm()
}
\arguments{
  \item{a}{lower value}
  \item{b}{upper value}
  \item{rate}{rate parameter of exponential RV}
}
\details{
For didactic purposes, we show how to build custom random number generators from the basic uniform distribution between 0 and 1. The following function defines a uniform generator between a and b using the scaling and shifting of U(0,1).

}
\value{
One value of a RV according to distribution selected
}
\references{
Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.
}
\author{
Miguel F. Acevedo \email{Acevedo@unt.edu}
}
\note{
Didactic purposes. Use runif,rexp, rnorm for your real work.
}

\seealso{ Random number generators \code{\link{runif}}, \code{\link{rexp}}, \code{\link{rnorm}}
}
\examples{
yunif <- array(); a <- 5;b <- 15
for(i in 1:10) yunif[i] <- ran.unif(a,b)
hist(yunif,prob=TRUE); title(sub="Uniform in 5,15") 

yexp <- array(); rate <- 1/10
for(i in 1:10) yexp[i] <- ran.exp(rate)
hist(yexp,prob=TRUE); title(sub="Exponential of rate 0.1") 

ynorm <- array(); 
for(i in 1:10) ynorm[i] <- ran.norm()
hist(ynorm,prob=TRUE); title(sub="Std Normal") 
}

\keyword{ distribution }
\keyword{ datagen }
