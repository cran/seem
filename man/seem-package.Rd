\name{seem-package}
\alias{seem-package}
\alias{seem}

\docType{package}

\title{Simulation of Ecological and Environmental Models}

\description{Used in Acevedo, M.F. 2012. "Simulation of Ecological and Environmental Models", CRC press.}

\details{

\tabular{ll}{
Package: \tab seem\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-06-05\cr
License: \tab GPL-2\cr
}

Package to perform simulations, including sensitivity analysis. It provides a suite of simulator functions: sim.comp, sim.mruns, sim.rnum, sim, simd, simr. Model functions are employed to define the ODE to be simulated by the simulation functions. Nominal parameter values are defined in input files. Variation of param values are defined in lists.

Input files are in 'datafiles.zip' in directory 'datafiles' and organized by chapters of Acevedo (2012). Input files are required to run the examples below.

In addition, seem includes wrappers to use several simulation programs written in Fortran and C. These include cerio.F, semi.F, zelig.F, river.C, and forsucc.C. The zelig simulator used here is based on the version 2.3 by D.L. Urban.   
}

\author{
M.F. Acevedo <Acevedo@unt.edu>
}
\references{
Acevedo, M.F. 2012. "Simulation of Ecological and Environmental Models", CRC press. 
}

\keyword{package}

\seealso{
Model functions \code{\link{expon}}, \code{\link{expon.z}}, \code{\link{expon.g}}, \code{\link{logistic}}, \code{\link{two.stage.cont}}, \code{\link{two.stage.cont.delay}}, \code{\link{two.stage.x0}}, and many others; Methods  \code{\link{RK4}}, \code{\link{euler}}, \code{\link{RK4D}}, \code{\link{ramos}}, Simulator output functions \code{\link{onerun.out}}, \code{\link{mruns.out}}, \code{\link{rnum.out}}, \code{\link{vars.out}}
}

\examples{
\dontrun{
model <- list(f=expon);file<-"chp4/exp-pop-inp.csv"
out.eu <- sim.comp(model,file,method="euler",lab.out="Euler")

model <- list(f=expon);file<-"chp4/exp-pop-inp.csv"
param <- list(plab="r", pval = seq(-0.02,0.02,0.01))
out.r <- sim.mruns(model,file, param, pdfout=TRUE)

model<-list(f=expon.rand); file<-"chp4/exp-rnum-inp.csv"
out.rk <- sim.rnum(model,file)

logis <- list(f=logistic)
t.X <- sim(logis,"chp6/logistic-inp.csv")
param <- list(plab="r", pval = seq(0.2,0.6,0.2))
t.X <- sim(logis,"chp6/logistic-inp.csv",param)

exp.sud <- list(f=expon,z=expon.z,g=expon.g)
t.X <- simd(exp.sud,file="chp7/exp-sud-inp.csv")

logis.sud <- list(f=logistic,z=logistic.z,g=logistic.g)
param <- list(plab="r",pval=c(0.1,0.2))
out.f <- simd(logis.sud,file="chp7/logis-sud-inp.csv",param)

two.stage.delay <-list(f=two.stage.cont.delay,x0=two.stage.x0)
param <- list(plab="nd", pval = c(2,5,6,7))
y.s <- simt(two.stage.delay,"chp10/two-stage-delay-inp.csv", param)

ga <- list(f=green.ampt.ramos);param <- list(plab="Ks",pval=c(10,20,30,40))
t.X <- simr(ga,"chp15/ga-inp.csv", param, pdfout=TRUE)

spp <- c("Role1", "Role2", "Role3", "Role4")
fileprefix <- "rolesm"; label <- "Semi-Markov Order 1" 
fileout <- semi.F("chp16",fileprefix)
x <- read.plot.semi.out(fileout, spp, label)

}
}
