\name{Sim functions}
\alias{sim.comp}
\alias{sim.mruns}
\alias{sim.rnum}
\alias{sim}
\alias{simd}
\alias{simt}
\alias{simr}

\title{
Simulation functions  
}
\description{

Simulates a single run (sim.comp), multiple runs (sim.mruns), Monte Carlo run (sim.rnum), and either single or multiple runs (sim). Function sim.comp can use the euler method (for didactic purposes). Function sim.rnum uses function rnum for integration. Function sim is a more generic function.

Function simd is a modification of sim that allows for discontinuities defined in model functions *.z and *.g that are required as components of the model list argument.

Function simt is a modification of sim that allows for time delays defined in model functions that are required as components of the model list argument.

Function simr uses the ramos function.
}
\usage{
sim.comp(model, file, method = "RK4", pdfout = F, lab.out = "")
sim.mruns(model, file, param, pdfout = F, lab.out = "")
sim.rnum(model, file, pdfout = F, lab.out = "")
sim(model, file, param=NULL, pdfout=F, lab.out="")
simd(model, file, param=NULL, pdfout=F, lab.out="")
simt(model, file, param=NULL, pdfout=F, lab.out="")
simr(model, file, param=NULL, pdfout=F, lab.out="")

}
\arguments{
  \item{model}{model to simulate as a list (required)}
  \item{file}{path/name of input file (required)}
  \item{method}{name of function for num method, (optional, default RK4)}
  \item{pdfout}{switch to produce pdf files, (optional, default F)}
  \item{lab.out}{label to identify output files, (optional, default none}
  \item{param}{param set for multiple runs; label of the parameter to be varied and its value}

}
\details{
Functions with at least two arguments (except sim.mrus which requires 3): a model defined as a list and an input filename. Optionally, we specify a numerical method, a param list for multiple runs, a logical variable to generate pdf output, and a label to identify output files. These functions have two parts: 1) reading and organizing input, and 2) performing iteration and writing output.

First part: Here the filename variable file is an argument and has extension .cvs for comma-separated values; for example look at file chp4/exp-pop-inp.csv. Here this function assigns the information read from file to a data frame object named input. We read names "Label", "Value", etc from the header (first row), and then read the records in all rows. Thus, we form vectors with the values of the column input$Val and input$Lab. Finally the units for time t and state X are extracted from the input data frame using $Uni and stored as a vector in order to facilitate referring to the units. Second part. First, we call integration function passing as arguments the initial condition, the time sequence, the function f to integrate as component model$f of list model, the parameter and the time step. Then we format the output according to number of significant figures, organize the results in a data frame, and call output functions *.out to do the heavy work of plotting graphs and writing the files.

Function sim.mruns is an extension of sim.comp that includes an outer loop around the call to the integration. This outer loop varies according to the condition of each run. It includes an argument to pass the information about th. Variation in a parameter; it contains the label of the parameter to be varied and its values. The last part of the function prepares the output, to identify the units and the run number. It uses another output function because there are several columns for the output, one for each run. This function is mruns.out.

The function sim.rnum is similar to sim.comp. The major differences are that  function rnum is called to perform the integration and function rnum.out is called to write the output files.

The first part of the sim function is like the sim.mruns function, but it includes the
rates and generalize using positions for the parameters and initial conditions, since we can have
more than one. This function is the workhorse for most cases.

In sim, we collect the two parameters in an array, p, which will be an argument to the model
function f. A key to the functionality of sim is detecting the position iNS of the initial condition in the input file. For this purpose, it assumes that the initial condition is labeled X0 for one-dimensional systems or starts with X1.0 for multidimensional systems. In addition, it is assumed that the parameter starts listing in position 5 of the file. Now we know that the parameters are listed from 5 to iNS-1 and the initial conditions are from iNS to the next to last position. This information is used in the parametric loop when making a decision on whether to vary the initial condition or a parameter and which parameter to vary.

For the sake of generality, the value of NS is tested to decide whether to include the rate dX.dt for one-dimensional systems using the model function or ignore it for the multidimensional systems.
Also for the sake of generality, we have made the param argument optional with default value NULL. Therefore, the null condition of param is tested with is.null(param) to decide on making one run or many according to the values of param.

Function simd is a modification of sim that allows for discontinuities. The only changes with respect to sim are in the parameter loop, where just before the call to integration, we calculate tz, invoking function z; we call a modified Runge-Kutta function RK4D, which has new arguments g and tz; and the rate calculation for one-dimensional systems uses function g.

Function simt is an extension of function sim to include intermediate stages of a delay and the total for all classes. These changes only affect the lines within the loop for parameter p in the following manner. When there is a delay, (1) call model$x0 component to compose extended initial conditions, (2) declare a 3-D array, Y, to store the intermediate stages, and (3) sum over the intermediate stages to get the total for class 1. These three modifications are executed only when model$x0 is not NULL.

Function simr is the same as sim but uses the ramos function instead of RK4. It is applied mainly to simulate Green-Ampt soil water dynamics. The ramos function implements the nonstandard explicit integration algorithm (EIA) that requires df/dx in addition to f.
  
}
\value{
\item{input}{Echo of simulation input file, just to make sure we associate the execution results to the input actually employed.}
\item{param}{For multiple runs only: Echo of parameter label and values varied for multiple runs, just to make sure we associate the execution results to the input actually employed. }
\item{output}{Contains the execution results in columns with labels that include the units for each variable as extracted from the input file.}
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
Model functions \code{\link{expon}}, \code{\link{expon.z}}, \code{\link{expon.g}}, \code{\link{logistic}}, \code{\link{two.stage.cont}}, \code{\link{two.stage.cont.delay}},\code{\link{two.stage.x0}}, and many others; Methods  \code{\link{RK4}}, \code{\link{euler}}, \code{\link{RK4D}}, \code{\link{ramos}}, Simulator output functions \code{\link{onerun.out}}, \code{\link{mruns.out}}, \code{\link{rnum.out}}, \code{\link{vars.out}}
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

}

}
\keyword{ math }
\keyword{ iteration }
\keyword{ programming }
