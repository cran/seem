\name{Trees functions}

\alias{DH.allom}
\alias{D.H.est}
\alias{Dl.allom}
\alias{dD.rate}
\alias{dD.dt.analysis}
\alias{dD.dt.calibra}
\alias{parab1}
\alias{env.factors.analysis}
\alias{ge.growth}
\alias{expon.alf}
\alias{gauss.smf}

\title{Functions to define tree responses}

\description{Support and didactic functions for forest simulators} 

\usage{
DH.allom(D,param)
D.H.est(D.H,D.Hmax,species.name,pdfout=F,linregout=F)
Dl.allom(D,param)
dD.rate(D,G,Dmax,param.H,param.l)
dD.dt.analysis(D,G,Dmax,param.H,param.l,D0,t,pdfout=F)
dD.dt.calibra(D,dD,Dmax,param.H,param.l)
parab1(E,Eopt,p)
env.factors.analysis(nsp,splab,FE,Eopt,param)
ge.growth(nsp,splab,FE,Eopt,param.E,Ef,param.Ef,Dmax,Gmax,param.H,param.l,D0,nt)
expon.alf(E,p)
gauss.smf(x,p)
}

\arguments{
 \item{D}{diameter or set of diameter values}
 \item{param}{For DH.allom: Hmax, b3, and b4. For Dl.allom: c, hc, and d. Also list of parameters for response functions}
 \item{D.H}{data frame with D and H values}
 \item{D.Hmax}{array with Dmax and Hmax}
 \item{species.name}{char string for label}
 \item{pdfout}{logical to decide to generate PDF}
 \item{linregout}{logical to decide to output the linear regression results}
 \item{G}{set of values of max growth rate}
 \item{Dmax}{Max diameter}
 \item{param.H}{parameters for function DH.allom}
 \item{param.l}{parameters for function Dl.allom}
 \item{D0}{initial diameter}
 \item{t}{time sequence}
 \item{dD}{set of diameter increment values}
 \item{E}{Environmental factors}
 \item{Eopt}{Optimum values of E}
 \item{p}{parameter of response function}
 \item{nsp}{number of species}
 \item{splab}{labels of species}
 \item{FE}{list of response functions}
 \item{param.E}{parameters of response functions by species}
 \item{Ef}{functions Ef to generate random samples of environmental variables}
 \item{param.Ef}{parameter to generate random samples of environmental variables}
 \item{Gmax}{max or optimum growth rate coeff}
 \item{nt}{length of time sequence}
 \item{x}{soil moisture indices}
}

\details{
DH.allom: calculates tree height, given the diameter and parameters Hmax, b3, and b4 according to Equation 16.5 of Acevedo (2012).

D.H.est: function performs a set of linear regressions on a log-transformed H versus D allometric relation in order to determine the best coefficient b3 value for each one of a set of values of b4. Then the b3, b4 pair with the smallest residual sum-of-squares error is selected to perform nonlinear regression. Arguments are a  and an array with Dmax and Hmax values. This function uses the previous function DH.allom.

Dl.allom: calculates the leaf area by applying Equation 16.7 of Acevedo(2012).

dD.rate: calculates the increment in the diameter using Equations 16.12 and 16.14 of Acevedo (2012) together with allometric functions DH.allom and Dl.allom already defined. It requires a value or a set of diameter and specification of G and seven parameters. Besides Dmax, three of these are for the H versus D relation and three others for the l versus D relation.

dD.dt.analysis: This function allows for analyzing the dD/dt increment for a given set of values of G. It uses the function dD.rate given above. Besides the values of G, it requires a set of values of D, seven parameters as described above in the function dD.rate, initial D0 condition for the numerical integration, and the time sequence t.

dD.dt.calibra: This function performs an estimation of Gmax, given a set of D and dD values together with Dmax and allometric coefficients. It uses Equation 16.26 in Acevedo (2012) and plots the corresponding growth curve for reference.

parab1: This is a simple function to evaluate a simple parabolic response to an environmental variable according to Equation 16.21 (Acevedo 2012); it takes the value 1 when E = Eopt and uses the argument p for a. Negative values are assigned zero. Although the function is not supposed to receive values of E leading to F greater than 1, we assign one to F values above 1 for caution.

env.factors.analysis: This function is designed to plot the environmental factors for several species. Their number and labels are arguments nsp and splab. Then four response functions (e.g., parab1) are passed as a list FE, Eopt as an array, and param for the functions as a list of arrays.

ge.growth: This function performs a simulation of random environmental fluctuation and uses tree response functions to calculate g(E) according to Equation 16.20. Then Gmaxg(E) is used in conjunction with the allometric parameters and Dmax to simulate tree diameter growth using dD.rate to implement the Euler method. The output includes graphs of g(E) and D versus time for each species.

expon.alf: This function implements the response given by Equation 16.31 (Acevedo, 2012). Arguments are a value of a set of values of environmental variable E and an array p containing the three parameters.

gauss.smf: This function applies Equations 16.32 and 16.33 (Acevedo, 2012) to two moisture indices and calculates the most limiting factor. Arguments are the indices and the tolerance parameters. Both are given as arrays with two components.
}

\value{
DH.allom:
H=height

D.H.est:
\item{b3}{b3 estimate}
\item{b4}{b4 estimate}
\item{H.est}{H values estimates}
\item{RSS.est}{RSS estimate}

Dl.allom
l=leaf area

dD.rate:
dD=diameter increment

dD.dt.analysis:
\item{D.dD}{data.frame(D,dD.dt) diameter and derivative}
\item{t.Xt}{data.frame(t,Dt,Ht,lt), time, diameter, height, and leaf area as function of t}


dD.dt.calibra:
Gmax= max value of G

parab1:
F= response to env factor

env.factors.analysis:
Does not return an object. It plots the results.

ge.growth:
\item{E}{env factors}
\item{F}{response}
\item{X}{data.frame(t,E,G,D)}

expon.alf:
F= response

gauss.smf:
Returns array of response and the most limiting factor
}
\references{
Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.

Phipps, R.L. 1979. Simulation of wetlands forest vegetation dynamics. Ecological Modelling 7:257-288.

Urban, D.L., and H.H. Shugart. 1992. Individual-based models of forest succession. In Plant Succession: Theory and Prediction, eds. D.C. Glenn-Lewin, R.K. Peet, and T.T. Veblan, 249-292. New York: Chapman and Hall.

}
\author{
Miguel F. Acevedo \email{Acevedo@unt.edu}
}
\note{
Nominal parameter values are defined in input files. Variation of param values are defined in lists.

Input files are in 'datafiles.zip' in directory 'datafiles' and organized by chapters of Acevedo (2012). Input files are required to run the examples below.    
}

\seealso{Forest simulation \code{\link{zelig.F}}, \code{\link{forsucc.C}}, \code{\link{semi.F}}  
}
\examples{
\dontrun{
# ---------- DH allometric
  filedata <- "chp16/greenash.txt"
  species.name <- "Green ash"
  # Input file is in two columns, H and DBH , with labels in the header 
  # listed in ascending order of DBH. Maximum H and DBH on the second row  
  D.Hmax <- scan(filedata, skip=1, nlines=1)
  D.H <- read.table(filedata, skip=2); names(D.H) <- c("D","H")
  d.h <- D.H.est(D.H,D.Hmax,species.name)

  Dmax <-100; D <- seq(0,Dmax,0.1)
  G <- c(600,1000)
  param.H <- c(20, -0.02,0.6); param.l <- c(0.16,0.8,2) 
  D0 <- 5; t <- seq(0,300)
  rates.var <- dD.dt.analysis(D,G,Dmax,param.H,param.l,D0,t)

  nsp <- 2; splab <- c("Intolerant","Tolerant")
  FE <- list(parab1,parab1,parab1,parab1)
  Eopt <- c(1,1,1,1)
  param<- list(c(3,2,1.5,1.1), rep(1,4))
  env.factors.analysis(nsp,splab,FE,Eopt,param)

  Dmax <-rep(100,nsp); Gmax <- rep(1000,nsp)
  param.H <- list(c(20, -0.02,0.6),c(20, -0.02,0.6))
  param.l <- list(c(0.16,0.8,2),c(0.16,0.8,2)) 
  D0 <- rep(5,nsp); t <- seq(0,50); nt <- length(t)
  param.E<- list(c(3,2,1.5,1.1), rep(1,4))
  Ef <- list(runif,runif,runif,runif)
  param.Ef <- list(c(Eopt[1],Eopt[1]),c(0,Eopt[2]),c(Eopt[3],Eopt[3]),c(0,Eopt[4]))
  t.GD <- ge.growth(nsp,splab,FE,Eopt,param.E,Ef,param.Ef,Dmax,Gmax,param.H,param.l,D0,nt)

 # calibration Licania densiflora
 dD <- c(0.5,0.6,0.7,0.9,0.8,1.0,1.1,1.3,1.2,1.0,1.0,0.8,0.8)
 D <- seq(15,75,5)
 Dmax <-100
 param.H <- c(40, -0.015,0.6); param.l <- c(0.16,0.6,2.2) 
 dD.dt.calibra(D,dD,Dmax,param.H,param.l)

 # light factor
 clab <- c("Very Int","Int","Med","Tol","Very Tol")
 c1<- c(1.58,1.26,1.13,1.05,1.02)
 c2<- c(1.19,1.79,2.44,3.29,4.16)
 c3<- c(0.15,0.12,0.09,0.06,0.03)
 L <- seq(0,1,0.01); nl <- length(L)
 ALF <- matrix(ncol=5,nrow=nl)
 for(i in 1:5)ALF[,i] <- expon.alf(L, c(c1[i],c2[i],c3[i]))

 # soil moisture factor
 Elab <- c("F(1-E_WD)","F(1-E_DD)")
 x <- seq(0,1,0.01)
 E <- matrix(nrow=length(x),ncol=3)
 y <- matrix(nrow=length(x),ncol=4)
 for(i in 1:2) E[,i] <- x
 p <- c(0.2,0.30)
 for(i in 1:length(x)) 
 y[i,] <- gauss.smf(E[i,],p)

}
}
\keyword{ math }
\keyword{ iteration }
