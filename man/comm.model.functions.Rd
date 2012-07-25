\name{Community model functions}

\alias{LVint.2sp}
\alias{holling.tanner}
\alias{LVint.3sp}
\alias{succession}
\alias{succession.z}
\alias{succession.g}
\alias{eq.2sp}
\alias{diversity}

\title{Functions to define community dynamics models for simulation, to determine equilibrium points, and calculate diversity from simulation output.}

\description{
Model functions are used as part of argument model to call simulation functions. Model functions are: LVint.2sp, holling.tanner, LVint.3sp, succession, succession.z and succession.g. Function eq.2sp and diversity support other calculations.
}

\usage{
LVint.2sp(t,p,x)
holling.tanner(t,p,x) 
LVint.3sp(t,p,x)
succession(t,p,x)
succession.z(t,p,x)
succession.g(t,p,x,tz)
eq.2sp(A,r)
diversity(output,param,fileout,pdfout) 
}

\arguments{
  \item{t}{time}
  \item{p}{parameters, an array}
  \item{x}{state variable}
  \item{tz}{times of sudden disturbance}
  \item{A}{matrix of LV system}
  \item{r}{growth rate coefficients}
  \item{param}{replicated content of input file}
  \item{output}{simulation output}
  \item{fileout}{output filename}
  \item{pdfout}{logical variable to decide whether to generate pdf output file}
}

\details{
State variables are population densities of species in the community.

LVint.2sp and LVint.3sp are for systems of two and three species in interaction. It is similar to the function used for compartments.

The succession model function is based on compartment models; the value in the diagonal is computed, taking into account that the sum of elements in a column should be 0.

Functions succession.z and succession.g define discontinuous disturbance regime. Function succession.z defines times of discontinuities, and succession.g applies a linear disturbance at tz times defined by succession.z. In function succession.g, we impact stages 2-5 negatively by the proportion given in u and then add all the cover taken to stage 1 to reset succession.

Function eq.2sp is a tool to check the condition for a nontrivial coexistence equilibrium point and gives the community matrix C=diag(xeq)A properties if the equilibrium exists.

Function diversity calculates proportions of population density and Shannon Index H. 

}

\value{
Model functions: Rate of change or derivative of model.

Function eq.2sp:
\item{detpos}{logical variable; TRUE if determinant is positive}
\item{coex}{logical variable; TRUE if there is nontrivial coexistence}
\item{xeq}{equilibrium point when there is coexistence}
\item{Det.Tra.Disc}{determinant, trace and discriminant of community matrix C}
\item{eval}{eigenvalues of community matrix C}

Function diversity: data frame with t and H(t)
}
\references{
Abrams, P.A., and J.D. Roth. 1994. The effects of enrichment of three species food chains with nonlinear functional responses. Ecology 75:1118-1130.

Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.

Acevedo, M. F. 1980. Electrical network simulation of tropical forests successional dynamics. In Progress in Ecological Engineering and Management by Mathematical Models, ed. D. Dubois, 883-892. Liege, Belgium: Centre Belge d'etudes et de documentation (Belgian Center of Studies and Documentation).

Bledsoe, L.J., and G.M. Van Dyne. 1971. A compartment model simulation of secondary succession. In Systems Analysis and Simulation in Ecology, Vol. 1, ed. B.C. Patten, 479-511. New York: Academic Press.

Dambacher, J.M., H.-K. Luh, H.W. Li, and P.A. Rossignol. 2003. Qualitative stability and ambiguity in model ecosystems. American Naturalist 161:876-888.

Hallam, T.G. 1986. Community dynamics in homogeneous environments. In Mathematical Ecology, ed. T. G. Hallam and S. A. Levin, 241-285. New York: Springer-Verlag.

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
# competition co-existence
A = matrix(c(-0.1,-0.01,-0.01,-0.1),ncol=2,byrow=T)
r = c(0.15,0.1)
eq.2sp(A,r)

LV.2sp <- list(f=LVint.2sp)
param <- list(plab=c("X1.0","X2.0"), pval = cbind(c(0.7,0.7,2.5,2.5),c(0.5,1.2,1.0,1.2)))
t.X <- sim(LV.2sp,"chp12/comp-inp.csv", param)

# holling-tanner
holl.tann <- list(f=holling.tanner)
param <- list(plab="X2.0", pval = seq(0.1,0.2,0.05))
t.X <- sim(holl.tann,"chp12/holl-tann-inp.csv", param,pdfout=TRUE)

#  three species competition co-existence
LV.3sp <- list(f=LVint.3sp)
param <- list(plab="h2", pval=c(0,-0.2,-0.5)) 
t.X <- sim(LV.3sp,"chp12/comp-3sp-inp.csv",param,pdfout=TRUE)
t.H <- diversity(t.X$output,param,fileout="chp12/comp3sp-div") 

# succession
succ <- list(f=succession,z=succession.z,g=succession.g)
param <- list(plab="u", pval=c(0,-0.4)) 
t.X <- simd(succ,"chp12/succession-inp.csv", param)
t.H <- diversity(t.X$output,param,fileout="chp12/succession-2u-div",pdfout=TRUE) 
}
}
\keyword{ math }
\keyword{ iteration }
