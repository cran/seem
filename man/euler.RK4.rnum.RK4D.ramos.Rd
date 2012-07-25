\name{Numerical integration}
\alias{euler}
\alias{RK4}
\alias{rnum}
\alias{RK4D}
\alias{ramos}

\title{
Numerical integration methods: Euler, Runge-Kutta 4th order, Monte Carlo, RK4 with discontinuities, and Ramos method.
}

\description{
Function euler is used only for didactic purposes. Uses the Euler method to numerically solve for X in the model dX/dt = f(t,p,X). 

Function RK4 uses the RK 4th order method to numerically solve for X in the model dX/dt = f(t,p,X)

Function rnum is based on the euler function but includes many realizations based on the mean and standard deviation of the parameter and calculations of mean and standard deviation of the result of all realizations.

Function RK4D explicitly includes the discontinuity function g in the integration loop.

Function ramos implements the Ramos nonstandard explicit integration algorithm (EIA) that requires df/dx in addition to f.
}

\usage{
euler(x0, t, f, p, dt)
RK4(x0, t, f, p, dt)
rnum(x0, t, f, p, dt, n)
RK4D(x0, t, f, p, dt, g, tz)
ramos(x0, t, f, p, dt)
}
\arguments{
  \item{x0}{initial condition}
  \item{t}{times for output}
  \item{f}{model function, right hand side of ODE}
  \item{p}{parameter}
  \item{dt}{time step}
  \item{n}{number of realizations}
  \item{g}{discontinuity function, passed as component of model list}
  \item{tz}{times of discontinuities, calculated from component *.z of model list}
}
\details{
Function euler: declare a matrix to store X, based on the number of entries of time and the number of entries in the variable as given by the initial conditions. Then initialize X and time, and execute two nested for loops. The outer loop controls the times to save output, whereas the inner loop runs the calculation in steps of dt updating X according to equation 4.6 (Acevedo 2012). Note that these loops include more than one line and therefore their scope requires curly braces. For the purposes of models covered in Acevedo (2012), all variables X will be positive or zero, so we force X to zero when negative.

Function RK4: essentially the same function as euler, but implements the Runge Kutta 4th order method instead of Euler inside the inner loop. First declare a matrix to store X, based on the number of entries of time and the number of entries in the variable as given by the initial conditions. Then initialize X and time, and execute two nested for loops. The outer loop controls the times to save output, whereas the inner loop runs the calculation in steps of dt updating X according to equations 4.8 and 4.9 (Acevedo 2012). For the purposes of models covered in Acevedo (2012), all variables X will be positive or zero, so we force X to zero when negative.

Function rnum is based on the euler function but includes the use of a three-dimensional array X, defined using structure, and the calculations of mean and standard deviation of all realizations. The third dimension of array X is for the realizations.

Function RK4D explicitly includes the discontinuity function g in the integration loop. This strategy separates the change of rate due to discontinuities from the continuous part, thus allowing more generality when we work with different sets of discontinuities.

Function ramos implements the Ramos nonstandard explicit integration algorithm (EIA) that requires df/dx in addition to f.
}
\value{
Results of integrated state variable at times specified by argument t. Function rnum returns mean and standard deviation of the realizations.
}
\references{
Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.

Ramos, H. 2007. A non-standard explicit integration scheme for initial-value problems. Applied Mathematics and Computation 189:710-718.

Swartzman, G.L., and S. Kaluzny. 1987. Ecological Simulation Primer. New York: MacMillan.
}
\author{
Miguel F. Acevedo \email{Acevedo@unt.edu}
}
\note{
These functions are employed by the simulation functions, therefore typically the end user will not work directly with the numerical functions but with the simulators. 
}

\seealso{
Simulation functions \code{\link{sim.comp}}, \code{\link{sim.rnum}}, \code{\link{sim.mruns}}, \code{\link{sim}}, \code{\link{simd}}, \code{\link{simr}}, Model functions \code{\link{expon}}, \code{\link{expon.z}}, \code{\link{expon.g}}, \code{\link{logistic}}, and many others.
}
\examples{
 # Euler
 model <- list(f=expon)
 t <- seq(0,10,1); dt <- 0.001
 p <- 0.1; X0 <- 1
 X <- euler(X0, t, model$f, p, dt)
 # Runge-Kutta
 model <- list(f=expon)
 t <- seq(0,10,1); dt <- 0.001
 p <- 0.1; X0 <- 1
 X <- RK4(X0, t, model$f, p, dt)
 # Stochastic
 model <- list(f=expon)
 t <- seq(0,10,1); dt <- 0.001
 p <- c(0.1,0.01); X0 <- 1
 X <- rnum(X0, t, model$f, p, dt,n=20)
 # RK4 with discontinuities
 model <- list(f=expon,z=expon.z,g=expon.g)
 t <- seq(0,100,1); dt <- 0.01
 p <- c(0.02,10,0,-10); X0 <- 100
 X <- RK4D(X0, t, model$f, p, dt, model$g, model$z(t,p,X))
}
\keyword{ math }
\keyword{ iteration }
\keyword{ programming }
