\name{Sun functions}

\alias{sun.dec}
\alias{sun.elev.max}
\alias{sun.rad.yr}
\alias{sun.elev.hr}
\alias{sun.rad.hr}
\alias{sun.rad.hr.mult}
\alias{sun.path}
\alias{sun.atmos}

\title{Solar radiation functions}
\description{Functions to calculate sun elevation angle and solar radiation by hour, and maximum daily. Applies to one day or sequence of multiple consecutive days. Also included is azimuth by hour and atmospheric effects.}

\usage{
sun.dec(nday)
sun.elev.max(nday,lat)
sun.rad.yr(nday,lat,Lm)
sun.elev.hr(nday, lat,hr.noon)
sun.rad.hr(nday,lat,hr.noon,Lm,sdr=0)
sun.rad.hr.mult(nday,lat,Lm,sdr=0,sw.plot=T)
sun.path(nday,lat)
sun.atmos(nday,lat,rho)
}

\arguments{
  \item{nday}{day of the year}
  \item{lat}{latitude}
  \item{Lm}{max radiation in the day}
  \item{hr.noon}{hour with respect to noon, negative for hours before noon, and positive for hours after noon}
  \item{sdr}{An optional argument with default value zero is the standard deviation of noise or variability to be imposed on the hourly radiation}
 \item{sw.plot}{optional logical variable to decide to plot the results}
 \item{rho}{reflectivity}
}

\details{
sun.dec: calculates declination for a given day of the year
sun.elev.max: maximum sun elevation angle for a given day
sun.rad.yr: uses the annual average of maximum radiation in the day Lm to convert maximum elevation into radiation.
sun.elev.hr: Sun elevation by hour during the day requires an argument for the hour with respect to noon. Elevation is made zero when the calculation yields
negative values.
sun.rad.hr: uses the annual average of maximum radiation in the day Lm to convert elevation by hour into radiation by hour (optionally can impose variability).
sun.rad.hr.mult:calculates hourly radiation for multiple consecutive days. 
sun.path: calculates sun elevation and azimuth by hour of the day.
sun.atmos: calculates atmospheric effects. 
}

\value{
\item{dec}{declination}
\item{elev}{maximum sun elevation angle for the day or for each hr of the day}
\item{rad}{radiation for the day, for each hour of teh day, or for the running hour in multiple days}
\item{hr.cum}{total running hr for multiple days}
\item{hr.noon}{hour with respect to noon}
\item{azi}{azimuth}
\item{I0}{ET solar radiation kW/m2}
\item{A}{atmospheric effect}
\item{opt.depth}{optical depth}
\item{air.mass}{air mass}
\item{I}{matrix containing IB=direct baeam, ID=diffuse, IR=reflected, and IT=total} 
}

\references{
Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.

Masters G. 2004. Renewable and Efficient Electric Power Systems, Wiley-IEEE Press.
}
\author{
Miguel F. Acevedo \email{Acevedo@unt.edu}
}

\note{
Support for ecosystem simulation
}

\seealso{
Primary productivity functions \code{\link{PPrates.depth}}, \code{\link{PP.Steele}}, \code{\link{DO.PP.pond}}
}

\examples{
\dontrun{
# day of year
nday <- seq(1,365); nd <- length(nday)
dec <- sun.dec(nday)

lat = c(-32,-23.45,23.45,32); nl <- length(lat)
elev <- matrix(nrow=nd,ncol=nl)
for(i in 1:nl) elev[,i] <- sun.elev.max(nday,lat[i])

nday <- seq(1,365); nd <- length(nday)
Lm = c(500,600); nl <- length(Lm)
lat <- 32.90
rad <- matrix(nrow=nd,ncol=nl)
for(i in 1:nl) rad[,i] <- sun.rad.yr(nday,lat,Lm[i])

# latitude tropic, and dallas (DFW)
lat = c(23.45,32.90); nl <- length(lat)
nday = c(15,81,180); nd <- length(nday)
hr.noon <- seq(-12,+12,0.1); nh <- length(hr.noon)
elev <- matrix(nrow=nh,ncol=nl);rad <- matrix(nrow=nh,ncol=nl)
j=1
for(i in 1:nl){
elev[,i] <- sun.elev.hr(nday[j],lat[i],hr.noon)
rad[,i] <- sun.rad.hr(nday[j],lat[i],hr.noon,Lm=600)
}
elev <- matrix(nrow=nh,ncol=nd);rad <- matrix(nrow=nh,ncol=nd)
j=2
for(i in 1:nd){
 elev[,i] <- sun.elev.hr(nday[i],lat[j],hr.noon)
 rad[,i] <- sun.rad.hr(nday[i],lat[j],hr.noon,Lm=600)
}
# latitude dallas (DFW)
lat = 32.90; nday = c(39:41); Lm=600
X <- sun.rad.hr.mult(nday,lat,Lm,sdr=0,sw.plot=TRUE)

}
}

\keyword{ math }
\keyword{ iteration }
