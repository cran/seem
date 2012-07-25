\name{Simulator output functions}
\alias{onerun.out}
\alias{mruns.out}
\alias{rnum.out}
\alias{vars.out}

\title{
Simulator output functions 
}
\description{
Produces output files for single runs and one variable (onerun.out),and multiple runs and one variable (mruns.out), or single/multiple runs plus rates (vars.out). Function rnum.out is similar to the onerun.out but plots the mean and the upper and lower bounds according to standard deviation of the realizations.
}

\usage{
onerun.out(prefix, lab.out = "", input, output, pdfout = F)
mruns.out(prefix, lab.out = "", input, param, output, tXlab, pdfout = F)
rnum.out(prefix, lab.out = "", input, output, tXlab, pdfout = F)
vars.out(prefix, lab.out="", input, param, output, tXlab, pdfout=F)
}

\arguments{
  \item{prefix}{character string to form the names of the output files}
  \item{lab.out}{character string to follow the prefix when forming the names of the output files}
  \item{input}{replicated content of input file}
  \item{param}{replicated content of input file}
  \item{output}{simulation output}
  \item{tXlab}{labels for graphs including units}
  \item{pdfout}{logical variable to decide whether to generate pdf output file}
}
\details{
Function onerun.out, it extracts t and X and their labels from the output data frame. Then, plots X vs. t, storing the graph in a pdf file when pdfout is T; the pdf filename is formed by concatenating the prefix, the output label passed as an argument, the label -out, and the extension pdf. Subsequently, it writes input and output data to a csv output file named in the same fashion as the pdf file; it uses function write.table with a variety of options; and transposes the names of the output to a row before writing them.

Output function mruns.out is similar to onerun.out but has additional arguments tXlab for labeling graphs with units, and param for the varying parameter. It uses matplot instead of plot, in order to draw several lines on one graph, one line per run. The output file name includes one more line in order to provide information about the varying parameter.

Function rnum.out, used by sim.rnum, is similar to the onerun.out but plots the mean and the upper and lower bounds. It requires the standard deviation provided by sim.rnum.

Function vars.out is long but that is the price we pay so that it can be generic and reusable for many model functions with two or more variables and simulations of two or more runs. In addition, it plots Xi+1 versus Xi to describe phase-plane portraits of the system. To make the legend more readable, it expands the horizontal axes to 50 percent more than the maximum value of t and start the legend at the top right. For more details, the interested user can examine the function on the console.
}
\value{
These functions do not return objects, they are called by simulator functions to write output files.
}
\references{
Acevedo M.F. 2012. Simulation of Ecological and Environmental Models. CRC Press.
}
\author{
Miguel F. Acevedo \email{Acevedo@unt.edu}
}
\note{
These functions are used by sim.comp, sim.mruns, sim.rnum, and sim.
}
\seealso{
Simulation functions \code{\link{sim.comp}}, \code{\link{sim.mruns}}, \code{\link{sim.rnum}}, \code{\link{sim}}, \code{\link{simd}}, \code{\link{simr}}
}
\examples{
\dontrun{ The end user would not typically have to call the function. Thus, no examples are provided.}
}

\keyword{ IO }
\keyword{ file }
\keyword{ dplot }
