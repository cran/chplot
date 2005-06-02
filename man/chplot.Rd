\name{chplot}
\alias{chplot}


\title{Augmented Convex Hull Plot}
\description{
  Plots 2D convex hulls or bivariate density contours, one for each
  group of data. Descriptive statistics are plotted as error bars or
  confidence ellipses within convex hulls. Marginal distributions
  as well as a special legend are added by default. Additionally, area
  of convex hull per point can be displayed.
}
\usage{
chplot(formula,data,chull=TRUE,clevel=0.95,band.power=.2,
mar.den=FALSE,descriptives="mean.sd",dlevel=0.68,bw=FALSE,ratio=.75,
plot.points=FALSE,log="",xlab,ylab,col,lty,legend,...)
}
\arguments{
\item{formula}{
a formula describing the form of conditioning plot. The
formula is generally of the form \code{y ~ x | g1}, indicating that
 plots of \code{y} (on the y axis) versus \code{x} (on the x axis)
 should be produced conditional on the variable \code{g1}.
 However, the conditioning variable \code{g1} may be omitted.
 The names of the \code{y}, \code{x} and \code{g1} variables are used
 for the axes and legend title.
}
\item{data}{
a data frame containing values for any variables in the formula.
}
\item{chull}{
 logical; if \code{TRUE} (default), convex hulls are plotted,
 otherwise density contours are drawn.
}
\item{clevel}{
the confidence level for the density plot if
\code{chull=FALSE} (default is 0.95).
}
\item{band.power}{
applies if \code{chull=FALSE}; bandwidth for bivariate density estimation
is calculated as \code{(group size)^(-band.power)} for each group
(default is 0.2).
}
\item{mar.den}{ logical; defines the nature of marginal distribution
plots. If \code{FALSE} (default), relative frequency polygons are
plotted, otherwise density plots are drawn.
}
\item{descriptives}{
the option to be used for depicting descriptive statistics. The default value
is \code{mean.sd}, which produces a cross with the lines intersecting at
the mean of \code{x} and \code{y} for each group and depicting the 68-percent tolerance interval
(i.e., stretching one standard deviation in each direction) with default \code{dlevel} setting.
Option \code{mean.se} does the same with standard errors of the means,
while \code{median} makes the lines one quartile long in each direction and
cross at the median. Option \code{ellipse} plots confidence ellipses.
}
\item{dlevel}{
the tolerance/confidence level applied if \code{descriptives=mean.sd} or
\code{mean.sd}.
}
\item{bw}{
logical; if \code{TRUE}, the plot is produced in black-and-white. The
default is \code{FALSE}, which plots in colour.
}
\item{ratio}{
the ratio of the main plot to the whole figure region. The default
value is 0.75. If equal to 1, the marginal distributions are not plotted.
If the ratio is less or equal to 0.75, the default legend position is
outside the main plot (i.e., in the top right corner), otherwise it is
within the main plot (i.e., to be selected by the user with mouse).
}
\item{plot.points}{
logical; if \code{true}, the points are drawn in addition to the convex
hulls or density contours.
}
\item{log}{
 the character strings \code{"x"}, \code{"y"} or \code{"xy"} makes a
 specific (or both) axes logarithmic; the default, which does nothing, is \code{""}.
}
\item{xlab}{
a title for the x axis in the main plot; the default is the name of
the x variable.
}
\item{ylab}{
a title for the y axis in the main plot; the default is the name of
the y variable.
}
\item{col}{the plotting colors; vector of length equalling the number of groups.}
\item{lty}{the type of line; vector of length equalling the number of groups.}
\item{legend}{
object of class \code{legend.control} specifying various parameters for the legend.
 Default is \code{legend.control(...)}. 
}
\item{...}{other arguments will be passed to \code{legend.control}}.
}

\details{
The relative frequency polygons chosen with \code{mar.den=FALSE}
connect the points in the middle of the intervals and the starting and
ending interval with zero frequency in each group. All the
frequencies are rescaled in order to make the plots immediately
comparable (the same is true for the density curves).

The default density contour might not fit into the plot - this can be
avoided by decreasing the \code{clevel}. In such cases the
area does not provide a sensible measure.
}

\value{
A list with components:
\item{usrc}{the limits of the central plotting region.}
\item{usru}{the limits of the top left plotting region.}
\item{usrr}{the limits of the bottom right plotting region.}
\item{ratio}{the ratio used.}
\item{is.xlog}{logical, denoting whether either \code{log="x"}
 or  \code{log="xy"} was used.}
\item{is.ylog}{logical, denoting whether either \code{log="y"}
 or  \code{log="xy"} was used.}
}


\references{
Vidmar, G., and Pohar, M. Augmented convex hull plots: rationale, implementation in R and
biomedical applications. \emph{Computer Methods and Programs in Biomedicine}, 2005, 78, 69-74.
}

\examples{
# the hdr dataset
data(hdr)
chplot(age~income|gender,data=hdr,log="x")
# the iris dataset
data(iris)
chplot(Sepal.Length~Sepal.Width|Species,data=iris,bw=TRUE,
       legend=legend.control(cex=.6))
}

\seealso{\code{\link{chadd}}, \code{\link{chull}}, \code{\link{bkde2D}}, \code{\link{legend.control}}. }

\keyword{hplot}
