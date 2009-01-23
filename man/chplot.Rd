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
plot.points=FALSE,log="",xlab,ylab,col,lty,lwd,legend=TRUE,...)
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
\item{mar.den}{
logical; defines the type of marginal distribution
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
\code{mean.se}.
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
\item{plot.points}{are the points added to the original plot? Default is \code{FALSE}.}

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
\item{lwd}{line width; vector of length equalling the number of groups.}
\item{legend}{
either logical, in which case the default is \code{TRUE} and a legend is drawn with the default settings,
or a list of legend parameters. For legend parameters, see the \code{legend} function; two additional
parameters are allowed: \code{area.in} (logical; specifies whether the area of convex hull per point
is displayed; default is \code{TRUE} if \code{chull=TRUE} and \code{FALSE} otherwise), and \code{pos} (\code{in} or \code{out};
default depends on \code{ratio} value; if \code{in} is chosen, the legend is positioned by the user with mouse,
while \code{out} places the legend in the top right corner of the plotting area).
Note that the default legend parameters are not the same as in the \code{legend} function:
default \code{title} is the name of the grouping variable, default \code{bty} depends on legend position
(\code{"o"} if inside and \code{"n"} if outside), while default \code{cex} is calculated on the basis of \code{ratio},
and the \code{cex} option is the multiplier of that default.
}
\item{...}{other arguments will be passed to the main plotting region and will affect points
if \code{plot.points=TRUE}.}
}

\details{
The relative frequency polygons chosen with \code{mar.den=FALSE}
connect the points in the middle of the intervals and the starting and
ending interval with zero frequency in each group. All the
frequencies are rescaled in order to make the plots immediately
comparable (the same is true for the density curves).

The default density contour might not fit into the plot - this can be
avoided by decreasing the \code{clevel}. In case of broken contours, the
area per point is not a sensible measure.
}

\value{
A list with components:
\item{area}{the area of convex hull per point for each group.}
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
postscript("file.ps",fonts="mono")         # For screen output, skip this
data(hdr)
chplot(age~income|gender,data=hdr,log="x")
# the iris dataset
data(iris)
chplot(Sepal.Length~Sepal.Width|Species,data=iris,bw=TRUE,
       legend=list(cex=.6))
chplot(Petal.Length ~ Petal.Width | Species, data = iris,
      legend = list(cex = 0.6),plot.points=TRUE,pch = 18, cex = 0.5)
dev.off()		# For screen output, skip this
}

\seealso{\code{\link{chadd}}, \code{\link{chull}}, \code{\link{bkde2D}}. }

\keyword{hplot}
