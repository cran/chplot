\name{legend.control}
\alias{legend.control}


\title{Options for legend in chplot}
\description{Checks and packages the legend options \code{chplot}
}
\usage{
legend.control(include=TRUE,area.in,pos,cex,bty,title)
}
\arguments{
\item{include}{
logical; if \code{TRUE} (default), a legend is created.
}
\item{pos}{
forces the position of the legend. The available options are \code{in} and \code{out}; the default depends on \code{ratio} value.
If \code{in} is chosen, the legend is positioned by the user with mouse; \code{out} places the legend in the top right corner of
the plotting area.
}
\item{area.in}{
logical; if \code{TRUE}, the area per point is calculated for each group and included in
the legend. The default is \code{TRUE} if \code{chull=TRUE} and \code{FALSE} otherwise.
}
\item{cex}{
a numerical value giving the amount by which text in the legend should be scaled relative to the default.
}
\item{bty}{
the type of box to be drawn around the legend. See \code{par}.
}
\item{title}{
the title of the legend; the default is the name of the grouping variable.
}
}

\value{
A list with the same elements as the input.
}

\seealso{\code{\link{chplot}}.}

\keyword{hplot}
