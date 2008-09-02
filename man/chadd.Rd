\name{chadd}
\alias{chadd}
\title{Add elements to a chplot}
\description{
  \code{chadd} is a function that adds any further elements to a plot
 produced with the function \code{chplot}.
}
\usage{
chadd(param,pos,add.fun,...)
}
\arguments{
\item{param}{the parameters specifying the plotting regions. As obtained from
\code{chplot}.
}
\item{pos}{
the plotting region to which an element is to be added: \code{1} for the main plot,
\code{2} for top left (marginal distribution plot of the x-variable),
\code{3} for top right (the legend),
\code{4} for bottom right (marginal distribution plot of the y-variable).
}
\item{add.fun}{the function to be applied.}
\item{...}{optional parameters to \code{add.fun}, separated by commas.}
}

\examples{
data(hdr)
# hdr dataset
param<-chplot(age~income|gender,data=hdr,log="x")
# box around the legend region
chadd(param,3,box,"figure")
# tickmark for overall mean in marginal distribution plots
chadd(param,2,lines,c(mean(hdr$income),mean(hdr$income)),c(0,.05))
chadd(param,4,lines,c(0,.025),c(mean(hdr$age),mean(hdr$age)))
}
\seealso{\code{\link{chplot}}.}

\keyword{aplot}
