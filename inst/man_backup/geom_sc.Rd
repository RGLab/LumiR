\name{geom_sc}
\alias{geom_sc}
\alias{geom_sc,slum-method}

\title{Plot the standard curves}

\description{
Create a geom_line object to add the standard curves to a ggplot.
}

\usage{
geom_sc(object, n=100, data=NULL,  mapping = aes(x=concentration, y=mfi),
                   stat = "identity", position = "identity", 
                   na.rm = FALSE, ...)
}

\arguments{
\item{object}{An object of class \code{slum}.}
\item{n}{A \code{numeric}. The number of points used to draw the curves.}
  \item{mapping}{The aesthetic mapping, usually constructed
  with \code{aes} or \code{aes_string}. Only
  needs to be set at the layer level if you are overriding
  the plot defaults.}

  \item{data}{A layer specific dataset - only needed if you want to
  override the plot defaults.}

  \item{stat}{The statistical transformation to use on the 
  data for this layer.}

  \item{position}{The position adjustment to use for 
  overlappling points on this layer}

  \item{na.rm}{Remove 'NA'. Argument passed to \code{geom_line}.}

  \item{...}{other arguments passed on to
  \code{layer}. This can include aesthetics whose
  values you want to set, not map. See \code{layer}
  for more details.}

}

\note{
See ggplot2's documentation for \code{geom_line} for more information.
}

\author{Renan Sauteraud}


\seealso{
\code{ggplot2::geom_line}
\code{\linkS4class{slum}}
}

\keyword{method}
