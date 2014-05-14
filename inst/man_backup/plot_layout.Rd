\name{plot_layout}
\alias{plot_layout}
\alias{plot_layout,blumORslum-method}

\title{Plot a plate layout}

\description{
Visualisation of a plate's sample information per well.
}

\usage{
plot_layout(object, plate_name=NULL, ...)
}

\arguments{
\item{object}{An object of class \code{blum} or \code{slum}.}
\item{plate_name}{A \code{character}. The name of the plate to visualize.}
\item{...}{Additional arguments to be passed to the underlying \code{geom_polygon} function.}
}

\value{
An invisible object of class \code{ggplot}.
}

\note{
The information passed in the additional arguments should be available in the phenoData.
}

\author{Renan Sauteraud}


\seealso{
\code{ggplot2::geom_polygon}
\code{\linkS4class{blum}}
\code{\linkS4class{slum}}
}

\keyword{method}
