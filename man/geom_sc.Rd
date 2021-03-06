% Generated by roxygen2 (4.0.0): do not edit by hand
\name{geom_sc}
\alias{geom_sc}
\alias{geom_sc,slum-method}
\title{Plot the standard curves}
\usage{
geom_sc(object, n = 100, data = NULL, mapping = aes(x = concentration, y =
  mfi), stat = "identity", position = "identity", na.rm = FALSE, ...)
}
\arguments{
\item{object}{A \code{slum}. The object that contains the experiment to plot.}

\item{n}{A \code{numeric}. The number of points used to draw the curves.}

\item{na.rm}{A \code{logical}. Set to TRUE if some of the calculated
concentrations are NA.}

\item{mapping,data,stat,position,...}{Additional parameters to be passed
 to geom_line}
}
\description{
Create a geom_line object to add the standard curves to a plot.
}
\details{
See ggplot2's documentation for more information regarding the additional
arguments.
}
\author{
Renan Sauteraud
}
\seealso{
\code{\link{geom_line}}
}

