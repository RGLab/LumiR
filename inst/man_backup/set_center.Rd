\name{set_center}
\alias{set_center}
\alias{set_center,blum,character-method}
\alias{set_center,slum,character-method}

\title{Add center information}

\description{
Adds center information to the phenoData and modifies the sample_id of a blum or slum object.
}

\usage{
set_center(object, center_name)
}

\arguments{
\item{object}{An object of class \code{blum} or \code{slum}.}
\item{center_name}{A \code{character}. The name of the center who generated the data.}
}

\value{
An object of the same class as the input object.
}

\author{Renan Sauteraud}


\seealso{
\code{\linkS4class{blum}}
\code{\linkS4class{slum}}
\code{writeMBAA}
}

\keyword{method}
