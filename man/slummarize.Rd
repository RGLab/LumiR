\name{slummarize}
\alias{slummarize}

\title{Summarize a blum object}

\description{
Summarize a blum object into a slum object. Fit the standard curves and calculate MFI and the corresponding concentrations.
}

\usage{
slummarize(from, type="MFI")
}

\arguments{
\item{from}{A \code{blum} object. The object to summarize. The object must have information about the standard location and concentration in its \code{phenoData} slot. See the User guide for more information on the requirements.}
\item{type}{A \code{character}. "MFI"}
}

\value{
An object of class \code{slum}.
}

\author{Renan Sauteraud}


\seealso{
\code{\linkS4class{blum}}
\code{\linkS4class{slum}}
}

\keyword{function}
