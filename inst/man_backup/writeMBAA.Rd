\name{writeMBAA}
\alias{writeMBAA}

\title{
Write a slum object as MBAA_results.xls
}

\description{
Format a \code{slum} object into a table as specified by ImmPort: 'https://immport.niaid.nih.gov/example_submission_packages/MBAA_Results.xls'.
}

\usage{
writeMBAA(object, outfile="./MBAA_results", type="csv", concentration_unit="pg/mL")
}

\arguments{
\item{object}{A \code{slum} object. The summarized experiment to export.}
\item{outfile}{A \code{character}. The name of the output file without extension.}
\item{type}{A \code{character}. Can be "csv" or "xls". Determines the type of file to be written. This will also be appended as the extension after the filename given in outfile.}
\item{concentration_unit}{A \code{character}. The unit of the concentrations to be reported in the document.}
}

\author{Renan Sauteraud}


\seealso{
\code{\linkS4class{slum}}
}

\keyword{function}
