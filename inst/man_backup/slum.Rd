\name{slum-class}
\alias{slum-class}
\alias{slum}
\alias{exprs,slum-method}
\alias{fit}
\alias{fit,slum-method}
\alias{concentration}
\alias{concentration,slum-method}
\alias{melt,slum-method}
\alias{formula<-}
\alias{formula<-,slum,formula-method}

\title{slum class}

\description{
The \code{slum} class is a summarized version of \code{blum}. The individual bead fluorescence measurements are replaced by the mfi of each analyte for each sample.
}

\section{Slots}{
  \describe{
    \item{phenodData:}{An object of class \code{AnnotatedDataFrame}. Contains the information regarding the samples (e.g: sample_type, sample_name, well, filename..).}
    \item{featureData:}{An object of class \code{AnnotatedDataFrame}. Contains the information regarding the analytes: ID and name.}
    \item{assayData:}{An \code{environment}. Contains two identically dimensioned matrices with analyte names as rows and sample_id as columns.} 
    \describe{
      \item{	exprs:}{The \code{matrix} of the MFI.}
      \item{	concentration:}{The \code{matrix} of the concentration calculated for the MFIs in exprs.}
    }
    \item{fit:}{A \code{data.frame}. Contains the information regarding the standard curve fitting: The values of the parameters to use in \code{inv} for each couple analyte / sample_id.}
    \item{formula:}{The formula used for the standard curve fitting. By default a 5-Parameters Logistic.}
    \item{inv:}{A \code{function}. The inverse function of the formula.}
    \item{unit:}{A \code{character}. "MFI"}
    \item{protocolData:}{An optional \code{AnnotatedDataFrame}, inherited from \code{ExpressionSet}.}
    \item{experimentData:}{An optional \code{MIAME}, inherited from \code{ExpressionSet}.}
  }
}

\section{Extends}{
  Class \code{\linkS4class{ExpressionSet}}, directly.
}

\section{Methods}{
  \describe{
    \item{exprs}{: \code{signature(object = "slum")}: Print the \code{matrix} of the MFI.}
    \item{concentration}{: \code{signature(object = "slum")}: Print the \code{matrix} of the concentrations.}
    \item{fit}{: \code{signature(object = "slum")}: Print the \code{data.frame} of the standard curves parameters.}
    \item{melt}{: \code{signature(x = "slum")}: Convert a `slum' into a data.frame.}
    \item{formula<-}{: \code{signature(object = "slum", value = "formula")}: Setter for the formula.}
  }
}

\author{Renan Sauteraud}

\seealso{
\code{slummarize}
}

\keyword{class}
