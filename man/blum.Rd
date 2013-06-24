\name{blum-class}
\alias{blum-class}
\alias{blum}
\alias{show,blum-method}
\alias{fData,blum-method}
\alias{fData<-,blum,data.frame-method}
\alias{pData,blum-method}
\alias{pData<-,blum,data.frame-method}
\alias{exprs,blum-method}
\alias{merge,blum,blum-method}
\alias{melt}
\alias{melt,blum-method}

\title{blum class}

\description{
The \code{blum} class is a container for the raw bead level information of a multiplex Luminex experiment. This is the base structure used in \code{LumiR}.
}

\section{Slots}{
  \describe{
    \item{phenodData:}{An object of class \code{AnnotatedDataFrame}. Contains the information regarding the samples (e.g: sample_type, sample_name, well, filename..).}
    \item{featureData:}{An object of class \code{AnnotatedDataFrame}. Contains the information regarding the analytes: ID and name.}
    \item{exprs}{A \code{data.table}. Contains the intensities measured for each bead.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "blum")}: Print informations about the object.}
    \item{fData}{\code{signature(object = "blum")}: Accessor for the phenoData slot.}
    \item{pData}{\code{signature(object = "blum")}: Accessor for the featureData slot.}
    \item{exprs}{\code{signature(object = "blum")}: Print the \code{data.table} with the bead events.}
    \item{melt}{\code{signature(x = "blum")}: Convert a `blum' into a data.frame.}
    \item{merge}{\code{signature(x = c("blum","blum"))}: Merge two \code{blum} objects if they can be merged (centers are set and experiment design is identical).}
  }
}

\author{Renan Sauteraud}

\seealso{
\code{read.experiment}
}

\keyword{class}
