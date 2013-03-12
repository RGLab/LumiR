\name{blum-class}
\alias{blum-class}
\alias{blum}

\title{blum class}

\description{
The \code{blum} class is a container for the raw bead level information of a multiplex Luminex experiment. This is the base structure used in \code{LumiR}.
}

\section{Slots}{
  \describe{
    \item{phenodData}{An object of class \code{AnnotatedDataFrame}. Contains the information regarding the samples (e.g: sample_type, sample_name, well, filename..).}
    \item{featureData}{An object of class \code{AnnotatedDataFrame}. Contains the information regarding the analytes: ID and name.}
    \item{exprs}{A two level \code{list}. Contains the intensities measured for each bead. List(1) of lists(2), with the first level being the samples, and the second being the analytes.}
  }
}

\author{Renan Sauteraud}

\seealso{
\code{read.experiment}
}

\keyword{class}
