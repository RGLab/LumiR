\name{read.experiment}
\alias{read.experiment}

\title{Read experiment folder}

\description{
Parse the root folder of an experiment to store bead level information in a \code{blum} object.
}

\usage{
read.experiment(path="./")
}

\arguments{
\item{path}{A \code{character}. The pathname of an experiment folder.}
}

\value{
An object of class \code{\linkS4class{blum}} that is the base structure for this package.
}

\note{
The folder passed in path argument must be structured in a specific way and contain mapping files describing the experiment. See the vignette LumiR user guide for more information.
}

\author{Renan Sauteraud}


\seealso{
\code{\linkS4class{blum}}
}

\keyword{function}
