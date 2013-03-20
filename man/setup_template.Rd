\name{setup_templates}
\alias{setup_templates}

\title{Read experiment folder}

\description{
Creates mapping file to be parsed in `read.experiment'.
}

\usage{
setup_templates<-function(path, templates=c("layout", "analyte", "phenotype"))
}

\arguments{
\item{path}{A \code{character}. The pathname of an experiment folder. This is the same path as the one that will be used in \code{read.experiment}.}
\item{templates}{A \code{character} vector. The mapping files to create.}
}

\value{
An invisible list of \code{data.frame} to be used internally by \code{read.experiment} if one or more mapping files are missing from the given path.
}

\note{
By if a mapping file specified by the templates list already exists, the function will not overwrite the existing file(s). The files should be removed before running the function in order to create new mapping file.
}

\author{Renan Sauteraud}


\seealso{
\code{read.experiment}
\code{\linkS4class{blum}}
}

\keyword{function}
