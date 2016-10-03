\name{aarhus_mplex}
\alias{aarhus_mplex}
\docType{data}
\title{
  Dataset of Computer Science Department employees at Aarhus University
}
\description{
  This multiplex social network consists of five kinds of online and offline relationships (Facebook, Leisure, Work, Co-authorship, Lunch) among 61 employees of the Computer Science Department at Aarhus.

  The original dataset has been re-arranged into a \code{list} object for an easily usage in the examples of the functions contained in this package.
}
\usage{data(aarhus_mplex)}
\format{
  Data are organized into a list of 7 objects:
  \itemize{
    \item{$nodes: a \code{data.frame} object with the information about the employees (nodes) of the CS Aarhus Department.}
    \item{$layerNames: a \code{vector} with names of the 5 levels of the multiplex network.}
    \item{$L1: a \code{data.frame} object with the edges of the first relationship (lunch) considered.}
    \item{$L2: a \code{data.frame} object with the edges of the second relationship (facebook) considered.}
    \item{$L3: a \code{data.frame} object with the edges of the third relationship (coauthor) considered.}
    \item{$L4: a \code{data.frame} object with the edges of the fourth relationship (leisure) considered.}
    \item{$L5: a \code{data.frame} object with the edges of the fifth relationship (work) considered.}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
 The original dataset is publicly available at \url{http://deim.urv.cat/~manlio.dedomenico/data/CS-Aarhus_Multiplex_Social.zip}.
 See \url{http://larica.uniurb.it/sigsna/data/} for further details.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(aarhus_mplex)
# Loads the dataset into the current environment.

head(aarhus_mplex$nodes)
# data.frame object with information about the node (employees) of the multiplex network.
# We see, for example, that employee named 'U33' has number 3 as univoque ID.

aarhus_mplex$layerNames
# list object with characters describing the layers (relationships) taken into consideration
# in the multiplex network.

head(aarhus_mplex$L3)
# data.frame object with edges for coautorship (third level) relationship.
# We see, for example, that employees with IDs 18 and 46 are linked into a coautorship
# relation in CS Aarhus Department.
}
\keyword{datasets}
