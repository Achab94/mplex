\name{varianceDegree.multiplex}
\alias{varianceDegree.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Layers' variance degree of the multiplex network
}
\description{
  This function calculates the variance of the degrees on each (intra)layer of the multiplex network. For (intra)layers with \emph{directed} relationships, the function considers the 'total degree' descriptor.
}
\usage{
varianceDegree.multiplex(obj,
                        indexNodeVar = 1:length(nodes.multiplex(obj)),
                        verbose = FALSE
                        )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
  }
  \item{indexNodeVar}{
  A vector of IDs (or labels) for the selected nodes on which to calculate the variance degree. Default selects all the nodes of the network.
  }
  \item{verbose}{
  Default is \code{FALSE}. If \code{TRUE}, the function prints on console the name of the selected nodes on which the variances are evaluated.
  }
}

\value{
  A \code{list} containing the variances of the degrees for each (intra)layer of the multiplex network.
}
\references{
Wassermann and Faust (1994), \emph{Social network analysis.} Cambridge University Press.
}

\seealso{
\code{\link{create.multiplex}}, \code{\link{degree.multiplex}}
}
\examples{
# Loading Aarhus CS Department dataset.

data(aarhus_mplex)

# Creating the multiplex object using the dataset loaded into aarhus_mplex object.

mplexObj <- create.multiplex(nodes = aarhus_mplex$nodes,
                            layersNames = aarhus_mplex$layerNames,
                            layer1 = aarhus_mplex$L1,
                            type1 = "undirected",
                            aarhus_mplex$L2,
                            aarhus_mplex$L3,
                            aarhus_mplex$L4,
                            aarhus_mplex$L5
                            )

# Calculating the variance degrees on the 5 intralayers of the multiplex network.
# We see, for example, that "facebook" is the more varied relationship (in terms of degree
# distribution), while "coauthor" is the less varied.

varianceDegree.multiplex(mplexObj)

# Calculating the variance degrees for the first 10 employees:

varianceDegree.multiplex(mplexObj, indexNodeVar = 1:10)

# Using the 'verbose = TRUE' option, the function prints the name of the nodes on which the variance
# degree is performed:

varianceDegree.multiplex(mplexObj, indexNodeVar = 1:10, verbose = TRUE)
}

\keyword{network, multiplex}
