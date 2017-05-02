\name{C1Global.multiplex}
\alias{C1Global.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  C1 global clustering coefficients for multiplex networks
}
\description{
  Computes Battiston et al.(2014)'s \emph{C1} global clustering coefficient of a multiplex network.
}
\details{
Instead of the classical definition of \emph{transitivity} and \emph{local clustering coefficient}, this index accounts for all the possible tradic structures that take form considering interlayers relationships. The index is defined as the ratio between the number of 2-triangles of the multiplex network, and (L - 1) times the number of 1-triads of the multiplex network, where L is the number of layers.

Function use \emph{permutations()} function from \emph{gtools} package. More informations available at \url{https://cran.r-project.org/web/packages/gtools/index.html}.
}
\usage{
C1Global.multiplex(obj, indexLayer = 1:length(layers.multiplex(obj)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
}
  \item{indexLayer}{
  A vector of IDs (or labels) for the selected (intra)layers on which to calculate the C1 global clustering coefficient.
}
}

\section{Warnings }{
  Note that the C1 global clustering coefficient can be defined for a multiplex network with, at least, 2 (intra)layers. So, if \code{indexLayer} argument has length 1, the function produces an error and halts.
}

\value{
  A vector of length one and value corresponding to the C1 global clustering coefficient.
}
\references{
Battiston et al. (2014) \emph{Structural measures for multiplex networks.} Phys. Rev. E 89, 032804.
}


\seealso{
  \code{\link{create.multiplex}}, \code{\link{c1Local.multiplex}}
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

# Computing the so-called C1 global clustering index for the nodes of the multiplex network:

C1Global.multiplex(mplexObj)

# It can be useful to select just some layers on which to calculate the
# C1 global clustering index. This can be easily done with 'indexLayer' argument:

C1Global.multiplex(mplexObj, indexLayer = c(2, 4, 5))
}

\keyword{network}
\keyword{multiplex}
