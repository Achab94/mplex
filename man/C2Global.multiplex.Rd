\name{C2Global.multiplex}
\alias{C2Global.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  C2 global clustering coefficients for multiplex networks
}
\description{
  Computes Battiston et al. (2014)'s \emph{C2} global clustering coefficient of a multiplex network. Instead of the classical definition of \emph{transitivity} and \emph{local clustering coefficient}, this index takes into consideration all the possible tradics structures that takes form even considering interlayers relationships. The index is defined as the ratio between the number of 3-triangles of the multiplex network, and (L - 2) times the number of 2-triads of the multiplex network, where L is the number of layers.
}
\usage{
C2Global.multiplex(obj, indexLayer = 1:length(layers.multiplex(obj)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
}
  \item{indexLayer}{
  A vector of IDs (or labels) for the selected (intra)layers on which to calculate the C2 global clustering coefficient.
}
}

\section{Warnings }{
  Note that the C2 global clustering coefficient can be defined for a multiplex network with, at least, 3 (intra)layers. So, if \code{indexLayer} argument has length < 3, the function produces an error and halts.
}

\value{
  A vector of length one and value corresponding to the C2 global clustering coefficient.
}
\references{
Battiston et al. (2014) \emph{Structural measures for multiplex networks.} Phys. Rev. E 89, 032804.
}
}

\seealso{
  \code{\link{create.multiplex}}, \code{\link{c2Local.multiplex}}
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

# Computing the so-called C2 global clustering index for the nodes of the multiplex network:

C2Global.multiplex(mplexObj)

# It can be useful to select just some layers on which to calculate the
# C2 global clustering index. This can be easily done with 'indexLayer' argument:

C2Global.multiplex(mplexObj, indexLayer = c(1, 2, 3))
}
