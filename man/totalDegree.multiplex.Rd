\name{totalDegree.multiplex}
\alias{totalDegree.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Total degrees of the multiplex network
}
\description{
  This function calculates the total degrees (as a sum of the degrees on the single layers) of the multiplex network.
}
\usage{
  totalDegree.multiplex(obj,
                        indexNode = 1:length(nodes.multiplex(obj)),
                        indexLayer = 1:length(layers.multiplex(obj)),
                        verbose = FALSE
                        )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
  }
  \item{indexNode}{
  A vector of IDs (or labels) for the selected nodes on which to calculate the total degree. Default selects all the nodes of the network.
  }
  \item{indexLayer}{
  A vector of IDs (or labels) for the selected (intra)layers on which to calculate the total degree. Default selects all the (intra)layers of the network.
  }
  \item{verbose}{
  Default is \code{FALSE}. If \code{TRUE}, the function prints on console the name of the selected layers on which the sum is evaluated.
  }
}


\value{
  A numeric \code{vector}: each cell contains the total degree for a specific node of the multiplex network.
}
\references{
Wassermann and Faust (1994), \emph{Social network analysis.} http://dx.doi.org/10.1017/CBO9780511815478
}
\author{
Emanuele Degani \email{emanuele.achab@gmail.com}
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

# Calculating the total degrees for the 61 employees of the multiplex network:

totalDegree.multiplex(mplexObj)

# Calculating the total degrees for 10 employees randomly choosen. Note the 'verbose = TRUE'
# option that prints on console the names of the layers used (they can eventually be selected
# with the 'indexLayer' argument):

totalDegree.multiplex(mplexObj,
                      indexNode = sample(1:length(nodes.multiplex(mplexObj)), 10),
                      verbose = TRUE
                      )
}
