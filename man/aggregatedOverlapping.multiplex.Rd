\name{aggregatedOverlapping.multiplex}
\alias{aggregatedOverlapping.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Overlapped adjacency matrix of the multiplex network
}
\description{
  Creates the \emph{overlapped adjacency matrix} (see References), which is the adjacency matrix of the \emph{overlay monoplex network} based on a specific multiplex network. It is mostly useful for detecting how many times two different nodes are linked each other on the (intra)layers of the multiplex network (eventually selected with the \code{indexLayer} argument).
}
\usage{
aggregatedOverlapping.multiplex(obj,
                                indexNode = 1:length(nodes.multiplex(obj)),
                                indexLayer = 1:length(layers.multiplex(obj)),
                                verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
}
  \item{indexNode}{
  A vector of IDs (or labels) for the selected nodes on which to calculate the overlapped adjacency matrix. Default selects all the nodes of the network.
}
  \item{indexLayer}{
  A vector of IDs (or labels) for the selected (intra)layers on which to calculate the overlapped adjacency matrix. Default selects all the (intra)layers of the network.
}
  \item{verbose}{
  Default is \code{FALSE}. If \code{TRUE}, the functions prints on console the names of the layers used to create the overlapped adjacency matrix.
}
}

\value{
A N x N \code{matrix}, the so-called overlapped adjacency matrix, where N is the number of nodes of the network. It is simply obtained summing all the adjacency matrices of the selected (intra)layers of the multiplex network.
}
\references{
De Domenico et al (2014). \emph{Mathematical formulation of multilayer networks.} http://doi.org/10.1103/PhysRevX.3.041022
}

\author{
Emanuele Degani \email{emanuele.achab@gmail.com}
}

\seealso{
\code{\link{create.multiplex}}
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

# Creating the overlapped adjacency matrix of the multiplex network:

aggregatedOverlapping.multiplex(mplexObj)

# Creating the overlapped adjacency matrix, casually selecting some nodes (employees) and
# 3 layers (relationships) between the 5 of the multiplex network.
# The command also prints on console (with option 'verbose = TRUE') the name of the layers
# used (and selected) for the construction of the aggregated topological matrix:

aggregatedOverlapping.multiplex(mplexObj,
                                indexNode = sample(1:length(nodes.multiplex(mplexObj)), 10),
                                indexLayer = c(1, 2, 3),
                                verbose = TRUE
                                )
}
