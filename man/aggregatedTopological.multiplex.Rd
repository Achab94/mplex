\name{aggregatedTopological.multiplex}
\alias{aggregatedTopological.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Aggregated topological matrix of the multiplex network
}
\description{
Creates the so-called \emph{aggregated topological matrix} of the multiplex network (see References for further information). This matrix is mostly used to see if one tie exists between two different nodes (eventually selected with the \code{indexNode} argument) in at least one (intra)layer of the multiplex network (eventually selected with the \code{indexLayer} argument).
}
\usage{
aggregatedTopological.multiplex(obj,
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
  A vector of IDs (or labels) for the selected nodes on which to calculate the aggregated topological matrix. Default selects all the nodes of the network.
}
  \item{indexLayer}{
  A vector of IDs (or labels) for the selected (intra)layers on which to calculate the aggregated topological matrix. Default selects all the (intra)layers of the network.
}
  \item{verbose}{
  Default is \code{FALSE}. If \code{TRUE}, the function prints on console the names of the layers used to create the aggregated topological matrix.
}
}

\value{
A N x N \code{matrix}, the so-called aggregated topological matrix, where N is the number of nodes of the network.
}
\references{
Battiston et al. (2014) \emph{Structural measures for multiplex networks.} Phys. Rev. E 89, 032804.
}

\note{
This function could also be used to \emph{dichotomize} a weighted adjacency matrix referred to a particular layer with weighted relationship: select the interested layer with its ID in the \emph{indexLayer} argument.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

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

# Creating the aggregated topological matrix of the multiplex network:

aggregatedTopological.multiplex(mplexObj)

# Creating the aggregated topological matrix, selecting some nodes (employees) and
# 3 layers (relationships) between the 5 of the multiplex network.
# The command also prints on console (with option 'verbose = TRUE') the name of the layers
# used (and selected) for the construction of the aggregated topological matrix:

aggregatedTopological.multiplex(mplexObj,
                                indexNode = 1:10,
                                indexLayer = c(1, 3, 5),
                                verbose = TRUE
                                )
}
