\name{nodes.multiplex}
\alias{nodes.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extraction of the nodes' vector
}
\description{
  The function allows a user-friendly method for extracting the information about selected nodes of the multiplex network.
}
\usage{
nodes.multiplex(obj, index = 1:nrow(obj$nodes), label = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
}
  \item{index}{
  A vector of the nodes' ID selected for the extraction. By default, the vector contains the IDs of all the nodes of the network.
  }
  \item{label}{
  \code{FALSE} by default. If \code{TRUE}, the function extracts only the labels of the selected  nodes.
  }
}

\value{
A \code{vector} with the nodes selected by the \code{index} argument.
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

# Extracting the information about the nodes (employees) of the network:

nodes.multiplex(mplexObj)

# Information can be extracted even for some specific nodes (their selection is
# randomized in the following example), selected with the 'indexNode' argument,
# on their label selected with 'label' argument:

nodes.multiplex(mplexObj, index = sample(1:length(nodes.multiplex(mplexObj)), 10), label = TRUE)
}
