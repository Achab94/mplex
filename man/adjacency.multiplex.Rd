\name{adjacency.multiplex}
\alias{adjacency.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extraction of the adjacency matrices of the selected (intra)layers
}
\description{
Extracts the adjacency matrices of selected layers in the multiplex network.
}
\usage{
adjacency.multiplex(obj, index = 1:length(obj$adjacency))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \samp{multiplex}.
}
  \item{index}{
   A vector of layer positions whose adjacency has to be computed. Default is all the layers.
}
}

\value{
  A \samp{list} having the adjacency matrices of the selected layers as components.
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

# Extracting a list of the adjacency matrices of the layers (relationships):

adjacency.multiplex(mplexObj)

# Adjacency matrices can be extracted even for some specific layers (for example,
# the first), selected with the 'index' argument:

adjacency.multiplex(mplexObj, index = 1)

# or, easily, referring to their names:

adjacency.multiplex(mplexObj, index = "lunch")
}

\keyword{network, multiplex}
