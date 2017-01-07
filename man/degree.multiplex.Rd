\name{degree.multiplex}
\alias{degree.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Layers' nodes' degree of the multiplex network
}
\description{
  This function calculates the degree of the nodes on each (intra)layer of the multiplex network, combining the \code{graph.multiplex} and \code{igraph}'s \code{degree} functions. For (intra)layers with \emph{directed} relationships, the \code{modeDirected} argument allows to evaluate 'in-degree', 'out-degree' and 'total degree' descriptors, though by default the function returns only the total degree.
}
\usage{
  degree.multiplex(obj, indexNode = 1:length(nodes.multiplex(obj)), modeDirected = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
  }
  \item{indexNode}{
  A vector of IDs (or labels) for the selected nodes on which to calculate the degree. Default selects all the nodes of the network.
  }
  \item{modeDirected}{
  Default is \code{FALSE}. If \code{TRUE}, the function evaluates 'in-degree' and 'out-degree' besides the usual 'total-degree' descriptor, for directed (intra)layers of the multiplex network.
  }
}

\value{
  A list with, for each (intra)layer, a vector of degrees for the nodes selected with the \code{indexNode} argument. If \code{modeDirected = TRUE}, it adds a list with 'in-degree' and 'out-degree' besides 'total degree' in the directed (intra)layers' position inside the output.
}
\references{
Wassermann and Faust (1994), \emph{Social network analysis.} Cambridge University Press.
}

\seealso{
\code{\link{create.multiplex}}, \code{\link[igraph]{degree}}
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

# Calculating the degrees of the 61 nodes of the multiplex network, for each layer:

degree.multiplex(mplexObj)

# Refining the evaluation only on 10 employees chosen, for each layer:

degree.multiplex(mplexObj, indexNode = 1:10, 10))

# Supposing one (or more) layers of the multiplex network have a directed relationship (for
# example, the third), we can obtain both the in-degree, out-degree and total-degree measures
# using the option 'modeDirected = TRUE':

# 1) Creating the multiplex network with a directed relationship (coauthor, L3)
mplexObjDirected <- create.multiplex(nodes = aarhus_mplex$nodes,
                                    layersNames = aarhus_mplex$layerNames,
                                    layer1 = aarhus_mplex$L1,
                                    type1 = "undirected",
                                    aarhus_mplex$L2,
                                    aarhus_mplex$L3,
                                    "directed",
                                    aarhus_mplex$L4,
                                    aarhus_mplex$L5
                                    )

# 2) Using the 'modeDirected = TRUE' option. Note the output and the list
# nested in $coauthor argument, containing the 3 measures:

degree.multiplex(mplexObjDirected, modeDirected = TRUE)
degree.multiplex(mplexObjDirected, modeDirected = TRUE)$coauthor
}
