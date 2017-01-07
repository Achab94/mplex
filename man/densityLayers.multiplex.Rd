\name{densityLayers.multiplex}
\alias{densityLayers.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Layers' density of the multiplex network
}
\description{
 Calculates the density of the (intra)layers of the multiplex network, combining the \code{graph.multiplex} and \code{igraph}'s \code{density} functions.
}
\usage{
densityLayers.multiplex(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
  }
}

\value{
  A \code{list} with the densities of the nodes for each (intra)layer of the multiplex network.
}
\references{
Wassermann and Faust (1994), \emph{Social network analysis.} Cambridge University Press.
}

\seealso{
\code{\link{create.multiplex}}, \code{\link[igraph]{graph.density}}}

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

# Calculating the density on the 5 intralayers of the multiplex network.
# We see, for example, that "lunch" and "work" are the denser intralayers, while "coauthor"
# is the less dense.

densityLayers.multiplex(mplexObj)

}
