\name{degreeCentrality.multiplex}
\alias{degreeCentrality.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Multiplex degree centrality indexes of the multiplex network
}
\description{
  Calculates the \emph{multiplex degree centrality} coefficients of the nodes (possibly selected by \code{indexNode} argument) of the multiplex network, defined by De Domenico et al. (see References for further information).
}
\details{The single index referred to a generic node v is defined as its degree on the \emph{projected monoplex network}, whose adjacency matrix is obtained as a sum of all the (intra)layer adjacency matrices (i.e., the \emph{aggregated overlapping multiplex} obtainable with \code{aggregated.overlapping.multiplex}) and the interlayer adjacency matrices.
}
\usage{
degreeCentrality.multiplex(obj, indexNode = 1:length(nodes.multiplex(obj)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
}
  \item{indexNode}{
  A vector of IDs (or labels) for the selected nodes on which to calculate the \emph{multiplex degree centrality} coefficients.
}
}

\value{
  A numerical \code{vector} with all the multiplex degree centrality index for the nodes of the network (eventually selected with \code{indexNode} argument).
}
\references{
De Domenico et al (2014). \emph{Mathematical formulation of multilayer networks.} Phys. Rev. X 3, 041022.
}

\seealso{
  \code{\link{create.multiplex}}, \code{\link{aggregatedOverlapping.multiplex}}
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

# Calculating the degree centrality coefficient of the multiplex network:

degreeCentrality.multiplex(mplexObj)

# It could be useful to calculate this coefficient only on some selected nodes of the network.
# It can be done using the 'indexNode' argument and selecting the appropriate ID of the nodes:

degreeCentrality.multiplex(mplexObj,
                          indexNode = sample(1:length(nodes.multiplex(mplexObj)), 10)
                          )
}

\keyword{network}
\keyword{multiplex}
