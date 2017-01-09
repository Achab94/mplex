\name{participationDegree.multiplex}
\alias{participationDegree.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Multiplex participation coefficients of degrees' distribution
}
\description{
  This function calculates the \emph{multiple participation coefficient} as a measure used to evaluate if a node is \emph{truly multiplex} (if it belongs to all the (intra)layers of the multiplex network) or \emph{focused} (if it belongs to just one (intra)layers of the multiplex networks). The formula is proposed in Battiston et al. (2014) work (see References) and assumes values in a range between 0 (if the node is \emph{focused}) and 1 (if the node is \emph{truly multiplex}).
}
\usage{
participationDegree.multiplex(obj,
                              indexNode = 1:length(nodes.multiplex(obj)),
                              indexOverlappingLayer = 1:length(layers.multiplex(obj)),
                              display = FALSE
                              )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
}
  \item{indexNode}{
  A vector of IDs (or labels) for the selected nodes on which to calculate the multiple participation coefficient.
}
  \item{indexOverlappingLayer}{
  A vector of IDs (or labels) for the selected (intra)layers on which to calculate the \emph{aggregated overlapping matrix} with \code{aggregatedOverlapping.multiplex} function, used inside multiplex partitipation coefficient's formula. Default selects all the layers of the network.
}
  \item{display}{
  Default is \code{FALSE}. If \code{TRUE}, the function plots a user-friendly visualization of the multiplex partitipation coefficients distribution that helps doing comparisons with others multiplex networks.
}
}
\details{
  This function takes values in (0, 1). Defined P(v_i) the multiplex participation coefficient of node v_i, then v_i is said to be \emph{focused} if \eqn{0 < P(v_i) < 1/3}, \emph{mixed} if \eqn{1/3 < P(v_i) < 2/3} or \emph{truly multiplex} if \eqn{2/3 < P(v_i) < 1}.
}
\value{
  A numeric \code{vector} with all the multiplex participation coefficients evaluated on the nodes of the multiplex network (eventually selected with the \code{indexNode} argument).
}
\references{
Battiston et al. (2014) \emph{Structural measures for multiplex networks.} Phys. Rev. E 89, 032804.
}

\seealso{
  \code{\link{create.multiplex}}, \code{\link{aggregatedOverlapping.multiplex}}, \code{\link[igraph]{degree}}
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

# Calculating the multiplex participation coefficients' distribution:

participationDegree.multiplex(mplexObj)

# A better way to visualize the results is given by the 'display = TRUE' option:

participationDegree.multiplex(mplexObj, display = TRUE)

# It could be also useful to select just some nodes or levels on which to calculate the
# multiplex participation coefficient, in order to investigate which layer(s) is (are)
# causing a low measure for a certain node.
# This can be easily done with both 'indexNode' and 'indexOverlappingLayer' arguments:

participationDegree.multiplex(mplexObj,
                              indexNode = sample(1:length(nodes.multiplex(mplexObj)), 10),
                              indexOverlappingLayer = (1:5)[-4],
                              display = TRUE
                              )
}

\keyword{network}
\keyword{multiplex}
