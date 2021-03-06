\name{heterEigenvectorCentrality.multiplex}
\alias{heterEigenvectorCentrality.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Multiplex heterogeneous eigenvector centrality indices of the multiplex network
}
\description{
 The function calculates the \emph{multiplex heterogeneous eigenvector centrality indices} proposed by Sola' et al.
}
\details{
This measure refeirs to the eigenvector of a special N*L x N*L matrix obtained from the (inter)layer adjacency matrices and the so-called \emph{influence matrix} of the multiplex network, where N is the number of nodes and L the number of layers of the multiplex network taken into consideration. Its eigenvector has length N*L, thus it's breaked into L sub-vectors, each of them refers to a certain (intra)layer of the multiplex network.

Irreducibility is a required assumption to satisfy the Perron-Frobenius theorem, which ensures the positivity of the eigenvector assosicated to the maximum eigenvalue of the supra adjacency matrix of the multiplex network; nevertheless, results are usually good even if it is not strictly satisfied.
}
\usage{
heterEigenvectorCentrality.multiplex(obj,
                                    indexNode = 1:length(nodes.multiplex(obj)),
                                    W = matrix(1, length(layers.multiplex(obj)),
                                    length(layers.multiplex(obj))),
                                    rowStand = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
}
  \item{indexNode}{
  A vector of IDs (or labels) for the selected nodes on which to calculate the \emph{multiplex heterogeneous eigenvector centrality} coefficients.
}
  \item{W}{
  The \emph{influence matrix} of the multiplex network. For further information, see References.
}
  \item{rowStand}{
  Default is \code{TRUE}. If \code{FALSE}, the vectors given in output are not standardized (i.e. the sum of their components is != 1) and comparison with other (intra)layers or multiplex networks is more difficult.
}

}

\value{
A \code{list} with the L vectors of the multiplex heterogeneous eigenvector centrality indices of the nodes (eventually selected with \code{indexNode} argument), where L is the number of layers of the multiplex network.
}
\references{
  Sola' et al. (2013) \emph{Eigenvector centrality of nodes in multiplex networks.} Chaos Volume 23, Issue 3.
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

# Calculating the multiplex heterogeneous eigenvector centrality indices for the multiplex network:

heterEigenvectorCentrality.multiplex(mplexObj)

# It could be useful to select a set of nodes on which to calculate the index. This can be done
# using the 'indexNode' argument, as it follows:

heterEigenvectorCentrality.multiplex(mplexObj,
                                      indexNode = sample(1:length(nodes.multiplex(mplexObj)), 10)
                                      )

# The particularity of this index is strictly linked to the possibility to include a so-called
# influence matrix in the argument 'W'. This matrix rapresents the weights of the links (i.e. the
# interlayer edges) of the multiplex network. If we set, as an example, a random influence
# matrix, we see that the values of the index changes:

W <- matrix(rbinom(25, 5, .5), 5, 5)
diag(W) <- 0
heterEigenvectorCentrality.multiplex(mplexObj, W = W)

# Another way to visualize the results is to consider the standardized measures. In this case,
# comparisons between indices on different layers can be done, because the sum of the indices
# for each layer are forced to be 1:

heterEigenvectorCentrality.multiplex(mplexObj, rowStand = TRUE)
apply(heterEigenvectorCentrality.multiplex(mplexObj, rowStand = TRUE), 1, sum)
}

\keyword{network}
\keyword{multiplex}
