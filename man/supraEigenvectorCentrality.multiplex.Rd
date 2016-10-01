\name{supraEigenvectorCentrality.multiplex}
\alias{supraEigenvectorCentrality.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Multiplex eigenvector centrality indices of the multiplex network
}
\description{
  This function calculates the \emph{multiplex eigenvector centrality} indices of the multiplex network, for each of its (inter)layers. The operation is conducted calculating the eigenvector referred to the maximum (in absolute value) eigenvalue of the \emph{supra adjacency matrix} of the multiplex network, obtained with the \code{supra.adjacency.multiplex} function. Defined N the number of nodes and L the number of (inter)layers of the multiplex network, the supra-adjacency matrix is N*L x N*L, thus the eigenvector has length N*L. The L vectors given in output are simply obtained breaking the N*L eigenvector into L vectors, each of length N and referred to the multiplex eigenvector centrality measures of the nodes (eventually selected with \code{indexNode} argument) on a certain (intra)layer of the multiplex network.
}
\usage{
supraEigenvectorCentrality.multiplex(obj, indexNode = 1:length(nodes.multiplex(obj)), rowStand = TRUE, testIrreducibility = FALSE, maxPower = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
}
  \item{indexNode}{
  A vector of IDs (or labels) for the selected nodes on which to calculate the \emph{multiplex eigenvector centrality} coefficients.
}
  \item{rowStand}{
  Default is \code{TRUE}. If \code{FALSE}, the vectors given in output are not standardized (i.e. the sum of their components is != 1) and comparison with other (intra)layers or multiplex networks is more difficult.
}
  \item{testIrreducibility}{
  Default is \code{FALSE}. If \code{TRUE}, the function tests the irreducibility of the supra adjacency matrix and prints on console the result. Irreducibility is a needed hypotesy to satisfy the Perron-Frobenius theorem, which ensures the positivity of the eigenvector assosiacted to the maximum eigenvalue of the supra adjacency matrix of the multiplex network.
}
  \item{maxPower}{
  Default is \code{100}. It is the power on which the irreducibility test stops.
}
}
\details{
  Irreducibility tested with the \code{testIrreducibility} is a requirred hypotesy of the Perron-Frobenius theorem, but results are usually good even if it is not strictly satisfied.
}
\value{
  Output is a \code{list} with the L vectors of the multiplex eigenvector centrality indices of the nodes (eventually selected with \code{indexNode} argument).
}
\references{
  De Domenico et al (2014). \emph{Mathematical formulation of multilayer networks.} http://doi.org/10.1103/PhysRevX.3.041022
}
\author{
  Emanuele Degani \email{emanuele.achab@gmail.com}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{create.multiplex}}, \code{\link{supraAdjacency.multiplex}}
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

# Calculating the multiplex eigenvector centrality indices for the multiplex network.
# Sometimes, a round( , 5) could be useful to better visualize the results:

supraEigenvectorCentrality.multiplex(mplexObj)
round(supraEigenvectorCentrality.multiplex(mplexObj), 5)

# It can also be possible to select the IDs of the nodes on which to calculate the index,
# using 'indexNode' argument as in this case:

round(
supraEigenvectorCentrality.multiplex(mplexObj, 
                            indexNode = sample(1:length(nodes.multiplex(mplexObj)), 10))
     , 5)
     
# Another way to visualize the results is to consider the standardized measures. In this case,
# comparisons between indices on different layers can be done, because the sum of the indices
# for each layer are forced to be 1:

supraEigenvectorCentrality.multiplex(mplexObj, rowStand = TRUE)
apply(supraEigenvectorCentrality.multiplex(mplexObj, rowStand = TRUE), 1, sum)

# One last feature of this function, is the possibility to test the irreducibility condition
# (needed for satisfying the Perron-Frobenius Theorem) with the 'testIrreducibility = TRUE'
# option, and to fix the maximum number of iterations to do before stopping the test with
# 'maxPower' argument. If 'testIrreducibility = TRUE', then it is printed on console whether
# the irreducibility condition is satisfied or not on the supra-adjacency matrix:

supraEigenvectorCentrality.multiplex(mplexObj, testIrreducibility = TRUE, maxPower = 100)
}