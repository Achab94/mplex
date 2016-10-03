\name{localClustering.multiplex}
\alias{localClustering.multiplex}
  %- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Layers' local clustering index of the multiplex network
}
\description{
  The function calculates the local clustering index (sometimes called \emph{transitivity measure}, see References) of the nodes (eventually selected with \code{indexNode} argument) for each (intra)layer of the multiplex network. It combines the \code{graph.multiplex} and \code{igraph}'s \code{transitivity} function, so \code{igraph} package is mandatory.
}
\usage{
localClustering.multiplex(obj, indexNode = 1:length(nodes.multiplex(obj)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
}
  \item{indexNode}{
  A vector of IDs (or labels) for the selected nodes on which to calculate the local clustering coefficient.
}
}
\details{
  Local clustering index calculated layer by layer is a good measure for studying the transitivity of a certain node in a specific (intra)layer, but it doesn't take into consideration the triadic structures that take form considering the intralayer relationships (those between two same nodes belonging to different layers). In this case, \code{\link{c1Local.multiplex}} and \code{\link{c2Local.multiplex}} functions are indicated for a \emph{multiplex} analysis of the network.
}
\value{
  A \code{list} of numerical \code{vectors}, each with the local clustering coefficients of the nodes (eventually selecte with \code{indexNode} argument) on a certain (intra)layer of the multiplex network.
}
\references{
  Wassermann and Faust (1994), \emph{Social network analysis.} http://dx.doi.org/10.1017/CBO9780511815478
}
\author{
  Emanuele Degani \email{emanuele.achab@gmail.com}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{create.multiplex}}, \code{\link[igraph]{transitivity}}
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

# Calculating the local clustering index for the nodes of the multiplex network:

localClustering.multiplex(mplexObj)

# Nodes can be eventually selected using the 'indexNode' argument:

localClustering.multiplex(mplexObj,
                          indexNode = sample(1:length(nodes.multiplex(mplexObj)) , 10)
                          )
}
