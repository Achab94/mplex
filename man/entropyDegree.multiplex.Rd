\name{entropyDegree.multiplex}
\alias{entropyDegree.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Entropy of the degrees' distribution
}
\description{
  This function calculates the entropy measure used to evaluate the homogeneity/heterogeneity of the nodes' distribution on a certain (intra)layer of the multiplex network. The formula is inspired by Battiston et al. (2014) work (see References).
}
\usage{
entropyDegree.multiplex(obj,
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
  A vector of IDs (or labels) for the selected nodes on which to calculate the entropy measure. Default selects all the nodes of the network.
}
  \item{indexOverlappingLayer}{
  A vector of IDs (or labels) for the selected (intra)layers on which to calculate the \emph{aggregated overlapping matrix} with \code{\link{aggregatedOverlapping.multiplex}} function, used inside entropy's formula. Default selects all the layers of the network.
}
  \item{display}{
  Default is \code{FALSE}. If \code{TRUE}, the function plots a user-friendly visualization of the entropy measures' distribution that helps doing comparisons with others multiplex networks.
}
}
\details{
  Instead of the formula written in Battiston et al. (2014) work, the function has a small correction in the denominator, avoiding \code{NaN} output if the degree is zero.
}
\value{
  A numeric \code{vector} with all the entropy measures evaluated on the nodes of the multiplex network (eventually selected with the \code{indexNode} argument).
}
\references{
Battiston et al. (2014) \emph{Structural measures for multiplex networks.} http://dx.doi.org/10.1103/PhysRevE.89.032804
}
\author{
  Emanuele Degani \email{emanuele.achab@gmail.com}
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


# Calculating the entropy of the degrees' distribution:

entropyDegree.multiplex(mplexObj)

# A better way to visualize the results is given by the 'display = TRUE' option.
# We see, in this case, that quite all the nodes have an entropy measure > 0.80: this mean
# that they have a similar degree in all the intralayer of the multiplex network; nonetheless
# 3 employees have a really low entropy measure:

entropyDegree.multiplex(mplexObj, display = TRUE)

# It could be also useful to select just some nodes or levels on which to calculate the
# entropy measure, in order to investigate which layer(s) is (are) causing a low measure
# of entropy for a certain node.
# This can be easily done with both 'indexNode' and 'indexOverlappingLayer' arguments:

entropyDegree.multiplex(mplexObj,
                        indexNode = sample(1:length(nodes.multiplex(mplexObj)), 10),
                        indexOverlappingLayer = (1:5)[-3],
                        display = TRUE
                        )
}
