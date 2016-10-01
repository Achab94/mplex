\name{degreeDistribution.multiplex}
\alias{degreeDistribution.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Layers' degrees' distribution of the multiplex network
}
\description{
  This function calculates the degrees' distribution of the nodes on each (intra)layer of the multiplex network, combining the \code{graph.multiplex} and \code{igraph}'s \code{degree.distribution} functions. For (intra)layers with \emph{directed} relationships, the function considers 'total degree' descriptors.
Degrees' distribution shows, for each (intra)layer of the network, the maximum degree on each layer and the percentage of a certain degree in the total distribution of the (intra)layer.
}
\usage{
  degreeDistribution.multiplex(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
  }
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A list with, for each (intra)layer, the degrees' distribution.
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
  \code{\link{create.multiplex}}, \code{\link[igraph]{degree.distribution}}
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

# Calculating the degrees' distribution of the 61 nodes of the multiplex network, for each
# layer. We see, for example, that in the first intralayer (the "lunch" relationship), 
# 0.0327 (3.27\%) of the employees (2 =~ 61 * 0.0327) are linked to other 15 employees 
# of the CS Aarhus Department.

degreeDistribution.multiplex(mplexObj)
}