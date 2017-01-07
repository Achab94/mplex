\name{supraAdjacency.multiplex}
\alias{supraAdjacency.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Supra-adjacency matrix of the multiplex network
}
\description{
  This function creates the so-called \emph{supra-adjacency matrix} (see References), obtained contracting the 4th order tensor (rapresenting De Domenico's et al. mathematical formulation in References) into a 2nd order N*L x N*L tensor (\code{matrix}), given N the number of nodes and L the number of (intra)layers of the network. The supra-adjacency matrix is a useful method for 'visualizing' all the multiplex network into a 2-dimensional workframe (a \code{matrix}).
}
\usage{
  supraAdjacency.multiplex(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
    An object of class \code{multiplex}.
  }
}
\details{
  The \emph{supra-adjacency matrix} visualize all the interlayer and intralayer adjacency matrices of the multiplex network into a single large matrix. It is 'filled' with all the intralayer adjacency matrices on its diagonal, and with all the interlayer adjacency matrices elsewhere.
}
\value{
  A N*L x N*L \code{matrix}, composed as described in 'Details', where N is the number of nodes and L the number of layers of the network.
}
\references{
De Domenico et al (2014). \emph{Mathematical formulation of multilayer networks.} Phys. Rev. X 3, 041022.

Kivela, M. et al (2014). \emph{Multilayer Networks.} J. Complex Network. 2(3): 203-271.
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

# Creating the supra-adjacency matrix of the multiplex network:

supraAdjacency.multiplex(mplexObj)
}
