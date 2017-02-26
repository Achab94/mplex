\name{c2Local.multiplex}
\alias{c2Local.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  c2 local clustering coefficients for multiplex networks
}
\description{
  Computes Battiston et al.'s \emph{C2} local clustering coefficient of a multiplex network.
}
\details{
Instead of the classical definition of \emph{transitivity} and \emph{local clustering coefficient}, this index accounts for all the possible tradics structures that takes form even considering interlayers relationships. The index is defined, for a certain node v_i (possibly selected with \code{indexNode} argument) as the ratio between the number of 3-triangles that contains v_i, and (L - 2) times the number of 2-traids centred in v_i, defines L as the number of (intra)layers taken into consideration.
}
\usage{
c2Local.multiplex(obj,
                  indexNode = 1:length(nodes.multiplex(obj)),
                  indexLayer = 1:length(layers.multiplex(obj))
                  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
}
  \item{indexNode}{
  A vector of IDs (or labels) for the selected nodes on which to calculate the C2 local clustering coefficient.
}
  \item{indexLayer}{
  A vector of IDs (or labels) for the selected (intra)layers on which to calculate the C2 local clustering coefficient.
}
}

\section{Warnings }{
  Note that the C2 local clustering coefficient can be defined for a multiplex network with, at least, 3 (intra)layers. So, if \code{indexLayer} argument has length <3, the function produces an error and halts.
}

\details{
Function use \emph{permutations()} function from \emph{gtools} package. More informations available at \url{https://cran.r-project.org/web/packages/gtools/index.html}.
}
\value{
  A numerical \code{vector} with all the C2 local clustering coefficients for the nodes selected in \code{indexNode} argument.
}
\references{
Battiston et al. (2014) \emph{Structural measures for multiplex networks.} Phys. Rev. E 89, 032804.
}


\seealso{
  \code{\link{create.multiplex}}, \code{\link{C2Global.multiplex}}
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

# Computing the so-called c2 clustering index for the nodes of the multiplex network:

c2Local.multiplex(mplexObj)

# It can be useful to select just some nodes or layers on which to calculate the
# c2 clustering index. This can be easily done with 'indexNode' and 'indexLayer' arguments:

c2Local.multiplex(mplexObj,
                  indexNode = sample(1:length(nodes.multiplex(mplexObj)), 10),
                  indexLayer = c(1, 3, 4)
                  )
}

\keyword{network}
\keyword{multiplex}
