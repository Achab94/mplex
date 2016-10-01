\name{c1Local.multiplex}
\alias{c1Local.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  c1 local clustering coefficients for multiplex networks
}
\description{
  Calculates Battiston et al. (2014)'s \emph{C1} local clustering coefficient of a multiplex network. Instead of the classical definition of \emph{transitivity} and \emph{local clustering coefficient}, this index takes into consideration all the possible triadic structures that takes form even considering interlayers relationships. The index is defined, for a certain node v_i (eventually selected with \code{indexNode} argument) as the ratio between the number of 2-triangles that contains v_i, and (L - 1) times the number of 1-traids centred in v_i, defined L as the number of (intra)layers taken into consideration.
}
\usage{
c1Local.multiplex(obj, indexNode = 1:length(nodes.multiplex(obj)), indexLayer = 1:length(layers.multiplex(obj)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
}
  \item{indexNode}{
  A vector of IDs (or labels) for the selected nodes on which to calculate the C1 local clustering coefficient.
}
  \item{indexLayer}{
  A vector of IDs (or labels) for the selected (intra)layers on which to calculate the C1 local clustering coefficient.
}
}
\details{
  
}
\section{Warnings }{
  Note that the C1 local clustering coefficient can be defined for a multiplex network with, at least, 2 (intra)layers. So, if \code{indexLayer} argument has length 1, the function produces a warning in console.
}
\value{
  A numerical \code{vector} with all the C1 local clustering coefficients for the nodes selected in \code{indexNode} argument.
}
\references{
  Battiston et al. (2014) \emph{Structural measures for multiplex networks.} http://dx.doi.org/10.1103/PhysRevE.89.032804
}
\author{
  Emanuele Degani \email{emanuele.achab@gmail.com}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{create.multiplex}}, \code{\link{C1Global.multiplex}}
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

# Calculating the so-called c1 clustering index for the nodes of the multiplex network:

c1Local.multiplex(mplexObj)

# It can be useful to select just some nodes or layers on which to calculate the
# c1 clustering index. This can be easily done with 'indexNode' and 'indexLayer' arguments:

c1Local.multiplex(mplexObj, 
                  indexNode = sample(1:length(nodes.multiplex(mplexObj)), 10), 
                  indexLayer = c(1, 3, 4)
                  )
}