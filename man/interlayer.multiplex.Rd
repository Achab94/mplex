\name{interlayer.multiplex}
\alias{interlayer.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extraction of the interlayer adjacency matrix
}
\description{
  The function allows a user-friendly method for extracting the interlayer adjacency between two different layers in the multiplex network.
}
\usage{
interlayer.multiplex(obj, level1, level2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
}
  \item{level1}{
  The ID (or name) of the first layer.
}
  \item{level2}{
  The ID (or name) of the second layer.
}
}

\value{
  A N*N \code{matrix} with the interlayer (adjacency) matrix between the selected levels, where N is the number of nodes of the multiplex network.
}
\references{
De Domenico et al (2014). \emph{Mathematical formulation of multilayer networks.} Phys. Rev. X 3, 041022.

Kivela, M. et al (2014). \emph{Multilayer Networks.} J. Complex Network. 2(3): 203-271.
}

\section{Warning }{If arguments \code{level1} = \code{level2}, the function can't perform. In fact, the interlayer matrix between the same level is called \emph{intralayer matrix} and it can be extracted with the \code{adjacency.multiplex} function.}


\seealso{
\code{\link{create.multiplex}}, \code{\link{adjacency.multiplex}}
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

# Extracting the interlayer adjacency matrix referred to the interlayer
# relationship between (for example) layer 1 and layer 4.
# Note that all the entries outside the diagonal are 0s (that's consistent
# with the definition of multiplex network) and that not all the entries of the diagonal are 1s.

interlayer.multiplex(mplexObj, level1 = 1, level2 = 4)
diag(interlayer.multiplex(mplexObj, level1 = 1, level2 = 4))

# Example of warning if 'level1' and 'level2' are equals (remove # for running):

# interlayer.multiplex(mplexObj, 3, 3)
}

\keyword{network}
\keyword{multiplex}
