\name{type.multiplex}
\alias{type.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extraction of the relationship types' vector of the selected layers
}
\description{
  The function allows a user-friendly method for extracting the vector of relationship types for the selected layers of the multiplex network.
}
\usage{
  type.multiplex(obj, index = 1:length(obj$type))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
}
  \item{index}{
  A vector of the layers' IDs selected for the extraction of their relative relationship type string. By default, the vector contains the IDs of all the layers of the network.
}
}

\value{
  A \code{vector} with the layers' relationship type \code{strings} selected by the \code{index} argument.
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

# Extracting the type of relationships of the layers:

type.multiplex(mplexObj)

# or for some of them:

type.multiplex(mplexObj, index = c("leisure", "facebook"))
}

\keyword{network}
\keyword{multiplex}
