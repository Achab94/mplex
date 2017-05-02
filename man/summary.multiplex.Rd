\name{summary.multiplex}
\alias{summary.multiplex}
\alias{summary}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Summary of the multiplex network
}
\description{
This function returns the number of layers and nodes of the multiplex network.
}
\usage{
  summary(obj, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
  }
  \item{\dots}{
  Any other useful options.
  }
}

\value{
  A string, namely the number of nodes and layers of the multiplex network
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

# Summarize the multiplex network:

summary(mplexObj)
}

\keyword{network}
\keyword{multiplex}
