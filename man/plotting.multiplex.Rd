\name{plotting.multiplex}
\alias{plotting.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plotting layers of a multiplex network
}
\description{
This function allows to have a visual approach to the network, according to the number of layers the user wants to display. The function plots a single, interactive graph is just one layer is selected, using \code{\link[networkD3]{simpleNetwork}} function, or a sequence of graphs if >1 layers are selected, using \code{\link[igraph]{plot.igraph}} method. All graphs are easily exportable using RStudio buttons.
}
\usage{
  plotting.multiplex(obj, index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
  }
  \item{index}{
  A vector of IDs (or labels) for the selected layers on which to display the graph.
  }
}

\value{
  If length(index) == 1, then an interactive graph for the selected layer is displayed in the Viewer pane. If length(index) > 1, a sequence of graphs are displayed in the Plot pane, browsable with the 'Enter' key.
}

\seealso{
\code{\link{create.multiplex}}, \code{\link[igraph]{plot.igraph}}, \code{\link[networkD3]{simpleNetwork}}
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

# Visualize the interactive graph for the "facebook" layer:

plotting.multiplex(mplexObj, 2)

# Visualize the sequence of graphs for the "coauthor" and "leisure" layers:

plotting.multiplex(mplexObj, c(3, 4))
}

\keyword{network}
\keyword{plots}
\keyword{multiplex}
