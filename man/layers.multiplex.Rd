\name{layers.multiplex}
\alias{layers.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extraction of the layers' dataframes
}
\description{
  The function allows a user-friendly method for extracting the information about selected layers of the multiplex network.
}
\usage{
layers.multiplex(obj, index = 1:length(obj$layers), label = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
  }
  \item{index}{
  A vector of the layers' IDs selected for the extraction. By default, the vector contains the IDs of all the layers of the network.
  }
  \item{label}{
  \code{FALSE} by default. If \code{TRUE}, the function extracts only the labels of the selected layers.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A \code{list} with the layers' \code{data.frame} structures selected by the \code{index} argument.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Emanuele Degani \email{emanuele.achab@gmail.com}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

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

# Extracting the edges of the layers (relationships):

layers.multiplex(mplexObj)

# Edges can be extracted even for some specific layers (for example, the
# second and third layers), selected with the 'index' argument:

layers.multiplex(mplexObj, index = c(2, 3))

# or, easily, referring to their names:

layers.multiplex(mplexObj, index = c("facebook", "coauthor"))

# The 'label' argument, if TRUE, returns the names/labels of the selected layers:

layers.multiplex(mplexObj, index = c(2, 4), label = TRUE)
}
