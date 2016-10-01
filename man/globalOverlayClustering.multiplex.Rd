\name{globalOverlayClustering.multiplex}
\alias{globalOverlayClustering.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Global overlay clustering coefficient of the multiplex network
}
\description{
 This function calculates De Domenico et al.'s \emph{global overlay clustering coefficient} (see References for further information). It is defined as a kind of generalization of \emph{C1 global clustering coefficient} and \emph{C2 global clustering coefficient}, and it's simply obtained as the local clustering coefficient of the \emph{overlay monoplex network} of the multiplex network obtained with \code{aggregated.overlapping.multiplex} function. The index has a standardization that allows comparisons with other multiplex network: the standardizing factor can be obtained with the \code{verbose} argument.
}
\usage{
globalOverlayClustering.multiplex(obj, indexLayer = 1:length(layers.multiplex(obj)), verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
}
  \item{indexLayer}{
  A vector of IDs (or labels) for the selected (intra)layers on which to calculate the global overlay clustering coefficient.
}
  \item{verbose}{
  Default is \code{FALSE}. If \code{TRUE}, it prints on console the value of the standardizing coefficient of the index.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Output is a number, the global overlay clustering coefficient.
}
\references{
  De Domenico et al (2014). \emph{Mathematical formulation of multilayer networks.} http://doi.org/10.1103/PhysRevX.3.041022
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

# Calculating the global overlay clustering coefficient of the multiplex network:

globalOverlayClustering.multiplex(mplexObj)

# Selecting the layers on which to calculate the coefficient, using 'indexLayer' argument:

globalOverlayClustering.multiplex(mplexObj, indexLayer = c(1, 2))

# It could be useful to visualize the standardizing coefficient applied to the denominator
# of the index, mostly if the relationships are weighted and not binary like in this example. 
# Using the 'verbose = TRUE' argument, it is printed on console:

globalOverlayClustering.multiplex(mplexObj, verbose = TRUE)
}