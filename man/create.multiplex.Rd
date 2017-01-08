\name{create.multiplex}
\alias{create.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Creation of the \code{multiplex} structure
}
\description{
  This function creates the \code{multiplex} structure, which is used in all the other functions of the package. This function implements the tensor-based framework proposed by De Domenico (2014) creating a list with nodes, layers and layers' names, adjacency matrices and types of relationships within each layer.
}
\usage{
  create.multiplex(nodes, layersNames = FALSE, layer1, type1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nodes}{
  Any \code{data.frame} structure containing information about the nodes: for each row, the first column is the ID of the node, the second is the label with its name/description.
}
  \item{layersNames}{
  Allows to specify a vector of labels for the layers given in input by the next arguments. Default is \code{FALSE}: in this case, the i-th layer will be called \emph{Layer(i)}.
}
  \item{layer1}{
  Any \code{data.frame} or \code{matrix} structure containing the information about the first layer of the multiplex network: for each row, the first column is the ID of the starting node, the second is the ID of the ending node, the third is the weight of the relationship (1 if unweighted).
}
  \item{type1}{
  A string describing the type of relationship expected for the first layer (\code{layer1}). It could be \emph{directed}, \emph{undirected} or any abbreviation.
}
  \item{\dots}{
  This argument allows to input any other layer belonging to the multiplex network, in the same structure as \code{layer1} argument. It also allows to specify the types of relationships of the layer (as in \code{type1} argument), just after the layer's structure. If it is not specified, the relationship is assumed to be \emph{undirected}.
}
}

\value{
  A list of the following objects, which will be useful in all the other functions developed in this package:
  \item{nodes}{
  The same \code{nodes} object given in input, containing the information about the nodes of the multiplex structure;
}
  \item{layers}{
  A \code{list} with all the \code{data.frame} structures (or \code{matrix} converted in \code{data.frame}) taken as input, in \code{layer1} or \dots argument;
}
  \item{adjacency}{
  A \code{list} object with all the adjacency matrices, one for each layer of the multiplex network;
}
  \item{type}{
  A \code{vector} object with all the strings describind the types of relationship, one for each layer of the multiplex network;
}
  \item{layersNames}{
  A \code{vector} object containing all the strings for the labels of the layers. If given in input (with \code{layersNames} argument), it is the same vector.
}
  \item{interlayersMatrix}{
  A N*L x N*L \code{matrix} object, where N is the number of nodes and L the number of layers of the multiplex network, rapresenting the interlayer structure as matrix of the adjacency matrices of the interlayer relationships between each couple of intralayers.
  }
}
\references{
De Domenico et al (2014). \emph{Mathematical formulation of multilayer networks.} Phys. Rev. X 3, 041022.
}

\section{Warning }{
  It is suggested to specify \code{layersNames = FALSE} argument even if no vector of layers names is given in input. If not, execution could cause problem/errors.
}

\examples{
# Loading Aarhus CS Department dataset.

data(aarhus_mplex)

# Creating the multiplex object using the dataset loaded into aarhus_mplex object.
# Note that all the layers refers to undirected relationship (so, their 'type'
# arguments are not specified).

mplexObj <- create.multiplex(nodes = aarhus_mplex$nodes,
                            layersNames = aarhus_mplex$layerNames,
                            layer1 = aarhus_mplex$L1,
                            type1 = "undirected",
                            aarhus_mplex$L2,
                            aarhus_mplex$L3,
                            aarhus_mplex$L4,
                            aarhus_mplex$L5
                            )

# Display the structure of the multiplex object created:

str(mplexObj)
}

\keyword{network, multiplex}
