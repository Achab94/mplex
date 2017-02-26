\name{mplex-package}
\alias{mplex-package}
\alias{mplex}
\docType{package}
\title{
mplex: A package for creating and managing multiplex networks
}
\description{
This package allows implement descriptive tools and function to handle multiplex networks into the R framework. Multiplex networks are composed by a multiplicity of overlapping layers that capture different types of connection between nodes and account for both relationships within a same layer and relationships between a same node located in different layers of the overall network. Some of the most common descriptors, clustering indices and centrality measures can be calculated starting from a set of nodes and layers given in input by the user. Functions are provided to compute some of the most common descriptors, clustering indices and centrality measures starting from a set of nodes and layers given in input by the user. Definitions and notations used throughout the package make reference to De Domenico et al (2014).
}
\details{
\tabular{ll}{
Package: \tab mplex\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2017-02-26\cr
License: \tab GPL (>=2)\cr
}
Functions contained in this package make use of the class \code{multiplex}, created using the function \code{\link{create.multiplex}}.

Some functions to easily recover the elements of the multiplex network (the \code{multiplex}-class object) are provided, such as \code{\link{nodes.multiplex}} for retrieving the list of nodes, \code{\link{adjacency.multiplex}} for the list of all the adjacency matrices of the layers, etc...

Further functions as \code{\link{degree.multiplex}}, \code{\link{meanDegree.multiplex}}, or \code{\link{densityLayers.multiplex}} allow the user to perform some of the basic descriptors for each of the layers composing the multiplex structure of the network.

Functions such as \code{\link{supraAdjacency.multiplex}} or \code{\link{aggregatedOverlapping.multiplex}} perform basic transformations of the network to structures useful for a first-approach visualization of the multiplex network, e.g. the so-called supra-adjacency matrix, or the aggregated overlapping matrix (see Kivela, M. et al (2014)).

Multiplex-specific functions are provided for computing clustering and centrality measures specifically defined for multiplex networks. Examples are \code{\link{degreeCentrality.multiplex}}, \code{\link{supraEigenvectorCentrality.multiplex}} or \code{\link{heterEigenvectorCentrality.multiplex}} and respectively \code{\link{c1Local.multiplex}}, \code{\link{C2Global.multiplex}} or \code{\link{globalOverlayClustering.multiplex}}.
}

\author{
Emanuele Degani

Maintainer: Emanuele Degani <emanuele.achab@gmail.com>
}
\references{

De Domenico et al (2014). \emph{Mathematical formulation of multilayer networks.} Phys. Rev. X 3, 041022.

Kivela, M. et al (2014). \emph{Multilayer Networks.} J. Complex Network. 2(3): 203-271.
}
\seealso{
\code{\link{create.multiplex}}
}

\keyword{package}
\keyword{network}
\keyword{multiplex}
