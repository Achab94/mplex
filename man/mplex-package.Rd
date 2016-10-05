\name{mplex-package}
\alias{mplex-package}
\alias{mplex}
\docType{package}
\title{
mplex: A package for creating and managing multiplex networks
}
\description{
This package allows to implement multiplex networks into the R framework. Functions are provided to compute some of the most common descriptors, clustering indices and centrality measures starting from a set of nodes and layers given in input by the user. Definitions and notations used throughout the package make reference to De Domenico et al (2014).
}
\details{
\tabular{ll}{
Package: \tab mplex\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-09-27\cr
License: \tab GPL (>=2)\cr
}

Functions contained in this package make use of the \code{multiplex} class, which can be created using the \code{\link{create.multiplex}} function.

Some functions for easily recovering the elements of the multiplex network (the \code{multiplex}-class object) are provided, such as \code{\link{nodes.multiplex}} for retrieving the list of nodes, \code{\link{adjacency.multiplex}} for the list of all the adjacency matrices of the layers, etc...

Other functions such as \code{\link{degree.multiplex}}, \code{\link{meanDegree.multiplex}} or \code{\link{density.multiplex}} allows the user to perform some of the basic descriptors for each of the layers composing the multiplex structure of the network, specifying the nodes and layers on which to compute them.

Functions such as \code{\link{supraAdjacency.multiplex}} or \code{\link{aggregatedOverlapping.multiplex}} performs some basic transformation of the network into particular structures which can be sometimes useful for having a first-approach visualization of the multiplex network, such as the so-called Supra-adjacency matrix, or the Aggregated Overlapping matrix.

Multiplex-specific functions are provided for computing some of the most widespread clustering and centrality measures strictly defined for multiplex networks. In this sense, some centrality measures are calculated using \code{\link{degreeCentrality.multiplex}}, \code{\link{supraEigenvectorCentrality.multiplex}} or \code{\link{heterEigenvectorCentrality.multiplex}}; clustering indices can be computed instead using functions such as \code{\link{c1Local.multiplex}}, \code{\link{C2Global.multiplex}} or \code{\link{globalOverlayClustering.multiplex}}.
}

\author{
Emanuele Degani

Maintainer: Emanuele Degani <emanuele.achab@gmail.com>
}
\references{
De Domenico et al (2014). \emph{Mathematical formulation of multilayer networks.} http://doi.org/10.1103/PhysRevX.3.041022

Kivela, M. et al (2014). \emph{Multilayer Networks.} http://doi.org/10.1093/comnet/cnu016
}

\keyword{ package }
\seealso{
\code{\link{create.multiplex}}
}

