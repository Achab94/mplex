\name{meanDegree.multiplex}
\alias{meanDegree.multiplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Layers' mean degree of the multiplex network
}
\description{
  This function calculates the mean of the degrees on each (intra)layer of the multiplex network. For (intra)layers with \emph{directed} relationships, the function considers the 'total degree' descriptor.
}
\usage{
meanDegree.multiplex(obj, indexNodeMean = 1:length(nodes.multiplex(obj)), verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  An object of class \code{multiplex}.
}
  \item{indexNodeMean}{
  A vector of IDs (or labels) for the selected nodes on which to calculate the mean degree. Default selects all the nodes of the network.
}
  \item{verbose}{
  Default is \code{FALSE}. If \code{TRUE}, the function prints on console the name of the selected nodes on which the means are evaluated.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The output is a \code{list} containing the means of the degrees for each (intra)layer of the multiplex network.
}
\references{
  Wassermann and Faust (1994), \emph{Social network analysis.} http://dx.doi.org/10.1017/CBO9780511815478 
}
\author{
  Emanuele Degani \email{emanuele.achab@gmail.com}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{create.multiplex}}, \code{\link{degree.multiplex}}
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

# Calculating the mean degrees on the 5 intralayers of the multiplex network.
# We see, for example, that on average each employee is connected to other 6 collegues
# through "lunch" and "work" relationship, and 0-1 through the "coauthor" relationship.

meanDegree.multiplex(mplexObj)

# Calculating the mean degrees selecting 10 employees randomly on which to execute the calculation:

meanDegree.multiplex(mplexObj, indexNodeMean = sample(1:length(nodes.multiplex(mplexObj)), 10))

# Using the 'verbose = TRUE' option, the function prints the name of the nodes on which the mean
# degree is performed:

meanDegree.multiplex(mplexObj, 
                    indexNodeMean = sample(1:length(nodes.multiplex(mplexObj)), 10), 
                    verbose = TRUE
                    )
}