mplex: an R package for multiplex networks
--------

mplex is a package that implement descriptive tools and function to handle multiplex networks into the R framework. Multiplex networks are composed by a multiplicity of overlapping layers that capture different types of connection between nodes and account for both relationships within a same layer and relationships between a same node located in different layers of the overall network. Some of the most common descriptors, clustering indices and centrality measures can be calculated starting from a set of nodes and layers given in input by the user.

Installation
------------

``` r
# The easiest way to get mplex is to install directly from GitHub:

# install.packages("devtools")
devtools::install_github("achab94/mplex")
```

Creating a `multiplex` structure:
------------
``` r
data(aarhus_mplex)

mplexObj <- create.multiplex(nodes = aarhus_mplex$nodes,
                            layersNames = aarhus_mplex$layerNames,
                            layer1 = aarhus_mplex$L1,
                            type1 = "undirected",
                            aarhus_mplex$L2,
                            aarhus_mplex$L3,
                            aarhus_mplex$L4,
                            aarhus_mplex$L5
                            )

class(aarhus_mplex)
```

Contacts:
------------
**Emanuele Degani**: emanuele.achab [at] gmail [dot] com
