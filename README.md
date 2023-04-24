# scclust for R

[![CRAN Status](https://www.r-pkg.org/badges/version/scclust)](https://cran.r-project.org/package=scclust)

The `scclust` package is an R wrapper for the [scclust library](https://github.com/fsavje/scclust). The package provides functions to construct near-optimal size-constrained clusterings.

Most conventional clustering functions restrict the number of clusters, but do not impose restrictions on the content of the clusters (see, for example, [k-means](https://en.wikipedia.org/wiki/K-means_clustering)). scclust takes another route. It imposes conditions on the content of the clusters, but allow any number of them to be formed. Specifically, subject to user-specified constraints on the size and composition of the clusters, scclust constructs a clustering so that within-cluster pair-wise distances are minimized.

It is possible to impose an overall size constraint so that each cluster must contain at least a certain number of points in total. It is also possible to impose constraints on the composition of the clusters so that each cluster must contain a certain number of points of different types. For example, in a sample with "red" and "blue" data points, one can constrain the clustering so that each cluster must contain at least 10 points in total of which at least 3 must be "red" and at least 2 must be "blue".

scclust was made with large data sets in mind, and it can cluster tens of millions of data points within minutes on an ordinary desktop computer. 


## How to install

`scclust` is on CRAN and can be installed by running:

```R
install.packages("scclust")
```


## How to install development version

It is recommended to use the stable CRAN version, but the latest development version can be installed directly from Github using [devtools](https://github.com/hadley/devtools):

```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("fsavje/scclust-R")
```

The package contains compiled code, and you must have a development environment to install the development version. (Use `devtools::has_devel()` to check whether you do.) If no development environment exists, Windows users download and install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and macOS users download and install [Xcode](https://itunes.apple.com/us/app/xcode/id497799835).


## How to use scclust

The following snippet shows how scclust can be used to make clusters with both size and type constraints. See the package documentation for more details. 

```R
# Make example data
my_data <- data.frame(id = 1:100000,
                      type = factor(rbinom(100000, 3, 0.3),
                                    labels = c("A", "B", "C", "D")),
                      x1 = rnorm(100000),
                      x2 = rnorm(100000),
                      x3 = rnorm(100000))

# Construct distance metric
my_dist <- distances(my_data,
                     id_variable = "id",
                     dist_variables = c("x1", "x2", "x3"))

# Make clustering with at least 3 data points in each cluster
my_clustering <- sc_clustering(my_dist, 3)

# Check so clustering satisfies constraints
check_clustering(my_clustering, 3)
# > TRUE

# Get statistics about the clustering
get_clustering_stats(my_dist, my_clustering)
# > num_data_points        1.000000e+05
# > ...

# Make clustering with at least one point of each type in each cluster
my_clustering <- sc_clustering(my_dist,
                               type_labels = my_data$type,
                               type_constraints = c("A" = 1, "B" = 1,
                                                    "C" = 1, "D" = 1))

# Check so clustering satisfies constraints
check_clustering(my_clustering,
                 type_labels = my_data$type,
                 type_constraints = c("A" = 1, "B" = 1,
                                      "C" = 1, "D" = 1))
# > TRUE

# Make clustering with at least 8 points in total of which at least
# one must be "A", two must be "B" and five can be any type
my_clustering <- sc_clustering(my_dist,
                               size_constraint = 8,
                               type_labels = my_data$type,
                               type_constraints = c("A" = 1, "B" = 2))
```
