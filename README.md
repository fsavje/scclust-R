# scclust for R

[![Build Status](https://travis-ci.org/fsavje/scclust-R.svg?branch=master)](https://travis-ci.org/fsavje/scclust-R)
[![Build status](https://ci.appveyor.com/api/projects/status/27c35hhx7vpigs7k/branch/master?svg=true)](https://ci.appveyor.com/project/fsavje/scclust-r/branch/master)
[![codecov](https://codecov.io/gh/fsavje/scclust-R/branch/master/graph/badge.svg)](https://codecov.io/gh/fsavje/scclust-R)

This package is an R wrapper for the [scclust library](https://github.com/fsavje/scclust). The library provides functions to construct near-optimal size constrained clusterings. Subject to user-specified conditions on the minimum size and composition of the clusters, scclust derives a partition of a set of data points so that the dissimilarity of points assigned to the same cluster is minimized.

scclust is made with large data sets in mind, and it can cluster tens of millions of data points within minutes on an ordinary desktop computer. 


## How to install `scclust`

`scclust` is not yet on CRAN, but you can install the current development version using [devtools](https://github.com/hadley/devtools):

```{r}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("fsavje/scclust-R")
```

The package contains compiled code, and you must have a development environment installed. (Use `devtools::has_devel()` to check whether you do.) If no development environment exists, Windows users download and install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and macOS users download and install [Xcode](https://itunes.apple.com/us/app/xcode/id497799835).


## How to use scclust

The following snippet shows how scclust can be used to make clusterings with both size and type constraints:

```{r}
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
                     dist_variables = c("x1", "x2", "x3"),
                     normalize = "mahalanobize")

# Make clustering with at least 3 data points in each cluster
my_clustering <- sc_clustering(my_dist, 3)

# Check so clustering satisfies constraints
check_scclust(my_clustering, 3)
# > TRUE

# Get statistics about the clustering
get_scclust_stats(my_clustering, my_dist)
# > num_data_points        1.000000e+05
# > ...

# Make clustering with at least one point of each type in each cluster
my_clustering <- sc_clustering(my_dist,
                               type_labels = my_data$type,
                               type_constraints = c("A" = 1, "B" = 1,
                                                    "C" = 1, "D" = 1))

# Check so clustering satisfies constraints
check_scclust(my_clustering,
              type_labels = my_data$type,
              type_constraints = c("A" = 1, "B" = 1,
                                   "C" = 1, "D" = 1))
# > TRUE

# Make clustering with at least 8 points in total of which at least
# one must be "A" and 2 "Bs" and five can be any type
my_clustering <- sc_clustering(my_dist,
                               size_constraint = 8,
                               type_labels = my_data$type,
                               type_constraints = c("A" = 1, "B" = 2))

```
