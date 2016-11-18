# Rscclust

[![Build Status](https://travis-ci.org/fsavje/Rscclust.svg?branch=master)](https://travis-ci.org/fsavje/Rscclust)
[![Build status](https://ci.appveyor.com/api/projects/status/5fvj2nm3ssfmsjab/branch/master?svg=true)](https://ci.appveyor.com/project/fsavje/rscclust/branch/master)
[![codecov](https://codecov.io/gh/fsavje/Rscclust/branch/master/graph/badge.svg)](https://codecov.io/gh/fsavje/Rscclust)

Rscclust is an R wrapper for the [scclust library](https://github.com/fsavje/scclust). The package can be used to construct size constrained clusterings. Subject to user-specified conditions on the minimum size and composition of the clusters, Rscclust derives a partition of a set of data points so that the similarity of points assigned to the same cluster is maximized.

Rscclust is made with large data sets in mind, and it can cluster hundreds of millions of data points within minutes on an ordinary desktop computer. 

The package is currently in alpha, and breaking changes to the API might happen.


## How to install Rscclust

Rscclust is not yet on CRAN, but you can install the current development version using [devtools](https://github.com/hadley/devtools):

```R
devtools::install_github("fsavje/Rscclust")
```

Rscclust contains compiled code. You must, therefore, have a development environment installed. (Use `devtools::has_devel()` to check whether you do.) If no development environment exists, Windows users download and install [Rtools](http://cran.r-project.org/bin/windows/Rtools) and macOS users download and install [Xcode](https://itunes.apple.com/us/app/xcode/id497799835).


## How to use Rscclust

The following snippet shows how Rscclust can be used to make clusterings with both size and type constraints:

```R
library("Rscclust")

set.seed(123456)

# Make example data
my_data <- data.frame(id = 1:100000,
                      type = factor(rbinom(100000, 3, 0.3),
                                    labels = c("A", "B", "C", "D")),
                      x1 = rnorm(100000),
                      x2 = rnorm(100000),
                      x3 = rnorm(100000))

### Construct distance metric
my_dist <- make_distances(my_data,
                          id_variable = "id",
                          dist_variables = c("x1", "x2", "x3"),
                          normalize = "mahalanobize")

### Make clustering with one data point of each type in each cluster
clustering1 <- nng_clustering_types(distance_object = my_dist,
                                    type_labels = my_data$type,
                                    type_size_constraints = c("A" = 1,
                                                              "B" = 1,
                                                              "C" = 1,
                                                              "D" = 1))

### Check that clustering fulfills constraints
check_clustering_types(clustering = clustering1,
                       type_labels = my_data$type,
                       type_size_constraints = c("A" = 1,
                                                 "B" = 1,
                                                 "C" = 1,
                                                 "D" = 1))

### Get statistics about the clustering
get_clustering_stats(clustering = clustering1,
                     distance_object = my_dist)

### Clustering with at least 2 "Bs" in each cluster
clustering2 <- nng_clustering_types(distance_object = my_dist,
                                    type_labels = my_data$type,
                                    type_size_constraints = c("A" = 1,
                                                              "B" = 2,
                                                              "C" = 1,
                                                              "D" = 1))

### At least 8 data points in each cluster in total
clustering3 <- nng_clustering_types(distance_object = my_dist,
                                    type_labels = my_data$type,
                                    type_size_constraints = c("A" = 1,
                                                              "B" = 1,
                                                              "C" = 1,
                                                              "D" = 1),
                                    total_size_constraint = 8)

### Use a radius of 0.5 in the distance metric
clustering4 <- nng_clustering_types(distance_object = my_dist,
                                    type_labels = my_data$type,
                                    type_size_constraints = c("A" = 1,
                                                              "B" = 1,
                                                              "C" = 1,
                                                              "D" = 1),
                                    radius = 0.1)

### Use a different seed finding function
clustering5 <- nng_clustering_types(distance_object = my_dist,
                                    type_labels = my_data$type,
                                    type_size_constraints = c("A" = 1,
                                                              "B" = 1,
                                                              "C" = 1,
                                                              "D" = 1),
                                    seed_method = "lexical")
# All options for `seed_method`:
# "lexical"
# "inwards_order"
# "inwards_updating"
# "inwards_alt_updating"
# "exclusion_order"
# "exclusion_updating"

### Use a different assign function
clustering6 <- nng_clustering_types(distance_object = my_dist,
                                    type_labels = my_data$type,
                                    type_size_constraints = c("A" = 1,
                                                              "B" = 1,
                                                              "C" = 1,
                                                              "D" = 1),
                                    unassigned_method = "closest_assigned")
# All options for `unassigned_method`:
# "ignore"
# "by_nng"
# "closest_assigned"
# "closest_seed"
# "estimated_radius_closest_seed"
```

