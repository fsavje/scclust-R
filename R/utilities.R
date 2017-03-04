# ==============================================================================
# scclust for R -- R wrapper for the scclust library
# https://github.com/fsavje/scclust-R
#
# Copyright (C) 2016  Fredrik Savje -- http://fredriksavje.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/
# ==============================================================================


#' Check clustering constraints
#'
#' \code{check_scclust} checks whether the clusters in a clustering satisfies a set
#' of conditions on their size and composition.
#'
#' @param clustering a \code{\link{scclust}} object containing a non-empty clustering.
#' @param size_constraint an integer with the required minimum cluster size. If \code{NULL},
#'                              only the type constraints will be checked.
#' @param type_labels a factor or integer containing the type of each data point. May be \code{NULL} when
#'                    \code{type_constraints} is NULL.
#' @param type_constraints a named vector containing type-specific size constraints. If \code{NULL},
#'                              only the overall constraint will be checked.
#'
#' @return Returns \code{TRUE} if \code{clustering} satisfies the clustering constraints, and
#'         \code{FALSE} if it does not. Throws an error if \code{clustering} is an invalid
#'         instance of the \code{\link{scclust}} class.
#'
#' @seealso See \code{\link{sc_clustering}} for details on
#'          how to specific type labels and constraints.
#
#' @examples
#' # Each cluster contains at least two data
#' # points; it's valid for `size_constraint == 2`.
#' my_clust_obj1 <- scclust(c("A", "A", "B", "C",
#'                                   "B", "C", "C", "A",
#'                                   "B", "B"))
#' check_scclust(my_clust_obj1, 2)
#' # > TRUE
#'
#'
#' # One cluster contains only one point. This is
#' # an invalid clustering for `size_constraint == 2`.
#' my_clust_obj2 <- scclust(c("A", "A", "B", "C",
#'                                   "B", "C", "C", "A",
#'                                   "B", "B", "D"))
#' check_scclust(my_clust_obj2, 2)
#' # > FALSE
#'
#'
#' # `my_clust_obj1` doesn't contain clusters with
#' # at least 5 units, it's invalid for `size_constraint == 2`.
#' check_scclust(my_clust_obj1, 5)
#' # > FALSE
#'
#'
#' # Clustering
#' my_clust_obj3 <- scclust(c(1, 1, 2, 3, 2, 3, 3, 1, 2, 2))
#'
#' # Data point types
#' my_types <- factor(c("x", "y", "y", "z", "z", "x", "y", "z", "x", "x"))
#'
#'
#' # Check whether each cluster contains one point of each
#' # of type "x", "y" and "z".
#' check_scclust(my_clust_obj3,
#'                  NULL,
#'                  my_types,
#'                  c("x" = 1, "y" = 1, "z" = 1))
#' # > TRUE
#'
#'
#' # Check whether each cluster contains one point of each
#' # of type "x", "y" and "z" and at least three points in total.
#' check_scclust(my_clust_obj3,
#'                  3,
#'                  my_types,
#'                  c("x" = 1, "y" = 1, "z" = 1))
#' # > TRUE
#'
#'
#' # Check whether each cluster contains five data points of type "y".
#' check_scclust(my_clust_obj3,
#'                  NULL,
#'                  my_types,
#'                  c("y" = 5))
#' # > FALSE
#'
#'
#' # Check whether each cluster contains one data point of
#' # both "x" and "z" and at least three points in total.
#' check_scclust(my_clust_obj3,
#'                  3,
#'                  my_types,
#'                  c("x" = 1, "z" = 1))
#' # > TRUE
#'
#'
#' # Using integers as type labels.
#' my_int_types <- c(1L, 2L, 2L, 3L, 3L, 1L, 2L, 3L, 1L, 1L)
#' check_scclust(my_clust_obj3,
#'                  NULL,
#'                  my_int_types,
#'                  c("1" = 1, "2" = 1, "3" = 1))
#' # > TRUE
#'
#' @export
check_scclust <- function(clustering,
                          size_constraint = NULL,
                          type_labels = NULL,
                          type_constraints = NULL) {
  ensure_scclust(clustering)
  num_data_points <- length(clustering)
  if (is.null(type_constraints)) {
    type_labels <- NULL
    size_constraint <- coerce_size_constraint(size_constraint, num_data_points)
  } else {
    type_labels <- coerce_type_labels(type_labels, num_data_points)
    type_constraints <- coerce_type_constraints(type_constraints)
    type_constraints <- make_type_size_constraints(type_constraints,
                                                   type_labels)
    size_constraint <- coerce_total_size_constraint(size_constraint,
                                                    type_constraints,
                                                    num_data_points)
  }

  .Call(Rscc_check_scclust,
        clustering,
        size_constraint,
        unclass(type_labels),
        type_constraints)
}


#' Derive statistics about a clustering.
#'
#' \code{get_scclust_stats} calculates important statistics about a clustering, such as average
#' cluster size and average within-cluster distance.
#'
#' The function reports the following statistics:
#'
#' \tabular{ll}{
#'   \code{num_data_points} \tab The total number of data points. \cr
#'   \code{num_assigned} \tab The number of points assigned to a cluster. \cr
#'   \code{num_clusters} \tab The total number of clusters. \cr
#'   \code{num_populated_clusters} \tab The number of non-empty clusters. \cr
#'   \code{min_cluster_size} \tab The size of the smallest cluster. \cr
#'   \code{max_cluster_size} \tab The size of the largest cluster. \cr
#'   \code{avg_cluster_size} \tab The average cluster size. \cr
#'   \code{sum_dists} \tab The sum of all within-cluster distances. \cr
#'   \code{min_dist} \tab The overall smallest within-cluster distance. \cr
#'   \code{max_dist} \tab The overall largest within-cluster distance. \cr
#'   \code{cl_avg_min_dist} \tab The average of the clusters' smallest distance. \cr
#'   \code{cl_avg_max_dist} \tab The average of the clusters' largest distance. \cr
#'   \code{cl_avg_dist_weighted} \tab The average of the clusters' average distance weighed by cluster size. \cr
#'   \code{cl_avg_dist_unweighted} \tab The average of the clusters' average distance (unweighed). \cr
#' }
#'
#' Let \eqn{d(i,j)}{d(i,j)} denote the distance between data points \eqn{i}{i} and \eqn{j}{j}. Let \eqn{c}{c} be a cluster
#' containing all indices of data points assigned to the cluster. Let
#' \deqn{D(c) = \{d(i,j): i,j \in c \wedge i>j\}}{D(c) = { d(i,j) : i,j in c and i > j }}
#' be a set function returning all unique within-cluster distances of \eqn{c}{c}. Finally, let \eqn{C}{C} be a set containing
#' all clusters.
#'
#' \code{sum_dists} is defined as:
#' \deqn{\sum_{c\in C} sum(D(c))}{\sum_[c in C] sum(D(c))}
#'
#' \code{min_dist} is defined as:
#' \deqn{\min_{c\in C} \min(D(c))}{min_[c in C] min(D(c))}
#'
#' \code{max_dist} is defined as:
#' \deqn{\max_{c\in C} \max(D(c))}{max_[c in C] max(D(c))}
#'
#' \code{cl_avg_min_dist} is defined as:
#' \deqn{\sum_{c\in C} \frac{\min(D(c))}{|C|}}{\sum_[c in C] min(D(c)) / count(C)}
#'
#' \code{cl_avg_max_dist} is defined as:
#' \deqn{\sum_{c\in C} \frac{\max(D(c))}{|C|}}{\sum_[c in C] max(D(c)) / count(C)}
#'
#' Let:
#' \deqn{AD(c) = \frac{sum(D(c))}{|D(c)|}}{AD(c) = sum(D(c)) / count(D(c))}
#' be the average within-cluster distance in cluster \eqn{c}{c}.
#'
#' \code{cl_avg_dist_weighted} is defined as:
#' \deqn{\sum_{c\in C} \frac{|c| AD(c)}{num_assigned}}{\sum_[c in C] count(c) * AD(c) / num_assigned}
#' where \eqn{num_assigned}{num_assigned} is the number of assigned data points (see above).
#'
#' \code{cl_avg_dist_unweighted} is defined as:
#' \deqn{\sum_{c\in C} \frac{AD(c)}{|C|}}{\sum_[c in C] AD(c) / count(C)}
#'
#' @param clustering a \code{scclust} object containing an existing non-empty clustering.
#' @param distances a distance object as produced by \code{\link[distances]{distances}}.
#'
#' @return Returns a list of class "scclust_stats" containing important statistics
#'         about the inputted clustering.
#'
#' @examples
#' library(distances)
#'
#' my_clust_obj <- scclust(c("A", "A", "B", "C", "B",
#'                                  "C", "C", "A", "B", "B"))
#' my_data_points <- data.frame(x = c(0.1, 0.2, 0.3, 0.4, 0.5,
#'                                    0.6, 0.7, 0.8, 0.9, 1.0),
#'                              y = c(10, 9, 8, 7, 6,
#'                                    10, 9, 8, 7, 6))
#' my_distance_obj <- distances(my_data_points)
#'
#' get_scclust_stats(my_clust_obj, my_distance_obj)
#'
#' # >                        Value
#' # > num_data_points        10.0000000
#' # > num_assigned           10.0000000
#' # > num_clusters            3.0000000
#' # > num_populated_clusters  3.0000000
#' # > min_cluster_size        3.0000000
#' # > max_cluster_size        4.0000000
#' # > avg_cluster_size        3.3333333
#' # > sum_dists              18.2013097
#' # > min_dist                0.5000000
#' # > max_dist                3.0066593
#' # > cl_avg_min_dist         0.8366584
#' # > cl_avg_max_dist         2.4148611
#' # > cl_avg_dist_weighted    1.5575594
#' # > cl_avg_dist_unweighted  1.5847484
#'
#' @export
get_scclust_stats <- function(clustering,
                              distances) {
  ensure_scclust(clustering)
  num_data_points <- length(clustering)
  ensure_distances(distances, num_data_points)

  clust_stats <- .Call(Rscc_get_scclust_stats,
                       clustering,
                       distances)
  structure(clust_stats,
            class = c("scclust_stats"))
}


#' @export
print.scclust_stats <- function(x, ...) {
  tmp_table <- as.table(format(as.matrix(unlist(x))))
  colnames(tmp_table) <- "Value"
  print(tmp_table)
}
