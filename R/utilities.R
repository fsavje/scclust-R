# ==============================================================================
# scclust for R -- R wrapper for the scclust library
# https://github.com/fsavje/scclust-R
#
# Copyright (C) 2016-2017  Fredrik Savje -- http://fredriksavje.com
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
#' @param type_labels a vector containing the type of each data point. May be
#'                    \code{NULL} when \code{type_constraints} is \code{NULL}.
#' @param type_constraints a named integer vector containing type-specific size constraints.
#'                         If \code{NULL}, only the overall constraint will be checked.
#' @param primary_data_points a vector specifying primary data points, either by point indices
#'                            or with a logical vector of length equal to the number of points.
#'                            \code{NULL} indicates that all data points are "primary".
#'                            \code{check_scclust} checks so no primary data point is unassigned.
#'
#' @return Returns \code{TRUE} if \code{clustering} satisfies the clustering constraints, and
#'         \code{FALSE} if it does not. Throws an error if \code{clustering} is an invalid
#'         instance of the \code{\link{scclust}} class.
#'
#' @seealso See \code{\link{sc_clustering}} for details on how to specify the
#'          \code{type_labels} and \code{type_constraints} parameters.
#'
#' @examples
#' # Example scclust clustering
#' my_scclust <- scclust(c("A", "A", "B", "C", "B",
#'                         "C", "C", "A", "B", "B"))
#'
#'
#' # Check so each cluster contains at least two data points
#' check_scclust(my_scclust, 2)
#' # > TRUE
#'
#'
#' # Check so each cluster contains at least three data points
#' check_scclust(my_scclust, 3)
#' # > FALSE
#'
#'
#' # Data point types
#' my_types <- factor(c("x", "y", "y", "z", "z",
#'                      "x", "y", "z", "x", "x"))
#'
#'
#' # Check so each cluster contains at least one point of each type
#' check_scclust(my_scclust,
#'               NULL,
#'               my_types,
#'               c("x" = 1, "y" = 1, "z" = 1))
#' # > TRUE
#'
#'
#' # Check so each cluster contains one point of each type
#' # and at least three points in total
#' check_scclust(my_scclust,
#'               3,
#'               my_types,
#'               c("x" = 1, "y" = 1, "z" = 1))
#' # > TRUE
#'
#'
#' # Check so each cluster contains five data points of type "y"
#' check_scclust(my_scclust,
#'               NULL,
#'               my_types,
#'               c("y" = 5))
#' # > FALSE
#'
#'
#' # Check so each cluster contains one data point of both "x" and "z"
#' # and at least three points in total
#' check_scclust(my_scclust,
#'               3,
#'               my_types,
#'               c("x" = 1, "z" = 1))
#' # > TRUE
#'
#' @export
check_scclust <- function(clustering,
                          size_constraint = NULL,
                          type_labels = NULL,
                          type_constraints = NULL,
                          primary_data_points = NULL) {
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

  primary_data_points <- coerce_data_point_indices(primary_data_points, num_data_points)
  if (is.logical(primary_data_points)) {
    primary_data_points <- which(primary_data_points) - 1L
  } else if (is.integer(primary_data_points)) {
    primary_data_points <- primary_data_points - 1L
  }

  .Call(Rscc_check_scclust,
        clustering,
        size_constraint,
        unclass(type_labels),
        type_constraints,
        primary_data_points)
}


#' Derive statistics about a clustering.
#'
#' \code{get_scclust_stats} calculates statistics about a clustering.
#'
#' The function reports the following statistics:
#'
#' \tabular{ll}{
#'   \code{num_data_points} \tab The total number of data points. \cr
#'   \code{num_assigned} \tab The number of points assigned to a cluster. \cr
#'   \code{num_clusters} \tab The number of clusters. \cr
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
#' containing all indices of points assigned to the cluster. Let
#' \deqn{D(c) = \{d(i,j): i,j \in c \wedge i>j\}}{D(c) = { d(i,j) : i,j in c and i > j }}
#' be a set function returning all unique within-cluster distances in \eqn{c}{c}. Let \eqn{C}{C} be a set containing
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
#' @param clustering a \code{\link{scclust}} object containing a non-empty clustering.
#' @param distances a \code{\link[distances]{distances}} object describing the distances
#'                  between the data points in \code{clustering}.
#'
#' @return Returns a list of class \code{scclust_stats} containing the statistics.
#'
#' @examples
#' my_data_points <- data.frame(x = c(0.1, 0.2, 0.3, 0.4, 0.5,
#'                                    0.6, 0.7, 0.8, 0.9, 1.0),
#'                              y = c(10, 9, 8, 7, 6,
#'                                    10, 9, 8, 7, 6))
#'
#' my_distances <- distances(my_data_points)
#'
#' my_scclust <- scclust(c("A", "A", "B", "C", "B",
#'                         "C", "C", "A", "B", "B"))
#'
#' get_scclust_stats(my_distances, my_scclust)
#'
#' # >                        Value
#' # > num_data_points        10.0000000
#' # > num_assigned           10.0000000
#' # > num_clusters            3.0000000
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
get_scclust_stats <- function(distances,
                              clustering) {
  ensure_scclust(clustering)
  num_data_points <- length(clustering)
  ensure_distances(distances, num_data_points)

  clust_stats <- .Call(Rscc_get_scclust_stats,
                       distances,
                       clustering)
  structure(clust_stats,
            class = c("scclust_stats"))
}


#' @export
print.scclust_stats <- function(x, ...) {
  tmp_table <- as.table(format(as.matrix(unlist(x))))
  colnames(tmp_table) <- "Value"
  print(tmp_table)
}
