# ==============================================================================
# Rscclust -- R wrapper for the scclust library
# https://github.com/fsavje/Rscclust
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


#' Nearest neighbor graph based clustering function
#'
#' \code{nng_clustering} derives a clustering statisfying the specified
#' size constraints. It implements an algorithm that first summaries the
#' distance information between data points in a sparse graph and then
#' constructs the clusterings based on the graph. The function is intended
#' to run fast while ensuring near-optimal performance.
#'
#' Must be one of
#'                    "lexical", "inwards_order", "inwards_updating",
#'                    "inwards_alt_updating", "exclusion_order", "exclusion_updating".
#'
#'  Must be one of
#'                               "ignore", "by_nng", "closest_assigned",
#'                               "closest_seed", "estimated_radius_closest_seed".
#'
#' @param distance_object a distance object as produced by \code{\link{make_distances}}.
#' @param size_constraint an integer with the required minimum cluster size.
#' @param seed_method how seeds are found. See below for additional details.
#' @param unassigned_method how remaining (primary) data points are assigned
#'                               to clusters. See below for additional details.
#' @param radius restricts the maximum length of an edge in the graph that
#'                    summaries the distance information. \code{NULL} indicates
#'                    no restriction. See below for details on how this translates
#'                    to maximum distance within a cluster.
#' @param primary_data_points a logical vector of length equal to number of data points
#'                            indicating primary and secondary data points. \code{NULL}
#'                            indicates that all data points are "primary". See below for
#'                            details.
#' @param secondary_unassigned_method how remaining secondary data points are assigned
#'                               to clusters. See below for additional details.
#' @param secondary_radius restricts the maximum distance when assigning secondary data
#'                         points to clusters.
#'
#' @return Returns a Rscclust cluster object containing the derived clustering.
#'
#' @keywords cluster
#' @family clustering functions
#'
#' @seealso \code{\link{nng_clustering_batches}} and \code{\link{nng_clustering_types}} are
#'          other NNG-based clustering functions in the package.
#'
#'          \code{\link{hierarchical_clustering}} can be used to refine the clustering
#'          constructed by this function.
#'
#' @useDynLib Rscclust Rscc_nng_clustering
#' @export
nng_clustering <- function(distance_object,
                           size_constraint,
                           seed_method = "exclusion_updating",
                           unassigned_method = "closest_seed",
                           radius = NULL,
                           primary_data_points = NULL,
                           secondary_unassigned_method = "ignore",
                           secondary_radius = NULL) {
  ensure_distances(distance_object)
  num_data_points <- data_point_count.Rscc_distances(distance_object)
  size_constraint <- coerce_size_constraint(size_constraint, num_data_points)
  seed_method <- coerce_args(seed_method, all_seed_methods)
  unassigned_method <- coerce_args(unassigned_method,
                                   c("ignore",
                                     "by_nng",
                                     "closest_assigned",
                                     "closest_seed",
                                     "estimated_radius_closest_seed"))
  radius <- coerce_radius(radius)
  if (is.null(primary_data_points)) {
    secondary_unassigned_method <- "ignore"
  } else {
    ensure_indicators(primary_data_points, num_data_points, TRUE)
  }
  secondary_unassigned_method <- coerce_args(secondary_unassigned_method,
                                             c("ignore",
                                               "closest_assigned",
                                               "closest_seed",
                                               "estimated_radius_closest_seed"))
  secondary_radius <- coerce_radius(secondary_radius)

  clustering <- .Call("Rscc_nng_clustering",
                      distance_object,
                      size_constraint,
                      seed_method,
                      unassigned_method,
                      radius,
                      primary_data_points,
                      secondary_unassigned_method,
                      secondary_radius,
                      PACKAGE = "Rscclust")
  make_Rscc_clustering(clustering$cluster_labels,
                       clustering$cluster_count,
                       attr(distance_object, "ids", exact = TRUE))
}



#' Nearest neighbor graph based clustering function in batch mode
#'
#' Dummy documentation
#'
#' @param distance_object abc ...
#' @param size_constraint abc ...
#' @param unassigned_method abc ...
#' @param radius abc ...
#' @param primary_data_points abc ...
#' @param batch_size abc ...
#'
#' @useDynLib Rscclust Rscc_nng_clustering_batches
#' @export
nng_clustering_batches <- function(distance_object,
                                   size_constraint,
                                   unassigned_method = "by_nng",
                                   radius = NULL,
                                   primary_data_points = NULL,
                                   batch_size = 100L) {
  ensure_distances(distance_object)
  num_data_points <- data_point_count.Rscc_distances(distance_object)
  size_constraint <- coerce_size_constraint(size_constraint, num_data_points)
  unassigned_method <- coerce_args(unassigned_method,
                                   c("ignore", "by_nng"))
  radius <- coerce_radius(radius)
  if (!is.null(primary_data_points)) {
    ensure_indicators(primary_data_points, num_data_points, TRUE)
  }
  batch_size <- coerce_counts(batch_size, 1L)

  clustering <- .Call("Rscc_nng_clustering_batches",
                      distance_object,
                      size_constraint,
                      unassigned_method,
                      radius,
                      primary_data_points,
                      batch_size,
                      PACKAGE = "Rscclust")
  make_Rscc_clustering(clustering$cluster_labels,
                       clustering$cluster_count,
                       attr(distance_object, "ids", exact = TRUE))
}


#' Nearest neighbor graph based clustering function with type constraints
#'
#' Dummy documentation
#'
#' @param distance_object abc ...
#' @param type_labels abc ...
#' @param type_size_constraints abc ...
#' @param total_size_constraint abc ...
#' @param seed_method abc ...
#' @param unassigned_method abc ...
#' @param radius abc ...
#' @param primary_data_points abc ...
#' @param secondary_unassigned_method abc ...
#' @param secondary_radius abc ...
#'
#' @useDynLib Rscclust Rscc_nng_clustering_types
#' @export
nng_clustering_types <- function(distance_object,
                                 type_labels,
                                 type_size_constraints,
                                 total_size_constraint = NULL,
                                 seed_method = "exclusion_updating",
                                 unassigned_method = "closest_seed",
                                 radius = NULL,
                                 primary_data_points = NULL,
                                 secondary_unassigned_method = "ignore",
                                 secondary_radius = NULL) {
  ensure_distances(distance_object)
  num_data_points <- data_point_count.Rscc_distances(distance_object)
  type_labels <- coerce_type_labels(type_labels, num_data_points)
  type_size_constraints <- coerce_type_constraints(type_size_constraints)
  type_size_constraints <- make_type_size_constraints(type_size_constraints,
                                                      type_labels)
  total_size_constraint <- coerce_total_size_constraint(total_size_constraint,
                                                        type_size_constraints,
                                                        num_data_points)
  seed_method <- coerce_args(seed_method, all_seed_methods)
  unassigned_method <- coerce_args(unassigned_method,
                                   c("ignore",
                                     "by_nng",
                                     "closest_assigned",
                                     "closest_seed",
                                     "estimated_radius_closest_seed"))
  radius <- coerce_radius(radius)
  if (is.null(primary_data_points)) {
    secondary_unassigned_method <- "ignore"
  } else {
    ensure_indicators(primary_data_points, num_data_points, TRUE)
  }
  secondary_unassigned_method <- coerce_args(secondary_unassigned_method,
                                             c("ignore",
                                               "closest_assigned",
                                               "closest_seed",
                                               "estimated_radius_closest_seed"))
  secondary_radius <- coerce_radius(secondary_radius)

  clustering <- .Call("Rscc_nng_clustering_types",
                      distance_object,
                      unclass(type_labels),
                      type_size_constraints,
                      total_size_constraint,
                      seed_method,
                      unassigned_method,
                      radius,
                      primary_data_points,
                      secondary_unassigned_method,
                      secondary_radius,
                      PACKAGE = "Rscclust")
  make_Rscc_clustering(clustering$cluster_labels,
                       clustering$cluster_count,
                       attr(distance_object, "ids", exact = TRUE))
}
