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
#' \code{make_clustering} derives a clustering statisfying the specified
#' size constraints. It implements an algorithm that first summaries the
#' distance information between data points in a sparse graph and then
#' constructs the clusterings based on the graph. The function is intended
#' to run fast while ensuring near-optimal performance.
#'
#' @param distance_object a distance object as produced by \code{\link{make_distances}}.
#' @param size_constraint an integer with the required minimum cluster size.
#' @param type_labels ...
#' @param type_constraints ...
#' @param seed_method how seeds are found. See below for additional details.
#' @param primary_data_points a logical vector of length equal to number of data points
#'                            indicating primary and secondary data points. \code{NULL}
#'                            indicates that all data points are "primary". See below for
#'                            details.
#' @param primary_unassigned_method how remaining (primary) data points are assigned
#'                               to clusters. See below for additional details.
#' @param secondary_unassigned_method how remaining secondary data points are assigned
#'                               to clusters. See below for additional details.
#' @param seed_radius restricts the maximum length of an edge in the graph that
#'                    summaries the distance information. \code{NULL} indicates
#'                    no restriction. See below for details on how this translates
#'                    to maximum distance within a cluster.
#' @param primary_radius ...
#' @param secondary_radius restricts the maximum distance when assigning secondary data
#'                         points to clusters.
#' @param batch_size ...
#'
#' @return Returns a Rscclust cluster object containing the derived clustering.
#'
#' @keywords cluster
#' @family clustering functions
#'
#' @seealso \code{\link{hierarchical_clustering}} can be used to refine the clustering
#'          constructed by this function.
#'
#' @useDynLib Rscclust Rscc_make_clustering
#' @export
make_clustering <- function(distance_object,
                            size_constraint = NULL,
                            type_labels = NULL,
                            type_constraints = NULL,
                            seed_method = "exclusion_updating",
                            primary_data_points = NULL,
                            primary_unassigned_method = "closest_seed",
                            secondary_unassigned_method = "ignore",
                            seed_radius = NULL,
                            primary_radius = "seed_radius",
                            secondary_radius = NULL,
                            batch_size = 100L) {
  ensure_distances(distance_object)
  num_data_points <- data_point_count.Rscc_distances(distance_object)

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

  seed_method <- coerce_args(seed_method, all_seed_methods)

  if (is.null(primary_data_points)) {
    secondary_unassigned_method <- "ignore"
  } else {
    ensure_indicators(primary_data_points, num_data_points, TRUE)
  }

  primary_unassigned_method <- coerce_args(primary_unassigned_method,
                                           c("ignore",
                                             "by_nng",
                                             "closest_assigned",
                                             "closest_seed",
                                             "estimated_radius_closest_seed"))
  secondary_unassigned_method <- coerce_args(secondary_unassigned_method,
                                             c("ignore",
                                               "closest_assigned",
                                               "closest_seed",
                                               "estimated_radius_closest_seed"))

  seed_radius <- coerce_radius(seed_radius, TRUE)
  primary_radius <- coerce_radius(primary_radius, FALSE)
  secondary_radius <- coerce_radius(secondary_radius, FALSE)

  if (!is.null(batch_size)) {
    batch_size <- coerce_counts(batch_size, 1L)
  }

  clustering <- .Call("Rscc_make_clustering",
                      distance_object,
                      size_constraint,
                      unclass(type_labels),
                      type_constraints,
                      seed_method,
                      primary_data_points,
                      primary_unassigned_method,
                      secondary_unassigned_method,
                      seed_radius,
                      primary_radius,
                      secondary_radius,
                      batch_size,
                      PACKAGE = "Rscclust")

  make_Rscc_clustering(clustering$cluster_labels,
                       clustering$cluster_count,
                       attr(distance_object, "ids", exact = TRUE))
}
