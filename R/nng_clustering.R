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


#' Dummy documentation
#'
#' Dummy documentation
#'
#' @param distance_object abc ...
#' @param size_constraint abc ...
#' @param seed_method abc ...
#' @param main_unassigned_method abc ...
#' @param main_radius abc ...
#' @param main_data_points abc ...
#' @param secondary_unassigned_method abc ...
#' @param secondary_radius abc ...
#'
#' @keywords cluster
#' @family clustering functions
#' @useDynLib Rscclust Rsccwrap_nng_clustering
#' @export
nng_clustering <- function(distance_object,
                           size_constraint,
                           seed_method = "exclusion_updating",
                           main_unassigned_method = "closest_seed",
                           main_radius = NULL,
                           main_data_points = NULL,
                           secondary_unassigned_method = "ignore",
                           secondary_radius = NULL) {
  ensure_distances(distance_object)
  num_data_points <- data_point_count_distances(distance_object)
  size_constraint <- coerce_size_constraint(size_constraint, num_data_points)
  seed_method <- coerce_args(seed_method, all_seed_methods)
  main_unassigned_method <- coerce_args(main_unassigned_method,
                                        c("ignore",
                                          "by_nng",
                                          "closest_assigned",
                                          "closest_seed",
                                          "estimated_radius_closest_seed"))
  main_radius <- coerce_radius(main_radius)
  if (is.null(main_data_points)) {
    secondary_unassigned_method <- "ignore"
  } else {
    ensure_indicators(main_data_points, num_data_points, TRUE)
  }
  secondary_unassigned_method <- coerce_args(secondary_unassigned_method,
                                             c("ignore",
                                               "closest_assigned",
                                               "closest_seed",
                                               "estimated_radius_closest_seed"))
  secondary_radius <- coerce_radius(secondary_radius)

  clustering <- .Call("Rsccwrap_nng_clustering",
                      distance_object,
                      size_constraint,
                      seed_method,
                      main_unassigned_method,
                      main_radius,
                      main_data_points,
                      secondary_unassigned_method,
                      secondary_radius,
                      PACKAGE = "Rscclust")
  make_Rscc_clustering(clustering$cluster_labels,
                       clustering$cluster_count,
                       attr(distance_object, "ids", exact = TRUE))
}


#' Dummy documentation
#'
#' Dummy documentation
#'
#' @param distance_object abc ...
#' @param size_constraint abc ...
#' @param main_unassigned_method abc ...
#' @param main_radius abc ...
#' @param main_data_points abc ...
#' @param batch_size abc ...
#'
#' @useDynLib Rscclust Rsccwrap_nng_clustering_batches
#' @export
nng_clustering_batches <- function(distance_object,
                                   size_constraint,
                                   main_unassigned_method = "by_nng",
                                   main_radius = NULL,
                                   main_data_points = NULL,
                                   batch_size = 100L) {
  ensure_distances(distance_object)
  num_data_points <- data_point_count_distances(distance_object)
  size_constraint <- coerce_size_constraint(size_constraint, num_data_points)
  main_unassigned_method <- coerce_args(main_unassigned_method,
                                        c("ignore", "by_nng"))
  main_radius <- coerce_radius(main_radius)
  if (!is.null(main_data_points)) {
    ensure_indicators(main_data_points, num_data_points, TRUE)
  }
  batch_size <- coerce_counts(batch_size, 1L)

  clustering <- .Call("Rsccwrap_nng_clustering_batches",
                      distance_object,
                      size_constraint,
                      main_unassigned_method,
                      main_radius,
                      main_data_points,
                      batch_size,
                      PACKAGE = "Rscclust")
  make_Rscc_clustering(clustering$cluster_labels,
                       clustering$cluster_count,
                       attr(distance_object, "ids", exact = TRUE))
}


#' Dummy documentation
#'
#' Dummy documentation
#'
#' @param distance_object abc ...
#' @param type_labels abc ...
#' @param type_size_constraints abc ...
#' @param total_size_constraint abc ...
#' @param seed_method abc ...
#' @param main_unassigned_method abc ...
#' @param main_radius abc ...
#' @param main_data_points abc ...
#' @param secondary_unassigned_method abc ...
#' @param secondary_radius abc ...
#'
#' @useDynLib Rscclust Rsccwrap_nng_clustering_types
#' @export
nng_clustering_types <- function(distance_object,
                                 type_labels,
                                 type_size_constraints,
                                 total_size_constraint = NULL,
                                 seed_method = "exclusion_updating",
                                 main_unassigned_method = "closest_seed",
                                 main_radius = NULL,
                                 main_data_points = NULL,
                                 secondary_unassigned_method = "ignore",
                                 secondary_radius = NULL) {
  ensure_distances(distance_object)
  num_data_points <- data_point_count_distances(distance_object)
  type_labels <- coerce_type_labels(type_labels, num_data_points)
  type_size_constraints <- coerce_type_constraints(type_size_constraints)
  #ensure_type_labels_exist(names(type_size_constraints), get_all_types(type_labels))
  type_size_constraints <- make_type_size_constraints(type_size_constraints,
                                                      type_labels)
  total_size_constraint <- coerce_total_size_constraint(total_size_constraint,
                                                        type_size_constraints,
                                                        num_data_points)
  seed_method <- coerce_args(seed_method, all_seed_methods)
  main_unassigned_method <- coerce_args(main_unassigned_method,
                                        c("ignore",
                                          "by_nng",
                                          "closest_assigned",
                                          "closest_seed",
                                          "estimated_radius_closest_seed"))
  main_radius <- coerce_radius(main_radius)
  if (is.null(main_data_points)) {
    secondary_unassigned_method <- "ignore"
  } else {
    ensure_indicators(main_data_points, num_data_points, TRUE)
  }
  secondary_unassigned_method <- coerce_args(secondary_unassigned_method,
                                             c("ignore",
                                               "closest_assigned",
                                               "closest_seed",
                                               "estimated_radius_closest_seed"))
  secondary_radius <- coerce_radius(secondary_radius)

  clustering <- .Call("Rsccwrap_nng_clustering_types",
                      distance_object,
                      total_size_constraint,
                      type_size_constraints,
                      unclass(type_labels),
                      seed_method,
                      main_unassigned_method,
                      main_radius,
                      main_data_points,
                      secondary_unassigned_method,
                      secondary_radius,
                      PACKAGE = "Rscclust")
  make_Rscc_clustering(clustering$cluster_labels,
                       clustering$cluster_count,
                       attr(distance_object, "ids", exact = TRUE))
}
