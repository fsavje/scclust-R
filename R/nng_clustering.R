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
  check_Rscc_distances(distance_object)
  num_data_points <- get_num_data_points(distance_object)
  size_constraint <- get_size_constraint(size_constraint)
  stopifnot(size_constraint <= num_data_points)
  seed_method <- get_seed_method(seed_method)
  main_unassigned_method <- match.arg(main_unassigned_method, c("ignore",
                                                                "by_nng",
                                                                "closest_assigned",
                                                                "closest_seed",
                                                                "estimated_radius_closest_seed"))
  main_radius <- get_radius(main_radius)
  check_main_data_points(main_data_points, num_data_points)
  if (is.null(main_data_points)) secondary_unassigned_method <- "ignore"
  secondary_unassigned_method <- match.arg(secondary_unassigned_method, c("ignore",
                                                                          "closest_assigned",
                                                                          "closest_seed",
                                                                          "estimated_radius_closest_seed"))
  secondary_radius <- get_radius(secondary_radius)

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


#' @useDynLib Rscclust Rsccwrap_nng_clustering_batches
#' @export
nng_clustering_batches <- function(distance_object,
                                   size_constraint,
                                   main_unassigned_method = "by_nng",
                                   main_radius = NULL,
                                   main_data_points = NULL,
                                   batch_size = 100L) {
  check_Rscc_distances(distance_object)
  num_data_points <- get_num_data_points(distance_object)
  size_constraint <- get_size_constraint(size_constraint)
  stopifnot(size_constraint <= num_data_points)
  main_unassigned_method <- match.arg(main_unassigned_method, c("ignore", "by_nng"))
  main_radius <- get_radius(main_radius)
  check_main_data_points(main_data_points, num_data_points)
  batch_size <- suppressWarnings(as.integer(batch_size)[1])
  stopifnot(!is.na(batch_size),
            batch_size >= 0L)

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
  check_Rscc_distances(distance_object)
  num_data_points <- get_num_data_points(distance_object)
  check_type_labels(type_labels)
  type_size_constraints <- get_type_size_constraints(type_size_constraints,
                                                     type_labels)
  total_size_constraint <- get_total_size_constraint(total_size_constraint,
                                                     type_size_constraints)
  stopifnot(total_size_constraint <= num_data_points)
  seed_method <- get_seed_method(seed_method)
  main_unassigned_method <- match.arg(main_unassigned_method, c("ignore",
                                                                "by_nng",
                                                                "closest_assigned",
                                                                "closest_seed",
                                                                "estimated_radius_closest_seed"))
  main_radius <- get_radius(main_radius)
  check_main_data_points(main_data_points, num_data_points)
  if (is.null(main_data_points)) secondary_unassigned_method <- "ignore"
  secondary_unassigned_method <- match.arg(secondary_unassigned_method, c("ignore",
                                                                          "closest_assigned",
                                                                          "closest_seed",
                                                                          "estimated_radius_closest_seed"))
  secondary_radius <- get_radius(secondary_radius)

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
