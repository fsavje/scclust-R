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

  stopifnot(inherits(distance_object, "Rscc_distances"),
            is.numeric(distance_object))
  num_data_points <- ncol(distance_object)

  size_constraint <- as.integer(size_constraint)[1]
  seed_method <- match.arg(seed_method, c("lexical",
                                          "inwards_order",
                                          "inwards_updating",
                                          "exclusion_order",
                                          "exclusion_updating"))
  main_unassigned_method <- match.arg(main_unassigned_method, c("ignore",
                                                                "by_nng",
                                                                "closest_assigned",
                                                                "closest_seed",
                                                                "estimated_radius_closest_seed"))
  if (!is.null(main_radius)) main_radius <- as.numeric(main_radius)[1]
  if (is.null(main_data_points)) secondary_unassigned_method <- "ignore"
  secondary_unassigned_method <- match.arg(secondary_unassigned_method, c("ignore",
                                                                          "closest_assigned",
                                                                          "closest_seed",
                                                                          "estimated_radius_closest_seed"))
  if (!is.null(secondary_radius)) secondary_radius <- as.numeric(secondary_radius)[1]

  stopifnot(size_constraint >= 2,
            size_constraint <= num_data_points,
            is.null(main_radius) || (main_radius > 0.0),
            is.null(main_data_points) || is.logical(main_data_points),
            !is.logical(main_data_points) || (length(x) < num_data_points),
            is.null(secondary_radius) || (secondary_radius > 0.0))

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

  structure(clustering$cluster_labels,
            cluster_count = clustering$cluster_count,
            ids = attr(distance_object, "ids", exact = TRUE),
            class = c("Rscc_clustering"))
}


#' @useDynLib Rscclust Rsccwrap_nng_clustering_batches
#' @export
nng_clustering_batches <- function(distance_object,
                                   size_constraint,
                                   main_unassigned_method = "by_nng",
                                   main_radius = NULL,
                                   main_data_points = NULL,
                                   batch_size = 100) {

  stopifnot(inherits(distance_object, "Rscc_distances"),
            is.numeric(distance_object))
  num_data_points <- ncol(distance_object)

  size_constraint <- as.integer(size_constraint)[1]
  main_unassigned_method <- match.arg(main_unassigned_method, c("ignore", "by_nng"))
  if (!is.null(main_radius)) main_radius <- as.numeric(main_radius)[1]
  batch_size <- as.integer(batch_size)[1]

  stopifnot(size_constraint >= 2,
            size_constraint <= num_data_points,
            is.null(main_radius) || (main_radius > 0.0),
            is.null(main_data_points) || is.logical(main_data_points),
            !is.logical(main_data_points) || (length(x) < num_data_points),
            batch_size >= 0)

  clustering <- .Call("Rsccwrap_nng_clustering_batches",
                      distance_object,
                      size_constraint,
                      main_unassigned_method,
                      main_radius,
                      main_data_points,
                      batch_size,
                      PACKAGE = "Rscclust")

  structure(clustering$cluster_labels,
            cluster_count = clustering$cluster_count,
            ids = attr(distance_object, "ids", exact = TRUE),
            class = c("Rscc_clustering"))
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

  stopifnot(inherits(distance_object, "Rscc_distances"),
            is.numeric(distance_object))
  num_data_points <- ncol(distance_object)

  stopifnot(is.factor(type_labels) || is.integer(type_labels),
            all(!is.na(type_labels)))

  if (is.factor(type_labels)) {
    stopifnot(!is.null(names(type_size_constraints)),
              !anyDuplicated(names(type_size_constraints)),
              all(names(type_size_constraints) %in% levels(type_labels)))
    save_type_size_constraints <- type_size_constraints
    type_size_constraints <- rep(0L, nlevels(type_labels))
    names(type_size_constraints) <- levels(type_labels)
    type_size_constraints[names(save_type_size_constraints)] <- as.integer(save_type_size_constraints)
    type_size_constraints <- c(0L, type_size_constraints)
    rm(save_type_size_constraints)
  } else {
    stopifnot(min(type_labels) >= 0,
              !is.null(names(type_size_constraints)),
              !anyDuplicated(names(type_size_constraints)))
    type_max <- max(type_labels)
    save_type_size_constraints <- type_size_constraints
    type_size_constraints <- rep(0L, type_max + 1)
    names(type_size_constraints) <- as.character(0:type_max)
    stopifnot(all(names(save_type_size_constraints) %in% names(type_size_constraints)))
    type_size_constraints[names(save_type_size_constraints)] <- as.integer(save_type_size_constraints)
    rm(save_type_size_constraints)
  }

  if (is.null(total_size_constraint)) {
    total_size_constraint <- as.integer(sum(type_size_constraints))
  } else {
    total_size_constraint <- as.integer(total_size_constraint)[1]
  }

  seed_method <- match.arg(seed_method, c("lexical",
                                          "inwards_order",
                                          "inwards_updating",
                                          "exclusion_order",
                                          "exclusion_updating"))
  main_unassigned_method <- match.arg(main_unassigned_method, c("ignore",
                                                                "by_nng",
                                                                "closest_assigned",
                                                                "closest_seed",
                                                                "estimated_radius_closest_seed"))
  if (!is.null(main_radius)) main_radius <- as.numeric(main_radius)[1]
  if (is.null(main_data_points)) secondary_unassigned_method <- "ignore"
  secondary_unassigned_method <- match.arg(secondary_unassigned_method, c("ignore",
                                                                          "closest_assigned",
                                                                          "closest_seed",
                                                                          "estimated_radius_closest_seed"))
  if (!is.null(secondary_radius)) secondary_radius <- as.numeric(secondary_radius)[1]

  stopifnot(min(type_size_constraints) >= 0,
            total_size_constraint >= 2,
            total_size_constraint <= num_data_points,
            is.null(main_radius) || (main_radius > 0.0),
            is.null(main_data_points) || is.logical(main_data_points),
            !is.logical(main_data_points) || (length(x) < num_data_points),
            is.null(secondary_radius) || (secondary_radius > 0.0))

  clustering <- .Call("Rsccwrap_nng_clustering_types",
                      distance_object,
                      total_size_constraint,
                      type_size_constraints,
                      type_labels,
                      seed_method,
                      main_unassigned_method,
                      main_radius,
                      main_data_points,
                      secondary_unassigned_method,
                      secondary_radius,
                      PACKAGE = "Rscclust")

  structure(clustering$cluster_labels,
            cluster_count = clustering$cluster_count,
            ids = attr(distance_object, "ids", exact = TRUE),
            class = c("Rscc_clustering"))
}
