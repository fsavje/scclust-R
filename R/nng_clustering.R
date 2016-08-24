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

  stopifnot(inherits(distance_object, "Rscc_distances"))
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
