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


# Check Rscclust objects =======================================================

check_Rscc_distances <- function(distance_object) {
  stopifnot(inherits(distance_object, "Rscc_distances"),
            is.matrix(distance_object),
            is.numeric(distance_object))
}


check_Rscc_clustering <- function(clustering) {
  stopifnot(inherits(clustering, "Rscc_clustering"),
            is.integer(clustering),
            length(clustering) > 0L,
            "cluster_count" %in% names(attributes(clustering)),
            as.integer(attr(clustering, "cluster_count", exact = TRUE))[1] > 0L)
}


# Make Rscclust objects ========================================================

make_Rscc_clustering <- function(cluster_labels,
                                 cluster_count,
                                 ids) {
  structure(cluster_labels,
            cluster_count = cluster_count,
            ids = ids,
            class = c("Rscc_clustering"))
}


# Check inputs =================================================================

check_main_data_points <- function(main_data_points,
                                   num_data_points) {
  if (!is.null(main_data_points)) {
    stopifnot(is.logical(main_data_points),
              length(main_data_points) >= num_data_points)
  }
}


check_type_labels <- function(type_labels) {
  stopifnot(is.factor(type_labels) || is.integer(type_labels),
            !is.factor(type_labels) || all(!is.na(type_labels)),
            !is.integer(type_labels) || isTRUE(min(type_labels, na.rm = FALSE) >= 0L))
}


# Get parameters ===============================================================

get_num_data_points <- function(distance_object) {
  ncol(distance_object)
}


get_size_constraint <- function(size_constraint) {
  stopifnot(is.numeric(size_constraint))
  stopifnot(as.integer(size_constraint)[1] >= 2L)
  as.integer(size_constraint)[1]
}


get_seed_method <- function(seed_method) {
  match.arg(seed_method, c("lexical",
                           "inwards_order",
                           "inwards_updating",
                           "exclusion_order",
                           "exclusion_updating"))
}


get_radius <- function(radius) {
  if (!is.null(radius)) {
    stopifnot(is.numeric(radius))
    radius <- as.numeric(radius)[1]
    stopifnot(radius > 0.0)
  }
  radius
}


get_num_clusters <- function(clustering) {
  as.integer(attr(clustering, "cluster_count", exact = TRUE))[1]
}


get_type_size_constraints <- function(type_size_constraints,
                                      type_labels) {
  stopifnot(!is.null(names(type_size_constraints)),
            !anyDuplicated(names(type_size_constraints)))

  out_type_size_constraints <- NULL

  if (is.factor(type_labels)) {
    out_type_size_constraints <- rep(0L, nlevels(type_labels))
    names(out_type_size_constraints) <- levels(type_labels)
    stopifnot(all(names(type_size_constraints) %in% names(out_type_size_constraints)))
    out_type_size_constraints[names(type_size_constraints)] <- as.integer(type_size_constraints)
    out_type_size_constraints <- c(0L, out_type_size_constraints)
  } else if (is.integer(type_labels)) {
    max_label <- max(type_labels)
    out_type_size_constraints <- rep(0L, max_label + 1L)
    names(out_type_size_constraints) <- as.character(0L:max_label)
    stopifnot(all(names(type_size_constraints) %in% names(out_type_size_constraints)))
    out_type_size_constraints[names(type_size_constraints)] <- as.integer(type_size_constraints)
  }

  stopifnot(min(out_type_size_constraints) >= 0L)

  unname(out_type_size_constraints)
}


get_total_size_constraint <- function(total_size_constraint,
                                      type_size_constraints) {
  if (is.null(total_size_constraint)) {
    total_size_constraint <- as.integer(sum(type_size_constraints))
  } else {
    total_size_constraint <- as.integer(total_size_constraint)[1]
    stopifnot(total_size_constraint >= as.integer(sum(type_size_constraints)))
  }

  stopifnot(total_size_constraint >= 2L)

  total_size_constraint
}
