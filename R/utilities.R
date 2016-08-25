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


#' @export
Rscc_clustering <- function(cluster_labels,
                            ids = NULL) {

  factor_labels <- factor(as.vector(cluster_labels))

  stopifnot(length(factor_labels) == 0,
            is.null(ids) || is.vector(ids),
            is.null(ids) || (length(factor_labels) == length(ids)))

  structure(as.integer(factor_labels) - 1,
            cluster_count = nlevels(factor_labels),
            ids = ids,
            class = c("Rscc_clustering"))
}


#' @useDynLib Rscclust Rsccwrap_check_clustering
#' @export
check_clustering <- function(clustering,
                             size_constraint) {

  stopifnot(inherits(clustering, "Rscc_clustering"),
            is.integer(clustering),
            "cluster_count" %in% names(attributes(clustering)))

  num_clusters <- as.integer(attr(clustering, "cluster_count", exact = TRUE))[1]
  size_constraint <- as.integer(size_constraint)[1]

  stopifnot(num_clusters > 0,
            size_constraint >= 2)

  .Call("Rsccwrap_check_clustering",
        clustering,
        num_clusters,
        size_constraint,
        PACKAGE = "Rscclust")
}


#' @useDynLib Rscclust Rsccwrap_check_clustering_types
#' @export
check_clustering_types <- function(clustering,
                                   type_labels,
                                   type_size_constraints,
                                   total_size_constraint = NULL) {

  stopifnot(inherits(clustering, "Rscc_clustering"),
            is.integer(clustering),
            "cluster_count" %in% names(attributes(clustering)),
            is.factor(type_labels) || is.integer(type_labels),
            all(!is.na(type_labels)))

  num_clusters <- as.integer(attr(clustering, "cluster_count", exact = TRUE))[1]

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

  stopifnot(num_clusters > 0,
            min(type_size_constraints) >= 0,
            total_size_constraint >= 2)

  .Call("Rsccwrap_check_clustering_types",
        clustering,
        num_clusters,
        total_size_constraint,
        type_size_constraints,
        type_labels,
        PACKAGE = "Rscclust")
}


#' @useDynLib Rscclust Rsccwrap_get_clustering_stats
#' @export
get_clustering_stats <- function(clustering,
                                 distance_object) {

  stopifnot(inherits(distance_object, "Rscc_distances"),
            is.numeric(distance_object),
            inherits(clustering, "Rscc_clustering"),
            is.integer(clustering),
            length(clustering) == ncol(distance_object),
            "cluster_count" %in% names(attributes(clustering)))

  num_clusters <- as.integer(attr(clustering, "cluster_count", exact = TRUE))[1]

  stopifnot(num_clusters > 0)

  clust_stats <- .Call("Rsccwrap_get_clustering_stats",
                       clustering,
                       num_clusters,
                       distance_object,
                       PACKAGE = "Rscclust")

  structure(clust_stats,
            class = c("Rscc_clustering_stats"))
}
