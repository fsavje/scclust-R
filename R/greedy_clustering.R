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
top_down_greedy_clustering <- function(distance_object,
                                       size_constraint,
                                       batch_assign = TRUE,
                                       existing_clustering = NULL) {
  top_down_greedy_clustering_internal(distance_object,
                                      size_constraint,
                                      batch_assign,
                                      existing_clustering,
                                      deep_copy = TRUE)
}


#' @useDynLib Rscclust Rsccwrap_top_down_greedy_clustering
top_down_greedy_clustering_internal <- function(distance_object,
                                                size_constraint,
                                                batch_assign = TRUE,
                                                existing_clustering = NULL,
                                                deep_copy = TRUE) {

  stopifnot(inherits(distance_object, "Rscc_distances"))
  num_data_points <- ncol(distance_object)
  existing_num_clusters <- 0L
  if (!is.null(existing_clustering)) {
    stopifnot(inherits(existing_clustering, "Rscc_clustering"),
              "cluster_count" %in% names(attributes(cl3)))
    existing_num_clusters <- as.integer(attr(existing_clustering, "cluster_count", exact = TRUE))[1]
  }

  size_constraint <- as.integer(size_constraint)[1]
  batch_assign <- as.logical(batch_assign)[1]
  deep_copy <- as.logical(deep_copy)[1]

  stopifnot(is.null(existing_clustering) || inherits(existing_clustering, "Rscc_clustering"),
            is.null(existing_clustering) || (existing_num_clusters > 0),
            !is.null(existing_clustering) || (existing_num_clusters == 0),
            size_constraint >= 2,
            size_constraint <= num_data_points)

  clustering <- .Call("Rsccwrap_top_down_greedy_clustering",
                      distance_object,
                      size_constraint,
                      batch_assign,
                      unclass(existing_clustering),
                      existing_num_clusters,
                      deep_copy,
                      PACKAGE = "Rscclust")

  structure(clustering$cluster_labels,
            cluster_count = clustering$cluster_count,
            ids = attr(distance_object, "ids", exact = TRUE),
            class = c("Rscc_clustering"))
}
