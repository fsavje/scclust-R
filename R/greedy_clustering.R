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
#' @export
top_down_greedy_clustering_internal <- function(distance_object,
                                                size_constraint,
                                                batch_assign = TRUE,
                                                existing_clustering = NULL,
                                                deep_copy = TRUE) {
  check_Rscc_distances(distance_object)
  num_data_points <- get_num_data_points(distance_object)
  existing_num_clusters <- 0L
  if (!is.null(existing_clustering)) {
    check_Rscc_clustering(existing_clustering)
    stopifnot(length(existing_clustering) == num_data_points)
    existing_num_clusters <- get_num_clusters(existing_clustering)
  }
  size_constraint <- get_size_constraint(size_constraint)
  stopifnot(size_constraint <= num_data_points)
  batch_assign <- get_bool_scalar(batch_assign)
  deep_copy <- get_bool_scalar(deep_copy)

  clustering <- .Call("Rsccwrap_top_down_greedy_clustering",
                      distance_object,
                      size_constraint,
                      batch_assign,
                      existing_clustering,
                      existing_num_clusters,
                      deep_copy,
                      PACKAGE = "Rscclust")
  make_Rscc_clustering(clustering$cluster_labels,
                       clustering$cluster_count,
                       attr(distance_object, "ids", exact = TRUE))
}
