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


#' @useDynLib Rscclust Rsccwrap_check_clustering
#' @export
check_clustering <- function(clustering,
                             size_constraint) {
  check_Rscc_clustering(clustering)
  num_clusters <- get_num_clusters(clustering)
  size_constraint <- get_size_constraint(size_constraint)

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
  check_Rscc_clustering(clustering)
  check_type_labels(type_labels)
  num_clusters <- get_num_clusters(clustering)
  type_size_constraints <- get_type_size_constraints(type_size_constraints,
                                                     type_labels)
  total_size_constraint <- get_total_size_constraint(total_size_constraint,
                                                     type_size_constraints)

  .Call("Rsccwrap_check_clustering_types",
        clustering,
        num_clusters,
        total_size_constraint,
        type_size_constraints,
        unclass(type_labels),
        PACKAGE = "Rscclust")
}


#' @useDynLib Rscclust Rsccwrap_get_clustering_stats
#' @export
get_clustering_stats <- function(clustering,
                                 distance_object) {
  check_Rscc_clustering(clustering)
  num_clusters <- get_num_clusters(clustering)
  check_Rscc_distances(distance_object)
  num_data_points <- get_num_data_points(distance_object)
  stopifnot(length(clustering) == num_data_points)

  clust_stats <- .Call("Rsccwrap_get_clustering_stats",
                       clustering,
                       num_clusters,
                       distance_object,
                       PACKAGE = "Rscclust")
  structure(clust_stats,
            class = c("Rscc_clustering_stats"))
}


#' @export
print.Rscc_clustering_stats <- function(x, ...) {
  tmp_table <- as.table(format(as.matrix(unlist(x))))
  colnames(tmp_table) <- "Value"
  print(tmp_table)
}
