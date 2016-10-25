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


#' Hierarchical clustering function
#'
#' \code{hierarchical_clustering} derives a clustering statisfying a specified
#' size-constraint using a hierarchical clustering algorithm. The primary purpose of
#' this function is to refine clusterings produced by the nearest neighbor graph
#' algorithms in the packages.
#'
#' To be written...
#'
#' @param distance_object a distance object as produced by \code{\link{make_distances}}.
#' @param size_constraint an integer with the required minimum cluster size.
#' @param batch_assign a bool indicating whether data points should be assigned in batches when
#'                            spliting clusters.
#' @param existing_clustering \code{NULL} or a \code{Rscc_clustering} object containing an existing
#'                            non-empty clustering.
#'
#' @return Returns a Rscclust cluster object containing the derived clustering.
#'
#' @keywords cluster
#' @family clustering functions
#'
#' @seealso \code{\link{nng_clustering}}, \code{\link{nng_clustering_batches}} and
#'          \code{\link{nng_clustering_types}} are other clustering functions in the package.
#'
#'          Use \code{\link{Rscc_clustering}} to create Rscclust cluster objects from external
#'          clusterings.
#'
#' @export
hierarchical_clustering <- function(distance_object,
                                    size_constraint,
                                    batch_assign = TRUE,
                                    existing_clustering = NULL) {
  hierarchical_clustering_internal(distance_object,
                                   size_constraint,
                                   batch_assign,
                                   existing_clustering,
                                   deep_copy = TRUE)
}


#' Internal hierarchical clustering function
#'
#' ATTENTION! Using this function has side-effects on your R environment,
#' ensure you properly understand what this does before using it.
#' In most cases, \code{\link{hierarchical_clustering}} is the preferred
#' choice.
#'
#' \code{hierarchical_clustering_internal} is identical to
#' \code{\link{hierarchical_clustering}} except for the
#' \code{deep_copy} parameter. When \code{deep_copy == TRUE}, the
#' two functions are identical. However, when
#' \code{deep_copy == FALSE}, \code{hierarchical_clustering_internal}
#' will \strong{not} make a deep copy of \code{existing_clustering} before
#' refining the clustering. That is, any changes will be made in place, and
#' \code{existing_clustering} will be invalidated by calling
#' this function. This is useful when refining an
#' existing clustering and the inputted clustering is not of interest.
#'
#' See \code{\link{hierarchical_clustering}} for detailed documentation.
#'
#' @param deep_copy a bool indicating whether a deep copy of
#'                  \code{existing_clustering} should be made before
#'                  refining the clustering.
#' @inheritParams hierarchical_clustering
#'
#' @return Returns a Rscclust cluster object containing the derived clustering.
#'
#' @keywords internal
#'
#' @useDynLib Rscclust Rsccwrap_hierarchical_clustering
#' @export
hierarchical_clustering_internal <- function(distance_object,
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

  clustering <- .Call("Rsccwrap_hierarchical_clustering",
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
