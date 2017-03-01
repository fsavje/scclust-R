# ==============================================================================
# scclust for R -- R wrapper for the scclust library
# https://github.com/fsavje/scclust-R
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
#' functions in the package.
#'
#' While \code{hierarchical_clustering} can be used alone to derive size constrained
#' clusterings, its main purpose is to be used with the nearest neighbor graph (NNG)
#' functions in the package (i.e., \code{\link{make_clustering}}). The clusterings
#' produced by the NNG functions tend produce large clusters in regions with many
#' data points. In some cases, it is beneficial to divide these clusters into smaller
#' groups. \code{hierarchical_clustering} can be use to achieve that.
#'
#' \code{hierarchical_clustering} implements a divisive hierarchical clustering
#' algorithm that respect size constraints. Starting from any clustering satisfying
#' the size constraints (which may be a clustering with a single cluster containing
#' all data points), \code{hierarchical_clustering} seraches for clusters that can
#' be broken into two or more new clusters without violating the constraints. It
#' continues in this fashion until all remaining clusters are unbreakable.
#'
#' \code{hierarchical_clustering} breaks a cluster in three stages. First,
#' it tries to find two data points as far as possible from each other in
#' the cluster it is about to break. The two points are called \emph{centers},
#' and they are the starting points for the two new clusters. The remaining
#' data points in the old cluster will be assigned to one of the centers.
#' In the second stage, each center picks the clostest data points so that
#' the new cluster satisfies the size constraint and assigns them to its
#' cluster. In the last stage, data points still not assigned are assigned to
#' their clostest centers either one-by-one or in batches (see below for
#' discussion). When all data points are assigned to a cluster, the old
#' cluster is removed and the two new are added to the clustering. If a new
#' cluster is breakable, it will go through the same procedure again.
#'
#' In some applications, it is desireable to avoid clusters that contain a number
#' of data points that are not multiples of \code{size_constraint}. After the
#' second stage described in the previous paragraph, both partial clusters are
#' exact multiples of the size constraint. By assigning remaining data points
#' in the third stage in batches of \code{size_constraint}, this property is
#' retained to the greatest extent possible.
#'
#'
#' @param distance_object a distance object as produced by \code{\link{make_distances}}.
#' @param size_constraint an integer with the required minimum cluster size.
#' @param batch_assign a bool indicating whether data points should be assigned in batches when
#'                            spliting clusters.
#' @param existing_clustering \code{NULL} or a \code{scc_clustering} object containing an existing
#'                            non-empty clustering. If \code{NULL}, the function will start with a
#'                            single cluster containing all data points (i.e., it derives a clustering
#'                            from scratch).
#'
#' @return Returns a scclust cluster object containing the derived clustering.
#'
#' @keywords cluster
#' @family clustering functions
#'
#' @seealso \code{\link{make_clustering}} is the main clustering function in the package.
#'
#'          Use \code{\link{scc_clustering}} to create scclust cluster objects from external
#'          clusterings.
#'
#' @export
hierarchical_clustering <- function(distance_object,
                                    size_constraint,
                                    batch_assign = TRUE,
                                    existing_clustering = NULL) {
  ensure_distances(distance_object)
  num_data_points <- data_point_count.scc_distances(distance_object)
  size_constraint <- coerce_size_constraint(size_constraint, num_data_points)
  ensure_indicators(batch_assign, 1L)
  if (!is.null(existing_clustering)) {
    ensure_scc_clustering(existing_clustering, num_data_points)
  }

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
#' @return Returns a scclust cluster object containing the derived clustering.
#'
#' @keywords internal
#'
#' @useDynLib scclust Rscc_hierarchical_clustering
hierarchical_clustering_internal <- function(distance_object,
                                             size_constraint,
                                             batch_assign = TRUE,
                                             existing_clustering = NULL,
                                             deep_copy = TRUE) {
  ensure_distances(distance_object)
  num_data_points <- data_point_count.scc_distances(distance_object)
  size_constraint <- coerce_size_constraint(size_constraint, num_data_points)
  ensure_indicators(batch_assign, 1L)
  if (!is.null(existing_clustering)) {
    ensure_scc_clustering(existing_clustering, num_data_points)
  }
  ensure_indicators(deep_copy, 1L)

  clustering <- .Call("Rscc_hierarchical_clustering",
                      distance_object,
                      size_constraint,
                      batch_assign,
                      existing_clustering,
                      deep_copy,
                      PACKAGE = "scclust")
  make_scc_clustering(clustering$cluster_labels,
                      clustering$cluster_count,
                      attr(distance_object, "ids", exact = TRUE))
}
