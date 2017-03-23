# ==============================================================================
# scclust for R -- R wrapper for the scclust library
# https://github.com/fsavje/scclust-R
#
# Copyright (C) 2016-2017  Fredrik Savje -- http://fredriksavje.com
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


#' Hierarchical size-constrained clustering
#'
#' \code{hierarchical_clustering} serves two purposes. Its primary use
#' is to refine existing clusterings subject to clustering
#' constraints. That is, given a (non-optimal) clustering satisfying
#' some constraints, the function splits clusters so to decrease within-cluster
#' distances without violating the constraints. The function can also be used to
#' derive size-constrained clusterings from scratch. In both cases, it
#' uses a hierarchical clustering algorithm.
#'
#' While \code{hierarchical_clustering} can be used alone to derive size-constrained
#' clusterings, its main purpose is to be used together with \code{\link{sc_clustering}}.
#' The clusterings produced by the main clustering function are guaranteed to be close to
#' optimal (in particular, within a constant factor of the optimal solution). However,
#' as \code{\link{sc_clustering}}'s solution is not exactly optimal, it is occasionally
#' possible to refine it. In particular, the main function tends produce large clusters
#' in regions of the metric space with many data points. In some cases, it is beneficial
#' to divide these clusters into smaller groups. \code{hierarchical_clustering} splits
#' these larger clusters so that all within-cluster distances weakly decreases while
#' respecting the overall the size constraint.
#'
#' \code{hierarchical_clustering} implements a divisive hierarchical clustering
#' algorithm that respect size constraints. Starting from any clustering satisfying
#' the size constraints (which may be a clustering with a single cluster containing
#' all data points), the function searches for clusters that can
#' be broken into two or more new clusters without violating the constraints. When
#' such a cluster is found, it breaks the cluster into two new cluster. It continues
#' in this fashion until all remaining clusters are unbreakable.
#'
#' Breakable clusters are broken in three stages. First,
#' it tries to find two data points as far as possible from each other. The
#' two points are called \emph{centers}, and they are the starting points
#' for the new clusters. The remaining
#' data points in the old cluster will be assigned to one of the centers.
#' In the second stage, each center picks the clostest data points so that
#' the new two clusters satisfy the size constraint. In the last stage,
#' data points that still are in the old cluster are assigned to the cluster containing
#' the closest center. The final stage is done either for each point one-by-one
#' or in batches (see below). When all data points are assigned to a cluster, the old
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
#' @param distances a \code{\link[distances]{distances}} object with distances
#'                  between the data points.
#' @param size_constraint an integer with the required minimum cluster size.
#' @param batch_assign a bool indicating whether data points should be assigned in batches when
#'                     spliting clusters (see below for details).
#' @param existing_clustering \code{NULL} or a \code{\link{scclust}} object containing a non-empty
#'                            clustering. If \code{NULL}, the function will start with a
#'                            single cluster containing all data points (i.e., it derives a clustering
#'                            from scratch).
#'
#' @return Returns a \code{\link{scclust}} object with the derived clustering.
#'
#' @seealso \code{\link{sc_clustering}} is the main clustering function in the package.
#'
#' @keywords cluster
#' @export
hierarchical_clustering <- function(distances,
                                    size_constraint,
                                    batch_assign = TRUE,
                                    existing_clustering = NULL) {
  ensure_distances(distances)
  num_data_points <- length(distances)
  size_constraint <- coerce_size_constraint(size_constraint, num_data_points)
  batch_assign <- coerce_scalar_indicator(batch_assign)
  if (!is.null(existing_clustering)) {
    ensure_scclust(existing_clustering, num_data_points)
  }

  clustering <- .Call(Rscc_hierarchical_clustering,
                      distances,
                      size_constraint,
                      batch_assign,
                      existing_clustering)

  make_scclust(clustering$cluster_labels,
               clustering$cluster_count,
               attr(distances, "ids", exact = TRUE))
}
