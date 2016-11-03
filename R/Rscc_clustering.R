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


#' Constructor for Rscclust clustering objects.
#'
#' \code{Rscc_clustering} constructs a Rscclust clustering object from existing cluster
#' labels. It is useful when using Rscclust's clustering functions to refine an
#' existing clustering.
#'
#' @param cluster_labels a vector containing each data point's cluster label.
#' @param unassigned_labels a vector containing cluster labels that denote unassigned data points.
#' @param ids optional IDs of the data points. Should be a vector of the same length as
#'            \code{cluster_labels} or \code{NULL}. If \code{NULL}, the IDs are set to
#'            \code{1:length(cluster_labels)}.
#'
#' @return Returns a Rscclust cluster object to be used when calling the
#'         clustering functions.
#'
#' @seealso
#' Functions that accept the Rscclust cluster object as input are
#' \code{\link{hierarchical_clustering}}, \code{\link{check_clustering}},
#' \code{\link{check_clustering_types}} and \code{\link{get_clustering_stats}} for
#'
#' \code{Rscc_clustering} does not derive clusterings from data sets; for that
#' functionality see \code{\link{nng_clustering}}, \code{\link{nng_clustering_batches}},
#' \code{\link{nng_clustering_types}} and \code{\link{hierarchical_clustering}}.
#'
#' @examples
#' # 10 data points in 3 clusters
#' my_labels1 <- c("A", "A", "B", "C", "B", "C", "C", "A", "B", "B")
#' my_clust_obj1 <- Rscc_clustering(my_labels1)
#'
#' # 10 data points in 3 clusters, integer labels
#' my_labels2 <- c(1, 1, 2, 3, 2, 3, 3, 1, 2, 2)
#' my_clust_obj2 <- Rscc_clustering(my_labels2)
#'
#' # 8 data points in 3 clusters, and 2 points unassigned
#' my_labels3 <- c("A", "A", "B", "C", "NONE", "C", "C", "NONE", "B", "B")
#' my_clust_obj3 <- Rscc_clustering(my_labels3, "NONE")
#'
#' # Different unassiged labels
#' my_labels4 <- c("A", "A", "B", "C", "NONE", "C", "C", "0", "B", "B")
#' my_clust_obj4 <- Rscc_clustering(my_labels4, c("NONE", "0"))
#'
#' # Custom data point IDs
#' my_labels5 <- c("A", "A", "B", "C", "B", "C", "C", "A", "B", "B")
#' my_ids <- letters[1:10]
#' my_clust_obj5 <- Rscc_clustering(my_labels5, ids = my_ids)
#'
#' @export
Rscc_clustering <- function(cluster_labels,
                            unassigned_labels = NA,
                            ids = NULL) {
  cluster_labels <- as.vector(cluster_labels)
  cluster_labels[cluster_labels %in% as.vector(unassigned_labels)] <- NA
  cluster_labels <- factor(cluster_labels)

  if (!is.null(ids)) {
    ids <- coerce_character(ids, length(cluster_labels))
  }

  make_Rscc_clustering(as.integer(cluster_labels) - 1L,
                       nlevels(cluster_labels),
                       ids)
}


#' Check \code{Rscc_clustering} object
#'
#' \code{is.Rscc_clustering} checks whether the provided object
#' is a valid instance of the \code{Rscc_clustering} class. It does
#' not whether the clustering itself is sensible.
#'
#' @param obj  object to check.
#'
#' @return Returns \code{TRUE} if \code{obj} is a valid
#'         \code{Rscc_clustering} object, otherwise \code{FALSE}.
#'
#' @export
is.Rscc_clustering <- function(obj) {
  is.integer(obj) &&
    inherits(obj, "Rscc_clustering") &&
    is.integer(attr(obj, "cluster_count", exact = TRUE)) &&
    length(attr(obj, "cluster_count", exact = TRUE)) == 1 &&
    attr(obj, "cluster_count", exact = TRUE) >= 0L
}


#' @export
cluster_count <- function(clustering) {
  ensure_Rscc_clustering(clustering)
  attr(clustering, "cluster_count", exact = TRUE)
}


#' @export
data_point_count.Rscc_clustering <- function(x) {
  stopifnot(is.integer(x))
  length(x)
}


#' @export
as.data.frame.Rscc_clustering <- function(x,
                                          row.names = NULL,
                                          ...) {
  stopifnot(is.Rscc_clustering(x))
  ids <- attr(x, "ids", exact = TRUE)
  if (is.null(ids)) ids <- 1:length(x)
  data.frame(id = ids,
             cluster_label = as.integer(x),
             row.names = row.names,
             ...)
}


#' @export
print.Rscc_clustering <- function(x,
                                  ...) {
  stopifnot(is.Rscc_clustering(x))
  ids <- attr(x, "ids", exact = TRUE)
  if (is.null(ids)) ids <- as.character(1:length(x))
  stopifnot(is.character(ids))
  xx <- as.integer(x)
  names(xx) <- ids
  print(xx)
}
