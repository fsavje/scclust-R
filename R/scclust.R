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


#' Constructor for scclust objects
#'
#' The \code{scclust} function constructs a \code{scclust} object
#' from existing cluster labels.
#'
#' \code{scclust} does not derive clusterings from sets of data points; see
#' \code{\link{sc_clustering}} and \code{\link{hierarchical_clustering}} for
#' that functionality.
#'
#' @param cluster_labels a vector containing each data point's cluster label.
#' @param unassigned_labels labels that denote unassigned data points. If \code{NULL},
#'                          \code{NA} values in \code{cluster_labels} are used to
#'                          denote unassigned points.
#' @param ids IDs of the data points. Should be a vector of the same length as
#'            \code{cluster_labels} or \code{NULL}. If \code{NULL}, the IDs are set to
#'            \code{1:length(cluster_labels)}.
#'
#' @return Returns a \code{scclust} object with the clustering described by
#'         the provided labels.
#'
#' @examples
#' # 10 data points in 3 clusters
#' my_labels1 <- c("A", "A", "B", "C", "B",
#'                 "C", "C", "A", "B", "B")
#' my_scclust1 <- scclust(my_labels1)
#'
#' # 8 data points in 3 clusters, 2 points unassigned, integer labels
#' my_labels2 <- c(1, 1, 2, 3, 2, NA, 3, 1, NA, 2)
#' my_scclust2 <- scclust(my_labels2)
#'
#' # Custom labels indicating unassiged points
#' my_labels3 <- c("A", "A", "B", "C", "NONE",
#'                 "C", "C", "NONE", "B", "B")
#' my_scclust3 <- scclust(my_labels3, unassigned_labels = "NONE")
#'
#' # Two different labels indicating unassiged points
#' my_labels4 <- c("A", "A", "B", "C", "NONE",
#'                 "C", "C", "0", "B", "B")
#' my_scclust4 <- scclust(my_labels4, unassigned_labels = c("NONE", "0"))
#'
#' # Custom data point IDs
#' my_ids <- letters[1:10]
#' my_labels1_ids <- scclust(my_labels1, ids = my_ids)
#'
#' @export
scclust <- function(cluster_labels,
                    unassigned_labels = NULL,
                    ids = NULL) {
  cluster_labels <- coerce_cluster_labels(cluster_labels, unassigned_labels)
  if (!is.null(ids)) {
    ids <- coerce_character(ids, length(cluster_labels))
  }

  make_scclust(as.integer(cluster_labels) - 1L,
               nlevels(cluster_labels),
               ids)
}


#' Check \code{scclust} object
#'
#' \code{is.scclust} checks whether the provided object
#' is a valid instance of the \code{\link{scclust}} class.
#' It does not check whether the clustering itself is sensible or
#' whether the clustering satisfies some set of constraints. See
#' \code{\link{check_clustering}} for that functionality.
#'
#'
#' @param x  object to check.
#'
#' @return Returns \code{TRUE} if \code{x} is a valid
#'         \code{\link{scclust}} object, otherwise \code{FALSE}.
#'
#' @export
is.scclust <- function(x) {
  is.integer(x) &&
    inherits(x, "scclust") &&
    is.integer(attr(x, "cluster_count", exact = TRUE)) &&
    length(attr(x, "cluster_count", exact = TRUE)) == 1 &&
    attr(x, "cluster_count", exact = TRUE) >= 0L
}


#' Count the number of clusters in a clustering
#'
#' \code{cluster_count} returns the number of clusters in a clustering.
#'
#' @param clustering a \code{\link{scclust}} object containing a non-empty clustering.
#'
#' @return Returns an integer with the number of clusters in \code{clustering}.
#'
#' @examples
#' # Example scclust clustering
#' my_scclust <- scclust(c("A", "A", "B", "C", "B",
#'                         "C", "C", "A", "B", "B"))
#'
#' cluster_count(my_scclust)
#' # > 3
#'
#' @export
cluster_count <- function(clustering) {
  ensure_scclust(clustering)
  attr(clustering, "cluster_count", exact = TRUE)
}


#' @export
as.data.frame.scclust <- function(x,
                                  row.names = NULL,
                                  ...) {
  stopifnot(is.scclust(x))
  ids <- attr(x, "ids", exact = TRUE)
  if (is.null(ids)) ids <- 1:length(x)
  data.frame(id = ids,
             cluster_label = as.integer(x),
             row.names = row.names,
             ...)
}


#' @export
print.scclust <- function(x,
                          ...) {
  stopifnot(is.scclust(x))
  ids <- attr(x, "ids", exact = TRUE)
  if (is.null(ids)) ids <- as.character(1:length(x))
  stopifnot(is.character(ids))
  xx <- as.integer(x)
  names(xx) <- ids
  print(xx)
}
