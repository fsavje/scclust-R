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


#' scclust: Size-constrained clustering
#'
#' This package is an R wrapper for the \code{scclust} library. The library
#' provides functions to construct near-optimal size-constrained clusterings.
#' Subject to user-specified constraints on the minimum size and composition
#' of the clusters, \code{scclust} derives a partition of a set of data points
#' so that the dissimilarity of points assigned to the same cluster is minimized.
#'
#' The main clustering function is \code{\link{sc_clustering}}. Statistics about
#' clusterings can be derived with the \code{\link{get_clustering_stats}} function.
#' To check if a clustering satisfies some set of clustering constraints use
#' \code{\link{check_clustering}}. To construct a \code{scclust} clustering object from
#' an existing clustering use \code{\link{scclust}}.
#'
#' Clusterings can also be constructed with \code{\link{hierarchical_clustering}}.
#' However, this function only support a subset of constraint options and
#' does not provide optimality guarantees. Its main use is to refine clusterings
#' produced by the \code{\link{sc_clustering}} function.
#'
#' See the package's website for more information:
#' \url{https://github.com/fsavje/scclust-R}.
#'
#' More information about the \code{scclust} library is found here:
#' \url{https://github.com/fsavje/scclust}.
#'
#' Bug reports and suggestions are greatly appreciated. They
#' are best reported here:
#' \url{https://github.com/fsavje/scclust-R/issues}.
#'
#' @docType package
#' @name scclust-package
#'
#' @import distances
NULL

#' @useDynLib scclust, .registration = TRUE
.onUnload <- function (libpath) {
  library.dynam.unload("scclust", libpath)
}
