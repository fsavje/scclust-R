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


#' scclust: Size-Constrained Clustering
#'
#' The \code{scclust} package is an R wrapper for the \code{scclust} library.
#' The package provides functions to construct near-optimal size-constrained
#' clusterings. Subject to user-specified constraints on the size and composition
#' of the clusters, \code{scclust} constructs a clustering so that within-cluster
#' pair-wise distances are minimized.
#'
#' The main clustering function is \code{\link{sc_clustering}}. Statistics about
#' clusters can be derived with the \code{\link{get_clustering_stats}}
#' function. To check if a clustering satisfies some set of
#' constraints, use \code{\link{check_clustering}}. Use \code{\link{scclust}} to
#' construct a \code{scclust} object from an existing clustering.
#'
#' Clusters can also be constructed with \code{\link{hierarchical_clustering}}.
#' However, this function does not support type constraints and does not provide
#' optimality guarantees. Its main use is to refine clusterings constructed with
#' the \code{\link{sc_clustering}} function.
#'
#' \code{scclust} was made with large data sets in mind, and it can cluster tens
#' of millions of data points within minutes on an ordinary desktop computer.
#'
#' See the package's website for more information:
#' \url{https://github.com/fsavje/scclust-R}.
#'
#' More information about the \code{scclust} library is found here:
#' \url{https://github.com/fsavje/scclust}.
#'
#' Bug reports and suggestions are greatly appreciated. They are best reported
#' here: \url{https://github.com/fsavje/scclust-R/issues}.
#'
#' @references
#' Higgins, Michael J., Fredrik Sävje and Jasjeet S. Sekhon (2016),
#' \sQuote{Improving massive experiments with threshold blocking},
#' \emph{Proceedings of the National Academy of Sciences}, \bold{113:27}, 7369--7376.
#'
#' Sävje, Fredrik and Michael J. Higgins and Jasjeet S. Sekhon (2017),
#' \sQuote{Generalized Full Matching}, arXiv 1703.03882.
#' \url{https://arxiv.org/abs/1703.03882}
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
