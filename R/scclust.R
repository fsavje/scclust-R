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


#' scclust: R wrapper for the scclust library
#'
#' This package wraps the \code{scclust} C library for R. \code{scclust}
#' provides functions for size constrained clustering. Besides direct wrapper
#' function, the package provides convenience and utility functions.
#'
#' The clustering functions are \code{\link{make_clustering}} and
#' \code{\link{hierarchical_clustering}}. To get statistics about a clustering
#' use \code{\link{get_clustering_stats}}. To construct a clustering object from
#' an existing clustering see \code{\link{scc_clustering}}. To check if a clustering
#' satisfies some clustering constraints use \code{\link{check_clustering}}.
#'
#' This package is under development, please exercise caution when using it.
#'
#' More information and the latest version is found here:
#' \url{https://github.com/fsavje/scclust-R}.
#'
#' More information about the underlying \code{scclust} library is found here:
#' \url{https://github.com/fsavje/scclust}.
#'
#' Bug reports and suggestions are greatly appreciated. They
#' are best reported here:
#' \url{https://github.com/fsavje/scclust-R/issues}.
#'
#' @docType package
#' @name scclust-package
#' @aliases scclust
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("The `scclust` package is under development. Please exercise caution when using it.")
  packageStartupMessage("Bug reports and suggestions are greatly appreciated. They are best reported here: https://github.com/fsavje/scclust-R/issues")
}

.onLoad <- function(libname, pkgname) {
  set_dist_functions()
}
