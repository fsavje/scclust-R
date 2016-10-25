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


#' Rscclust: R wrapper for the scclust library
#'
#' This package wraps the \code{scclust} C library for R. \code{scclust}
#' provides functions for size constrained clustering.
#'
#' This package is under heavy development, please exercise great caution when using it.
#'
#' More information and the latest version is found here:
#' \url{https://github.com/fsavje/Rscclust}.
#'
#' More information about the underlying \code{scclust} library is found here:
#' \url{https://github.com/fsavje/scclust}.
#'
#' Bug reports and suggestions are greatly appreciated. They
#' are best reported here:
#' \url{https://github.com/fsavje/Rscclust/issues}.
#'
#' @docType package
#' @name Rscclust-package
#' @aliases Rscclust
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("The `Rscclust` package is under heavy development. Please exercise great caution when using it.")
  packageStartupMessage("Bug reports and suggestions are greatly appreciated. They are best reported here: https://github.com/fsavje/Rscclust/issues")
}
