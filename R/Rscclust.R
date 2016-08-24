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
