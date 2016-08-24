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


#' Construct distance metric object.
#'
#' \code{make_distances} constructs Euclidean distances to be used in clustering
#' problems (including blocking and matching problems). It optionally normalize and
#' weights the inputted data so to achieve, e.g., Mahalanobis distances or
#' normalized Euclidean distances.
#'
#' Let \eqn{x} and \eqn{y} be two data points in \code{data} described by two vectors. \code{make_distances}
#' uses the following metric to derive the distance between \eqn{x} and \eqn{y}:
#'
#' \deqn{\sqrt{(x - y) N^{-0.5} W N^{-0.5} (x - y)}}{\sqrt((x - y) * N^-0.5 * W * N^-0.5 * (x - y))}
#'
#' where \eqn{N^{-0.5}}{N^-0.5} is the Cholesky decomposition of the inverse of the matrix speficied by \code{normalize}, and \eqn{W}
#' is matrix speficied by \code{weights}.
#'
#' When \code{normalize} is \code{var(data)} (i.e., using the \code{"mahalanobize"} option), this
#' derives the (weighted) Mahalanobis distances. When \code{normalize} is \code{diag(var(data))} (i.e., using
#' the \code{"studentize"} option), this will effective divide each column by its variance leading to (weighted) normalized
#' Euclidean distances. If \code{normalize} is the identity matrix (i.e., using the \code{"none"} or \code{NULL} option), this
#' derives the (weighted) Euclidean distances.
#'
#' @param data a matrix or data frame containing the data points between distances should be derived.
#' @param id_variable optional IDs of the data points.
#'                    If \code{id_variable} is a single string and \code{data} is a data frame, the
#'                    corresponding column in \code{data} will be taken as IDs. That column will be
#'                    excluded from \code{data} when constructing distances. If \code{id_variable} is
#'                    \code{NULL}, the IDs are set to \code{1:nrow(data)}. Otherwise, \code{id_variable}
#'                    must be of length \code{nrow(data)} and will be used directly as IDs.
#' @param dist_variables optional names of the columns in \code{data} that should
#'                       be used when constructing distances. If \code{dist_variables} is \code{NULL},
#'                       all columns will be used (net of eventual column specified by \code{id_variable}).
#'                       If \code{data} is a matrix, \code{dist_variables} must be \code{NULL}.
#' @param normalize optional normalization of the data prior to distance construction. If \code{normalize}
#'                  is \code{NULL} or \code{"none"}, no normalization will be done (effectively setting \code{normalize}
#'                  to the identity matrix). If \code{normalize} is \code{"mahalanobize"}, normalization will be
#'                  done with \code{var(data)} (i.e., resulting in Mahalanobis distances). If \code{normalize} is
#'                  \code{"studentize"}, normalization is done with the diagonal of \code{var(data)}. If \code{normalize}
#'                  is a matrix, it will be used in the normalization. If \code{normalize} is a vector, a diagonal matrix
#'                  with the supplied vector as its diagonal will be used. The matrix used for normalization must be
#'                  positive-semidefinite.
#' @param weights optional weighting of the data prior to distance construction. If \code{normalize} is \code{NULL}
#'                no weighting will be done (effectively setting \code{weights} to the identity matrix). If \code{weights}
#'                is a matrix, that will be used in the weighting. If \code{normalize} is a vector, a diagonal matrix
#'                with the supplied vector as its diagonal will be used. The matrix used for weighting must be
#'                  positive-semidefinite.
#'
#' @return \code{make_distances} returns a distance object to be used when calling the clustering, blocking or matching
#'         functions.
#'
#' @export
make_distances <- function(data,
                           id_variable = NULL,
                           dist_variables = NULL,
                           normalize = NULL,
                           weights = NULL) {

   stopifnot((is.data.frame(data) || is.matrix(data)),
             (is.data.frame(data) || is.null(dist_variables)),
             (nrow(data) >= 2))

   num_data_points <- nrow(data)

   if (is.data.frame(data)) {
      if (!is.null(id_variable) && (length(id_variable) == 1)) {
         stopifnot(is.character(id_variable),
                   (id_variable %in% colnames(data)))
         if (is.null(dist_variables)) dist_variables <- colnames(data)
         dist_variables <- setdiff(dist_variables, id_variable)
         id_variable <- data[, id_variable, drop = TRUE]
      }
      if (!is.null(dist_variables)) {
         stopifnot(is.character(dist_variables),
                   (dist_variables %in% colnames(data)))
         data <- data[, dist_variables, drop = FALSE]
      }
      data <- as.matrix(data)
   }

   stopifnot(is.matrix(data),
             is.numeric(data),
             all(!is.na(data)),
             (is.null(id_variable) ||
                 (length(id_variable) == num_data_points)))

   if (is.character(normalize)) {
      if (normalize == "mahalanobis") normalize <- "mahalanobize"
      normalize <- match.arg(normalize, c("none",
                                          "mahalanobize",
                                          "studentize"))
      normalize <- switch(normalize,
                          none = NULL,
                          mahalanobize = var(data),
                          studentize = diag(diag(var(data))))
   }

   if (!is.null(normalize)) {
      if (is.data.frame(normalize)) normalize <- as.matrix(normalize)
      if (is.vector(normalize)) normalize <- diag(normalize)

      stopifnot(is.matrix(normalize),
                is.numeric(normalize),
                all(!is.na(normalize)),
                isSymmetric(normalize),
                ncol(normalize) == ncol(data))
      if (any(eigen(normalize, symmetric = TRUE, only.values = TRUE)$values <= 2 * .Machine$double.eps)) {
         stop("`normalize` must be positive-semidefinite.")
      }

      data <- tcrossprod(data, chol(solve(normalize)))
   } else {
      normalize <- diag(ncol(data))
   }

   if (!is.null(weights)) {
      if (is.data.frame(weights)) weights <- as.matrix(weights)
      if (is.vector(weights)) weights <- diag(weights)
      stopifnot(is.matrix(weights),
                is.numeric(weights),
                all(!is.na(weights)),
                isSymmetric(weights),
                ncol(weights) == ncol(data))
      if (any(eigen(weights, symmetric = TRUE, only.values = TRUE)$values <= 2 * .Machine$double.eps)) {
         stop("`weights` must be positive-semidefinite.")
      }

      data <- tcrossprod(data, chol(weights))
   } else {
      weights <- diag(ncol(data))
   }

   structure(t(data),
             ids = id_variable,
             normalization = normalize,
             weights = weights,
             class = c("Rscc_distances"))
}


#' @export
print.Rscc_distances <- function(x) {
   cat("This is a Rscc_distances object. It's only purpose is to be used in the Rscclust package.",
       "\nUse `as.matrix` to generate an R matrix with the distances (a very slow operation for large distance matrices).")
}


#' @export
as.matrix.Rscc_distances <- function(x) {
   as.matrix(dist(t(unclass(x))))
}
