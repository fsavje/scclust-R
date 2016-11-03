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


# ==============================================================================
# Helper functions
# ==============================================================================

# Throw error
new_error <- function(...) {
  stop(structure(list(message = paste0(...),
                      call = match.call(definition = sys.function(-2),
                                        call = sys.call(which = -2),
                                        expand.dots = TRUE,
                                        envir = sys.frame(-3))),
                 class = c("error", "condition")))
}


# Throw warning
new_warning <- function(...) {
  warning(structure(list(message = paste0(...),
                         call = match.call(definition = sys.function(-2),
                                           call = sys.call(which = -2),
                                           expand.dots = TRUE,
                                           envir = sys.frame(-3))),
                    class = c("warning", "condition")))
}


# Is `x` a numeric that can be coerced into integer without loss of information?
is.numeric_integer <- function(x) {
  is.numeric(x) &&
    !any(is.nan(x)) &&
    !any(is.infinite(x)) &&
    all(is.na(x) | as.integer(x) == x)
}


# ==============================================================================
# Ensure functions
# ==============================================================================

# Ensure that `distances` is `Rscc_distances` object
ensure_distances <- function(distances,
                             req_length = NULL) {
  if (!is.Rscc_distances(distances)) {
    new_error("`", match.call()$distances, "` is not a `Rscc_distances` object.")
  }
  if (!is.null(req_length) && (data_point_count_distances(distances) != req_length)) {
    new_error("`", match.call()$distances, "` does not contain `", match.call()$req_length, "` data points.")
  }
}


# Ensure that `indicators` are non-NA logicals
ensure_indicators <- function(indicators,
                              req_length = NULL,
                              any_true = FALSE) {
  if (!is.logical(indicators)) {
    new_error("`", match.call()$indicators, "` must be logical.")
  }
  if (any(is.na(indicators))) {
    new_error("`", match.call()$indicators, "` may not contain NAs.")
  }
  if (!is.null(req_length) && (length(indicators) != req_length)) {
    new_error("`", match.call()$indicators, "` is not of length `", match.call()$req_length, "`.")
  }
  if (any_true) {
    if (!any(indicators)) {
      new_error("`", match.call()$indicators, "` cannot be all `FALSE`.")
    }
  }
}


# Ensure that `clustering` is a `Rscc_clustering` object
ensure_Rscc_clustering <- function(clustering,
                                   req_length = NULL) {
  if (!is.Rscc_clustering(clustering)) {
    new_error("`", match.call()$clustering, "` is not a `Rscc_clustering` object.")
  }
  if (!is.null(req_length) && (data_point_count_clustering(clustering) != req_length)) {
    new_error("`", match.call()$clustering, "` does not contain `", match.call()$req_length, "` data points.")
  }
}


# ==============================================================================
# Coerce functions
# ==============================================================================

# Similar to `match.arg` but with custom error message
coerce_args <- function(arg,
                        choices) {
  stopifnot(is.character(choices),
            length(choices) > 0)
  if (!is.character(arg) || (length(arg) != 1L)) {
    new_error("`", match.call()$arg, "` must be character scalar.")
  }
  i <- pmatch(arg, choices, nomatch = 0L)
  if (i == 0) {
    new_error("`", match.call()$arg, "` must be one of ", paste0(paste0("\"", choices, "\""), collapse = ", "), ".")
  }
  choices[i]
}


# Coerce `x` to character vector
coerce_character <- function(x,
                             req_length = NULL) {
  x <- as.character(x)
  if (!is.null(req_length) && (length(x) != req_length)) {
    new_error("`", match.call()$x, "` is not of length `", match.call()$req_length, "`.")
  }
  x
}


# Coerce `counts` to non-NA, non-negative integers
coerce_counts <- function(counts,
                          req_length = NULL) {
  if (!is.integer(counts)) {
    if (is.numeric_integer(counts)) {
      storage.mode(counts) <- "integer"
    } else {
      new_error("`", match.call()$counts, "` must be integer.")
    }
  }
  if (any(is.na(counts))) {
    new_error("`", match.call()$counts, "` may not contain NAs.")
  }
  if (any(counts < 0L)) {
    new_error("`", match.call()$counts, "` must be non-negative.")
  }
  if (!is.null(req_length) && (length(counts) != req_length)) {
    new_error("`", match.call()$counts, "` is not of length `", match.call()$req_length, "`.")
  }
  counts
}


# Coerce `data` to non-NA, numeric matrix and extract `id_variable`
coerce_distance_data <- function(data,
                                 id_variable,
                                 dist_variables) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    new_error("`", match.call()$data, "` must be matrix or data frame.")
  }
  if (nrow(data) < 2L) {
    new_error("`", match.call()$data, "` must contain at least two data points.")
  }
  if (is.matrix(data)) {
    if (!is.null(dist_variables)) {
      new_error("`", match.call()$dist_variables, "` must be NULL when `", match.call()$data, "` is matrix.")
    }
  } else {
    # is.data.frame(data) == TRUE
    if (!is.null(id_variable) && (length(id_variable) == 1)) {
      if (!(as.character(id_variable) %in% colnames(data))) {
        new_error("`", match.call()$id_variable, "` does not exists as column in `", match.call()$data, "`.")
      }
      if (is.null(dist_variables)) {
        dist_variables <- setdiff(colnames(data), as.character(id_variable))
      }
      id_variable <- data[, as.character(id_variable), drop = TRUE]
    }
    if (!is.null(dist_variables)) {
      if (!all(as.character(dist_variables) %in% colnames(data))) {
        new_error("Some entries in `", match.call()$dist_variables, "` do not exist as columns in `", match.call()$data, "`.")
      }
      data <- data[, as.character(dist_variables), drop = FALSE]
    }
    for (col in names(data)) {
      if (!is.double(data[, col])) {
        if (is.numeric(data[, col])) {
          data[, col] <- as.double(data[, col])
        } else if (is.factor(data[, col])) {
          new_warning("Factor columns in `", match.call()$data, "` are coerced to numeric.")
          data[, col] <- as.double(data[, col])
        } else {
          new_error("Cannot coerce all data columns in `", match.call()$data, "` to numeric.")
        }
      }
    }
    data <- as.matrix(data)
  }
  if (!is.double(data)) {
    if (is.numeric(data)) {
      storage.mode(data) <- "double"
    } else {
      new_error("`", match.call()$data, "` must be numeric.")
    }
  }
  if (any(is.na(data))) {
    new_error("`", match.call()$data, "` may not contain NAs.")
  }
  if (!is.null(id_variable) && (length(id_variable) != nrow(data))) {
    new_error("`", match.call()$id_variable, "` does not match `", match.call()$data, "`.")
  }
  list(data = unname(data),
       id_variable = id_variable)
}


# Coerce `mat` to symmetric, positive-semidefinite, numeric matrix
coerce_norm_matrix <- function(mat,
                               num_cov) {
  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  } else if (is.vector(mat)) {
    mat <- diag(mat)
  }
  mat <- unname(mat)
  if (!is.matrix(mat)) {
    new_error("`", match.call()$mat, "` must be matrix, data.frame or vector.")
  }
  if (!is.numeric(mat)) {
    new_error("`", match.call()$mat, "` must be numeric.")
  }
  if (any(is.na(mat))) {
    new_error("`", match.call()$mat, "` may not contain NAs.")
  }
  if (!isSymmetric(mat)) {
    new_error("`", match.call()$mat, "` must be symmetric.")
  }
  if (ncol(mat) != num_cov) {
    new_error("The dimensions of `", match.call()$mat, "` do not correspond to `", match.call()$num_cov, "`.")
  }
  if (any(eigen(mat, symmetric = TRUE, only.values = TRUE)$values <= 2 * .Machine$double.eps)) {
    new_error("`", match.call()$mat, "` must be positive-semidefinite.")
  }
  mat
}


# Coerce `radius` to NULL or a scalar, positive, non-na, numeric
coerce_radius <- function(radius) {
  if (!is.null(radius)) {
    if (!is.numeric(radius)) {
      new_error("`", match.call()$radius, "` must be numeric or `NULL`.")
    }
    if (length(radius) != 1L) {
      new_error("`", match.call()$radius, "` must be scalar.")
    }
    if (is.na(radius)) {
      new_error("`", match.call()$radius, "` may not be NA.")
    }
    if (radius <= 0.0) {
      new_error("`", match.call()$radius, "` must be positive or `NULL`.")
    }
    # If `radius` is integer
    radius <- as.numeric(radius)
  }
  radius
}


# Coerce `size_constraint` to scalar, non-NA integer
coerce_size_constraint <- function(size_constraint,
                                   num_data_points) {
  if (length(size_constraint) != 1L) {
    new_error("`", match.call()$size_constraint, "` must be scalar.")
  }
  if (!is.integer(size_constraint)) {
    if (is.numeric_integer(size_constraint)) {
      storage.mode(size_constraint) <- "integer"
    } else {
      new_error("`", match.call()$size_constraint, "` must be integer.")
    }
  }
  if (is.na(size_constraint)) {
    new_error("`", match.call()$size_constraint, "` may not be NA.")
  }
  if (size_constraint < 2L) {
    new_error("`", match.call()$size_constraint, "` must be greater or equal to two.")
  }
  if (size_constraint > num_data_points) {
    new_error("`", match.call()$size_constraint, "` may not be great than the number of data points.")
  }
  size_constraint
}


# Coerce `total_size_constraint` to scalar, non-NA integer with default as `sum(type_constraints)`
coerce_total_size_constraint <- function(total_size_constraint,
                                         type_constraints,
                                         num_data_points) {
  sum_type_constraints <- sum(type_constraints)
  stopifnot(is.integer(type_constraints),
            sum_type_constraints >= 0L)

  if (is.null(total_size_constraint)) {
    total_size_constraint <- sum_type_constraints
  }
  if (length(total_size_constraint) != 1L) {
    new_error("`", match.call()$total_size_constraint, "` must be scalar.")
  }
  if (!is.integer(total_size_constraint)) {
    if (is.numeric_integer(total_size_constraint)) {
      storage.mode(total_size_constraint) <- "integer"
    } else {
      new_error("`", match.call()$total_size_constraint, "` must be integer.")
    }
  }
  if (is.na(total_size_constraint)) {
    new_error("`", match.call()$total_size_constraint, "` may not be NA.")
  }
  if (total_size_constraint < 2L) {
    new_error("`", match.call()$total_size_constraint, "` must be greater or equal to two.")
  }
  if (total_size_constraint < sum_type_constraints) {
    new_error("`", match.call()$total_size_constraint, "` must be greater or equal to the sum of the type constraints.")
  }
  if (total_size_constraint > num_data_points) {
    new_error("`", match.call()$total_size_constraint, "` may not be great than the number of data points.")
  }
  total_size_constraint
}


# Coerce `type_constraints` to valid type constraints
coerce_type_constraints <- function(type_constraints) {
  if (is.null(names(type_constraints))) {
    new_error("`", match.call()$type_constraints, "` must be named.")
  }
  if (anyDuplicated(names(type_constraints))) {
    new_error("`", match.call()$type_constraints, "` may not contain duplicate names.")
  }
  if (!is.integer(type_constraints)) {
    if (is.numeric_integer(type_constraints)) {
      storage.mode(type_constraints) <- "integer"
    } else {
      new_error("`", match.call()$type_constraints, "` must be integer.")
    }
  }
  if (any(is.na(type_constraints))) {
    new_error("`", match.call()$type_constraints, "` may not contain NAs.")
  }
  if (any(type_constraints < 0L)) {
    new_error("`", match.call()$type_constraints, "` must be non-negative.")
  }
  type_constraints
}


# Coerce `type_labels` to non-NA factor or non-NA, non-negative integer
coerce_type_labels <- function(type_labels,
                               req_length = NULL) {
  if (!is.factor(type_labels) && !is.integer(type_labels)) {
    if (is.numeric(type_labels)) {
      if (is.numeric_integer(type_labels)) {
        type_labels <- as.integer(type_labels)
      } else {
        new_error("`", match.call()$type_labels, "` must be integer or factor.")
      }
    } else if (is.character(type_labels)) {
      type_labels <- as.factor(type_labels)
    } else {
      new_warning("Coercing `", match.call()$type_labels, "` to factor.")
      type_labels <- as.factor(type_labels)
    }
  }
  if (any(is.na(type_labels))) {
    new_error("`", match.call()$type_labels, "` may not contain NAs.")
  }
  if (is.integer(type_labels) && any(type_labels < 0L)) {
    new_error("`", match.call()$type_labels, "` may not contain negtive entries.")
  }
  if (!is.null(req_length) && (length(type_labels) != req_length)) {
    new_error("`", match.call()$type_labels, "` is not of length `", match.call()$req_length, "`.")
  }
  type_labels
}
