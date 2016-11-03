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

library(Rscclust)
context("Input checking in R code")


# ==============================================================================
# new_error & new_warning
# ==============================================================================

t_new_error <- function(...) {
  temp_func <- function(...) {
    new_error(...)
  }
  temp_func(...)
}

t_new_warning <- function(...) {
  temp_func <- function(...) {
    new_warning(...)
  }
  temp_func(...)
}

test_that("`new_error` & `new_warning` make warnings and errors.", {
  expect_error(t_new_error("This is an error."),
               regexp = "This is an error.")
  expect_error(t_new_error("This is", " also ", "an error."),
               regexp = "This is also an error.")
  expect_warning(t_new_warning("This is a warning."),
                 regexp = "This is a warning.")
  expect_warning(t_new_warning("This is", " also ", "a warning."),
                 regexp = "This is also a warning.")
})


# ==============================================================================
# is.numeric_integer
# ==============================================================================

test_that("`is.numeric_integer` makes correct output.", {
  expect_true(is.numeric_integer(c(1, 2, 3, 4, 5)))
  expect_true(is.numeric_integer(1:5))
  expect_true(is.numeric_integer(c(1, 2, NA, 4, 5)))
  expect_false(is.numeric_integer(c(1, 2, NaN, 4, 5)))
  expect_false(is.numeric_integer(c(1, 2, 3, Inf, 5)))
  expect_false(is.numeric_integer(c(1, 2.5, 3, 4, 5)))
})


# ==============================================================================
# ensure_distances
# ==============================================================================

t_ensure_distances <- function(t_distances = make_distances(matrix(1:10, nrow = 5)),
                               t_req_length = NULL) {
  ensure_distances(t_distances, t_req_length)
}

test_that("`ensure_distances` checks input.", {
  expect_silent(t_ensure_distances())
  expect_silent(t_ensure_distances(t_req_length = 5))
  expect_error(t_ensure_distances(t_distances = "a"),
               regexp = "`t_distances` is not a `Rscc_distances` object.")
  expect_error(t_ensure_distances(t_req_length = 4),
               regexp = "`t_distances` does not contain `t_req_length` data points.")
})


# ==============================================================================
# ensure_indicators
# ==============================================================================

t_ensure_indicators <- function(t_indicators = c(TRUE, FALSE, TRUE, TRUE, FALSE),
                                t_req_length = NULL,
                                t_any_true = FALSE) {
  ensure_indicators(t_indicators, t_req_length, t_any_true)
}

test_that("`ensure_indicators` checks input.", {
  expect_silent(t_ensure_indicators())
  expect_silent(t_ensure_indicators(t_indicators = rep(FALSE, 5)))
  expect_silent(t_ensure_indicators(t_req_length = 5))
  expect_silent(t_ensure_indicators(t_any_true = TRUE))
  expect_silent(t_ensure_indicators(t_req_length = 5, t_any_true = TRUE))
  expect_error(t_ensure_indicators(t_indicators = letters[1:5]),
               regexp = "`t_indicators` must be logical.")
  expect_error(t_ensure_indicators(t_indicators = c(TRUE, FALSE, NA, TRUE, FALSE)),
               regexp = "`t_indicators` may not contain NAs.")
  expect_error(t_ensure_indicators(t_req_length = 4),
               regexp = "`t_indicators` is not of length `t_req_length`.")
  expect_error(t_ensure_indicators(t_indicators = rep(FALSE, 5), t_any_true = TRUE),
               regexp = "`t_indicators` cannot be all `FALSE`.")
})


# ==============================================================================
# ensure_Rscc_clustering
# ==============================================================================

t_ensure_Rscc_clustering <- function(t_clustering = Rscc_clustering(rep(letters[1:5], 2)),
                                     t_req_length = NULL) {
  ensure_Rscc_clustering(t_clustering, t_req_length)
}

test_that("`ensure_Rscc_clustering` checks input.", {
  expect_silent(t_ensure_Rscc_clustering())
  expect_silent(t_ensure_Rscc_clustering(t_req_length = 10))
  expect_error(t_ensure_Rscc_clustering(t_clustering = "a"),
               regexp = "`t_clustering` is not a `Rscc_clustering` object.")
  expect_error(t_ensure_Rscc_clustering(t_req_length = 4),
               regexp = "`t_clustering` does not contain `t_req_length` data points.")
})


# ==============================================================================
# coerce_args
# ==============================================================================

t_coerce_args <- function(t_arg = "abc",
                          t_choices = c("abcdef", "123456", "xzy", "amb")) {
  coerce_args(t_arg, t_choices)
}

test_that("`coerce_args` checks input.", {
  expect_silent(t_coerce_args())
  expect_error(t_coerce_args(t_choices = 1L))
  expect_error(t_coerce_args(t_choices = character()))
  expect_error(t_coerce_args(t_arg = 1L),
               regexp = "`t_arg` must be character scalar.")
  expect_error(t_coerce_args(t_arg = c("a", "z")),
               regexp = "`t_arg` must be character scalar.")
  expect_error(t_coerce_args(t_arg = "nonexist"),
               regexp = "`t_arg` must be one of \"abcdef\", \"123456\", \"xzy\", \"amb\".")
  expect_error(t_coerce_args(t_arg = "a"),
               regexp = "`t_arg` must be one of \"abcdef\", \"123456\", \"xzy\", \"amb\".")
})

test_that("`coerce_args` coerces correctly.", {
  expect_identical(t_coerce_args(), "abcdef")
  expect_identical(t_coerce_args(t_arg = "123456"), "123456")
  expect_identical(t_coerce_args(t_arg = "x"), "xzy")
})


# ==============================================================================
# coerce_character
# ==============================================================================

t_coerce_character <- function(t_x = letters[1:10],
                               t_req_length = NULL) {
  coerce_character(t_x, t_req_length)
}

test_that("`coerce_character` checks input.", {
  expect_silent(t_coerce_character())
  expect_silent(t_coerce_character(t_req_length = 10))
  expect_error(t_coerce_character(t_req_length = 8),
               regexp = "`t_x` is not of length `t_req_length`.")
})

test_that("`coerce_character` coerces correctly.", {
  expect_identical(t_coerce_character(), letters[1:10])
  expect_identical(t_coerce_character(t_x = "123456"), "123456")
  expect_identical(t_coerce_character(t_x = c(1, 2, 3, 4)), c("1", "2", "3", "4"))
})


# ==============================================================================
# coerce_cluster_labels
# ==============================================================================

t_coerce_cluster_labels <- function(t_cluster_labels = 1:10,
                                    t_unassigned_labels = NULL) {
  coerce_cluster_labels(t_cluster_labels, t_unassigned_labels)
}

test_that("`coerce_cluster_labels` checks input.", {
  expect_silent(t_coerce_cluster_labels())
  expect_silent(t_coerce_cluster_labels(t_cluster_labels = factor(1:10)))
  expect_silent(t_coerce_cluster_labels(t_unassigned_labels = 1L))
  expect_silent(t_coerce_cluster_labels(t_cluster_labels = factor(1:10),
                                        t_unassigned_labels = "3"))
  expect_error(t_coerce_cluster_labels(t_cluster_labels = dist(1:10)),
               regexp = "`t_cluster_labels` must be factor or vector.")
  expect_error(t_coerce_cluster_labels(t_unassigned_labels = 20L),
               regexp = "`t_unassigned_labels` contains entries not in `t_cluster_labels`.")
})

test_that("`coerce_cluster_labels` coerces correctly.", {
  expect_identical(t_coerce_cluster_labels(),
                   factor(1:10))
  expect_identical(t_coerce_cluster_labels(t_cluster_labels = factor(1:10)),
                   factor(1:10))
  expect_identical(t_coerce_cluster_labels(t_unassigned_labels = 1L),
                   factor(c(NA, 2:10)))
  expect_identical(t_coerce_cluster_labels(t_cluster_labels = factor(1:10),
                                        t_unassigned_labels = "3"),
                   factor(c(1:2, NA, 4:10)))
})


# ==============================================================================
# coerce_counts
# ==============================================================================

t_coerce_counts <- function(t_counts = 1:10,
                            t_req_length = NULL) {
  coerce_counts(t_counts, t_req_length)
}

test_that("`coerce_counts` checks input.", {
  expect_silent(t_coerce_counts())
  expect_silent(t_coerce_counts(t_counts = as.numeric(1:10)))
  expect_silent(t_coerce_counts(t_req_length = 10))
  expect_error(t_coerce_counts(t_counts = c(1, 1.5, 2)),
               regexp = "`t_counts` must be integer.")
  expect_error(t_coerce_counts(t_counts = c(1, NA, 2)),
               regexp = "`t_counts` may not contain NAs.")
  expect_error(t_coerce_counts(t_counts = c(1, -4, 2)),
               regexp = "`t_counts` must be non-negative.")
  expect_error(t_coerce_counts(t_req_length = 8),
               regexp = "`t_counts` is not of length `t_req_length`.")
})

test_that("`coerce_counts` coerces correctly.", {
  expect_identical(t_coerce_counts(), 1:10)
  expect_identical(t_coerce_counts(t_counts = as.numeric(1:10)), 1:10)
})


# ==============================================================================
# coerce_distance_data
# ==============================================================================

t_distance_data_frame <- data.frame(x = 1:10,
                                    y = 11:20,
                                    z = 21:30,
                                    id = letters[1:10],
                                    stringsAsFactors = FALSE)

t_coerce_distance_data <- function(t_data = matrix(as.numeric(1:10), nrow = 5),
                                   t_id_variable = NULL,
                                   t_dist_variables = NULL) {
  coerce_distance_data(t_data, t_id_variable, t_dist_variables)
}

test_that("`coerce_distance_data` checks input.", {
  expect_silent(t_coerce_distance_data())
  expect_silent(t_coerce_distance_data(t_data = matrix(1:10, nrow = 5)))
  expect_silent(t_coerce_distance_data(t_data = data.frame(matrix(1:10, nrow = 5))))
  expect_silent(t_coerce_distance_data(t_id_variable = letters[1:5]))
  expect_silent(t_coerce_distance_data(t_data = t_distance_data_frame,
                                       t_id_variable = "id"))
  expect_silent(t_coerce_distance_data(t_data = t_distance_data_frame,
                                       t_id_variable = "id",
                                       t_dist_variables = c("x", "y")))
  expect_error(t_coerce_distance_data(t_data = 1:10),
               regexp = "`t_data` must be matrix or data frame.")
  expect_error(t_coerce_distance_data(t_data = matrix(1:2, nrow = 1)),
               regexp = "`t_data` must contain at least two data points.")
  expect_error(t_coerce_distance_data(t_dist_variables = "X1"),
               regexp = "`t_dist_variables` must be NULL when `t_data` is matrix.")
  expect_error(t_coerce_distance_data(t_data = t_distance_data_frame,
                                      t_id_variable = "xxx"),
               regexp = "`t_id_variable` does not exists as column in `t_data`.")
  expect_error(t_coerce_distance_data(t_data = t_distance_data_frame,
                                      t_id_variable = "id",
                                      t_dist_variables = c("x", "nkjsne")),
               regexp = "Some entries in `t_dist_variables` do not exist as columns in `t_data`.")
  expect_warning(t_coerce_distance_data(t_data = data.frame(x = 1:10,
                                                            y = factor(1:10))),
                 regexp = "Factor columns in `t_data` are coerced to numeric.")
  expect_error(t_coerce_distance_data(t_data = data.frame(x = 1:10,
                                                          y = letters[1:10],
                                                          stringsAsFactors = FALSE)),
               regexp = "Cannot coerce all data columns in `t_data` to numeric.")
  expect_error(t_coerce_distance_data(t_data = matrix(letters[1:10], nrow = 5)),
               regexp = "`t_data` must be numeric.")
  expect_error(t_coerce_distance_data(t_data = matrix(c(1:9, NA), nrow = 5)),
               regexp = "`t_data` may not contain NAs.")
  expect_error(t_coerce_distance_data(t_id_variable = 1:4),
               regexp = "`t_id_variable` does not match `t_data`.")
})

t_dist_test1 <- data.frame(x = as.numeric(1:5),
                           y = 6:10)
t_dist_test2 <- data.frame(x = 1:5,
                           y = 6:10,
                           id = letters[1:5],
                           stringsAsFactors = FALSE)
t_dist_test3 <- data.frame(x = 1:5,
                           y = 6:10,
                           z = 20:24,
                           id = letters[1:5],
                           stringsAsFactors = FALSE)
t_dist_test4 <- data.frame(x = factor(letters[1:5]),
                           y = 6:10)
t_dist_test5 <- matrix(1:10, nrow = 5)
dimnames(t_dist_test5) <- list(1:5, c("a", "b"))

t_dist_ref1 <- list(data = matrix(as.numeric(1:10), nrow = 5),
                  id_variable = NULL)
t_dist_ref2 <- list(data = matrix(as.numeric(1:10), nrow = 5),
                  id_variable = letters[1:5])

test_that("`coerce_distance_data` coerces correctly.", {
  expect_equal(t_coerce_distance_data(), t_dist_ref1)
  expect_equal(t_coerce_distance_data(t_data = matrix(1:10, nrow = 5)), t_dist_ref1)
  expect_equal(t_coerce_distance_data(t_id_variable = letters[1:5]), t_dist_ref2)
  expect_equal(t_coerce_distance_data(t_data = t_dist_test1), t_dist_ref1)
  expect_equal(t_coerce_distance_data(t_data = t_dist_test2,
                                      t_id_variable = "id"), t_dist_ref2)
  expect_equal(t_coerce_distance_data(t_data = t_dist_test3,
                                      t_dist_variables = c("x", "y")), t_dist_ref1)
  expect_equal(t_coerce_distance_data(t_data = t_dist_test3,
                                      t_dist_variables = c("x", "y"),
                                      t_id_variable = "id"), t_dist_ref2)
  expect_warning(expect_equal(t_coerce_distance_data(t_data = t_dist_test4), t_dist_ref1))
  expect_equal(t_coerce_distance_data(t_data = t_dist_test5), t_dist_ref1)
})


# ==============================================================================
# coerce_norm_matrix
# ==============================================================================

t_coerce_norm_matrix <- function(t_mat = rep(1, 4),
                                 t_num_cov = 4) {
  coerce_norm_matrix(t_mat, t_num_cov)
}

test_that("`coerce_norm_matrix` checks input.", {
  expect_silent(t_coerce_norm_matrix())
  expect_silent(t_coerce_norm_matrix(t_mat = diag(rep(1, 4))))
  expect_silent(t_coerce_norm_matrix(t_mat = data.frame(diag(rep(1, 4)))))
  expect_silent(t_coerce_norm_matrix(t_mat = matrix(diag(rep(1, 4)), ncol = 4,
                                                    dimnames = list(c(1:4), letters[1:4]))))
  expect_error(t_coerce_norm_matrix(t_mat = dist(1:4)),
               regexp = "`t_mat` must be matrix, data.frame or vector.")
  expect_error(t_coerce_norm_matrix(t_mat = matrix(letters[1:16], ncol = 4)),
               regexp = "`t_mat` must be numeric.")
  expect_error(t_coerce_norm_matrix(t_mat = c(1, NA, 1, 1)),
               regexp = "`t_mat` may not contain NAs.")
  expect_error(t_coerce_norm_matrix(t_mat = matrix(1:16, ncol = 4)),
               regexp = "`t_mat` must be symmetric.")
  expect_error(t_coerce_norm_matrix(t_num_cov = 3),
               regexp = "The dimensions of `t_mat` do not correspond to `t_num_cov`.")
  expect_error(t_coerce_norm_matrix(t_mat = matrix(c(-1, 2, 2, 3), ncol = 2),
                                    t_num_cov = 2),
               regexp = "`t_mat` must be positive-semidefinite.")
})

test_that("`coerce_norm_matrix` coerces correctly.", {
  expect_equal(t_coerce_norm_matrix(),
               diag(rep(1, 4)))
  expect_equal(t_coerce_norm_matrix(t_mat = diag(rep(1, 4))),
               diag(rep(1, 4)))
  expect_equal(t_coerce_norm_matrix(t_mat = data.frame(diag(rep(1, 4)))),
               diag(rep(1, 4)))
  expect_equal(t_coerce_norm_matrix(t_mat = matrix(diag(rep(1, 4)), ncol = 4,
                                                    dimnames = list(c(1:4), letters[1:4]))),
               diag(rep(1, 4)))
})


# ==============================================================================
# coerce_radius
# ==============================================================================

t_coerce_radius <- function(t_radius = 0.5) {
  coerce_radius(t_radius)
}

test_that("`coerce_radius` checks input.", {
  expect_silent(t_coerce_radius())
  expect_silent(t_coerce_radius(t_radius = 1L))
  expect_silent(t_coerce_radius(t_radius = NULL))
  expect_error(t_coerce_radius(t_radius = "a"),
               regexp = "`t_radius` must be numeric or `NULL`.")
  expect_error(t_coerce_radius(t_radius = c(1.4, 2.4)),
               regexp = "`t_radius` must be scalar.")
  expect_error(t_coerce_radius(t_radius = as.numeric(NA)),
               regexp = "`t_radius` may not be NA.")
  expect_error(t_coerce_radius(t_radius = -0.5),
               regexp = "`t_radius` must be positive or `NULL`.")
})

test_that("`coerce_radius` coerces correctly.", {
  expect_equal(t_coerce_radius(), 0.5)
  expect_type(t_coerce_radius(t_radius = 1L), "double")
  expect_equal(t_coerce_radius(t_radius = 1L), 1)
  expect_null(t_coerce_radius(t_radius = NULL))
})


# ==============================================================================
# coerce_size_constraint
# ==============================================================================

t_coerce_size_constraint <- function(t_size_constraint = 2L,
                                     t_num_data_points = 1000L) {
  coerce_size_constraint(t_size_constraint, t_num_data_points)
}

test_that("`coerce_size_constraint` checks input.", {
  expect_silent(t_coerce_size_constraint())
  expect_silent(t_coerce_size_constraint(t_size_constraint = 2.0))
  expect_error(t_coerce_size_constraint(t_size_constraint = 1:2),
               regexp = "`t_size_constraint` must be scalar.")
  expect_error(t_coerce_size_constraint(t_size_constraint = "a"),
               regexp = "`t_size_constraint` must be integer.")
  expect_error(t_coerce_size_constraint(t_size_constraint = as.integer(NA)),
               regexp = "`t_size_constraint` may not be NA.")
  expect_error(t_coerce_size_constraint(t_size_constraint = 1L),
               regexp = "`t_size_constraint` must be greater or equal to two.")
  expect_error(t_coerce_size_constraint(t_size_constraint = 2000L),
               regexp = "`t_size_constraint` may not be great than the number of data points.")
})

test_that("`coerce_size_constraint` coerces correctly.", {
  expect_identical(t_coerce_size_constraint(), 2L)
  expect_identical(t_coerce_size_constraint(t_size_constraint = 2.0), 2L)
})


# ==============================================================================
# coerce_total_size_constraint
# ==============================================================================

t_coerce_total_size_constraint <- function(t_total_size_constraint = 4L,
                                           t_type_constraints = 1:2,
                                           t_num_data_points = 1000L) {
  coerce_total_size_constraint(t_total_size_constraint, t_type_constraints, t_num_data_points)
}

test_that("`coerce_total_size_constraint` checks input.", {
  expect_silent(t_coerce_total_size_constraint())
  expect_silent(t_coerce_total_size_constraint(t_total_size_constraint = 4.0))
  expect_silent(t_coerce_total_size_constraint(t_total_size_constraint = NULL))
  expect_error(t_coerce_total_size_constraint(t_type_constraints = "s"))
  expect_error(t_coerce_total_size_constraint(t_type_constraints = -5:2))
  expect_error(t_coerce_total_size_constraint(t_total_size_constraint = 1:2),
               regexp = "`t_total_size_constraint` must be scalar.")
  expect_error(t_coerce_total_size_constraint(t_total_size_constraint = "a"),
               regexp = "`t_total_size_constraint` must be integer.")
  expect_error(t_coerce_total_size_constraint(t_total_size_constraint = as.integer(NA)),
               regexp = "`t_total_size_constraint` may not be NA.")
  expect_error(t_coerce_total_size_constraint(t_total_size_constraint = 1L),
               regexp = "`t_total_size_constraint` must be greater or equal to two.")
  expect_error(t_coerce_total_size_constraint(t_total_size_constraint = 2L),
               regexp = "`t_total_size_constraint` must be greater or equal to the sum of the type constraints.")
  expect_error(t_coerce_total_size_constraint(t_total_size_constraint = 2000L),
               regexp = "`t_total_size_constraint` may not be great than the number of data points.")
})

test_that("`coerce_total_size_constraint` coerces correctly.", {
  expect_identical(t_coerce_total_size_constraint(), 4L)
  expect_identical(t_coerce_total_size_constraint(t_total_size_constraint = 4.0), 4L)
  expect_identical(t_coerce_total_size_constraint(t_total_size_constraint = NULL), 3L)
})


# ==============================================================================
# coerce_type_constraints
# ==============================================================================

t_coerce_type_constraints <- function(t_type_constraints = c("1" = 2L, "a" = 4L)) {
  coerce_type_constraints(t_type_constraints)
}

test_that("`coerce_type_constraints` checks input.", {
  expect_silent(t_coerce_type_constraints())
  expect_silent(t_coerce_type_constraints(t_type_constraints = c("1" = 2.0, "a" = 4.0)))

  expect_error(t_coerce_type_constraints(t_type_constraints = c(2L, 4L)),
               regexp = "`t_type_constraints` must be named.")
  expect_error(t_coerce_type_constraints(t_type_constraints = c("1" = 2L, "1" = 4L)),
               regexp = "`t_type_constraints` may not contain duplicate names.")
  expect_error(t_coerce_type_constraints(t_type_constraints = c("1" = "x", "a" = "y")),
               regexp = "`t_type_constraints` must be integer.")
  expect_error(t_coerce_type_constraints(t_type_constraints = c("1" = 2L, "a" = NA)),
               regexp = "`t_type_constraints` may not contain NAs.")
  expect_error(t_coerce_type_constraints(t_type_constraints = c("1" = 2L, "a" = -3L)),
               regexp = "`t_type_constraints` must be non-negative.")
})

test_that("`coerce_type_constraints` coerces correctly.", {
  expect_identical(t_coerce_type_constraints(),
                   c("1" = 2L, "a" = 4L))
  expect_identical(t_coerce_type_constraints(t_type_constraints = c("1" = 2.0, "a" = 4.0)),
                   c("1" = 2L, "a" = 4L))
})



# ==============================================================================
# coerce_type_labels
# ==============================================================================


t_coerce_type_labels <- function(t_type_labels = 1:10,
                                 t_req_length = NULL) {
  coerce_type_labels(t_type_labels, t_req_length)
}

test_that("`coerce_type_labels` checks input.", {
  expect_silent(t_coerce_type_labels())
  expect_silent(t_coerce_type_labels(t_type_labels = factor(letters[1:10])))
  expect_silent(t_coerce_type_labels(t_type_labels = as.numeric(1:10)))
  expect_silent(t_coerce_type_labels(t_type_labels = letters[1:10]))
  expect_error(t_coerce_type_labels(t_type_labels = c(1, 5, 3.5, 3, 5)),
                regexp = "`t_type_labels` must be integer or factor.")
  expect_warning(t_coerce_type_labels(t_type_labels = c(TRUE, FALSE, TRUE, FALSE)),
                 regexp = "Coercing `t_type_labels` to factor.")
  expect_error(t_coerce_type_labels(t_type_labels = c(1L, NA, 4L, 3L)),
               regexp = "`t_type_labels` may not contain NAs.")
  expect_error(t_coerce_type_labels(t_type_labels = c(1L, 5L, -4L, 3L)),
               regexp = "`t_type_labels` may not contain negtive entries.")
  expect_error(t_coerce_type_labels(t_req_length = 5L),
               regexp = "`t_type_labels` is not of length `t_req_length`.")
})

test_that("`coerce_type_labels` coerces correctly.", {
  expect_identical(t_coerce_type_labels(),
                   1:10)
  expect_identical(t_coerce_type_labels(t_type_labels = factor(letters[1:10])),
                   factor(letters[1:10]))
  expect_identical(t_coerce_type_labels(t_type_labels = as.numeric(1:10)),
                   1:10)
  expect_identical(t_coerce_type_labels(t_type_labels = letters[1:10]),
                   factor(letters[1:10]))
  expect_warning(expect_identical(t_coerce_type_labels(t_type_labels = c(TRUE, FALSE, TRUE, FALSE)),
                   factor(c(TRUE, FALSE, TRUE, FALSE))))
})
