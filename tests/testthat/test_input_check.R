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

library(scclust)
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

t_ensure_distances <- function(t_distances = distances::distances(matrix(1:10, nrow = 5)),
                               t_req_length = NULL) {
  ensure_distances(t_distances, t_req_length)
}

test_that("`ensure_distances` checks input.", {
  expect_silent(t_ensure_distances())
  expect_silent(t_ensure_distances(t_req_length = 5))
  expect_error(t_ensure_distances(t_distances = "a"),
               regexp = "`t_distances` is not a `scc_distances` object.")
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
# ensure_scclust
# ==============================================================================

t_ensure_scclust <- function(t_clustering = scclust(rep(letters[1:5], 2)),
                             t_req_length = NULL) {
  ensure_scclust(t_clustering, t_req_length)
}

test_that("`ensure_scclust` checks input.", {
  expect_silent(t_ensure_scclust())
  expect_silent(t_ensure_scclust(t_req_length = 10))
  expect_error(t_ensure_scclust(t_clustering = "a"),
               regexp = "`t_clustering` is not a `scclust` object.")
  expect_error(t_ensure_scclust(t_req_length = 4),
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
# coerce_data_point_indices
# ==============================================================================

t_coerce_data_point_indices <- function(t_indices = c(TRUE, FALSE, TRUE, TRUE, FALSE),
                                        t_num_data_points = 5L) {
  coerce_data_point_indices(t_indices, t_num_data_points)
}

test_that("`coerce_data_point_indices` checks input.", {
  expect_silent(t_coerce_data_point_indices())
  expect_silent(t_coerce_data_point_indices(t_indices = NULL))
  expect_silent(t_coerce_data_point_indices(t_indices = c(1L, 4L)))

  expect_error(t_coerce_data_point_indices(t_indices = c(TRUE, FALSE, NA, TRUE, FALSE)),
               regexp = "`t_indices` may not contain NAs.")
  expect_error(t_coerce_data_point_indices(t_indices = rep(FALSE, 5)),
               regexp = "`t_indices` cannot be all `FALSE`.")
  expect_error(t_coerce_data_point_indices(t_indices = c(TRUE, FALSE, FALSE, TRUE)),
               regexp = "`t_indices` is not of length `t_num_data_points`.")
  expect_error(t_coerce_data_point_indices(t_indices = letters[1:5]),
               regexp = "`t_indices` must be integer, logical or NULL.")

  expect_error(t_coerce_data_point_indices(t_indices = c(1L, NA, 4L)),
               regexp = "`t_indices` may not contain NAs.")
  expect_error(t_coerce_data_point_indices(t_indices = c(-1L, 2L, 4L)),
               regexp = "`t_indices` must be positive.")
  expect_error(t_coerce_data_point_indices(t_indices = integer()),
               regexp = "`t_indices` cannot be empty.")
})

test_that("`coerce_data_point_indices` coerces correctly.", {
  expect_equal(t_coerce_data_point_indices(),
               c(TRUE, FALSE, TRUE, TRUE, FALSE))
  expect_equal(t_coerce_data_point_indices(t_indices = NULL),
               NULL)
  expect_equal(t_coerce_data_point_indices(t_indices = c(1L, 4L)),
               c(1L, 4L))
  expect_equal(t_coerce_data_point_indices(t_indices = as.numeric(c(1L, 4L))),
               c(1L, 4L))
})


# ==============================================================================
# coerce_scalar_indicator
# ==============================================================================

t_coerce_scalar_indicator <- function(t_x = TRUE) {
  coerce_scalar_indicator(t_x)
}

test_that("`coerce_scalar_indicator` checks input.", {
  expect_silent(t_coerce_scalar_indicator())
  expect_error(t_coerce_scalar_indicator(t_x = "A"),
               regexp = "`t_x` must be TRUE or FALSE.")
})

test_that("`coerce_scalar_indicator` coerces correctly.", {
  expect_equal(t_coerce_scalar_indicator(), TRUE)
  expect_equal(t_coerce_scalar_indicator(FALSE), FALSE)
  expect_equal(t_coerce_scalar_indicator("TRUE"), TRUE)
  expect_equal(t_coerce_scalar_indicator("T"), TRUE)
  expect_equal(t_coerce_scalar_indicator(NULL), FALSE)
})


# ==============================================================================
# coerce_radius
# ==============================================================================

t_coerce_radius <- function(t_radius = 0.5, t_is_seed = FALSE) {
  coerce_radius(t_radius, t_is_seed)
}

test_that("`coerce_radius` checks input.", {
  expect_silent(t_coerce_radius())
  expect_silent(t_coerce_radius(t_radius = 1L))
  expect_silent(t_coerce_radius(t_radius = "seed_radius"))
  expect_silent(t_coerce_radius(t_radius = NULL))
  expect_error(t_coerce_radius(t_radius = c(1.4, 2.4)),
               regexp = "`t_radius` must be scalar.")
  expect_error(t_coerce_radius(t_radius = as.numeric(NA)),
               regexp = "`t_radius` may not be NA.")
  expect_error(t_coerce_radius(t_radius = "invalid"),
               regexp = "`t_radius` must be one of \"no_radius\", \"seed_radius\", \"estimated_radius\".")
  expect_error(t_coerce_radius(t_radius = "seed_radius", t_is_seed = TRUE),
               regexp = "`t_radius` must be numeric or `NULL`.")
  expect_error(t_coerce_radius(t_radius = -0.5),
               regexp = "`t_radius` must be positive.")
  expect_error(t_coerce_radius(t_radius = TRUE),
               regexp = "`t_radius` must be numeric, character or `NULL`.")
})

test_that("`coerce_radius` coerces correctly.", {
  expect_equal(t_coerce_radius(), 0.5)
  expect_type(t_coerce_radius(t_radius = 1L), "double")
  expect_equal(t_coerce_radius(t_radius = 1L), 1)
  expect_equal(t_coerce_radius(t_radius = "no_radius"), "no_radius")
  expect_equal(t_coerce_radius(t_radius = "seed_radius"), "seed_radius")
  expect_equal(t_coerce_radius(t_radius = "estimated_radius"), "estimated_radius")
  expect_equal(t_coerce_radius(t_radius = "no_rad"), "no_radius")
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
