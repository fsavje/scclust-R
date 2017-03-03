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

library(scclust)
context("Input checking in C code")


# ==============================================================================
# scc_hierarchical.c
# ==============================================================================

c_hierarchical_clustering <- function(distance_object = distances::distances(matrix(as.numeric(1:16), ncol = 2)),
                                      size_constraint = 2L,
                                      batch_assign = FALSE,
                                      existing_clustering = NULL) {
  .Call(Rscc_hierarchical_clustering,
        distance_object,
        size_constraint,
        batch_assign,
        existing_clustering)
}

temp_existing_clustering1 <- 1:6
attr(temp_existing_clustering1, "cluster_count") <- 2L
temp_existing_clustering2 <- 1:8
attr(temp_existing_clustering2, "cluster_count") <- 0L

test_that("`Rscc_hierarchical_clustering` checks input.", {
  expect_silent(c_hierarchical_clustering())
  expect_error(c_hierarchical_clustering(distance_object = as.numeric(1:16)),
               regexp = "`R_distances` is not a valid distance object.")
  expect_error(c_hierarchical_clustering(distance_object = matrix(1:16, ncol = 8)),
               regexp = "`R_distances` is not a valid distance object.")
  expect_error(c_hierarchical_clustering(size_constraint = 2.5),
               regexp = "`R_size_constraint` must be integer.")
  expect_error(c_hierarchical_clustering(batch_assign = 1),
               regexp = "`R_batch_assign` must be logical.")
  expect_error(c_hierarchical_clustering(existing_clustering = letters[1:8]),
               regexp = "`R_existing_clustering` is not a valid clustering object.")
  expect_error(c_hierarchical_clustering(existing_clustering = 1:8),
               regexp = "`R_existing_clustering` is not a valid clustering object.")
  expect_error(c_hierarchical_clustering(existing_clustering = temp_existing_clustering1),
               regexp = "`R_existing_clustering` does not match `R_distances`.")
  expect_error(c_hierarchical_clustering(existing_clustering = temp_existing_clustering2),
               regexp = "`R_existing_clustering` is empty.")
})


# ==============================================================================
# scc_nng.c
# ==============================================================================

c_make_clustering <- function(distance_object = distances::distances(matrix(as.numeric(1:16), ncol = 2)),
                              size_constraint = 2L,
                              type_labels = c(1L, 1L, 2L, 1L, 1L, 2L, 2L, 2L),
                              type_constraints = c("0" = 0L, "1" = 1L, "2" = 1L),
                              seed_method = "exclusion_updating",
                              primary_data_points = NULL,
                              primary_unassigned_method = "closest_seed",
                              secondary_unassigned_method = "ignore",
                              seed_radius = NULL,
                              primary_radius = NULL,
                              secondary_radius = NULL,
                              batch_size = NULL) {
  .Call(Rscc_make_clustering,
        distance_object,
        size_constraint,
        type_labels,
        type_constraints,
        seed_method,
        primary_data_points,
        primary_unassigned_method,
        secondary_unassigned_method,
        seed_radius,
        primary_radius,
        secondary_radius,
        batch_size)
}

test_that("`Rscc_make_clustering` checks input.", {
  expect_silent(c_make_clustering())
  expect_error(c_make_clustering(distance_object = as.numeric(1:16)),
               regexp = "`R_distances` is not a valid distance object.")
  expect_error(c_make_clustering(distance_object = matrix(1:16, ncol = 8)),
               regexp = "`R_distances` is not a valid distance object.")
  expect_error(c_make_clustering(size_constraint = 2.5),
               regexp = "`R_size_constraint` must be integer.")
  expect_error(c_make_clustering(type_labels = letters[1:8]),
               regexp = "`R_type_labels` must be factor, integer or NULL.")
  expect_error(c_make_clustering(type_labels = NULL),
               regexp = "`R_type_constraints` must be NULL when no types are supplied.")
  expect_error(c_make_clustering(type_labels = c(1L, 1L, 2L, 1L, 1L, 2L)),
               regexp = "`R_type_labels` does not match `R_distances`.")
  expect_error(c_make_clustering(type_constraints = c("0", "1", "2")),
               regexp = "`R_type_constraints` must be integer.")
  expect_error(c_make_clustering(type_constraints = c("0" = 0L, "1" = -1L, "2" = 1L)),
               regexp = "Negative type size constraint.")
  expect_error(c_make_clustering(seed_method = 1L),
               regexp = "`R_seed_method` must be string.")
  expect_error(c_make_clustering(seed_method = "invalid"),
               regexp = "Not a valid seed method.")
  expect_error(c_make_clustering(primary_data_points = letters[1:8]),
               regexp = "`R_primary_data_points` must be NULL or integer.")
  expect_error(c_make_clustering(primary_unassigned_method = 1L),
               regexp = "`R_primary_unassigned_method` must be string.")
  expect_error(c_make_clustering(primary_unassigned_method = "invalid"),
               regexp = "Not a valid unassigned method.")
  expect_error(c_make_clustering(secondary_unassigned_method = 1L),
               regexp = "`R_secondary_unassigned_method` must be string.")
  expect_error(c_make_clustering(secondary_unassigned_method = "invalid"),
               regexp = "Not a valid unassigned method.")
  expect_error(c_make_clustering(seed_radius = "a"),
               regexp = "`R_seed_radius` must be NULL or double.")
  expect_error(c_make_clustering(primary_radius = FALSE),
               regexp = "`R_primary_radius` must be NULL, string or double.")
  expect_error(c_make_clustering(primary_radius = "invalid"),
               regexp = "Not a valid radius method.")
  expect_error(c_make_clustering(secondary_radius = FALSE),
               regexp = "`R_secondary_radius` must be NULL, string or double.")
  expect_error(c_make_clustering(secondary_radius = "invalid"),
               regexp = "Not a valid radius method.")
  expect_error(c_make_clustering(batch_size = "invalid"),
               regexp = "`R_batch_size` must be NULL or integer.")
})


# ==============================================================================
# scc_utilities.c
# ==============================================================================

temp_clustering1 <- c(1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L)
attr(temp_clustering1, "cluster_count") <- 2L
temp_clustering2 <- c(1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L)
attr(temp_clustering2, "cluster_count") <- 0L

c_check_clustering <- function(clustering = temp_clustering1,
                               size_constraint = 2L,
                               type_labels = c(1L, 2L, 1L, 1L, 2L, 1L, 2L, 2L),
                               type_constraints = c("0" = 0L, "1" = 1L, "2" = 1L)) {
  .Call(Rscc_check_clustering,
        clustering,
        size_constraint,
        unclass(type_labels),
        type_constraints)
}

test_that("`Rscc_check_clustering` checks input.", {
  expect_silent(c_check_clustering())
  expect_error(c_check_clustering(clustering = letters[1:8]),
               regexp = "`R_clustering` is not a valid clustering object.")
  expect_error(c_check_clustering(clustering = 1:8),
               regexp = "`R_clustering` is not a valid clustering object.")
  expect_error(c_check_clustering(clustering = temp_clustering2),
               regexp = "`R_clustering` is empty.")
  expect_error(c_check_clustering(size_constraint = 2.5),
               regexp = "`R_size_constraint` must be integer.")
  expect_error(c_check_clustering(type_labels = letters[1:8]),
               regexp = "`R_type_labels` must be factor, integer or NULL.")
  expect_error(c_make_clustering(type_labels = NULL),
               regexp = "`R_type_constraints` must be NULL when no types are supplied.")
  expect_error(c_check_clustering(type_labels = c(1L, 1L, 2L, 1L, 1L, 2L)),
               regexp = "`R_type_labels` does not match `R_clustering`.")
  expect_error(c_check_clustering(type_constraints = c("0", "1", "2")),
               regexp = "`R_type_constraints` must be integer.")
  expect_error(c_check_clustering(type_constraints = c("0" = 0L, "1" = -1L, "2" = 1L)),
               regexp = "Negative type size constraint.")
})


c_get_clustering_stats <- function(clustering = temp_clustering1,
                                   distance_object = distances::distances(matrix(as.numeric(1:16), ncol = 2))) {
  .Call(Rscc_get_clustering_stats,
        clustering,
        distance_object)
}

test_that("`Rscc_get_clustering_stats` checks input.", {
  expect_silent(c_get_clustering_stats())
  expect_error(c_get_clustering_stats(clustering = letters[1:8]),
               regexp = "`R_clustering` is not a valid clustering object.")
  expect_error(c_get_clustering_stats(clustering = 1:8),
               regexp = "`R_clustering` is not a valid clustering object.")
  expect_error(c_get_clustering_stats(clustering = temp_clustering2),
               regexp = "`R_clustering` is empty.")
  expect_error(c_get_clustering_stats(distance_object = as.numeric(1:16)),
               regexp = "`R_distances` is not a valid distance object.")
  expect_error(c_get_clustering_stats(distance_object = matrix(1:16, ncol = 8)),
               regexp = "`R_distances` is not a valid distance object.")
  expect_error(c_get_clustering_stats(distance_object = matrix(as.numeric(1:14), ncol = 7)),
               regexp = "`R_distances` does not match `R_clustering`.")
})


# ==============================================================================
# Check scclust error
# ==============================================================================

test_that("scclust returns errors correctly.", {
  expect_silent(c_hierarchical_clustering())
  expect_error(c_hierarchical_clustering(size_constraint = 1L),
               regexp = "[(]scclust:src/hierarchical_clustering.c")
  expect_error(make_clustering(distances::distances(matrix(c(0.1, 0.2, 0.3), ncol = 1)), size_constraint = 2L, seed_radius = 0.001),
               "Infeasible radius constraint.")
})
