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
context("Input checking in C code")


# ==============================================================================
# Rscc_hierarchical.c
# ==============================================================================

c_hierarchical_clustering <- function(distance_object = matrix(as.numeric(1:16), ncol = 8),
                                      size_constraint = 2L,
                                      batch_assign = FALSE,
                                      existing_clustering = NULL,
                                      deep_copy = TRUE) {
  .Call("Rscc_hierarchical_clustering",
        distance_object,
        size_constraint,
        batch_assign,
        existing_clustering,
        deep_copy,
        PACKAGE = "Rscclust")
}

temp_existing_clustering1 <- 1:6
attr(temp_existing_clustering1, "cluster_count") <- 2L
temp_existing_clustering2 <- 1:8
attr(temp_existing_clustering2, "cluster_count") <- 0L

test_that("`Rscc_hierarchical_clustering` checks input.", {
  expect_silent(c_hierarchical_clustering())
  expect_error(c_hierarchical_clustering(distance_object = as.numeric(1:16)),
               regexp = "`R_distance_object` is not a valid distance object.")
  expect_error(c_hierarchical_clustering(distance_object = matrix(1:16, ncol = 8)),
               regexp = "`R_distance_object` is not a valid distance object.")
  expect_error(c_hierarchical_clustering(size_constraint = 2.5),
               regexp = "`R_size_constraint` must be integer.")
  expect_error(c_hierarchical_clustering(batch_assign = 1),
               regexp = "`R_batch_assign` must be logical.")
  expect_error(c_hierarchical_clustering(existing_clustering = letters[1:8]),
               regexp = "`R_existing_clustering` is not a valid clustering object.")
  expect_error(c_hierarchical_clustering(existing_clustering = 1:8),
               regexp = "`R_existing_clustering` is not a valid clustering object.")
  expect_error(c_hierarchical_clustering(existing_clustering = temp_existing_clustering1),
               regexp = "`R_existing_clustering` does not match `R_distance_object`.")
  expect_error(c_hierarchical_clustering(existing_clustering = temp_existing_clustering2),
               regexp = "`R_existing_clustering` is empty.")
  expect_error(c_hierarchical_clustering(deep_copy = 1),
               regexp = "`R_deep_copy` must be logical.")
})


# ==============================================================================
# Rscc_nng.c
# ==============================================================================

c_nng_clustering <- function(distance_object = matrix(as.numeric(1:16), ncol = 8),
                             size_constraint = 2L,
                             seed_method = "exclusion_updating",
                             main_unassigned_method = "closest_seed",
                             main_radius = NULL,
                             main_data_points = NULL,
                             secondary_unassigned_method = "ignore",
                             secondary_radius = NULL) {
  .Call("Rscc_nng_clustering",
        distance_object,
        size_constraint,
        seed_method,
        main_unassigned_method,
        main_radius,
        main_data_points,
        secondary_unassigned_method,
        secondary_radius,
        PACKAGE = "Rscclust")
}

test_that("`Rscc_nng_clustering` checks input.", {
  expect_silent(c_nng_clustering())
  expect_error(c_nng_clustering(distance_object = as.numeric(1:16)),
               regexp = "`R_distance_object` is not a valid distance object.")
  expect_error(c_nng_clustering(distance_object = matrix(1:16, ncol = 8)),
               regexp = "`R_distance_object` is not a valid distance object.")
  expect_error(c_nng_clustering(size_constraint = 2.5),
               regexp = "`R_size_constraint` must be integer.")
  expect_error(c_nng_clustering(seed_method = 1L),
               regexp = "`R_seed_method` must be string.")
  expect_error(c_nng_clustering(seed_method = "invalid"),
               regexp = "Not a valid seed method.")
  expect_error(c_nng_clustering(main_unassigned_method = 1L),
               regexp = "`R_main_unassigned_method` must be string.")
  expect_error(c_nng_clustering(main_unassigned_method = "invalid"),
               regexp = "Not a valid main unassigned method.")
  expect_error(c_nng_clustering(main_radius = "a"),
               regexp = "`R_main_radius` must be NULL or double.")
  expect_error(c_nng_clustering(main_data_points = letters[1:8]),
               regexp = "`R_main_data_points` must be NULL or logical.")
  expect_error(c_nng_clustering(main_data_points = rep(TRUE, 6L)),
               regexp = "Invalid `R_main_data_points`.")
  expect_error(c_nng_clustering(secondary_unassigned_method = 1L),
               regexp = "`R_secondary_unassigned_method` must be string.")
  expect_error(c_nng_clustering(secondary_unassigned_method = "invalid"),
               regexp = "Not a valid main unassigned method.")
  expect_error(c_nng_clustering(secondary_radius = "a"),
               regexp = "`R_secondary_radius` must be NULL or double.")
})


c_nng_clustering_batches <- function(distance_object = matrix(as.numeric(1:16), ncol = 8),
                                     size_constraint = 2L,
                                     main_unassigned_method = "by_nng",
                                     main_radius = NULL,
                                     main_data_points = NULL,
                                     batch_size = 100L) {
  .Call("Rscc_nng_clustering_batches",
        distance_object,
        size_constraint,
        main_unassigned_method,
        main_radius,
        main_data_points,
        batch_size,
        PACKAGE = "Rscclust")
}

test_that("`Rscc_nng_clustering_batches` checks input.", {
  expect_silent(c_nng_clustering_batches())
  expect_error(c_nng_clustering_batches(distance_object = as.numeric(1:16)),
               regexp = "`R_distance_object` is not a valid distance object.")
  expect_error(c_nng_clustering_batches(distance_object = matrix(1:16, ncol = 8)),
               regexp = "`R_distance_object` is not a valid distance object.")
  expect_error(c_nng_clustering_batches(size_constraint = 2.5),
               regexp = "`R_size_constraint` must be integer.")
  expect_error(c_nng_clustering_batches(main_unassigned_method = 1L),
               regexp = "`R_main_unassigned_method` must be string.")
  expect_error(c_nng_clustering_batches(main_unassigned_method = "invalid"),
               regexp = "Not a valid main unassigned method.")
  expect_error(c_nng_clustering_batches(main_radius = "a"),
               regexp = "`R_main_radius` must be NULL or double.")
  expect_error(c_nng_clustering_batches(main_data_points = letters[1:8]),
               regexp = "`R_main_data_points` must be NULL or logical.")
  expect_error(c_nng_clustering_batches(main_data_points = rep(TRUE, 6L)),
               regexp = "Invalid `R_main_data_points`.")
  expect_error(c_nng_clustering_batches(batch_size = "a"),
               regexp = "`R_batch_size` must be integer.")
})


c_nng_clustering_types <- function(distance_object = matrix(as.numeric(1:16), ncol = 8),
                                   type_labels = c(1L, 1L, 2L, 1L, 1L, 2L, 2L, 2L),
                                   type_size_constraints = c("0" = 0L, "1" = 1L, "2" = 1L),
                                   total_size_constraint = 2L,
                                   seed_method = "exclusion_updating",
                                   main_unassigned_method = "closest_seed",
                                   main_radius = NULL,
                                   main_data_points = NULL,
                                   secondary_unassigned_method = "ignore",
                                   secondary_radius = NULL) {
  .Call("Rscc_nng_clustering_types",
        distance_object,
        unclass(type_labels),
        type_size_constraints,
        total_size_constraint,
        seed_method,
        main_unassigned_method,
        main_radius,
        main_data_points,
        secondary_unassigned_method,
        secondary_radius,
        PACKAGE = "Rscclust")
}

test_that("`Rscc_nng_clustering_types` checks input.", {
  expect_silent(c_nng_clustering_types())
  expect_error(c_nng_clustering_types(distance_object = as.numeric(1:16)),
               regexp = "`R_distance_object` is not a valid distance object.")
  expect_error(c_nng_clustering_types(distance_object = matrix(1:16, ncol = 8)),
               regexp = "`R_distance_object` is not a valid distance object.")
  expect_error(c_nng_clustering_types(type_labels = letters[1:8]),
               regexp = "`R_type_labels` must be factor or integer.")
  expect_error(c_nng_clustering_types(type_labels = c(1L, 1L, 2L, 1L, 1L, 2L)),
               regexp = "`R_type_labels` does not match `R_distance_object`.")
  expect_error(c_nng_clustering_types(type_size_constraints = c("0", "1", "2")),
               regexp = "`R_type_size_constraints` must be integer.")
  expect_error(c_nng_clustering_types(type_size_constraints = c("0" = 0L, "1" = -1L, "2" = 1L)),
               regexp = "Negative type size constraint.")
  expect_error(c_nng_clustering_types(total_size_constraint = 2.5),
               regexp = "`R_total_size_constraint` must be integer.")
  expect_error(c_nng_clustering_types(seed_method = 1L),
               regexp = "`R_seed_method` must be string.")
  expect_error(c_nng_clustering_types(seed_method = "invalid"),
               regexp = "Not a valid seed method.")
  expect_error(c_nng_clustering_types(main_unassigned_method = 1L),
               regexp = "`R_main_unassigned_method` must be string.")
  expect_error(c_nng_clustering_types(main_unassigned_method = "invalid"),
               regexp = "Not a valid main unassigned method.")
  expect_error(c_nng_clustering_types(main_radius = "a"),
               regexp = "`R_main_radius` must be NULL or double.")
  expect_error(c_nng_clustering_types(main_data_points = letters[1:8]),
               regexp = "`R_main_data_points` must be NULL or logical.")
  expect_error(c_nng_clustering_types(main_data_points = rep(TRUE, 6L)),
               regexp = "Invalid `R_main_data_points`.")
  expect_error(c_nng_clustering_types(secondary_unassigned_method = 1L),
               regexp = "`R_secondary_unassigned_method` must be string.")
  expect_error(c_nng_clustering_types(secondary_unassigned_method = "invalid"),
               regexp = "Not a valid main unassigned method.")
  expect_error(c_nng_clustering_types(secondary_radius = "a"),
               regexp = "`R_secondary_radius` must be NULL or double.")
})


# ==============================================================================
# Rscc_utilities.c
# ==============================================================================

temp_clustering1 <- c(1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L)
attr(temp_clustering1, "cluster_count") <- 2L
temp_clustering2 <- c(1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L)
attr(temp_clustering2, "cluster_count") <- 0L

c_check_clustering <- function(clustering = temp_clustering1,
                               size_constraint = 2L) {
  .Call("Rscc_check_clustering",
        clustering,
        size_constraint,
        PACKAGE = "Rscclust")
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
})


c_check_clustering_types <- function(clustering = temp_clustering1,
                                     type_labels = c(1L, 2L, 1L, 1L, 2L, 1L, 2L, 2L),
                                     type_size_constraints = c("0" = 0L, "1" = 1L, "2" = 1L),
                                     total_size_constraint = 2L) {
  .Call("Rscc_check_clustering_types",
        clustering,
        unclass(type_labels),
        type_size_constraints,
        total_size_constraint,
        PACKAGE = "Rscclust")
}

test_that("`Rscc_check_clustering_types` checks input.", {
  expect_silent(c_check_clustering_types())
  expect_error(c_check_clustering_types(clustering = letters[1:8]),
               regexp = "`R_clustering` is not a valid clustering object.")
  expect_error(c_check_clustering_types(clustering = 1:8),
               regexp = "`R_clustering` is not a valid clustering object.")
  expect_error(c_check_clustering_types(clustering = temp_clustering2),
               regexp = "`R_clustering` is empty.")
  expect_error(c_check_clustering_types(type_labels = letters[1:8]),
               regexp = "`R_type_labels` must be factor or integer.")
  expect_error(c_check_clustering_types(type_labels = c(1L, 1L, 2L, 1L, 1L, 2L)),
               regexp = "`R_type_labels` does not match `R_clustering`.")
  expect_error(c_check_clustering_types(type_size_constraints = c("0", "1", "2")),
               regexp = "`R_type_size_constraints` must be integer.")
  expect_error(c_check_clustering_types(type_size_constraints = c("0" = 0L, "1" = -1L, "2" = 1L)),
               regexp = "Negative type size constraint.")
  expect_error(c_check_clustering_types(total_size_constraint = 2.5),
               regexp = "`R_total_size_constraint` must be integer.")
})


c_get_clustering_stats <- function(clustering = temp_clustering1,
                                   distance_object = matrix(as.numeric(1:16), ncol = 8)) {
  .Call("Rscc_get_clustering_stats",
        clustering,
        distance_object,
        PACKAGE = "Rscclust")
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
               regexp = "`R_distance_object` is not a valid distance object.")
  expect_error(c_get_clustering_stats(distance_object = matrix(1:16, ncol = 8)),
               regexp = "`R_distance_object` is not a valid distance object.")
  expect_error(c_get_clustering_stats(distance_object = matrix(as.numeric(1:14), ncol = 7)),
               regexp = "`R_distance_object` does not match `R_clustering`.")
})


# ==============================================================================
# Check scclust error
# ==============================================================================

test_that("scclust returns errors correctly.", {
  expect_silent(c_hierarchical_clustering())
  expect_error(c_hierarchical_clustering(size_constraint = 1L),
               regexp = "[(]scclust:src/hierarchical_clustering.c")
  expect_error(nng_clustering(make_distances(matrix(c(0.1, 0.2, 0.3), ncol = 1)), size_constraint = 2L, main_radius = 0.001),
               "No clustering satisfying the specified radius constraints exists.")
})
