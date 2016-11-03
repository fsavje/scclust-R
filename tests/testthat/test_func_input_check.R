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
context("Input checking in external functions")


# ==============================================================================
# Shared objects
# ==============================================================================

sound_distance_object <- make_distances(matrix(c(1, 4, 3, 2, 45, 6, 3, 2, 6, 5,
                                                 34, 2, 4, 6, 4, 6, 4, 2, 7, 8), nrow = 10))
unsound_distance_object <- letters[1:10]
sound_size_constraint <- 2L
unsound_size_constraint <- 0L
sound_clustering <- Rscc_clustering(c(rep("a", 5), rep("b", 5)))
unsound_clustering <- c(rep("a", 5), rep("b", 5))
sound_bool <- TRUE
unsound_bool <- "a"
sound_seed_method <- "lexical"
unsound_seed_method <- "unknown"
sound_unassigned_method <- "ignore"
unsound_unassigned_method <- "unknown"
sound_radius <- NULL
unsound_radius <- "a"
sound_main_data_points <- rep(TRUE, 10)
unsound_main_data_points <- rep(TRUE, 5)
sound_type_labels <- rep(c(1L, 2L), 5)
unsound_type_labels <- rep(c(1L, 2L), 4)
sound_type_size_constraints <- c("1" = 1L, "2" = 1L)
unsound_type_size_constraints <- c(1L, 1L)
sound_total_size_constraint <- 3L
unsound_total_size_constraint <- 100L


# ==============================================================================
# make_distances
# ==============================================================================

sound_data <- matrix(c(1, 4, 3, 2, 45, 6, 3, 2, 6, 5), nrow = 5)
unsound_data <- matrix(letters[1:10], nrow = 5)
sound_id_variable <- letters[1:5]
unsound_id_variable <- letters[1:3]
sound_dist_variables <- NULL
unsound_dist_variables <- dist(1:10)
sound_normalize <- "mahalanobize"
unsound_normalize <- matrix(c(-1, 2, 2, 3), ncol = 2)
sound_weights <- NULL
unsound_weights <- matrix(letters[1:4], ncol = 2)

test_that("`make_distances` checks input.", {
  expect_silent(make_distances(data = sound_data,
                               id_variable = sound_id_variable,
                               dist_variables = sound_dist_variables,
                               normalize = sound_normalize,
                               weights = sound_weights))
  expect_error(make_distances(data = unsound_data,
                              id_variable = sound_id_variable,
                              dist_variables = sound_dist_variables,
                              normalize = sound_normalize,
                              weights = sound_weights))
  expect_error(make_distances(data = sound_data,
                              id_variable = unsound_id_variable,
                              dist_variables = sound_dist_variables,
                              normalize = sound_normalize,
                              weights = sound_weights))
  expect_error(make_distances(data = sound_data,
                              id_variable = sound_id_variable,
                              dist_variables = unsound_dist_variables,
                              normalize = sound_normalize,
                              weights = sound_weights))
  expect_error(make_distances(data = sound_data,
                              id_variable = sound_id_variable,
                              dist_variables = sound_dist_variables,
                              normalize = unsound_normalize,
                              weights = sound_weights))
  expect_error(make_distances(data = sound_data,
                              id_variable = sound_id_variable,
                              dist_variables = sound_dist_variables,
                              normalize = sound_normalize,
                              weights = unsound_weights))
})


# ==============================================================================
# hierarchical_clustering
# ==============================================================================

test_that("`hierarchical_clustering` checks input.", {
  expect_silent(hierarchical_clustering(distance_object = sound_distance_object,
                                        size_constraint = sound_size_constraint,
                                        batch_assign = sound_bool,
                                        existing_clustering = sound_clustering))
  expect_error(hierarchical_clustering(distance_object = unsound_distance_object,
                                       size_constraint = sound_size_constraint,
                                       batch_assign = sound_bool,
                                       existing_clustering = sound_clustering))
  expect_error(hierarchical_clustering(distance_object = sound_distance_object,
                                       size_constraint = unsound_size_constraint,
                                       batch_assign = sound_bool,
                                       existing_clustering = sound_clustering))
  expect_error(hierarchical_clustering(distance_object = sound_distance_object,
                                       size_constraint = sound_size_constraint,
                                       batch_assign = unsound_bool,
                                       existing_clustering = sound_clustering))
  expect_error(hierarchical_clustering(distance_object = sound_distance_object,
                                       size_constraint = sound_size_constraint,
                                       batch_assign = sound_bool,
                                       existing_clustering = unsound_clustering))
})


# ==============================================================================
# hierarchical_clustering_internal
# ==============================================================================

test_that("`hierarchical_clustering_internal` checks input.", {
  expect_silent(hierarchical_clustering_internal(distance_object = sound_distance_object,
                                                 size_constraint = sound_size_constraint,
                                                 batch_assign = sound_bool,
                                                 existing_clustering = sound_clustering,
                                                 deep_copy = sound_bool))
  expect_error(hierarchical_clustering_internal(distance_object = unsound_distance_object,
                                                size_constraint = sound_size_constraint,
                                                batch_assign = sound_bool,
                                                existing_clustering = sound_clustering,
                                                deep_copy = sound_bool))
  expect_error(hierarchical_clustering_internal(distance_object = sound_distance_object,
                                                size_constraint = unsound_size_constraint,
                                                batch_assign = sound_bool,
                                                existing_clustering = sound_clustering,
                                                deep_copy = sound_bool))
  expect_error(hierarchical_clustering_internal(distance_object = sound_distance_object,
                                                size_constraint = sound_size_constraint,
                                                batch_assign = unsound_bool,
                                                existing_clustering = sound_clustering,
                                                deep_copy = sound_bool))
  expect_error(hierarchical_clustering_internal(distance_object = sound_distance_object,
                                                size_constraint = sound_size_constraint,
                                                batch_assign = sound_bool,
                                                existing_clustering = unsound_clustering,
                                                deep_copy = sound_bool))
  expect_error(hierarchical_clustering_internal(distance_object = sound_distance_object,
                                                size_constraint = sound_size_constraint,
                                                batch_assign = sound_bool,
                                                existing_clustering = sound_clustering,
                                                deep_copy = unsound_bool))
})


# ==============================================================================
# nng_clustering
# ==============================================================================

test_that("`nng_clustering` checks input.", {
  expect_silent(nng_clustering(distance_object = sound_distance_object,
                               size_constraint = sound_size_constraint,
                               seed_method = sound_seed_method,
                               main_unassigned_method = sound_unassigned_method,
                               main_radius = sound_radius,
                               main_data_points = sound_main_data_points,
                               secondary_unassigned_method = sound_unassigned_method,
                               secondary_radius = sound_radius))
  expect_error(nng_clustering(distance_object = unsound_distance_object,
                              size_constraint = sound_size_constraint,
                              seed_method = sound_seed_method,
                              main_unassigned_method = sound_unassigned_method,
                              main_radius = sound_radius,
                              main_data_points = sound_main_data_points,
                              secondary_unassigned_method = sound_unassigned_method,
                              secondary_radius = sound_radius))
  expect_error(nng_clustering(distance_object = sound_distance_object,
                              size_constraint = unsound_size_constraint,
                              seed_method = sound_seed_method,
                              main_unassigned_method = sound_unassigned_method,
                              main_radius = sound_radius,
                              main_data_points = sound_main_data_points,
                              secondary_unassigned_method = sound_unassigned_method,
                              secondary_radius = sound_radius))
  expect_error(nng_clustering(distance_object = sound_distance_object,
                              size_constraint = sound_size_constraint,
                              seed_method = unsound_seed_method,
                              main_unassigned_method = sound_unassigned_method,
                              main_radius = sound_radius,
                              main_data_points = sound_main_data_points,
                              secondary_unassigned_method = sound_unassigned_method,
                              secondary_radius = sound_radius))
  expect_error(nng_clustering(distance_object = sound_distance_object,
                              size_constraint = sound_size_constraint,
                              seed_method = sound_seed_method,
                              main_unassigned_method = unsound_unassigned_method,
                              main_radius = sound_radius,
                              main_data_points = sound_main_data_points,
                              secondary_unassigned_method = sound_unassigned_method,
                              secondary_radius = sound_radius))
  expect_error(nng_clustering(distance_object = sound_distance_object,
                              size_constraint = sound_size_constraint,
                              seed_method = sound_seed_method,
                              main_unassigned_method = sound_unassigned_method,
                              main_radius = unsound_radius,
                              main_data_points = sound_main_data_points,
                              secondary_unassigned_method = sound_unassigned_method,
                              secondary_radius = sound_radius))
  expect_error(nng_clustering(distance_object = sound_distance_object,
                              size_constraint = sound_size_constraint,
                              seed_method = sound_seed_method,
                              main_unassigned_method = sound_unassigned_method,
                              main_radius = sound_radius,
                              main_data_points = unsound_main_data_points,
                              secondary_unassigned_method = sound_unassigned_method,
                              secondary_radius = sound_radius))
  expect_error(nng_clustering(distance_object = sound_distance_object,
                              size_constraint = sound_size_constraint,
                              seed_method = sound_seed_method,
                              main_unassigned_method = sound_unassigned_method,
                              main_radius = sound_radius,
                              main_data_points = sound_main_data_points,
                              secondary_unassigned_method = unsound_unassigned_method,
                              secondary_radius = sound_radius))
  expect_error(nng_clustering(distance_object = sound_distance_object,
                              size_constraint = sound_size_constraint,
                              seed_method = sound_seed_method,
                              main_unassigned_method = sound_unassigned_method,
                              main_radius = sound_radius,
                              main_data_points = sound_main_data_points,
                              secondary_unassigned_method = sound_unassigned_method,
                              secondary_radius = unsound_radius))
})


# ==============================================================================
# nng_clustering_batches
# ==============================================================================

sound_batch_size <- 100L
unsound_batch_size <- -100L

test_that("`nng_clustering_batches` checks input.", {
  expect_silent(nng_clustering_batches(distance_object = sound_distance_object,
                                       size_constraint = sound_size_constraint,
                                       main_unassigned_method = sound_unassigned_method,
                                       main_radius = sound_radius,
                                       main_data_points = sound_main_data_points,
                                       batch_size = sound_batch_size))
  expect_error(nng_clustering_batches(distance_object = unsound_distance_object,
                                      size_constraint = sound_size_constraint,
                                      main_unassigned_method = sound_unassigned_method,
                                      main_radius = sound_radius,
                                      main_data_points = sound_main_data_points,
                                      batch_size = sound_batch_size))
  expect_error(nng_clustering_batches(distance_object = sound_distance_object,
                                      size_constraint = unsound_size_constraint,
                                      main_unassigned_method = sound_unassigned_method,
                                      main_radius = sound_radius,
                                      main_data_points = sound_main_data_points,
                                      batch_size = sound_batch_size))
  expect_error(nng_clustering_batches(distance_object = sound_distance_object,
                                      size_constraint = sound_size_constraint,
                                      main_unassigned_method = unsound_unassigned_method,
                                      main_radius = sound_radius,
                                      main_data_points = sound_main_data_points,
                                      batch_size = sound_batch_size))
  expect_error(nng_clustering_batches(distance_object = sound_distance_object,
                                      size_constraint = sound_size_constraint,
                                      main_unassigned_method = sound_unassigned_method,
                                      main_radius = unsound_radius,
                                      main_data_points = sound_main_data_points,
                                      batch_size = sound_batch_size))
  expect_error(nng_clustering_batches(distance_object = sound_distance_object,
                                      size_constraint = sound_size_constraint,
                                      main_unassigned_method = sound_unassigned_method,
                                      main_radius = sound_radius,
                                      main_data_points = unsound_main_data_points,
                                      batch_size = sound_batch_size))
  expect_error(nng_clustering_batches(distance_object = sound_distance_object,
                                      size_constraint = sound_size_constraint,
                                      main_unassigned_method = sound_unassigned_method,
                                      main_radius = sound_radius,
                                      main_data_points = sound_main_data_points,
                                      batch_size = unsound_batch_size))
})


# ==============================================================================
# nng_clustering_types
# ==============================================================================

test_that("`nng_clustering_types` checks input.", {
  expect_silent(nng_clustering_types(distance_object = sound_distance_object,
                                     type_labels = sound_type_labels,
                                     type_size_constraints = sound_type_size_constraints,
                                     total_size_constraint = sound_total_size_constraint,
                                     seed_method = sound_seed_method,
                                     main_unassigned_method = sound_unassigned_method,
                                     main_radius = sound_radius,
                                     main_data_points = sound_main_data_points,
                                     secondary_unassigned_method = sound_unassigned_method,
                                     secondary_radius = sound_radius))
  expect_error(nng_clustering_types(distance_object = unsound_distance_object,
                                    type_labels = sound_type_labels,
                                    type_size_constraints = sound_type_size_constraints,
                                    total_size_constraint = sound_total_size_constraint,
                                    seed_method = sound_seed_method,
                                    main_unassigned_method = sound_unassigned_method,
                                    main_radius = sound_radius,
                                    main_data_points = sound_main_data_points,
                                    secondary_unassigned_method = sound_unassigned_method,
                                    secondary_radius = sound_radius))
  expect_error(nng_clustering_types(distance_object = sound_distance_object,
                                    type_labels = unsound_type_labels,
                                    type_size_constraints = sound_type_size_constraints,
                                    total_size_constraint = sound_total_size_constraint,
                                    seed_method = sound_seed_method,
                                    main_unassigned_method = sound_unassigned_method,
                                    main_radius = sound_radius,
                                    main_data_points = sound_main_data_points,
                                    secondary_unassigned_method = sound_unassigned_method,
                                    secondary_radius = sound_radius))
  expect_error(nng_clustering_types(distance_object = sound_distance_object,
                                    type_labels = sound_type_labels,
                                    type_size_constraints = unsound_type_size_constraints,
                                    total_size_constraint = sound_total_size_constraint,
                                    seed_method = sound_seed_method,
                                    main_unassigned_method = sound_unassigned_method,
                                    main_radius = sound_radius,
                                    main_data_points = sound_main_data_points,
                                    secondary_unassigned_method = sound_unassigned_method,
                                    secondary_radius = sound_radius))
  expect_error(nng_clustering_types(distance_object = sound_distance_object,
                                    type_labels = sound_type_labels,
                                    type_size_constraints = sound_type_size_constraints,
                                    total_size_constraint = unsound_total_size_constraint,
                                    seed_method = sound_seed_method,
                                    main_unassigned_method = sound_unassigned_method,
                                    main_radius = sound_radius,
                                    main_data_points = sound_main_data_points,
                                    secondary_unassigned_method = sound_unassigned_method,
                                    secondary_radius = sound_radius))
  expect_error(nng_clustering_types(distance_object = sound_distance_object,
                                    type_labels = sound_type_labels,
                                    type_size_constraints = sound_type_size_constraints,
                                    total_size_constraint = sound_total_size_constraint,
                                    seed_method = unsound_seed_method,
                                    main_unassigned_method = sound_unassigned_method,
                                    main_radius = sound_radius,
                                    main_data_points = sound_main_data_points,
                                    secondary_unassigned_method = sound_unassigned_method,
                                    secondary_radius = sound_radius))
  expect_error(nng_clustering_types(distance_object = sound_distance_object,
                                    type_labels = sound_type_labels,
                                    type_size_constraints = sound_type_size_constraints,
                                    total_size_constraint = sound_total_size_constraint,
                                    seed_method = sound_seed_method,
                                    main_unassigned_method = unsound_unassigned_method,
                                    main_radius = sound_radius,
                                    main_data_points = sound_main_data_points,
                                    secondary_unassigned_method = sound_unassigned_method,
                                    secondary_radius = sound_radius))
  expect_error(nng_clustering_types(distance_object = sound_distance_object,
                                    type_labels = sound_type_labels,
                                    type_size_constraints = sound_type_size_constraints,
                                    total_size_constraint = sound_total_size_constraint,
                                    seed_method = sound_seed_method,
                                    main_unassigned_method = sound_unassigned_method,
                                    main_radius = unsound_radius,
                                    main_data_points = sound_main_data_points,
                                    secondary_unassigned_method = sound_unassigned_method,
                                    secondary_radius = sound_radius))
  expect_error(nng_clustering_types(distance_object = sound_distance_object,
                                    type_labels = sound_type_labels,
                                    type_size_constraints = sound_type_size_constraints,
                                    total_size_constraint = sound_total_size_constraint,
                                    seed_method = sound_seed_method,
                                    main_unassigned_method = sound_unassigned_method,
                                    main_radius = sound_radius,
                                    main_data_points = unsound_main_data_points,
                                    secondary_unassigned_method = sound_unassigned_method,
                                    secondary_radius = sound_radius))
  expect_error(nng_clustering_types(distance_object = sound_distance_object,
                                    type_labels = sound_type_labels,
                                    type_size_constraints = sound_type_size_constraints,
                                    total_size_constraint = sound_total_size_constraint,
                                    seed_method = sound_seed_method,
                                    main_unassigned_method = sound_unassigned_method,
                                    main_radius = sound_radius,
                                    main_data_points = sound_main_data_points,
                                    secondary_unassigned_method = unsound_unassigned_method,
                                    secondary_radius = sound_radius))
  expect_error(nng_clustering_types(distance_object = sound_distance_object,
                                    type_labels = sound_type_labels,
                                    type_size_constraints = sound_type_size_constraints,
                                    total_size_constraint = sound_total_size_constraint,
                                    seed_method = sound_seed_method,
                                    main_unassigned_method = sound_unassigned_method,
                                    main_radius = sound_radius,
                                    main_data_points = sound_main_data_points,
                                    secondary_unassigned_method = sound_unassigned_method,
                                    secondary_radius = unsound_radius))
})


# ==============================================================================
# check_clustering
# ==============================================================================

test_that("`check_clustering` checks input.", {
  expect_silent(check_clustering(clustering = sound_clustering,
                                 size_constraint = sound_size_constraint))
  expect_error(check_clustering(clustering = unsound_clustering,
                                size_constraint = sound_size_constraint))
  expect_error(check_clustering(clustering = sound_clustering,
                                size_constraint = unsound_size_constraint))
})


# ==============================================================================
# check_clustering_types
# ==============================================================================

test_that("`check_clustering_types` checks input.", {
  expect_silent(check_clustering_types(clustering = sound_clustering,
                                       type_labels = sound_type_labels,
                                       type_size_constraints = sound_type_size_constraints,
                                       total_size_constraint = sound_total_size_constraint))
  expect_error(check_clustering_types(clustering = unsound_clustering,
                                      type_labels = sound_type_labels,
                                      type_size_constraints = sound_type_size_constraints,
                                      total_size_constraint = sound_total_size_constraint))
  expect_error(check_clustering_types(clustering = sound_clustering,
                                      type_labels = unsound_type_labels,
                                      type_size_constraints = sound_type_size_constraints,
                                      total_size_constraint = sound_total_size_constraint))
  expect_error(check_clustering_types(clustering = sound_clustering,
                                      type_labels = sound_type_labels,
                                      type_size_constraints = unsound_type_size_constraints,
                                      total_size_constraint = sound_total_size_constraint))
  expect_error(check_clustering_types(clustering = sound_clustering,
                                      type_labels = sound_type_labels,
                                      type_size_constraints = sound_type_size_constraints,
                                      total_size_constraint = unsound_total_size_constraint))
})


# ==============================================================================
# get_clustering_stats
# ==============================================================================

test_that("`get_clustering_stats` checks input.", {
  expect_silent(get_clustering_stats(clustering = sound_clustering,
                                     distance_object = sound_distance_object))
  expect_error(get_clustering_stats(clustering = unsound_clustering,
                                    distance_object = sound_distance_object))
  expect_error(get_clustering_stats(clustering = sound_clustering,
                                    distance_object = unsound_distance_object))
})
