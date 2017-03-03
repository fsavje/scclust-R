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
context("Input checking in exported functions")


# ==============================================================================
# Shared objects
# ==============================================================================

sound_distance_object <- distances::distances(matrix(c(1, 4, 3, 2, 45, 6, 3, 2, 6, 5,
                                                 34, 2, 4, 6, 4, 6, 4, 2, 7, 8), nrow = 10))
unsound_distance_object <- letters[1:10]
sound_size_constraint <- 2L
unsound_size_constraint <- 0L
sound_clustering <- scc_clustering(c(rep("a", 5), rep("b", 5)))
unsound_clustering <- c(rep("a", 5), rep("b", 5))
sound_bool <- TRUE
unsound_bool <- "a"
sound_seed_method <- "lexical"
unsound_seed_method <- "unknown"
sound_unassigned_method <- "ignore"
unsound_unassigned_method <- "unknown"
sound_radius <- NULL
unsound_radius <- TRUE
sound_primary_data_points <- rep(TRUE, 10)
unsound_primary_data_points <- rep(TRUE, 5)
sound_type_labels <- rep(c(1L, 2L), 5)
unsound_type_labels <- rep(c(1L, 2L), 4)
sound_type_size_constraints <- c("1" = 1L, "2" = 1L)
unsound_type_size_constraints <- c(1L, 1L)
sound_total_size_constraint <- 3L
unsound_total_size_constraint <- 100L


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
# make_clustering
# ==============================================================================

sound_batch_size <- 100L
unsound_batch_size <- -100L

test_make_clustering <- function(distance_object = sound_distance_object,
                                 size_constraint = sound_total_size_constraint,
                                 type_labels = sound_type_labels,
                                 type_constraints = sound_type_size_constraints,
                                 seed_method = sound_seed_method,
                                 primary_data_points = sound_primary_data_points,
                                 primary_unassigned_method = sound_unassigned_method,
                                 secondary_unassigned_method = sound_unassigned_method,
                                 seed_radius = sound_radius,
                                 primary_radius = sound_radius,
                                 secondary_radius = sound_radius,
                                 batch_size = sound_batch_size) {
  make_clustering(distance_object = distance_object,
                  size_constraint = size_constraint,
                  type_labels = type_labels,
                  type_constraints = type_constraints,
                  seed_method = seed_method,
                  primary_data_points = primary_data_points,
                  primary_unassigned_method = primary_unassigned_method,
                  secondary_unassigned_method = secondary_unassigned_method,
                  seed_radius = seed_radius,
                  primary_radius = primary_radius,
                  secondary_radius = secondary_radius,
                  batch_size = batch_size)
}


test_that("`nng_clustering` checks input.", {
  expect_silent(test_make_clustering())
  expect_error(test_make_clustering(distance_object = unsound_distance_object))
  expect_error(test_make_clustering(size_constraint = unsound_total_size_constraint))
  expect_error(test_make_clustering(type_labels = unsound_type_labels))
  expect_error(test_make_clustering(type_constraints = unsound_type_size_constraints))
  expect_error(test_make_clustering(seed_method = unsound_seed_method))
  expect_error(test_make_clustering(primary_data_points = unsound_primary_data_points))
  expect_error(test_make_clustering(primary_unassigned_method = unsound_unassigned_method))
  expect_error(test_make_clustering(secondary_unassigned_method = unsound_unassigned_method))
  expect_error(test_make_clustering(seed_radius = unsound_radius))
  expect_error(test_make_clustering(primary_radius = unsound_radius))
  expect_error(test_make_clustering(secondary_radius = unsound_radius))
  expect_error(test_make_clustering(batch_size = unsound_batch_size))
})


# ==============================================================================
# scc_clustering methods
# ==============================================================================

sound_cluster_labels <- 1:10
unsound_cluster_labels <- dist(1:10)
sound_unassigned_labels <- c(1L, 3L)
unsound_unassigned_labels <- c(1L, "a")
sound_ids <- letters[1:10]
unsound_ids <- letters[1:5]

test_that("`scc_clustering` checks input.", {
  expect_silent(scc_clustering(cluster_labels = sound_cluster_labels,
                               unassigned_labels = sound_unassigned_labels,
                               ids = sound_ids))
  expect_error(scc_clustering(cluster_labels = unsound_cluster_labels,
                              unassigned_labels = sound_unassigned_labels,
                              ids = sound_ids))
  expect_error(scc_clustering(cluster_labels = sound_cluster_labels,
                              unassigned_labels = unsound_unassigned_labels,
                              ids = sound_ids))
  expect_error(scc_clustering(cluster_labels = sound_cluster_labels,
                              unassigned_labels = sound_unassigned_labels,
                              ids = unsound_ids))
})

test_that("`cluster_count` checks input.", {
  expect_silent(cluster_count(sound_clustering))
  expect_error(cluster_count(unsound_clustering))
})

test_that("`as.data.frame.scc_clustering` checks input.", {
  expect_silent(as.data.frame.scc_clustering(sound_clustering))
  expect_error(as.data.frame.scc_clustering(unsound_clustering))
})

test_that("`print.scc_clustering` checks input.", {
  expect_output(print.scc_clustering(sound_clustering))
  expect_error(print.scc_clustering(unsound_clustering))
})


# ==============================================================================
# check_clustering
# ==============================================================================

test_that("`check_clustering` checks input.", {
  expect_silent(check_clustering(clustering = sound_clustering,
                                 size_constraint = sound_total_size_constraint,
                                 type_labels = sound_type_labels,
                                 type_constraints = sound_type_size_constraints))
  expect_error(check_clustering(clustering = unsound_clustering,
                                size_constraint = sound_total_size_constraint,
                                type_labels = sound_type_labels,
                                type_constraints = sound_type_size_constraints))
  expect_error(check_clustering(clustering = sound_clustering,
                                size_constraint = sound_total_size_constraint,
                                type_labels = unsound_type_labels,
                                type_constraints = sound_type_size_constraints))
  expect_error(check_clustering(clustering = sound_clustering,
                                size_constraint = sound_total_size_constraint,
                                type_labels = sound_type_labels,
                                type_constraints = unsound_type_size_constraints))
  expect_error(check_clustering(clustering = sound_clustering,
                                size_constraint = unsound_total_size_constraint,
                                type_labels = sound_type_labels,
                                type_constraints = sound_type_size_constraints))
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
