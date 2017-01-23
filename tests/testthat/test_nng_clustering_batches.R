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
context("nng_clustering_batches")

source("config.R", local = TRUE)
source("../replica/replica_make_nng.R", local = TRUE)
source("../replica/replica_nng_batches.R", local = TRUE)
source("utils_nng.R", local = TRUE)


test_that("`nng_clustering_batches` returns correct output", {

  test_nng_batch_against_replica(test_distances1,
                                 2L,
                                 "ignore",
                                 NULL,
                                 NULL,
                                 10L)
  test_nng_batch_against_replica(test_distances1,
                                 3L,
                                 "ignore",
                                 NULL,
                                 NULL,
                                 100L)
  test_nng_batch_against_replica(test_distances1,
                                 6L,
                                 "ignore",
                                 NULL,
                                 NULL,
                                 10L)
  test_nng_batch_against_replica(test_distances1,
                                 10L,
                                 "ignore",
                                 NULL,
                                 NULL,
                                 100L)
  test_nng_batch_against_replica(test_distances1,
                                 2L,
                                 "any_neighbor",
                                 NULL,
                                 NULL,
                                 10L)
  test_nng_batch_against_replica(test_distances1,
                                 2L,
                                 "ignore",
                                 NULL,
                                 primary_data_points,
                                 10L)
  test_nng_batch_against_replica(test_distances1,
                                 3L,
                                 "ignore",
                                 NULL,
                                 primary_data_points,
                                 100L)
  test_nng_batch_against_replica(test_distances1,
                                 6L,
                                 "ignore",
                                 NULL,
                                 primary_data_points,
                                 10L)
  test_nng_batch_against_replica(test_distances1,
                                 10L,
                                 "ignore",
                                 NULL,
                                 primary_data_points,
                                 100L)
  test_nng_batch_against_replica(test_distances1,
                                 2L,
                                 "any_neighbor",
                                 NULL,
                                 primary_data_points,
                                 10L)
  test_nng_batch_against_replica(test_distances1,
                                 2L,
                                 "ignore",
                                 test_radius,
                                 NULL,
                                 10L)
  test_nng_batch_against_replica(test_distances1,
                                 3L,
                                 "ignore",
                                 test_radius,
                                 NULL,
                                 100L)
  test_nng_batch_against_replica(test_distances1,
                                 2L,
                                 "any_neighbor",
                                 test_radius,
                                 NULL,
                                 10L)
  test_nng_batch_against_replica(test_distances1,
                                 2L,
                                 "ignore",
                                 test_radius,
                                 primary_data_points,
                                 10L)
  test_nng_batch_against_replica(test_distances1,
                                 3L,
                                 "ignore",
                                 test_radius,
                                 primary_data_points,
                                 100L)
  test_nng_batch_against_replica(test_distances1,
                                 2L,
                                 "any_neighbor",
                                 test_radius,
                                 primary_data_points,
                                 10L)
})


test_that("`nng_clustering_batches` returns correct output (combinations)", {
  skip_on_cran()
  skip_if_not(run_slow_tests)

  for (unassigned_method in c("ignore", "any_neighbor")) {
    for (radius in c(0, 0.1, 0.2)) {
      for (do_primary_data_points in c("N", "Y")) {
        for (batch_size in c(1L, 10L, 100L, 1000L)) {
          use_primary_data_points <- NULL
          if (do_primary_data_points == "Y") use_primary_data_points <- primary_data_points
          use_radius <- radius
          if (use_radius == 0) use_radius <- NULL

          #cat(c(unassigned_method,
          #      use_radius,
          #      is.null(use_primary_data_points),
          #      batch_size), "\n")
          test_nng_batch_against_replica(test_distances1,
                                         2L,
                                         unassigned_method,
                                         use_radius,
                                         use_primary_data_points,
                                         batch_size)
        }
      }
    }
  }
})
