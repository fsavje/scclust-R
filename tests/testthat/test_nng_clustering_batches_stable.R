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
context("nng_clustering_batches (stable)")

source("config.R", local = TRUE)
source("../replica/replica_make_nng.R", local = TRUE)
source("../replica/replica_nng_batches.R", local = TRUE)
source("utils_nng.R", local = TRUE)


test_that("`nng_clustering_batches` (stable) returns correct output", {
  skip_on_cran()
  skip_if_not(run_slow_tests)
  skip_if_not(compiled_with_stable_nng, "Only run this when scclust is compiled with the -DSCC_STABLE_NNG flag.")

  for (size_constraint in c(2L, 3L, 6L)) {
    for (unassigned_method in c("ignore", "any_neighbor")) {
      radius_to_use <- c(0, 0.1, 0.2, 0.3)
      if (size_constraint == 6L) {
        radius_to_use <- c(0, 0.3, 0.4)
      }
      for (radius in radius_to_use) {
        for (do_primary_data_points in c("N", "Y")) {
          for (batch_size in c(1L, 10L, 100L, 1000L)) {
            use_primary_data_points <- NULL
            if (do_primary_data_points == "Y") use_primary_data_points <- primary_data_points
            use_radius <- radius
            if (use_radius == 0) use_radius <- NULL

            #cat(c(size_constraint,
            #      unassigned_method,
            #      use_radius,
            #      is.null(use_primary_data_points),
            #      batch_size), "\n")
            test_nng_batch_against_replica(test_distances1,
                                           size_constraint,
                                           unassigned_method,
                                           use_radius,
                                           use_primary_data_points,
                                           batch_size)
          }
        }
      }
    }
  }
})
