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
context("nng_clustering (stable)")

source("config.R", local = TRUE)
source("../replica/replica_findseed.R", local = TRUE)
source("../replica/replica_make_nng.R", local = TRUE)
source("../replica/replica_nng.R", local = TRUE)
source("../replica/replica_nng_stable.R", local = TRUE)
source("utils_nng.R", local = TRUE)


test_that("`nng_clustering` (stable) returns correct output", {
  skip_on_cran()
  skip_if_not(run_slow_tests)
  skip_if_not(compiled_with_stable_nng, "Only run this when scclust is compiled with the -DSCC_STABLE_NNG flag.")
  skip_if_not(!compiled_with_stable_findseed, "Only run this when scclust is *not* compiled with the -DSCC_STABLE_FINDSEED flag.")

  for (seed_method in c("lexical",
                        "inwards_order",
                        "inwards_updating",
                        "inwards_alt_updating",
                        "exclusion_order",
                        "exclusion_updating")) {
    for (unassigned_method in c("ignore",
                                "any_neighbor",
                                "closest_assigned",
                                "closest_seed",
                                "estimated_radius_closest_seed")) {
      for (size_constraint in c(2L, 3L, 6L)) {

        radius_to_use <- c(0, 0.1, 0.2, 0.3)
        if (size_constraint == 6L) {
          radius_to_use <- c(0, 0.3, 0.4)
        }

        for (radius in radius_to_use) {
          for (secondary_unassigned_method in c("ignore",
                                                "ignore_withmain",
                                                "closest_assigned",
                                                "closest_seed",
                                                "estimated_radius_closest_seed")) {
            for (secondary_radius in radius_to_use) {

              use_primary_data_points <- NULL
              if (secondary_unassigned_method != "ignore") use_primary_data_points <- primary_data_points

              use_secondary_unassigned_method <- secondary_unassigned_method
              if (use_secondary_unassigned_method == "ignore_withmain") use_secondary_unassigned_method <- "ignore"

              use_radius <- radius
              if (use_radius == 0) use_radius <- NULL
              use_secondary_radius <- secondary_radius
              if (use_secondary_radius == 0) use_secondary_radius <- NULL

              #cat(c(size_constraint,
              #      seed_method,
              #      unassigned_method,
              #      use_radius,
              #      is.null(use_primary_data_points),
              #      use_secondary_unassigned_method,
              #      use_secondary_radius), "\n")
              test_nng_against_replica(test_distances1,
                                       size_constraint,
                                       seed_method,
                                       unassigned_method,
                                       use_radius,
                                       use_primary_data_points,
                                       use_secondary_unassigned_method,
                                       use_secondary_radius)
            }
          }
        }
      }
    }
  }
})
