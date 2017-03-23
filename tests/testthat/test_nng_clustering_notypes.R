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
context("nng_clustering")

source("config.R", local = TRUE)
source("../replica/replica_findseed.R", local = TRUE)
source("../replica/replica_make_nng.R", local = TRUE)
source("../replica/replica_nng.R", local = TRUE)
source("utils_nng.R", local = TRUE)


test_that("`nng_clustering` returns correct output", {
  skip_if_not(!compiled_with_stable_nng, "Only run this when scclust is *not* compiled with the -DSCC_STABLE_NNG flag.")
  skip_if_not(!compiled_with_stable_findseed, "Only run this when scclust is *not* compiled with the -DSCC_STABLE_FINDSEED flag.")

  test_nng_against_replica(test_distances1,
                           2L,
                           "lexical",
                           "ignore",
                           NULL,
                           NULL,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           NULL,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           10L,
                           "lexical",
                           "ignore",
                           NULL,
                           NULL,
                           "ignore",
                           NULL)

  test_nng_against_replica(test_distances1,
                           3L,
                           "inwards_order",
                           "ignore",
                           NULL,
                           NULL,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           2L,
                           "inwards_updating",
                           "ignore",
                           NULL,
                           NULL,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           2L,
                           "inwards_updating",
                           "ignore",
                           NULL,
                           NULL,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "exclusion_order",
                           "ignore",
                           NULL,
                           NULL,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           2L,
                           "exclusion_updating",
                           "ignore",
                           NULL,
                           NULL,
                           "ignore",
                           NULL)

  test_nng_against_replica(test_distances1,
                           3L,
                           "inwards_order",
                           "ignore",
                           NULL,
                           which(primary_data_points),
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           2L,
                           "inwards_updating",
                           "ignore",
                           NULL,
                           which(primary_data_points),
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "inwards_order",
                           "ignore",
                           NULL,
                           primary_data_points,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           2L,
                           "inwards_updating",
                           "ignore",
                           NULL,
                           primary_data_points,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           2L,
                           "inwards_updating",
                           "ignore",
                           NULL,
                           primary_data_points,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "exclusion_order",
                           "ignore",
                           NULL,
                           primary_data_points,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           2L,
                           "exclusion_updating",
                           "ignore",
                           NULL,
                           primary_data_points,
                           "ignore",
                           NULL)

  test_nng_against_replica(test_distances1,
                           2L,
                           "lexical",
                           "any_neighbor",
                           NULL,
                           NULL,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "closest_assigned",
                           NULL,
                           NULL,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "closest_seed",
                           NULL,
                           NULL,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "estimated_radius_closest_seed",
                           NULL,
                           NULL,
                           "ignore",
                           NULL)

  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           test_radius,
                           NULL,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           2L,
                           "lexical",
                           "any_neighbor",
                           test_radius,
                           NULL,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "closest_assigned",
                           test_radius,
                           NULL,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "closest_seed",
                           test_radius,
                           NULL,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "estimated_radius_closest_seed",
                           test_radius,
                           NULL,
                           "ignore",
                           NULL)

  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           primary_data_points,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           primary_data_points,
                           "closest_assigned",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           primary_data_points,
                           "closest_seed",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           primary_data_points,
                           "estimated_radius_closest_seed",
                           NULL)

  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           primary_data_points,
                           "ignore",
                           test_radius)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           primary_data_points,
                           "closest_assigned",
                           test_radius)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           primary_data_points,
                           "closest_seed",
                           test_radius)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           primary_data_points,
                           "estimated_radius_closest_seed",
                           test_radius)
})


test_that("`nng_clustering` returns correct output (combinations)", {
  skip_on_cran()
  skip_if_not(run_slow_tests)
  skip_if_not(!compiled_with_stable_nng, "Only run this when scclust is *not* compiled with the -DSCC_STABLE_NNG flag.")
  skip_if_not(!compiled_with_stable_findseed, "Only run this when scclust is *not* compiled with the -DSCC_STABLE_FINDSEED flag.")

  for (seed_method in c("lexical",
                        "inwards_order",
                        "inwards_updating",
                        "exclusion_order",
                        "exclusion_updating")) {
    for (unassigned_method in c("ignore",
                                "any_neighbor",
                                "closest_assigned",
                                "closest_seed",
                                "estimated_radius_closest_seed")) {

      sc_to_test <- c(2L, 3L, 6L)
      # scclust and the replica sorts the arcs for a vertex differently,
      # this causes problems with all updating seed finding functions
      # when there's more than 1 arc (i.e., `size_constraint` >= 2).
      # Here we only test options that produces the same results. Compiling
      # scclust with the -DSCC_STABLE_NNG flag make the arc order the same as
      # in R. See <test_nng_clustering_notypes_stable.R>.
      if (seed_method %in% c("inwards_updating",
                             "exclusion_updating")) {
        sc_to_test <- 2L
      }

      for (size_constraint in sc_to_test) {

        radius_to_use <- c(0, 0.1, 0.2)
        if (size_constraint == 6L) {
          radius_to_use <- c(0, 0.3)
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
