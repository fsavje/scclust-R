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
context("nng_clustering_types (stable)")

source("config.R", local = TRUE)
source("../replica/replica_findseed.R", local = TRUE)
source("../replica/replica_make_nng.R", local = TRUE)
source("../replica/replica_nng.R", local = TRUE)
source("utils_nng.R", local = TRUE)


by_nng_or_closest_assigned <- function(unassigned_method,
                                       unassigned,
                                       nng,
                                       assigned,
                                       radius,
                                       cl_label,
                                       distances) {
  if (unassigned_method == "by_nng") {
    for(i in which(unassigned)) {
      pick <- intersect(which(nng[, i]), which(assigned))
      if (length(pick) > 0) {
        cl_label[i] <- cl_label[pick[1]]
        unassigned[i] <- FALSE
      }
    }
  } else if (any(unassigned)) {
    cl_label <- match_n_assign(cl_label, which(assigned), unassigned, radius, distances)
  }

  cl_label
}


test_that("`nng_clustering_types` (stable) returns correct output", {
  skip_on_cran()
  skip_if_not(run_slow_tests)
  skip_if_not(compiled_with_stable_nng, "Only run this when scclust is compiled with the -DSCC_STABLE_NNG flag.")
  skip_if_not(!compiled_with_stable_findseed, "Only run this when scclust is *not* compiled with the -DSCC_STABLE_FINDSEED flag.")

  for (tc in list(list(tc = "tc1",
                       type_labels = types1,
                       type_size_constraints = c("0" = 1L, "1" = 1L, "2" = 1L, "3" = 1L),
                       total_size_constraint = NULL),
                  list(tc = "tc2",
                       type_labels = types1,
                       type_size_constraints = c("0" = 1L, "1" = 2L, "2" = 2L, "3" = 1L),
                       total_size_constraint = 8L),
                  list(tc = "tc3",
                       type_labels = types1,
                       type_size_constraints = c("0" = 3L, "1" = 3L, "2" = 3L, "3" = 3L),
                       total_size_constraint = 12L),
                  list(tc = "tc4",
                       type_labels = types1,
                       type_size_constraints = c("0" = 1L, "1" = 0L, "2" = 1L, "3" = 0L),
                       total_size_constraint = 4L),
                  list(tc = "tc5",
                       type_labels = types2,
                       type_size_constraints = c("0" = 1L, "1" = 1L),
                       total_size_constraint = 2L),
                  list(tc = "tc6",
                       type_labels = types2,
                       type_size_constraints = c("0" = 2L, "1" = 0L),
                       total_size_constraint = 4L),
                  list(tc = "tc7",
                       type_labels = types2,
                       type_size_constraints = c("0" = 2L, "1" = 2L),
                       total_size_constraint = NULL))) {

    for (seed_method in c("lexical",
                          "inwards_order",
                          "inwards_updating",
                          "inwards_alt_updating",
                          "exclusion_order",
                          "exclusion_updating")) {

      for (unassigned_method in c("ignore",
                                  "by_nng",
                                  "closest_assigned",
                                  "closest_seed",
                                  "estimated_radius_closest_seed")) {

        radius_to_use <- c(0, 0.3, 0.4)
        if (isTRUE(tc$total_size_constraint > 6)) {
          radius_to_use <- c(0, 0.8)
        }
        if (isTRUE(tc$total_size_constraint > 10)) {
          radius_to_use <- c(0, 1.2)
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

              #cat(c(tc$tc,
              #      seed_method,
              #      unassigned_method,
              #      use_radius,
              #      is.null(use_primary_data_points),
              #      use_secondary_unassigned_method,
              #      use_secondary_radius), "\n")
              test_nng_types_against_replica(test_distances1,
                                             tc$type_labels,
                                             tc$type_size_constraints,
                                             tc$total_size_constraint,
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
