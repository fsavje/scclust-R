library(Rscclust)
context("nng_clustering.R (no types, stable)")

source("../replica/replica_findseed.R", local = TRUE)
source("../replica/replica_nng.R", local = TRUE)
source("../replica/replica_nng_stable.R", local = TRUE)
source("utils_nng.R", local = TRUE)

test_that("non-type nng clustering function cluster correctly all combinations", {
  skip("Only run this when scclust is compiled with the -DSCC_STABLE_NNG flag.")
  for (seed_method in c("lexical",
                        "inwards_order",
                        "inwards_updating",
                        "inwards_alt_updating",
                        "exclusion_order",
                        "exclusion_updating")) {
    for (main_unassigned_method in c("ignore",
                                     "by_nng",
                                     "closest_assigned",
                                     "closest_seed",
                                     "estimated_radius_closest_seed")) {
      for (size_constraint in c(2L, 3L, 6L)) {

        radius_to_use <- c(0, 0.1, 0.2, 0.3)
        if (size_constraint == 6L) {
          radius_to_use <- c(0, 0.3, 0.4)
        }

        for (main_radius in radius_to_use) {
          for (secondary_unassigned_method in c("ignore",
                                                "ignore_withmain",
                                                "closest_assigned",
                                                "closest_seed",
                                                "estimated_radius_closest_seed")) {
            for (secondary_radius in radius_to_use) {

              use_main_data_points <- NULL
              if (secondary_unassigned_method != "ignore") use_main_data_points <- main_data_points

              use_secondary_unassigned_method <- secondary_unassigned_method
              if (use_secondary_unassigned_method == "ignore_withmain") use_secondary_unassigned_method <- "ignore"

              use_main_radius <- main_radius
              if (use_main_radius == 0) use_main_radius <- NULL
              use_secondary_radius <- secondary_radius
              if (use_secondary_radius == 0) use_secondary_radius <- NULL

              cat(c(size_constraint,
                    seed_method,
                    main_unassigned_method,
                    use_main_radius,
                    is.null(use_main_data_points),
                    use_secondary_unassigned_method,
                    use_secondary_radius), "\n")
              test_nng_against_replica(test_distances1,
                                       size_constraint,
                                       seed_method,
                                       main_unassigned_method,
                                       use_main_radius,
                                       use_main_data_points,
                                       use_secondary_unassigned_method,
                                       use_secondary_radius)
            }
          }
        }
      }
    }
  }
})
