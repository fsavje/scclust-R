library(Rscclust)
context("nng_clustering.R (no types)")

source("../replica/replica_findseed.R", local = TRUE)
source("../replica/replica_nng.R", local = TRUE)
source("utils_nng.R", local = TRUE)

test_that("non-type nng clustering function cluster correctly.", {

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
                           "inwards_alt_updating",
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
                           main_data_points,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           2L,
                           "inwards_updating",
                           "ignore",
                           NULL,
                           main_data_points,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           2L,
                           "inwards_alt_updating",
                           "ignore",
                           NULL,
                           main_data_points,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "exclusion_order",
                           "ignore",
                           NULL,
                           main_data_points,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           2L,
                           "exclusion_updating",
                           "ignore",
                           NULL,
                           main_data_points,
                           "ignore",
                           NULL)

  test_nng_against_replica(test_distances1,
                           2L,
                           "lexical",
                           "by_nng",
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
                           "by_nng",
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
                           main_data_points,
                           "ignore",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           main_data_points,
                           "closest_assigned",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           main_data_points,
                           "closest_seed",
                           NULL)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           main_data_points,
                           "estimated_radius_closest_seed",
                           NULL)

  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           main_data_points,
                           "ignore",
                           test_radius)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           main_data_points,
                           "closest_assigned",
                           test_radius)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           main_data_points,
                           "closest_seed",
                           test_radius)
  test_nng_against_replica(test_distances1,
                           3L,
                           "lexical",
                           "ignore",
                           NULL,
                           main_data_points,
                           "estimated_radius_closest_seed",
                           test_radius)
})


test_that("non-type nng clustering function cluster correctly all combinations", {
  skip_on_cran()
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

      sc_to_test <- c(2L, 3L, 6L)
      # scclust and the replica sorts the arcs for a vertex differently,
      # this causes problems with all updating seed finding functions
      # when there's more than 1 arc (i.e., `size_constraint` >= 2).
      # Here we only test options that produces the same results. Compiling
      # scclust with the -DSCC_STABLE_NNG flag make the arc order the same as
      # in R. See <test_nng_clustering_notypes_stable.R>.
      if (seed_method %in% c("inwards_updating",
                             "inwards_alt_updating",
                             "exclusion_updating")) {
        sc_to_test <- 2L
      }

      for (size_constraint in sc_to_test) {

        radius_to_use <- c(0, 0.1, 0.2)
        if (size_constraint == 6L) {
          radius_to_use <- c(0, 0.3)
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
