library(Rscclust)
context("nng_clustering.R (batches, stable)")

source("config.R", local = TRUE)
source("../replica/replica_make_nng.R", local = TRUE)
source("../replica/replica_nng_batches.R", local = TRUE)
source("utils_nng.R", local = TRUE)


test_that("batch nng clustering function cluster correctly all combinations", {
  skip_on_cran()
  skip_if_not(run_slow_tests)
  skip_if_not(compiled_with_stable_nng, "Only run this when scclust is compiled with the -DSCC_STABLE_NNG flag.")

  for (size_constraint in c(2L, 3L, 6L)) {
    for (main_unassigned_method in c("ignore", "by_nng")) {
      radius_to_use <- c(0, 0.1, 0.2, 0.3)
      if (size_constraint == 6L) {
        radius_to_use <- c(0, 0.3, 0.4)
      }
      for (main_radius in radius_to_use) {
        for (do_main_data_points in c("N", "Y")) {
          for (batch_size in c(1L, 10L, 100L, 1000L)) {
            use_main_data_points <- NULL
            if (do_main_data_points == "Y") use_main_data_points <- main_data_points
            use_main_radius <- main_radius
            if (use_main_radius == 0) use_main_radius <- NULL

            cat(c(size_constraint,
                  main_unassigned_method,
                  use_main_radius,
                  is.null(use_main_data_points),
                  batch_size), "\n")
            test_nng_batch_against_replica(test_distances1,
                                           size_constraint,
                                           main_unassigned_method,
                                           use_main_radius,
                                           use_main_data_points,
                                           batch_size)
          }
        }
      }
    }
  }
})
