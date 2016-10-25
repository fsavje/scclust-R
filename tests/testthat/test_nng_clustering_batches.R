library(Rscclust)
context("nng_clustering.R (batches)")

source("config.R", local = TRUE)
source("../replica/replica_make_nng.R", local = TRUE)
source("../replica/replica_nng_batches.R", local = TRUE)
source("utils_nng.R", local = TRUE)

test_that("batch nng clustering function cluster correctly.", {

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
                                 "by_nng",
                                 NULL,
                                 NULL,
                                 10L)
  test_nng_batch_against_replica(test_distances1,
                                 2L,
                                 "ignore",
                                 NULL,
                                 main_data_points,
                                 10L)
  test_nng_batch_against_replica(test_distances1,
                                 3L,
                                 "ignore",
                                 NULL,
                                 main_data_points,
                                 100L)
  test_nng_batch_against_replica(test_distances1,
                                 6L,
                                 "ignore",
                                 NULL,
                                 main_data_points,
                                 10L)
  test_nng_batch_against_replica(test_distances1,
                                 10L,
                                 "ignore",
                                 NULL,
                                 main_data_points,
                                 100L)
  test_nng_batch_against_replica(test_distances1,
                                 2L,
                                 "by_nng",
                                 NULL,
                                 main_data_points,
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
                                 "by_nng",
                                 test_radius,
                                 NULL,
                                 10L)
  test_nng_batch_against_replica(test_distances1,
                                 2L,
                                 "ignore",
                                 test_radius,
                                 main_data_points,
                                 10L)
  test_nng_batch_against_replica(test_distances1,
                                 3L,
                                 "ignore",
                                 test_radius,
                                 main_data_points,
                                 100L)
  test_nng_batch_against_replica(test_distances1,
                                 2L,
                                 "by_nng",
                                 test_radius,
                                 main_data_points,
                                 10L)
})


test_that("non-type nng clustering function cluster correctly all combinations", {
  skip_on_cran()
  skip_if_not(run_slow_tests)

  for (main_unassigned_method in c("ignore", "by_nng")) {
    for (main_radius in c(0, 0.1, 0.2)) {
      for (do_main_data_points in c("N", "Y")) {
        for (batch_size in c(1L, 10L, 100L, 1000L)) {
          use_main_data_points <- NULL
          if (do_main_data_points == "Y") use_main_data_points <- main_data_points
          use_main_radius <- main_radius
          if (use_main_radius == 0) use_main_radius <- NULL

          #cat(c(main_unassigned_method,
          #      use_main_radius,
          #      is.null(use_main_data_points),
          #      batch_size), "\n")
          test_nng_batch_against_replica(test_distances1,
                                         2L,
                                         main_unassigned_method,
                                         use_main_radius,
                                         use_main_data_points,
                                         batch_size)
        }
      }
    }
  }
})
