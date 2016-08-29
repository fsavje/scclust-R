library(Rscclust)
context("utilities.R")

sound_clustering <- Rscc_clustering(c("a", "b", "c", "a", "b", "c", "a"))
unsound_clustering <- make_Rscc_clustering(c("a", "b", "c", "a", "b", "c", "a"), 3, NULL)
sound_size_constraint <- 3
unsound_size_constraint <- 0
sound_type_labels <- factor(c("a", "b", "c", "a", "b", "c", "a"))
unsound_type_labels <- c("a", "b", "c", "a", "b", "c", "a")
sound_type_size_constraints <- c("a" = 1, "b" = 2)
unsound_type_size_constraints <- c("a" = 1, "bbb" = 2)
sound_distance_obj <- make_distances(matrix(c(0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1), ncol = 1))
unsound_distance_obj <- structure(matrix(letters[1:9], ncol = 3),
                                  ids = NULL,
                                  normalization = diag(3),
                                  weights = diag(3),
                                  class = c("Rscc_distances"))

test_that("utility functions check input.", {
  expect_equal(check_clustering(sound_clustering, sound_size_constraint), FALSE)
  expect_error(check_clustering(unsound_clustering, sound_size_constraint))
  expect_error(check_clustering(sound_clustering, unsound_size_constraint))
  expect_equal(check_clustering_types(sound_clustering,
                                      sound_type_labels,
                                      sound_type_size_constraints,
                                      sound_size_constraint), FALSE)
  expect_error(check_clustering_types(unsound_clustering,
                                      sound_type_labels,
                                      sound_type_size_constraints,
                                      sound_size_constraint))
  expect_error(check_clustering_types(sound_clustering,
                                      unsound_type_labels,
                                      sound_type_size_constraints,
                                      sound_size_constraint))
  expect_error(check_clustering_types(sound_clustering,
                                      sound_type_labels,
                                      unsound_type_size_constraints,
                                      sound_size_constraint))
  expect_error(check_clustering_types(sound_clustering,
                                      sound_type_labels,
                                      sound_type_size_constraints,
                                      unsound_size_constraint))
  expect_is(get_clustering_stats(sound_clustering, sound_distance_obj), "Rscc_clustering_stats")
  expect_error(get_clustering_stats(unsound_clustering, sound_distance_obj))
  expect_error(get_clustering_stats(sound_clustering, unsound_distance_obj))
})



cl1 <- Rscc_clustering(c(1, 1, 2, 3, 2, 3, 3, 1, 2, 2))
cl2 <- Rscc_clustering(c(1, 1, 2, 3, 2, 3, 3, 1, 2, 2), 1)
cl3 <- Rscc_clustering(c(1, 1, 2, 3, 2, 3, 3, 1, 2, 2), ids = letters[1:10])
cl4 <- Rscc_clustering(c(1, 1, 2, 3, 2, 3, 3, 1, 2, 2), 1, ids = letters[1:10])

test_that("`check_clustering` returns correctly.", {
  expect_equal(check_clustering(cl1, 2), TRUE)
  expect_equal(check_clustering(cl1, 3), TRUE)
  expect_equal(check_clustering(cl1, 4), FALSE)
  expect_equal(check_clustering(cl2, 2), TRUE)
  expect_equal(check_clustering(cl2, 3), TRUE)
  expect_equal(check_clustering(cl2, 4), FALSE)
  expect_equal(check_clustering(cl3, 2), TRUE)
  expect_equal(check_clustering(cl3, 3), TRUE)
  expect_equal(check_clustering(cl3, 4), FALSE)
  expect_equal(check_clustering(cl4, 2), TRUE)
  expect_equal(check_clustering(cl4, 3), TRUE)
  expect_equal(check_clustering(cl4, 4), FALSE)
})


cl1 <- Rscc_clustering(c(1, 1, 2, 3, 2, 3, 3, 1, 2, 2))
cl2 <- Rscc_clustering(c(1, 1, 2, 3, 2, 3, 3, 1, 2, 2), 1)
cl3 <- Rscc_clustering(c(1, 1, 2, 3, 2, 3, 3, 1, 2, 2), ids = letters[1:10])
cl4 <- Rscc_clustering(c(1, 1, 2, 3, 2, 3, 3, 1, 2, 2), 1, ids = letters[1:10])
dp_types <- factor(c("x", "y", "y", "z", "z", "x", "y", "z", "x", "x"))

test_that("`check_clustering_types` returns correctly.", {
  expect_equal(check_clustering_types(cl1, dp_types, c("x" = 1, "y" = 1, "z" = 1)), TRUE)
  expect_equal(check_clustering_types(cl1, dp_types, c("x" = 1, "y" = 1, "z" = 1), 3), TRUE)
  expect_equal(check_clustering_types(cl1, dp_types, c("x" = 1, "z" = 1)), TRUE)
  expect_equal(check_clustering_types(cl1, dp_types, c("y" = 1, "z" = 1), 3), TRUE)
  expect_equal(check_clustering_types(cl1, dp_types, c("x" = 1, "y" = 1, "z" = 1), 4), FALSE)
  expect_equal(check_clustering_types(cl1, dp_types, c("x" = 3, "y" = 1, "z" = 1)), FALSE)

  expect_equal(check_clustering_types(cl2, dp_types, c("x" = 1, "y" = 1, "z" = 1)), TRUE)
  expect_equal(check_clustering_types(cl2, dp_types, c("x" = 1, "y" = 1, "z" = 1), 3), TRUE)
  expect_equal(check_clustering_types(cl2, dp_types, c("x" = 1, "z" = 1)), TRUE)
  expect_equal(check_clustering_types(cl2, dp_types, c("y" = 1, "z" = 1), 3), TRUE)
  expect_equal(check_clustering_types(cl2, dp_types, c("x" = 1, "y" = 1, "z" = 1), 4), FALSE)
  expect_equal(check_clustering_types(cl2, dp_types, c("x" = 3, "y" = 1, "z" = 1)), FALSE)

  expect_equal(check_clustering_types(cl3, dp_types, c("x" = 1, "y" = 1, "z" = 1)), TRUE)
  expect_equal(check_clustering_types(cl3, dp_types, c("x" = 1, "y" = 1, "z" = 1), 3), TRUE)
  expect_equal(check_clustering_types(cl3, dp_types, c("x" = 1, "z" = 1)), TRUE)
  expect_equal(check_clustering_types(cl3, dp_types, c("y" = 1, "z" = 1), 3), TRUE)
  expect_equal(check_clustering_types(cl3, dp_types, c("x" = 1, "y" = 1, "z" = 1), 4), FALSE)
  expect_equal(check_clustering_types(cl3, dp_types, c("x" = 3, "y" = 1, "z" = 1)), FALSE)

  expect_equal(check_clustering_types(cl4, dp_types, c("x" = 1, "y" = 1, "z" = 1)), TRUE)
  expect_equal(check_clustering_types(cl4, dp_types, c("x" = 1, "y" = 1, "z" = 1), 3), TRUE)
  expect_equal(check_clustering_types(cl4, dp_types, c("x" = 1, "z" = 1)), TRUE)
  expect_equal(check_clustering_types(cl4, dp_types, c("y" = 1, "z" = 1), 3), TRUE)
  expect_equal(check_clustering_types(cl4, dp_types, c("x" = 1, "y" = 1, "z" = 1), 4), FALSE)
  expect_equal(check_clustering_types(cl4, dp_types, c("x" = 3, "y" = 1, "z" = 1)), FALSE)
})


cl1 <- Rscc_clustering(c(1, 1, 2, 3, 2, 3, 3, 1, 2, 2))
cl2 <- Rscc_clustering(c(1, 1, 2, 3, 2, 3, 3, 1, 2, 2), 1)
cl3 <- Rscc_clustering(c(1, 1, 2, 3, 2, 3, 3, 1, 2, 2), ids = letters[1:10])
cl4 <- Rscc_clustering(c(1, 1, 2, 3, 2, 3, 3, 1, 2, 2), 1, ids = letters[1:10])
cl5 <- make_Rscc_clustering(c(1L, 1L, 3L, 4L, 3L, 4L, 4L, 1L, 3L, 3L),
                            5L, NULL)
cl6 <- make_Rscc_clustering(c(NA, NA, 3L, 4L, 3L, 4L, 4L, NA, 3L, 3L),
                            5L, NULL)
cl7 <- make_Rscc_clustering(c(1L, 1L, 3L, 4L, 3L, 4L, 4L, 1L, 3L, 3L),
                            5L, ids = letters[1:10])
cl8 <- make_Rscc_clustering(c(NA, NA, 3L, 4L, 3L, 4L, 4L, NA, 3L, 3L),
                            5L, ids = letters[1:10])
dp_data <- matrix(c(1.43247031975424, -0.0743079629128504, 0.664288699457021, 0.041257931075698, -0.745804833010271, -0.603411650921577, -0.705441137672159, 1.1906572821574, 1.25508544074284, -0.732099707199829, 1.18248873752832, 0.0593032947441447, 0.547530248329181, -1.30988291210263, -1.92018121606691, 0.212216066820278, -1.32505024286282, -1.27119614395428, -0.207545316189576, -1.51865610892823),
                  ncol = 2)

all_assigned_stats <- structure(list(
  num_data_points = 10L,
  num_assigned = 10L,
  num_clusters = 3L,
  num_populated_clusters = 3L,
  min_cluster_size = 3L,
  max_cluster_size = 4L,
  avg_cluster_size = 10/3,
  sum_dists = 21.8322314142849,
  min_dist = 0.401758935353289,
  max_dist = 2.84217586398484,
  cl_avg_min_dist = 0.994822547780713,
  cl_avg_max_dist = 2.32024663290364,
  cl_avg_dist_weighted = 1.79285753468519,
  cl_avg_dist_unweighted = 1.77519414590393
), class = c("Rscc_clustering_stats"))

some_assigned_stats <- structure(list(
  num_data_points = 10L,
  num_assigned = 7L,
  num_clusters = 2L,
  num_populated_clusters = 2L,
  min_cluster_size = 3,
  max_cluster_size = 4,
  avg_cluster_size = 3.5,
  sum_dists = 15.6514622631497,
  min_dist = 0.401758935353289,
  max_dist = 2.84217586398484,
  cl_avg_min_dist = 0.57430601573908,
  cl_avg_max_dist = 2.24758417693106,
  cl_avg_dist_weighted = 1.67825802795953,
  cl_avg_dist_unweighted = 1.63266302700004
), class = c("Rscc_clustering_stats"))

all_assigned_emptycl_stats <- structure(list(
  num_data_points = 10L,
  num_assigned = 10L,
  num_clusters = 5L,
  num_populated_clusters = 3L,
  min_cluster_size = 3L,
  max_cluster_size = 4L,
  avg_cluster_size = 10/3,
  sum_dists = 21.8322314142849,
  min_dist = 0.401758935353289,
  max_dist = 2.84217586398484,
  cl_avg_min_dist = 0.994822547780713,
  cl_avg_max_dist = 2.32024663290364,
  cl_avg_dist_weighted = 1.79285753468519,
  cl_avg_dist_unweighted = 1.77519414590393
), class = c("Rscc_clustering_stats"))

some_assigned_emptycl_stats <- structure(list(
  num_data_points = 10L,
  num_assigned = 7L,
  num_clusters = 5L,
  num_populated_clusters = 2L,
  min_cluster_size = 3,
  max_cluster_size = 4,
  avg_cluster_size = 3.5,
  sum_dists = 15.6514622631497,
  min_dist = 0.401758935353289,
  max_dist = 2.84217586398484,
  cl_avg_min_dist = 0.57430601573908,
  cl_avg_max_dist = 2.24758417693106,
  cl_avg_dist_weighted = 1.67825802795953,
  cl_avg_dist_unweighted = 1.63266302700004
), class = c("Rscc_clustering_stats"))

test_that("`check_clustering_types` returns correctly.", {
  expect_equal(get_clustering_stats(cl1, make_distances(dp_data)), all_assigned_stats)
  expect_equal(get_clustering_stats(cl2, make_distances(dp_data)), some_assigned_stats)
  expect_equal(get_clustering_stats(cl3, make_distances(dp_data)), all_assigned_stats)
  expect_equal(get_clustering_stats(cl4, make_distances(dp_data)), some_assigned_stats)
  expect_equal(get_clustering_stats(cl5, make_distances(dp_data)), all_assigned_emptycl_stats)
  expect_equal(get_clustering_stats(cl6, make_distances(dp_data)), some_assigned_emptycl_stats)
  expect_equal(get_clustering_stats(cl7, make_distances(dp_data)), all_assigned_emptycl_stats)
  expect_equal(get_clustering_stats(cl8, make_distances(dp_data)), some_assigned_emptycl_stats)
})
