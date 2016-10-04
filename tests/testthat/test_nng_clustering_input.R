library(Rscclust)
context("nng_clustering.R input check")

sound_distance_object <- make_distances(matrix(c(0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1), ncol = 1))
unsound_distance_object <- structure(matrix(letters[1:9], ncol = 3),
                                     ids = NULL,
                                     normalization = diag(3),
                                     weights = diag(3),
                                     class = c("Rscc_distances"))
sound_size_constraint <- 3
unsound_size_constraint <- NULL
sound_seed_method <- "lexical"
unsound_seed_method <- "dontexist"
sound_radius <- NULL
unsound_radius <- -1.0
sound_main_data_points <- NULL
sound_main_data_points2 <- c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE)
unsound_main_data_points <- "a"
sound_type_labels <- factor(c("a", "b", "c", "a", "b", "c", "a"))
unsound_type_labels <- c("a", "b", "c", "a", "b", "c", "a")
sound_type_size_constraints <- c("a" = 1, "b" = 2)
unsound_type_size_constraints <- c("a" = 1, "bbb" = 2)
sound_total_size_constraint <- NULL
unsound_total_size_constraint <- 0

test_that("nng clustering functions check input.", {
  expect_match(class(nng_clustering(sound_distance_object,
                                    sound_size_constraint,
                                    sound_seed_method,
                                    "closest_seed",
                                    sound_radius,
                                    sound_main_data_points,
                                    "ignore",
                                    sound_radius)), "Rscc_clustering", fixed = TRUE)
  expect_error(nng_clustering(unsound_distance_object,
                              sound_size_constraint,
                              sound_seed_method,
                              "closest_seed",
                              sound_radius,
                              sound_main_data_points,
                              "ignore",
                              sound_radius))
  expect_error(nng_clustering(sound_distance_object,
                              unsound_size_constraint,
                              sound_seed_method,
                              "closest_seed",
                              sound_radius,
                              sound_main_data_points,
                              "ignore",
                              sound_radius))
  expect_error(nng_clustering(sound_distance_object,
                              sound_size_constraint,
                              unsound_seed_method,
                              "closest_seed",
                              sound_radius,
                              sound_main_data_points,
                              "ignore",
                              sound_radius))
  expect_error(nng_clustering(sound_distance_object,
                              sound_size_constraint,
                              sound_seed_method,
                              "closest_seed",
                              unsound_radius,
                              sound_main_data_points,
                              "ignore",
                              sound_radius))
  expect_error(nng_clustering(sound_distance_object,
                              sound_size_constraint,
                              sound_seed_method,
                              "closest_seed",
                              sound_radius,
                              unsound_main_data_points,
                              "ignore",
                              sound_radius))
  expect_error(nng_clustering(sound_distance_object,
                              sound_size_constraint,
                              sound_seed_method,
                              "closest_seed",
                              sound_radius,
                              sound_main_data_points,
                              "ignore",
                              unsound_radius))
  expect_error(nng_clustering(sound_distance_object,
                              sound_size_constraint,
                              sound_seed_method,
                              "dontexist",
                              sound_radius,
                              sound_main_data_points,
                              "ignore",
                              sound_radius))
  expect_error(nng_clustering(sound_distance_object,
                              sound_size_constraint,
                              sound_seed_method,
                              "closest_seed",
                              sound_radius,
                              sound_main_data_points2,
                              "by_nng",
                              sound_radius))


  expect_match(class(nng_clustering_batches(sound_distance_object,
                                            sound_size_constraint,
                                            "ignore",
                                            sound_radius,
                                            sound_main_data_points,
                                            10)), "Rscc_clustering", fixed = TRUE)
  expect_error(nng_clustering_batches(unsound_distance_object,
                                      sound_size_constraint,
                                      "ignore",
                                      sound_radius,
                                      sound_main_data_points,
                                      10))
  expect_error(nng_clustering_batches(sound_distance_object,
                                      unsound_size_constraint,
                                      "ignore",
                                      sound_radius,
                                      sound_main_data_points,
                                      10))
  expect_error(nng_clustering_batches(sound_distance_object,
                                      sound_size_constraint,
                                      "ignore",
                                      unsound_radius,
                                      sound_main_data_points,
                                      10))
  expect_error(nng_clustering_batches(sound_distance_object,
                                      sound_size_constraint,
                                      "ignore",
                                      sound_radius,
                                      unsound_main_data_points,
                                      10))
  expect_error(nng_clustering_batches(sound_distance_object,
                                      sound_size_constraint,
                                      "closest_assigned",
                                      sound_radius,
                                      sound_main_data_points,
                                      10))
  expect_error(nng_clustering_batches(sound_distance_object,
                                      sound_size_constraint,
                                      "ignore",
                                      sound_radius,
                                      sound_main_data_points,
                                      -1))
  expect_error(nng_clustering_batches(sound_distance_object,
                                      sound_size_constraint,
                                      "ignore",
                                      sound_radius,
                                      sound_main_data_points,
                                      "a"))


  expect_match(class(nng_clustering_types(sound_distance_object,
                                          sound_type_labels,
                                          sound_type_size_constraints,
                                          sound_total_size_constraint,
                                          sound_seed_method,
                                          "closest_seed",
                                          sound_radius,
                                          sound_main_data_points,
                                          "ignore",
                                          sound_radius)), "Rscc_clustering", fixed = TRUE)
  expect_error(nng_clustering_types(unsound_distance_object,
                                    sound_type_labels,
                                    sound_type_size_constraints,
                                    sound_total_size_constraint,
                                    sound_seed_method,
                                    "closest_seed",
                                    sound_radius,
                                    sound_main_data_points,
                                    "ignore",
                                    sound_radius))
  expect_error(nng_clustering_types(sound_distance_object,
                                    unsound_type_labels,
                                    sound_type_size_constraints,
                                    sound_total_size_constraint,
                                    sound_seed_method,
                                    "closest_seed",
                                    sound_radius,
                                    sound_main_data_points,
                                    "ignore",
                                    sound_radius))
  expect_error(nng_clustering_types(sound_distance_object,
                                    sound_type_labels,
                                    unsound_type_size_constraints,
                                    sound_total_size_constraint,
                                    sound_seed_method,
                                    "closest_seed",
                                    sound_radius,
                                    sound_main_data_points,
                                    "ignore",
                                    sound_radius))
  expect_error(nng_clustering_types(sound_distance_object,
                                    sound_type_labels,
                                    sound_type_size_constraints,
                                    unsound_total_size_constraint,
                                    sound_seed_method,
                                    "closest_seed",
                                    sound_radius,
                                    sound_main_data_points,
                                    "ignore",
                                    sound_radius))
  expect_error(nng_clustering_types(sound_distance_object,
                                    sound_type_labels,
                                    sound_type_size_constraints,
                                    sound_total_size_constraint,
                                    unsound_seed_method,
                                    "closest_seed",
                                    sound_radius,
                                    sound_main_data_points,
                                    "ignore",
                                    sound_radius))
  expect_error(nng_clustering_types(sound_distance_object,
                                    sound_type_labels,
                                    sound_type_size_constraints,
                                    sound_total_size_constraint,
                                    sound_seed_method,
                                    "closest_seed",
                                    unsound_radius,
                                    sound_main_data_points,
                                    "ignore",
                                    sound_radius))
  expect_error(nng_clustering_types(sound_distance_object,
                                    sound_type_labels,
                                    sound_type_size_constraints,
                                    sound_total_size_constraint,
                                    sound_seed_method,
                                    "closest_seed",
                                    sound_radius,
                                    unsound_main_data_points,
                                    "ignore",
                                    sound_radius))
  expect_error(nng_clustering_types(sound_distance_object,
                                    sound_type_labels,
                                    sound_type_size_constraints,
                                    sound_total_size_constraint,
                                    sound_seed_method,
                                    "closest_seed",
                                    sound_radius,
                                    sound_main_data_points,
                                    "ignore",
                                    unsound_radius))
  expect_error(nng_clustering_types(sound_distance_object,
                                    sound_type_labels,
                                    sound_type_size_constraints,
                                    sound_total_size_constraint,
                                    sound_seed_method,
                                    "dontexist",
                                    sound_radius,
                                    sound_main_data_points,
                                    "ignore",
                                    sound_radius))
  expect_error(nng_clustering_types(sound_distance_object,
                                    sound_type_labels,
                                    sound_type_size_constraints,
                                    sound_total_size_constraint,
                                    sound_seed_method,
                                    "closest_seed",
                                    sound_radius,
                                    sound_main_data_points2,
                                    "by_nng",
                                    sound_radius))
})

test_that("C error functions work.", {
  expect_error(nng_clustering(make_distances(matrix(c(0.1, 0.2, 0.3), ncol = 1)), size_constraint = 2L, main_radius = 0.001),
               "No clustering satisfying the specified radius constraints exists.")
  expect_error(.Call("Rsccwrap_nng_clustering",
                     NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                     PACKAGE = "Rscclust"),
               "Invalid distance object.")
})
