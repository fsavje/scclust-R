library(Rscclust)
context("misc.R")

distance_obj <- make_distances(matrix(c(0.1, 0.2, 0.3, 0.1, 0.2, 0.3), ncol = 2))
distance_obj_fake1 <- structure(c(0.1, 0.2, 0.3),
                                ids = NULL,
                                normalization = diag(3),
                                weights = diag(3),
                                class = c("Rscc_distances"))
distance_obj_fake2 <- structure(matrix(letters[1:9], ncol = 3),
                                ids = NULL,
                                normalization = diag(3),
                                weights = diag(3),
                                class = c("Rscc_distances"))

test_that("`check_Rscc_distances` checks input.", {
  expect_silent(check_Rscc_distances(distance_obj))
  expect_error(check_Rscc_distances(1:100))
  expect_error(check_Rscc_distances(distance_obj_fake1))
  expect_error(check_Rscc_distances(distance_obj_fake2))
})


cl1 <- Rscc_clustering(c("a", "b", "c", "a", "b", "c", "a"))
cl_fake1 <- make_Rscc_clustering(c("a", "b", "c", "a", "b", "c", "a"), 3, NULL)
cl_fake2 <- make_Rscc_clustering(integer(), 3, NULL)
cl_fake3 <- structure(c(0L, 1L, 2L, 0L, 1L, 2L, 0L),
                      ids = NULL,
                      class = c("Rscc_clustering"))
cl_fake4 <- make_Rscc_clustering(c(0L, 1L, 2L, 0L, 1L, 2L, 0L), 0, NULL)

test_that("`check_Rscc_clustering` checks input.", {
  expect_silent(check_Rscc_clustering(cl1))
  expect_error(check_Rscc_clustering(1:100))
  expect_error(check_Rscc_clustering(cl_fake1))
  expect_error(check_Rscc_clustering(cl_fake2))
  expect_error(check_Rscc_clustering(cl_fake3))
  expect_error(check_Rscc_clustering(cl_fake4))
})


test_that("`check_main_data_points` checks input.", {
  expect_silent(check_main_data_points(NULL, 4))
  expect_silent(check_main_data_points(c(TRUE, FALSE, TRUE, TRUE), 4))
  expect_error(check_main_data_points(1:100, 4))
  expect_error(check_main_data_points(c(TRUE, FALSE, TRUE, TRUE), 5))
})


type_label1 <- factor(c("a", "b", "c", "a", "a", "b", "c", "a", "b", "c", "a", "a"))
type_label2 <- c(0L, 1L, 2L, 2L, 1L, 0L, 2L, 1L, 2L)
type_label_fake1 <- factor(c("a", "b", "c", "a", NA, "b", "c", "a", NA, "c", "a", "a"))
type_label_fake2 <- c(0L, 1L, NA, 2L, 1L, NA, 2L, 1L, 2L)
type_label_fake3 <- c(0L, 1L, 2L, -2L, 1L, 0L, 2L, 1L, 2L)

test_that("`check_type_labels` checks input.", {
  expect_silent(check_type_labels(type_label1))
  expect_silent(check_type_labels(type_label2))
  expect_error(check_type_labels(c("a", "b", "c", "a", "a", "b", "c", "a", "b", "c", "a", "a")))
  expect_error(check_type_labels(type_label_fake1))
  expect_error(check_type_labels(type_label_fake2))
  expect_error(check_type_labels(type_label_fake3))
})


distance_obj1 <- make_distances(matrix(c(0.1, 0.2, 0.3, 0.1, 0.2, 0.3), ncol = 2))
distance_obj2 <- make_distances(matrix(c(0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.3), ncol = 2))

test_that("`get_num_data_points` produces correct output.", {
  expect_identical(get_num_data_points(distance_obj1), 3L)
  expect_identical(get_num_data_points(distance_obj2), 5L)
})


test_that("`get_size_constraint` produces correct output.", {
  expect_identical(get_size_constraint(4), 4L)
  expect_identical(get_size_constraint(100L), 100L)
  expect_error(get_size_constraint(1))
  expect_error(get_size_constraint(0))
  expect_error(get_size_constraint(-100L))
  expect_error(get_size_constraint("a"))
})


test_that("`get_seed_method` produces correct output.", {
  expect_match(get_seed_method(NULL), "lexical", fixed = TRUE)
  expect_match(get_seed_method("lexical"), "lexical", fixed = TRUE)
  expect_match(get_seed_method("lex"), "lexical", fixed = TRUE)
  expect_match(get_seed_method("inwards_order"), "inwards_order", fixed = TRUE)
  expect_match(get_seed_method("inwards_o"), "inwards_order", fixed = TRUE)
  expect_match(get_seed_method("inwards_updating"), "inwards_updating", fixed = TRUE)
  expect_match(get_seed_method("inwards_u"), "inwards_updating", fixed = TRUE)
  expect_match(get_seed_method("exclusion_order"), "exclusion_order", fixed = TRUE)
  expect_match(get_seed_method("exclusion_o"), "exclusion_order", fixed = TRUE)
  expect_match(get_seed_method("exclusion_updating"), "exclusion_updating", fixed = TRUE)
  expect_match(get_seed_method("exclusion_u"), "exclusion_updating", fixed = TRUE)
  expect_error(get_seed_method("inwards"))
  expect_error(get_seed_method("exclusion"))
  expect_error(get_seed_method("exclusion_nonexist"))
  expect_error(get_seed_method("abc"))
  expect_error(get_seed_method(1))
})


test_that("`get_radius` produces correct output.", {
  expect_identical(get_radius(NULL), NULL)
  expect_identical(get_radius(0.4), 0.4)
  expect_error(get_radius("a"))
  expect_error(get_radius(-0.4))
})


cl1 <- Rscc_clustering(c("a", "b", "c", "a", "b", "c", "a"))
cl2 <- Rscc_clustering(c("a", "b", "c", "a", "b", "c", "a", "d", "e", "d"))

test_that("`get_num_clusters` produces correct output.", {
  expect_identical(get_num_clusters(cl1), 3L)
  expect_identical(get_num_clusters(cl2), 5L)
})


factor1 <- factor(c("a", "b", "c", "a", "a", "b", "c", "a", "b", "c", "a", "a"))
factor2 <- factor(c("a", "b", "c", "a", " ", "a", " ", "b", "c", " ", "a", "b", "c", "a", "a", " "))
int1 <- c(0L, 1L, 2L, 2L, 1L, 0L, 2L, 1L, 2L)
int2 <- c(3L, 1L, 2L, 2L, 1L, 3L, 2L, 1L, 2L)
int3 <- c(3L, 5L, 2L, 2L, 5L, 3L, 2L, 5L, 2L)

test_that("`get_type_size_constraints` produces correct output.", {
  expect_equal(get_type_size_constraints(c(a = 1, b = 2, c = 3), factor1),
               c(0L, 1L, 2L, 3L))
  expect_equal(get_type_size_constraints(c(a = 0, b = 2, c = 3), factor1),
               c(0L, 0L, 2L, 3L))
  expect_equal(get_type_size_constraints(c(b = 2), factor1),
               c(0L, 0L, 2L, 0L))
  expect_equal(get_type_size_constraints(c(b = 2, " " = 3), factor2),
               c(0L, 3L, 0L, 2L, 0L))
  expect_equal(get_type_size_constraints(c("0" = 5, "1" = 6, "2" = 7), int1),
               c(5L, 6L, 7L))
  expect_equal(get_type_size_constraints(c("0" = 5, "2" = 7), int1),
               c(5L, 0L, 7L))
  expect_equal(get_type_size_constraints(c("1" = 5, "2" = 7), int2),
               c(0L, 5L, 7L, 0L))
  expect_equal(get_type_size_constraints(c("0" = 4, "1" = 5, "2" = 7), int2),
               c(4L, 5L, 7L, 0L))
  expect_equal(get_type_size_constraints(c("2" = 4, "3" = 5), int3),
               c(0L, 0L, 4L, 5L, 0L, 0L))
  expect_error(get_type_size_constraints(c(5, 7), factor1))
  expect_error(get_type_size_constraints(c("a" = 5, "b" = 7, "a" = 4), factor1))
  expect_error(get_type_size_constraints(c("a" = 5, "b" = 7, "none" = 4), factor1))
  expect_error(get_type_size_constraints(c("0" = 5, "1" = 6, "77" = 7), int1))
  expect_error(get_type_size_constraints(c(a = 1, b = -2, c = 3), factor1))
})


test_that("`get_total_size_constraint` produces correct output.", {
  expect_identical(get_total_size_constraint(NULL, 1L:3L), 6L)
  expect_identical(get_total_size_constraint(8L, 1L:3L), 8L)
  expect_error(get_total_size_constraint(NULL, c(0L, 0L, 0L)))
  expect_error(get_total_size_constraint(2L, 1L:3L))
})
