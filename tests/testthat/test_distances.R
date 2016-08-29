library(Rscclust)
context("distance.R")

cov1 <- c(-0.152, -0.619, 0.231, -0.524, 0.538, 0.414, 0.842, 0.173, -0.897, -0.963, 2.3, 0.206, -0.417, 0.462, 0.139, 1.262, -1.564, -0.385, 0.525, 0.158, 1.103, 0.793, -0.901, 0.744, -1.254, 0.045, -0.891, -0.154, 0.872, -0.118, 1.307, 0.242, -1.242, -0.47, 1.106, 1.126, 0.258, 0.943, 0.89, -0.869, 1.929, 0.115, -0.046, -0.162, 0.776, 0.903, -0.273, -0.982, 1.175, -0.519, -0.296, 0.709, -0.51, 0.969, -1.902, -0.936, -0.738, 0.582, 0.023, 0.578, 0.295, -0.005, -1.268, -0.038, 1.023, -0.733, -1.469, 0.545, 1.1, 0.962, 0.292, 0.317, -0.066, -0.394, -1.008, 1.621, 0.725, 0.963, -2.635, 0.092, -1.744, 1.083, 0.299, 0.605, 1.324, 0.083, 0.467, 0.288, 0.52, 0.258, -0.208, 0.207, -1.501, 0.54, -2.589, -0.276, -0.825, 0.67, 0.173, 2.45)
cov2 <- c(1.354, 1.321, 0.775, 0.779, 0.192, -0.027, -0.347, 0.472, 0.551, 0.797, -0.892, -0.684, 0.99, 0.525, -1.696, -0.063, 0.512, 0.023, 1.059, -0.786, -0.638, 0.35, 0.418, 0.767, -1.044, 1.962, -1.504, 1.19, 0.485, 0.464, -0.139, -1.571, 0.25, -0.228, 1.435, -2.201, 1.218, 2.444, 0.746, 0.523, 0.039, 1.422, -0.63, 0.46, 0.086, 1.242, 1.576, 0.313, -1.243, 0.706, 0.712, -0.066, 0.249, 0.359, 0.12, -0.872, -0.609, 1.115, 0.558, -1.743, 0.096, -0.336, 1.372, 0.298, 0.026, 1.012, 0.211, 1.058, 0.566, -1.473, -0.257, 0.835, 1.039, -0.884, -0.923, -0.066, 0.074, -1.923, 0.653, -0.636, 1.089, 1.433, -1.732, -0.544, -1.44, 2.669, 0.755, 0.206, -0.927, 0.046, -1.917, 0.256, -0.469, 0.312, 0.058, -0.202, -0.541, -0.942, 1.569, -1.76)
cov3 <- c(-1.134, 0.499, -0.551, 2.068, -0.004, -0.648, -1.139, 0.707, 0.666, -0.151, -0.193, 0.29, -0.113, 0.225, 0.161, 0.201, -0.753, 2.078, 0.507, -1.478, 1.126, -0.064, 0.437, -0.871, -1.364, 0.326, -0.807, -0.671, -0.424, -0.21, -0.771, 0.056, -0.022, 1.07, 1.072, 0.795, 1.751, -0.106, 0.05, -0.932, -1.04, -0.95, -0.339, -0.783, -1.276, 0.86, 1.303, -0.572, 0.976, 0.976, 0.742, 0.353, -0.582, -0.545, 0.154, -0.244, -1.32, -0.527, 0.805, -2.125, 0.23, 0.17, 0.971, -1.719, 0.543, -0.057, 0.487, 0.7, 0.091, -0.41, 0.509, -0.564, 0.149, 0.557, -0.061, 0.306, -1.357, 0.624, 0.142, -0.462, -0.992, 0.48, -1.247, -0.688, 0.036, -0.669, -0.123, 0.367, 1.649, 0.02, -0.484, -0.547, 0.969, 2.443, -0.465, 1.632, 0.067, -0.826, -0.29, 0.771)
idvar <- paste0("a", 1:100)

test_data_matrix <- matrix(c(cov1, cov2, cov3), nrow = 100)
test_data_matrix_single <- matrix(cov1, nrow = 100)
test_data_df1 <- data.frame(cov1 = cov1, cov2 = cov2, cov3 = cov3)
test_data_df2 <- data.frame(cov1 = cov1, cov2 = cov2, cov3 = cov3, extracol = 101:200)
test_data_df3 <- data.frame(idcol = idvar, cov1 = cov1, cov2 = cov2, cov3 = cov3, stringsAsFactors = FALSE)
test_data_df4 <- data.frame(idcol = idvar, cov1 = cov1, cov2 = cov2, cov3 = cov3, extracol = 101:200, stringsAsFactors = FALSE)
test_data_df1_single <- data.frame(cov1 = cov1)
test_data_df2_single <- data.frame(cov1 = cov1, extracol = 101:200)
test_data_df3_single <- data.frame(idcol = idvar, cov1 = cov1, stringsAsFactors = FALSE)
test_data_df4_single <- data.frame(idcol = idvar, cov1 = cov1, extracol = 101:200, stringsAsFactors = FALSE)

test_data_matrix_wNA <- test_data_matrix
test_data_matrix_wNA[55, 2] <- NA
test_data_df1_wNA <- test_data_df1
test_data_df1_wNA$cov2[34] <- NA


ref_out_vanilla <- structure(t(test_data_matrix),
                             normalization = diag(rep(1, 3)),
                             weights = diag(rep(1, 3)),
                             class = c("Rscc_distances"))
ref_out_vanilla_single <- structure(t(test_data_matrix_single),
                                    normalization = diag(1),
                                    weights = diag(1),
                                    class = c("Rscc_distances"))
ref_out_ids <- structure(t(test_data_matrix),
                         ids = idvar,
                         normalization = diag(rep(1, 3)),
                         weights = diag(rep(1, 3)),
                         class = c("Rscc_distances"))
ref_out_ids_single <- structure(t(test_data_matrix_single),
                                ids = idvar,
                                normalization = diag(1),
                                weights = diag(1),
                                class = c("Rscc_distances"))

ref_dist_mat_simple <- as.matrix(dist(test_data_matrix))
ref_dist_mat_simple_single <- as.matrix(dist(test_data_matrix_single))


test_that("`make_distances` checks data input.", {
  expect_is(make_distances(test_data_matrix), "Rscc_distances")
  expect_error(make_distances(NULL))
  expect_error(make_distances(1:100))
  expect_error(make_distances(matrix(1:100, nrow = 1)))
})


test_that("`make_distances` prints warning text.", {
  expect_message(print(make_distances(test_data_matrix)), "This is a Rscc_distances object. It's only purpose is to be used in the Rscclust package. Use `as.matrix` to generate an R matrix with the distances (a very slow operation for large distance matrices).\n", fixed = TRUE)
})


test_that("`make_distances` accepts matrix input.", {
  expect_is(make_distances(test_data_matrix), "Rscc_distances")
  expect_equal(make_distances(test_data_matrix), ref_out_vanilla)
  expect_equal(as.matrix(make_distances(test_data_matrix)), ref_dist_mat_simple)
  expect_error(make_distances(test_data_matrix_wNA))
  expect_is(make_distances(test_data_matrix, id_variable = idvar), "Rscc_distances")
  expect_equal(make_distances(test_data_matrix, id_variable = idvar), ref_out_ids)
  expect_equal(unname(as.matrix(make_distances(test_data_matrix, id_variable = idvar))), unname(ref_dist_mat_simple))
  expect_error(make_distances(test_data_matrix, id_variable = idvar[1:50]))
  expect_error(make_distances(test_data_matrix, id_variable = "idcol"))
  expect_error(make_distances(test_data_matrix, dist_variables = c("cov1", "cov2", "cov3")))
})


test_that("`make_distances` accepts 1D matrix input.", {
  expect_is(make_distances(test_data_matrix_single), "Rscc_distances")
  expect_equal(make_distances(test_data_matrix_single), ref_out_vanilla_single)
  expect_equal(as.matrix(make_distances(test_data_matrix_single)), ref_dist_mat_simple_single)
  expect_is(make_distances(test_data_matrix_single, id_variable = idvar), "Rscc_distances")
  expect_equal(make_distances(test_data_matrix_single, id_variable = idvar), ref_out_ids_single)
  expect_equal(unname(as.matrix(make_distances(test_data_matrix_single, id_variable = idvar))), unname(ref_dist_mat_simple_single))
  expect_error(make_distances(test_data_matrix_single, id_variable = idvar[1:50]))
  expect_error(make_distances(test_data_matrix_single, id_variable = "idcol"))
  expect_error(make_distances(test_data_matrix_single, dist_variables = c("cov1", "cov2", "cov3")))
})


test_that("`make_distances` accepts data.frame input.", {
  expect_is(make_distances(test_data_df1), "Rscc_distances")
  expect_equal(make_distances(test_data_df1), ref_out_vanilla)
  expect_equal(as.matrix(make_distances(test_data_df1)), ref_dist_mat_simple)
  expect_error(make_distances(test_data_df1_wNA))
  expect_is(make_distances(test_data_df1, id_variable = idvar), "Rscc_distances")
  expect_equal(make_distances(test_data_df1, id_variable = idvar), ref_out_ids)
  expect_equal(unname(as.matrix(make_distances(test_data_df1, id_variable = idvar))), unname(ref_dist_mat_simple))
  expect_error(make_distances(test_data_df1, id_variable = idvar[1:50]))

  expect_is(make_distances(test_data_df3, id_variable = "idcol"), "Rscc_distances")
  expect_equal(make_distances(test_data_df3, id_variable = "idcol"), ref_out_ids)
  expect_equal(unname(as.matrix(make_distances(test_data_df3, id_variable = "idcol"))), unname(ref_dist_mat_simple))
  expect_error(make_distances(test_data_df3, id_variable = "nonexisting"))
  expect_error(make_distances(test_data_df3, id_variable = 6L))

  expect_is(make_distances(test_data_df2, dist_variables = c("cov1", "cov2", "cov3")), "Rscc_distances")
  expect_equal(make_distances(test_data_df2, dist_variables = c("cov1", "cov2", "cov3")), ref_out_vanilla)
  expect_equal(as.matrix(make_distances(test_data_df2, dist_variables = c("cov1", "cov2", "cov3"))), ref_dist_mat_simple)
  expect_error(make_distances(test_data_df2, dist_variables = c("cov1", "nonexisting", "cov3")))
  expect_error(make_distances(test_data_df2, dist_variables = 1:3))

  expect_is(make_distances(test_data_df4, id_variable = "idcol", dist_variables = c("cov1", "cov2", "cov3")), "Rscc_distances")
  expect_equal(make_distances(test_data_df4, id_variable = "idcol", dist_variables = c("cov1", "cov2", "cov3")), ref_out_ids)
  expect_equal(unname(as.matrix(make_distances(test_data_df4, id_variable = "idcol", dist_variables = c("cov1", "cov2", "cov3")))), unname(ref_dist_mat_simple))
  expect_error(make_distances(test_data_df4, id_variable = "idcol", dist_variables = c("cov1", "nonexisting", "cov3")))
  expect_error(make_distances(test_data_df4, id_variable = "idcol", dist_variables = 1:3))
  expect_error(make_distances(test_data_df4, id_variable = "nonexisting", dist_variables = c("cov1", "cov2", "cov3")))
  expect_error(make_distances(test_data_df4, id_variable = 6L, dist_variables = c("cov1", "cov2", "cov3")))

  expect_equal(unname(as.matrix(make_distances(test_data_df4, id_variable = "cov1", dist_variables = c("cov1", "cov2", "cov3")))), unname(ref_dist_mat_simple))
})


test_that("`make_distances` accepts 1D data.frame input.", {
  expect_is(make_distances(test_data_df1_single), "Rscc_distances")
  expect_equal(make_distances(test_data_df1_single), ref_out_vanilla_single)
  expect_equal(as.matrix(make_distances(test_data_df1_single)), ref_dist_mat_simple_single)
  expect_is(make_distances(test_data_df1_single, id_variable = idvar), "Rscc_distances")
  expect_equal(make_distances(test_data_df1_single, id_variable = idvar), ref_out_ids_single)
  expect_equal(unname(as.matrix(make_distances(test_data_df1_single, id_variable = idvar))), unname(ref_dist_mat_simple_single))
  expect_error(make_distances(test_data_df1_single, id_variable = idvar[1:50]))

  expect_is(make_distances(test_data_df3_single, id_variable = "idcol"), "Rscc_distances")
  expect_equal(make_distances(test_data_df3_single, id_variable = "idcol"), ref_out_ids_single)
  expect_equal(unname(as.matrix(make_distances(test_data_df3_single, id_variable = "idcol"))), unname(ref_dist_mat_simple_single))
  expect_error(make_distances(test_data_df3_single, id_variable = "nonexisting"))
  expect_error(make_distances(test_data_df3_single, id_variable = 6L))

  expect_is(make_distances(test_data_df2_single, dist_variables = c("cov1")), "Rscc_distances")
  expect_equal(make_distances(test_data_df2_single, dist_variables = c("cov1")), ref_out_vanilla_single)
  expect_equal(as.matrix(make_distances(test_data_df2_single, dist_variables = c("cov1"))), ref_dist_mat_simple_single)
  expect_error(make_distances(test_data_df2_single, dist_variables = c("nonexisting")))
  expect_error(make_distances(test_data_df2_single, dist_variables = 1L))

  expect_is(make_distances(test_data_df4_single, id_variable = "idcol", dist_variables = c("cov1")), "Rscc_distances")
  expect_equal(make_distances(test_data_df4_single, id_variable = "idcol", dist_variables = c("cov1")), ref_out_ids_single)
  expect_equal(unname(as.matrix(make_distances(test_data_df4_single, id_variable = "idcol", dist_variables = c("cov1")))), unname(ref_dist_mat_simple_single))
  expect_error(make_distances(test_data_df4_single, id_variable = "idcol", dist_variables = c("nonexisting")))
  expect_error(make_distances(test_data_df4_single, id_variable = "idcol", dist_variables = 1L))
  expect_error(make_distances(test_data_df4_single, id_variable = "nonexisting", dist_variables = c("cov1")))
  expect_error(make_distances(test_data_df4_single, id_variable = 6L, dist_variables = c("cov1")))

  expect_equal(unname(as.matrix(make_distances(test_data_df4_single, id_variable = "cov1", dist_variables = c("cov1")))), unname(ref_dist_mat_simple_single))
})


ref_make_distance <- function(mat, covmat, inverted = FALSE) {
  tmp <- sqrt(apply(mat, 1, function(x) { mahalanobis(mat, x, covmat, inverted) }))
  dimnames(tmp) <- list(as.character(1:nrow(mat)), as.character(1:nrow(mat)))
  tmp
}

ref_dist_mat_simple2 <- ref_make_distance(test_data_matrix, diag(3))
ref_dist_mat_mahalanobis <- ref_make_distance(test_data_matrix, var(test_data_matrix))
ref_dist_mat_student <- ref_make_distance(test_data_matrix, diag(diag(var(test_data_matrix))))
ref_dist_mat_custom <- ref_make_distance(test_data_matrix, diag(c(1, 2, 3)))

test_that("`normalize` works.", {
  expect_is(make_distances(test_data_matrix), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix), "normalization"), diag(3))
  expect_equal(as.matrix(make_distances(test_data_matrix)), ref_dist_mat_simple2)

  expect_is(make_distances(test_data_matrix, normalize = NULL), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, normalize = NULL), "normalization"), diag(3))
  expect_equal(as.matrix(make_distances(test_data_matrix, normalize = NULL)), ref_dist_mat_simple2)

  expect_is(make_distances(test_data_matrix, normalize = "none"), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, normalize = "none"), "normalization"), diag(3))
  expect_equal(as.matrix(make_distances(test_data_matrix, normalize = "none")), ref_dist_mat_simple2)

  expect_is(make_distances(test_data_matrix, normalize = "no"), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, normalize = "no"), "normalization"), diag(3))
  expect_equal(as.matrix(make_distances(test_data_matrix, normalize = "no")), ref_dist_mat_simple2)

  expect_is(make_distances(test_data_matrix, normalize = "mahalanobis"), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, normalize = "mahalanobis"), "normalization"), var(test_data_matrix))
  expect_equal(as.matrix(make_distances(test_data_matrix, normalize = "mahalanobis")), ref_dist_mat_mahalanobis)

  expect_is(make_distances(test_data_matrix, normalize = "mahalanobize"), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, normalize = "mahalanobize"), "normalization"), var(test_data_matrix))
  expect_equal(as.matrix(make_distances(test_data_matrix, normalize = "mahalanobize")), ref_dist_mat_mahalanobis)

  expect_is(make_distances(test_data_matrix, normalize = "maha"), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, normalize = "maha"), "normalization"), var(test_data_matrix))
  expect_equal(as.matrix(make_distances(test_data_matrix, normalize = "maha")), ref_dist_mat_mahalanobis)

  expect_is(make_distances(test_data_matrix, normalize = "studentize"), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, normalize = "studentize"), "normalization"), diag(diag(var(test_data_matrix))))
  expect_equal(as.matrix(make_distances(test_data_matrix, normalize = "studentize")), ref_dist_mat_student)

  expect_is(make_distances(test_data_matrix, normalize = "student"), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, normalize = "student"), "normalization"), diag(diag(var(test_data_matrix))))
  expect_equal(as.matrix(make_distances(test_data_matrix, normalize = "student")), ref_dist_mat_student)

  expect_is(make_distances(test_data_matrix, normalize = matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3)), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, normalize = matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3)), "normalization"),
               matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3))
  expect_equal(as.matrix(make_distances(test_data_matrix, normalize = matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3))), ref_dist_mat_custom)

  expect_is(make_distances(test_data_matrix, normalize = data.frame(matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3))), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, normalize = data.frame(matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3))), "normalization"),
               matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3))
  expect_equal(as.matrix(make_distances(test_data_matrix, normalize = data.frame(matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3)))), ref_dist_mat_custom)

  expect_is(make_distances(test_data_matrix, normalize = c(1, 2, 3)), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, normalize = c(1, 2, 3)), "normalization"),
               matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3))
  expect_equal(as.matrix(make_distances(test_data_matrix, normalize = c(1, 2, 3))), ref_dist_mat_custom)

  expect_error(make_distances(test_data_matrix, normalize = "abc"))
  expect_error(make_distances(test_data_matrix, normalize = dist(1:4)))
  expect_error(make_distances(test_data_matrix, normalize = matrix(letters[1:9], nrow = 3)))
  expect_error(make_distances(test_data_matrix, normalize = matrix(c(1, 0, NA, 0, 2, 0, 0, 0, 3), nrow = 3)))
  expect_error(make_distances(test_data_matrix, normalize = matrix(c(1, 0, 1, 0, 2, 0, 0, 0, 3), nrow = 3)))
  expect_error(make_distances(test_data_matrix, normalize = matrix(c(1, 2, 0, 2, 1, 2, 0, 2, 1), nrow = 3)))
  expect_error(make_distances(test_data_matrix, normalize = matrix(c(1, 0, 0, 2), nrow = 2)))
})


ref_dist_mat_weights <- ref_make_distance(test_data_matrix, diag(c(1, 2, 3)), inverted = TRUE)

test_that("`weights` works.", {
  expect_is(make_distances(test_data_matrix), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix), "weights"), diag(3))
  expect_equal(as.matrix(make_distances(test_data_matrix)), ref_dist_mat_simple2)

  expect_is(make_distances(test_data_matrix, weights = NULL), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, weights = NULL), "weights"), diag(3))
  expect_equal(as.matrix(make_distances(test_data_matrix, weights = NULL)), ref_dist_mat_simple2)

  expect_is(make_distances(test_data_matrix, weights = matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3)), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, weights = matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3)), "weights"),
               matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3))
  expect_equal(as.matrix(make_distances(test_data_matrix, weights = matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3))), ref_dist_mat_weights)

  expect_is(make_distances(test_data_matrix, weights = data.frame(matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3))), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, weights = data.frame(matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3))), "weights"),
               matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3))
  expect_equal(as.matrix(make_distances(test_data_matrix, weights = data.frame(matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3)))), ref_dist_mat_weights)

  expect_is(make_distances(test_data_matrix, weights = c(1, 2, 3)), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, weights = c(1, 2, 3)), "weights"),
               matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3))
  expect_equal(as.matrix(make_distances(test_data_matrix, weights = c(1, 2, 3))), ref_dist_mat_weights)

  expect_error(make_distances(test_data_matrix, weights = dist(1:4)))
  expect_error(make_distances(test_data_matrix, weights = matrix(letters[1:9], nrow = 3)))
  expect_error(make_distances(test_data_matrix, weights = matrix(c(1, 0, NA, 0, 2, 0, 0, 0, 3), nrow = 3)))
  expect_error(make_distances(test_data_matrix, weights = matrix(c(1, 0, 1, 0, 2, 0, 0, 0, 3), nrow = 3)))
  expect_error(make_distances(test_data_matrix, weights = matrix(c(1, 2, 0, 2, 1, 2, 0, 2, 1), nrow = 3)))
  expect_error(make_distances(test_data_matrix, weights = matrix(c(1, 0, 0, 2), nrow = 2)))
})


replica_make_distances <- function(mat, normalize, weights) {
  covmat <- t(chol(solve(normalize))) %*% weights %*% chol(solve(normalize))
  tmp <- sqrt(apply(mat, 1, function(x) { mahalanobis(mat, x, covmat, TRUE) }))
  dimnames(tmp) <- list(as.character(1:nrow(mat)), as.character(1:nrow(mat)))
  tmp
}

ref_dist_mat_simple2_wweights <- replica_make_distances(test_data_matrix, diag(3), diag(c(4, 5, 6)))
ref_dist_mat_mahalanobis_wweights <- replica_make_distances(test_data_matrix, var(test_data_matrix), diag(c(4, 5, 6)))
ref_dist_mat_student_wweights <- replica_make_distances(test_data_matrix, diag(diag(var(test_data_matrix))), diag(c(4, 5, 6)))
ref_dist_mat_custom_wweights <- replica_make_distances(test_data_matrix, diag(c(1, 2, 3)), diag(c(4, 5, 6)))

test_that("`normalize` and `weights` works.", {
  expect_is(make_distances(test_data_matrix, normalize = "none", weights = c(4, 5, 6)), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, normalize = "none", weights = c(4, 5, 6)), "normalization"), diag(3))
  expect_equal(attr(make_distances(test_data_matrix, normalize = "none", weights = c(4, 5, 6)), "weights"), diag(c(4, 5, 6)))
  expect_equal(as.matrix(make_distances(test_data_matrix, normalize = "none", weights = c(4, 5, 6))), ref_dist_mat_simple2_wweights)

  expect_is(make_distances(test_data_matrix, normalize = "mahalanobize", weights = c(4, 5, 6)), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, normalize = "mahalanobize", weights = c(4, 5, 6)), "normalization"), var(test_data_matrix))
  expect_equal(attr(make_distances(test_data_matrix, normalize = "mahalanobize", weights = c(4, 5, 6)), "weights"), diag(c(4, 5, 6)))
  expect_equal(as.matrix(make_distances(test_data_matrix, normalize = "mahalanobize", weights = c(4, 5, 6))), ref_dist_mat_mahalanobis_wweights)

  expect_is(make_distances(test_data_matrix, normalize = "studentize", weights = c(4, 5, 6)), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix, normalize = "studentize", weights = c(4, 5, 6)), "normalization"), diag(diag(var(test_data_matrix))))
  expect_equal(attr(make_distances(test_data_matrix, normalize = "studentize", weights = c(4, 5, 6)), "weights"), diag(c(4, 5, 6)))
  expect_equal(as.matrix(make_distances(test_data_matrix, normalize = "studentize", weights = c(4, 5, 6))), ref_dist_mat_student_wweights)

  expect_is(make_distances(test_data_matrix,
                           normalize = matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3),
                           weights = c(4, 5, 6)), "Rscc_distances")
  expect_equal(attr(make_distances(test_data_matrix,
                                   normalize = matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3),
                                   weights = c(4, 5, 6)), "normalization"), matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3))
  expect_equal(attr(make_distances(test_data_matrix,
                                   normalize = matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3),
                                   weights = c(4, 5, 6)), "weights"), diag(c(4, 5, 6)))
  expect_equal(as.matrix(make_distances(test_data_matrix,
                                        normalize = matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3),
                                        weights = c(4, 5, 6))), ref_dist_mat_custom_wweights)
})
