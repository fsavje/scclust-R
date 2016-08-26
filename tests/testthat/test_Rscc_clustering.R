library(Rscclust)
context("Rscc_clustering.R")

cl_labels <- c(0, 0, 1, 2, 1, 2, 2, 0, 1, 1)
cl_labels_wNA <- c(0, 0, 1, 2, NA, 2, 2, NA, 1, 1)

cl_obj_simple <- structure(cl_labels,
                           cluster_count = 3L,
                           ids = NULL,
                           class = c("Rscc_clustering"))

cl_obj_wNA <- structure(cl_labels_wNA,
                        cluster_count = 3L,
                        ids = NULL,
                        class = c("Rscc_clustering"))

cl_obj_wID <- structure(cl_labels,
                        cluster_count = 3L,
                        ids = letters[1:10],
                        class = c("Rscc_clustering"))

cl_obj_wIDnNA <- structure(cl_labels_wNA,
                        cluster_count = 3L,
                        ids = letters[1:10],
                        class = c("Rscc_clustering"))

test_that("`Rscc_clustering` makes correct output.", {
  expect_equal(Rscc_clustering(c("A", "A", "B", "C", "B", "C", "C", "A", "B", "B")),
               cl_obj_simple)
  expect_equal(Rscc_clustering(c(1, 1, 2, 3, 2, 3, 3, 1, 2, 2)),
               cl_obj_simple)
  expect_equal(Rscc_clustering(c("A", "A", "B", "C", "NONE", "C", "C", "NONE", "B", "B"), "NONE"),
               cl_obj_wNA)
  expect_equal(Rscc_clustering(c("A", "A", "B", "C", "NONE", "C", "C", "0", "B", "B"), c("NONE", "0")),
               cl_obj_wNA)
  expect_equal(Rscc_clustering(c("A", "A", "B", "C", "B", "C", "C", "A", "B", "B"), ids = letters[1:10]),
               cl_obj_wID)
  expect_equal(Rscc_clustering(c("A", "A", "B", "C", "NONE", "C", "C", "0", "B", "B"), c("NONE", "0", "nnn"), ids = letters[1:10]),
               cl_obj_wIDnNA)
})

test_that("`Rscc_clustering` checks input.", {
  expect_error(Rscc_clustering(NULL))
  expect_error(Rscc_clustering(c("A", "A", "B", "C", "B", "C", "C", "A", "B", "B"), ids = dist(1:10)))
  expect_error(Rscc_clustering(c("A", "A", "B", "C", "B", "C", "C", "A", "B", "B"), ids = letters[1:9]))
})

test_that("`Rscc_clustering` converts to data.frame", {
  expect_equal(as.data.frame(cl_obj_simple),
               data.frame(id = 1:10, cluster_label = cl_labels))
  expect_equal(as.data.frame(cl_obj_wNA),
               data.frame(id = 1:10, cluster_label = cl_labels_wNA))
  expect_equal(as.data.frame(cl_obj_wID),
               data.frame(id = letters[1:10], cluster_label = cl_labels))
  expect_equal(as.data.frame(cl_obj_wIDnNA),
               data.frame(id = letters[1:10], cluster_label = cl_labels_wNA))
})
