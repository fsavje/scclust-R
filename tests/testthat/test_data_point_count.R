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
context("data_point_count")


test_that("`data_point_count` dispatches correctly", {
  expect_identical(data_point_count(make_distances(matrix(c(0.1, 0.2, 0.3, 0.1, 0.2, 0.3), ncol = 2))),
                   3L)
  expect_identical(data_point_count(make_distances(matrix(c(0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.3), ncol = 2))),
                   5L)
  expect_identical(data_point_count(Rscc_clustering(c("a", "b", "c", "a", "b", "c", "a"))),
                   7L)
  expect_identical(data_point_count(Rscc_clustering(c("a", "b", "c", "a", "b", "c", "a", "d", "e", "d"))),
                   10L)
  expect_error(data_point_count("invalid"),
               regexp = "Unknown class")
})



