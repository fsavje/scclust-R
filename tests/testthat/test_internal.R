# ==============================================================================
# scclust for R -- R wrapper for the scclust library
# https://github.com/fsavje/scclust-R
#
# Copyright (C) 2016-2017  Fredrik Savje -- http://fredriksavje.com
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

library(scclust)
context("Internal functions")


# ==============================================================================
# Constants
# ==============================================================================

test_that("Constants are correct", {
  expect_identical(all_seed_methods,
                   c("lexical",
                     "batches",
                     "inwards_order",
                     "inwards_updating",
                     "exclusion_order",
                     "exclusion_updating"))
})


# ==============================================================================
# make_scclust
# ==============================================================================

test_that("`make_scclust` constructs correct object", {
  expect_identical(make_scclust(1:10, 10L, letters[1:10]),
                   structure(1:10,
                             cluster_count = 10L,
                             ids = letters[1:10],
                             class = c("scclust")))
})


# ==============================================================================
# make_type_indicators
# ==============================================================================

test_that("`make_type_indicators` checks input.", {
  expect_silent(make_type_indicators(c(1L, 3L), rep(1:3, 100L)))
  expect_silent(make_type_indicators(c("1", "3"), factor(rep(1:3, 100L))))
  expect_error(make_type_indicators(c(1L, 1L, 3L), rep(1:3, 100L)))
  expect_error(make_type_indicators(c(1L, 3L), letters))
  expect_error(make_type_indicators(c(1L, 3L, 4L), rep(1:3, 100L)))
  expect_error(make_type_indicators(c("1", "3", "a"), factor(rep(1:3, 100L))))
})

test_that("`make_type_indicators` returns correct output.", {
  expect_identical(make_type_indicators(c(1L, 3L), rep(1:3, 100L)),
                   c("0" = FALSE, "1" = TRUE, "2" = FALSE, "3" = TRUE))
  expect_identical(make_type_indicators(c("1", "3"), factor(rep(1:3, 100L))),
                   c("0" = FALSE, "1" = TRUE, "2" = FALSE, "3" = TRUE))
})


# ==============================================================================
# make_type_size_constraints
# ==============================================================================

test_that("`make_type_size_constraints` checks input.", {
  expect_silent(make_type_size_constraints(c("1" = 3L, "3" = 2L), rep(1:3, 100L)))
  expect_silent(make_type_size_constraints(c("1" = 3L, "3" = 2L), factor(rep(1:3, 100L))))
  expect_error(make_type_size_constraints(c("1" = "a", "3" = "b"), rep(1:3, 100L)))
  expect_error(make_type_size_constraints(c("1" = 3L, "3" = NA), rep(1:3, 100L)))
  expect_error(make_type_size_constraints(c("1" = -3L, "3" = 2L), rep(1:3, 100L)))
  expect_error(make_type_size_constraints(c(3L, 2L), rep(1:3, 100L)))
  expect_error(make_type_size_constraints(c("1" = 3L, "3" = 2L, "1" = 2L), rep(1:3, 100L)))
  expect_error(make_type_size_constraints(c("1" = 3L, "3" = 2L), letters))
  expect_error(make_type_size_constraints(c("1" = 3L, "3" = 2L, "4" = 1L), rep(1:3, 100L)))
  expect_error(make_type_size_constraints(c("1" = 3L, "3" = 2L, "4" = 1L), factor(rep(1:3, 100L))))
})

test_that("`make_type_size_constraints` returns correct output.", {
  expect_identical(make_type_size_constraints(c("1" = 3L, "3" = 2L), rep(1:3, 100L)),
                   c("0" = 0L, "1" = 3L, "2" = 0L, "3" = 2L))
  expect_identical(make_type_size_constraints(c("1" = 3L, "3" = 2L), factor(rep(1:3, 100L))),
                   c("0" = 0L, "1" = 3L, "2" = 0L, "3" = 2L))
})
