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


# ==============================================================================
# Constants
# ==============================================================================

all_seed_methods <-  c("lexical",
                       "inwards_order",
                       "inwards_updating",
                       "inwards_alt_updating",
                       "exclusion_order",
                       "exclusion_updating")


# ==============================================================================
# Class constructors
# ==============================================================================

make_Rscc_clustering <- function(cluster_labels,
                                 cluster_count,
                                 ids) {
  structure(cluster_labels,
            cluster_count = cluster_count,
            ids = ids,
            class = c("Rscc_clustering"))
}


# ==============================================================================
# Property functions
# ==============================================================================

data_point_count_distances <- function(distance_obj) {
  stopifnot(is.matrix(distance_obj))
  ncol(distance_obj)
}


data_point_count_clustering <- function(clustering_obj) {
  length(clustering_obj)
}


cluster_count <- function(clustering) {
  stopifnot(is.integer(attr(clustering, "cluster_count", exact = TRUE)))
  attr(clustering, "cluster_count", exact = TRUE)
}


get_all_types <- function(type_labels) {
  stopifnot(is.factor(type_labels) || is.integer(type_labels))
  if (is.factor(type_labels)) {
    levels(type_labels)
  } else if (is.integer(type_labels)) {
    sort(unique(type_labels))
  }
}


# ==============================================================================
# Translation functions (from R to C)
# ==============================================================================

make_type_indicators <- function(targets,
                                 type_labels) {
  stopifnot(is.character(names(targets)),
            !anyDuplicated(names(targets)),
            is.factor(type_labels) || is.integer(type_labels))
  if (is.factor(type_labels)) {
    out_indicators <- rep(FALSE, nlevels(type_labels))
    names(out_indicators) <- levels(type_labels)
    stopifnot(all(as.character(targets) %in% names(out_indicators)))
    out_indicators[as.character(targets)] <- TRUE
    out_indicators <- c(FALSE, out_indicators)
  } else if (is.integer(type_labels)) {
    max_label <- max(type_labels)
    out_indicators <- rep(FALSE, max_label + 1L)
    names(out_indicators) <- as.character(0L:max_label)
    stopifnot(all(as.character(targets) %in% names(out_indicators)))
    out_indicators[as.character(targets)] <- TRUE
  }
  out_indicators
}


make_type_size_constraints <- function(type_constraints,
                                       type_labels) {
  stopifnot(is.integer(type_constraints),
            !any(is.na(type_constraints)),
            all(type_constraints >= 0L),
            is.character(names(type_constraints)),
            !anyDuplicated(names(type_constraints)),
            is.factor(type_labels) || is.integer(type_labels))
  if (is.factor(type_labels)) {
    out_type_size_constraints <- rep(0L, nlevels(type_labels))
    names(out_type_size_constraints) <- levels(type_labels)
    stopifnot(all(names(type_constraints) %in% names(out_type_size_constraints)))
    out_type_size_constraints[names(type_constraints)] <- as.integer(type_constraints)
    out_type_size_constraints <- c(0L, out_type_size_constraints)
  } else if (is.integer(type_labels)) {
    max_label <- max(type_labels)
    out_type_size_constraints <- rep(0L, max_label + 1L)
    names(out_type_size_constraints) <- as.character(0L:max_label)
    stopifnot(all(names(type_constraints) %in% names(out_type_size_constraints)))
    out_type_size_constraints[names(type_constraints)] <- as.integer(type_constraints)
  }
  unname(out_type_size_constraints)
}
