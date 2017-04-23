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


# ==============================================================================
# Constants
# ==============================================================================

all_seed_methods <-  c("lexical",
                       "batches",
                       "inwards_order",
                       "inwards_updating",
                       "exclusion_order",
                       "exclusion_updating")


# ==============================================================================
# Class constructors
# ==============================================================================

make_scclust <- function(cluster_labels,
                         cluster_count,
                         ids) {
  structure(cluster_labels,
            cluster_count = cluster_count,
            ids = ids,
            class = c("scclust"))
}


# ==============================================================================
# Translation functions (from R to C)
# ==============================================================================

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
    out_type_size_constraints <- c("0" = 0L, out_type_size_constraints)
  } else if (is.integer(type_labels)) {
    max_label <- max(type_labels)
    out_type_size_constraints <- rep(0L, max_label + 1L)
    names(out_type_size_constraints) <- as.character(0L:max_label)
    stopifnot(all(names(type_constraints) %in% names(out_type_size_constraints)))
    out_type_size_constraints[names(type_constraints)] <- as.integer(type_constraints)
  }
  out_type_size_constraints
}
