/* =============================================================================
 * scclust for R -- R wrapper for the scclust library
 * https://github.com/fsavje/scclust-R
 *
 * Copyright (C) 2016  Fredrik Savje -- http://fredriksavje.com
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/
 * ========================================================================== */

#ifndef RSCC_MAKE_CLUSTERING_HG
#define RSCC_MAKE_CLUSTERING_HG

#include <R.h>
#include <Rinternals.h>

SEXP Rscc_make_clustering(SEXP R_distances,
                          SEXP R_size_constraint,
                          SEXP R_type_labels,
                          SEXP R_type_constraints,
                          SEXP R_seed_method,
                          SEXP R_primary_data_points,
                          SEXP R_primary_unassigned_method,
                          SEXP R_secondary_unassigned_method,
                          SEXP R_seed_radius,
                          SEXP R_primary_radius,
                          SEXP R_secondary_radius,
                          SEXP R_batch_size);

#endif // ifndef RSCC_MAKE_CLUSTERING_HG
