/* =============================================================================
 * Rscclust -- R wrapper for the scclust library
 * https://github.com/fsavje/Rscclust
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

#ifndef RSCC_NNG_HG
#define RSCC_NNG_HG

#include <R.h>
#include <Rinternals.h>

SEXP Rscc_nng_clustering(SEXP R_distance_object,
                         SEXP R_size_constraint,
                         SEXP R_seed_method,
                         SEXP R_main_unassigned_method,
                         SEXP R_main_radius,
                         SEXP R_main_data_points,
                         SEXP R_secondary_unassigned_method,
                         SEXP R_secondary_radius);

SEXP Rscc_nng_clustering_batches(SEXP R_distance_object,
                                 SEXP R_size_constraint,
                                 SEXP R_main_unassigned_method,
                                 SEXP R_main_radius,
                                 SEXP R_main_data_points,
                                 SEXP R_batch_size);

SEXP Rscc_nng_clustering_types(SEXP R_distance_object,
                               SEXP R_type_labels,
                               SEXP R_type_size_constraints,
                               SEXP R_total_size_constraint,
                               SEXP R_seed_method,
                               SEXP R_main_unassigned_method,
                               SEXP R_main_radius,
                               SEXP R_main_data_points,
                               SEXP R_secondary_unassigned_method,
                               SEXP R_secondary_radius);

#endif // ifndef RSCC_NNG_HG
