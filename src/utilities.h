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

#ifndef RSCC_UTILITIES_HG
#define RSCC_UTILITIES_HG

#include <R.h>
#include <Rinternals.h>


#define set_dist_functions)



bool dist_functions_are_set = false;




SEXP Rscc_set_dist_functions(SEXP R_dist_functions);

SEXP Rscc_check_clustering(SEXP R_clustering,
                           SEXP R_size_constraint,
                           SEXP R_type_labels,
                           SEXP R_type_constraints);

SEXP Rscc_get_clustering_stats(SEXP R_clustering,
                               SEXP R_distance_object);

#endif // ifndef RSCC_UTILITIES_HG
