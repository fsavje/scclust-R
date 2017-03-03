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

#ifndef RSCC_INTERNAL_HG
#define RSCC_INTERNAL_HG

#include <stdbool.h>
#include <R.h>
#include <Rinternals.h>

#define Rscc_set_dist_functions() (void)((Rscc_dist_functions_are_set) || (Rscc_set_dist_functions__(), 0))

#define Rscc_get_distances_pointer(distances) ((void*) distances)

extern bool Rscc_dist_functions_are_set;

void Rscc_set_dist_functions__(void);

bool idist_check_distance_object(SEXP R_distances);

int idist_num_data_points(SEXP R_distances);

#endif // ifndef RSCC_INTERNAL_HG
