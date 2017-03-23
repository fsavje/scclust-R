/* =============================================================================
 * scclust for R -- R wrapper for the scclust library
 * https://github.com/fsavje/scclust-R
 *
 * Copyright (C) 2016-2017  Fredrik Savje -- http://fredriksavje.com
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

#include "internal.h"
#include <stdbool.h>
#include <stddef.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <scclust_spi.h>
#include "error.h"


// =============================================================================
// Extern variables
// =============================================================================

bool Rscc_dist_functions_are_set = false;


// =============================================================================
// External function implementations
// =============================================================================

void Rscc_set_dist_functions__(void)
{
	if (!scc_set_dist_functions((scc_check_data_set) R_GetCCallable("distances", "idist_check_distance_object"),
	                            (scc_num_data_points) idist_num_data_points,
	                            (scc_get_dist_matrix) R_GetCCallable("distances", "idist_get_dist_matrix"),
	                            (scc_get_dist_rows) R_GetCCallable("distances", "idist_get_dist_columns"),
	                            (scc_init_max_dist_object) R_GetCCallable("distances", "idist_init_max_distance_search"),
	                            (scc_get_max_dist) R_GetCCallable("distances", "idist_max_distance_search"),
	                            (scc_close_max_dist_object) R_GetCCallable("distances", "idist_close_max_distance_search"),
	                            (scc_init_nn_search_object) R_GetCCallable("distances", "idist_init_nearest_neighbor_search"),
	                            (scc_nearest_neighbor_search) R_GetCCallable("distances", "idist_nearest_neighbor_search"),
	                            (scc_close_nn_search_object) R_GetCCallable("distances", "idist_close_nearest_neighbor_search"))) {
		iRscc_error("Could not set distance search functions in scclust.");
	}
	Rscc_dist_functions_are_set = true;
}


bool idist_check_distance_object(const SEXP R_distances)
{
	static bool(*func)(SEXP) = NULL;
	if (func == NULL) {
		func = (bool(*)(SEXP)) R_GetCCallable("distances", "idist_check_distance_object");
	}
	return func(R_distances);
}


size_t idist_num_data_points(const SEXP R_distances)
{
	static int(*func)(SEXP) = NULL;
	if (func == NULL) {
		func = (int(*)(SEXP)) R_GetCCallable("distances", "idist_num_data_points");
	}
	return (size_t) func(R_distances);
}
