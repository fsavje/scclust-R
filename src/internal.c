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

#include "internal.h"
#include <stdbool.h>
#include <stddef.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <scclust_spi.h>
#include "error.h"


bool Rscc_dist_functions_are_set = false;

bool Rscc_check_data_set(void* data_set,
                         size_t num_data_points);


void Rscc_set_dist_functions__(void)
{
	if (!scc_set_dist_functions(Rscc_check_data_set,
                               (scc_get_dist_matrix) R_GetCCallable("distances", "idist_get_dist_matrix"),
                               (scc_get_dist_rows) R_GetCCallable("distances", "idist_get_dist_rows"),
                               (scc_init_max_dist_object) R_GetCCallable("distances", "idist_init_max_distance_search"),
                               (scc_get_max_dist) R_GetCCallable("distances", "idist_max_distance_search"),
                               (scc_close_max_dist_object) R_GetCCallable("distances", "idist_close_max_distance_search"),
                               (scc_init_nn_search_object) R_GetCCallable("distances", "idist_init_nearest_neighbor_search"),
                               (scc_nearest_neighbor_search) R_GetCCallable("distances", "idist_nearest_neighbor_search"),
                               (scc_close_nn_search_object) R_GetCCallable("distances", "idist_close_nearest_neighbor_search"))) {
		iRscc_scc_error();
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


int idist_num_data_points(const SEXP R_distances)
{
	static int(*func)(SEXP) = NULL;
	if (func == NULL) {
		func = (int(*)(SEXP)) R_GetCCallable("distances", "idist_num_data_points");
	}
	return func(R_distances);
}


bool Rscc_check_data_set(void* const data_set,
                         const size_t num_data_points)
{
	if (data_set == NULL) return false;
	SEXP R_distances = (SEXP) data_set;
	if (!idist_check_distance_object(R_distances)) return false;
	const size_t num_dp = (size_t) idist_num_data_points(R_distances);
	if (num_dp < num_data_points) return false;
	return true;
}
