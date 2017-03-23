/* =============================================================================
 * scclust -- A C library for size constrained clustering
 * https://github.com/fsavje/scclust
 *
 * Copyright (C) 2015-2017  Fredrik Savje -- http://fredriksavje.com
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library. If not, see http://www.gnu.org/licenses/
 * ========================================================================== */

#include "../include/scclust_spi.h"

#include <stddef.h>
#include "dist_search.h"
#include "dist_search_imp.h"


// =============================================================================
// External variable initialization
// =============================================================================

// See "dist_search.h" for definition
iscc_dist_functions_struct iscc_dist_functions = {
	.check_data_set = iscc_imp_check_data_set,
	.num_data_points = iscc_imp_num_data_points,
	.get_dist_matrix = iscc_imp_get_dist_matrix,
	.get_dist_rows = iscc_imp_get_dist_rows,
	.init_max_dist_object = iscc_imp_init_max_dist_object,
	.get_max_dist = iscc_imp_get_max_dist,
	.close_max_dist_object = iscc_imp_close_max_dist_object,
	.init_nn_search_object = iscc_imp_init_nn_search_object,
	.nearest_neighbor_search = iscc_imp_nearest_neighbor_search,
	.close_nn_search_object = iscc_imp_close_nn_search_object,
};


// =============================================================================
// Public function implementations
// =============================================================================

bool scc_reset_dist_functions(void)
{
	iscc_dist_functions = (iscc_dist_functions_struct) {
		.check_data_set = iscc_imp_check_data_set,
		.num_data_points = iscc_imp_num_data_points,
		.get_dist_matrix = iscc_imp_get_dist_matrix,
		.get_dist_rows = iscc_imp_get_dist_rows,
		.init_max_dist_object = iscc_imp_init_max_dist_object,
		.get_max_dist = iscc_imp_get_max_dist,
		.close_max_dist_object = iscc_imp_close_max_dist_object,
		.init_nn_search_object = iscc_imp_init_nn_search_object,
		.nearest_neighbor_search = iscc_imp_nearest_neighbor_search,
		.close_nn_search_object = iscc_imp_close_nn_search_object,
	};

	return true;
}


bool scc_set_dist_functions(scc_check_data_set check_data_set,
                            scc_num_data_points num_data_points,
                            scc_get_dist_matrix get_dist_matrix,
                            scc_get_dist_rows get_dist_rows,
                            scc_init_max_dist_object init_max_dist_object,
                            scc_get_max_dist get_max_dist,
                            scc_close_max_dist_object close_max_dist_object,
                            scc_init_nn_search_object init_nn_search_object,
                            scc_nearest_neighbor_search nearest_neighbor_search,
                            scc_close_nn_search_object close_nn_search_object)
{
	if (check_data_set != NULL) {
		iscc_dist_functions.check_data_set = check_data_set;
	}

	if (num_data_points != NULL) {
		iscc_dist_functions.num_data_points = num_data_points;
	}

	if (get_dist_matrix != NULL) {
		iscc_dist_functions.get_dist_matrix = get_dist_matrix;
	}

	if (get_dist_rows != NULL) {
		iscc_dist_functions.get_dist_rows = get_dist_rows;
	}

	if (init_max_dist_object != NULL &&
			get_max_dist != NULL &&
			close_max_dist_object != NULL) {
		iscc_dist_functions.init_max_dist_object = init_max_dist_object;
		iscc_dist_functions.get_max_dist = get_max_dist;
		iscc_dist_functions.close_max_dist_object = close_max_dist_object;
	} else if (init_max_dist_object != NULL ||
			get_max_dist != NULL ||
			close_max_dist_object != NULL) {
		return false;
	}

	if (init_nn_search_object != NULL &&
			nearest_neighbor_search != NULL &&
			close_nn_search_object != NULL) {
		iscc_dist_functions.init_nn_search_object = init_nn_search_object;
		iscc_dist_functions.nearest_neighbor_search = nearest_neighbor_search;
		iscc_dist_functions.close_nn_search_object = close_nn_search_object;
	} else if (init_nn_search_object != NULL ||
			nearest_neighbor_search != NULL ||
			close_nn_search_object != NULL) {
		return false;
	}

	return true;
}
