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

#include "../include/scclust.h"

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "error.h"
#include "data_set_struct.h"
#include "scclust_types.h"


// =============================================================================
// Public function implementations
// =============================================================================

scc_ErrorCode scc_init_data_set(const uint64_t num_data_points,
                                const uint32_t num_dimensions,
                                const size_t len_data_matrix,
                                const double data_matrix[const],
                                scc_DataSet** const out_data_set)
{
	if (out_data_set == NULL) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Output parameter may not be NULL.");
	}
	// Initialize to null, so subsequent functions detect invalid clustering
	// if user doesn't check for errors.
	*out_data_set = NULL;

	if (num_data_points == 0) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Data set must have positive number of data points.");
	}
	if (num_data_points > ISCC_POINTINDEX_MAX) {
		return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many data points (adjust the `scc_PointIndex` type).");
	}
	if (num_data_points > SIZE_MAX - 1) {
		return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many data points.");
	}
	if (num_dimensions == 0) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Data set must have positive number of dimensions.");
	}
	if (num_dimensions > UINT16_MAX) {
		return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many data dimensions.");
	}
	if (len_data_matrix < num_data_points * num_dimensions) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid data matrix.");
	}
	if (data_matrix == NULL) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid data matrix.");
	}

	scc_DataSet* tmp_dso = malloc(sizeof(scc_DataSet));
	if (tmp_dso == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);

	*tmp_dso = (scc_DataSet) {
		.data_set_version = ISCC_DATASET_STRUCT_VERSION,
		.num_data_points = (size_t) num_data_points,
		.num_dimensions = (uint_fast16_t) num_dimensions,
		.data_matrix = data_matrix,
	};

	*out_data_set = tmp_dso;

	return iscc_no_error();
}


void scc_free_data_set(scc_DataSet** const data_set)
{
	if ((data_set != NULL) && (*data_set != NULL)) {
		free(*data_set);
		*data_set = NULL;
	}
}


bool scc_is_initialized_data_set(const scc_DataSet* const data_set)
{
	if (data_set == NULL) return false;
	if (data_set->data_set_version != ISCC_DATASET_STRUCT_VERSION) return false;
	if (data_set->num_data_points == 0) return false;
	if (data_set->num_dimensions == 0) return false;
	if (data_set->data_matrix == NULL) return false;
	return true;
}
