/* =============================================================================
 * scclust -- A C library for size constrained clustering
 * https://github.com/fsavje/scclust
 *
 * Copyright (C) 2015-2016  Fredrik Savje -- http://fredriksavje.com
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
#include "scc_data_set_struct.h"
#include "scclust_internal.h"


// =============================================================================
// External function implementations
// =============================================================================

scc_ErrorCode scc_init_data_set(const uintmax_t num_data_points,
                                const uintmax_t num_dimensions,
                                const size_t len_data_matrix,
                                double data_matrix[const],
                                const bool deep_matrix_copy,
                                scc_DataSet** const out_data_set)
{
	if (out_data_set == NULL) return iscc_make_error(SCC_ER_NULL_INPUT);
	// Initialize to null, so subsequent functions detect invalid clustering
	// if user doesn't check for errors.
	*out_data_set = NULL;

	if (num_data_points == 0) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (num_data_points > ISCC_DPID_MAX) return iscc_make_error(SCC_ER_TOO_LARGE_PROBLEM);
	if (num_data_points > SIZE_MAX - 1) return iscc_make_error(SCC_ER_TOO_LARGE_PROBLEM);
	if (num_dimensions == 0) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (num_dimensions > UINT16_MAX) return iscc_make_error(SCC_ER_TOO_LARGE_PROBLEM);
	if (len_data_matrix < num_data_points * num_dimensions) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (data_matrix == NULL) return iscc_make_error(SCC_ER_NULL_INPUT);

	scc_DataSet* tmp_dso = malloc(sizeof(scc_DataSet));
	if (tmp_dso == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);

	*tmp_dso = (scc_DataSet) {
		.data_set_version = ISCC_DATASET_STRUCT_VERSION,
		.num_data_points = (size_t) num_data_points,
		.num_dimensions = (uint_fast16_t) num_dimensions,
		.data_matrix = NULL,
		.external_matrix = !deep_matrix_copy,
	};

	if (deep_matrix_copy) {
		tmp_dso->data_matrix = malloc(sizeof(double[tmp_dso->num_data_points * tmp_dso->num_dimensions]));
		if (tmp_dso->data_matrix == NULL) {
			free(tmp_dso);
			return iscc_make_error(SCC_ER_NO_MEMORY);
		}
		memcpy(tmp_dso->data_matrix, data_matrix, tmp_dso->num_data_points * tmp_dso->num_dimensions * sizeof(double));
	} else {
		tmp_dso->data_matrix = data_matrix;
	}

	assert(tmp_dso->data_matrix != NULL);

	*out_data_set = tmp_dso;

	return iscc_no_error();
}


void scc_free_data_set(scc_DataSet** const data_set)
{
	if ((data_set != NULL) && (*data_set != NULL)) {
		if (!((*data_set)->external_matrix)) free((void*) (*data_set)->data_matrix);
		free(*data_set);
		*data_set = NULL;
	}
}


bool scc_is_initialized_data_set(const scc_DataSet* const data_set,
                                 const uintmax_t num_data_points)
{
	if (data_set == NULL) return false;
	if (data_set->data_set_version != ISCC_DATASET_STRUCT_VERSION) return false;
	if (data_set->num_data_points < num_data_points) return false;
	if (data_set->num_dimensions == 0) return false;
	if (data_set->data_matrix == NULL) return false;
	return true;
}
