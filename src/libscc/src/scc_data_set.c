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

void scc_free_data_set_object(scc_DataSetObject** const out_data_set_object)
{
	if ((out_data_set_object != NULL) && (*out_data_set_object != NULL)) {
		if (!((*out_data_set_object)->external_matrix)) free((void*) (*out_data_set_object)->data_matrix);
		free(*out_data_set_object);
		*out_data_set_object = NULL;
	}
}

scc_ErrorCode scc_get_data_set_object(const uintmax_t num_data_points,
                                      const uintmax_t num_dimensions,
                                      const size_t len_data_matrix,
                                      double data_matrix[const],
                                      const bool deep_matrix_copy,
                                      scc_DataSetObject** const out_data_set_object)
{
	if (out_data_set_object == NULL) return iscc_make_error(SCC_ER_NULL_INPUT);
	// Initialize to null, so subsequent functions detect invalid clustering
	// if user doesn't check for errors.
	*out_data_set_object = NULL;

	if (num_data_points == 0) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (num_data_points > ISCC_DPID_MAX) return iscc_make_error(SCC_ER_TOO_LARGE_PROBLEM);
	if (num_data_points > SIZE_MAX - 1) return iscc_make_error(SCC_ER_TOO_LARGE_PROBLEM);
	if (num_dimensions == 0) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (num_dimensions > UINT16_MAX) return iscc_make_error(SCC_ER_TOO_LARGE_PROBLEM);
	if (len_data_matrix < num_data_points * num_dimensions) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (data_matrix == NULL) return iscc_make_error(SCC_ER_NULL_INPUT);

	scc_DataSetObject* tmp_dso = malloc(sizeof(scc_DataSetObject));
	if (tmp_dso == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);

	*tmp_dso = (scc_DataSetObject) {
		.data_set_object_version = ISCC_CURRENT_DATASETOBJ_VERSION,
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

	*out_data_set_object = tmp_dso;

	return iscc_no_error();
}
