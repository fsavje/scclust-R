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

#include "dist_search.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include "../include/scclust.h"
#include "dist_inline_sqdist.h"
#include "scc_data_set_struct.h"
#include "scclust_internal.h"


// =============================================================================
// Miscellaneous function implementations
// =============================================================================

bool iscc_get_dist_matrix(void* const data_set,
                          const size_t len_point_indices,
                          const iscc_Dpid point_indices[const],
                          double output_dists[])
{
	if (!scc_is_initialized_data_set(data_set, 1)) return false;
	assert(len_point_indices > 1);
	assert(output_dists != NULL);

	if (point_indices == NULL) {
		for (size_t p1 = 0; p1 < len_point_indices; ++p1) {
			for (size_t p2 = p1 + 1; p2 < len_point_indices; ++p2) {
				*output_dists = sqrt(iscc_get_sq_dist(data_set, p1, p2));
				++output_dists;
			}
		}
	} else {
		for (size_t p1 = 0; p1 < len_point_indices; ++p1) {
			for (size_t p2 = p1 + 1; p2 < len_point_indices; ++p2) {
				*output_dists = sqrt(iscc_get_sq_dist(data_set, (size_t) point_indices[p1], (size_t) point_indices[p2]));
				++output_dists;
			}
		}
	}

	return true;
}


bool iscc_get_dist_rows(void* const data_set,
                        const size_t len_query_indices,
                        const iscc_Dpid query_indices[const],
                        const size_t len_column_indices,
                        const iscc_Dpid column_indices[const],
                        double output_dists[])
{
	if (!scc_is_initialized_data_set(data_set, 1)) return false;
	assert(len_query_indices > 0);
	assert(len_column_indices > 0);
	assert(output_dists != NULL);

	if ((query_indices != NULL) && (column_indices != NULL)) {
		for (size_t q = 0; q < len_query_indices; ++q) {
			for (size_t c = 0; c < len_column_indices; ++c) {
				*output_dists = sqrt(iscc_get_sq_dist(data_set, (size_t) query_indices[q], (size_t) column_indices[c]));
				++output_dists;
			}
		}

	} else if ((query_indices == NULL) && (column_indices != NULL)) {
		for (size_t q = 0; q < len_query_indices; ++q) {
			for (size_t c = 0; c < len_column_indices; ++c) {
				*output_dists = sqrt(iscc_get_sq_dist(data_set, q, (size_t) column_indices[c]));
				++output_dists;
			}
		}

	} else if ((query_indices != NULL) && (column_indices == NULL)) {
		for (size_t q = 0; q < len_query_indices; ++q) {
			for (size_t c = 0; c < len_column_indices; ++c) {
				*output_dists = sqrt(iscc_get_sq_dist(data_set, (size_t) query_indices[q], c));
				++output_dists;
			}
		}

	} else if ((query_indices == NULL) && (column_indices == NULL)) {
		for (size_t q = 0; q < len_query_indices; ++q) {
			for (size_t c = 0; c < len_column_indices; ++c) {
				*output_dists = sqrt(iscc_get_sq_dist(data_set, q, c));
				++output_dists;
			}
		}
	}

	return true;
}


// =============================================================================
// Max dist search function implementations
// =============================================================================

struct iscc_MaxDistObject {
	scc_DataSet* data_set;
	size_t len_search_indices;
	const iscc_Dpid* search_indices;
};


bool iscc_init_max_dist_object(void* const data_set,
                               const size_t len_search_indices,
                               const iscc_Dpid search_indices[const],
                               iscc_MaxDistObject** const out_max_dist_object)
{
	if (!scc_is_initialized_data_set(data_set, 1)) return false;
	assert(len_search_indices > 0);
	assert(out_max_dist_object != NULL);

	*out_max_dist_object = malloc(sizeof(iscc_MaxDistObject));
	if (*out_max_dist_object == NULL) return false;

	**out_max_dist_object = (iscc_MaxDistObject) {
		.data_set = data_set,
		.len_search_indices = len_search_indices,
		.search_indices = search_indices,
	};

	return true;
}


bool iscc_get_max_dist(iscc_MaxDistObject* const max_dist_object,
                       const size_t len_query_indices,
                       const iscc_Dpid query_indices[const],
                       iscc_Dpid out_max_indices[const],
                       double out_max_dists[const])
{
	assert(max_dist_object != NULL);
	scc_DataSet* const data_set = max_dist_object->data_set;
	const size_t len_search_indices = max_dist_object->len_search_indices;
	const iscc_Dpid* const search_indices = max_dist_object->search_indices;

	if (!scc_is_initialized_data_set(data_set, 1)) return false;
	assert(len_search_indices > 0);
	assert(len_query_indices > 0);
	assert(out_max_indices != NULL);
	assert(out_max_dists != NULL);

	double tmp_dist;
	double max_dist;

	if ((query_indices != NULL) && (search_indices != NULL)) {
		for (size_t q = 0; q < len_query_indices; ++q) {
			max_dist = -1.0;
			for (size_t s = 0; s < len_search_indices; ++s) {
				tmp_dist = iscc_get_sq_dist(data_set, (size_t) query_indices[q], (size_t) search_indices[s]);
				if (max_dist < tmp_dist) {
					max_dist = tmp_dist;
					out_max_indices[q] = search_indices[s];
				}
			}
			out_max_dists[q] = sqrt(max_dist);
		}

	} else if ((query_indices == NULL) && (search_indices != NULL)) {
		for (size_t q = 0; q < len_query_indices; ++q) {
			max_dist = -1.0;
			for (size_t s = 0; s < len_search_indices; ++s) {
				tmp_dist = iscc_get_sq_dist(data_set, q, (size_t) search_indices[s]);
				if (max_dist < tmp_dist) {
					max_dist = tmp_dist;
					out_max_indices[q] = search_indices[s];
				}
			}
			out_max_dists[q] = sqrt(max_dist);
		}

	} else if ((query_indices != NULL) && (search_indices == NULL)) {
		for (size_t q = 0; q < len_query_indices; ++q) {
			max_dist = -1.0;
			for (size_t s = 0; s < len_search_indices; ++s) {
				tmp_dist = iscc_get_sq_dist(data_set, (size_t) query_indices[q], s);
				if (max_dist < tmp_dist) {
					max_dist = tmp_dist;
					out_max_indices[q] = (iscc_Dpid) s;
				}
			}
			out_max_dists[q] = sqrt(max_dist);
		}

	} else if ((query_indices == NULL) && (search_indices == NULL)) {
		for (size_t q = 0; q < len_query_indices; ++q) {
			max_dist = -1.0;
			for (size_t s = 0; s < len_search_indices; ++s) {
				tmp_dist = iscc_get_sq_dist(data_set, q, s);
				if (max_dist < tmp_dist) {
					max_dist = tmp_dist;
					out_max_indices[q] = (iscc_Dpid) s;
				}
			}
			out_max_dists[q] = sqrt(max_dist);
		}
	}

	return true;
}


bool iscc_close_max_dist_object(iscc_MaxDistObject** const max_dist_object)
{
	if (max_dist_object != NULL) {
		free(*max_dist_object);
		*max_dist_object = NULL;
	}
	return true;
}
