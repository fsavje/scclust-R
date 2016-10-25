/* ==============================================================================
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
 * ============================================================================== */

#include "dist_search.h"

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include "../include/scc_data_obj.h"
#include "dist_inline_sqdist.h"
#include "scclust_int.h"


// ==============================================================================
// Internal function prototypes
// ==============================================================================

static inline void iscc_add_dist_to_list(double add_dist,
                                         iscc_Dpid add_index,
                                         double* dist_list,
                                         iscc_Dpid* index_list,
                                         const double* dist_list_start);


// ==============================================================================
// External function implementations
// ==============================================================================

struct iscc_NNSearchObject {
	scc_DataSetObject* data_set_object;
	size_t len_search_indices;
	const iscc_Dpid* search_indices;
};


bool iscc_init_nn_search_object(void* const data_set_object,
                                const size_t len_search_indices,
                                const iscc_Dpid search_indices[const],
                                iscc_NNSearchObject** const out_nn_search_object)
{
	assert(iscc_check_data_set_object(data_set_object, 1));
	assert(len_search_indices > 0);
	assert(out_nn_search_object != NULL);

	*out_nn_search_object = malloc(sizeof(iscc_NNSearchObject));
	if (*out_nn_search_object == NULL) return false;

	**out_nn_search_object = (iscc_NNSearchObject) {
		.data_set_object = data_set_object,
		.len_search_indices = len_search_indices,
		.search_indices = search_indices,
	};

	return true;
}


bool iscc_nearest_neighbor_search_digraph(iscc_NNSearchObject* const nn_search_object,
                                          const size_t len_query_indicators,
                                          const bool query_indicators[const],
                                          bool out_query_indicators[const],
                                          const uint32_t k,
                                          const bool radius_search,
                                          const double radius,
                                          const bool accept_partial,
                                          iscc_Arci out_nn_ref[const],
                                          iscc_Dpid out_nn_indices[const])
{
	assert(nn_search_object != NULL);
	scc_DataSetObject* const data_set_object = nn_search_object->data_set_object;
	const size_t len_search_indices = nn_search_object->len_search_indices;
	const iscc_Dpid* const search_indices = nn_search_object->search_indices;

	assert(iscc_check_data_set_object(data_set_object, len_query_indicators));
	assert(len_search_indices > 0);
	assert(len_query_indicators > 0);
	assert(k > 0);
	assert(k <= len_search_indices);
	assert(!radius_search || (radius > 0.0));
	assert(out_nn_ref != NULL);
	assert(out_nn_indices != NULL);

	double tmp_dist;
	iscc_Dpid* index_write = out_nn_indices;
	double* const sort_scratch = malloc(sizeof(double[k]));
	if (sort_scratch == NULL) return false;
	double* const sort_scratch_end = sort_scratch + k - 1;

	out_nn_ref[0] = 0;
	if (search_indices == NULL) {
		for (size_t q = 0; q < len_query_indicators; ++q) {
			if ((query_indicators == NULL) || query_indicators[q]) {
				size_t s = 0;
				uint32_t found;
				iscc_Dpid* const index_write_end = index_write + k - 1;

				if (radius_search) {
					const double radius_sq = radius * radius;
					found = 0;
					for (; (s < len_search_indices) && (found < k); ++s) {
						tmp_dist = iscc_get_sq_dist(data_set_object, q, s);
						if (tmp_dist > radius_sq) continue;
						iscc_add_dist_to_list(tmp_dist, (iscc_Dpid) s, sort_scratch + found, index_write + found, sort_scratch);
						++found;
					}
				} else {
					found = k;
					for (; s < k; ++s) {
						tmp_dist = iscc_get_sq_dist(data_set_object, q, s);
						iscc_add_dist_to_list(tmp_dist, (iscc_Dpid) s, sort_scratch + s, index_write + s, sort_scratch);
					}
				}

				for (; s < len_search_indices; ++s) {
					assert(found == k);
					tmp_dist = iscc_get_sq_dist(data_set_object, q, s);
					if (tmp_dist >= *sort_scratch_end) continue;
					iscc_add_dist_to_list(tmp_dist, (iscc_Dpid) s, sort_scratch_end, index_write_end, sort_scratch);
				}

				if (radius_search && !accept_partial && (found < k)) {
					found = 0;
					if (out_query_indicators != NULL) out_query_indicators[q] = false;
				}
				out_nn_ref[q + 1] = out_nn_ref[q] + found;
				index_write += found;
			} else {
				out_nn_ref[q + 1] = out_nn_ref[q];
			}
		}

	} else if (search_indices != NULL) {
		for (size_t q = 0; q < len_query_indicators; ++q) {
			if ((query_indicators == NULL) || query_indicators[q]) {
				size_t s = 0;
				uint32_t found;
				iscc_Dpid* const index_write_end = index_write + k - 1;

				if (radius_search) {
					const double radius_sq = radius * radius;
					found = 0;
					for (; (s < len_search_indices) && (found < k); ++s) {
						tmp_dist = iscc_get_sq_dist(data_set_object, q, (size_t) search_indices[s]);
						if (tmp_dist > radius_sq) continue;
						iscc_add_dist_to_list(tmp_dist, search_indices[s], sort_scratch + found, index_write + found, sort_scratch);
						++found;
					}
				} else {
					found = k;
					for (; s < k; ++s) {
						tmp_dist = iscc_get_sq_dist(data_set_object, q, (size_t) search_indices[s]);
						iscc_add_dist_to_list(tmp_dist, search_indices[s], sort_scratch + s, index_write + s, sort_scratch);
					}
				}

				for (; s < len_search_indices; ++s) {
					assert(found == k);
					tmp_dist = iscc_get_sq_dist(data_set_object, q, (size_t) search_indices[s]);
					if (tmp_dist >= *sort_scratch_end) continue;
					iscc_add_dist_to_list(tmp_dist, search_indices[s], sort_scratch_end, index_write_end, sort_scratch);
				}

				if (radius_search && !accept_partial && (found < k)) {
					found = 0;
					if (out_query_indicators != NULL) out_query_indicators[q] = false;
				}
				out_nn_ref[q + 1] = out_nn_ref[q] + found;
				index_write += found;
			} else {
				out_nn_ref[q + 1] = out_nn_ref[q];
			}
		}
	}

	free(sort_scratch);

	return true;
}


bool iscc_nearest_neighbor_search_index(iscc_NNSearchObject* const nn_search_object,
                                        const size_t len_query_indices,
                                        const iscc_Dpid query_indices[const],
                                        const uint32_t k,
                                        const bool radius_search,
                                        const double radius,
                                        iscc_Dpid out_nn_indices[const])
{
	assert(nn_search_object != NULL);
	scc_DataSetObject* const data_set_object = nn_search_object->data_set_object;
	const size_t len_search_indices = nn_search_object->len_search_indices;
	const iscc_Dpid* const search_indices = nn_search_object->search_indices;

	assert(iscc_check_data_set_object(data_set_object, 1));
	assert(len_search_indices > 0);
	assert(len_query_indices > 0);
	assert(query_indices != NULL);
	assert(k > 0);
	assert(k <= len_search_indices);
	assert(!radius_search || (radius > 0.0));
	assert(out_nn_indices != NULL);

	double tmp_dist;
	iscc_Dpid* index_write = out_nn_indices;
	double* const sort_scratch = malloc(sizeof(double[k]));
	if (sort_scratch == NULL) return false;
	double* const sort_scratch_end = sort_scratch + k - 1;

	if (search_indices == NULL) {
		for (size_t q = 0; q < len_query_indices; ++q) {
			size_t s = 0;
			iscc_Dpid* const index_write_end = index_write + k - 1;

			if (radius_search) {
				const double radius_sq = radius * radius;
				uint32_t found = 0;
				for (; (s < len_search_indices) && (found < k); ++s) {
					tmp_dist = iscc_get_sq_dist(data_set_object, (size_t) query_indices[q], s);
					if (tmp_dist > radius_sq) continue;
					iscc_add_dist_to_list(tmp_dist, (iscc_Dpid) s, sort_scratch + found, index_write + found, sort_scratch);
					++found;
				}
				for (; found < k; ++found) {
					index_write[found] = ISCC_DPID_NA;
				}
			} else {
				for (; s < k; ++s) {
					tmp_dist = iscc_get_sq_dist(data_set_object, (size_t) query_indices[q], s);
					iscc_add_dist_to_list(tmp_dist, (iscc_Dpid) s, sort_scratch + s, index_write + s, sort_scratch);
				}
			}

			for (; s < len_search_indices; ++s) {
				tmp_dist = iscc_get_sq_dist(data_set_object, (size_t) query_indices[q], s);
				if (tmp_dist >= *sort_scratch_end) continue;
				iscc_add_dist_to_list(tmp_dist, (iscc_Dpid) s, sort_scratch_end, index_write_end, sort_scratch);
			}
			index_write += k;
		}
	} else if (search_indices != NULL) {
		for (size_t q = 0; q < len_query_indices; ++q) {
			size_t s = 0;
			iscc_Dpid* const index_write_end = index_write + k - 1;

			if (radius_search) {
				const double radius_sq = radius * radius;
				uint32_t found = 0;
				for (; (s < len_search_indices) && (found < k); ++s) {
					tmp_dist = iscc_get_sq_dist(data_set_object, (size_t) query_indices[q], (size_t) search_indices[s]);
					if (tmp_dist > radius_sq) continue;
					iscc_add_dist_to_list(tmp_dist, search_indices[s], sort_scratch + found, index_write + found, sort_scratch);
					++found;
				}
				for (; found < k; ++found) {
					index_write[found] = ISCC_DPID_NA;
				}
			} else {
				for (; s < k; ++s) {
					tmp_dist = iscc_get_sq_dist(data_set_object, (size_t) query_indices[q], (size_t) search_indices[s]);
					iscc_add_dist_to_list(tmp_dist, search_indices[s], sort_scratch + s, index_write + s, sort_scratch);
				}
			}

			for (; s < len_search_indices; ++s) {
				tmp_dist = iscc_get_sq_dist(data_set_object, (size_t) query_indices[q], (size_t) search_indices[s]);
				if (tmp_dist >= *sort_scratch_end) continue;
				iscc_add_dist_to_list(tmp_dist, search_indices[s], sort_scratch_end, index_write_end, sort_scratch);
			}
			index_write += k;
		}
	}

	free(sort_scratch);

	return true;
}


bool iscc_close_nn_search_object(iscc_NNSearchObject** const nn_search_object)
{
	if (nn_search_object != NULL) {
		free(*nn_search_object);
		*nn_search_object = NULL;
	}
	return true;
}

// ==============================================================================
// Internal function implementations
// ==============================================================================

static inline void iscc_add_dist_to_list(const double add_dist,
                                         const iscc_Dpid add_index,
                                         double* dist_list,
                                         iscc_Dpid* index_list,
                                         const double* const dist_list_start)
{
	assert(dist_list != NULL);
	assert(index_list != NULL);
	assert(dist_list_start != NULL);

	for (; (dist_list != dist_list_start) && (add_dist < dist_list[-1]); --dist_list, --index_list) {
		dist_list[0] = dist_list[-1];
		index_list[0] = index_list[-1];
	}
	dist_list[0] = add_dist;
	index_list[0] = add_index;
}
