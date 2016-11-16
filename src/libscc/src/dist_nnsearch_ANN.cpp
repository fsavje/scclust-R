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

// So g++ defines integer limits
#define __STDC_LIMIT_MACROS

#include "dist_search.h"

#include <cassert>
#include <climits>
#include <cstddef>
#include "../exlib/libANN/include/ANN/ANN.h"
#include "../include/scclust.h"
#include "scc_data_set_struct.h"
#include "scclust_internal.h"

#ifdef SCC_ANN_BDTREE
	#define ANNpointSetConstructor ANNbd_tree
#else
	#define ANNpointSetConstructor ANNkd_tree
#endif

#ifndef SCC_ANN_EPS
	#define SCC_ANN_EPS 0.0
#endif

// =============================================================================
// Internal structs and variables
// =============================================================================

static int iscc_open_search_objects = 0;

struct iscc_NNSearchObject {
	scc_DataSetObject* data_set_object;
	size_t len_search_indices;
	const iscc_Dpid* search_indices;
	ANNpoint* search_points;
	ANNpointSet* search_tree;
};


// =============================================================================
// External function implementations
// =============================================================================

bool iscc_init_nn_search_object(void* const data_set_object,
                                const size_t len_search_indices,
                                const iscc_Dpid* const search_indices,
                                iscc_NNSearchObject** const out_nn_search_object)
{
	assert(iscc_open_search_objects >= 0);
	assert(iscc_check_data_set_object(data_set_object, 1));
	assert(len_search_indices > 0);
	assert(out_nn_search_object != NULL);

	if (len_search_indices > INT_MAX) return false;

	scc_DataSetObject* const data_set_object_cast = static_cast<scc_DataSetObject*>(data_set_object);

	ANNpoint* search_points;
	try {
		search_points = new ANNpoint[len_search_indices];
		*out_nn_search_object = new iscc_NNSearchObject;
	} catch (...) {
		return false;
	}

	if (search_indices == NULL) {
		assert(len_search_indices <= data_set_object_cast->num_data_points);
		double* search_point = data_set_object_cast->data_matrix;
		for (size_t i = 0; i < len_search_indices; ++i, search_point += data_set_object_cast->num_dimensions) {
			search_points[i] = search_point;
		}
	} else if (search_indices != NULL) {
		for (size_t i = 0; i < len_search_indices; ++i) {
			assert(static_cast<size_t>(search_indices[i]) < data_set_object_cast->num_data_points);
			search_points[i] = data_set_object_cast->data_matrix + search_indices[i] * data_set_object_cast->num_dimensions;
		}
	}

	ANNpointSet* search_tree;
	try {
		search_tree = new ANNpointSetConstructor(search_points,
		                                         static_cast<int>(len_search_indices),
		                                         static_cast<int>(data_set_object_cast->num_dimensions));
	} catch (...) {
		delete[] search_points;
		delete *out_nn_search_object;
		*out_nn_search_object = NULL;
		return false;
	}

	(*out_nn_search_object)->data_set_object = data_set_object_cast;
	(*out_nn_search_object)->len_search_indices = len_search_indices;
	(*out_nn_search_object)->search_indices = search_indices;
	(*out_nn_search_object)->search_points = search_points;
	(*out_nn_search_object)->search_tree = search_tree;

	++iscc_open_search_objects;
	return true;
}


bool iscc_nearest_neighbor_search_digraph(iscc_NNSearchObject* const nn_search_object,
                                          const size_t len_query_indicators,
                                          const bool* const query_indicators,
                                          bool* const out_query_indicators,
                                          const uint32_t k,
                                          const bool radius_search,
                                          const double radius,
                                          const bool accept_partial,
                                          iscc_Arci* const out_nn_ref,
                                          iscc_Dpid* const out_nn_indices)
{
	assert(nn_search_object != NULL);
	scc_DataSetObject* const data_set_object = nn_search_object->data_set_object;
	#ifndef NDEBUG
		const size_t len_search_indices = nn_search_object->len_search_indices;
		assert(len_search_indices <= INT_MAX);
	#endif
	const iscc_Dpid* const search_indices = nn_search_object->search_indices;
	ANNpointSet* const search_tree = nn_search_object->search_tree;

	assert(iscc_open_search_objects > 0);
	assert(iscc_check_data_set_object(data_set_object, len_query_indicators));
	assert(len_search_indices > 0);
	assert(search_tree != NULL);
	assert(len_query_indicators > 0);
	assert(k > 0);
	assert(k <= len_search_indices);
	assert(!radius_search || (radius > 0.0));
	assert(out_nn_ref != NULL);
	assert(out_nn_indices != NULL);

	if (k > INT_MAX) return false;
	const int k_int = static_cast<int>(k);

	out_nn_ref[0] = 0;

	#if defined(SCC_DPID_INT) || defined(SCC_ANN_INT_OVERRIDE)

		if (search_indices == NULL) {

			// If `iscc_Dpid` is `int` we can send `out_nn_indices` directly to ANN.
			// When searching on sequential indices (i.e., `search_indices == NULL`),
			// ANN produces the desired result and we do not need to do translating and/or casting

			ANNdist* dist_scratch;
			try {
				dist_scratch = new ANNdist[k];
			} catch (...) {
				return false;
			}

			if (!radius_search) {
				for (size_t q = 0; q < len_query_indicators; ++q) {
					out_nn_ref[q + 1] = out_nn_ref[q];
					if ((query_indicators == NULL) || query_indicators[q]) {
						const ANNpoint query_point = data_set_object->data_matrix + q * data_set_object->num_dimensions;
						search_tree->annkSearch(query_point,    // pointer to query point
						                        k_int,          // number of neighbors
						                        out_nn_indices + out_nn_ref[q],    // pointer to start of index result
						                        dist_scratch,   // pointer to start of distance result
						                        SCC_ANN_EPS);   // error margin
						out_nn_ref[q + 1] += k;
					}
				}

			} else {
				assert(radius_search);
				const double radius_sq = radius * radius;
				for (size_t q = 0; q < len_query_indicators; ++q) {
					out_nn_ref[q + 1] = out_nn_ref[q];
					if ((query_indicators == NULL) || query_indicators[q]) {
						const ANNpoint query_point = data_set_object->data_matrix + q * data_set_object->num_dimensions;
						const int num_found = search_tree->annkFRSearch(query_point,              // pointer to query point
						                                                radius_sq,                // squared caliper
						                                                k_int,                    // number of neighbors
						                                                out_nn_indices + out_nn_ref[q], // pointer to start of index result
						                                                dist_scratch,             // pointer to start of distance result
						                                                SCC_ANN_EPS);             // error margin
						assert(num_found >= 0);
						if (num_found >= k_int) {
							out_nn_ref[q + 1] += k;
						} else if (accept_partial) {
							out_nn_ref[q + 1] += static_cast<iscc_Arci>(num_found);
						} else if (out_query_indicators != NULL) {
							out_query_indicators[q] = false;
						}
					}
				}
			}

			delete[] dist_scratch;

			return true;
		}

	#endif // #if defined(SCC_DPID_INT) || defined(SCC_ANN_INT_OVERRIDE)

	ANNidx* idx_scratch;
	ANNdist* dist_scratch;
	try {
		idx_scratch = new ANNidx[k];
		dist_scratch = new ANNdist[k];
	} catch (...) {
		return false;
	}

	if (!radius_search) {
		iscc_Dpid* write_nnidx = out_nn_indices;
		const ANNidx* const idx_stop = idx_scratch + k;
		for (size_t q = 0; q < len_query_indicators; ++q) {
			out_nn_ref[q + 1] = out_nn_ref[q];
			if ((query_indicators == NULL) || query_indicators[q]) {
				const ANNpoint query_point = data_set_object->data_matrix + q * data_set_object->num_dimensions;
				search_tree->annkSearch(query_point,    // pointer to query point
				                        k_int,          // number of neighbors
				                        idx_scratch,    // pointer to start of index result
				                        dist_scratch,   // pointer to start of distance result
				                        SCC_ANN_EPS);   // error margin
				if (search_indices == NULL) {
					// Sequential indices, just do casting
					for (const ANNidx* idx_tmp = idx_scratch; idx_tmp != idx_stop; ++idx_tmp, ++write_nnidx) {
						*write_nnidx = static_cast<iscc_Dpid>(*idx_tmp);
					}
				} else {
					// Not sequential indices, translate to original indices
					for (const ANNidx* idx_tmp = idx_scratch; idx_tmp != idx_stop; ++idx_tmp, ++write_nnidx) {
						*write_nnidx = search_indices[*idx_tmp];
					}
				}
				out_nn_ref[q + 1] += k;
			}
		}

	} else { // radius_search
		const double radius_sq = radius * radius;
		iscc_Dpid* write_nnidx = out_nn_indices;
		for (size_t q = 0; q < len_query_indicators; ++q) {
			out_nn_ref[q + 1] = out_nn_ref[q];
			if ((query_indicators == NULL) || query_indicators[q]) {
				const ANNpoint query_point = data_set_object->data_matrix + q * data_set_object->num_dimensions;
				int num_found = search_tree->annkFRSearch(query_point,     // pointer to query point
				                                          radius_sq,       // squared caliper
				                                          k_int,           // number of neighbors
				                                          idx_scratch,     // pointer to start of index result
				                                          dist_scratch,    // pointer to start of distance result
				                                          SCC_ANN_EPS);    // error margin

				assert(num_found >= 0);
				if (accept_partial || (num_found >= k_int)) {
					if (num_found >= k_int) num_found = k_int;
					const ANNidx* const idx_stop = idx_scratch + num_found;
					if (search_indices == NULL) {
						// Sequential indices, just do casting
						for (const ANNidx* idx_tmp = idx_scratch; idx_tmp != idx_stop; ++idx_tmp, ++write_nnidx) {
							*write_nnidx = static_cast<iscc_Dpid>(*idx_tmp);
						}
					} else {
						// Not sequential indices, translate to original indices
						for (const ANNidx* idx_tmp = idx_scratch; idx_tmp != idx_stop; ++idx_tmp, ++write_nnidx) {
							*write_nnidx = search_indices[*idx_tmp];
						}
					}
					out_nn_ref[q + 1] += static_cast<iscc_Arci>(num_found);
				} else if (out_query_indicators != NULL) {
					out_query_indicators[q] = false;
				}
			}
		}
	}

	delete[] idx_scratch;
	delete[] dist_scratch;

	return true;
}


bool iscc_nearest_neighbor_search_index(iscc_NNSearchObject* const nn_search_object,
                                        const size_t len_query_indices,
                                        const iscc_Dpid* const query_indices,
                                        const uint32_t k,
                                        const bool radius_search,
                                        const double radius,
                                        iscc_Dpid* const out_nn_indices)
{
	assert(nn_search_object != NULL);
	scc_DataSetObject* const data_set_object = nn_search_object->data_set_object;
	#ifndef NDEBUG
		const size_t len_search_indices = nn_search_object->len_search_indices;
		assert(len_search_indices <= INT_MAX);
	#endif
	const iscc_Dpid* const search_indices = nn_search_object->search_indices;
	ANNpointSet* const search_tree = nn_search_object->search_tree;

	assert(iscc_open_search_objects > 0);
	assert(iscc_check_data_set_object(data_set_object, 1));
	assert(len_search_indices > 0);
	assert(search_tree != NULL);
	assert(len_query_indices > 0);
	assert(query_indices != NULL);
	assert(k > 0);
	assert(k <= len_search_indices);
	assert(!radius_search || (radius > 0.0));
	assert(out_nn_indices != NULL);

	if (k > INT_MAX) return false;
	const int k_int = static_cast<int>(k);

	#if defined(SCC_DPID_INT) || defined(SCC_ANN_INT_OVERRIDE)

		if (search_indices == NULL) {

			iscc_Dpid* write_nnidx = out_nn_indices;

			// If `iscc_Dpid` is `int` we can send `out_nn_indices` directly to ANN.
			// When searching on sequential indices (i.e., `search_indices == NULL`),
			// ANN produces the desired result and we do not need to do translating and/or casting

			ANNdist* dist_scratch;
			try {
				dist_scratch = new ANNdist[k];
			} catch (...) {
				return false;
			}

			if (!radius_search) {
				for (size_t q = 0; q < len_query_indices; ++q) {
					const ANNpoint query_point = data_set_object->data_matrix + query_indices[q] * data_set_object->num_dimensions;
					search_tree->annkSearch(query_point,    // pointer to query point
					                        k_int,          // number of neighbors
					                        write_nnidx,    // pointer to start of index result
					                        dist_scratch,   // pointer to start of distance result
					                        SCC_ANN_EPS);   // error margin
					write_nnidx += k;
				}

			} else {
				assert(radius_search);
				const double radius_sq = radius * radius;
				for (size_t q = 0; q < len_query_indices; ++q) {
					const ANNpoint query_point = data_set_object->data_matrix + query_indices[q] * data_set_object->num_dimensions;
					int num_found = search_tree->annkFRSearch(query_point,     // pointer to query point
					                                          radius_sq,       // squared caliper
					                                          k_int,           // number of neighbors
					                                          write_nnidx,     // pointer to start of index result
					                                          dist_scratch,    // pointer to start of distance result
					                                          SCC_ANN_EPS);    // error margin
					assert(num_found >= 0);
					for (; num_found < k_int; ++num_found) {
						write_nnidx[num_found] = ISCC_DPID_NA;
					}
					write_nnidx += k;
				}
			}

			delete[] dist_scratch;

			return true;
		}

	#endif // #if defined(SCC_DPID_INT) || defined(SCC_ANN_INT_OVERRIDE)

	ANNidx* idx_scratch;
	ANNdist* dist_scratch;
	try {
		idx_scratch = new ANNidx[k];
		dist_scratch = new ANNdist[k];
	} catch (...) {
		return false;
	}

	if (!radius_search) {
		iscc_Dpid* write_nnidx = out_nn_indices;
		const ANNidx* const idx_stop = idx_scratch + k;
		for (size_t q = 0; q < len_query_indices; ++q) {
			const ANNpoint query_point = data_set_object->data_matrix + query_indices[q] * data_set_object->num_dimensions;
			search_tree->annkSearch(query_point,    // pointer to query point
			                        k_int,          // number of neighbors
			                        idx_scratch,    // pointer to start of index result
			                        dist_scratch,   // pointer to start of distance result
			                        SCC_ANN_EPS);   // error margin
			if (search_indices == NULL) {
				// Sequential indices, just do casting
				for (const ANNidx* idx_tmp = idx_scratch; idx_tmp != idx_stop; ++idx_tmp, ++write_nnidx) {
					*write_nnidx = static_cast<iscc_Dpid>(*idx_tmp);
				}
			} else {
				// Not sequential indices, translate to original indices
				for (const ANNidx* idx_tmp = idx_scratch; idx_tmp != idx_stop; ++idx_tmp, ++write_nnidx) {
					*write_nnidx = search_indices[*idx_tmp];
				}
			}
		}

	} else { // radius_search
		const double radius_sq = radius * radius;
		iscc_Dpid* write_nnidx = out_nn_indices;
		const ANNidx* const idx_stop = idx_scratch + k;
		for (size_t q = 0; q < len_query_indices; ++q) {
			const ANNpoint query_point = data_set_object->data_matrix + query_indices[q] * data_set_object->num_dimensions;
			int num_found = search_tree->annkFRSearch(query_point,     // pointer to query point
			                                          radius_sq,       // squared caliper
			                                          k_int,           // number of neighbors
			                                          idx_scratch,     // pointer to start of index result
			                                          dist_scratch,    // pointer to start of distance result
			                                          SCC_ANN_EPS);    // error margin

			assert(num_found >= 0);
			const ANNidx* idx_tmp = idx_scratch;
			if (num_found >= k_int) num_found = k_int;
			const ANNidx* const idx_found_stop = idx_scratch + num_found;
			if (search_indices == NULL) {
				// Sequential indices, just do casting
				for (; idx_tmp != idx_found_stop; ++idx_tmp, ++write_nnidx) {
					*write_nnidx = static_cast<iscc_Dpid>(*idx_tmp);
				}
			} else {
				// Not sequential indices, translate to original indices
				for (; idx_tmp != idx_found_stop; ++idx_tmp, ++write_nnidx) {
					*write_nnidx = search_indices[*idx_tmp];
				}
			}
			for (; idx_tmp != idx_stop; ++idx_tmp, ++write_nnidx) {
				*write_nnidx = ISCC_DPID_NA;
			}
		}
	}

	delete[] idx_scratch;
	delete[] dist_scratch;

	return true;
}


bool iscc_close_nn_search_object(iscc_NNSearchObject** const nn_search_object)
{
	assert(iscc_open_search_objects >= 0);

	if (nn_search_object != NULL) {
		delete (*nn_search_object)->search_tree;
		delete[] (*nn_search_object)->search_points;
		delete *nn_search_object;
		*nn_search_object = NULL;
	}

	if (iscc_open_search_objects <= 0) {
		annClose();
		return false;
	}

	--iscc_open_search_objects;
	if (iscc_open_search_objects == 0) annClose();

	return true;
}
