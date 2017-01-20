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

#include "ann_wrapper.h"

#include <cassert>
#include <climits>
#include <cstddef>
#include <include/ANN/ANN.h>
#include <include/scclust.h>
#include <include/scclust_spi.h>
#include <src/data_set_struct.h>
#include <src/dist_search_imp.h>


#ifdef SCC_ANN_BDTREE
	#define ANNpointSetConstructor ANNbd_tree
#else
	#define ANNpointSetConstructor ANNkd_tree
#endif

#ifndef SCC_ANN_EPS
	#define SCC_ANN_EPS 0.0
#endif


// =============================================================================
// Internal functions
// =============================================================================

extern "C" {

	bool iscc_ann_init_nn_search_object(void* data_set,
	                                    size_t len_search_indices,
	                                    const scc_PointIndex search_indices[],
	                                    iscc_NNSearchObject** out_nn_search_object);

	bool iscc_ann_nearest_neighbor_search(iscc_NNSearchObject* nn_search_object,
	                                      size_t len_query_indices,
	                                      const scc_PointIndex query_indices[],
	                                      uint32_t k,
	                                      bool radius_search,
	                                      double radius,
	                                      size_t* out_num_ok_queries,
	                                      scc_PointIndex out_query_indices[],
	                                      scc_PointIndex out_nn_indices[]);

	bool iscc_ann_close_nn_search_object(iscc_NNSearchObject** nn_search_object);

}


// =============================================================================
// Internal structs and variables
// =============================================================================

static int iscc_ann_open_search_objects = 0;

static const int32_t ISCC_ANN_NN_SEARCH_STRUCT_VERSION = 155294001;

struct iscc_NNSearchObject {
	int32_t nn_search_version;
	scc_DataSet* data_set;
	size_t len_search_indices;
	const scc_PointIndex* search_indices;
	ANNpoint* search_points;
	ANNpointSet* search_tree;
};


// =============================================================================
// External function implementations
// =============================================================================

bool scc_set_ann_dist_search()
{
	return scc_set_dist_functions(NULL, NULL, NULL, NULL, NULL, NULL,
	                              iscc_ann_init_nn_search_object,
	                              iscc_ann_nearest_neighbor_search,
	                              iscc_ann_close_nn_search_object);
}


// =============================================================================
// Internal function implementations
// =============================================================================

bool iscc_ann_init_nn_search_object(void* const data_set,
                                    const size_t len_search_indices,
                                    const scc_PointIndex* const search_indices,
                                    iscc_NNSearchObject** const out_nn_search_object)
{
	assert(iscc_ann_open_search_objects >= 0);
	assert(len_search_indices > 0);
	assert(out_nn_search_object != NULL);
	assert(iscc_imp_check_data_set(data_set, len_search_indices));

	if (len_search_indices > INT_MAX) return false;

	scc_DataSet* const data_set_cast = static_cast<scc_DataSet*>(data_set);

	ANNpoint* search_points;
	try {
		search_points = new ANNpoint[len_search_indices];
		*out_nn_search_object = new iscc_NNSearchObject;
	} catch (...) {
		return false;
	}

	if (search_indices == NULL) {
		assert(len_search_indices <= data_set_cast->num_data_points);
		double* search_point = const_cast<double*>(data_set_cast->data_matrix);
		for (size_t i = 0; i < len_search_indices; ++i, search_point += data_set_cast->num_dimensions) {
			search_points[i] = search_point;
		}
	} else if (search_indices != NULL) {
		for (size_t i = 0; i < len_search_indices; ++i) {
			assert(static_cast<size_t>(search_indices[i]) < data_set_cast->num_data_points);
			search_points[i] = const_cast<double*>(data_set_cast->data_matrix) + search_indices[i] * data_set_cast->num_dimensions;
		}
	}

	ANNpointSet* search_tree;
	try {
		search_tree = new ANNpointSetConstructor(search_points,
		                                         static_cast<int>(len_search_indices),
		                                         static_cast<int>(data_set_cast->num_dimensions));
	} catch (...) {
		delete[] search_points;
		delete *out_nn_search_object;
		*out_nn_search_object = NULL;
		return false;
	}

	(*out_nn_search_object)->nn_search_version = ISCC_ANN_NN_SEARCH_STRUCT_VERSION;
	(*out_nn_search_object)->data_set = data_set_cast;
	(*out_nn_search_object)->len_search_indices = len_search_indices;
	(*out_nn_search_object)->search_indices = search_indices;
	(*out_nn_search_object)->search_points = search_points;
	(*out_nn_search_object)->search_tree = search_tree;

	++iscc_ann_open_search_objects;
	return true;
}


bool iscc_ann_nearest_neighbor_search(iscc_NNSearchObject* const nn_search_object,
                                      const size_t len_query_indices,
                                      const scc_PointIndex* const query_indices,
                                      const uint32_t k,
                                      const bool radius_search,
                                      const double radius,
                                      size_t* const out_num_ok_queries,
                                      scc_PointIndex* const out_query_indices,
                                      scc_PointIndex* const out_nn_indices)
{

	assert(nn_search_object != NULL);
	assert(nn_search_object->nn_search_version == ISCC_ANN_NN_SEARCH_STRUCT_VERSION);
	scc_DataSet* const data_set = nn_search_object->data_set;
	#ifndef NDEBUG
		const size_t len_search_indices = nn_search_object->len_search_indices;
		assert(len_search_indices <= INT_MAX);
	#endif
	const scc_PointIndex* const search_indices = nn_search_object->search_indices;
	ANNpointSet* const search_tree = nn_search_object->search_tree;

	assert(iscc_ann_open_search_objects > 0);
	assert(iscc_imp_check_data_set(data_set, 0));
	assert(len_search_indices > 0);
	assert(search_tree != NULL);
	assert(len_query_indices > 0);
	assert(k > 0);
	assert(k <= len_search_indices);
	assert(!radius_search || (radius > 0.0));
	assert(out_num_ok_queries != NULL);
	assert(out_nn_indices != NULL);

	if (k > INT_MAX) return false;
	const int k_int = static_cast<int>(k);

	#ifdef SCC_M_POINTINDEX_TYPE_int

		if (search_indices == NULL) {

			// If `scc_PointIndex` is `int` we can send `out_nn_indices` directly to ANN.
			// When searching on sequential indices (i.e., `search_indices == NULL`),
			// ANN produces the desired result and we do not need to do translating and/or casting

			ANNdist* dist_scratch;
			try {
				dist_scratch = new ANNdist[k];
			} catch (...) {
				return false;
			}

			size_t num_ok_queries = 0;

			if (!radius_search) {
				int* write_nnidx = out_nn_indices;
				for (size_t q = 0; q < len_query_indices; ++q) {
					size_t query = q;
					if (query_indices != NULL) {
						query = (size_t) query_indices[q];
					}
					const ANNpoint query_point = const_cast<double*>(data_set->data_matrix) + query * data_set->num_dimensions;
					search_tree->annkSearch(query_point,    // pointer to query point
					                        k_int,          // number of neighbors
					                        write_nnidx,    // pointer to start of index result
					                        dist_scratch,   // pointer to start of distance result
					                        SCC_ANN_EPS);   // error margin
					if (out_query_indices != NULL) {
						out_query_indices[num_ok_queries] = (scc_PointIndex) query;
					}
					++num_ok_queries;
					write_nnidx += k_int;
				}
			} else {
				assert(radius_search);
				const double radius_sq = radius * radius;
				int* write_nnidx = out_nn_indices;
				for (size_t q = 0; q < len_query_indices; ++q) {
					size_t query = q;
					if (query_indices != NULL) {
						query = (size_t) query_indices[q];
					}
					const ANNpoint query_point = const_cast<double*>(data_set->data_matrix) + query * data_set->num_dimensions;
					const int num_found = search_tree->annkFRSearch(query_point,              // pointer to query point
					                                                radius_sq,                // squared caliper
					                                                k_int,                    // number of neighbors
					                                                write_nnidx,              // pointer to start of index result
					                                                dist_scratch,             // pointer to start of distance result
					                                                SCC_ANN_EPS);             // error margin
					assert(num_found >= 0);
					if (num_found >= k_int) {
						if (out_query_indices != NULL) {
							out_query_indices[num_ok_queries] = (scc_PointIndex) query;
						}
						++num_ok_queries;
						write_nnidx += k_int;
					}
			}

			delete[] dist_scratch;

			*out_num_ok_queries = num_ok_queries;

			return true;
		}
	}

	#endif // #ifdef SCC_M_POINTINDEX_TYPE_int


	ANNidx* idx_scratch;
	ANNdist* dist_scratch;
	try {
		idx_scratch = new ANNidx[k];
		dist_scratch = new ANNdist[k];
	} catch (...) {
		return false;
	}

	size_t num_ok_queries = 0;

	if (!radius_search) {
		scc_PointIndex* write_nnidx = out_nn_indices;
		const ANNidx* const idx_stop = idx_scratch + k;
		for (size_t q = 0; q < len_query_indices; ++q) {
			size_t query = q;
			if (query_indices != NULL) {
				query = (size_t) query_indices[q];
			}
			const ANNpoint query_point = const_cast<double*>(data_set->data_matrix) + query * data_set->num_dimensions;
			search_tree->annkSearch(query_point,    // pointer to query point
			                        k_int,          // number of neighbors
			                        idx_scratch,    // pointer to start of index result
			                        dist_scratch,   // pointer to start of distance result
			                        SCC_ANN_EPS);   // error margin
			if (search_indices == NULL) {
				// Sequential indices, just do casting
				for (const ANNidx* idx_tmp = idx_scratch; idx_tmp != idx_stop; ++idx_tmp, ++write_nnidx) {
					*write_nnidx = static_cast<scc_PointIndex>(*idx_tmp);
				}
			} else {
				// Not sequential indices, translate to original indices
				for (const ANNidx* idx_tmp = idx_scratch; idx_tmp != idx_stop; ++idx_tmp, ++write_nnidx) {
					*write_nnidx = search_indices[*idx_tmp];
				}
			}
			if (out_query_indices != NULL) {
				out_query_indices[num_ok_queries] = (scc_PointIndex) query;
			}
			++num_ok_queries;
		}

	} else {
		assert(radius_search);
		const double radius_sq = radius * radius;
		scc_PointIndex* write_nnidx = out_nn_indices;
		for (size_t q = 0; q < len_query_indices; ++q) {
			size_t query = q;
			if (query_indices != NULL) {
				query = (size_t) query_indices[q];
			}
			const ANNpoint query_point = const_cast<double*>(data_set->data_matrix) + query * data_set->num_dimensions;
			int num_found = search_tree->annkFRSearch(query_point,     // pointer to query point
			                                          radius_sq,       // squared caliper
			                                          k_int,           // number of neighbors
			                                          idx_scratch,     // pointer to start of index result
			                                          dist_scratch,    // pointer to start of distance result
			                                          SCC_ANN_EPS);    // error margin
			assert(num_found >= 0);
			if (num_found >= k_int) {
				const ANNidx* const idx_stop = idx_scratch + k;
				if (search_indices == NULL) {
					// Sequential indices, just do casting
					for (const ANNidx* idx_tmp = idx_scratch; idx_tmp != idx_stop; ++idx_tmp, ++write_nnidx) {
						*write_nnidx = static_cast<scc_PointIndex>(*idx_tmp);
					}
				} else {
					// Not sequential indices, translate to original indices
					for (const ANNidx* idx_tmp = idx_scratch; idx_tmp != idx_stop; ++idx_tmp, ++write_nnidx) {
						*write_nnidx = search_indices[*idx_tmp];
					}
				}
				if (out_query_indices != NULL) {
					out_query_indices[num_ok_queries] = (scc_PointIndex) query;
				}
				++num_ok_queries;
			}
		}
	}

	delete[] idx_scratch;
	delete[] dist_scratch;

	*out_num_ok_queries = num_ok_queries;

	return true;
}


bool iscc_ann_close_nn_search_object(iscc_NNSearchObject** const nn_search_object)
{
	assert(iscc_ann_open_search_objects >= 0);

	if (nn_search_object != NULL && *nn_search_object != NULL) {
		assert((*nn_search_object)->nn_search_version == ISCC_ANN_NN_SEARCH_STRUCT_VERSION);
		delete (*nn_search_object)->search_tree;
		delete[] (*nn_search_object)->search_points;
		delete *nn_search_object;
		*nn_search_object = NULL;
	}

	if (iscc_ann_open_search_objects <= 0) {
		annClose();
		return false;
	}

	--iscc_ann_open_search_objects;
	if (iscc_ann_open_search_objects == 0) {
		annClose();
	}

	return true;
}
