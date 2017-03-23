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

#include "dist_search_imp.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include "../include/scclust.h"
#include "data_set_struct.h"
#include "scclust_types.h"


// =============================================================================
// Distance calculations
// =============================================================================

static inline double iscc_get_sq_dist(const scc_DataSet* const data_set,
                                      const size_t index1,
                                      const size_t index2)
{
	assert(index1 < data_set->num_data_points);
	assert(index2 < data_set->num_data_points);

	const double* data1 = &data_set->data_matrix[index1 * data_set->num_dimensions];
	const double* const data1_stop = data1 + data_set->num_dimensions;
	const double* data2 = &data_set->data_matrix[index2 * data_set->num_dimensions];

	double tmp_dist = 0.0;
	while (data1 != data1_stop) {
		const double value_diff = (*data1 - *data2);
		++data1;
		++data2;
		tmp_dist += value_diff * value_diff;
	}
	return tmp_dist;
}


// =============================================================================
// Miscellaneous functions implementations
// =============================================================================

bool iscc_imp_check_data_set(void* const data_set)
{
	if (data_set == NULL) return false;
	const scc_DataSet* const data_set_cast = (const scc_DataSet*) data_set;
	if (!scc_is_initialized_data_set(data_set_cast)) return false;
	return true;
}


size_t iscc_imp_num_data_points(void* const data_set)
{
	assert(iscc_imp_check_data_set(data_set));
	if (data_set == NULL) return 0;
	const scc_DataSet* const data_set_cast = (const scc_DataSet*) data_set;
	return data_set_cast->num_data_points;
}


bool iscc_imp_get_dist_matrix(void* const data_set,
                              const size_t len_point_indices,
                              const scc_PointIndex point_indices[const],
                              double output_dists[])
{
	assert(iscc_imp_check_data_set(data_set));
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


bool iscc_imp_get_dist_rows(void* const data_set,
                            const size_t len_query_indices,
                            const scc_PointIndex query_indices[const],
                            const size_t len_column_indices,
                            const scc_PointIndex column_indices[const],
                            double output_dists[])
{
	assert(iscc_imp_check_data_set(data_set));
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
// Max dist functions implementations
// =============================================================================

struct iscc_MaxDistObject {
	int32_t max_dist_version;
	scc_DataSet* data_set;
	size_t len_search_indices;
	const scc_PointIndex* search_indices;
};


static const int32_t ISCC_MAXDIST_STRUCT_VERSION = 722439001;


bool iscc_imp_init_max_dist_object(void* const data_set,
                                   const size_t len_search_indices,
                                   const scc_PointIndex search_indices[const],
                                   iscc_MaxDistObject** const out_max_dist_object)
{
	assert(iscc_imp_check_data_set(data_set));
	assert(len_search_indices > 0);
	assert(out_max_dist_object != NULL);

	*out_max_dist_object = malloc(sizeof(iscc_MaxDistObject));
	if (*out_max_dist_object == NULL) return false;

	**out_max_dist_object = (iscc_MaxDistObject) {
		.max_dist_version = ISCC_MAXDIST_STRUCT_VERSION,
		.data_set = data_set,
		.len_search_indices = len_search_indices,
		.search_indices = search_indices,
	};

	return true;
}


bool iscc_imp_get_max_dist(iscc_MaxDistObject* const max_dist_object,
                           const size_t len_query_indices,
                           const scc_PointIndex query_indices[const],
                           scc_PointIndex out_max_indices[const],
                           double out_max_dists[const])
{
	assert(max_dist_object != NULL);
	assert(max_dist_object->max_dist_version == ISCC_MAXDIST_STRUCT_VERSION);
	scc_DataSet* const data_set = max_dist_object->data_set;
	const size_t len_search_indices = max_dist_object->len_search_indices;
	const scc_PointIndex* const search_indices = max_dist_object->search_indices;

	assert(iscc_imp_check_data_set(data_set));
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
					out_max_indices[q] = (scc_PointIndex) s;
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
					out_max_indices[q] = (scc_PointIndex) s;
				}
			}
			out_max_dists[q] = sqrt(max_dist);
		}
	}

	return true;
}


bool iscc_imp_close_max_dist_object(iscc_MaxDistObject** const max_dist_object)
{
	if (max_dist_object != NULL && *max_dist_object != NULL) {
		assert((*max_dist_object)->max_dist_version == ISCC_MAXDIST_STRUCT_VERSION);
		free(*max_dist_object);
		*max_dist_object = NULL;
	}
	return true;
}


// =============================================================================
// Nearest neighbor search functions implementations
// =============================================================================

struct iscc_NNSearchObject {
	int32_t nn_search_version;
	scc_DataSet* data_set;
	size_t len_search_indices;
	const scc_PointIndex* search_indices;
};


static const int32_t ISCC_NN_SEARCH_STRUCT_VERSION = 722294001;


static inline void iscc_add_dist_to_list(const double add_dist,
                                         const scc_PointIndex add_index,
                                         double* dist_list,
                                         scc_PointIndex* index_list,
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


bool iscc_imp_init_nn_search_object(void* const data_set,
                                    const size_t len_search_indices,
                                    const scc_PointIndex search_indices[const],
                                    iscc_NNSearchObject** const out_nn_search_object)
{
	assert(iscc_imp_check_data_set(data_set));
	assert(len_search_indices > 0);
	assert(out_nn_search_object != NULL);

	*out_nn_search_object = malloc(sizeof(iscc_NNSearchObject));
	if (*out_nn_search_object == NULL) return false;

	**out_nn_search_object = (iscc_NNSearchObject) {
		.nn_search_version = ISCC_NN_SEARCH_STRUCT_VERSION,
		.data_set = data_set,
		.len_search_indices = len_search_indices,
		.search_indices = search_indices,
	};

	return true;
}


bool iscc_imp_nearest_neighbor_search(iscc_NNSearchObject* const nn_search_object,
                                      const size_t len_query_indices,
                                      const scc_PointIndex query_indices[const],
                                      const uint32_t k,
                                      const bool radius_search,
                                      const double radius,
                                      size_t* const out_num_ok_queries,
                                      scc_PointIndex out_query_indices[const],
                                      scc_PointIndex out_nn_indices[const])
{
	assert(nn_search_object != NULL);
	assert(nn_search_object->nn_search_version == ISCC_NN_SEARCH_STRUCT_VERSION);
	scc_DataSet* const data_set = nn_search_object->data_set;
	const size_t len_search_indices = nn_search_object->len_search_indices;
	const scc_PointIndex* const search_indices = nn_search_object->search_indices;

	assert(iscc_imp_check_data_set(data_set));
	assert(len_search_indices > 0);
	assert(len_query_indices > 0);
	assert(k > 0);
	assert(k <= len_search_indices);
	assert(!radius_search || (radius > 0.0));
	assert(out_num_ok_queries != NULL);
	assert(out_nn_indices != NULL);

	double tmp_dist;
	size_t num_ok_queries = 0;
	scc_PointIndex* index_write = out_nn_indices;
	double* const sort_scratch = malloc(sizeof(double[k]));
	if (sort_scratch == NULL) return false;
	double* const sort_scratch_end = sort_scratch + k - 1;
	const double radius_sq = radius * radius;

	if (search_indices == NULL) {
		for (size_t q = 0; q < len_query_indices; ++q) {
			size_t query = q;
			if (query_indices != NULL) {
				query = (size_t) query_indices[q];
			}
			size_t s = 0;
			uint32_t found;
			scc_PointIndex* const index_write_end = index_write + k - 1;

			if (radius_search) {
				found = 0;
				for (; (s < len_search_indices) && (found < k); ++s) {
					tmp_dist = iscc_get_sq_dist(data_set, query, s);
					if (tmp_dist > radius_sq) continue;
					iscc_add_dist_to_list(tmp_dist, (scc_PointIndex) s, sort_scratch + found, index_write + found, sort_scratch);
					++found;
				}
			} else {
				for (; s < k; ++s) {
					tmp_dist = iscc_get_sq_dist(data_set, query, s);
					iscc_add_dist_to_list(tmp_dist, (scc_PointIndex) s, sort_scratch + s, index_write + s, sort_scratch);
				}
				found = k;
			}

			for (; s < len_search_indices; ++s) {
				assert(found == k);
				tmp_dist = iscc_get_sq_dist(data_set, query, s);
				if (tmp_dist >= *sort_scratch_end) continue;
				iscc_add_dist_to_list(tmp_dist, (scc_PointIndex) s, sort_scratch_end, index_write_end, sort_scratch);
			}

			assert(found == k || out_query_indices != NULL);
			if (found == k) {
				if (out_query_indices != NULL) {
					out_query_indices[num_ok_queries] = (scc_PointIndex) query;
				}
				++num_ok_queries;
				index_write += k;
			}
		}
	} else {
		assert(search_indices != NULL);
		for (size_t q = 0; q < len_query_indices; ++q) {
			size_t query = q;
			if (query_indices != NULL) {
				query = (size_t) query_indices[q];
			}
			size_t s = 0;
			uint32_t found;
			scc_PointIndex* const index_write_end = index_write + k - 1;

			if (radius_search) {
				found = 0;
				for (; (s < len_search_indices) && (found < k); ++s) {
					tmp_dist = iscc_get_sq_dist(data_set, query, (size_t) search_indices[s]);
					if (tmp_dist > radius_sq) continue;
					iscc_add_dist_to_list(tmp_dist, search_indices[s], sort_scratch + found, index_write + found, sort_scratch);
					++found;
				}
			} else {
				for (; s < k; ++s) {
					tmp_dist = iscc_get_sq_dist(data_set, query, (size_t) search_indices[s]);
					iscc_add_dist_to_list(tmp_dist, search_indices[s], sort_scratch + s, index_write + s, sort_scratch);
				}
				found = k;
			}

			for (; s < len_search_indices; ++s) {
				assert(found == k);
				tmp_dist = iscc_get_sq_dist(data_set, query, (size_t) search_indices[s]);
				if (tmp_dist >= *sort_scratch_end) continue;
				iscc_add_dist_to_list(tmp_dist, search_indices[s], sort_scratch_end, index_write_end, sort_scratch);
			}

			assert(found == k || out_query_indices != NULL);
			if (found == k) {
				if (out_query_indices != NULL) {
					out_query_indices[num_ok_queries] = (scc_PointIndex) query;
				}
				++num_ok_queries;
				index_write += k;
			}
		}
	}

	*out_num_ok_queries = num_ok_queries;

	free(sort_scratch);

	return true;
}


bool iscc_imp_close_nn_search_object(iscc_NNSearchObject** const nn_search_object)
{
	if (nn_search_object != NULL && *nn_search_object != NULL) {
		assert((*nn_search_object)->nn_search_version == ISCC_NN_SEARCH_STRUCT_VERSION);
		free(*nn_search_object);
		*nn_search_object = NULL;
	}
	return true;
}
