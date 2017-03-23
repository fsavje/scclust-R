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

#include "nng_core.h"

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "../include/scclust.h"
#include "clustering_struct.h"
#include "digraph_core.h"
#include "digraph_operations.h"
#include "dist_search.h"
#include "error.h"
#include "nng_findseeds.h"
#include "scclust_types.h"


// =============================================================================
// Internal structs & variables
// =============================================================================

typedef struct iscc_TypeCount {
	uint32_t sum_type_constraints;
	size_t* type_group_size;
	scc_PointIndex* point_store;
	scc_PointIndex** type_groups;
} iscc_TypeCount;


static const size_t ISCC_ESTIMATE_AVG_MAX_SAMPLE = 1000;


// =============================================================================
// Static function prototypes
// =============================================================================

static scc_ErrorCode iscc_make_nng(void* data_set,
                                   size_t num_data_points,
                                   size_t len_search_indices,
                                   const scc_PointIndex search_indices[],
                                   size_t len_query_indices,
                                   const scc_PointIndex query_indices[],
                                   uint32_t k,
                                   bool radius_search,
                                   double radius,
                                   size_t* out_len_query_indices,
                                   scc_PointIndex out_query_indices[],
                                   iscc_Digraph* out_nng);


static scc_ErrorCode iscc_make_nng_from_search_object(iscc_NNSearchObject* nn_search_object,
                                                      size_t num_data_points,
                                                      size_t len_query_indices,
                                                      const scc_PointIndex query_indices[],
                                                      uint32_t k,
                                                      bool radius_search,
                                                      double radius,
                                                      size_t* out_len_query_indices,
                                                      scc_PointIndex out_query_indices[],
                                                      iscc_Digraph* out_nng);


static inline void iscc_ensure_self_match(iscc_Digraph* nng,
                                          size_t len_search_indices,
                                          const scc_PointIndex search_indices[]);


static scc_ErrorCode iscc_type_count(size_t num_data_points,
                                     uint32_t size_constraint,
                                     uint_fast16_t num_types,
                                     const uint32_t type_constraints[static num_types],
                                     const scc_TypeLabel type_labels[static num_data_points],
                                     iscc_TypeCount* out_type_result);


static size_t iscc_assign_seeds_and_neighbors(scc_Clustering* clustering,
                                              const iscc_SeedResult* seed_result,
                                              iscc_Digraph* nng);


static size_t iscc_assign_by_nng(scc_Clustering* clustering,
                                 iscc_Digraph* nng);


static scc_ErrorCode iscc_assign_by_nn_search(scc_Clustering* clustering,
                                              iscc_NNSearchObject* nn_search_object,
                                              size_t num_to_assign,
                                              scc_PointIndex to_assign[restrict static num_to_assign],
                                              bool radius_constraint,
                                              double radius);


#ifdef SCC_STABLE_NNG

static void iscc_sort_nng(iscc_Digraph* nng);

#endif // ifdef SCC_STABLE_NNG


// =============================================================================
// External function implementations
// =============================================================================

scc_ErrorCode iscc_get_nng_with_size_constraint(void* const data_set,
                                                const size_t num_data_points,
                                                const uint32_t size_constraint,
                                                size_t len_primary_data_points,
                                                const scc_PointIndex primary_data_points[],
                                                const bool radius_constraint,
                                                const double radius,
                                                iscc_Digraph* const out_nng)
{
	assert(iscc_check_data_set(data_set));
	assert(iscc_num_data_points(data_set) == num_data_points);
	assert(num_data_points >= 2);
	assert(size_constraint <= num_data_points);
	assert(size_constraint >= 2);
	assert(!radius_constraint || (radius > 0.0));
	assert(out_nng != NULL);

	size_t num_queries;
	if (primary_data_points == NULL) {
		num_queries = num_data_points;
	} else {
		num_queries = len_primary_data_points;
	}

	scc_ErrorCode ec;
	if ((ec = iscc_make_nng(data_set,
	                        num_data_points,
	                        num_data_points,
	                        NULL,
	                        num_queries,
	                        primary_data_points,
	                        size_constraint,
	                        radius_constraint,
	                        radius,
	                        NULL,
	                        NULL,
	                        out_nng)) != SCC_ER_OK) {
		return ec;
	}

	if (iscc_digraph_is_empty(out_nng)) {
		iscc_free_digraph(out_nng);
		return iscc_make_error_msg(SCC_ER_NO_SOLUTION, "Infeasible radius constraint.");
	}

	iscc_ensure_self_match(out_nng, num_data_points, NULL);

	if ((ec = iscc_delete_loops(out_nng)) != SCC_ER_OK) {
		iscc_free_digraph(out_nng);
		return ec;
	}

	#ifdef SCC_STABLE_NNG
		iscc_sort_nng(out_nng);
	#endif // ifdef SCC_STABLE_NNG

	return iscc_no_error();
}


scc_ErrorCode iscc_get_nng_with_type_constraint(void* const data_set,
                                                const size_t num_data_points,
                                                const uint32_t size_constraint,
                                                const uint_fast16_t num_types,
                                                const uint32_t type_constraints[const static num_types],
                                                const scc_TypeLabel type_labels[const static num_data_points],
                                                size_t len_primary_data_points,
                                                const scc_PointIndex primary_data_points[],
                                                const bool radius_constraint,
                                                const double radius,
                                                iscc_Digraph* const out_nng)
{
	assert(iscc_check_data_set(data_set));
	assert(iscc_num_data_points(data_set) == num_data_points);
	assert(num_data_points >= 2);
	assert(size_constraint <= num_data_points);
	assert(size_constraint >= 2);
	assert(num_types >= 2);
	assert(num_types <= ISCC_TYPELABEL_MAX);
	assert(type_constraints != NULL);
	assert(type_labels != NULL);
	assert(!radius_constraint || (radius > 0.0));
	assert(out_nng != NULL);

	size_t num_queries;
	if (primary_data_points == NULL) {
		num_queries = num_data_points;
	} else {
		num_queries = len_primary_data_points;
	}

	scc_PointIndex* seedable;
	const scc_PointIndex* seedable_const;
	if (radius_constraint) {
		seedable = malloc(sizeof(scc_PointIndex[num_queries]));
		if (seedable == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);
		seedable_const = seedable;
		if (primary_data_points == NULL) {
			for (scc_PointIndex i = 0; i < (scc_PointIndex) num_data_points; ++i) {
				seedable[i] = i;
			}
		} else {
			memcpy(seedable, primary_data_points, sizeof(scc_PointIndex[num_queries]));
		}
	} else {
		seedable = NULL;
		seedable_const = primary_data_points;
	}

	iscc_Digraph* const nng_by_type = malloc(sizeof(iscc_Digraph[num_types]));
	if (nng_by_type == NULL) {
		free(seedable);
		return iscc_make_error(SCC_ER_NO_MEMORY);
	}

	scc_ErrorCode ec;
	iscc_TypeCount tc;
	if ((ec = iscc_type_count(num_data_points,
	                          size_constraint,
	                          num_types,
	                          type_constraints,
	                          type_labels,
	                          &tc)) != SCC_ER_OK) {
		free(seedable);
		free(nng_by_type);
		return ec;
	}

	uint_fast16_t num_non_zero_type_constraints = 0;
	for (uint_fast16_t i = 0; i < num_types; ++i) {
		if (type_constraints[i] > 0) {
			if ((ec = iscc_make_nng(data_set,
			                        num_data_points,
			                        tc.type_group_size[i],
			                        tc.type_groups[i],
			                        num_queries,
			                        seedable_const,
			                        type_constraints[i],
			                        radius_constraint,
			                        radius,
			                        &num_queries,
			                        seedable,
			                        &nng_by_type[num_non_zero_type_constraints])) != SCC_ER_OK) {
				break;
			}
			++num_non_zero_type_constraints;
			if (iscc_digraph_is_empty(&nng_by_type[num_non_zero_type_constraints - 1])) {
				ec = iscc_make_error_msg(SCC_ER_NO_SOLUTION, "Infeasible radius constraint.");
				break;
			}
			iscc_ensure_self_match(&nng_by_type[num_non_zero_type_constraints - 1],
			                       tc.type_group_size[i],
			                       tc.type_groups[i]);
		}
	}

	free(tc.type_group_size);
	free(tc.point_store);
	free(tc.type_groups);

	if (ec == SCC_ER_OK) {
		if (size_constraint > tc.sum_type_constraints) {
			// If general size constaint (besides type constraints), we need to keep self-loops
			ec = iscc_digraph_union_and_delete(num_non_zero_type_constraints, nng_by_type, num_queries, seedable_const, true, out_nng);
		} else {
			ec = iscc_digraph_union_and_delete(num_non_zero_type_constraints, nng_by_type, num_queries, seedable_const, false, out_nng);
		}
	}

	for (uint_fast16_t i = 0; i < num_non_zero_type_constraints; ++i) {
		iscc_free_digraph(&nng_by_type[i]);
	}
	free(nng_by_type);

	if (ec != SCC_ER_OK) {
		// When `ec != SCC_ER_OK`, error is from `iscc_digraph_union_and_delete` so `out_nng` is already freed
		free(seedable);
		return ec;
	}

	if (size_constraint > tc.sum_type_constraints) {
		uint32_t additional_nn_needed = size_constraint - tc.sum_type_constraints;
		iscc_Digraph nng_sum[2];
		nng_sum[0] = *out_nng;

		if ((ec = iscc_make_nng(data_set,
		                        num_data_points,
		                        num_data_points,
		                        NULL,
		                        num_queries,
		                        seedable_const,
		                        size_constraint,
		                        radius_constraint,
		                        radius,
		                        &num_queries,
		                        seedable,
		                        &nng_sum[1])) != SCC_ER_OK) {
			free(seedable);
			iscc_free_digraph(&nng_sum[0]);
			return ec;
		}

		if (iscc_digraph_is_empty(&nng_sum[1])) {
			ec = iscc_make_error_msg(SCC_ER_NO_SOLUTION, "Infeasible radius constraint.");
		} else {
			ec = iscc_digraph_difference(&nng_sum[1], &nng_sum[0], additional_nn_needed);
		}

		if (ec == SCC_ER_OK) {
			ec = iscc_digraph_union_and_delete(2, nng_sum, num_queries, seedable_const, false, out_nng);
		}

		iscc_free_digraph(&nng_sum[0]);
		iscc_free_digraph(&nng_sum[1]);

		if (ec != SCC_ER_OK) {
			free(seedable);
			return ec;
		}
	}

	free(seedable);

	#ifdef SCC_STABLE_NNG
		iscc_sort_nng(out_nng);
	#endif // ifdef SCC_STABLE_NNG

	return iscc_no_error();
}


scc_ErrorCode iscc_estimate_avg_seed_dist(void* const data_set,
                                          const iscc_SeedResult* const seed_result,
                                          const iscc_Digraph* const nng,
                                          const uint32_t size_constraint,
                                          double* const out_avg_seed_dist)
{
	assert(iscc_check_data_set(data_set));
	assert(iscc_num_data_points(data_set) == nng->vertices);
	assert(seed_result->count > 0);
	assert(seed_result->seeds != NULL);
	assert(iscc_digraph_is_valid(nng));
	assert(!iscc_digraph_is_empty(nng));
	assert(size_constraint >= 2);
	assert(out_avg_seed_dist != NULL);

	const size_t step = (seed_result->count > ISCC_ESTIMATE_AVG_MAX_SAMPLE) ? (seed_result->count / ISCC_ESTIMATE_AVG_MAX_SAMPLE) : 1;
	assert(step > 0);

	size_t sampled = 0;
	double sum_dist = 0.0;
	double* const dist_scratch = malloc(sizeof(double[size_constraint]));
	if (dist_scratch == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);

	for (size_t s = 0; s < seed_result->count; s += step) {
		const scc_PointIndex seed = seed_result->seeds[s];
		const size_t num_neighbors = (nng->tail_ptr[seed + 1] - nng->tail_ptr[seed]);
		const scc_PointIndex* const neighbors = nng->head + nng->tail_ptr[seed];

		// Either zero or one self-loops
		assert((num_neighbors == size_constraint) ||
		       (num_neighbors == size_constraint - 1));

		if (!iscc_get_dist_rows(data_set,
		                        1,
		                        &seed,
		                        num_neighbors,
		                        neighbors,
		                        dist_scratch)) {
			free(dist_scratch);
			return iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
		}

		double tmp_dist = 0.0;
		size_t num_non_self_loops = 0;
		for (size_t i = 0; i < num_neighbors; ++i) {
			if (neighbors[i] != seed) {
				tmp_dist += dist_scratch[i];
				++num_non_self_loops;
			}
		}
		assert((num_non_self_loops == size_constraint) ||
		       (num_non_self_loops == size_constraint - 1));
		assert(num_non_self_loops > 0);
		++sampled;
		sum_dist += tmp_dist / ((double) num_non_self_loops);
	}

	free(dist_scratch);

	*out_avg_seed_dist = sum_dist / ((double) sampled);

	return iscc_no_error();
}


scc_ErrorCode iscc_make_nng_clusters_from_seeds(scc_Clustering* const clustering,
                                                void* const data_set,
                                                const iscc_SeedResult* const seed_result,
                                                iscc_Digraph* const nng,
                                                const bool nng_is_ordered,
                                                scc_UnassignedMethod unassigned_method,
                                                const bool radius_constraint,
                                                const double radius,
                                                size_t len_primary_data_points,
                                                const scc_PointIndex primary_data_points[],
                                                scc_UnassignedMethod secondary_unassigned_method,
                                                const bool secondary_radius_constraint,
                                                const double secondary_radius)
{
	assert(iscc_check_input_clustering(clustering));
	assert(iscc_check_data_set(data_set));
	assert(iscc_num_data_points(data_set) == clustering->num_data_points);
	assert(seed_result->count > 0);
	assert(seed_result->seeds != NULL);
	assert(iscc_digraph_is_valid(nng));
	assert(!iscc_digraph_is_empty(nng));
	assert((unassigned_method == SCC_UM_IGNORE) ||
	       (unassigned_method == SCC_UM_ANY_NEIGHBOR) ||
	       (unassigned_method == SCC_UM_CLOSEST_ASSIGNED) ||
	       (unassigned_method == SCC_UM_CLOSEST_SEED));
	assert(!radius_constraint || (radius > 0.0));
	assert((secondary_unassigned_method == SCC_UM_IGNORE) ||
	       (secondary_unassigned_method == SCC_UM_CLOSEST_ASSIGNED) ||
	       (secondary_unassigned_method == SCC_UM_CLOSEST_SEED));
	assert(!secondary_radius_constraint || (secondary_radius > 0.0));

	// Assign seeds and their neighbors
	const size_t num_assigned_as_seed_or_neighbor = iscc_assign_seeds_and_neighbors(clustering, seed_result, nng);
	size_t total_assigned = num_assigned_as_seed_or_neighbor;

	// Are we done?
	if ((total_assigned == clustering->num_data_points) ||
	        ((unassigned_method == SCC_UM_IGNORE) && (secondary_unassigned_method == SCC_UM_IGNORE))) {
		return iscc_no_error();
	}

	// If SCC_UM_CLOSEST_ASSIGNED, construct array with seeds and their neighbors for the nn search object
	scc_PointIndex* seed_or_neighbor = NULL;
	if ((unassigned_method == SCC_UM_CLOSEST_ASSIGNED) ||
	        (secondary_unassigned_method == SCC_UM_CLOSEST_ASSIGNED)) {
		seed_or_neighbor = malloc(sizeof(scc_PointIndex[num_assigned_as_seed_or_neighbor]));
		if (seed_or_neighbor == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);

		scc_PointIndex* write_seed_or_neighbor = seed_or_neighbor;
		assert(clustering->num_data_points <= ISCC_POINTINDEX_MAX);
		const scc_PointIndex num_data_points_pi = (scc_PointIndex) clustering->num_data_points; // If `scc_PointIndex` is signed.
		for (scc_PointIndex i = 0; i < num_data_points_pi; ++i) {
			if (clustering->cluster_label[i] != SCC_CLABEL_NA) {
				*write_seed_or_neighbor = i;
				++write_seed_or_neighbor;
			}
		}
		assert(((size_t) (write_seed_or_neighbor - seed_or_neighbor)) == num_assigned_as_seed_or_neighbor);
	}

	// Run assignment by nng. When nng is ordered, we can use it for `SCC_UM_CLOSEST_ASSIGNED` as well.
	// (NNG already contains radius constraint.)
	if ((unassigned_method == SCC_UM_ANY_NEIGHBOR) ||
	        (nng_is_ordered && (unassigned_method == SCC_UM_CLOSEST_ASSIGNED))) {
		total_assigned += iscc_assign_by_nng(clustering, nng);

		// Ignore remaining points if SCC_UM_ANY_NEIGHBOR
		if (unassigned_method == SCC_UM_ANY_NEIGHBOR) {
			unassigned_method = SCC_UM_IGNORE;
		}

		// Are we done?
		if ((total_assigned == clustering->num_data_points) ||
		        ((unassigned_method == SCC_UM_IGNORE) && (secondary_unassigned_method == SCC_UM_IGNORE))) {
			free(seed_or_neighbor);
			return iscc_no_error();
		}
	}

	// No need for nng any more
	iscc_free_digraph(nng);

	scc_ErrorCode ec = SCC_ER_OK;
	iscc_NNSearchObject* nn_assigned_search_object = NULL;
	iscc_NNSearchObject* nn_seed_search_object = NULL;

	if ((unassigned_method == SCC_UM_CLOSEST_ASSIGNED) ||
	        (secondary_unassigned_method == SCC_UM_CLOSEST_ASSIGNED)) {
		assert(seed_or_neighbor != NULL);
		if (!iscc_init_nn_search_object(data_set,
		                                num_assigned_as_seed_or_neighbor,
		                                seed_or_neighbor,
		                                &nn_assigned_search_object)) {
			ec = iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
		}
	}

	if (ec != SCC_ER_OK) {
		free(seed_or_neighbor);
		return ec;
	}

	if ((unassigned_method == SCC_UM_CLOSEST_SEED) ||
	        (secondary_unassigned_method == SCC_UM_CLOSEST_SEED)) {
		if (!iscc_init_nn_search_object(data_set,
		                                seed_result->count,
		                                seed_result->seeds,
		                                &nn_seed_search_object)) {
			ec = iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
		}
	}

	if (ec != SCC_ER_OK) {
		free(seed_or_neighbor);
		if (nn_assigned_search_object != NULL) {
			iscc_close_nn_search_object(&nn_assigned_search_object);
		}
		return ec;
	}

	size_t num_to_assign = 0;
	scc_PointIndex* const to_assign = malloc(sizeof(scc_PointIndex[clustering->num_data_points - total_assigned + 1]));
	if (to_assign == NULL) {
		free(seed_or_neighbor);
		if (nn_assigned_search_object != NULL) {
			iscc_close_nn_search_object(&nn_assigned_search_object);
		}
		if (nn_seed_search_object != NULL) {
			iscc_close_nn_search_object(&nn_seed_search_object);
		}
		return iscc_make_error(SCC_ER_NO_MEMORY);
	}

	if (primary_data_points != NULL) {
		for (size_t i = 0; i < len_primary_data_points; ++i) {
			to_assign[num_to_assign] = primary_data_points[i];
			num_to_assign += (clustering->cluster_label[primary_data_points[i]] == SCC_CLABEL_NA);
		}
	} else {
		const scc_PointIndex num_data_points_pi = (scc_PointIndex) clustering->num_data_points;
		for (scc_PointIndex i = 0; i < num_data_points_pi; ++i) {
			to_assign[num_to_assign] = i;
			num_to_assign += (clustering->cluster_label[i] == SCC_CLABEL_NA);
		}
	}

	if (num_to_assign > 0) {
		if (unassigned_method == SCC_UM_CLOSEST_ASSIGNED) {
			ec = iscc_assign_by_nn_search(clustering,
			                              nn_assigned_search_object,
			                              num_to_assign,
			                              to_assign,
			                              radius_constraint,
			                              radius);
		} else if (unassigned_method == SCC_UM_CLOSEST_SEED) {
			ec = iscc_assign_by_nn_search(clustering,
			                              nn_seed_search_object,
			                              num_to_assign,
			                              to_assign,
			                              radius_constraint,
			                              radius);
		}
	}

	if (ec != SCC_ER_OK) {
		free(seed_or_neighbor);
		free(to_assign);
		if (nn_assigned_search_object != NULL) {
			iscc_close_nn_search_object(&nn_assigned_search_object);
		}
		if (nn_seed_search_object != NULL) {
			iscc_close_nn_search_object(&nn_seed_search_object);
		}
		return ec;
	}

	if (secondary_unassigned_method != SCC_UM_IGNORE) {
		size_t num_to_assign = 0;
		const scc_PointIndex num_data_points_pi = (scc_PointIndex) clustering->num_data_points;
		for (scc_PointIndex i = 0; i < num_data_points_pi; ++i) {
			to_assign[num_to_assign] = i;
			num_to_assign += (clustering->cluster_label[i] == SCC_CLABEL_NA);
		}

		if (num_to_assign > 0) {
			if (secondary_unassigned_method == SCC_UM_CLOSEST_ASSIGNED) {
				ec = iscc_assign_by_nn_search(clustering,
				                              nn_assigned_search_object,
				                              num_to_assign,
				                              to_assign,
				                              secondary_radius_constraint,
				                              secondary_radius);
			} else if (secondary_unassigned_method == SCC_UM_CLOSEST_SEED) {
				ec = iscc_assign_by_nn_search(clustering,
				                              nn_seed_search_object,
				                              num_to_assign,
				                              to_assign,
				                              secondary_radius_constraint,
				                              secondary_radius);
			}
		}
	}

	free(seed_or_neighbor);
	free(to_assign);
	if (nn_assigned_search_object != NULL) {
		iscc_close_nn_search_object(&nn_assigned_search_object);
	}
	if (nn_seed_search_object != NULL) {
		iscc_close_nn_search_object(&nn_seed_search_object);
	}

	return iscc_no_error();
}


// =============================================================================
// Static function implementations
// =============================================================================

static scc_ErrorCode iscc_make_nng(void* const data_set,
                                   const size_t num_data_points,
                                   const size_t len_search_indices,
                                   const scc_PointIndex search_indices[const],
                                   const size_t len_query_indices,
                                   const scc_PointIndex query_indices[const],
                                   const uint32_t k,
                                   const bool radius_search,
                                   const double radius,
                                   size_t* const out_len_query_indices,
                                   scc_PointIndex out_query_indices[const],
                                   iscc_Digraph* const out_nng)
{
	assert(iscc_check_data_set(data_set));
	assert(len_search_indices > 0);
	assert(len_query_indices > 0);
	assert(k > 0);
	assert(len_search_indices >= k);
	assert(!radius_search || (radius > 0.0));
	assert(out_nng != NULL);

	iscc_NNSearchObject* nn_search_object;
	if (!iscc_init_nn_search_object(data_set,
	                                len_search_indices,
	                                search_indices,
	                                &nn_search_object)) {
		return iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
	}

	scc_ErrorCode ec;
	if ((ec = iscc_make_nng_from_search_object(nn_search_object,
	                                           num_data_points,
	                                           len_query_indices,
	                                           query_indices,
	                                           k,
	                                           radius_search,
	                                           radius,
	                                           out_len_query_indices,
	                                           out_query_indices,
	                                           out_nng)) != SCC_ER_OK) {
		iscc_close_nn_search_object(&nn_search_object);
		return ec;
	}

	if (!iscc_close_nn_search_object(&nn_search_object)) {
		iscc_free_digraph(out_nng);
		return iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
	}

	return iscc_no_error();
}


static scc_ErrorCode iscc_make_nng_from_search_object(iscc_NNSearchObject* const nn_search_object,
                                                      const size_t num_data_points,
                                                      const size_t len_query_indices,
                                                      const scc_PointIndex query_indices[const],
                                                      const uint32_t k,
                                                      const bool radius_search,
                                                      const double radius,
                                                      size_t* const out_len_query_indices,
                                                      scc_PointIndex out_query_indices[const],
                                                      iscc_Digraph* const out_nng)
{
	assert(nn_search_object != NULL);
	assert(len_query_indices > 0);
	assert(k > 0);
	assert(!radius_search || (radius > 0.0));
	assert(out_nng != NULL);

	scc_PointIndex* internal_out_query_indices = NULL;
	scc_PointIndex* dist_out_query_indices = NULL;

	if (radius_search) {
		if (out_query_indices != NULL) {
			dist_out_query_indices = out_query_indices;
		} else {
			internal_out_query_indices = malloc(sizeof(scc_PointIndex[len_query_indices]));
			if (internal_out_query_indices == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);
			dist_out_query_indices = internal_out_query_indices;
		}
	}

	scc_ErrorCode ec;
	if ((ec = iscc_init_digraph(num_data_points,
	                            len_query_indices * k,
	                            out_nng)) != SCC_ER_OK) {
		free(internal_out_query_indices);
		return ec;
	}

	size_t num_ok_queries = 0;
	if (!iscc_nearest_neighbor_search(nn_search_object,
	                                  len_query_indices,
	                                  query_indices,
	                                  k,
	                                  radius_search,
	                                  radius,
	                                  &num_ok_queries,
	                                  dist_out_query_indices,
	                                  out_nng->head)) {
		free(internal_out_query_indices);
		iscc_free_digraph(out_nng);
		return iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
	}

	iscc_ArcIndex* write_tail_ptr = out_nng->tail_ptr;
	*write_tail_ptr = 0;
	++write_tail_ptr;
	const iscc_ArcIndex* const write_tail_ptr_stop = write_tail_ptr + num_data_points;

	if (radius_search || query_indices != NULL) {
		const scc_PointIndex* ok_q;
		if (radius_search) {
			assert(dist_out_query_indices != NULL);
			ok_q = dist_out_query_indices;
		} else {
			assert(len_query_indices == num_ok_queries);
			assert(query_indices != NULL);
			ok_q = query_indices;
		}

		scc_PointIndex i = 0;
		const scc_PointIndex* const ok_q_stop = ok_q + num_ok_queries;
		for (; ok_q < ok_q_stop; ++ok_q) {
			for (; i < *ok_q; ++i) {
				*write_tail_ptr = *(write_tail_ptr - 1);
				++write_tail_ptr;
			}
			*write_tail_ptr = *(write_tail_ptr - 1) + k;
			++write_tail_ptr;
			++i;
		}
	} else {
		assert(!radius_search && query_indices == NULL);
		assert(len_query_indices == num_ok_queries);
		for (size_t q = 0; q < len_query_indices; ++q) {
			*write_tail_ptr = *(write_tail_ptr - 1) + k;
			++write_tail_ptr;
		}
	}

	for (; write_tail_ptr < write_tail_ptr_stop; ++write_tail_ptr) {
		*write_tail_ptr = *(write_tail_ptr - 1);
	}

	if (internal_out_query_indices != NULL) {
		assert(radius_search);
		assert(out_query_indices == NULL);
		free(internal_out_query_indices);
	}

	if (len_query_indices > num_ok_queries) {
		assert(radius_search);
		if ((ec = iscc_change_arc_storage(out_nng, num_ok_queries * k)) != SCC_ER_OK) {
			iscc_free_digraph(out_nng);
			return ec;
		}
	}

	if (out_len_query_indices != NULL) {
		*out_len_query_indices = num_ok_queries;
	}

	return iscc_no_error();
}


static inline void iscc_ensure_self_match(iscc_Digraph* const nng,
                                          const size_t len_search_indices,
                                          const scc_PointIndex search_indices[const])
{
	assert(iscc_digraph_is_valid(nng));
	assert(!iscc_digraph_is_empty(nng));
	assert(len_search_indices > 0);

	/* When there's identical data points, `iscc_nearest_neighbor_search_digraph` may not
	 * return a self-loop when a query is a search point. The NNG clustering functions
	 * require this. However, if all data points are unique, or the query and search sets
	 * are disjoint, it's safe to call `iscc_make_nng` without `iscc_ensure_self_match`. */

	if (search_indices == NULL) {
		assert(len_search_indices <= ISCC_POINTINDEX_MAX);
		const scc_PointIndex len_search_indices_pi = (scc_PointIndex) len_search_indices; // If `scc_PointIndex` is signed.
		for (scc_PointIndex search_point = 0; search_point < len_search_indices_pi; ++search_point) {
			scc_PointIndex* v_arc = nng->head + nng->tail_ptr[search_point];
			const scc_PointIndex* const v_arc_stop = nng->head + nng->tail_ptr[search_point + 1];
			if ((*v_arc == search_point) || (v_arc == v_arc_stop)) continue;
			for (++v_arc; (*v_arc != search_point) && (v_arc != v_arc_stop); ++v_arc);
			if (v_arc == v_arc_stop) *(v_arc - 1) = search_point;
		}

	} else if (search_indices != NULL) {
		for (size_t s = 0; s < len_search_indices; ++s) {
			const scc_PointIndex search_point = search_indices[s];
			scc_PointIndex* v_arc = nng->head + nng->tail_ptr[search_point];
			const scc_PointIndex* const v_arc_stop = nng->head + nng->tail_ptr[search_point + 1];
			if ((*v_arc == search_point) || (v_arc == v_arc_stop)) continue;
			for (++v_arc; (*v_arc != search_point) && (v_arc != v_arc_stop); ++v_arc);
			if (v_arc == v_arc_stop) *(v_arc - 1) = search_point;
		}
	}
}


static scc_ErrorCode iscc_type_count(const size_t num_data_points,
                                     const uint32_t size_constraint,
                                     const uint_fast16_t num_types,
                                     const uint32_t type_constraints[const static num_types],
                                     const scc_TypeLabel type_labels[const static num_data_points],
                                     iscc_TypeCount* const out_type_result)
{
	assert(num_data_points > 1);
	assert(size_constraint >= 2);
	assert(num_types >= 2);
	assert(num_types <= ISCC_TYPELABEL_MAX);
	assert(type_constraints != NULL);
	assert(type_labels != NULL);
	assert(out_type_result != NULL);

	*out_type_result = (iscc_TypeCount) {
		.sum_type_constraints = 0,
		.type_group_size = calloc(num_types, sizeof(size_t)),
		.point_store = malloc(sizeof(scc_PointIndex[num_data_points])),
		.type_groups = malloc(sizeof(scc_PointIndex*[num_types])),
	};

	if ((out_type_result->type_group_size == NULL) || (out_type_result->point_store == NULL) || (out_type_result->type_groups == NULL)) {
		free(out_type_result->type_group_size);
		free(out_type_result->point_store);
		free(out_type_result->type_groups);
		return iscc_make_error(SCC_ER_NO_MEMORY);
	}

	for (size_t i = 0; i < num_data_points; ++i) {
		assert(type_labels[i] < (scc_TypeLabel) num_types);
		++out_type_result->type_group_size[type_labels[i]];
	}

	for (uint_fast16_t i = 0; i < num_types; ++i) {
		if (out_type_result->type_group_size[i] < type_constraints[i]) {
			free(out_type_result->type_group_size);
			free(out_type_result->point_store);
			free(out_type_result->type_groups);
			return iscc_make_error_msg(SCC_ER_NO_SOLUTION, "Fewer data points than type size constraint.");
		}
		out_type_result->sum_type_constraints += type_constraints[i];
	}

	if (out_type_result->sum_type_constraints > size_constraint) {
		free(out_type_result->type_group_size);
		free(out_type_result->point_store);
		free(out_type_result->type_groups);
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Type constraint cannot be larger than overall size constraint.");
	}

	out_type_result->type_groups[0] = out_type_result->point_store + out_type_result->type_group_size[0];
	for (uint_fast16_t i = 1; i < num_types; ++i) {
		out_type_result->type_groups[i] = out_type_result->type_groups[i - 1] + out_type_result->type_group_size[i];
	}

	assert(num_data_points <= ISCC_POINTINDEX_MAX);
	const scc_PointIndex num_data_points_pi = (scc_PointIndex) num_data_points; // if case `scc_PointIndex` is signed.
	for (scc_PointIndex i = 0; i < num_data_points_pi; ++i) {
		--(out_type_result->type_groups[type_labels[i]]);
		*(out_type_result->type_groups[type_labels[i]]) = i;
	}

	return iscc_no_error();
}


static size_t iscc_assign_seeds_and_neighbors(scc_Clustering* const clustering,
                                              const iscc_SeedResult* const seed_result,
                                              iscc_Digraph* const nng)
{
	assert(iscc_check_input_clustering(clustering));
	assert(clustering->cluster_label != NULL);
	assert(seed_result->count > 0);
	assert(seed_result->seeds != NULL);
	assert(iscc_digraph_is_valid(nng));
	assert(!iscc_digraph_is_empty(nng));

	clustering->num_clusters = seed_result->count;

	for (size_t i = 0; i < clustering->num_data_points; ++i) {
		clustering->cluster_label[i] = SCC_CLABEL_NA;
	}

	size_t num_assigned = 0;
	scc_Clabel clabel = 0;
	const scc_PointIndex* const seed_stop = seed_result->seeds + seed_result->count;
	for (const scc_PointIndex* seed = seed_result->seeds;
	        seed != seed_stop; ++seed, ++clabel) {
		assert(clabel != SCC_CLABEL_NA);
		assert(clabel < SCC_CLABEL_MAX);
		assert(clustering->cluster_label[*seed] == SCC_CLABEL_NA);

		const scc_PointIndex* const s_arc_stop = nng->head + nng->tail_ptr[*seed + 1];
		for (const scc_PointIndex* s_arc = nng->head + nng->tail_ptr[*seed];
		        s_arc != s_arc_stop; ++s_arc) {
			assert(clustering->cluster_label[*s_arc] == SCC_CLABEL_NA);
			clustering->cluster_label[*s_arc] = clabel;
		}
		num_assigned += (nng->tail_ptr[*seed + 1] - nng->tail_ptr[*seed]) + // Number of arcs from seed
		                    (clustering->cluster_label[*seed] == SCC_CLABEL_NA); // In the case of no seed self-loop
		clustering->cluster_label[*seed] = clabel; // Assign seed last so seed `assert` work also in case of self-loops
	}

	assert(clabel == (scc_Clabel) clustering->num_clusters);

	return num_assigned;
}


static size_t iscc_assign_by_nng(scc_Clustering* const clustering,
                                 iscc_Digraph* const nng)
{
	assert(iscc_check_input_clustering(clustering));
	assert(iscc_digraph_is_valid(nng));
	assert(!iscc_digraph_is_empty(nng));

	bool* const scratch = malloc(sizeof(bool[clustering->num_data_points]));
	if (scratch == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);
	for (size_t i = 0; i < clustering->num_data_points; ++i) {
		scratch[i] = (clustering->cluster_label[i] == SCC_CLABEL_NA);
	}

	size_t num_assigned_by_nng = 0;
	for (size_t i = 0; i < clustering->num_data_points; ++i) {
		if (scratch[i]) {
			assert(clustering->cluster_label[i] == SCC_CLABEL_NA);
			const scc_PointIndex* const v_arc_stop = nng->head + nng->tail_ptr[i + 1];
			for (const scc_PointIndex* v_arc = nng->head + nng->tail_ptr[i];
			        v_arc != v_arc_stop; ++v_arc) {
				if (!scratch[*v_arc]) {
					assert(clustering->cluster_label[*v_arc] != SCC_CLABEL_NA);
					clustering->cluster_label[i] = clustering->cluster_label[*v_arc];
					++num_assigned_by_nng;
					break;
				}
			}
		}
	}

	free(scratch);

	return num_assigned_by_nng;
}


static scc_ErrorCode iscc_assign_by_nn_search(scc_Clustering* const clustering,
                                              iscc_NNSearchObject* const nn_search_object,
                                              const size_t num_to_assign,
                                              scc_PointIndex to_assign[restrict const static num_to_assign],
                                              const bool radius_constraint,
                                              const double radius)
{
	assert(iscc_check_input_clustering(clustering));
	assert(nn_search_object != NULL);
	assert(num_to_assign > 0);
	assert(to_assign != NULL);
	assert(!radius_constraint || (radius > 0.0));

	size_t num_ok_queries = 0;
	scc_PointIndex* out_ok_query = NULL;
	if (radius_constraint) {
		out_ok_query = to_assign;
	}
	scc_PointIndex* const out_nn_indices = malloc(sizeof(scc_PointIndex[num_to_assign]));

	if (!iscc_nearest_neighbor_search(nn_search_object,
	                                  num_to_assign,
	                                  to_assign,
	                                  1,
	                                  radius_constraint,
	                                  radius,
	                                  &num_ok_queries,
	                                  out_ok_query,
	                                  out_nn_indices)) {
		free(out_nn_indices);
		return iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
	}

	if (!radius_constraint) {
		assert(num_ok_queries == num_to_assign);
		out_ok_query = to_assign;
	}

	for (size_t i = 0; i < num_ok_queries; ++i) {
		assert(clustering->cluster_label[out_ok_query[i]] == SCC_CLABEL_NA);
		assert(clustering->cluster_label[out_nn_indices[i]] != SCC_CLABEL_NA);
		clustering->cluster_label[out_ok_query[i]] = clustering->cluster_label[out_nn_indices[i]];
	}

	free(out_nn_indices);

	return iscc_no_error();
}


#ifdef SCC_STABLE_NNG

static int iscc_compare_PointIndex(const void* const a, const void* const b)
{
    const scc_PointIndex arg1 = *(const scc_PointIndex* const)a;
    const scc_PointIndex arg2 = *(const scc_PointIndex* const)b;
    return (arg1 > arg2) - (arg1 < arg2);
}


static void iscc_sort_nng(iscc_Digraph* const nng)
{
	for (size_t v = 0; v < nng->vertices; ++v) {
		const size_t count = nng->tail_ptr[v + 1] - nng->tail_ptr[v];
		if (count > 1) {
			qsort(nng->head + nng->tail_ptr[v], count, sizeof(scc_PointIndex), iscc_compare_PointIndex);
		}
	}
}

#endif // ifdef SCC_STABLE_NNG
