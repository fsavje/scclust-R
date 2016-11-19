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

#include "nng_core.h"

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "../include/scclust.h"
#include "digraph_core.h"
#include "digraph_operations.h"
#include "dist_search.h"
#include "error.h"
#include "nng_findseeds.h"
#include "scclust_internal.h"


// =============================================================================
// Internal structs & variables
// =============================================================================

typedef struct iscc_TypeCount iscc_TypeCount;
struct iscc_TypeCount {
	uint32_t sum_type_constraints;
	size_t* type_group_size;
	iscc_Dpid* point_store;
	iscc_Dpid** type_groups;
};


static const size_t ISCC_ESTIMATE_AVG_MAX_SAMPLE = 1000;


// =============================================================================
// Internal function prototypes
// =============================================================================

static scc_ErrorCode iscc_make_nng(void* data_set,
                                   size_t len_search_indices,
                                   const iscc_Dpid search_indices[],
                                   size_t len_query_indicators,
                                   const bool query_indicators[],
                                   bool out_query_indicators[],
                                   uint32_t k,
                                   bool radius_search,
                                   double radius,
                                   bool accept_partial,
                                   uintmax_t max_arcs,
                                   iscc_Digraph* out_nng);

static scc_ErrorCode iscc_make_nng_from_search_object(iscc_NNSearchObject* nn_search_object,
                                                      size_t len_query_indicators,
                                                      const bool query_indicators[],
                                                      bool out_query_indicators[],
                                                      uint32_t k,
                                                      bool radius_search,
                                                      double radius,
                                                      bool accept_partial,
                                                      uintmax_t max_arcs,
                                                      iscc_Digraph* out_nng);

static inline void iscc_ensure_self_match(iscc_Digraph* nng,
                                          size_t len_search_indices,
                                          const iscc_Dpid search_indices[]);

static scc_ErrorCode iscc_type_count(size_t num_data_points,
                                     uint32_t size_constraint,
                                     uint_fast16_t num_types,
                                     const uint32_t type_size_constraints[static num_types],
                                     const scc_TypeLabel type_labels[static num_data_points],
                                     iscc_TypeCount* out_type_result);

static size_t iscc_assign_seeds_and_neighbors(scc_Clustering* clustering,
                                              const iscc_SeedResult* seed_result,
                                              iscc_Digraph* nng);

static size_t iscc_assign_by_nng(scc_Clustering* clustering,
                                 iscc_Digraph* nng,
                                 bool scratch[restrict static clustering->num_data_points]);

static scc_ErrorCode iscc_assign_by_nn_search(scc_Clustering* clustering,
                                              iscc_NNSearchObject* nn_search_object,
                                              size_t num_to_assign,
                                              const bool to_assign[restrict static clustering->num_data_points],
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
                                                const bool primary_data_points[const],
                                                const bool radius_constraint,
                                                const double radius,
                                                iscc_Digraph* const out_nng)
{
	assert(data_set != NULL);
	assert(num_data_points >= 2);
	assert(size_constraint <= num_data_points);
	assert(size_constraint >= 2);
	assert(!radius_constraint || (radius > 0.0));
	assert(out_nng != NULL);

	size_t num_queries;
	if (primary_data_points == NULL) {
		num_queries = num_data_points;
	} else {
		num_queries = 0;
		for (size_t i = 0; i < num_data_points; ++i) {
			num_queries += primary_data_points[i];
		}
		if (num_queries == 0) return iscc_make_error(SCC_ER_NO_CLUST_EXIST_CONSTRAINT);
	}

	scc_ErrorCode ec;
	if ((ec = iscc_make_nng(data_set,
	                        num_data_points,
	                        NULL,
	                        num_data_points,
	                        primary_data_points,
	                        NULL,
	                        size_constraint,
	                        radius_constraint,
	                        radius,
	                        false,
	                        (size_constraint * num_queries),
	                        out_nng)) != SCC_ER_OK) {
		return ec;
	}

	if (iscc_digraph_is_empty(out_nng)) {
		iscc_free_digraph(out_nng);
		return iscc_make_error(SCC_ER_NO_CLUST_EXIST_RADIUS);
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
                                                const uint32_t type_size_constraints[const static num_types],
                                                const scc_TypeLabel type_labels[const static num_data_points],
                                                const bool primary_data_points[const],
                                                const bool radius_constraint,
                                                const double radius,
                                                iscc_Digraph* const out_nng)
{
	assert(data_set != NULL);
	assert(num_data_points >= 2);
	assert(size_constraint <= num_data_points);
	assert(size_constraint >= 2);
	assert(num_types >= 2);
	assert(num_types <= ISCC_TYPELABEL_MAX);
	assert(type_size_constraints != NULL);
	assert(type_labels != NULL);
	assert(!radius_constraint || (radius > 0.0));
	assert(out_nng != NULL);

	size_t num_queries;
	if (primary_data_points == NULL) {
		num_queries = num_data_points;
	} else {
		num_queries = 0;
		for (size_t i = 0; i < num_data_points; ++i) {
			num_queries += primary_data_points[i];
		}
		if (num_queries == 0) return iscc_make_error(SCC_ER_NO_CLUST_EXIST_CONSTRAINT);
	}

	bool* seedable;
	const bool* seedable_const;
	if (radius_constraint) {
		seedable = malloc(sizeof(bool[num_data_points]));
		if (seedable == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);
		seedable_const = seedable;
		if (primary_data_points == NULL) {
			for (size_t i = 0; i < num_data_points; ++i) {
				seedable[i] = true;
			}
		} else {
			memcpy(seedable, primary_data_points, sizeof(bool[num_data_points]));
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
	                          type_size_constraints,
	                          type_labels,
	                          &tc)) != SCC_ER_OK) {
		free(seedable);
		free(nng_by_type);
		return ec;
	}

	uint_fast16_t num_non_zero_type_constraints = 0;
	for (uint_fast16_t i = 0; i < num_types; ++i) {
		if (type_size_constraints[i] > 0) {
			if ((ec = iscc_make_nng(data_set,
			                        tc.type_group_size[i],
			                        tc.type_groups[i],
			                        num_data_points,
			                        seedable_const,
			                        seedable,
			                        type_size_constraints[i],
			                        radius_constraint,
			                        radius,
			                        false,
			                        (type_size_constraints[i] * num_queries),
			                        &nng_by_type[num_non_zero_type_constraints])) != SCC_ER_OK) {
				break;
			}
			++num_non_zero_type_constraints;
			if (iscc_digraph_is_empty(&nng_by_type[num_non_zero_type_constraints - 1])) {
				ec = iscc_make_error(SCC_ER_NO_CLUST_EXIST_RADIUS);
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
			ec = iscc_digraph_union_and_delete(num_non_zero_type_constraints, nng_by_type, seedable_const, true, out_nng);
		} else {
			ec = iscc_digraph_union_and_delete(num_non_zero_type_constraints, nng_by_type, seedable_const, false, out_nng);
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
		                        NULL,
		                        num_data_points,
		                        seedable_const,
		                        seedable,
		                        size_constraint,
		                        radius_constraint,
		                        radius,
		                        false,
		                        (size_constraint * num_queries),
		                        &nng_sum[1])) != SCC_ER_OK) {
			free(seedable);
			iscc_free_digraph(&nng_sum[0]);
			return ec;
		}

		if (iscc_digraph_is_empty(&nng_sum[1])) {
			ec = iscc_make_error(SCC_ER_NO_CLUST_EXIST_RADIUS);
		} else {
			ec = iscc_digraph_difference(&nng_sum[1], &nng_sum[0], additional_nn_needed);
		}

		if (ec == SCC_ER_OK) {
			ec = iscc_digraph_union_and_delete(2, nng_sum, seedable_const, false, out_nng);
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
	assert(data_set != NULL);
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
		const iscc_Dpid seed = seed_result->seeds[s];
		const size_t num_neighbors = (nng->tail_ptr[seed + 1] - nng->tail_ptr[seed]);
		const iscc_Dpid* const neighbors = nng->head + nng->tail_ptr[seed];

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
                                                const bool primary_data_points[const],
                                                scc_UnassignedMethod secondary_unassigned_method,
                                                const bool secondary_radius_constraint,
                                                const double secondary_radius)
{
	assert(iscc_check_input_clustering(clustering));
	assert(data_set != NULL);
	assert(seed_result->count > 0);
	assert(seed_result->seeds != NULL);
	assert(iscc_digraph_is_valid(nng));
	assert(!iscc_digraph_is_empty(nng));
	assert((unassigned_method == SCC_UM_IGNORE) ||
	       (unassigned_method == SCC_UM_ASSIGN_BY_NNG) ||
	       (unassigned_method == SCC_UM_CLOSEST_ASSIGNED) ||
	       (unassigned_method == SCC_UM_CLOSEST_SEED));
	assert(!radius_constraint || (radius > 0.0));
	assert((primary_data_points != NULL) || (secondary_unassigned_method == SCC_UM_IGNORE));
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
	iscc_Dpid* seed_or_neighbor = NULL;
	if ((unassigned_method == SCC_UM_CLOSEST_ASSIGNED) ||
	        (secondary_unassigned_method == SCC_UM_CLOSEST_ASSIGNED)) {
		seed_or_neighbor = malloc(sizeof(iscc_Dpid[num_assigned_as_seed_or_neighbor]));
		if (seed_or_neighbor == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);

		iscc_Dpid* write_seed_or_neighbor = seed_or_neighbor;
		assert(clustering->num_data_points <= ISCC_DPID_MAX);
		const iscc_Dpid num_data_points_dpid = (iscc_Dpid) clustering->num_data_points; // If `iscc_Dpid` is signed.
		for (iscc_Dpid i = 0; i < num_data_points_dpid; ++i) {
			if (clustering->cluster_label[i] != SCC_CLABEL_NA) {
				*write_seed_or_neighbor = i;
				++write_seed_or_neighbor;
			}
		}
		assert(((size_t) (write_seed_or_neighbor - seed_or_neighbor)) == num_assigned_as_seed_or_neighbor);
	}

	// Indicators of main data points that are unassigned (first used as scratch)
	bool* main_assign = malloc(sizeof(bool[clustering->num_data_points]));
	if (main_assign == NULL) {
		free(seed_or_neighbor);
		return iscc_make_error(SCC_ER_NO_MEMORY);
	}

	// Run assignment by nng. When nng is ordered, we can use it for `SCC_UM_CLOSEST_ASSIGNED` as well.
	// (NNG already contains radius constraint.)
	if ((unassigned_method == SCC_UM_ASSIGN_BY_NNG) ||
	        (nng_is_ordered && (unassigned_method == SCC_UM_CLOSEST_ASSIGNED))) {
		total_assigned += iscc_assign_by_nng(clustering, nng, main_assign); // Use `main_assign` as scratch

		// Ignore remaining points if SCC_UM_ASSIGN_BY_NNG
		if (unassigned_method == SCC_UM_ASSIGN_BY_NNG) {
			unassigned_method = SCC_UM_IGNORE;
		}

		// Are we done?
		if ((total_assigned == clustering->num_data_points) ||
		        ((unassigned_method == SCC_UM_IGNORE) && (secondary_unassigned_method == SCC_UM_IGNORE))) {
			free(seed_or_neighbor);
			free(main_assign);
			return iscc_no_error();
		}
	}

	// No need for nng any more
	iscc_free_digraph(nng);

	// Derive which data points to assign
	bool* secondary_assign = NULL;
	size_t num_main_assign = 0;
	size_t num_secondary_assign = 0;
	if (primary_data_points == NULL) {
		// All data points are in main
		assert(secondary_unassigned_method == SCC_UM_IGNORE);
		for (size_t i = 0; i < clustering->num_data_points; ++i) {
			main_assign[i] = (clustering->cluster_label[i] == SCC_CLABEL_NA);
		}
		num_main_assign = clustering->num_data_points - total_assigned;

		#ifndef NDEBUG
			size_t dbg_main_count = 0;
			for (size_t i = 0; i < clustering->num_data_points; ++i) {
				dbg_main_count += main_assign[i];
			}
			assert(dbg_main_count == num_main_assign);
			assert(total_assigned + dbg_main_count == clustering->num_data_points);
		#endif

	} else if (secondary_unassigned_method == SCC_UM_IGNORE) {
		// Assign only data points in main
		for (size_t i = 0; i < clustering->num_data_points; ++i) {
			main_assign[i] = primary_data_points[i] && (clustering->cluster_label[i] == SCC_CLABEL_NA);
			num_main_assign += main_assign[i];
		}

		#ifndef NDEBUG
			size_t dbg_secondary_count = 0;
			for (size_t i = 0; i < clustering->num_data_points; ++i) {
				dbg_secondary_count += !primary_data_points[i] && (clustering->cluster_label[i] == SCC_CLABEL_NA);
			}
			assert(total_assigned + num_main_assign + dbg_secondary_count == clustering->num_data_points);
		#endif

	} else {
		// Assign both main and secondary
		secondary_assign = malloc(sizeof(bool[clustering->num_data_points]));
		if (secondary_assign == NULL) {
			free(seed_or_neighbor);
			free(main_assign);
			return iscc_make_error(SCC_ER_NO_MEMORY);
		}
		for (size_t i = 0; i < clustering->num_data_points; ++i) {
			main_assign[i] = primary_data_points[i] && (clustering->cluster_label[i] == SCC_CLABEL_NA);
			num_main_assign += main_assign[i];
			secondary_assign[i] = !primary_data_points[i] && (clustering->cluster_label[i] == SCC_CLABEL_NA);
			num_secondary_assign += secondary_assign[i];
		}
		assert(total_assigned + num_main_assign + num_secondary_assign == clustering->num_data_points);
	}

	if ((num_main_assign == 0) || (unassigned_method == SCC_UM_IGNORE)) {
		assert(main_assign != NULL);
		free(main_assign);
		main_assign = NULL;
		unassigned_method = SCC_UM_IGNORE;
	}

	if ((num_secondary_assign == 0) && (secondary_unassigned_method != SCC_UM_IGNORE)) {
		// If `secondary_unassigned_method == SCC_UM_IGNORE`, then `secondary_assign` should already be NULL
		assert(secondary_assign != NULL);
		free(secondary_assign);
		secondary_assign = NULL;
		secondary_unassigned_method = SCC_UM_IGNORE;
	}

	assert(((unassigned_method == SCC_UM_IGNORE) && (main_assign == NULL)) ||
	       ((unassigned_method != SCC_UM_IGNORE) && (main_assign != NULL)));
	assert(((secondary_unassigned_method == SCC_UM_IGNORE) && (secondary_assign == NULL)) ||
	       ((secondary_unassigned_method != SCC_UM_IGNORE) && (secondary_assign != NULL)));

	scc_ErrorCode ec = SCC_ER_OK;
	// Run assign to closest assigned
	if ((unassigned_method == SCC_UM_CLOSEST_ASSIGNED) ||
	        (secondary_unassigned_method == SCC_UM_CLOSEST_ASSIGNED)) {
		assert(seed_or_neighbor != NULL);
		iscc_NNSearchObject* nn_search_object;
		if (!iscc_init_nn_search_object(data_set,
		                                num_assigned_as_seed_or_neighbor,
		                                seed_or_neighbor,
		                                &nn_search_object)) {
			ec = iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
		}

		if ((ec == SCC_ER_OK) && (unassigned_method == SCC_UM_CLOSEST_ASSIGNED)) {
			assert(num_main_assign > 0);
			assert(main_assign != NULL);
			ec = iscc_assign_by_nn_search(clustering,
			                              nn_search_object,
			                              num_main_assign,
			                              main_assign,
			                              radius_constraint,
			                              radius);
			free(main_assign);
			main_assign = NULL;
		}

		if ((ec == SCC_ER_OK) && (secondary_unassigned_method == SCC_UM_CLOSEST_ASSIGNED)) {
			assert(num_secondary_assign > 0);
			assert(secondary_assign != NULL);
			ec = iscc_assign_by_nn_search(clustering,
			                              nn_search_object,
			                              num_secondary_assign,
			                              secondary_assign,
			                              secondary_radius_constraint,
			                              secondary_radius);
			free(secondary_assign);
			secondary_assign = NULL;
		}

		if (!iscc_close_nn_search_object(&nn_search_object)) {
			if (ec == SCC_ER_OK) ec = iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
		}

		if (ec != SCC_ER_OK) {
			free(seed_or_neighbor);
			free(main_assign);
			free(secondary_assign);
			return ec;
		}
	}

	free(seed_or_neighbor);
	seed_or_neighbor = NULL;

	// Run assign to closest seed
	if ((unassigned_method == SCC_UM_CLOSEST_SEED) ||
	        (secondary_unassigned_method == SCC_UM_CLOSEST_SEED)) {
		iscc_NNSearchObject* nn_search_object;
		if (!iscc_init_nn_search_object(data_set,
		                                seed_result->count,
		                                seed_result->seeds,
		                                &nn_search_object)) {
			ec = iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
		}

		if ((ec == SCC_ER_OK) && (unassigned_method == SCC_UM_CLOSEST_SEED)) {
			assert(num_main_assign > 0);
			assert(main_assign != NULL);
			ec = iscc_assign_by_nn_search(clustering,
			                              nn_search_object,
			                              num_main_assign,
			                              main_assign,
			                              radius_constraint,
			                              radius);
			free(main_assign);
			main_assign = NULL;
		}

		if ((ec == SCC_ER_OK) && (secondary_unassigned_method == SCC_UM_CLOSEST_SEED)) {
			assert(num_secondary_assign > 0);
			assert(secondary_assign != NULL);
			ec = iscc_assign_by_nn_search(clustering,
			                              nn_search_object,
			                              num_secondary_assign,
			                              secondary_assign,
			                              secondary_radius_constraint,
			                              secondary_radius);
			free(secondary_assign);
			secondary_assign = NULL;
		}

		if (!iscc_close_nn_search_object(&nn_search_object)) {
			if (ec == SCC_ER_OK) ec = iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
		}

		if (ec != SCC_ER_OK) {
			free(main_assign);
			free(secondary_assign);
			return ec;
		}
	}

	assert(main_assign == NULL);
	assert(secondary_assign == NULL);

	return iscc_no_error();
}


// =============================================================================
// Internal function implementations
// =============================================================================

static scc_ErrorCode iscc_make_nng(void* const data_set,
                                   const size_t len_search_indices,
                                   const iscc_Dpid search_indices[const],
                                   const size_t len_query_indicators,
                                   const bool query_indicators[const],
                                   bool out_query_indicators[const],
                                   const uint32_t k,
                                   const bool radius_search,
                                   const double radius,
                                   const bool accept_partial,
                                   const uintmax_t max_arcs,
                                   iscc_Digraph* const out_nng)
{
	assert(data_set != NULL);
	assert(len_search_indices > 0);
	assert(len_query_indicators > 0);
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
                                               len_query_indicators,
                                               query_indicators,
                                               out_query_indicators,
                                               k,
                                               radius_search,
                                               radius,
                                               accept_partial,
                                               max_arcs,
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
                                                      const size_t len_query_indicators,
                                                      const bool query_indicators[const],
                                                      bool out_query_indicators[const],
                                                      const uint32_t k,
                                                      const bool radius_search,
                                                      const double radius,
                                                      const bool accept_partial,
                                                      const uintmax_t max_arcs,
                                                      iscc_Digraph* const out_nng)
{
	assert(nn_search_object != NULL);
	assert(len_query_indicators > 0);
	assert(k > 0);
	assert(!radius_search || (radius > 0.0));
	assert(max_arcs >= k);
	assert(out_nng != NULL);

	scc_ErrorCode ec;
	if ((ec = iscc_init_digraph(len_query_indicators, max_arcs, out_nng)) != SCC_ER_OK) return ec;

	if (!iscc_nearest_neighbor_search_digraph(nn_search_object,
	                                          len_query_indicators,
	                                          query_indicators,
	                                          out_query_indicators,
	                                          k,
	                                          radius_search,
	                                          radius,
	                                          accept_partial,
	                                          out_nng->tail_ptr,
	                                          out_nng->head)) {
		iscc_free_digraph(out_nng);
		return iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
	}

	if (max_arcs > out_nng->tail_ptr[out_nng->vertices]) {
		if ((ec = iscc_change_arc_storage(out_nng, out_nng->tail_ptr[out_nng->vertices])) != SCC_ER_OK) {
			iscc_free_digraph(out_nng);
			return ec;
		}
	}

	return iscc_no_error();
}


static inline void iscc_ensure_self_match(iscc_Digraph* const nng,
                                          const size_t len_search_indices,
                                          const iscc_Dpid search_indices[const])
{
	assert(iscc_digraph_is_valid(nng));
	assert(!iscc_digraph_is_empty(nng));
	assert(len_search_indices > 0);

	/* When there's identical data points, `iscc_nearest_neighbor_search_digraph` may not
	 * return a self-loop when a query is a search point. The NNG clustering functions
	 * require this. However, if all data points are unique, or the query and search sets
	 * are disjoint, it's safe to call `iscc_make_nng` without `iscc_ensure_self_match`. */

	if (search_indices == NULL) {
		assert(len_search_indices <= ISCC_DPID_MAX);
		const iscc_Dpid len_search_indices_dpid = (iscc_Dpid) len_search_indices; // If `iscc_Dpid` is signed.
		for (iscc_Dpid search_point = 0; search_point < len_search_indices_dpid; ++search_point) {
			iscc_Dpid* v_arc = nng->head + nng->tail_ptr[search_point];
			const iscc_Dpid* const v_arc_stop = nng->head + nng->tail_ptr[search_point + 1];
			if ((*v_arc == search_point) || (v_arc == v_arc_stop)) continue;
			for (++v_arc; (*v_arc != search_point) && (v_arc != v_arc_stop); ++v_arc);
			if (v_arc == v_arc_stop) *(v_arc - 1) = search_point;
		}

	} else if (search_indices != NULL) {
		for (size_t s = 0; s < len_search_indices; ++s) {
			const iscc_Dpid search_point = search_indices[s];
			iscc_Dpid* v_arc = nng->head + nng->tail_ptr[search_point];
			const iscc_Dpid* const v_arc_stop = nng->head + nng->tail_ptr[search_point + 1];
			if ((*v_arc == search_point) || (v_arc == v_arc_stop)) continue;
			for (++v_arc; (*v_arc != search_point) && (v_arc != v_arc_stop); ++v_arc);
			if (v_arc == v_arc_stop) *(v_arc - 1) = search_point;
		}
	}
}


static scc_ErrorCode iscc_type_count(const size_t num_data_points,
                                     const uint32_t size_constraint,
                                     const uint_fast16_t num_types,
                                     const uint32_t type_size_constraints[const static num_types],
                                     const scc_TypeLabel type_labels[const static num_data_points],
                                     iscc_TypeCount* const out_type_result)
{
	assert(num_data_points > 1);
	assert(size_constraint >= 2);
	assert(num_types >= 2);
	assert(num_types <= ISCC_TYPELABEL_MAX);
	assert(type_size_constraints != NULL);
	assert(type_labels != NULL);
	assert(out_type_result != NULL);

	*out_type_result = (iscc_TypeCount) {
		.sum_type_constraints = 0,
		.type_group_size = calloc(num_types, sizeof(size_t)),
		.point_store = malloc(sizeof(iscc_Dpid[num_data_points])),
		.type_groups = malloc(sizeof(iscc_Dpid*[num_types])),
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
		if (out_type_result->type_group_size[i] < type_size_constraints[i]) {
			free(out_type_result->type_group_size);
			free(out_type_result->point_store);
			free(out_type_result->type_groups);
			return iscc_make_error(SCC_ER_NO_CLUST_EXIST_CONSTRAINT);
		}
		out_type_result->sum_type_constraints += type_size_constraints[i];
	}

	if (out_type_result->sum_type_constraints > size_constraint) {
		free(out_type_result->type_group_size);
		free(out_type_result->point_store);
		free(out_type_result->type_groups);
		return iscc_make_error(SCC_ER_INVALID_INPUT);
	}

	out_type_result->type_groups[0] = out_type_result->point_store + out_type_result->type_group_size[0];
	for (uint_fast16_t i = 1; i < num_types; ++i) {
		out_type_result->type_groups[i] = out_type_result->type_groups[i - 1] + out_type_result->type_group_size[i];
	}

	assert(num_data_points <= ISCC_DPID_MAX);
	const iscc_Dpid num_data_points_dpid = (iscc_Dpid) num_data_points; // if case `iscc_Dpid` is signed.
	for (iscc_Dpid i = 0; i < num_data_points_dpid; ++i) {
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
	const iscc_Dpid* const seed_stop = seed_result->seeds + seed_result->count;
	for (const iscc_Dpid* seed = seed_result->seeds;
	        seed != seed_stop; ++seed, ++clabel) {
		assert(clabel != SCC_CLABEL_NA);
		assert(clabel < SCC_CLABEL_MAX);
		assert(clustering->cluster_label[*seed] == SCC_CLABEL_NA);

		const iscc_Dpid* const s_arc_stop = nng->head + nng->tail_ptr[*seed + 1];
		for (const iscc_Dpid* s_arc = nng->head + nng->tail_ptr[*seed];
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
                                 iscc_Digraph* const nng,
                                 bool scratch[restrict const static clustering->num_data_points])
{
	assert(iscc_check_input_clustering(clustering));
	assert(iscc_digraph_is_valid(nng));
	assert(!iscc_digraph_is_empty(nng));
	assert(scratch != NULL);

	for (size_t i = 0; i < clustering->num_data_points; ++i) {
		scratch[i] = (clustering->cluster_label[i] == SCC_CLABEL_NA);
	}

	size_t num_assigned_by_nng = 0;
	for (size_t i = 0; i < clustering->num_data_points; ++i) {
		if (scratch[i]) {
			assert(clustering->cluster_label[i] == SCC_CLABEL_NA);
			const iscc_Dpid* const v_arc_stop = nng->head + nng->tail_ptr[i + 1];
			for (const iscc_Dpid* v_arc = nng->head + nng->tail_ptr[i];
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

	return num_assigned_by_nng;
}


static scc_ErrorCode iscc_assign_by_nn_search(scc_Clustering* const clustering,
                                              iscc_NNSearchObject* const nn_search_object,
                                              const size_t num_to_assign,
                                              const bool to_assign[restrict const static clustering->num_data_points],
                                              const bool radius_constraint,
                                              const double radius)
{
	assert(iscc_check_input_clustering(clustering));
	assert(nn_search_object != NULL);
	assert(num_to_assign > 0);
	assert(to_assign != NULL);
	assert(!radius_constraint || (radius > 0.0));

	scc_ErrorCode ec;
	iscc_Digraph priority_graph;
	if ((ec = iscc_make_nng_from_search_object(nn_search_object,
	                                           clustering->num_data_points,
	                                           to_assign,
	                                           NULL,
	                                           1,
	                                           radius_constraint,
	                                           radius,
	                                           true,
	                                           num_to_assign,
	                                           &priority_graph)) != SCC_ER_OK) {
		return ec;
	}

	if (radius_constraint) {
		for (size_t i = 0; i < clustering->num_data_points; ++i) {
			if (to_assign[i] && (priority_graph.tail_ptr[i] < priority_graph.tail_ptr[i + 1])) {
				assert(clustering->cluster_label[i] == SCC_CLABEL_NA);
				assert(priority_graph.tail_ptr[i] + 1 == priority_graph.tail_ptr[i + 1]);
				assert(!to_assign[priority_graph.head[priority_graph.tail_ptr[i]]]);
				assert(clustering->cluster_label[priority_graph.head[priority_graph.tail_ptr[i]]] != SCC_CLABEL_NA);
				clustering->cluster_label[i] = clustering->cluster_label[priority_graph.head[priority_graph.tail_ptr[i]]];
			}
		}
	} else {
		for (size_t i = 0; i < clustering->num_data_points; ++i) {
			if (to_assign[i]) {
				assert(clustering->cluster_label[i] == SCC_CLABEL_NA);
				assert(priority_graph.tail_ptr[i] + 1 == priority_graph.tail_ptr[i + 1]);
				assert(!to_assign[priority_graph.head[priority_graph.tail_ptr[i]]]);
				assert(clustering->cluster_label[priority_graph.head[priority_graph.tail_ptr[i]]] != SCC_CLABEL_NA);
				clustering->cluster_label[i] = clustering->cluster_label[priority_graph.head[priority_graph.tail_ptr[i]]];
			}
		}
	}

	iscc_free_digraph(&priority_graph);

	return iscc_no_error();
}


#ifdef SCC_STABLE_NNG

static int iscc_compare_Dpid(const void* const a, const void* const b)
{
    const iscc_Dpid arg1 = *(const iscc_Dpid* const)a;
    const iscc_Dpid arg2 = *(const iscc_Dpid* const)b;
    return (arg1 > arg2) - (arg1 < arg2);
}

static void iscc_sort_nng(iscc_Digraph* const nng)
{
	for (size_t v = 0; v < nng->vertices; ++v) {
		const size_t count = nng->tail_ptr[v + 1] - nng->tail_ptr[v];
		if (count > 1) {
			qsort(nng->head + nng->tail_ptr[v], count, sizeof(iscc_Dpid), iscc_compare_Dpid);
		}
	}
}

#endif // ifdef SCC_STABLE_NNG
