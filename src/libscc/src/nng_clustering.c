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
#include "digraph_core.h"
#include "dist_search.h"
#include "error.h"
#include "nng_core.h"
#include "nng_findseeds.h"
#include "scclust_internal.h"


// =============================================================================
// Internal variables
// =============================================================================

static const scc_SeedMethod SCC_MAX_SEED_METHOD = SCC_SM_EXCLUSION_UPDATING;
static const scc_UnassignedMethod SCC_MAX_UNASSIGNED_METHOD = SCC_UM_CLOSEST_SEED_EST_RADIUS;


// =============================================================================
// Internal function prototypes
// =============================================================================

static scc_ErrorCode iscc_make_clustering_from_nng(scc_Clustering* clustering,
                                                   void* data_set_object,
                                                   iscc_Digraph* nng,
                                                   bool nng_is_ordered,
                                                   uint32_t size_constraint,
                                                   scc_SeedMethod seed_method,
                                                   scc_UnassignedMethod main_unassigned_method,
                                                   bool main_radius_constraint,
                                                   double main_radius,
                                                   const bool main_data_points[],
                                                   scc_UnassignedMethod secondary_unassigned_method,
                                                   bool secondary_radius_constraint,
                                                   double secondary_radius);


// =============================================================================
// External function implementations
// =============================================================================

scc_ErrorCode scc_nng_clustering(scc_Clustering* const clustering,
                                 void* const data_set_object,
                                 const uint32_t size_constraint,
                                 const scc_SeedMethod seed_method,
                                 const scc_UnassignedMethod main_unassigned_method,
                                 const bool main_radius_constraint,
                                 const double main_radius,
                                 const size_t len_main_data_points,
                                 const bool main_data_points[const],
                                 const scc_UnassignedMethod secondary_unassigned_method,
                                 const bool secondary_radius_constraint,
                                 const double secondary_radius)
{
	if (!iscc_check_input_clustering(clustering)) return iscc_make_error(SCC_ER_INVALID_CLUSTERING);
	if (!iscc_check_data_set_object(data_set_object, clustering->num_data_points)) return iscc_make_error(SCC_ER_INVALID_DATA_OBJ);
	if (size_constraint < 2) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (clustering->num_data_points < size_constraint) return iscc_make_error(SCC_ER_NO_CLUST_EXIST_CONSTRAINT);
	if (seed_method > SCC_MAX_SEED_METHOD) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (main_unassigned_method > SCC_MAX_UNASSIGNED_METHOD) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (main_radius_constraint && (main_radius <= 0.0)) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if ((main_data_points != NULL) && (len_main_data_points < clustering->num_data_points)) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if ((main_data_points == NULL) && (secondary_unassigned_method != SCC_UM_IGNORE)) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (secondary_unassigned_method > SCC_MAX_UNASSIGNED_METHOD) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (secondary_unassigned_method == SCC_UM_ASSIGN_BY_NNG) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (secondary_radius_constraint && (secondary_radius <= 0.0)) return iscc_make_error(SCC_ER_INVALID_INPUT);

	if (clustering->num_clusters != 0) return iscc_make_error(SCC_ER_NOT_IMPLEMENTED);

	scc_ErrorCode ec;
	iscc_Digraph nng;
	if ((ec = iscc_get_nng_with_size_constraint(data_set_object,
	                                            clustering->num_data_points,
	                                            size_constraint,
	                                            main_data_points,
	                                            main_radius_constraint,
	                                            main_radius,
	                                            &nng)) != SCC_ER_OK) {
		return ec;
	}

	assert(!iscc_digraph_is_empty(&nng));

	ec = iscc_make_clustering_from_nng(clustering,
	                                   data_set_object,
	                                   &nng,
	                                   true,
	                                   size_constraint,
	                                   seed_method,
	                                   main_unassigned_method,
	                                   main_radius_constraint,
	                                   main_radius,
	                                   main_data_points,
	                                   secondary_unassigned_method,
	                                   secondary_radius_constraint,
	                                   secondary_radius);

	iscc_free_digraph(&nng);

	return ec;
}


scc_ErrorCode scc_nng_clustering_types(scc_Clustering* const clustering,
                                       void* const data_set_object,
                                       const uint32_t size_constraint,
                                       const uintmax_t num_types,
                                       const uint32_t type_size_constraints[const],
                                       const size_t len_type_labels,
                                       const scc_TypeLabel type_labels[const],
                                       const scc_SeedMethod seed_method,
                                       const scc_UnassignedMethod main_unassigned_method,
                                       const bool main_radius_constraint,
                                       const double main_radius,
                                       const size_t len_main_data_points,
                                       const bool main_data_points[const],
                                       const scc_UnassignedMethod secondary_unassigned_method,
                                       const bool secondary_radius_constraint,
                                       const double secondary_radius)
{
	if (!iscc_check_input_clustering(clustering)) return iscc_make_error(SCC_ER_INVALID_CLUSTERING);
	if (!iscc_check_data_set_object(data_set_object, clustering->num_data_points)) return iscc_make_error(SCC_ER_INVALID_DATA_OBJ);
	if (size_constraint < 2) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (clustering->num_data_points < size_constraint) return iscc_make_error(SCC_ER_NO_CLUST_EXIST_CONSTRAINT);
	if (num_types < 2) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (num_types > ISCC_TYPELABEL_MAX) return iscc_make_error(SCC_ER_TOO_LARGE_PROBLEM);
	if (num_types > UINT_FAST16_MAX) return iscc_make_error(SCC_ER_TOO_LARGE_PROBLEM);
	if (type_size_constraints == NULL) return iscc_make_error(SCC_ER_NULL_INPUT);
	if (len_type_labels < clustering->num_data_points) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (type_labels == NULL) return iscc_make_error(SCC_ER_NULL_INPUT);
	if (seed_method > SCC_MAX_SEED_METHOD) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (main_unassigned_method > SCC_MAX_UNASSIGNED_METHOD) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (main_radius_constraint && (main_radius <= 0.0)) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if ((main_data_points != NULL) && (len_main_data_points < clustering->num_data_points)) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if ((main_data_points == NULL) && (secondary_unassigned_method != SCC_UM_IGNORE)) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (secondary_unassigned_method > SCC_MAX_UNASSIGNED_METHOD) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (secondary_unassigned_method == SCC_UM_ASSIGN_BY_NNG) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (secondary_radius_constraint && (secondary_radius <= 0.0)) return iscc_make_error(SCC_ER_INVALID_INPUT);

	if (clustering->num_clusters != 0) return iscc_make_error(SCC_ER_NOT_IMPLEMENTED);

	assert(num_types <= UINT_FAST16_MAX);
	const uint_fast16_t num_types_f16 = (uint_fast16_t) num_types;

	scc_ErrorCode ec;
	iscc_Digraph nng;
	if ((ec = iscc_get_nng_with_type_constraint(data_set_object,
	                                            clustering->num_data_points,
	                                            size_constraint,
	                                            num_types_f16,
	                                            type_size_constraints,
	                                            type_labels,
	                                            main_data_points,
	                                            main_radius_constraint,
	                                            main_radius,
	                                            &nng)) != SCC_ER_OK) {
		return ec;
	}

	assert(!iscc_digraph_is_empty(&nng));

	ec = iscc_make_clustering_from_nng(clustering,
	                                   data_set_object,
	                                   &nng,
	                                   false,
	                                   size_constraint,
	                                   seed_method,
	                                   main_unassigned_method,
	                                   main_radius_constraint,
	                                   main_radius,
	                                   main_data_points,
	                                   secondary_unassigned_method,
	                                   secondary_radius_constraint,
	                                   secondary_radius);

	iscc_free_digraph(&nng);

	return ec;
}


// =============================================================================
// Internal function implementations
// =============================================================================

static scc_ErrorCode iscc_make_clustering_from_nng(scc_Clustering* const clustering,
                                                   void* const data_set_object,
                                                   iscc_Digraph* const nng,
                                                   const bool nng_is_ordered,
                                                   const uint32_t size_constraint,
                                                   const scc_SeedMethod seed_method,
                                                   scc_UnassignedMethod main_unassigned_method,
                                                   bool main_radius_constraint,
                                                   double main_radius,
                                                   const bool main_data_points[const],
                                                   scc_UnassignedMethod secondary_unassigned_method,
                                                   bool secondary_radius_constraint,
                                                   double secondary_radius)
{
	assert(iscc_check_input_clustering(clustering));
	assert(iscc_check_data_set_object(data_set_object, clustering->num_data_points));
	assert(iscc_digraph_is_valid(nng));
	assert(!iscc_digraph_is_empty(nng));
	assert(size_constraint >= 2);
	assert(seed_method <= SCC_MAX_SEED_METHOD);
	assert((main_unassigned_method == SCC_UM_IGNORE) ||
	       (main_unassigned_method == SCC_UM_ASSIGN_BY_NNG) ||
	       (main_unassigned_method == SCC_UM_CLOSEST_ASSIGNED) ||
	       (main_unassigned_method == SCC_UM_CLOSEST_SEED) ||
	       (main_unassigned_method == SCC_UM_CLOSEST_SEED_EST_RADIUS));
	assert(!main_radius_constraint || (main_radius > 0.0));
	assert((main_data_points != NULL) || (secondary_unassigned_method == SCC_UM_IGNORE));
	assert((secondary_unassigned_method == SCC_UM_IGNORE) ||
	       (secondary_unassigned_method == SCC_UM_CLOSEST_ASSIGNED) ||
	       (secondary_unassigned_method == SCC_UM_CLOSEST_SEED) ||
	       (secondary_unassigned_method == SCC_UM_CLOSEST_SEED_EST_RADIUS));
	assert(!secondary_radius_constraint || (secondary_radius > 0.0));

	iscc_SeedResult seed_result = {
		.capacity = 1 + (clustering->num_data_points / size_constraint),
		.count = 0,
		.seeds = NULL,
	};

	scc_ErrorCode ec;
	if ((ec = iscc_find_seeds(nng, seed_method, &seed_result)) != SCC_ER_OK) {
		return ec;
	}

	// Estimate assign radius if we need to, and modify options
	if ((main_unassigned_method == SCC_UM_CLOSEST_SEED_EST_RADIUS) ||
	        (secondary_unassigned_method == SCC_UM_CLOSEST_SEED_EST_RADIUS)) {
		double avg_seed_dist;
		if ((ec = iscc_estimate_avg_seed_dist(data_set_object,
		                                      &seed_result,
		                                      nng,
		                                      size_constraint,
		                                      &avg_seed_dist)) != SCC_ER_OK) {
			free(seed_result.seeds);
			return ec;
		}

		if (main_unassigned_method == SCC_UM_CLOSEST_SEED_EST_RADIUS) {
			if (avg_seed_dist > 0.0) {
				main_unassigned_method = SCC_UM_CLOSEST_SEED;
				main_radius_constraint = true;
				main_radius = avg_seed_dist;
			} else {
				free(seed_result.seeds);
				return iscc_make_error(SCC_ER_NOT_IMPLEMENTED);
			}
		}

		if (secondary_unassigned_method == SCC_UM_CLOSEST_SEED_EST_RADIUS) {
			if (avg_seed_dist > 0.0) {
				secondary_unassigned_method = SCC_UM_CLOSEST_SEED;
				secondary_radius_constraint = true;
				secondary_radius = avg_seed_dist;
			} else {
				free(seed_result.seeds);
				return iscc_make_error(SCC_ER_NOT_IMPLEMENTED);
			}
		}
	}

	// Initialize cluster labels
	if (clustering->cluster_label == NULL) {
		clustering->external_labels = false;
		clustering->cluster_label = malloc(sizeof(scc_Clabel[clustering->num_data_points]));
		if (clustering->cluster_label == NULL) {
			free(seed_result.seeds);
			return iscc_make_error(SCC_ER_NO_MEMORY);
		}
	}

	ec = iscc_make_nng_clusters_from_seeds(clustering,
	                                       data_set_object,
	                                       &seed_result,
	                                       nng,
	                                       nng_is_ordered,
	                                       main_unassigned_method,
	                                       main_radius_constraint,
	                                       main_radius,
	                                       main_data_points,
	                                       secondary_unassigned_method,
	                                       secondary_radius_constraint,
	                                       secondary_radius);

	free(seed_result.seeds);
	return ec;
}
