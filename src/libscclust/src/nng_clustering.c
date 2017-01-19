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
#include "clustering_struct.h"
#include "digraph_core.h"
#include "dist_search.h"
#include "error.h"
#include "nng_batch_clustering.h"
#include "nng_core.h"
#include "nng_findseeds.h"


// =============================================================================
// Internal variables
// =============================================================================

#define ISCC_M_OPTIONS_STRUCT_VERSION 722678001
static const int32_t ISCC_OPTIONS_STRUCT_VERSION = ISCC_M_OPTIONS_STRUCT_VERSION;

const scc_ClusterOptions scc_default_cluster_options = {
	.options_version = ISCC_M_OPTIONS_STRUCT_VERSION, // GCC error if not init with macro
	.size_constraint = 0,
	.num_types = 0,
	.type_constraints = NULL,
	.len_type_labels = 0,
	.type_labels = NULL,
	.seed_method = SCC_SM_LEXICAL,
	.len_primary_data_points = 0,
	.primary_data_points = NULL,
	.primary_unassigned_method = SCC_UM_ANY_NEIGHBOR,
	.secondary_unassigned_method = SCC_UM_IGNORE,
	.seed_radius = SCC_RM_NO_RADIUS,
	.seed_supplied_radius = 0.0,
	.primary_radius = SCC_RM_USE_SEED_RADIUS,
	.primary_supplied_radius = 0.0,
	.secondary_radius = SCC_RM_USE_SEED_RADIUS,
	.secondary_supplied_radius = 0.0,
	.batch_size = 0,
};


// =============================================================================
// Internal function prototypes
// =============================================================================

static scc_ErrorCode iscc_check_cluster_options(const scc_ClusterOptions* options,
                                                size_t num_data_points);

static scc_ErrorCode iscc_make_clustering_from_nng(scc_Clustering* clustering,
                                                   void* data_set,
                                                   iscc_Digraph* nng,
                                                   const scc_ClusterOptions* options);


// =============================================================================
// External function implementations
// =============================================================================

scc_ErrorCode scc_make_clustering(void* const data_set,
                                  scc_Clustering* const clustering,
                                  const scc_ClusterOptions* const options)
{
	if (!iscc_check_input_clustering(clustering)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid clustering object.");
	}
	if (!iscc_check_data_set(data_set, clustering->num_data_points)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid data set object.");
	}
	scc_ErrorCode ec;
	if ((ec = iscc_check_cluster_options(options, clustering->num_data_points)) != SCC_ER_OK) {
		return ec;
	}
	if (clustering->num_clusters != 0) {
		return iscc_make_error_msg(SCC_ER_NOT_IMPLEMENTED, "Cannot refine existing clusterings.");
	}

	if (options->seed_method == SCC_SM_BATCHES) {
		return scc_nng_clustering_batches(clustering,
		                                  data_set,
		                                  options->size_constraint,
		                                  options->primary_unassigned_method,
		                                  (options->seed_radius == SCC_RM_USE_SUPPLIED),
		                                  options->seed_supplied_radius,
		                                  options->len_primary_data_points,
		                                  options->primary_data_points,
		                                  options->batch_size);
	}

	iscc_Digraph nng;
	if (options->num_types < 2) {
		if ((ec = iscc_get_nng_with_size_constraint(data_set,
		                                            clustering->num_data_points,
		                                            options->size_constraint,
		                                            options->len_primary_data_points,
		                                            options->primary_data_points,
		                                            (options->seed_radius == SCC_RM_USE_SUPPLIED),
		                                            options->seed_supplied_radius,
		                                            &nng)) != SCC_ER_OK) {
			return ec;
		}
	} else {
		assert(options->num_types <= UINT_FAST16_MAX);
		if ((ec = iscc_get_nng_with_type_constraint(data_set,
		                                            clustering->num_data_points,
		                                            options->size_constraint,
		                                            (uint_fast16_t) options->num_types,
		                                            options->type_constraints,
		                                            options->type_labels,
		                                            options->len_primary_data_points,
		                                            options->primary_data_points,
		                                            (options->seed_radius == SCC_RM_USE_SUPPLIED),
		                                            options->seed_supplied_radius,
		                                            &nng)) != SCC_ER_OK) {
			return ec;
		}
	}

	assert(!iscc_digraph_is_empty(&nng));

	ec = iscc_make_clustering_from_nng(clustering,
	                                   data_set,
	                                   &nng,
	                                   options);

	iscc_free_digraph(&nng);

	return ec;
}


// =============================================================================
// Internal function implementations
// =============================================================================

static scc_ErrorCode iscc_check_cluster_options(const scc_ClusterOptions* const options,
                                                const size_t num_data_points)
{
	if (options->options_version != ISCC_OPTIONS_STRUCT_VERSION) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Incompatible scc_ClusterOptions version.");
	}
	if (options->size_constraint < 2) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Size constraint must be 2 or greater.");
	}
	if (num_data_points < options->size_constraint) {
		return iscc_make_error_msg(SCC_ER_NO_SOLUTION, "Fewer data points than size constraint.");
	}

	if (options->num_types < 2) {
		if (options->type_constraints != NULL) {
			return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid type constraints.");
		}
		if (options->len_type_labels != 0) {
			return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid type labels.");
		}
		if (options->type_labels != NULL) {
			return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid type labels.");
		}
	} else {
		if (options->num_types > ISCC_TYPELABEL_MAX) {
			return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many data point types.");
		}
		if (options->num_types > UINT_FAST16_MAX) {
			return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many data point types.");
		}
		if (options->type_constraints == NULL) {
			return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid type constraints.");
		}
		if (options->len_type_labels < num_data_points) {
			return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid type labels.");
		}
		if (options->type_labels == NULL) {
			return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid type labels.");
		}
	}

	if ((options->seed_method != SCC_SM_LEXICAL) &&
			(options->seed_method != SCC_SM_INWARDS_ORDER) &&
			(options->seed_method != SCC_SM_INWARDS_UPDATING) &&
			(options->seed_method != SCC_SM_INWARDS_ALT_UPDATING) &&
			(options->seed_method != SCC_SM_EXCLUSION_ORDER) &&
			(options->seed_method != SCC_SM_EXCLUSION_UPDATING) &&
			(options->seed_method != SCC_SM_BATCHES)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Unknown seed method.");
	}
	if ((options->primary_data_points != NULL) && (options->len_primary_data_points == 0)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid primary data points.");
	}
	if (options->primary_data_points != NULL) {
		for (size_t i = 1; i < options->len_primary_data_points; ++i) {
			if (options->primary_data_points[i - 1] >= options->primary_data_points[i]) {
				return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "`primary_data_points` is not sorted.");
			}
		}
	}
	if ((options->primary_data_points == NULL) && (options->len_primary_data_points > 0)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid primary data points.");
	}
	if ((options->primary_data_points == NULL) && (options->secondary_unassigned_method != SCC_UM_IGNORE)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid unassigned method.");
	}

	if ((options->primary_unassigned_method != SCC_UM_IGNORE) &&
			(options->primary_unassigned_method != SCC_UM_ANY_NEIGHBOR) &&
			(options->primary_unassigned_method != SCC_UM_CLOSEST_ASSIGNED) &&
			(options->primary_unassigned_method != SCC_UM_CLOSEST_SEED)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Unknown unassigned method.");
	}
	if (options->secondary_unassigned_method == SCC_UM_ANY_NEIGHBOR) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid unassigned method.");
	}
	if ((options->secondary_unassigned_method != SCC_UM_IGNORE) &&
			(options->secondary_unassigned_method != SCC_UM_CLOSEST_ASSIGNED) &&
			(options->secondary_unassigned_method != SCC_UM_CLOSEST_SEED)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Unknown unassigned method.");
	}
	if ((options->seed_radius != SCC_RM_NO_RADIUS) &&
			(options->seed_radius != SCC_RM_USE_SUPPLIED)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid radius method.");
	}
	if ((options->seed_radius == SCC_RM_USE_SUPPLIED) && (options->seed_supplied_radius <= 0.0)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid radius.");
	}
	if ((options->primary_radius != SCC_RM_NO_RADIUS) &&
			(options->primary_radius != SCC_RM_USE_SUPPLIED) &&
			(options->primary_radius != SCC_RM_USE_SEED_RADIUS) &&
			(options->primary_radius != SCC_RM_USE_ESTIMATED)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid radius method.");
	}
	if ((options->primary_radius == SCC_RM_USE_SUPPLIED) && (options->primary_supplied_radius <= 0.0)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid radius.");
	}
	if ((options->secondary_radius != SCC_RM_NO_RADIUS) &&
			(options->secondary_radius != SCC_RM_USE_SUPPLIED) &&
			(options->secondary_radius != SCC_RM_USE_SEED_RADIUS) &&
			(options->secondary_radius != SCC_RM_USE_ESTIMATED)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid radius method.");
	}
	if ((options->secondary_radius == SCC_RM_USE_SUPPLIED) && (options->secondary_supplied_radius <= 0.0)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid radius.");
	}

	if (options->seed_method == SCC_SM_BATCHES) {
		if (options->num_types >= 2) {
			return iscc_make_error_msg(SCC_ER_NOT_IMPLEMENTED, "SCC_SM_BATCHES cannot be used with type constraints.");
		}
		if (options->secondary_unassigned_method != SCC_UM_IGNORE) {
			return iscc_make_error_msg(SCC_ER_NOT_IMPLEMENTED, "SCC_SM_BATCHES must be used with `secondary_unassigned_method = SCC_UM_IGNORE`.");
		}
		if (options->primary_radius != SCC_RM_USE_SEED_RADIUS) {
			return iscc_make_error_msg(SCC_ER_NOT_IMPLEMENTED, "SCC_SM_BATCHES must be used with `primary_radius = SCC_RM_USE_SEED_RADIUS`.");
		}
	}

	return iscc_no_error();
}


static scc_ErrorCode iscc_make_clustering_from_nng(scc_Clustering* const clustering,
                                                   void* const data_set,
                                                   iscc_Digraph* const nng,
                                                   const scc_ClusterOptions* options)
{
	assert(iscc_check_input_clustering(clustering));
	assert(iscc_check_data_set(data_set, clustering->num_data_points));
	assert(iscc_digraph_is_valid(nng));
	assert(!iscc_digraph_is_empty(nng));

	const bool nng_is_ordered = (options->num_types < 2);
	const uint32_t size_constraint = options->size_constraint;
	const scc_SeedMethod seed_method = options->seed_method;
	bool* primary_data_points = NULL;
	if (options->primary_data_points != NULL) {
		primary_data_points = calloc(clustering->num_data_points, sizeof(bool));
		for (size_t i = 0; i < options->len_primary_data_points; ++i) {
			primary_data_points[options->primary_data_points[i]] = true;
		}
	}
	scc_UnassignedMethod primary_unassigned_method = options->primary_unassigned_method;
	scc_UnassignedMethod secondary_unassigned_method = options->secondary_unassigned_method;
	scc_RadiusMethod seed_radius = options->seed_radius;
	double seed_supplied_radius = options->seed_supplied_radius;
	scc_RadiusMethod primary_radius = options->primary_radius;
	double primary_supplied_radius = options->primary_supplied_radius;
	scc_RadiusMethod secondary_radius = options->secondary_radius;
	double secondary_supplied_radius = options->secondary_supplied_radius;

	iscc_SeedResult seed_result = {
		.capacity = 1 + (clustering->num_data_points / size_constraint),
		.count = 0,
		.seeds = NULL,
	};

	scc_ErrorCode ec;
	if ((ec = iscc_find_seeds(nng, seed_method, &seed_result)) != SCC_ER_OK) {
		free(primary_data_points);
		return ec;
	}

	// Estimate assign radius if we need to, and modify options
	if ((primary_radius == SCC_RM_USE_ESTIMATED) ||
	        (secondary_radius == SCC_RM_USE_ESTIMATED)) {
		double avg_seed_dist;
		if ((ec = iscc_estimate_avg_seed_dist(data_set,
		                                      &seed_result,
		                                      nng,
		                                      size_constraint,
		                                      &avg_seed_dist)) != SCC_ER_OK) {
			free(seed_result.seeds);
			free(primary_data_points);
			return ec;
		}

		if (primary_radius == SCC_RM_USE_ESTIMATED) {
			if (avg_seed_dist > 0.0) {
				primary_radius = SCC_RM_USE_SUPPLIED;
				primary_supplied_radius = avg_seed_dist;
			} else {
				free(seed_result.seeds);
				free(primary_data_points);
				return iscc_make_error_msg(SCC_ER_NO_SOLUTION, "Infeasible radius constraint.");
			}
		}

		if (secondary_radius == SCC_RM_USE_ESTIMATED) {
			if (avg_seed_dist > 0.0) {
				secondary_radius = SCC_RM_USE_SUPPLIED;
				secondary_supplied_radius = avg_seed_dist;
			} else {
				free(seed_result.seeds);
				free(primary_data_points);
				return iscc_make_error_msg(SCC_ER_NO_SOLUTION, "Infeasible radius constraint.");
			}
		}
	}

	if (primary_radius == SCC_RM_USE_SEED_RADIUS) {
		primary_radius = seed_radius;
		primary_supplied_radius = seed_supplied_radius;
	}

	if (secondary_radius == SCC_RM_USE_SEED_RADIUS) {
		secondary_radius = seed_radius;
		secondary_supplied_radius = seed_supplied_radius;
	}

	assert((primary_radius == SCC_RM_NO_RADIUS) || (primary_radius == SCC_RM_USE_SUPPLIED));
	assert((secondary_radius == SCC_RM_NO_RADIUS) || (secondary_radius == SCC_RM_USE_SUPPLIED));

	// Initialize cluster labels
	if (clustering->cluster_label == NULL) {
		clustering->external_labels = false;
		clustering->cluster_label = malloc(sizeof(scc_Clabel[clustering->num_data_points]));
		if (clustering->cluster_label == NULL) {
			free(seed_result.seeds);
			free(primary_data_points);
			return iscc_make_error(SCC_ER_NO_MEMORY);
		}
	}

	ec = iscc_make_nng_clusters_from_seeds(clustering,
	                                       data_set,
	                                       &seed_result,
	                                       nng,
	                                       nng_is_ordered,
	                                       primary_unassigned_method,
	                                       (primary_radius == SCC_RM_USE_SUPPLIED),
	                                       primary_supplied_radius,
	                                       primary_data_points,
	                                       secondary_unassigned_method,
	                                       (secondary_radius == SCC_RM_USE_SUPPLIED),
	                                       secondary_supplied_radius);

	free(seed_result.seeds);
	free(primary_data_points);
	return ec;
}
