/* =============================================================================
 * scclust -- A C library for size-constrained clustering
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
#include "utilities.h"


// =============================================================================
// Static function prototypes
// =============================================================================

static scc_ErrorCode iscc_make_clustering_from_nng(scc_Clustering* clustering,
                                                   void* data_set,
                                                   iscc_Digraph* nng,
                                                   const scc_ClusterOptions* options);


// =============================================================================
// Public function implementations
// =============================================================================

scc_ErrorCode scc_sc_clustering(void* const data_set,
                                const scc_ClusterOptions* const options,
                                scc_Clustering* const out_clustering)
{
	if (!iscc_check_input_clustering(out_clustering)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid clustering object.");
	}
	if (!iscc_check_data_set(data_set)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid data set object.");
	}
	if (iscc_num_data_points(data_set) != out_clustering->num_data_points) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Number of data points in data set does not match clustering object.");
	}
	scc_ErrorCode ec;
	if ((ec = iscc_check_cluster_options(options, out_clustering->num_data_points)) != SCC_ER_OK) {
		return ec;
	}
	if (out_clustering->num_clusters != 0) {
		return iscc_make_error_msg(SCC_ER_NOT_IMPLEMENTED, "Cannot refine existing clusterings.");
	}

	if (options->seed_method == SCC_SM_BATCHES) {
		return scc_nng_clustering_batches(out_clustering,
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
		                                            out_clustering->num_data_points,
		                                            options->size_constraint,
		                                            options->len_primary_data_points,
		                                            options->primary_data_points,
		                                            (options->seed_radius == SCC_RM_USE_SUPPLIED),
		                                            options->seed_supplied_radius,
		                                            &nng)) != SCC_ER_OK) {
			return ec;
		}
	} else {
		assert(options->num_types <= UINT16_MAX);
		if ((ec = iscc_get_nng_with_type_constraint(data_set,
		                                            out_clustering->num_data_points,
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

	ec = iscc_make_clustering_from_nng(out_clustering,
	                                   data_set,
	                                   &nng,
	                                   options);

	iscc_free_digraph(&nng);

	return ec;
}


// =============================================================================
// Static function implementations
// =============================================================================

static scc_ErrorCode iscc_make_clustering_from_nng(scc_Clustering* const clustering,
                                                   void* const data_set,
                                                   iscc_Digraph* const nng,
                                                   const scc_ClusterOptions* options)
{
	assert(iscc_check_input_clustering(clustering));
	assert(iscc_check_data_set(data_set));
	assert(iscc_num_data_points(data_set) == clustering->num_data_points);
	assert(iscc_digraph_is_valid(nng));
	assert(!iscc_digraph_is_empty(nng));

	iscc_SeedResult seed_result = {
		.capacity = 1 + (clustering->num_data_points / options->size_constraint),
		.count = 0,
		.seeds = NULL,
	};

	scc_ErrorCode ec;
	if ((ec = iscc_find_seeds(nng, options->seed_method, &seed_result)) != SCC_ER_OK) {
		return ec;
	}

	scc_RadiusMethod primary_radius = options->primary_radius;
	double primary_supplied_radius = options->primary_supplied_radius;
	scc_RadiusMethod secondary_radius = options->secondary_radius;
	double secondary_supplied_radius = options->secondary_supplied_radius;

	if (options->primary_unassigned_method == SCC_UM_IGNORE) {
		primary_radius = SCC_RM_NO_RADIUS;
		primary_supplied_radius = 0.0;
	}

	if (options->secondary_unassigned_method == SCC_UM_IGNORE) {
		secondary_radius = SCC_RM_NO_RADIUS;
		secondary_supplied_radius = 0.0;
	}

	if (primary_radius == SCC_RM_USE_SEED_RADIUS) {
		primary_radius = options->seed_radius;
		primary_supplied_radius = options->seed_supplied_radius;
	}

	if (secondary_radius == SCC_RM_USE_SEED_RADIUS) {
		secondary_radius = options->seed_radius;
		secondary_supplied_radius = options->seed_supplied_radius;
	}

	// Estimate assign radius if we need to, and modify options
	if ((primary_radius == SCC_RM_USE_ESTIMATED) ||
	        (secondary_radius == SCC_RM_USE_ESTIMATED)) {
		double avg_seed_dist;
		if ((ec = iscc_estimate_avg_seed_dist(data_set,
		                                      &seed_result,
		                                      nng,
		                                      options->size_constraint,
		                                      &avg_seed_dist)) != SCC_ER_OK) {
			free(seed_result.seeds);
			return ec;
		}

		if (primary_radius == SCC_RM_USE_ESTIMATED) {
			if (avg_seed_dist > 0.0) {
				primary_radius = SCC_RM_USE_SUPPLIED;
				primary_supplied_radius = avg_seed_dist;
			} else {
				free(seed_result.seeds);
				return iscc_make_error_msg(SCC_ER_NO_SOLUTION, "Infeasible radius constraint.");
			}
		}

		if (secondary_radius == SCC_RM_USE_ESTIMATED) {
			if (avg_seed_dist > 0.0) {
				secondary_radius = SCC_RM_USE_SUPPLIED;
				secondary_supplied_radius = avg_seed_dist;
			} else {
				free(seed_result.seeds);
				return iscc_make_error_msg(SCC_ER_NO_SOLUTION, "Infeasible radius constraint.");
			}
		}
	}

	assert((primary_radius == SCC_RM_NO_RADIUS) || (primary_radius == SCC_RM_USE_SUPPLIED));
	assert((secondary_radius == SCC_RM_NO_RADIUS) || (secondary_radius == SCC_RM_USE_SUPPLIED));

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
	                                       data_set,
	                                       &seed_result,
	                                       nng,
	                                       (options->num_types < 2),
	                                       options->primary_unassigned_method,
	                                       (primary_radius == SCC_RM_USE_SUPPLIED),
	                                       primary_supplied_radius,
	                                       options->len_primary_data_points,
	                                       options->primary_data_points,
	                                       options->secondary_unassigned_method,
	                                       (secondary_radius == SCC_RM_USE_SUPPLIED),
	                                       secondary_supplied_radius);

	free(seed_result.seeds);
	return ec;
}
