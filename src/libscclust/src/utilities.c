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

#include "../include/scclust.h"
#include "utilities.h"

#include <assert.h>
#include <float.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include "clustering_struct.h"
#include "dist_search.h"
#include "error.h"
#include "scclust_types.h"


// =============================================================================
// Internal variables
// =============================================================================

/** The null clustering statistics struct.
 *
 *  This is an easily detectable invalid struct used as return value on errors.
 */
static const scc_ClusteringStats ISCC_NULL_CLUSTERING_STATS = { 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

static const int32_t ISCC_OPTIONS_STRUCT_VERSION = 722678001;


// =============================================================================
// Public function implementations
// =============================================================================

scc_ClusterOptions scc_get_default_options(void)
{
	return (scc_ClusterOptions) {
		.options_version = ISCC_OPTIONS_STRUCT_VERSION,
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
}


scc_ErrorCode scc_check_clustering(const scc_Clustering* const clustering,
                                   const scc_ClusterOptions* const options,
                                   bool* const out_is_OK)
{
	if (out_is_OK == NULL) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Output parameter may not be NULL.");
	}
	*out_is_OK = false;
	if (!iscc_check_input_clustering(clustering)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid clustering object.");
	}
	if (clustering->num_clusters == 0) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Empty clustering.");
	}
	scc_ErrorCode ec;
	if ((ec = iscc_check_cluster_options(options, clustering->num_data_points)) != SCC_ER_OK) {
		return ec;
	}

	const uint32_t size_constraint = options->size_constraint;
	const uintmax_t num_types = options->num_types;
	const uint32_t* const type_constraints = options->type_constraints;
	const scc_TypeLabel* const type_labels = options->type_labels;

	assert(clustering->num_clusters <= ((uintmax_t) SCC_CLABEL_MAX));
	const scc_Clabel max_cluster = (scc_Clabel) clustering->num_clusters;
	for (size_t i = 0; i < clustering->num_data_points; ++i) {
		if ((clustering->cluster_label[i] > 0) && (clustering->cluster_label[i] < max_cluster)) continue;
		if (clustering->cluster_label[i] == 0) continue; // Since `scc_Clabel` can be unsigned
		if (clustering->cluster_label[i] == SCC_CLABEL_NA) continue;
		return iscc_no_error(); // Error found, return. (`out_is_OK` is set to false)
	}

	if (options->primary_data_points != NULL) {
		for (size_t i = 0; i < options->len_primary_data_points; ++i) {
			if (clustering->cluster_label[options->primary_data_points[i]] == SCC_CLABEL_NA) {
				return iscc_no_error(); // Error found, return. (`out_is_OK` is set to false)
			}
		}
	}

	if (num_types < 2) {

		size_t* const cluster_sizes = calloc(clustering->num_clusters, sizeof(size_t));
		if (cluster_sizes == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);

		for (size_t i = 0; i < clustering->num_data_points; ++i) {
			if (clustering->cluster_label[i] != SCC_CLABEL_NA) {
				++cluster_sizes[clustering->cluster_label[i]];
			}
		}

		for (size_t i = 0; i < clustering->num_clusters; ++i) {
			if (cluster_sizes[i] < size_constraint) {
				free(cluster_sizes);
				return iscc_no_error(); // Error found, return. (`out_is_OK` is set to false)
			}
		}

		free(cluster_sizes);

	} else { // num_types >= 2

		size_t* const cluster_type_sizes = calloc(num_types * clustering->num_clusters, sizeof(size_t));
		if (cluster_type_sizes == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);

		for (size_t i = 0; i < clustering->num_data_points; ++i) {
			if (clustering->cluster_label[i] != SCC_CLABEL_NA) {
				++cluster_type_sizes[(((size_t) clustering->cluster_label[i]) * num_types) + ((size_t) type_labels[i])];
			}
		}

		for (size_t i = 0; i < clustering->num_clusters; ++i) {
			size_t tmp_total_size = 0;
			for (size_t t = 0; t < num_types; ++t) {
				tmp_total_size += cluster_type_sizes[(i * num_types) + t];
				if (cluster_type_sizes[(i * num_types) + t] < type_constraints[t]) {
					free(cluster_type_sizes);
					return iscc_no_error(); // Error found, return. (`out_is_OK` is set to false)
				}
			}
			if (tmp_total_size < size_constraint) {
				free(cluster_type_sizes);
				return iscc_no_error(); // Error found, return. (`out_is_OK` is set to false)
			}
		}

		free(cluster_type_sizes);

	}

	*out_is_OK = true;
	return iscc_no_error();
}


scc_ErrorCode scc_get_clustering_stats(void* const data_set,
                                       const scc_Clustering* const clustering,
                                       scc_ClusteringStats* const out_stats)
{
	if (out_stats == NULL) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Output parameter may not be NULL.");
	}
	*out_stats = ISCC_NULL_CLUSTERING_STATS;
	if (!iscc_check_input_clustering(clustering)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid clustering object.");
	}
	if (clustering->num_clusters == 0) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Empty clustering.");
	}
	if (!iscc_check_data_set(data_set)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid data set object.");
	}
	if (iscc_num_data_points(data_set) != clustering->num_data_points) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Number of data points in data set does not match clustering object.");
	}

	size_t* const cluster_size = calloc(clustering->num_clusters, sizeof(size_t));
	if (cluster_size == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);

	for (size_t i = 0; i < clustering->num_data_points; ++i) {
		if (clustering->cluster_label[i] != SCC_CLABEL_NA) {
			++cluster_size[clustering->cluster_label[i]];
		}
	}

	scc_ClusteringStats tmp_stats = {
		.num_data_points = clustering->num_data_points,
		.num_assigned = 0,
		.num_clusters = clustering->num_clusters,
		.num_populated_clusters = 0,
		.min_cluster_size = UINT64_MAX,
		.max_cluster_size = 0,
		.avg_cluster_size = 0.0,
		.sum_dists = 0.0,
		.min_dist = DBL_MAX,
		.max_dist = 0.0,
		.cl_avg_min_dist = 0.0,
		.cl_avg_max_dist = 0.0,
		.cl_avg_dist_weighted = 0.0,
		.cl_avg_dist_unweighted = 0.0,
	};

	for (size_t c = 0; c < clustering->num_clusters; ++c) {
		if (cluster_size[c] == 0) continue;
		++tmp_stats.num_populated_clusters;
		tmp_stats.num_assigned += cluster_size[c];
		if (tmp_stats.min_cluster_size > cluster_size[c]) {
			tmp_stats.min_cluster_size = cluster_size[c];
		}
		if (tmp_stats.max_cluster_size < cluster_size[c]) {
			tmp_stats.max_cluster_size = cluster_size[c];
		}
	}

	if (tmp_stats.num_populated_clusters == 0) {
		free(cluster_size);
		*out_stats = tmp_stats;
		return iscc_no_error();
	}

	const size_t largest_dist_matrix = (tmp_stats.max_cluster_size * (tmp_stats.max_cluster_size - 1)) / 2;
	scc_PointIndex* const id_store = malloc(sizeof(scc_PointIndex[tmp_stats.num_assigned]));
	scc_PointIndex** const cl_members = malloc(sizeof(scc_PointIndex*[clustering->num_clusters]));
	double* const dist_scratch = malloc(sizeof(double[largest_dist_matrix]));
	if ((id_store == NULL) || (cl_members == NULL) || (dist_scratch == NULL)) {
		free(cluster_size);
		free(id_store);
		free(cl_members);
		free(dist_scratch);
		return iscc_make_error(SCC_ER_NO_MEMORY);
	}

	cl_members[0] = id_store + cluster_size[0];
	for (size_t c = 1; c < clustering->num_clusters; ++c) {
		cl_members[c] = cl_members[c - 1] + cluster_size[c];
	}

	assert(clustering->num_data_points <= ISCC_POINTINDEX_MAX);
	const scc_PointIndex num_data_points = (scc_PointIndex) clustering->num_data_points; // If `scc_PointIndex` is signed
	for (scc_PointIndex i = 0; i < num_data_points; ++i) {
		if (clustering->cluster_label[i] != SCC_CLABEL_NA) {
			--cl_members[clustering->cluster_label[i]];
			*(cl_members[clustering->cluster_label[i]]) = i;
		}
	}

	for (size_t c = 0; c < clustering->num_clusters; ++c) {
		if (cluster_size[c] < 2) {
			if (cluster_size[c] == 1) tmp_stats.min_dist = 0.0;
			continue;
		}

		const size_t size_dist_matrix = (cluster_size[c] * (cluster_size[c] - 1)) / 2;
		if (!iscc_get_dist_matrix(data_set, cluster_size[c], cl_members[c], dist_scratch)) {
			free(cluster_size);
			free(id_store);
			free(cl_members);
			free(dist_scratch);
			return iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
		}

		double cluster_sum_dists = dist_scratch[0];
		double cluster_min = dist_scratch[0];
		double cluster_max = dist_scratch[0];

		for (size_t d = 1; d < size_dist_matrix; ++d) {
			cluster_sum_dists += dist_scratch[d];
			if (cluster_min > dist_scratch[d]) {
				cluster_min = dist_scratch[d];
			}
			if (cluster_max < dist_scratch[d]) {
				cluster_max = dist_scratch[d];
			}
		}

		tmp_stats.sum_dists += cluster_sum_dists;

		if (tmp_stats.min_dist > cluster_min) {
			tmp_stats.min_dist = cluster_min;
		}
		if (tmp_stats.max_dist < cluster_max) {
			tmp_stats.max_dist = cluster_max;
		}
		tmp_stats.cl_avg_min_dist += cluster_min;
		tmp_stats.cl_avg_max_dist += cluster_max;

		tmp_stats.cl_avg_dist_weighted += ((double) cluster_size[c]) * cluster_sum_dists / ((double) size_dist_matrix);
		tmp_stats.cl_avg_dist_unweighted += cluster_sum_dists / ((double) size_dist_matrix);
	}

	tmp_stats.avg_cluster_size = ((double) tmp_stats.num_assigned) / ((double) tmp_stats.num_populated_clusters);
	tmp_stats.cl_avg_min_dist = tmp_stats.cl_avg_min_dist / ((double) tmp_stats.num_populated_clusters);
	tmp_stats.cl_avg_max_dist = tmp_stats.cl_avg_max_dist / ((double) tmp_stats.num_populated_clusters);
	tmp_stats.cl_avg_dist_weighted = tmp_stats.cl_avg_dist_weighted / ((double) tmp_stats.num_assigned);
	tmp_stats.cl_avg_dist_unweighted = tmp_stats.cl_avg_dist_unweighted / ((double) tmp_stats.num_populated_clusters);

	free(cluster_size);
	free(id_store);
	free(cl_members);
	free(dist_scratch);

	*out_stats = tmp_stats;

	return iscc_no_error();
}


scc_ErrorCode scc_get_cluster_seeds(void* const data_set,
                                    const scc_ClusterOptions* const options,
                                    scc_SeedVector* const out_seed_vector)
{
	if (!iscc_check_data_set(data_set)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid data set object.");
	}
	scc_ErrorCode ec;
	if ((ec = iscc_check_cluster_options(options, iscc_num_data_points(data_set))) != SCC_ER_OK) {
		return ec;
	}
	if (out_seed_vector->seeds != NULL) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "`out_seed_vector->seeds` must be NULL.");
	}
	return iscc_make_error(SCC_ER_NOT_IMPLEMENTED);
}


void scc_free_seed_vector(scc_SeedVector* const seed_vector)
{
	if (seed_vector != NULL) {
		free(seed_vector->seeds);
		*seed_vector = (scc_SeedVector) {
			.num_seeds = 0,
			.seeds = NULL,
		};
	}
}


// =============================================================================
// External function implementations
// =============================================================================

scc_ErrorCode iscc_check_cluster_options(const scc_ClusterOptions* const options,
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
		if (options->num_types > UINT16_MAX) {
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
			(options->seed_method != SCC_SM_BATCHES) &&
			(options->seed_method != SCC_SM_INWARDS_ORDER) &&
			(options->seed_method != SCC_SM_INWARDS_UPDATING) &&
			(options->seed_method != SCC_SM_EXCLUSION_ORDER) &&
			(options->seed_method != SCC_SM_EXCLUSION_UPDATING)) {
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
