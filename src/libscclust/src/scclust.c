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
#include <float.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
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


// =============================================================================
// External function implementations
// =============================================================================

void scc_get_compiled_version(uint32_t* const out_major,
                              uint32_t* const out_minor,
                              uint32_t* const out_patch)
{
	if (out_major != NULL) *out_major = SCC_SCCLUST_MAJOR_VERSION;
	if (out_minor != NULL) *out_minor = SCC_SCCLUST_MINOR_VERSION;
	if (out_patch != NULL) *out_patch = SCC_SCCLUST_PATCH_VERSION;
}


scc_ErrorCode scc_init_empty_clustering(const uintmax_t num_data_points,
                                        scc_Clabel external_cluster_labels[const],
                                        scc_Clustering** const out_clustering)
{
	if (out_clustering == NULL) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Output parameter may not be NULL.");
	}
	// Initialize to null, so subsequent functions detect invalid clustering
	// if user doesn't check for errors.
	*out_clustering = NULL;

	if (num_data_points == 0) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Clustering must have positive number of data points.");
	}
	if (num_data_points > ISCC_POINTINDEX_MAX) {
		return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many data points (adjust the 'scc_PointIndex' type).");
	}
	if (num_data_points > SIZE_MAX - 1) {
		return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many data points.");
	}

	scc_Clustering* tmp_cl = malloc(sizeof(scc_Clustering));
	if (tmp_cl == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);

	*tmp_cl = (scc_Clustering) {
		.clustering_version = ISCC_CLUSTERING_STRUCT_VERSION,
		.num_data_points = (size_t) num_data_points,
		.num_clusters = 0,
		.cluster_label = external_cluster_labels,
		.external_labels = (external_cluster_labels != NULL),
	};

	assert(iscc_check_input_clustering(tmp_cl));

	*out_clustering = tmp_cl;

	return iscc_no_error();
}


scc_ErrorCode scc_init_existing_clustering(const uintmax_t num_data_points,
                                           const uintmax_t num_clusters,
                                           scc_Clabel current_cluster_labels[const],
                                           const bool deep_label_copy,
                                           scc_Clustering** const out_clustering)
{
	if (out_clustering == NULL) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Output parameter may not be NULL.");
	}
	// Initialize to null, so subsequent functions detect invalid clustering
	// if user doesn't check for errors.
	*out_clustering = NULL;

	if (num_data_points == 0) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Clustering must have positive number of data points.");
	}
	if (num_data_points > ISCC_POINTINDEX_MAX) {
		return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many data points (adjust the `scc_PointIndex` type).");
	}
	if (num_data_points > SIZE_MAX - 1) {
		return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many data points.");
	}
	if (num_clusters == 0) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Empty clustering.");
	}
	if (num_clusters > ((uintmax_t) SCC_CLABEL_MAX)) {
		return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many clusters (adjust the `scc_Clabel` type).");
	}
	if (num_clusters > SIZE_MAX) {
		return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many clusters.");
	}
	if (current_cluster_labels == NULL) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid cluster labels.");
	}

	const size_t num_data_points_st = (size_t) num_data_points;

	scc_Clustering* tmp_cl = malloc(sizeof(scc_Clustering));
	if (tmp_cl == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);

	*tmp_cl = (scc_Clustering) {
		.clustering_version = ISCC_CLUSTERING_STRUCT_VERSION,
		.num_data_points = num_data_points_st,
		.num_clusters = (size_t) num_clusters,
		.cluster_label = NULL,
		.external_labels = !deep_label_copy,
	};

	if (deep_label_copy) {
		tmp_cl->cluster_label = malloc(sizeof(scc_Clabel[num_data_points_st]));
		if (tmp_cl->cluster_label == NULL) {
			free(tmp_cl);
			return iscc_make_error(SCC_ER_NO_MEMORY);
		}
		memcpy(tmp_cl->cluster_label, current_cluster_labels, num_data_points_st * sizeof(scc_Clabel));
	} else {
		tmp_cl->cluster_label = current_cluster_labels;
	}

	assert(iscc_check_input_clustering(tmp_cl));

	*out_clustering = tmp_cl;

	return iscc_no_error();
}


void scc_free_clustering(scc_Clustering** const clustering)
{
	if ((clustering != NULL) && (*clustering != NULL)) {
		if (!((*clustering)->external_labels)) free((*clustering)->cluster_label);
		free(*clustering);
		*clustering = NULL;
	}
}


bool scc_is_initialized_clustering(const scc_Clustering* const clustering)
{
	if (clustering == NULL) return false;
	if (clustering->clustering_version != ISCC_CLUSTERING_STRUCT_VERSION) return false;
	if (clustering->num_data_points == 0) return false;
	if (clustering->num_data_points > ISCC_POINTINDEX_MAX) return false;
	if (clustering->num_clusters > ((uintmax_t) SCC_CLABEL_MAX)) return false;
	if ((clustering->num_clusters > 0) && (clustering->cluster_label == NULL)) return false;

	return true;
}


scc_ErrorCode scc_copy_clustering(const scc_Clustering* const in_clustering,
                                  scc_Clustering** const out_clustering)
{
	if (out_clustering == NULL) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Output parameter may not be NULL.");
	}
	// Initialize to null, so subsequent functions detect invalid clustering
	// if user doesn't check for errors.
	*out_clustering = NULL;

	if (!iscc_check_input_clustering(in_clustering)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid clustering object.");
	}

	scc_Clustering* tmp_cl = malloc(sizeof(scc_Clustering));
	if (tmp_cl == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);

	*tmp_cl = (scc_Clustering) {
		.clustering_version = ISCC_CLUSTERING_STRUCT_VERSION,
		.num_data_points = in_clustering->num_data_points,
		.num_clusters = in_clustering->num_clusters,
		.cluster_label = NULL,
		.external_labels = false,
	};

	if (in_clustering->num_clusters > 0) {
		tmp_cl->cluster_label = malloc(sizeof(scc_Clabel[in_clustering->num_data_points]));
		if (tmp_cl->cluster_label == NULL) {
			free(tmp_cl);
			return iscc_make_error(SCC_ER_NO_MEMORY);
		}
		memcpy(tmp_cl->cluster_label, in_clustering->cluster_label, in_clustering->num_data_points * sizeof(scc_Clabel));
	}

	assert(iscc_check_input_clustering(tmp_cl));

	*out_clustering = tmp_cl;

	return iscc_no_error();
}


scc_ErrorCode scc_check_clustering(const scc_Clustering* const clustering,
                                   const uint32_t size_constraint,
                                   const uintmax_t num_types,
                                   const uint32_t type_constraints[const],
                                   const size_t len_type_labels,
                                   const scc_TypeLabel type_labels[const],
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

	if (num_types < 2) {
		if (type_constraints != NULL) {
			return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid type constraints.");
		}
		if (len_type_labels != 0) {
			return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid type labels.");
		}
		if (type_labels != NULL) {
			return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid type labels.");
		}
	} else {
		if (num_types > ISCC_TYPELABEL_MAX) {
			return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many data point types.");
		}
		if (num_types > UINT_FAST16_MAX) {
			return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many data point types.");
		}
		if (type_constraints == NULL) {
			return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid type constraints.");
		}
		if (len_type_labels < clustering->num_data_points) {
			return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid type labels.");
		}
		if (type_labels == NULL) {
			return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid type labels.");
		}
	}

	assert(clustering->num_clusters <= ((uintmax_t) SCC_CLABEL_MAX));
	const scc_Clabel max_cluster = (scc_Clabel) clustering->num_clusters;
	for (size_t i = 0; i < clustering->num_data_points; ++i) {
		if ((clustering->cluster_label[i] > 0) && (clustering->cluster_label[i] < max_cluster)) continue;
		if (clustering->cluster_label[i] == 0) continue; // Since `scc_Clabel` can be unsigned
		if (clustering->cluster_label[i] == SCC_CLABEL_NA) continue;
		return iscc_no_error(); // Error found, return. (`out_is_OK` is set to false)
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


scc_ErrorCode scc_get_clustering_info(const scc_Clustering* const clustering,
                                      uintmax_t* const out_num_data_points,
                                      uintmax_t* const out_num_clusters)
{
	if (!iscc_check_input_clustering(clustering)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid clustering object.");
	}

	if (out_num_data_points != NULL) *out_num_data_points = clustering->num_data_points;
	if (out_num_clusters != NULL) *out_num_clusters = clustering->num_clusters;

	return iscc_no_error();
}


scc_ErrorCode scc_get_cluster_labels(const scc_Clustering* const clustering,
                                     const size_t len_out_label_buffer,
                                     scc_Clabel out_label_buffer[const])
{
	if (!iscc_check_input_clustering(clustering)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid clustering object.");
	}
	if (clustering->num_clusters == 0) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Empty clustering.");
	}
	if (len_out_label_buffer == 0) {
		return iscc_make_error(SCC_ER_INVALID_INPUT);
	}
	if (out_label_buffer == NULL) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Output parameter may not be NULL.");
	}

	size_t write = 0;
	for (; (write < len_out_label_buffer) && (write < clustering->num_data_points); ++write) {
		out_label_buffer[write] = clustering->cluster_label[write];
	}
	for (; write < len_out_label_buffer; ++write) {
		out_label_buffer[write] = SCC_CLABEL_NA;
	}

	return iscc_no_error();
}


scc_ErrorCode scc_get_clustering_stats(const scc_Clustering* const clustering,
                                       void* const data_set,
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
	if (!iscc_check_data_set(data_set, clustering->num_data_points)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid data set object.");
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
		.min_cluster_size = UINTMAX_MAX,
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
