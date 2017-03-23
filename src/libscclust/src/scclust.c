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
#include <string.h>
#include "clustering_struct.h"
#include "error.h"
#include "scclust_types.h"


// =============================================================================
// Public function implementations
// =============================================================================

void scc_get_compiled_version(uint32_t* const out_major,
                              uint32_t* const out_minor,
                              uint32_t* const out_patch)
{
	if (out_major != NULL) *out_major = SCC_SCCLUST_MAJOR_VERSION;
	if (out_minor != NULL) *out_minor = SCC_SCCLUST_MINOR_VERSION;
	if (out_patch != NULL) *out_patch = SCC_SCCLUST_PATCH_VERSION;
}


scc_ErrorCode scc_init_empty_clustering(const uint64_t num_data_points,
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


scc_ErrorCode scc_init_existing_clustering(const uint64_t num_data_points,
                                           const uint64_t num_clusters,
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
	if (num_clusters > ((uint64_t) SCC_CLABEL_MAX)) {
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


scc_ErrorCode scc_get_clustering_info(const scc_Clustering* const clustering,
                                      uint64_t* const out_num_data_points,
                                      uint64_t* const out_num_clusters)
{
	if (!iscc_check_input_clustering(clustering)) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Invalid clustering object.");
	}

	if (out_num_data_points != NULL) *out_num_data_points = (uint64_t) clustering->num_data_points;
	if (out_num_clusters != NULL) *out_num_clusters = (uint64_t) clustering->num_clusters;

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
