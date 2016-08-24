/* ==============================================================================
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
 * ============================================================================== */

#include "../include/scclust.h"

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include "dist_search.h"
#include "error.h"
#include "scclust_int.h"


// ==============================================================================
// Internal function prototypes
// ==============================================================================

scc_ErrorCode iscc_run_nng_batches(scc_Clustering* clustering,
                                   iscc_NNSearchObject* nn_search_object,
                                   uint32_t size_constraint,
                                   bool ignore_unassigned,
                                   bool main_radius_constraint,
                                   double main_radius,
                                   const bool main_data_points[],
                                   uint32_t batch_size,
                                   iscc_Dpid* batch_indices,
                                   iscc_Dpid* out_indices,
                                   bool* assigned);


// ==============================================================================
// External function implementations
// ==============================================================================

scc_ErrorCode scc_nng_clustering_batches(scc_Clustering* const clustering,
                                         void* const data_set_object,
                                         const uint32_t size_constraint,
                                         const scc_UnassignedMethod main_unassigned_method,
                                         const bool main_radius_constraint,
                                         const double main_radius,
                                         const size_t len_main_data_points,
                                         const bool main_data_points[const],
                                         uint32_t batch_size)
{
	if (!iscc_check_input_clustering(clustering)) return iscc_make_error(SCC_ER_INVALID_CLUSTERING);
	if (!iscc_check_data_set_object(data_set_object, clustering->num_data_points)) return iscc_make_error(SCC_ER_INVALID_DATA_OBJ);
	if (size_constraint < 2) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (clustering->num_data_points < size_constraint) return iscc_make_error(SCC_ER_NO_CLUST_EXIST_CONSTRAINT);
	if ((main_unassigned_method != SCC_UM_IGNORE) && (main_unassigned_method != SCC_UM_ASSIGN_BY_NNG)) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if (main_radius_constraint && (main_radius <= 0.0)) return iscc_make_error(SCC_ER_INVALID_INPUT);
	if ((main_data_points != NULL) && (len_main_data_points < clustering->num_data_points)) return iscc_make_error(SCC_ER_INVALID_INPUT);

	if (clustering->num_clusters != 0) return iscc_make_error(SCC_ER_NOT_IMPLEMENTED);

	if (batch_size == 0) batch_size = UINT32_MAX;
	if (batch_size > clustering->num_data_points) {
		batch_size = (uint32_t) clustering->num_data_points;
	}

	iscc_NNSearchObject* nn_search_object;
	if (!iscc_init_nn_search_object(data_set_object,
	                                clustering->num_data_points,
	                                NULL,
	                                &nn_search_object)) {
		return iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
	}

	iscc_Dpid* const batch_indices = malloc(sizeof(iscc_Dpid[batch_size]));
	iscc_Dpid* const out_indices = malloc(sizeof(iscc_Dpid[size_constraint * batch_size]));
	bool* const assigned = calloc(clustering->num_data_points, sizeof(bool));
	if ((batch_indices == NULL) || (out_indices == NULL) || (assigned == NULL)) {
		free(batch_indices);
		free(out_indices);
		free(assigned);
		iscc_close_nn_search_object(&nn_search_object);
		return iscc_make_error(SCC_ER_NO_MEMORY);
	}

	// Initialize cluster labels
	if (clustering->cluster_label == NULL) {
		clustering->external_labels = false;
		clustering->cluster_label = malloc(sizeof(scc_Clabel[clustering->num_data_points]));
		if (clustering->cluster_label == NULL) {
			free(batch_indices);
			free(out_indices);
			free(assigned);
			iscc_close_nn_search_object(&nn_search_object);
			return iscc_make_error(SCC_ER_NO_MEMORY);
		}
	}

	scc_ErrorCode ec = iscc_run_nng_batches(clustering,
	                                        nn_search_object,
	                                        size_constraint,
	                                        (main_unassigned_method == SCC_UM_IGNORE),
	                                        main_radius_constraint,
	                                        main_radius,
	                                        main_data_points,
	                                        batch_size,
	                                        batch_indices,
	                                        out_indices,
	                                        assigned);

	free(batch_indices);
	free(out_indices);
	free(assigned);
	iscc_close_nn_search_object(&nn_search_object);

	return ec;
}


// ==============================================================================
// Internal function implementations 
// ==============================================================================

scc_ErrorCode iscc_run_nng_batches(scc_Clustering* const clustering,
                                   iscc_NNSearchObject* const nn_search_object,
                                   const uint32_t size_constraint,
                                   const bool ignore_unassigned,
                                   const bool main_radius_constraint,
                                   const double main_radius,
                                   const bool main_data_points[const],
                                   const uint32_t batch_size,
                                   iscc_Dpid* const batch_indices,
                                   iscc_Dpid* const out_indices,
                                   bool* const assigned)
{
	assert(iscc_check_input_clustering(clustering));
	assert(clustering->cluster_label != NULL);
	assert(clustering->num_clusters == 0);
	assert(nn_search_object != NULL);
	assert(size_constraint >= 2);
	assert(clustering->num_data_points >= size_constraint);
	assert(!main_radius_constraint || (main_radius > 0.0));
	assert(batch_size > 0);
	assert(batch_indices != NULL);
	assert(out_indices != NULL);
	assert(assigned != NULL);

	bool search_done = false;
	scc_Clabel next_cluster_label = 0;
	assert(clustering->num_data_points <= ISCC_DPID_MAX);
	const iscc_Dpid num_data_points = (iscc_Dpid) clustering->num_data_points; // If `iscc_Dpid` is signed

	for (iscc_Dpid curr_point = 0; curr_point < num_data_points; ) {

		size_t in_batch = 0;
		if (main_data_points == NULL) {
			for (; (in_batch < batch_size) && (curr_point < num_data_points); ++curr_point) {
				if (!assigned[curr_point]) {
					batch_indices[in_batch] = curr_point;
					++in_batch;
				}
			}
		} else {
			for (; (in_batch < batch_size) && (curr_point < num_data_points); ++curr_point) {
				if (!assigned[curr_point]) {
					if (main_data_points[curr_point]) {
						batch_indices[in_batch] = curr_point;
						++in_batch;
					} else {
						assert(!assigned[curr_point]);
						clustering->cluster_label[curr_point] = SCC_CLABEL_NA;
					}
				}
			}
		}

		if (in_batch == 0) {
			assert(curr_point == num_data_points);
			break;
		}

		search_done = true;
		if (!iscc_nearest_neighbor_search_index(nn_search_object,
		                                        in_batch,
		                                        batch_indices,
		                                        size_constraint,
		                                        main_radius_constraint,
		                                        main_radius,
		                                        out_indices)) {
			return iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
		}

		const iscc_Dpid* check_indices = out_indices;
		for (size_t i = 0; i < in_batch; ++i) {
			const iscc_Dpid* const stop_check_indices = check_indices + size_constraint;
			if (!assigned[batch_indices[i]]) {
				if (stop_check_indices[-1] == ISCC_DPID_NA) {
					// Radius constraint binding
					assert(!assigned[batch_indices[i]]);
					clustering->cluster_label[batch_indices[i]] = SCC_CLABEL_NA;
				} else {
					for (; (check_indices != stop_check_indices) && !assigned[*check_indices]; ++check_indices) {}
					if (check_indices == stop_check_indices) {
						// `i` has no assigned neighbors and can be seed
						if (next_cluster_label == SCC_CLABEL_MAX) {
							return iscc_make_error(SCC_ER_TOO_LARGE_PROBLEM);
						}

						assert(!assigned[batch_indices[i]]);
						const iscc_Dpid* const stop_assign_indices = stop_check_indices - 1;
						for (check_indices -= size_constraint; check_indices != stop_assign_indices; ++check_indices) {
							assert(!assigned[*check_indices]);
							assigned[*check_indices] = true;
							clustering->cluster_label[*check_indices] = next_cluster_label;
						}
						if (assigned[batch_indices[i]]) {
							// Self-loop from `batch_indices[i]` to `batch_indices[i]` existed among NN
							assert(!assigned[*check_indices]);
							assigned[*check_indices] = true;
							clustering->cluster_label[*check_indices] = next_cluster_label;
						} else {
							// Self-loop did not exist
							assert(!assigned[batch_indices[i]]);
							assigned[batch_indices[i]] = true;
							clustering->cluster_label[batch_indices[i]] = next_cluster_label;
						}

						assert(clustering->cluster_label[batch_indices[i]] == next_cluster_label);
						++next_cluster_label;
					} else {
						// `i` has assigned neighbors and cannot be seed
						if (ignore_unassigned) {
							assert(!assigned[batch_indices[i]]);
							clustering->cluster_label[batch_indices[i]] = SCC_CLABEL_NA;
						} else {
							// Assign `batch_indices[i]` to a preliminary cluster.
							// If a future seed wants it as neighbor, it switches cluster.
							assert(assigned[*check_indices]);
							assert(clustering->cluster_label[*check_indices] != SCC_CLABEL_NA);
							assert(!assigned[batch_indices[i]]);
							clustering->cluster_label[batch_indices[i]] = clustering->cluster_label[*check_indices];
						}
					}
				}
			}
			check_indices = stop_check_indices;
		} // Loop in batch
	} // Loop between batches

	if (next_cluster_label == 0) {
		if (!search_done) {
			// Never did search, i.e., main_data_points are all false
			assert(main_data_points != NULL);
			return iscc_make_error(SCC_ER_NO_CLUST_EXIST_CONSTRAINT);
		} else {
			// Did search but still no clusters, i.e., too tight radius constraint
			assert(main_radius_constraint);
			return iscc_make_error(SCC_ER_NO_CLUST_EXIST_RADIUS);
		}
	}

	clustering->num_clusters = (size_t) next_cluster_label;

	return iscc_no_error();
}
