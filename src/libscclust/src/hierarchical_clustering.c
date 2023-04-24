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
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include "dist_search.h"
#include "clustering_struct.h"
#include "error.h"
#include "scclust_types.h"

// Maximum number of data points to check when finding centers.
static const uint_fast16_t ISCC_HI_NUM_TO_CHECK = 100;


// =============================================================================
// Internal structs
// =============================================================================

typedef struct iscc_hi_DistanceEdge iscc_hi_DistanceEdge;
struct iscc_hi_DistanceEdge {
	scc_PointIndex head;
	double distance;
	iscc_hi_DistanceEdge* next_dist;
};


typedef struct iscc_hi_ClusterItem {
	size_t size;
	uint_fast16_t marker;
	scc_PointIndex* members;
} iscc_hi_ClusterItem;


typedef struct iscc_hi_ClusterStack {
	size_t capacity;
	size_t items;
	iscc_hi_ClusterItem* clusters;
	scc_PointIndex* pointindex_store;
} iscc_hi_ClusterStack;


typedef struct iscc_hi_WorkArea {
	scc_PointIndex* const pointindex_array1;
	scc_PointIndex* const pointindex_array2;
	double* const dist_array;
	uint_fast16_t* const vertex_markers;
	iscc_hi_DistanceEdge* const edge_store1;
	iscc_hi_DistanceEdge* const edge_store2;
} iscc_hi_WorkArea;


// =============================================================================
// Static function prototypes
// =============================================================================

static scc_ErrorCode iscc_hi_empty_cl_stack(size_t num_data_points,
                                            iscc_hi_ClusterStack* out_cl_stack);


static scc_ErrorCode iscc_hi_init_cl_stack(const scc_Clustering* in_cl,
                                           iscc_hi_ClusterStack* out_cl_stack,
                                           size_t* out_size_largest_cluster);


static scc_ErrorCode iscc_hi_run_hierarchical_clustering(iscc_hi_ClusterStack* cl_stack,
                                                         scc_Clustering* cl,
                                                         void* data_set,
                                                         iscc_hi_WorkArea* work_area,
                                                         uint32_t size_constraint,
                                                         bool batch_assign);


static scc_ErrorCode iscc_hi_check_capacity(iscc_hi_ClusterStack* cl_stack);


static inline void iscc_hi_push_to_stack(iscc_hi_ClusterStack* cl_stack,
                                         iscc_hi_ClusterItem** cl);


static scc_ErrorCode iscc_hi_break_cluster_into_two(iscc_hi_ClusterItem* cluster_to_break,
                                                    void* data_set,
                                                    iscc_hi_WorkArea* work_area,
                                                    uint32_t size_constraint,
                                                    bool batch_assign,
                                                    iscc_hi_ClusterItem* out_new_cluster);


static inline uint_fast16_t iscc_hi_get_next_marker(iscc_hi_ClusterItem* cl,
                                                    uint_fast16_t vertex_markers[]);


static inline iscc_hi_DistanceEdge* iscc_hi_get_next_k_nn(iscc_hi_DistanceEdge* prev_dist,
                                                          uint32_t k,
                                                          const uint_fast16_t vertex_markers[],
                                                          uint_fast16_t curr_marker,
                                                          scc_PointIndex out_dist_array[static k]);


static inline iscc_hi_DistanceEdge* iscc_hi_get_next_dist(iscc_hi_DistanceEdge* prev_dist,
                                                          const uint_fast16_t vertex_markers[],
                                                          uint_fast16_t curr_marker);


static inline void iscc_hi_move_point_to_cluster1(scc_PointIndex id,
                                                  iscc_hi_ClusterItem* cl,
                                                  uint_fast16_t vertex_markers[],
                                                  uint_fast16_t curr_marker);


static inline void iscc_hi_move_point_to_cluster2(scc_PointIndex id,
                                                  iscc_hi_ClusterItem* cl,
                                                  uint_fast16_t vertex_markers[],
                                                  uint_fast16_t curr_marker);


static inline void iscc_hi_move_array_to_cluster1(uint32_t len_ids,
                                                  const scc_PointIndex ids[static len_ids],
                                                  iscc_hi_ClusterItem* cl,
                                                  uint_fast16_t vertex_markers[],
                                                  uint_fast16_t curr_marker);


static inline void iscc_hi_move_array_to_cluster2(uint32_t len_ids,
                                                  const scc_PointIndex ids[static len_ids],
                                                  iscc_hi_ClusterItem* cl,
                                                  uint_fast16_t vertex_markers[],
                                                  uint_fast16_t curr_marker);


static scc_ErrorCode iscc_hi_find_centers(iscc_hi_ClusterItem* cl,
                                          void* data_set,
                                          iscc_hi_WorkArea* work_area,
                                          scc_PointIndex* out_center1,
                                          scc_PointIndex* out_center2);


static scc_ErrorCode iscc_hi_populate_edge_lists(const iscc_hi_ClusterItem* cl,
                                                 void* data_set,
                                                 scc_PointIndex center1,
                                                 scc_PointIndex center2,
                                                 iscc_hi_WorkArea* work_area);


static inline void iscc_hi_sort_edge_list(const iscc_hi_ClusterItem* cl,
                                          scc_PointIndex center,
                                          const double row_dists[static cl->size],
                                          iscc_hi_DistanceEdge edge_store[static cl->size]);


static int iscc_hi_compare_dist_edges(const void* a,
                                      const void* b);


// =============================================================================
// Public function implementations
// =============================================================================

scc_ErrorCode scc_hierarchical_clustering(void* const data_set,
                                          const uint32_t size_constraint,
                                          const bool batch_assign,
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
	if (size_constraint < 2) {
		return iscc_make_error_msg(SCC_ER_INVALID_INPUT, "Size constraint must be 2 or greater.");
	}
	if (out_clustering->num_data_points < size_constraint) {
		return iscc_make_error_msg(SCC_ER_NO_SOLUTION, "Fewer data points than size constraint.");
	}

	scc_ErrorCode ec;
	size_t size_largest_cluster = 0; // Initialize to avoid gcc warning
	iscc_hi_ClusterStack cl_stack;
	if (out_clustering->num_clusters == 0) {
		if (out_clustering->cluster_label == NULL) {
			out_clustering->external_labels = false;
			out_clustering->cluster_label = malloc(sizeof(scc_Clabel[out_clustering->num_data_points]));
			if (out_clustering->cluster_label == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);
		}

		size_largest_cluster = out_clustering->num_data_points;
		if ((ec = iscc_hi_empty_cl_stack(out_clustering->num_data_points, &cl_stack)) != SCC_ER_OK) {
			return ec;
		}
	} else {
		if ((ec = iscc_hi_init_cl_stack(out_clustering, &cl_stack, &size_largest_cluster)) != SCC_ER_OK) {
			return ec;
		}
	}

	assert(cl_stack.items > 0);
	assert(cl_stack.clusters != NULL);
	assert(cl_stack.pointindex_store != NULL);

	const size_t size_pointindex_array = (size_constraint > ISCC_HI_NUM_TO_CHECK) ? size_constraint : ISCC_HI_NUM_TO_CHECK;
	const size_t size_dist_array = ((2 * size_largest_cluster) > ISCC_HI_NUM_TO_CHECK) ? (2 * size_largest_cluster) : ISCC_HI_NUM_TO_CHECK;
	iscc_hi_WorkArea work_area = {
		.pointindex_array1 = malloc(sizeof(scc_PointIndex[size_pointindex_array])),
		.pointindex_array2 = malloc(sizeof(scc_PointIndex[size_pointindex_array])),
		.dist_array = malloc(sizeof(double[size_dist_array])),
		.vertex_markers = calloc(out_clustering->num_data_points, sizeof(uint_fast16_t)),
		.edge_store1 = malloc(sizeof(iscc_hi_DistanceEdge[size_largest_cluster])),
		.edge_store2 = malloc(sizeof(iscc_hi_DistanceEdge[size_largest_cluster])),
	};

	if ((work_area.pointindex_array1 == NULL) || (work_area.pointindex_array2 == NULL) ||
	        (work_area.dist_array == NULL) || (work_area.vertex_markers == NULL) ||
	        (work_area.edge_store1 == NULL) || (work_area.edge_store2 == NULL)) {
		ec = iscc_make_error(SCC_ER_NO_MEMORY);
	}

	if (ec == SCC_ER_OK) {
		ec = iscc_hi_run_hierarchical_clustering(&cl_stack,
		                                         out_clustering,
		                                         data_set,
		                                         &work_area,
		                                         size_constraint,
		                                         batch_assign);
	}

	free(work_area.pointindex_array1);
	free(work_area.pointindex_array2);
	free(work_area.dist_array);
	free(work_area.vertex_markers);
	free(work_area.edge_store1);
	free(work_area.edge_store2);
	free(cl_stack.clusters);
	free(cl_stack.pointindex_store);

	return ec;
}


// =============================================================================
// Static function implementations
// =============================================================================

static scc_ErrorCode iscc_hi_empty_cl_stack(const size_t num_data_points,
                                            iscc_hi_ClusterStack* const out_cl_stack)
{
	assert(num_data_points >= 2);
	assert(out_cl_stack != NULL);

	const size_t tmp_capacity = 1 + ((size_t) (20 * log2((double) num_data_points)));
	*out_cl_stack = (iscc_hi_ClusterStack) {
		.capacity = tmp_capacity,
		.items = 1,
		.clusters = malloc(sizeof(iscc_hi_ClusterItem[tmp_capacity])),
		.pointindex_store = malloc(sizeof(scc_PointIndex[num_data_points])),
	};
	if ((out_cl_stack->clusters == NULL) || (out_cl_stack->pointindex_store == NULL)) {
		free(out_cl_stack->clusters);
		free(out_cl_stack->pointindex_store);
		return iscc_make_error(SCC_ER_NO_MEMORY);
	}

	assert(num_data_points <= ISCC_POINTINDEX_MAX);
	const scc_PointIndex num_data_points_pi = (scc_PointIndex) num_data_points; // If `scc_PointIndex` is signed
	for (scc_PointIndex i = 0; i < num_data_points_pi; ++i) {
		out_cl_stack->pointindex_store[i] = i;
	}

	out_cl_stack->clusters[0] = (iscc_hi_ClusterItem) {
		.size = num_data_points,
		.marker = 0,
		.members = out_cl_stack->pointindex_store,
	};

	return iscc_no_error();
}


static scc_ErrorCode iscc_hi_init_cl_stack(const scc_Clustering* const in_cl,
                                           iscc_hi_ClusterStack* const out_cl_stack,
                                           size_t* const out_size_largest_cluster)
{
	assert(iscc_check_input_clustering(in_cl));
	assert(in_cl->num_data_points >= 2);
	assert(out_cl_stack != NULL);
	assert(out_size_largest_cluster != NULL);

	const uint64_t tmp_capacity = in_cl->num_clusters + 1 + ((uint64_t) (10 * log2((double) in_cl->num_data_points)));
	if (tmp_capacity > SIZE_MAX) {
		return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many clusters.");
	}
	*out_cl_stack = (iscc_hi_ClusterStack) {
		.capacity = (size_t) tmp_capacity,
		.items = in_cl->num_clusters,
		.clusters = calloc((size_t) tmp_capacity, sizeof(iscc_hi_ClusterItem)),
		.pointindex_store = malloc(sizeof(scc_PointIndex[in_cl->num_data_points])),
	};
	if ((out_cl_stack->clusters == NULL) || (out_cl_stack->pointindex_store == NULL)) {
		free(out_cl_stack->clusters);
		free(out_cl_stack->pointindex_store);
		return iscc_make_error(SCC_ER_NO_MEMORY);
	}

	const scc_Clabel* const cluster_label = in_cl->cluster_label;
	iscc_hi_ClusterItem* const clusters = out_cl_stack->clusters;

	for (size_t i = 0; i < in_cl->num_data_points; ++i) {
		if (cluster_label[i] != SCC_CLABEL_NA) {
			++(clusters[cluster_label[i]].size);
		}
	}

	size_t size_largest_cluster = 0;
	clusters[0].members = out_cl_stack->pointindex_store + clusters[0].size;
	for (size_t c = 1; c < in_cl->num_clusters; ++c) {
		clusters[c].members = clusters[c - 1].members + clusters[c].size;
		size_largest_cluster = (clusters[c].size > size_largest_cluster) ? clusters[c].size : size_largest_cluster;
	}
	*out_size_largest_cluster = size_largest_cluster;

	assert(in_cl->num_data_points <= ISCC_POINTINDEX_MAX);
	const scc_PointIndex num_data_points = (scc_PointIndex) in_cl->num_data_points; // If `scc_PointIndex` is signed
	for (scc_PointIndex i = 0; i < num_data_points; ++i) {
		if (cluster_label[i] != SCC_CLABEL_NA) {
			--clusters[cluster_label[i]].members;
			*(clusters[cluster_label[i]].members) = i;
		}
	}

	return iscc_no_error();
}


static scc_ErrorCode iscc_hi_run_hierarchical_clustering(iscc_hi_ClusterStack* const cl_stack,
                                                         scc_Clustering* const cl,
                                                         void* const data_set,
                                                         iscc_hi_WorkArea* const work_area,
                                                         const uint32_t size_constraint,
                                                         const bool batch_assign)
{
	assert(cl_stack != NULL);
	assert(cl_stack->items > 0);
	assert(cl_stack->items <= cl_stack->capacity);
	assert(cl_stack->clusters != NULL);
	assert(cl_stack->pointindex_store != NULL);
	assert(iscc_check_input_clustering(cl));
	assert(iscc_check_data_set(data_set));
	assert(iscc_num_data_points(data_set) == cl->num_data_points);
	assert(work_area != NULL);
	assert(size_constraint >= 2);

	scc_ErrorCode ec;
	scc_Clabel current_label = 0;
	while (cl_stack->items > 0) {

		if ((ec = iscc_hi_check_capacity(cl_stack)) != SCC_ER_OK) {
			return ec;
		}

		iscc_hi_ClusterItem* current_cluster = &cl_stack->clusters[cl_stack->items - 1];

		if (current_cluster->size < (2 * size_constraint)) {
			if (current_cluster->size > 0) {
				if (current_label == SCC_CLABEL_MAX) {
					return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many clusters (adjust the `scc_Clabel` type).");
				}
				for (size_t v = 0; v < current_cluster->size; ++v) {
					cl->cluster_label[current_cluster->members[v]] = current_label;
				}
				++current_label;
			}
			--(cl_stack->items);
		} else {
			iscc_hi_ClusterItem* new_cluster = NULL; // Initialize to avoid gcc warning
			iscc_hi_push_to_stack(cl_stack, &new_cluster);
			if ((ec = iscc_hi_break_cluster_into_two(current_cluster,
			                                         data_set,
			                                         work_area,
			                                         size_constraint,
			                                         batch_assign,
			                                         new_cluster)) != SCC_ER_OK) {
				return ec;
			}
		}
	}

	cl->num_clusters = (size_t) current_label;

	assert(cl_stack->items == 0);

	return iscc_no_error();
}


static scc_ErrorCode iscc_hi_check_capacity(iscc_hi_ClusterStack* const cl_stack)
{
	assert(cl_stack != NULL);
	assert(cl_stack->items <= cl_stack->capacity);
	assert(cl_stack->clusters != NULL);

	if (cl_stack->items == cl_stack->capacity) {
		const uintmax_t capacity_tmp = cl_stack->capacity + 16 + (cl_stack->capacity >> 4);
		if ((capacity_tmp > SIZE_MAX) || (capacity_tmp < cl_stack->capacity)) {
			return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many clusters.");
		}
		iscc_hi_ClusterItem* const clusters_tmp = realloc(cl_stack->clusters, sizeof(iscc_hi_ClusterItem[(size_t) capacity_tmp]));
		if (clusters_tmp == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);
		cl_stack->clusters = clusters_tmp;
		cl_stack->capacity = (size_t) capacity_tmp;
	}

	return iscc_no_error();
}


static inline void iscc_hi_push_to_stack(iscc_hi_ClusterStack* const cl_stack,
                                         iscc_hi_ClusterItem** const cl)
{
	assert(cl_stack != NULL);
	assert(cl_stack->items <= cl_stack->capacity);
	assert(cl_stack->clusters != NULL);
	assert(cl != NULL);

	*cl = &cl_stack->clusters[cl_stack->items];
	++(cl_stack->items);
}


static scc_ErrorCode iscc_hi_break_cluster_into_two(iscc_hi_ClusterItem* const cluster_to_break,
                                                    void* const data_set,
                                                    iscc_hi_WorkArea* const work_area,
                                                    const uint32_t size_constraint,
                                                    const bool batch_assign,
                                                    iscc_hi_ClusterItem* const out_new_cluster)
{
	assert(cluster_to_break != NULL);
	assert(cluster_to_break->size >= 2 * size_constraint);
	assert(cluster_to_break->members != NULL);
	assert(iscc_check_data_set(data_set));
	assert(work_area != NULL);
	assert(work_area->pointindex_array1 != NULL);
	assert(work_area->pointindex_array2 != NULL);
	assert(work_area->vertex_markers != NULL);
	assert(work_area->edge_store1 != NULL);
	assert(work_area->edge_store2 != NULL);
	assert(size_constraint >= 2);
	assert(out_new_cluster != NULL);

	scc_ErrorCode ec;
	scc_PointIndex center1 = ISCC_POINTINDEX_MAX_PI, center2 = ISCC_POINTINDEX_MAX_PI; // Initialize these to avoid gcc warning
	// `iscc_hi_find_centers` must be before `iscc_hi_get_next_marker`
	// since the marker becomes invalid after `iscc_hi_find_centers`
	if ((ec = iscc_hi_find_centers(cluster_to_break,
	                               data_set,
	                               work_area,
	                               &center1,
	                               &center2)) != SCC_ER_OK) {
		return ec;
	}

	if ((ec = iscc_hi_populate_edge_lists(cluster_to_break,
	                                      data_set,
	                                      center1,
	                                      center2,
	                                      work_area)) != SCC_ER_OK) {
		return ec;
	}

	scc_PointIndex* const k_nn_array1 = work_area->pointindex_array1;
	scc_PointIndex* const k_nn_array2 = work_area->pointindex_array2;
	uint_fast16_t* const vertex_markers = work_area->vertex_markers;

	// `edge_store1` and `edge_store2` have been populated by `iscc_hi_populate_edge_lists`
	iscc_hi_DistanceEdge* last_assigned_edge1 = work_area->edge_store1;
	iscc_hi_DistanceEdge* last_assigned_edge2 = work_area->edge_store2;

	iscc_hi_DistanceEdge* temp_edge1;
	iscc_hi_DistanceEdge* temp_edge2;

	size_t num_unassigned = cluster_to_break->size;
	const uint_fast16_t curr_marker = iscc_hi_get_next_marker(cluster_to_break, vertex_markers);

	*out_new_cluster = (iscc_hi_ClusterItem) {
		.size = 0,
		.marker = curr_marker,
		.members = cluster_to_break->members + num_unassigned,
	};

	iscc_hi_ClusterItem* const cluster1 = cluster_to_break;
	cluster1->size = 0;
	iscc_hi_ClusterItem* const cluster2 = out_new_cluster;

	iscc_hi_move_point_to_cluster1(center1, cluster1, vertex_markers, curr_marker);
	iscc_hi_move_point_to_cluster2(center2, cluster2, vertex_markers, curr_marker);

	temp_edge1 = iscc_hi_get_next_k_nn(last_assigned_edge1, size_constraint - 1, vertex_markers, curr_marker, k_nn_array1);
	temp_edge2 = iscc_hi_get_next_k_nn(last_assigned_edge2, size_constraint - 1, vertex_markers, curr_marker, k_nn_array2);

	if (temp_edge1->distance >= temp_edge2->distance) {
		iscc_hi_move_array_to_cluster1(size_constraint - 1, k_nn_array1, cluster1, vertex_markers, curr_marker);
		last_assigned_edge1 = temp_edge1;

		last_assigned_edge2 = iscc_hi_get_next_k_nn(last_assigned_edge2, size_constraint - 1, vertex_markers, curr_marker, k_nn_array2);
		iscc_hi_move_array_to_cluster2(size_constraint - 1, k_nn_array2, cluster2, vertex_markers, curr_marker);
	} else {
		iscc_hi_move_array_to_cluster2(size_constraint - 1, k_nn_array2, cluster2, vertex_markers, curr_marker);
		last_assigned_edge2 = temp_edge2;

		last_assigned_edge1 = iscc_hi_get_next_k_nn(last_assigned_edge1, size_constraint - 1, vertex_markers, curr_marker, k_nn_array1);
		iscc_hi_move_array_to_cluster1(size_constraint - 1, k_nn_array1, cluster1, vertex_markers, curr_marker);
	}

	assert((cluster1->size == size_constraint) && (cluster2->size == size_constraint));

	num_unassigned -= 2 * size_constraint;

	if (batch_assign) {
		uint32_t num_assign_in_batch = size_constraint;
		for (; num_unassigned > 0; num_unassigned -= num_assign_in_batch) {

			if (num_assign_in_batch > num_unassigned) num_assign_in_batch = (uint32_t) num_unassigned;

			temp_edge1 = iscc_hi_get_next_k_nn(last_assigned_edge1, num_assign_in_batch, vertex_markers, curr_marker, k_nn_array1);
			temp_edge2 = iscc_hi_get_next_k_nn(last_assigned_edge2, num_assign_in_batch, vertex_markers, curr_marker, k_nn_array2);

			if (temp_edge1->distance <= temp_edge2->distance) {
				iscc_hi_move_array_to_cluster1(num_assign_in_batch, k_nn_array1, cluster1, vertex_markers, curr_marker);
				last_assigned_edge1 = temp_edge1;
			} else {
				iscc_hi_move_array_to_cluster2(num_assign_in_batch, k_nn_array2, cluster2, vertex_markers, curr_marker);
				last_assigned_edge2 = temp_edge2;
			}
		}

	} else {
		for (; num_unassigned > 0; --num_unassigned) {
			temp_edge1 = iscc_hi_get_next_dist(last_assigned_edge1, vertex_markers, curr_marker);
			temp_edge2 = iscc_hi_get_next_dist(last_assigned_edge2, vertex_markers, curr_marker);

			if (temp_edge1->distance <= temp_edge2->distance) {
				iscc_hi_move_point_to_cluster1(temp_edge1->head, cluster1, vertex_markers, curr_marker);
				last_assigned_edge1 = temp_edge1;
			} else {
				iscc_hi_move_point_to_cluster2(temp_edge2->head, cluster2, vertex_markers, curr_marker);
				last_assigned_edge2 = temp_edge2;
			}
		}
	}

	assert(num_unassigned == 0);
	assert(cluster1->size >= size_constraint);
	assert(cluster2->size >= size_constraint);
	assert(cluster1->members + cluster1->size == cluster2->members);

	return iscc_no_error();
}


static inline uint_fast16_t iscc_hi_get_next_marker(iscc_hi_ClusterItem* const cl,
                                                    uint_fast16_t vertex_markers[const])
{
	assert(cl != NULL);
	assert(cl->size > 0);
	assert(cl->members != NULL);
	assert(vertex_markers != NULL);

	if (cl->marker == UINT_FAST16_MAX) {
		cl->marker = 0;
		for (size_t i = 0; i < cl->size; ++i) {
			vertex_markers[cl->members[i]] = 0;
		}
	}

	++(cl->marker);
	return cl->marker;
}


static inline iscc_hi_DistanceEdge* iscc_hi_get_next_k_nn(iscc_hi_DistanceEdge* prev_dist,
                                                          const uint32_t k,
                                                          const uint_fast16_t vertex_markers[const],
                                                          const uint_fast16_t curr_marker,
                                                          scc_PointIndex out_dist_array[const static k])
{
	assert(prev_dist != NULL);
	assert(prev_dist->next_dist != NULL); // We should never reach the end!
	assert(k > 0);
	assert(vertex_markers != NULL);
	assert(out_dist_array != NULL);

	for (uint32_t found = 0; found < k; ++found) {
		prev_dist = iscc_hi_get_next_dist(prev_dist, vertex_markers, curr_marker);
		out_dist_array[found] = prev_dist->head;
	}

	return prev_dist;
}


static inline iscc_hi_DistanceEdge* iscc_hi_get_next_dist(iscc_hi_DistanceEdge* const prev_dist,
                                                          const uint_fast16_t vertex_markers[const],
                                                          const uint_fast16_t curr_marker)
{
	assert(prev_dist != NULL);
	assert(prev_dist->next_dist != NULL); // We should never reach the end!
	assert(vertex_markers != NULL);

	while(vertex_markers[prev_dist->next_dist->head] == curr_marker) {
		// Vertex has already been assigned to a new cluster, skip it
		prev_dist->next_dist = prev_dist->next_dist->next_dist;
		assert(prev_dist->next_dist != NULL); // We should never reach the end!
	}

	return prev_dist->next_dist;
}


static inline void iscc_hi_move_point_to_cluster1(const scc_PointIndex id,
                                                  iscc_hi_ClusterItem* const cl,
                                                  uint_fast16_t vertex_markers[const],
                                                  const uint_fast16_t curr_marker)
{
	assert(cl != NULL);
	assert(cl->members != NULL);
	assert(vertex_markers != NULL);
	assert(vertex_markers[id] != curr_marker);
	assert(cl->marker == curr_marker);

	vertex_markers[id] = curr_marker;
	cl->members[cl->size] = id;
	cl->size += 1;
}


static inline void iscc_hi_move_point_to_cluster2(const scc_PointIndex id,
                                                  iscc_hi_ClusterItem* const cl,
                                                  uint_fast16_t vertex_markers[const],
                                                  const uint_fast16_t curr_marker)
{
	assert(cl != NULL);
	assert(cl->members != NULL);
	assert(vertex_markers != NULL);
	assert(vertex_markers[id] != curr_marker);
	assert(cl->marker == curr_marker);

	--(cl->members);
	vertex_markers[id] = curr_marker;
	*cl->members = id;
	cl->size += 1;
}


static inline void iscc_hi_move_array_to_cluster1(const uint32_t len_ids,
                                                  const scc_PointIndex ids[static len_ids],
                                                  iscc_hi_ClusterItem* const cl,
                                                  uint_fast16_t vertex_markers[const],
                                                  const uint_fast16_t curr_marker)
{
	assert(len_ids > 0);
	assert(ids != NULL);
	assert(cl != NULL);
	assert(cl->members != NULL);
	assert(vertex_markers != NULL);
	assert(cl->marker == curr_marker);

	const scc_PointIndex* const ids_stop = ids + len_ids;
	scc_PointIndex* cl_mem = cl->members + cl->size;
	for (; ids != ids_stop; ++ids, ++cl_mem) {
		assert(vertex_markers[*ids] != curr_marker);
		vertex_markers[*ids] = curr_marker;
		*cl_mem = *ids;
	}

	cl->size += len_ids;
}


static inline void iscc_hi_move_array_to_cluster2(const uint32_t len_ids,
                                                  const scc_PointIndex ids[static len_ids],
                                                  iscc_hi_ClusterItem* const cl,
                                                  uint_fast16_t vertex_markers[const],
                                                  const uint_fast16_t curr_marker)
{
	assert(len_ids > 0);
	assert(ids != NULL);
	assert(cl != NULL);
	assert(cl->members != NULL);
	assert(vertex_markers != NULL);
	assert(cl->marker == curr_marker);

	const scc_PointIndex* const ids_stop = ids + len_ids;
	for (; ids != ids_stop; ++ids) {
		assert(vertex_markers[*ids] != curr_marker);
		--(cl->members);
		vertex_markers[*ids] = curr_marker;
		*cl->members = *ids;
	}

	cl->size += len_ids;
}


static scc_ErrorCode iscc_hi_find_centers(iscc_hi_ClusterItem* const cl,
                                          void* const data_set,
                                          iscc_hi_WorkArea* const work_area,
                                          scc_PointIndex* const out_center1,
                                          scc_PointIndex* const out_center2)
{
	assert(cl != NULL);
	assert(cl->size >= 4);
	assert(cl->members != NULL);
	assert(iscc_check_data_set(data_set));
	assert(work_area != NULL);
	assert(work_area->pointindex_array1 != NULL);
	assert(work_area->pointindex_array2 != NULL);
	assert(work_area->dist_array != NULL);
	assert(work_area->vertex_markers != NULL);
	assert(out_center1 != NULL);
	assert(out_center2 != NULL);

	scc_PointIndex* const to_check = work_area->pointindex_array1;
	scc_PointIndex* const max_indices = work_area->pointindex_array2;
	double* const max_dists = work_area->dist_array;
	uint_fast16_t* const vertex_markers = work_area->vertex_markers;

	const uint_fast16_t curr_marker = iscc_hi_get_next_marker(cl, vertex_markers);

	size_t step = cl->size / ISCC_HI_NUM_TO_CHECK;
	if (step < 2) step = 2;
	// num_to_check = ceil(size / step) = floor((size + step - 1) / step) = 1 + floor((size - 1) / step)
	uint_fast16_t num_to_check = ((uint_fast16_t) (1 + (cl->size - 1) / step));
	num_to_check = (ISCC_HI_NUM_TO_CHECK < num_to_check) ? ISCC_HI_NUM_TO_CHECK : num_to_check;
	assert(num_to_check <= ISCC_HI_NUM_TO_CHECK);

	for (size_t i = 0; i < num_to_check; ++i) {
		to_check[i] = cl->members[i * step];
		vertex_markers[to_check[i]] = curr_marker;
	}

	iscc_MaxDistObject* max_dist_object;
	if (!iscc_init_max_dist_object(data_set, cl->size, cl->members, &max_dist_object)) {
		return iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
	}

	double max_dist = -1.0;
	while (num_to_check > 0) {
		if (!iscc_get_max_dist(max_dist_object, num_to_check, to_check, max_indices, max_dists)) {
			iscc_close_max_dist_object(&max_dist_object);
			return iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
		}

		uint_fast16_t write_in_to_check = 0;
		for (uint_fast16_t i = 0; i < num_to_check; ++i) {
			if (max_dists[i] > max_dist) {
				max_dist = max_dists[i];
				*out_center1 = to_check[i];
				*out_center2 = max_indices[i];
			}

			if (vertex_markers[max_indices[i]] != curr_marker) {
				vertex_markers[max_indices[i]] = curr_marker;
				to_check[write_in_to_check] = max_indices[i];
				++write_in_to_check;
			}
		}

		num_to_check = write_in_to_check;
	}

	if (!iscc_close_max_dist_object(&max_dist_object)) {
		return iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
	}

	return iscc_no_error();
}


static scc_ErrorCode iscc_hi_populate_edge_lists(const iscc_hi_ClusterItem* const cl,
                                                 void* const data_set,
                                                 const scc_PointIndex center1,
                                                 const scc_PointIndex center2,
                                                 iscc_hi_WorkArea* const work_area)
{
	assert(cl != NULL);
	assert(cl->size >= 4);
	assert(cl->members != NULL);
	assert(center1 != center2);
	assert(iscc_check_data_set(data_set));
	assert(work_area != NULL);
	assert(work_area->dist_array != NULL);
	assert(work_area->edge_store1 != NULL);
	assert(work_area->edge_store2 != NULL);

	double* const row_dists = work_area->dist_array;
	const scc_PointIndex query_indices[2] = { center1, center2 };

	if (!iscc_get_dist_rows(data_set,
	                        2,
	                        query_indices,
	                        cl->size,
	                        cl->members,
	                        row_dists)) {
		return iscc_make_error(SCC_ER_DIST_SEARCH_ERROR);
	}

	iscc_hi_sort_edge_list(cl, center1, row_dists, work_area->edge_store1);
	iscc_hi_sort_edge_list(cl, center2, row_dists + cl->size, work_area->edge_store2);

	return iscc_no_error();
}


static inline void iscc_hi_sort_edge_list(const iscc_hi_ClusterItem* const cl,
                                          const scc_PointIndex center,
                                          const double row_dists[const static cl->size],
                                          iscc_hi_DistanceEdge edge_store[const static cl->size])
{
	assert(cl != NULL);
	assert(cl->size >= 4);
	assert(cl->members != NULL);
	assert(row_dists != NULL);
	assert(edge_store != NULL);

	iscc_hi_DistanceEdge* write_edge = edge_store + 1;
	for (size_t i = 0; i < cl->size; ++i) {
		if (cl->members[i] == center) continue;
		write_edge->head = cl->members[i];
		write_edge->distance = row_dists[i];
		++write_edge;
	}

	assert(write_edge == (edge_store + cl->size));

	qsort(edge_store + 1, cl->size - 1, sizeof(iscc_hi_DistanceEdge), iscc_hi_compare_dist_edges);

	iscc_hi_DistanceEdge* const edge_stop = edge_store + cl->size - 1;
	for (iscc_hi_DistanceEdge* edge = edge_store; edge != edge_stop; ++edge) {
		edge->next_dist = edge + 1;
	}
	edge_stop->next_dist = NULL;
}


static int iscc_hi_compare_dist_edges(const void* const a,
                                      const void* const b)
{
	const double dist_a = ((const iscc_hi_DistanceEdge*)a)->distance;
	const double dist_b = ((const iscc_hi_DistanceEdge*)b)->distance;

	if (dist_a < dist_b) return -1;
	if (dist_a > dist_b) return 1;
	return 0;
}
