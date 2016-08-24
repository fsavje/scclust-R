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

/** @file
 *
 *  The scclust library...
 */

#ifndef SCC_SCCLUST_HG
#define SCC_SCCLUST_HG

#ifdef __cplusplus
// So g++ defines integer limits
#define __STDC_LIMIT_MACROS
#endif

#include <limits.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif


// ==============================================================================
// Library specific types, user-serviceable
// ==============================================================================

/** Type used for cluster labels. May be unsigned or signed.
 *  
 *  \note
 *  Number of clusters in any clustering problem must be strictly less
 *  than the maximum number that can be stored in #scc_Clabel. I.e., 
 *  cluster labels must be in the sequence `[0, 1, ..., SCC_CLABEL_MAX - 1]`, 
 *  and `SCC_CLABEL_NA` may not be in this sequence (but it may be `SCC_CLABEL_MAX`).
 *
 *  \note
 *  The library has been tested with #scc_Clabel set to `uint32_t`, `uint64_t` and `int`.
 */
typedef int scc_Clabel;

/// Maximum number that can be stored in #scc_Clabel. May not be greater than `SIZE_MAX`.
static const scc_Clabel SCC_CLABEL_MAX = INT_MAX;

/// Label given to unassigned vertices.
static const scc_Clabel SCC_CLABEL_NA = INT_MIN;

/** Type used to indicate data point type (for the NNG method). May be unsigned or signed.
 *  
 *  \note
 *  Type labels must be in the sequence `[0, 1, ..., 65534]`.
 *
 *  \note
 *  The library has been tested with #scc_TypeLabel set to `uint_fast16_t` and `int`.
 */
typedef int scc_TypeLabel;


// ==============================================================================
// Version information
// ==============================================================================

#define SCC_SCCLUST_MAJOR_VERSION 0
#define SCC_SCCLUST_MINOR_VERSION 1
#define SCC_SCCLUST_PATCH_VERSION 1
#define SCC_CHECK_VERSION(major, minor) ((major == SCC_SCCLUST_MAJOR_VERSION) && (minor <= SCC_SCCLUST_MINOR_VERSION))

void scc_get_compiled_version(uint32_t* out_major,
                              uint32_t* out_minor,
                              uint32_t* out_patch);


// ==============================================================================
// Error handling
// ==============================================================================

enum scc_ErrorCode {
	SCC_ER_OK,
	SCC_ER_NULL_INPUT,
	SCC_ER_INVALID_INPUT,
	SCC_ER_INVALID_CLUSTERING,
	SCC_ER_EMPTY_CLUSTERING,
	SCC_ER_INVALID_DATA_OBJ,
	SCC_ER_NO_MEMORY,
	SCC_ER_TOO_LARGE_PROBLEM,
	SCC_ER_TOO_LARGE_DIGRAPH,
	SCC_ER_DIST_SEARCH_ERROR,
	SCC_ER_NO_CLUST_EXIST_CONSTRAINT,
	SCC_ER_NO_CLUST_EXIST_RADIUS,
	SCC_ER_NOT_IMPLEMENTED,
};

/// Typedef for the scc_ErrorCode enum
typedef enum scc_ErrorCode scc_ErrorCode;

bool scc_get_latest_error(size_t len_error_message_buffer,
                          char error_message_buffer[]);


// ==============================================================================
// Utility functions
// ==============================================================================

/// Type used for clusterings
typedef struct scc_Clustering scc_Clustering;

/// Type used for clustering statistics
typedef struct scc_ClusteringStats scc_ClusteringStats;

/// Struct to store clustering statistics
struct scc_ClusteringStats {
	uintmax_t num_data_points;
	uintmax_t num_assigned;
	uintmax_t num_clusters;
	uintmax_t num_populated_clusters;
	uintmax_t min_cluster_size;
	uintmax_t max_cluster_size;
	double avg_cluster_size;
	double sum_dists;
	double min_dist;
	double max_dist;
	double cl_avg_min_dist;
	double cl_avg_max_dist;
	double cl_avg_dist_weighted;
	double cl_avg_dist_unweighted;
};


void scc_free_clustering(scc_Clustering** clustering);

scc_ErrorCode scc_init_empty_clustering(uintmax_t num_data_points,
                                        scc_Clabel external_cluster_labels[],
                                        scc_Clustering** out_clustering);

scc_ErrorCode scc_init_existing_clustering(uintmax_t num_data_points,
                                           uintmax_t num_clusters,
                                           scc_Clabel current_cluster_labels[],
                                           bool deep_label_copy,
                                           scc_Clustering** out_clustering);

scc_ErrorCode scc_copy_clustering(const scc_Clustering* in_clustering,
                                  scc_Clustering** out_clustering);

bool scc_is_initialized_clustering(const scc_Clustering* clustering);

scc_ErrorCode scc_check_clustering(const scc_Clustering* clustering,
                                   uint32_t size_constraint,
                                   bool* out_is_OK);

scc_ErrorCode scc_check_clustering_types(const scc_Clustering* clustering,
                                         uint32_t size_constraint,
                                         uintmax_t num_types,
                                         const uint32_t type_size_constraints[],
                                         size_t len_type_labels,
                                         const scc_TypeLabel type_labels[],
                                         bool* out_is_OK);

scc_ErrorCode scc_get_clustering_info(const scc_Clustering* clustering,
                                      uintmax_t* out_num_data_points,
                                      uintmax_t* out_num_clusters);

scc_ErrorCode scc_get_cluster_labels(const scc_Clustering* clustering,
                                     size_t len_out_label_buffer,
                                     scc_Clabel out_label_buffer[]);

scc_ErrorCode scc_get_clustering_stats(const scc_Clustering* clustering,
                                       void* data_set_object,
                                       scc_ClusteringStats* out_stats);


// ==============================================================================
// NNG-based clustering
// ==============================================================================

/** Enum to specify seed finding methods.
 *
 *  The NNG-based clustering algorithms find seeds to build the clustering on. This enum specifies which method is used to find the seeds.
 *
 *  In most settings, it is desired to find as many clusters (i.e., as many seeds) as possible. The choice between the methods is
 *  largely one between performance in this aspect and resource requirements.
 *
 *  See #scc_get_seed_clustering for more details. See also the appendix of \cite Higgins2016.
 */
enum scc_SeedMethod {

	/** Find seeds lexically by vertex ID. 
	 *
	 *  This method finds seed sequentially by checking whether adding the next seed satisfy the four conditions described in #scc_get_seed_clustering.
	 */
	SCC_SM_LEXICAL,

	/** Find seeds ordered by inwards pointing arcs. 
	 *
	 *  This method counts vertices' inwards pointing arcs and finds seeds in ascending order by the arc count. Vertices pointing to a seed cannot
	 *  themselves be seeds, thus a vertex with many inwards pointing arcs will exclude many vertices from being seeds. Heuristically, picking such
	 *  a vertex as seed will lead to fewer clusters. 
	 */
	SCC_SM_INWARDS_ORDER,

	/** Find seeds ordered by inwards pointing arcs from unassigned vertices. 
	 *
	 *  This method counts vertices' inwards pointing arcs from *unassigned* vertices and finds seeds by in ascending order by the arc count. Unlike
	 *  the #SCC_SM_INWARDS_ORDER, this method updates the arc count after finding a seed so that only arcs where the tails are unassigned are counted.
	 *  Vertices already assigned to a cluster cannot be a seed, thus it is inconsequential whether they are pointing to a seed.
	 *
	 *  If the desired size is two, this method ensures that the maximum distance between any two vertices in a common cluster in the
	 *  final clustering is bounded by twice the maximum distance in the NNG.
	 */
	SCC_SM_INWARDS_UPDATING,

	/** Find seeds ordered by edge count in the exclusion graph.
	 *
	 *  The exclusion graph is the undirected graph where an edge is drawn between two vertices if they cannot both be seeds.
	 *  Any maximal independent set in this graph is a valid set of seeds. This method counts the edges of each vertex in this
	 *  graph and find seeds in ascending order.
	 */
	SCC_SM_EXCLUSION_ORDER,

	/** Find seeds ordered by edge count in the exclusion graph from non-excluded vertices.
	 *
	 *  The exclusion graph is the undirected graph where an edge is draw between two vertices if they cannot both be seeds.
	 *  Any maximal independent set in this graph is a valid set of seeds. This method counts the edges of each that vertex is not already excluded
	 *  and find seeds in ascending order by this count. Unlike the #SCC_SM_EXCLUSION_ORDER, this method updates the edge count after finding a
	 *  seed so that only edges where the tails that still can become seeds are counted.
	 */
	SCC_SM_EXCLUSION_UPDATING,
};

/// Typedef for the scc_NNGMethod enum
typedef enum scc_SeedMethod scc_SeedMethod;

enum scc_UnassignedMethod {
	SCC_UM_IGNORE,
	SCC_UM_ASSIGN_BY_NNG,
	SCC_UM_CLOSEST_ASSIGNED,
	SCC_UM_CLOSEST_SEED,
	SCC_UM_CLOSEST_SEED_EST_RADIUS,
};

typedef enum scc_UnassignedMethod scc_UnassignedMethod;

scc_ErrorCode scc_nng_clustering(scc_Clustering* clustering,
                                 void* data_set_object,
                                 uint32_t size_constraint,
                                 scc_SeedMethod seed_method,
                                 scc_UnassignedMethod main_unassigned_method,
                                 bool main_radius_constraint,
                                 double main_radius,
                                 size_t len_main_data_points,
                                 const bool main_data_points[],
                                 scc_UnassignedMethod secondary_unassigned_method,
                                 bool secondary_radius_constraint,
                                 double secondary_radius);

scc_ErrorCode scc_nng_clustering_batches(scc_Clustering* clustering,
                                         void* data_set_object,
                                         uint32_t size_constraint,
                                         scc_UnassignedMethod main_unassigned_method,
                                         bool main_radius_constraint,
                                         double main_radius,
                                         size_t len_main_data_points,
                                         const bool main_data_points[],
                                         uint32_t batch_size);

scc_ErrorCode scc_nng_clustering_with_types(scc_Clustering* clustering,
                                            void* data_set_object,
                                            uint32_t size_constraint,
                                            uintmax_t num_types,
                                            const uint32_t type_size_constraints[],
                                            size_t len_type_labels,
                                            const scc_TypeLabel type_labels[],
                                            scc_SeedMethod seed_method,
                                            scc_UnassignedMethod main_unassigned_method,
                                            bool main_radius_constraint,
                                            double main_radius,
                                            size_t len_main_data_points,
                                            const bool main_data_points[],
                                            scc_UnassignedMethod secondary_unassigned_method,
                                            bool secondary_radius_constraint,
                                            double secondary_radius);


// ==============================================================================
// Greedy clustering function
// ==============================================================================

scc_ErrorCode scc_top_down_greedy_clustering(scc_Clustering* clustering,
                                             void* data_set_object,
                                             uint32_t size_constraint,
                                             bool batch_assign);

scc_ErrorCode scc_bottom_up_greedy_clustering(scc_Clustering* clustering,
                                              void* data_set_object,
                                              uint32_t size_constraint);


#ifdef __cplusplus
}
#endif

#endif // ifndef SCC_SCCLUST_HG
