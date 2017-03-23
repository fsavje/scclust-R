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


// =============================================================================
// Version information
// =============================================================================

#define SCC_SCCLUST_MAJOR_VERSION 0
#define SCC_SCCLUST_MINOR_VERSION 1
#define SCC_SCCLUST_PATCH_VERSION 1


/** Get the version the library was compiled with.
 *
 *  \param[out] out_major variable to write major version to.
 *  \param[out] out_minor variable to write minor version to.
 *  \param[out] out_patch variable to write patch version to.
 *
 */
void scc_get_compiled_version(uint32_t* out_major,
                              uint32_t* out_minor,
                              uint32_t* out_patch);


// =============================================================================
// Error handling
// =============================================================================

/// Error code returned by methods in the library.
typedef enum scc_ErrorCode {
	/// No error.
	SCC_ER_OK,

	/// Unkonwn error.
	SCC_ER_UNKNOWN_ERROR,

	/// Function parameters are invalid.
	SCC_ER_INVALID_INPUT,

	/// Cannot allocate required memory.
	SCC_ER_NO_MEMORY,

	/// Clustering problem has no solution.
	SCC_ER_NO_SOLUTION,

	/// Clustering problem is too large.
	SCC_ER_TOO_LARGE_PROBLEM,

	/// Failed to calculate distances.
	SCC_ER_DIST_SEARCH_ERROR,

	/// Functionality not yet implemented.
	SCC_ER_NOT_IMPLEMENTED

} scc_ErrorCode;


/** Get latest scclust error.
 *
 *  Writes a description of the latest error prroduced by the library to the
 *  supplied buffer.
 *
 *  \param[in] len_error_message_buffer the length of the buffer #error_message_buffer.
 *  \param[out] error_message_buffer the buffer to write to.
 *
 *  \return \c true if write was a success, otherwise \c false.
 */
bool scc_get_latest_error(size_t len_error_message_buffer,
                          char error_message_buffer[]);


// =============================================================================
// Library types
// =============================================================================

/** Type used for data point IDs.
 *
 *  \note
 *  Number of data points in any clustering problem must be strictly less
 *  than the maximum number that can be stored in #scc_PointIndex.
 */
typedef int scc_PointIndex;


/// Macro for data point ID type.
#define SCC_M_POINTINDEX_TYPE_int


/** Type used for cluster labels.
 *
 *  \note
 *  Number of clusters in any clustering problem must be strictly less
 *  than the maximum number that can be stored in #scc_Clabel.
 */
typedef int scc_Clabel;


/// Macro for cluster labels type.
#define SCC_M_CLABEL_TYPE_int


/// Label given to unassigned vertices.
static const scc_Clabel SCC_CLABEL_NA = INT_MIN;


/// Macro for unassigned label.
#define SCC_M_CLABEL_NA INT_MIN


/** Type used to indicate data point type.
 *
 *  \note
 *  Type labels must be in the sequence `[0, 1, ..., 65534]`.
 */
typedef int scc_TypeLabel;


/// Macro for cluster data point type.
#define SCC_M_TYPELABEL_TYPE_int


// =============================================================================
// Data set object
// =============================================================================

/// Typedef for struct containing data sets.
typedef struct scc_DataSet scc_DataSet;


/** Construct new data set from raw data.
 *
 *  Creates a #scc_DataSet based on supplied raw data.
 *
 *  \param[in] num_data_points the number of data points in the data set.
 *  \param[in] num_dimensions the number of dimensions for each data point.
 *  \param[in] len_data_matrix the length of #data_matrix.
 *  \param[in] data_matrix the raw data. The data should be ordered first by point,
 *                         then by dimension. With three units (A, B, C) in two
 *                         dimensions, #data_matrix should be `[A_1, A_2, B_1,
 *                         B_2, C_1, C_2]`.
 *  \param[out] out_data_set double pointer to where to write the data set reference.
 *
 *  \return #scc_ErrorCode describing eventual error.
 */
scc_ErrorCode scc_init_data_set(uint64_t num_data_points,
                                uint32_t num_dimensions,
                                size_t len_data_matrix,
                                const double data_matrix[],
                                scc_DataSet** out_data_set);


/** Free data set.
 *
 *  Frees a #scc_DataSet previously allocated by #scc_init_data_set.
 *
 *  \param[in,out] data_set double pointer to a #scc_DataSet objec to free.
 */
void scc_free_data_set(scc_DataSet** data_set);


/** Check data set.
 *
 *  Checks whether inputted data set is properly initialized.
 *
 *  \param[in] data_set pointer to a #scc_DataSet objec to check.
 *
 *  \return \c true if #data_set is initialized, otherwise \c false.
 *
 *  \note This function does not check whether #data_set is a sensible data_set.
 *        It only ensures that #data_set is initialized properly and is of a
 *        compatible version.
 */
bool scc_is_initialized_data_set(const scc_DataSet* data_set);


// =============================================================================
// Clustering object
// =============================================================================

/// Type used for clusterings
typedef struct scc_Clustering scc_Clustering;


scc_ErrorCode scc_init_empty_clustering(uint64_t num_data_points,
                                        scc_Clabel external_cluster_labels[],
                                        scc_Clustering** out_clustering);


scc_ErrorCode scc_init_existing_clustering(uint64_t num_data_points,
                                           uint64_t num_clusters,
                                           scc_Clabel current_cluster_labels[],
                                           bool deep_label_copy,
                                           scc_Clustering** out_clustering);


void scc_free_clustering(scc_Clustering** clustering);


bool scc_is_initialized_clustering(const scc_Clustering* clustering);


scc_ErrorCode scc_copy_clustering(const scc_Clustering* in_clustering,
                                  scc_Clustering** out_clustering);


scc_ErrorCode scc_get_clustering_info(const scc_Clustering* clustering,
                                      uint64_t* out_num_data_points,
                                      uint64_t* out_num_clusters);


scc_ErrorCode scc_get_cluster_labels(const scc_Clustering* clustering,
                                     size_t len_out_label_buffer,
                                     scc_Clabel out_label_buffer[]);


// =============================================================================
// Clustering functions
// =============================================================================

/** Enum to specify seed finding methods.
 *
 *  The NNG-based clustering algorithms find seeds to build the clustering on. This enum specifies which method is used to find the seeds.
 *
 *  In most settings, it is desired to find as many clusters (i.e., as many seeds) as possible. The choice between the methods is
 *  largely one between performance in this aspect and resource requirements.
 *
 *  See #scc_get_seed_clustering for more details. See also the appendix of \cite Higgins2016.
 */
typedef enum scc_SeedMethod {
	/** Find seeds lexically by vertex ID.
	 *
	 *  This method finds seed sequentially by checking whether adding the next seed satisfy the four conditions described in #scc_get_seed_clustering.
	 */
	SCC_SM_LEXICAL,


	SCC_SM_BATCHES,

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
	SCC_SM_EXCLUSION_UPDATING

} scc_SeedMethod;


typedef enum scc_UnassignedMethod {
	SCC_UM_IGNORE,
	SCC_UM_ANY_NEIGHBOR,
	SCC_UM_CLOSEST_ASSIGNED,
	SCC_UM_CLOSEST_SEED
} scc_UnassignedMethod;


typedef enum scc_RadiusMethod {
	SCC_RM_NO_RADIUS,
	SCC_RM_USE_SUPPLIED,
	SCC_RM_USE_SEED_RADIUS,
	SCC_RM_USE_ESTIMATED
} scc_RadiusMethod;


typedef struct scc_ClusterOptions {
	/** scc_ClusterOptions struct version
	 *
	 *  \note
	 *  This must be set to "722678001".
	 */
	int32_t options_version;
	uint32_t size_constraint;
	uint32_t num_types;
	const uint32_t* type_constraints;
	size_t len_type_labels;
	const scc_TypeLabel* type_labels;
	scc_SeedMethod seed_method;
	size_t len_primary_data_points;
	const scc_PointIndex* primary_data_points;
	scc_UnassignedMethod primary_unassigned_method;
	scc_UnassignedMethod secondary_unassigned_method;
	scc_RadiusMethod seed_radius;
	double seed_supplied_radius;
	scc_RadiusMethod primary_radius;
	double primary_supplied_radius;
	scc_RadiusMethod secondary_radius;
	double secondary_supplied_radius;
	uint32_t batch_size;
} scc_ClusterOptions;


scc_ClusterOptions scc_get_default_options(void);


scc_ErrorCode scc_sc_clustering(void* data_set,
                                const scc_ClusterOptions* options,
                                scc_Clustering* out_clustering);


scc_ErrorCode scc_hierarchical_clustering(void* data_set,
                                          uint32_t size_constraint,
                                          bool batch_assign,
                                          scc_Clustering* out_clustering);


// =============================================================================
// Utility functions
// =============================================================================

/** Check clustering
 *
 *  Checks whether supplied clustering satisfies clustering constraints in options.
 *  Currently it checks `size_constraint`, `type_constraints` and `primary_data_points`.
 *  I.e., it ignores all radius constraints.
 *
 *  \return #scc_ErrorCode describing eventual error.
 */
scc_ErrorCode scc_check_clustering(const scc_Clustering* clustering,
                                   const scc_ClusterOptions* options,
                                   bool* out_is_OK);


/// Struct to report clustering statistics
typedef struct scc_ClusteringStats {
	uint64_t num_data_points;
	uint64_t num_assigned;
	uint64_t num_clusters;
	uint64_t num_populated_clusters;
	uint64_t min_cluster_size;
	uint64_t max_cluster_size;
	double avg_cluster_size;
	double sum_dists;
	double min_dist;
	double max_dist;
	double avg_min_dist;
	double avg_max_dist;
	double avg_dist_weighted;
	double avg_dist_unweighted;
} scc_ClusteringStats;


scc_ErrorCode scc_get_clustering_stats(void* data_set,
                                       const scc_Clustering* clustering,
                                       scc_ClusteringStats* out_stats);


#ifdef __cplusplus
}
#endif

#endif // ifndef SCC_SCCLUST_HG
