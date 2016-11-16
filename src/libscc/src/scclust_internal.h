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

/** @file
 *
 *  Header with internal type definitions.
 */

#ifndef SCC_SCCLUST_INTERNAL_HG
#define SCC_SCCLUST_INTERNAL_HG

#ifdef __cplusplus
// So g++ defines integer limits
#define __STDC_LIMIT_MACROS
#endif

#include <limits.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "../include/scclust.h"

#ifdef __cplusplus
	extern "C" {
#endif


// =============================================================================
// Macros
// =============================================================================

/// Macro for internal clustering checks.
#define iscc_check_input_clustering(cl) scc_is_initialized_clustering(cl)


// =============================================================================
// Structs, types and variables
// =============================================================================

/** Type used for data point IDs. May be unsigned or signed.
 *
 *  \note
 *  Number of data points in any clustering problem must be strictly less
 *  than the maximum number that can be stored in #iscc_Dpid. I.e.,
 *  possible data point IDs must be in the sequence `[0, 1, ..., ISCC_DPID_MAX - 1]`.
 *  Independent of `ISCC_DPID_MAX`, number of data points may not be greater than `SIZE_MAX - 1`.
 */
#if defined(SCC_DPID_INT)
	typedef int iscc_Dpid;
	#define ISCC_DPID_MAX_MACRO INT_MAX
	static const uintmax_t ISCC_DPID_MAX = INT_MAX;
	static const iscc_Dpid ISCC_DPID_NA = INT_MAX;
#elif defined(SCC_DPID_UINT64)
	typedef uint64_t iscc_Dpid;
	#define ISCC_DPID_MAX_MACRO UINT64_MAX
	static const uintmax_t ISCC_DPID_MAX = UINT64_MAX;
	static const iscc_Dpid ISCC_DPID_NA = UINT64_MAX;
#else
	typedef uint32_t iscc_Dpid;
	#define ISCC_DPID_MAX_MACRO UINT32_MAX
	static const uintmax_t ISCC_DPID_MAX = UINT32_MAX;
	static const iscc_Dpid ISCC_DPID_NA = UINT32_MAX;
#endif


/** Type used for arc indices. Must be unsigned.
 *
 *  \note
 *  Number of arcs in any digraph must be less or equal to
 *  the maximum number that can be stored in #iscc_Arci.
 */
#ifdef SCC_ARC64
	typedef uint64_t iscc_Arci;
	#define ISCC_ARCI_MAX_MACRO UINT64_MAX
	static const uintmax_t ISCC_ARCI_MAX = UINT64_MAX;
#else
	typedef uint32_t iscc_Arci;
	#define ISCC_ARCI_MAX_MACRO UINT32_MAX
	static const uintmax_t ISCC_ARCI_MAX = UINT32_MAX;
#endif


/// Label given to unassigned vertices.
static const uintmax_t ISCC_TYPELABEL_MAX = 65535;


/** Clustering struct.
 *
 *  This struct describes clusterings by enumerating a cluster label for each vertex.
 */
struct scc_Clustering {
	/// Version of the struct.
	int32_t clustering_version;

	/// Number of data points in the clustering problem.
	size_t num_data_points;

	/// Number of clusters.
	size_t num_clusters;

	/** Array of length #num_data_points with cluster labels of the assigned data points.
	 *  Unassigned vertices have the value #SCC_CLABEL_NA.
	 */
	scc_Clabel* cluster_label;

	/** Indicator whether the #cluster_label array was assigned by the user (rather than the library). If #external_labels is \c true,
	 *  #cluster_label will not be freed when the instance of #scc_Clustering is destroyed.
	 */
	bool external_labels;
};

/// Current version of the clustering struct.
static const int32_t ISCC_CURRENT_CLUSTSTRUCT_VERSION = 1;

/** The null clustering.
 *
 *  The null clustering is an easily detectable invalid clustering. It is mainly used as return
 *  value when functions encounter errors.
 */
static const scc_Clustering ISCC_NULL_CLUSTERING = { 0, 0, 0, NULL, false };


#ifdef __cplusplus
	}
#endif

#endif // ifndef SCC_SCCLUST_INTERNAL_HG
