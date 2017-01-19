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
 *  The scclust library SPI...
 */

#ifndef SCC_SCCLUST_SPI_HG
#define SCC_SCCLUST_SPI_HG

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
// Internal types
// =============================================================================

/** Type used for data point IDs. May be unsigned or signed.
 *
 *  \note
 *  Number of data points in any clustering problem must be strictly less
 *  than the maximum number that can be stored in #scc_PointIndex.
 */
typedef int scc_PointIndex;

static const scc_PointIndex SCC_POINTINDEX_NA = INT_MAX;

/** Type used for arc indices. Must be unsigned.
 *
 *  \note
 *  Number of arcs in any digraph must be less or equal to
 *  the maximum number that can be stored in #iscc_ArcIndex.
 */
typedef uint32_t iscc_ArcIndex;

typedef struct iscc_MaxDistObject iscc_MaxDistObject;

typedef struct iscc_NNSearchObject iscc_NNSearchObject;

#define SCC_M_POINTINDEX_TYPE_int
#define SCC_M_POINTINDEX_NA INT_MAX
#define ISCC_M_ARCINDEX_TYPE_uint32_t


// =============================================================================
// Distance search functions
// =============================================================================

typedef bool (*scc_check_data_set) (void*,
                                    size_t);

typedef bool (*scc_get_dist_matrix) (void*,
                                     size_t,
                                     const scc_PointIndex*,
                                     double*);

typedef bool (*scc_get_dist_rows) (void*,
                                   size_t,
                                   const scc_PointIndex*,
                                   size_t,
                                   const scc_PointIndex*,
                                   double*);

typedef bool (*scc_init_max_dist_object) (void*,
                                          size_t,
                                          const scc_PointIndex*,
                                          iscc_MaxDistObject**);

typedef bool (*scc_get_max_dist) (iscc_MaxDistObject*,
                                  size_t,
                                  const scc_PointIndex*,
                                  scc_PointIndex*,
                                  double*);

typedef bool (*scc_close_max_dist_object) (iscc_MaxDistObject**);

typedef bool (*scc_init_nn_search_object) (void*,
                                           size_t,
                                           const scc_PointIndex*,
                                           iscc_NNSearchObject**);

typedef bool (*scc_nearest_neighbor_search_digraph) (iscc_NNSearchObject*,
                                                     size_t,
                                                     const bool*,
                                                     bool*,
                                                     uint32_t,
                                                     bool,
                                                     double,
                                                     iscc_ArcIndex*,
                                                     scc_PointIndex*);

typedef bool (*scc_nearest_neighbor_search_index) (iscc_NNSearchObject*,
                                                   size_t,
                                                   const scc_PointIndex*,
                                                   uint32_t,
                                                   bool,
                                                   double,
                                                   scc_PointIndex*);

typedef bool (*scc_close_nn_search_object) (iscc_NNSearchObject**);


// =============================================================================
// SPI functions
// =============================================================================

bool scc_reset_dist_functions(void);

bool scc_set_dist_functions(scc_check_data_set,
                            scc_get_dist_matrix,
                            scc_get_dist_rows,
                            scc_init_max_dist_object,
                            scc_get_max_dist,
                            scc_close_max_dist_object,
                            scc_init_nn_search_object,
                            scc_nearest_neighbor_search_digraph,
                            scc_nearest_neighbor_search_index,
                            scc_close_nn_search_object);


#ifdef __cplusplus
}
#endif

#endif // ifndef SCC_SCCLUST_HG
