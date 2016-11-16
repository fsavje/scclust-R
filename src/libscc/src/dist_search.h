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

#ifndef SCC_DIST_SEARCH_HG
#define SCC_DIST_SEARCH_HG

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "scclust_internal.h"

#ifdef __cplusplus
extern "C" {
#endif


// =============================================================================
// Miscellaneous functions
// =============================================================================

bool iscc_check_data_set_object(void* data_set_object,
                                size_t required_data_points);


// `output_dists` must be of length `(len_point_indices - 1) len_point_indices / 2`
bool iscc_get_dist_matrix(void* data_set_object,
                          size_t len_point_indices,
                          const iscc_Dpid point_indices[],
                          double output_dists[]);

// `output_dists` must be of length `len_query_indices * len_column_indices`
bool iscc_get_dist_rows(void* data_set_object,
                        size_t len_query_indices,
                        const iscc_Dpid query_indices[],
                        size_t len_column_indices,
                        const iscc_Dpid column_indices[],
                        double output_dists[]);


// =============================================================================
// Max dist functions
// =============================================================================

typedef struct iscc_MaxDistObject iscc_MaxDistObject;

bool iscc_init_max_dist_object(void* data_set_object,
                               size_t len_search_indices,
                               const iscc_Dpid search_indices[],
                               iscc_MaxDistObject** out_max_dist_object);

// `max_indices` and `max_dists` must be of length `n_query_points`
bool iscc_get_max_dist(iscc_MaxDistObject* max_dist_object,
                       size_t len_query_indices,
                       const iscc_Dpid query_indices[],
                       iscc_Dpid out_max_indices[],
                       double out_max_dists[]);

bool iscc_close_max_dist_object(iscc_MaxDistObject** max_dist_object);


// =============================================================================
// Nearest neighbor search functions
// =============================================================================

typedef struct iscc_NNSearchObject iscc_NNSearchObject;

bool iscc_init_nn_search_object(void* data_set_object,
                                size_t len_search_indices,
                                const iscc_Dpid search_indices[],
                                iscc_NNSearchObject** out_nn_search_object);

// `out_nn_ref` must be of length `len_query_indicators + 1`
// `out_nn_indices` must be of length `k * sum(num_query_indicators)`
bool iscc_nearest_neighbor_search_digraph(iscc_NNSearchObject* nn_search_object,
                                          size_t len_query_indicators,
                                          const bool query_indicators[],
                                          bool out_query_indicators[],
                                          uint32_t k,
                                          bool radius_search,
                                          double radius,
                                          bool accept_partial,
                                          iscc_Arci out_nn_ref[],
                                          iscc_Dpid out_nn_indices[]);

// `out_nn_indices` must be of length `k * len_query_indices`
bool iscc_nearest_neighbor_search_index(iscc_NNSearchObject* nn_search_object,
                                        size_t len_query_indices,
                                        const iscc_Dpid query_indices[],
                                        uint32_t k,
                                        bool radius_search,
                                        double radius,
                                        iscc_Dpid out_nn_indices[]);

bool iscc_close_nn_search_object(iscc_NNSearchObject** nn_search_object);


#ifdef __cplusplus
}
#endif

#endif // ifndef SCC_DIST_SEARCH_HG
