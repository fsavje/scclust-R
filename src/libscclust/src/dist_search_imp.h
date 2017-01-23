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

#ifndef SCC_DIST_SEARCH_IMP_HG
#define SCC_DIST_SEARCH_IMP_HG

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "../include/scclust.h"
#include "../include/scclust_spi.h"

#ifdef __cplusplus
extern "C" {
#endif


// =============================================================================
// Miscellaneous functions
// =============================================================================

bool iscc_imp_check_data_set(void* data_set,
                             size_t num_data_points);

// `output_dists` must be of length `(len_point_indices - 1) len_point_indices / 2`
bool iscc_imp_get_dist_matrix(void* data_set,
                              size_t len_point_indices,
                              const scc_PointIndex point_indices[],
                              double output_dists[]);

// `output_dists` must be of length `len_query_indices * len_column_indices`
bool iscc_imp_get_dist_rows(void* data_set,
                            size_t len_query_indices,
                            const scc_PointIndex query_indices[],
                            size_t len_column_indices,
                            const scc_PointIndex column_indices[],
                            double output_dists[]);


// =============================================================================
// Max dist functions
// =============================================================================

bool iscc_imp_init_max_dist_object(void* data_set,
                                   size_t len_search_indices,
                                   const scc_PointIndex search_indices[],
                                   iscc_MaxDistObject** out_max_dist_object);

// `max_indices` and `max_dists` must be of length `n_query_points`
bool iscc_imp_get_max_dist(iscc_MaxDistObject* max_dist_object,
                           size_t len_query_indices,
                           const scc_PointIndex query_indices[],
                           scc_PointIndex out_max_indices[],
                           double out_max_dists[]);

bool iscc_imp_close_max_dist_object(iscc_MaxDistObject** max_dist_object);


// =============================================================================
// Nearest neighbor search functions
// =============================================================================

bool iscc_imp_init_nn_search_object(void* data_set,
                                    size_t len_search_indices,
                                    const scc_PointIndex search_indices[],
                                    iscc_NNSearchObject** out_nn_search_object);

// `out_nn_indices` must be of length `k * len_query_indices`
bool iscc_imp_nearest_neighbor_search(iscc_NNSearchObject* nn_search_object,
                                      size_t len_query_indices,
                                      const scc_PointIndex query_indices[],
                                      uint32_t k,
                                      bool radius_search,
                                      double radius,
                                      size_t* out_num_ok_queries,
                                      scc_PointIndex out_query_indices[],
                                      scc_PointIndex out_nn_indices[]);

bool iscc_imp_close_nn_search_object(iscc_NNSearchObject** nn_search_object);


#ifdef __cplusplus
}
#endif

#endif // ifndef SCC_DIST_SEARCH_IMP_HG
