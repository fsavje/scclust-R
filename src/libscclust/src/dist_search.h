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
#include "../include/scclust_spi.h"


struct iscc_dist_functions_struct {
	scc_check_data_set check_data_set;
	scc_get_dist_matrix get_dist_matrix;
	scc_get_dist_rows get_dist_rows;
	scc_init_max_dist_object init_max_dist_object;
	scc_get_max_dist get_max_dist;
	scc_close_max_dist_object close_max_dist_object;
	scc_init_nn_search_object init_nn_search_object;
	scc_nearest_neighbor_search nearest_neighbor_search;
	scc_close_nn_search_object close_nn_search_object;
};

typedef struct iscc_dist_functions_struct iscc_dist_functions_struct;

extern iscc_dist_functions_struct iscc_dist_functions;


// =============================================================================
// Miscellaneous functions
// =============================================================================

static inline bool iscc_check_data_set(void* data_set,
                                       size_t num_data_points)
{
	return iscc_dist_functions.check_data_set(data_set,
	                                          num_data_points);
}


static inline bool iscc_get_dist_matrix(void* data_set,
                                        size_t len_point_indices,
                                        const scc_PointIndex point_indices[],
                                        double output_dists[])
{
	return iscc_dist_functions.get_dist_matrix(data_set,
	                                           len_point_indices,
	                                           point_indices,
	                                           output_dists);
}


static inline bool iscc_get_dist_rows(void* data_set,
                                      size_t len_query_indices,
                                      const scc_PointIndex query_indices[],
                                      size_t len_column_indices,
                                      const scc_PointIndex column_indices[],
                                      double output_dists[])
{
	return iscc_dist_functions.get_dist_rows(data_set,
	                                         len_query_indices,
	                                         query_indices,
	                                         len_column_indices,
	                                         column_indices,
	                                         output_dists);
}


// =============================================================================
// Max dist functions
// =============================================================================

static inline bool iscc_init_max_dist_object(void* data_set,
                                             size_t len_search_indices,
                                             const scc_PointIndex search_indices[],
                                             iscc_MaxDistObject** out_max_dist_object)
{
	return iscc_dist_functions.init_max_dist_object(data_set,
	                                                len_search_indices,
	                                                search_indices,
	                                                out_max_dist_object);
}


static inline bool iscc_get_max_dist(iscc_MaxDistObject* max_dist_object,
                                     size_t len_query_indices,
                                     const scc_PointIndex query_indices[],
                                     scc_PointIndex out_max_indices[],
                                     double out_max_dists[])
{
	return iscc_dist_functions.get_max_dist(max_dist_object,
	                                        len_query_indices,
	                                        query_indices,
	                                        out_max_indices,
	                                        out_max_dists);
}


static inline bool iscc_close_max_dist_object(iscc_MaxDistObject** max_dist_object)
{
	return iscc_dist_functions.close_max_dist_object(max_dist_object);
}


// =============================================================================
// Nearest neighbor search functions
// =============================================================================

static inline bool iscc_init_nn_search_object(void* data_set,
                                              size_t len_search_indices,
                                              const scc_PointIndex search_indices[],
                                              iscc_NNSearchObject** out_nn_search_object)
{
	return iscc_dist_functions.init_nn_search_object(data_set,
	                                                 len_search_indices,
	                                                 search_indices,
	                                                 out_nn_search_object);
}


static inline bool iscc_nearest_neighbor_search(iscc_NNSearchObject* nn_search_object,
                                                size_t len_query_indices,
                                                const scc_PointIndex query_indices[],
                                                uint32_t k,
                                                bool radius_search,
                                                double radius,
                                                size_t* out_num_ok_queries,
                                                scc_PointIndex out_query_indices[],
                                                scc_PointIndex out_nn_indices[])
{
	return iscc_dist_functions.nearest_neighbor_search(nn_search_object,
	                                                   len_query_indices,
	                                                   query_indices,
	                                                   k,
	                                                   radius_search,
	                                                   radius,
	                                                   out_num_ok_queries,
	                                                   out_query_indices,
	                                                   out_nn_indices);
}


static inline bool iscc_close_nn_search_object(iscc_NNSearchObject** nn_search_object)
{
	return iscc_dist_functions.close_nn_search_object(nn_search_object);
}

#endif // ifndef SCC_DIST_SEARCH_HG
