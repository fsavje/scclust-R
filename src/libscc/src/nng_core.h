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

#ifndef SCC_NNG_CORE_HG
#define SCC_NNG_CORE_HG

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "../include/scclust.h"
#include "digraph_core.h"
#include "nng_findseeds.h"
#include "scclust_internal.h"


// =============================================================================
// Function prototypes
// =============================================================================

scc_ErrorCode iscc_get_nng_with_size_constraint(void* data_set_object,
                                                size_t num_data_points,
                                                uint32_t size_constraint,
                                                const bool main_data_points[],
                                                bool radius_constraint,
                                                double radius,
                                                iscc_Digraph* out_nng);

scc_ErrorCode iscc_get_nng_with_type_constraint(void* data_set_object,
                                                size_t num_data_points,
                                                uint32_t size_constraint,
                                                uint_fast16_t num_types,
                                                const uint32_t type_size_constraints[static num_types],
                                                const scc_TypeLabel type_labels[static num_data_points],
                                                const bool main_data_points[],
                                                bool radius_constraint,
                                                double radius,
                                                iscc_Digraph* out_nng);

scc_ErrorCode iscc_estimate_avg_seed_dist(void* data_set_object,
                                          const iscc_SeedResult* seed_result,
                                          const iscc_Digraph* nng,
                                          uint32_t size_constraint,
                                          double* out_avg_seed_dist);

scc_ErrorCode iscc_make_nng_clusters_from_seeds(scc_Clustering* clustering,
                                                void* data_set_object,
                                                const iscc_SeedResult* seed_result,
                                                iscc_Digraph* nng,
                                                bool nng_is_ordered,
                                                scc_UnassignedMethod main_unassigned_method,
                                                bool main_radius_constraint,
                                                double main_radius,
                                                const bool main_data_points[],
                                                scc_UnassignedMethod secondary_unassigned_method,
                                                bool secondary_radius_constraint,
                                                double secondary_radius);


#endif // ifndef SCC_NNG_CORE_HG
