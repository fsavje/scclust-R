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

#ifndef SCC_BATCH_CLUSTERING_HG
#define SCC_BATCH_CLUSTERING_HG

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "../include/scclust.h"


// =============================================================================
// Function prototypes
// =============================================================================

scc_ErrorCode scc_nng_clustering_batches(scc_Clustering* clustering,
                                         void* data_set,
                                         uint32_t size_constraint,
                                         scc_UnassignedMethod unassigned_method,
                                         bool radius_constraint,
                                         double radius,
                                         size_t len_primary_data_points,
                                         const bool primary_data_points[],
                                         uint32_t batch_size);


#endif // ifndef SCC_BATCH_CLUSTERING_HG
