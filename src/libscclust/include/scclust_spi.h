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

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "scclust.h"

#ifdef __cplusplus
extern "C" {
#endif


// =============================================================================
// Internal types
// =============================================================================

typedef struct iscc_MaxDistObject iscc_MaxDistObject;

typedef struct iscc_NNSearchObject iscc_NNSearchObject;


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

typedef bool (*scc_nearest_neighbor_search) (iscc_NNSearchObject*,
                                             size_t,
                                             const scc_PointIndex*,
                                             uint32_t,
                                             bool,
                                             double,
                                             size_t*,
                                             scc_PointIndex*,
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
                            scc_nearest_neighbor_search,
                            scc_close_nn_search_object);


#ifdef __cplusplus
}
#endif

#endif // ifndef SCC_SCCLUST_HG
