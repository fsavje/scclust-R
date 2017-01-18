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

#ifndef SCC_NNG_FINDSEEDS_HG
#define SCC_NNG_FINDSEEDS_HG

#include <stddef.h>
#include "../include/scclust.h"
#include "digraph_core.h"
#include "scclust_types.h"


// =============================================================================
// Structs, types and variables
// =============================================================================

typedef struct iscc_SeedResult iscc_SeedResult;
struct iscc_SeedResult {
	size_t capacity;
	size_t count;
	scc_PointIndex* seeds;
};


// =============================================================================
// Function prototypes
// =============================================================================

scc_ErrorCode iscc_find_seeds(const iscc_Digraph* nng,
                              scc_SeedMethod seed_method,
                              iscc_SeedResult* out_seeds);


#endif // ifndef SCC_NNG_FINDSEEDS_HG
