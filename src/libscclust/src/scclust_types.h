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
#include "../include/scclust_spi.h"


/** Type used for arc indices. Must be unsigned.
 *
 *  \note
 *  Number of arcs in any digraph must be less or equal to
 *  the maximum number that can be stored in #iscc_ArcIndex.
 */
typedef uint32_t iscc_ArcIndex;

#define ISCC_M_ARCINDEX_TYPE_uint32_t

static const scc_Clabel SCC_CLABEL_MAX = INT_MAX;
static const scc_PointIndex ISCC_POINTINDEX_MAX_PI = INT_MAX;
static const uintmax_t ISCC_POINTINDEX_MAX = INT_MAX;
static const uintmax_t ISCC_ARCINDEX_MAX = UINT32_MAX;
static const uintmax_t ISCC_TYPELABEL_MAX = 65535;

#define ISCC_M_CLABEL_MAX INT_MAX
#define ISCC_M_POINTINDEX_MAX INT_MAX
#define ISCC_M_ARCINDEX_MAX UINT32_MAX
#define ISCC_M_TYPELABEL_MAX 65535


#endif // ifndef SCC_SCCLUST_INTERNAL_HG
