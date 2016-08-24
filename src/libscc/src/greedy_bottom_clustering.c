/* ==============================================================================
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
 * ============================================================================== */

#include "../include/scclust.h"
#include "error.h"

// ==============================================================================
// External function implementations
// ==============================================================================

scc_ErrorCode scc_bottom_up_greedy_clustering(scc_Clustering* const clustering,
                                              void* const data_set_object,
                                              const uint32_t size_constraint)
{
	(void) clustering;
	(void) data_set_object;
	(void) size_constraint;
	return iscc_make_error(SCC_ER_NOT_IMPLEMENTED);
}
