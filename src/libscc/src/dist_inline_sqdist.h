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

#ifndef SCC_DIST_INLINE_SQDIST_HG
#define SCC_DIST_INLINE_SQDIST_HG

#include <assert.h>
#include <stddef.h>
#include "../include/scc_data_obj.h"
#include "scc_data_obj_int.h"


// ==============================================================================
// Inline function implementations
// ==============================================================================

static inline double iscc_get_sq_dist(const scc_DataSetObject* const data_set_object,
                                      size_t index1,
                                      size_t index2)
{
	assert(index1 < data_set_object->num_data_points);
	assert(index2 < data_set_object->num_data_points);

	const double* data1 = &data_set_object->data_matrix[index1 * data_set_object->num_dimensions];
	const double* const data1_stop = data1 + data_set_object->num_dimensions;
	const double* data2 = &data_set_object->data_matrix[index2 * data_set_object->num_dimensions];

	double tmp_dist = 0.0;
	while (data1 != data1_stop) {
		const double value_diff = (*data1 - *data2);
		++data1;
		++data2;
		tmp_dist += value_diff * value_diff;
	}
	return tmp_dist;
}

#endif // ifndef SCC_DIST_INLINE_SQDIST_HG
