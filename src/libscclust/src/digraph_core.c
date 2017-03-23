/* =============================================================================
 * scclust -- A C library for size constrained clustering
 * https://github.com/fsavje/scclust
 *
 * Copyright (C) 2015-2017  Fredrik Savje -- http://fredriksavje.com
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

#include "digraph_core.h"

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include "../include/scclust.h"
#include "error.h"
#include "scclust_types.h"


// =============================================================================
// External function implementations
// =============================================================================

void iscc_free_digraph(iscc_Digraph* const dg)
{
	if (dg != NULL) {
		free(dg->head);
		free(dg->tail_ptr);
		*dg = ISCC_NULL_DIGRAPH;
	}
}


bool iscc_digraph_is_initialized(const iscc_Digraph* const dg)
{
	if ((dg == NULL) || (dg->tail_ptr == NULL)) return false;
	if ((dg->vertices > ISCC_POINTINDEX_MAX) || (dg->max_arcs > ISCC_ARCINDEX_MAX)) return false;
	if ((dg->max_arcs == 0) && (dg->head != NULL)) return false;
	if ((dg->max_arcs > 0) && (dg->head == NULL)) return false;
	return true;
}


bool iscc_digraph_is_valid(const iscc_Digraph* const dg)
{
	if (!iscc_digraph_is_initialized(dg)) return false;
	if (dg->tail_ptr[0] != 0) return false;
	if (dg->tail_ptr[dg->vertices] > dg->max_arcs) return false;
	for (size_t i = 0; i < dg->vertices; ++i) {
		if (dg->tail_ptr[i] > dg->tail_ptr[i + 1]) return false;
	}
	if (dg->tail_ptr[dg->vertices] > 0) {
		assert(dg->vertices <= ISCC_POINTINDEX_MAX);
		scc_PointIndex vertices = (scc_PointIndex) dg->vertices; // If `scc_PointIndex` is signed.
		const scc_PointIndex* const arc_stop = dg->head + dg->tail_ptr[dg->vertices];
		for (const scc_PointIndex* arc = dg->head; arc != arc_stop; ++arc) {
			if (*arc >= vertices) return false;
		}
	}
	return true;
}


bool iscc_digraph_is_empty(const iscc_Digraph* const dg)
{
	assert(iscc_digraph_is_initialized(dg));
	return (dg->tail_ptr[dg->vertices] == 0);
}


scc_ErrorCode iscc_init_digraph(const size_t vertices,
                                const uintmax_t max_arcs,
                                iscc_Digraph* const out_dg)
{
	assert(vertices > 0);
	assert(vertices <= ISCC_POINTINDEX_MAX);
	assert(vertices < SIZE_MAX);
	assert(out_dg != NULL);
	if ((max_arcs > ISCC_ARCINDEX_MAX) || (max_arcs > SIZE_MAX)) {
		return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many arcs in graph (adjust the `iscc_ArcIndex` type).");
	}

	*out_dg = (iscc_Digraph) {
		.vertices = vertices,
		.max_arcs = (size_t) max_arcs,
		.head = NULL,
		.tail_ptr = malloc(sizeof(iscc_ArcIndex[vertices + 1])),
	};
	if (out_dg->tail_ptr == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);

	if (max_arcs > 0) {
		out_dg->head = malloc(sizeof(scc_PointIndex[max_arcs]));
		if (out_dg->head == NULL) {
			iscc_free_digraph(out_dg);
			return iscc_make_error(SCC_ER_NO_MEMORY);
		}
	}

	assert(iscc_digraph_is_initialized(out_dg));

	return iscc_no_error();
}


scc_ErrorCode iscc_empty_digraph(const size_t vertices,
                                 const uintmax_t max_arcs,
                                 iscc_Digraph* const out_dg)
{
	assert(vertices > 0);
	assert(vertices <= ISCC_POINTINDEX_MAX);
	assert(vertices < SIZE_MAX);
	assert(out_dg != NULL);
	if ((max_arcs > ISCC_ARCINDEX_MAX) || (max_arcs > SIZE_MAX)) {
		return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many arcs in graph (adjust the `iscc_ArcIndex` type).");
	}

	*out_dg = (iscc_Digraph) {
		.vertices = vertices,
		.max_arcs = (size_t) max_arcs,
		.head = NULL,
		.tail_ptr = calloc(vertices + 1, sizeof(iscc_ArcIndex)),
	};
	if (out_dg->tail_ptr == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);

	if (max_arcs > 0) {
		out_dg->head = malloc(sizeof(scc_PointIndex[max_arcs]));
		if (out_dg->head == NULL) {
			iscc_free_digraph(out_dg);
			return iscc_make_error(SCC_ER_NO_MEMORY);
		}
	}

	assert(iscc_digraph_is_valid(out_dg));

	return iscc_no_error();
}


scc_ErrorCode iscc_change_arc_storage(iscc_Digraph* const dg,
                                      const uintmax_t new_max_arcs)
{
	assert(iscc_digraph_is_initialized(dg));
	assert(dg->tail_ptr[dg->vertices] <= new_max_arcs);
	if ((new_max_arcs > ISCC_ARCINDEX_MAX) || (new_max_arcs > SIZE_MAX)) {
		return iscc_make_error_msg(SCC_ER_TOO_LARGE_PROBLEM, "Too many arcs in graph (adjust the `iscc_ArcIndex` type).");
	}
	if (dg->max_arcs == new_max_arcs) return iscc_no_error();

	if (new_max_arcs == 0) {
		free(dg->head);
		dg->head = NULL;
		dg->max_arcs = 0;
	} else {
		scc_PointIndex* const tmp_ptr = realloc(dg->head, sizeof(scc_PointIndex[new_max_arcs]));
		if (tmp_ptr == NULL) return iscc_make_error(SCC_ER_NO_MEMORY);
		dg->head = tmp_ptr;
		dg->max_arcs = (size_t) new_max_arcs;
	}

	return iscc_no_error();
}
