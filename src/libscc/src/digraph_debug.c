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

#include "digraph_debug.h"

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "digraph_core.h"
#include "error.h"
#include "scclust_internal.h"


// =============================================================================
// External function implementations
// =============================================================================

bool iscc_is_balanced_digraph(const iscc_Digraph* const dg,
                              const iscc_Arci arcs_per_vertex)
{
	if (!iscc_digraph_is_valid(dg)) return false;

	for (size_t i = 0; i <= dg->vertices; ++i) {
		if (dg->tail_ptr[i] != i * arcs_per_vertex) return false;
	}

	return true;
}


bool iscc_digraphs_equal(const iscc_Digraph* const dg_a,
                         const iscc_Digraph* const dg_b)
{
	assert(iscc_digraph_is_valid(dg_a));
	assert(iscc_digraph_is_valid(dg_b));
	if (dg_a->vertices != dg_b->vertices) return false;
	if ((dg_a->tail_ptr[dg_a->vertices] == 0) && (dg_b->tail_ptr[dg_b->vertices] == 0)) return true;

	int_fast8_t* const single_row = calloc(dg_a->vertices, sizeof(int_fast8_t));

	for (size_t v = 0; v < dg_a->vertices; ++v) {
		const iscc_Dpid* const arc_a_stop = dg_a->head + dg_a->tail_ptr[v + 1];
		for (const iscc_Dpid* arc_a = dg_a->head + dg_a->tail_ptr[v];
		        arc_a != arc_a_stop; ++arc_a) {
			single_row[*arc_a] = 1;
		}

		const iscc_Dpid* const arc_b_stop = dg_b->head + dg_b->tail_ptr[v + 1];
		for (const iscc_Dpid* arc_b = dg_b->head + dg_b->tail_ptr[v];
		        arc_b != arc_b_stop; ++arc_b) {
			if (single_row[*arc_b] == 0) {
				free(single_row);
				return false;
			}
			single_row[*arc_b] = 2;
		}

		for (size_t i = 0; i < dg_a->vertices; ++i) {
			if (single_row[i] == 1) {
				free(single_row);
				return false;
			}
			single_row[i] = 0;
		}
	}

	free(single_row);

	return true;
}


scc_ErrorCode iscc_digraph_from_pieces(const size_t vertices,
                                       const uintmax_t max_arcs,
                                       const iscc_Arci tail_ptr[const static vertices + 1],
                                       const iscc_Dpid head[const static max_arcs],
                                       iscc_Digraph* const out_dg)
{
	assert(vertices > 0);
	assert(vertices <= ISCC_DPID_MAX);
	assert(max_arcs > 0);
	assert(max_arcs <= ISCC_ARCI_MAX);
	assert(max_arcs <= SIZE_MAX);
	assert(tail_ptr != NULL);
	assert(head != NULL);
	assert(out_dg != NULL);

	scc_ErrorCode ec;
	if ((ec = iscc_init_digraph(vertices, max_arcs, out_dg)) != SCC_ER_OK) return ec;

	memcpy(out_dg->tail_ptr, tail_ptr, (vertices + 1) * sizeof(iscc_Arci));
	memcpy(out_dg->head, head, max_arcs * sizeof(iscc_Dpid));

	return iscc_no_error();
}


scc_ErrorCode iscc_digraph_from_string(const char dg_str[const],
                                       iscc_Digraph* const out_dg)
{
	assert(dg_str != NULL);
	assert(out_dg != NULL);

	size_t vertices = 0;
	uintmax_t all_arcs = 0;
	uintmax_t max_arcs = 0;

	for (size_t c = 0; dg_str[c] != '\0'; ++c) {
		if (dg_str[c] == '#' || dg_str[c] == '.') ++all_arcs;
		if (dg_str[c] == '#') ++max_arcs;
		if (dg_str[c] == '/' && vertices == 0) vertices = all_arcs;
		if (dg_str[c] == '/' && (all_arcs % vertices) != 0) return iscc_make_error(SCC_ER_INVALID_INPUT);
	}

	scc_ErrorCode ec;
	if ((ec = iscc_init_digraph(vertices, max_arcs, out_dg)) != SCC_ER_OK) return ec;

	iscc_Arci curr_array_pos = 0;
	size_t curr_row = 0;
	iscc_Dpid curr_col = 0;
	out_dg->tail_ptr[0] = 0;

	for (size_t c = 0; dg_str[c] != '\0'; ++c) {
		if (dg_str[c] == '#') {
			out_dg->head[curr_array_pos] = curr_col;
			++curr_array_pos;
		}
		if (dg_str[c] == '#' || dg_str[c] == '.') ++curr_col;
		if (dg_str[c] == '/') {
			++curr_row;
			curr_col = 0;
			out_dg->tail_ptr[curr_row] = curr_array_pos;
		}
	}
	out_dg->tail_ptr[vertices] = curr_array_pos;

	assert(iscc_digraph_is_valid(out_dg));

	return iscc_no_error();
}


scc_ErrorCode iscc_copy_digraph(const iscc_Digraph* const in_dg,
                                iscc_Digraph* const out_dg)
{
	scc_ErrorCode ec;
	assert(iscc_digraph_is_initialized(in_dg));
	assert(out_dg != NULL);
	if (in_dg->vertices == 0) return iscc_empty_digraph(0, 0, out_dg);

	const size_t num_vertices = in_dg->vertices;
	const uintmax_t num_arcs = in_dg->tail_ptr[in_dg->vertices];

	if ((ec = iscc_init_digraph(num_vertices, num_arcs, out_dg)) != SCC_ER_OK) return ec;

	memcpy(out_dg->tail_ptr, in_dg->tail_ptr, (num_vertices + 1) * sizeof(iscc_Arci));
	if (num_arcs > 0) {
		memcpy(out_dg->head, in_dg->head, num_arcs * sizeof(iscc_Dpid));
	}

	return iscc_no_error();
}


void iscc_print_digraph(const iscc_Digraph* const dg)
{
	assert(iscc_digraph_is_initialized(dg));

	if (dg->vertices == 0) {
		printf("[]\n\n");
		return;
	}

	bool* const single_row = calloc(dg->vertices, sizeof(bool));
	if (single_row == NULL) {
		printf("Out of memory.\n\n");
		return;
	}

	for (size_t v = 0; v < dg->vertices; ++v) {
		const iscc_Dpid* const a_stop = dg->head + dg->tail_ptr[v + 1];
		for (const iscc_Dpid* a = dg->head + dg->tail_ptr[v];
		        a != a_stop; ++a) {
			single_row[*a] = true;
		}

		for (size_t i = 0; i < dg->vertices; ++i) {
			if (single_row[i]) {
				putchar('#');
			} else {
				putchar('.');
			}
			single_row[i] = false;
		}

		putchar('\n');
	}
	putchar('\n');

	free(single_row);
}
