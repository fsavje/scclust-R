/* =============================================================================
 * scclust -- A C library for size-constrained clustering
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

/** @file
 *
 * Digraph structs and miscellaneous functions.
 *
 * Sparse digraphs are the backbone of NNG clustering in the scclust library.
 * This header defines the digraph struct and functions to generate them.
 */

#ifndef SCC_DIGRAPH_CORE_HG
#define SCC_DIGRAPH_CORE_HG

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "../include/scclust.h"
#include "scclust_types.h"


// =============================================================================
// Structs and variables
// =============================================================================

/** Main digraph struct stored as sparse matrix.
 *
 *  Stores the digraph in Yale sparse matrix format. For any vertex `i` in the digraph,
 *  `#tail_ptr[i]` indicates the arc index in #head of the first arc for which `i` is the tail, and `#tail_ptr[i+1]-1`
 *  indicates the last arc. If `#tail_ptr[i] == #tail_ptr[i+1]`, there exists no arc for which `i` is the tail.
 *  Thus if `i` is a tail for at least one arc, `#head[#tail_ptr[i]]` is the head of the first arc for which
 *  `i` is the tail, and `#head[#tail_ptr[i+1]-1]` is the last.
 *
 *  In other words, if there is an arc `i` -> `j`, there exists some `k` such that `#tail_ptr[i] <= k < #tail_ptr[i+1]` and `#head[k]==j`.
 */
typedef struct iscc_Digraph {
	/** Number of vertices in the digraph. May not be greater than `ISCC_POINTINDEX_MAX`.
	 *
	 *  \note Valid vertices in this digraph is any `i` such that `0 <= i < #vertices`.
	 */
	size_t vertices;

	/// Maximum number of arcs in digraph.
	size_t max_arcs;

	/** Array of vertex IDs indicating arc heads.
	 *
	 *  If `#max_arcs == 0`, #head must be `NULL`. If `#max_arcs > 0`,
	 *  #head should point a memory area of length #max_arcs.
	 *
	 *  \note All used elements of #head must be less than #vertices.
	 */
	scc_PointIndex* head;

	/** Array of arc indices indicating arcs for which a vertex is the tail.
	 *
	 *  #tail_ptr may never be `NULL` and must point a memory area of length `#vertices + 1`.
	 *
	 *  The first element of #tail_ptr must be zero (`#tail_ptr[0] == 0`). For all `i < #vertices`,
	 *  we must have `#tail_ptr[i] <= #tail_ptr[i+1] <= #max_arcs`.
	 */
	iscc_ArcIndex* tail_ptr;
} iscc_Digraph;


/** The null digraph.
 *
 *  The null digraph is an easily detectable invalid digraph.
 */
static const iscc_Digraph ISCC_NULL_DIGRAPH = { 0, 0, NULL, NULL };


// =============================================================================
// Function prototypes
// =============================================================================

/** Destructor for digraphs.
 *
 *  Frees the memory allocated by the inputted digraph and writes the null digraph to it.
 *
 *  \param[in,out] dg digraph to destroy. When #scc_free_digraph returns, \p dg is set to #SCC_NULL_DIGRAPH.
 */
void iscc_free_digraph(iscc_Digraph* dg);


/** Checks whether provided digraph is initialized.
 *
 *  This function returns \c true if \p dg is initialized. That is, scc_Digraph::tail_ptr
 *  and scc_Digraph::head are allocated. If scc_Digraph::max_arcs is zero, it checks so
 *  scc_Digraph::head is \c NULL.
 *
 *  \param[in] dg digraph to check.
 *
 *  \return \c true if \p dg is correctly initialized, otherwise \c false.
 *
 *  \note This function does not check whether \p dg is a valid digraph, that is whether
 *        the information is sound.
 */
bool iscc_digraph_is_initialized(const iscc_Digraph* dg);


/** Checks whether provided digraph is valid.
 *
 *  This function returns \c true if \p dg is a valid scc_Digraph instance. That is,
 *  \p dg describes a valid digraph.
 *
 *  \param[in] dg digraph to check.
 *
 *  \return \c true if \p dg is valid, otherwise \c false.
 */
bool iscc_digraph_is_valid(const iscc_Digraph* dg);


/** Checks whether provided digraph is empty.
 *
 *  This function returns \c true if \p dg does not contain any arcs.
 *
 *  \param[in] dg digraph to check.
 *
 *  \return \c true if \p dg is empty, otherwise \c false.
 *
 *  \note This function does not check whether \p dg is a valid digraph, that is whether
 *        the information is sound.
 */
bool iscc_digraph_is_empty(const iscc_Digraph* dg);


/** Generic constructor for digraphs.
 *
 *  Initializes and allocates memory for specified digraph. The memory spaces
 *  (i.e., scc_Digraph::head and scc_Digraph::tail_ptr) are uninitialized, thus
 *  the produced digraph is in general invalid.
 *
 *  \param vertices number of vertices that can be represented in the digraph.å
 *  \param max_arcs memory space to be allocated for arcs.
 *  \param[out] out_dg a scc_Digraph with allocated memory.
 */
scc_ErrorCode iscc_init_digraph(size_t vertices,
                                uintmax_t max_arcs,
                                iscc_Digraph* out_dg);


/** Construct an empty digraph.
 *
 *  This function returns a digraph where all elements of scc_Digraph::tail_ptr are set to `0`.
 *  The memory space pointed to by scc_Digraph::head is left uninitialized.
 *
 *  \param vertices number of vertices that can be represented in the digraph.
 *  \param max_arcs memory space to be allocated for arcs.
 *  \param[out] out_dg a scc_Digraph with allocated memory.
 */
scc_ErrorCode iscc_empty_digraph(size_t vertices,
                                 uintmax_t max_arcs,
                                 iscc_Digraph* out_dg);


/** Reallocate arc memory.
 *
 *  Increases or decreases the memory space for arcs in \p dg to fit exactly \p new_max_arcs arcs.
 *  Requires that the number of arcs in \p dg is less or equally to \p new_max_arcs.
 *  If `new_max_arcs == 0`, the memory space is deallocated and scc_Digraph::head is set to `NULL`.
 *
 *  \param[in,out] dg digraph to reallocate arc memory for.
 *  \param         new_max_arcs new size of memory.
 */
scc_ErrorCode iscc_change_arc_storage(iscc_Digraph* dg,
                                      uintmax_t new_max_arcs);


#endif // ifndef SCC_DIGRAPH_CORE_HG
