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

/** @file
 *
 *  Functions for debugging digraphs.
 *
 *  These function should in general not be used in production, but they can
 *  be useful when debugging.
 */

#ifndef SCC_DIGRAPH_DEBUG_HG
#define SCC_DIGRAPH_DEBUG_HG

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "digraph_core.h"
#include "scclust_int.h"


// ==============================================================================
// Function prototypes
// ==============================================================================

/** Checks whether provided digraph is a balanced digraph.
 *
 *  This function returns \c true if \p dg is valid and all vertices have equally many
 *  outwards pointing arcs.
 *
 *  \param[in] dg digraph to check.
 *  \param arcs_per_vertex number of outwards pointing arcs each vertex is expected to have.
 *
 *  \return \c true if \p dg is balanced, otherwise \c false.
 */
bool iscc_is_balanced_digraph(const iscc_Digraph* dg,
                              iscc_Arci arcs_per_vertex);

/** Checks whether two digraphs are logically identical.
 *
 *  Two digraphs are considered logically identical if they contain
 *  equally many vertices and they contain the same set of arcs as judged
 *  by vertex IDs. Duplicate arcs are ignored. 
 *
 *  \param[in] dg_a first digraph to compare.
 *  \param[in] dg_b second digraph to compare.
 *
 *  \return \c true if \p dg_a and \p dg_b are identical.
 */
bool iscc_digraphs_equal(const iscc_Digraph* dg_a,
                         const iscc_Digraph* dg_b);

/** Constructs digraph from separate struct parts.
 *
 *  This function constructs a new scc_Digraph instance and initializes it
 *  with the provided values. It takes a deep copy of \p tail_ptr and \p head.
 *
 *  \param vertices number of vertices in the digraph.
 *  \param max_arcs length of the memory space for arcs.
 *  \param[in] tail_ptr array of arc indices of length `vertices + 1` (see scc_Digraph::tail_ptr).
 *  \param[in] head array of arc heads of length `max_arcs` (see scc_Digraph::head).
 *  \param[out] out_dg the constructed scc_Digraph.
 */
scc_ErrorCode iscc_digraph_from_pieces(size_t vertices,
                                       uintmax_t max_arcs,
                                       const iscc_Arci tail_ptr[static vertices + 1],
                                       const iscc_Dpid head[static max_arcs],
                                       iscc_Digraph* out_dg);

/** Constructs digraph from human readable strings.
 *
 *  This function builds a digraph from a string describing an adjacency matrix, where rows
 *  indicate tails and columns heads. `#` denotes an arc, `.` denotes absence of an arc and
 *  `/` the end of a row. All other characters are ignored.
 *  
 *  Thus, the string ".##./..#./...#/#..#/" describes the following digraph:
 *  \dot
 *  digraph example {
 *      0 -> 1;
 *      0 -> 2;
 *      1 -> 2;
 *      2 -> 3;
 *      3 -> 0;
 *      3 -> 3;
 *  }
 *  \enddot
 *
 *  \param[in] dg_str a null terminated string describing an adjacency matrix.
 *  \param[out] out_dg the digraph described by \p dg_str.
 */
scc_ErrorCode iscc_digraph_from_string(const char dg_str[],
                                       iscc_Digraph* out_dg);

/** Deep copy of a digraph.
 *
 *  This function produces a deep copy of the inputted digraph.
 *
 *  \param[in] dg digraph to copy.
 *  \param[out] out_dg a copy of \p dg that does not share memory space with it.
 *
 *  \note This function allocates memory space to fit the arcs actually in \p dg. If \p dg
 *        contains excess space, scc_Digraph::max_arcs will differ between the original and copy.
 */
scc_ErrorCode iscc_copy_digraph(const iscc_Digraph* in_dg,
                                iscc_Digraph* out_dg);

/** Print a digraph in human readable format.
 *
 *  This function prints the adjacency matrix of the provided digraph. Rows
 *  indicate tails and columns heads. `#` denotes an arc, `.` denotes absence of an arc.
 *  Digraphs without vertices are printed `[]`.
 *
 *  Thus, the following digraph:
 *  \dot
 *  digraph example {
 *      0 -> 1;
 *      0 -> 2;
 *      1 -> 2;
 *      2 -> 3;
 *      3 -> 0;
 *      3 -> 3;
 *  }
 *  \enddot
 *  would be printed as:
 *  \code
 *	.##.
 *	..#.
 *	...#
 *	#..#
 *  \endcode
 *
 *  \param[in] dg digraph to print.
 */
void iscc_print_digraph(const iscc_Digraph* dg);


#endif // ifndef SCC_DIGRAPH_DEBUG_HG
