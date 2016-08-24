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
 * Operations on digraphs.
 */

#ifndef SCC_DIGRAPH_OPERATIONS_HG
#define SCC_DIGRAPH_OPERATIONS_HG

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "digraph_core.h"


// ==============================================================================
// Function prototypes
// ==============================================================================

/** Delete all self-loops.
 *
 *  This function deletes all arcs where the tail and the head is the same vertex.
 *
 *  \param[in,out] dg digraph to delete self-loops from.
 *
 *  \note Arc memory space that is freed due to the deletion is deallocated.
 *
 *  \note The deletion is stable so that the internal ordering of remaining arcs in \p dg->head is unchanged.
 */
scc_ErrorCode iscc_delete_loops(iscc_Digraph* dg);

/** Calculates the union of arbitrary number of digraphs.
 *
 *  This function produces the union of the inputted digraphs. Optionally, the function can also delete
 *  arc as indicated by tails.
 *
 *  The union of the following first two digraphs is the third digraph:
 *  \dot
 *  digraph example {
 *      A -> C;
 *      B -> C;
 *      C -> D;
 *
 *      A2 [ label="A" ];
 *      B2 [ label="B" ];
 *      C2 [ label="C" ];
 *      D2 [ label="D" ];
 *      A2 -> C2;
 *      B2 -> D2;
 *      A2 -> D2;
 *
 *      A3 [ label="A" ];
 *      B3 [ label="B" ];
 *      C3 [ label="C" ];
 *      D3 [ label="D" ];
 *
 *      A3 -> C3;
 *      B3 -> C3;
 *      C3 -> D3;
 *      B3 -> D3;
 *      A3 -> D3;
 *  }
 *  \enddot
 *
 *  \param num_dgs number of digraph to calculate union for. Must be non-zero.
 *  \param[in] dgs the digraphs. Must be of length \p num_dgs.
 *  \param[in] tails_to_keep indicators of tails for which the arcs should be *kept*.
 *                           If `NULL` no arcs (except self-loops) are deleted.
 *                           If not `NULL`, it must be of the same length as the number of vertices in the digraphs.
 *  \param[out] out_dg the union of \p dgs.
 *
 *  \note All digraphs in \p dgs must contain equally many vertices.
 *  \note All self-loops in the digraphs will be ignored.
 */
scc_ErrorCode iscc_digraph_union_and_delete(uint_fast16_t num_in_dgs,
                                            const iscc_Digraph in_dgs[static num_in_dgs],
                                            const bool tails_to_keep[],
                                            bool keep_self_loops,
                                            iscc_Digraph* out_dg);

scc_ErrorCode iscc_digraph_difference(iscc_Digraph* minuend_dg,
                                      const iscc_Digraph* subtrahend_dg,
                                      uint32_t max_out_degree);

/** Derives the digraph transpose a digraph.
 *
 *  This function produces the transpose of the inputted digraph. That is, a digraph where
 *  all the arcs are reversed.
 *
 *  The transpose of the following first digraph is the second digraph:
 *  \dot
 *  digraph example {
 *      A -> C;
 *      B -> C;
 *      C -> D;
 *
 *      A2 [ label="A" ];
 *      B2 [ label="B" ];
 *      C2 [ label="C" ];
 *      D2 [ label="D" ];
 *      C2 -> A2;
 *      C2 -> B2;
 *      D2 -> C2;
 *  }
 *  \enddot
 *
 *  \param[in] dg digraph to transpose.
 *  \param[out] out_dg the transpose of \p dg.
 */
scc_ErrorCode iscc_digraph_transpose(const iscc_Digraph* in_dg,
                                     iscc_Digraph* out_dg);

/** Calculates the product of the adjacency matrices of two digraphs.
 *
 *  Digraphs are stored as sparse adjacency matrices. By multiplying the underlying
 *  matrices one can derive paths and powers of the graphs.
 *
 *  Let \f$A\f$ and \f$B\f$ be the adjacency matrices of two arbitrary digraph. #scc_adjacency_product
 *  returns \f$A B\f$.
 *
 *  The main purpose of this function is to derive paths and powers. Let \f$A\f$ be an adjacency matrix of some digraph.
 *  \f$A A\f$ then gives the adjacency matrix of the digraph that contains all path of length two in the original digraph.
 *  Let \f$I\f$ be the identity matrix, then \f$(A + I) A\f$ gives the second power of the digraph. That is,
 *  all possible paths of length two or less. Moreover, if \f$A_2 = A A\f$ then \f$A A_2\f$ gives
 *  all paths of length three and \f$A_2 A_2\f$ of length four. If \f$A_p = (A + I) A\f$, then \f$(A + I) A_p\f$
 *  gives the third power and \f$(A_p + I) A_p\f$ gives the fourth power.
 *
 *  \code
 *	// Dummy digraph generator
 *  scc_Digraph my_dg = some_digraph(); 
 *
 *  // All paths of length 2 in `my_dg`
 *  scc_Digraph my_dg_path2 = scc_adjacency_product(&my_dg, &my_dg, false);
 *  
 *  // All paths of length 3 in `my_dg`
 *  scc_Digraph my_dg_path3 = scc_adjacency_product(&my_dg, &my_dg_path2, false);
 *  
 *  // Second power of `my_dg`
 *  scc_Digraph my_dg_power2 = scc_adjacency_product(&my_dg, &my_dg, true);
 *  
 *  // Fourth power of `my_dg`
 *  scc_Digraph my_dg_power4 = scc_adjacency_product(&my_dg_power2, &my_dg_power2, true);
 *  
 *  // Free all digraphs
 *  scc_free_digraph(&my_dg); scc_free_digraph(&my_dg_path2); [...]
 *  \endcode
 *  
 *  \param[in] dg_a the first digraph of the product.
 *  \param[in] dg_b the second digraph of the product.
 *  \param     force_loops when \c true, forces self-loops in \p dg_a (i.e., all vertices have an arc to themselves).
 *  \param[out] out_dg the digraph described by the product of the adjacency matrices of \p dg_a and \p dg_b.
 *
 *  \note \p dg_a and \p dg_b must contain equally many vertices.
 *
 *  \note The output digraph will never have self-loops (independently of \p force_loops).
 */
scc_ErrorCode iscc_adjacency_product(const iscc_Digraph* in_dg_a,
                                     const iscc_Digraph* in_dg_b,
                                     bool force_loops,
                                     iscc_Digraph* out_dg);


#endif // ifndef SCC_DIGRAPH_OPERATIONS_HG
