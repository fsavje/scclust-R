/* ==============================================================================
 * Rscclust -- R wrapper for the scclust library
 * https://github.com/fsavje/Rscclust
 *
 * Copyright (C) 2016  Fredrik Savje -- http://fredriksavje.com
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/
 * ============================================================================== */

#ifndef RSCCWRAP_GREEDY_HG
#define RSCCWRAP_GREEDY_HG

#include <R.h>
#include <Rinternals.h>

SEXP Rsccwrap_top_down_greedy_clustering(SEXP R_distance_object,
                                         SEXP R_size_constraint,
                                         SEXP R_batch_assign,
                                         SEXP R_existing_clustering,
                                         SEXP R_existing_num_clusters,
                                         SEXP R_deep_copy);

#endif // ifndef RSCCWRAP_GREEDY_HG
