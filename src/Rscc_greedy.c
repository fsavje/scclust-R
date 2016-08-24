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

#include "Rscc_greedy.h"

#include <limits.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <R.h>
#include <Rinternals.h>
#include <scclust.h>
#include <scc_data_obj.h>
#include "Rscc_misc.h"


// ==============================================================================
// External function implementations
// ============================================================================== 

SEXP Rsccwrap_top_down_greedy_clustering(const SEXP R_distance_object,
                                         const SEXP R_size_constraint,
                                         const SEXP R_batch_assign,
                                         const SEXP R_existing_clustering,
                                         const SEXP R_existing_num_clusters,
                                         const SEXP R_deep_copy)
{
	if (!isMatrix(R_distance_object) || !isReal(R_distance_object)) iRsccwrap_error("Invalid distance object.");
	if (!isInteger(R_size_constraint)) iRsccwrap_error("`R_size_constraint` must be integer.");
	if (!isLogical(R_batch_assign)) iRsccwrap_error("`R_batch_assign` must be logical.");
	if (!isNull(R_existing_clustering) && !isInteger(R_existing_clustering)) iRsccwrap_error("`R_existing_clustering` must be integer.");
	if (!isInteger(R_existing_num_clusters)) iRsccwrap_error("`R_existing_num_clusters` must be integer.");
	if (!isLogical(R_deep_copy)) iRsccwrap_error("`R_deep_copy` must be logical.");

	const uintmax_t num_data_points = (uintmax_t) INTEGER(getAttrib(R_distance_object, R_DimSymbol))[1];
	const uintmax_t num_dimensions = (uintmax_t) INTEGER(getAttrib(R_distance_object, R_DimSymbol))[0];
	const uint32_t size_constraint = (uint32_t) asInteger(R_size_constraint);
	const bool batch_assign = (bool) asLogical(R_batch_assign);
	const uintmax_t existing_num_clusters = (uintmax_t) asInteger(R_existing_num_clusters);
	const bool deep_copy = (bool) asLogical(R_deep_copy);

	scc_ErrorCode ec;
	scc_DataSetObject* data_set_object;
	if ((ec = scc_get_data_set_object(num_data_points,
	                                  num_dimensions,
	                                  (size_t) xlength(R_distance_object),
	                                  REAL(R_distance_object),
	                                  false,
	                                  &data_set_object)) != SCC_ER_OK) {
		iRsccwrap_scc_error();
	}

	SEXP R_cluster_labels;
	scc_Clustering* clustering;
	if (isNull(R_existing_clustering)) {
		if (existing_num_clusters > 0) iRsccwrap_error("`R_existing_num_clusters` must be zero when creating new clustering.");
		R_cluster_labels = PROTECT(allocVector(INTSXP, (R_xlen_t) num_data_points));
		if ((ec = scc_init_empty_clustering(num_data_points,
		                                    INTEGER(R_cluster_labels),
		                                    &clustering)) != SCC_ER_OK) {
			scc_free_data_set_object(&data_set_object);
			UNPROTECT(1);
			iRsccwrap_scc_error();
		}
	} else {
		if (!isInteger(R_existing_clustering)) iRsccwrap_error("`R_existing_clustering` must be integer.");
		if (xlength(R_existing_clustering) != num_data_points) iRsccwrap_error("Existing clustering does not match distance object.");
		if (existing_num_clusters == 0) iRsccwrap_error("`R_existing_num_clusters` must be non-zero when using existing clustering.");

		if (deep_copy) {
			R_cluster_labels = PROTECT(duplicate(R_existing_clustering));
		} else {
			R_cluster_labels = PROTECT(R_existing_clustering);
		}

		if ((ec = scc_init_existing_clustering(num_data_points,
		                                       existing_num_clusters,
		                                       INTEGER(R_existing_clustering),
		                                       false,
		                                       &clustering)) != SCC_ER_OK) {
			scc_free_data_set_object(&data_set_object);
			UNPROTECT(1);
			iRsccwrap_scc_error();
		}
	}

	if ((ec = scc_top_down_greedy_clustering(clustering,
	                                         data_set_object,
	                                         size_constraint,
	                                         batch_assign)) != SCC_ER_OK) {
		scc_free_clustering(&clustering);
		scc_free_data_set_object(&data_set_object);
		UNPROTECT(1);
		iRsccwrap_scc_error();
	}

	scc_free_data_set_object(&data_set_object);

	uintmax_t num_clusters = 0;
	if ((ec = scc_get_clustering_info(clustering,
	                                  NULL,
	                                  &num_clusters)) != SCC_ER_OK) {
		scc_free_clustering(&clustering);
		UNPROTECT(1);
		iRsccwrap_scc_error();
	}

	scc_free_clustering(&clustering);

	if (num_clusters > INT_MAX) iRsccwrap_error("Too many clusters.");
	if (num_clusters < existing_num_clusters) iRsccwrap_error("Clustering failed: number of clusters decreased.");
	const int num_clusters_int = (int) num_clusters;

	const SEXP R_clustering_obj = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(R_clustering_obj, 0, R_cluster_labels);
	SET_VECTOR_ELT(R_clustering_obj, 1, ScalarInteger(num_clusters_int));

	const SEXP R_obj_elem_names = PROTECT(allocVector(STRSXP, 2));
	SET_STRING_ELT(R_obj_elem_names, 0, mkChar("cluster_labels"));
	SET_STRING_ELT(R_obj_elem_names, 1, mkChar("cluster_count"));
	setAttrib(R_clustering_obj, R_NamesSymbol, R_obj_elem_names);

	UNPROTECT(3);
	return R_clustering_obj;
}
