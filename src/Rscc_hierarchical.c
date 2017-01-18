/* =============================================================================
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
 * ========================================================================== */

#include "Rscc_hierarchical.h"

#include <limits.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <R.h>
#include <Rinternals.h>
#include <scclust.h>
#include <ann_wrapper.h>
#include "Rscc_error.h"


// =============================================================================
// External function implementations
// =============================================================================

SEXP Rscc_hierarchical_clustering(const SEXP R_distance_object,
                                  const SEXP R_size_constraint,
                                  const SEXP R_batch_assign,
                                  const SEXP R_existing_clustering,
                                  const SEXP R_deep_copy)
{
	if (!scc_set_ann_dist_search()) {
		iRscc_error("Cannot change NN search functions to ANN.");
	}
	if (!isMatrix(R_distance_object) || !isReal(R_distance_object)) {
		iRscc_error("`R_distance_object` is not a valid distance object.");
	}
	if (!isInteger(R_size_constraint)) {
		iRscc_error("`R_size_constraint` must be integer.");
	}
	if (!isLogical(R_batch_assign)) {
		iRscc_error("`R_batch_assign` must be logical.");
	}
	if (!isNull(R_existing_clustering) && !isInteger(R_existing_clustering)) {
		iRscc_error("`R_existing_clustering` is not a valid clustering object.");
	}
	if (!isLogical(R_deep_copy)) {
		iRscc_error("`R_deep_copy` must be logical.");
	}

	const uintmax_t num_data_points = (uintmax_t) INTEGER(getAttrib(R_distance_object, R_DimSymbol))[1];
	const uintmax_t num_dimensions = (uintmax_t) INTEGER(getAttrib(R_distance_object, R_DimSymbol))[0];
	const uint32_t size_constraint = (uint32_t) asInteger(R_size_constraint);
	const bool batch_assign = (bool) asLogical(R_batch_assign);
	const bool deep_copy = (bool) asLogical(R_deep_copy);

	scc_ErrorCode ec;
	scc_DataSet* data_set;
	if ((ec = scc_init_data_set(num_data_points,
	                            num_dimensions,
	                            (size_t) xlength(R_distance_object),
	                            REAL(R_distance_object),
	                            &data_set)) != SCC_ER_OK) {
		iRscc_scc_error();
	}

	SEXP R_cluster_labels;
	scc_Clustering* clustering;
	if (isNull(R_existing_clustering)) {
		R_cluster_labels = PROTECT(allocVector(INTSXP, (R_xlen_t) num_data_points));
		if ((ec = scc_init_empty_clustering(num_data_points,
		                                    INTEGER(R_cluster_labels),
		                                    &clustering)) != SCC_ER_OK) {
			scc_free_data_set(&data_set);
			UNPROTECT(1);
			iRscc_scc_error();
		}
	} else {
		if (!isInteger(getAttrib(R_existing_clustering, install("cluster_count")))) {
			iRscc_error("`R_existing_clustering` is not a valid clustering object.");
		}
		if (((uintmax_t) xlength(R_existing_clustering)) != num_data_points) {
			iRscc_error("`R_existing_clustering` does not match `R_distance_object`.");
		}
		const uintmax_t existing_num_clusters = (uintmax_t) asInteger(getAttrib(R_existing_clustering, install("cluster_count")));
		if (existing_num_clusters == 0) {
			iRscc_error("`R_existing_clustering` is empty.");
		}

		if (deep_copy) {
			R_cluster_labels = PROTECT(duplicate(R_existing_clustering));
		} else {
			R_cluster_labels = PROTECT(R_existing_clustering);
		}

		setAttrib(R_cluster_labels, install("class"), R_NilValue);
		setAttrib(R_cluster_labels, install("cluster_count"), R_NilValue);
		setAttrib(R_cluster_labels, install("ids"), R_NilValue);

		if ((ec = scc_init_existing_clustering(num_data_points,
		                                       existing_num_clusters,
		                                       INTEGER(R_cluster_labels),
		                                       false,
		                                       &clustering)) != SCC_ER_OK) {
			scc_free_data_set(&data_set);
			UNPROTECT(1);
			iRscc_scc_error();
		}
	}

	if ((ec = scc_hierarchical_clustering(data_set,
	                                      clustering,
	                                      size_constraint,
	                                      batch_assign)) != SCC_ER_OK) {
		scc_free_clustering(&clustering);
		scc_free_data_set(&data_set);
		UNPROTECT(1);
		iRscc_scc_error();
	}

	scc_free_data_set(&data_set);

	uintmax_t num_clusters = 0;
	if ((ec = scc_get_clustering_info(clustering,
	                                  NULL,
	                                  &num_clusters)) != SCC_ER_OK) {
		scc_free_clustering(&clustering);
		UNPROTECT(1);
		iRscc_scc_error();
	}

	scc_free_clustering(&clustering);

	if (num_clusters > INT_MAX) iRscc_error("Too many clusters.");
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
