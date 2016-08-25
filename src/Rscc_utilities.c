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

#include "Rscc_utilities.h"

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

SEXP Rsccwrap_check_clustering(const SEXP R_cluster_labels,
                               const SEXP R_num_clusters,
                               const SEXP R_size_constraint)
{
	if (!isInteger(R_cluster_labels)) iRsccwrap_error("`R_cluster_labels` must be integer.");
	if (!isInteger(R_num_clusters)) iRsccwrap_error("`R_num_clusters` must be integer.");
	if (!isInteger(R_size_constraint)) iRsccwrap_error("`R_size_constraint` must be integer.");

	const uintmax_t num_data_points = (uintmax_t) xlength(R_cluster_labels);
	const uintmax_t num_clusters = (uintmax_t) asInteger(R_num_clusters);
	const uint32_t size_constraint = (uint32_t) asInteger(R_size_constraint);

	if (num_clusters == 0) iRsccwrap_error("Clustering cannot be empty.");

	scc_ErrorCode ec;
	scc_Clustering* clustering;
	if ((ec = scc_init_existing_clustering(num_data_points,
	                                       num_clusters,
	                                       INTEGER(R_cluster_labels),
	                                       false,
	                                       &clustering)) != SCC_ER_OK) {
		iRsccwrap_scc_error();
	}

	bool is_OK = false;
	if ((ec = scc_check_clustering(clustering,
	                               size_constraint,
	                               &is_OK)) != SCC_ER_OK) {
		scc_free_clustering(&clustering);
		iRsccwrap_scc_error();
	}

	scc_free_clustering(&clustering);

	return ScalarLogical((int) is_OK);
}


SEXP Rsccwrap_check_clustering_types(const SEXP R_cluster_labels,
                                     const SEXP R_num_clusters,
                                     const SEXP R_size_constraint,
                                     const SEXP R_type_size_constraints,
                                     const SEXP R_type_labels)
{
	if (!isInteger(R_cluster_labels)) iRsccwrap_error("`R_cluster_labels` must be integer.");
	if (!isInteger(R_num_clusters)) iRsccwrap_error("`R_num_clusters` must be integer.");
	if (!isInteger(R_size_constraint)) iRsccwrap_error("`R_size_constraint` must be integer.");
	if (!isInteger(R_type_size_constraints)) iRsccwrap_error("`R_type_size_constraints` must be integer.");
	if (!isInteger(R_type_labels)) iRsccwrap_error("`R_type_labels` must be factor.");

	const uintmax_t num_data_points = (uintmax_t) xlength(R_cluster_labels);
	const uintmax_t num_clusters = (uintmax_t) asInteger(R_num_clusters);
	const uint32_t size_constraint = (uint32_t) asInteger(R_size_constraint);
	const size_t len_type_labels = (size_t) xlength(R_type_labels);
	const int* const type_labels = INTEGER(R_type_labels);

	if (num_clusters == 0) iRsccwrap_error("Clustering cannot be empty.");

	const uintmax_t num_types = (uintmax_t) xlength(R_type_size_constraints);
	uint32_t* const type_size_constraints = (uint32_t*) R_alloc(num_types, sizeof(uint32_t)); // Automatically freed by R on return
	if (type_size_constraints == NULL) iRsccwrap_error("Could not allocate memory.");
	const int* const tmp_type_size_constraints = INTEGER(R_type_size_constraints);
	for (size_t i = 0; i < num_types; ++i) {
		if (tmp_type_size_constraints[i] < 0) iRsccwrap_error("Negative type size constraint.");
		type_size_constraints[i] = (uint32_t) tmp_type_size_constraints[i];
	}

	scc_ErrorCode ec;
	scc_Clustering* clustering;
	if ((ec = scc_init_existing_clustering(num_data_points,
	                                       num_clusters,
	                                       INTEGER(R_cluster_labels),
	                                       false,
	                                       &clustering)) != SCC_ER_OK) {
		iRsccwrap_scc_error();
	}

	bool is_OK = false;
	if ((ec = scc_check_clustering_types(clustering,
	                                     size_constraint,
	                                     num_types,
	                                     type_size_constraints,
	                                     len_type_labels,
	                                     type_labels,
	                                     &is_OK)) != SCC_ER_OK) {
		scc_free_clustering(&clustering);
		iRsccwrap_scc_error();
	}

	scc_free_clustering(&clustering);

	return ScalarLogical((int) is_OK);
}


SEXP Rsccwrap_get_clustering_stats(const SEXP R_cluster_labels,
                                   const SEXP R_num_clusters,
                                   const SEXP R_distance_object)
{
	if (!isInteger(R_cluster_labels)) iRsccwrap_error("`R_cluster_labels` must be integer.");
	if (!isInteger(R_num_clusters)) iRsccwrap_error("`R_num_clusters` must be integer.");
	if (!isMatrix(R_distance_object) || !isReal(R_distance_object)) iRsccwrap_error("Invalid distance object.");

	const uintmax_t num_data_points = (uintmax_t) INTEGER(getAttrib(R_distance_object, R_DimSymbol))[1];
	const uintmax_t num_dimensions = (uintmax_t) INTEGER(getAttrib(R_distance_object, R_DimSymbol))[0];
	const uintmax_t num_clusters = (uintmax_t) asInteger(R_num_clusters);

	if (xlength(R_cluster_labels) != num_data_points) iRsccwrap_error("Clustering does not match distance object.");
	if (num_clusters == 0) iRsccwrap_error("Clustering cannot be empty.");

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

	scc_Clustering* clustering;
	if ((ec = scc_init_existing_clustering(num_data_points,
	                                       num_clusters,
	                                       INTEGER(R_cluster_labels),
	                                       false,
	                                       &clustering)) != SCC_ER_OK) {
		scc_free_data_set_object(&data_set_object);
		iRsccwrap_scc_error();
	}

	scc_ClusteringStats clust_stats;
	if ((ec = scc_get_clustering_stats(clustering,
	                                   data_set_object,
	                                   &clust_stats)) != SCC_ER_OK) {
		scc_free_clustering(&clustering);
		scc_free_data_set_object(&data_set_object);
		iRsccwrap_scc_error();
	}

	scc_free_clustering(&clustering);
	scc_free_data_set_object(&data_set_object);

	if (clust_stats.num_data_points > INT_MAX) iRsccwrap_error("Too many data points.");
	if (clust_stats.num_assigned > INT_MAX) iRsccwrap_error("Too many assigned data points.");
	if (clust_stats.num_clusters > INT_MAX) iRsccwrap_error("Too many clusters.");
	if (clust_stats.num_populated_clusters > INT_MAX) iRsccwrap_error("Too many populated clusters.");
	if (clust_stats.min_cluster_size > INT_MAX) iRsccwrap_error("Too large clusters.");
	if (clust_stats.max_cluster_size > INT_MAX) iRsccwrap_error("Too large clusters.");

	const SEXP R_clust_stats = PROTECT(allocVector(VECSXP, 14));
	SET_VECTOR_ELT(R_clust_stats, 0, ScalarInteger((int) clust_stats.num_data_points));
	SET_VECTOR_ELT(R_clust_stats, 1, ScalarInteger((int) clust_stats.num_assigned));
	SET_VECTOR_ELT(R_clust_stats, 2, ScalarInteger((int) clust_stats.num_clusters));
	SET_VECTOR_ELT(R_clust_stats, 3, ScalarInteger((int) clust_stats.num_populated_clusters));
	SET_VECTOR_ELT(R_clust_stats, 4, ScalarInteger((int) clust_stats.min_cluster_size));
	SET_VECTOR_ELT(R_clust_stats, 5, ScalarInteger((int) clust_stats.max_cluster_size));
	SET_VECTOR_ELT(R_clust_stats, 6, ScalarReal(clust_stats.avg_cluster_size));
	SET_VECTOR_ELT(R_clust_stats, 7, ScalarReal(clust_stats.sum_dists));
	SET_VECTOR_ELT(R_clust_stats, 8, ScalarReal(clust_stats.min_dist));
	SET_VECTOR_ELT(R_clust_stats, 9, ScalarReal(clust_stats.max_dist));
	SET_VECTOR_ELT(R_clust_stats, 10, ScalarReal(clust_stats.cl_avg_min_dist));
	SET_VECTOR_ELT(R_clust_stats, 11, ScalarReal(clust_stats.cl_avg_max_dist));
	SET_VECTOR_ELT(R_clust_stats, 12, ScalarReal(clust_stats.cl_avg_dist_weighted));
	SET_VECTOR_ELT(R_clust_stats, 13, ScalarReal(clust_stats.cl_avg_dist_unweighted));

	const SEXP R_clust_stats_names = PROTECT(allocVector(STRSXP, 14));
	SET_STRING_ELT(R_clust_stats_names, 0, mkChar("num_data_points"));
	SET_STRING_ELT(R_clust_stats_names, 1, mkChar("num_assigned"));
	SET_STRING_ELT(R_clust_stats_names, 2, mkChar("num_clusters"));
	SET_STRING_ELT(R_clust_stats_names, 3, mkChar("num_populated_clusters"));
	SET_STRING_ELT(R_clust_stats_names, 4, mkChar("min_cluster_size"));
	SET_STRING_ELT(R_clust_stats_names, 5, mkChar("max_cluster_size"));
	SET_STRING_ELT(R_clust_stats_names, 6, mkChar("avg_cluster_size"));
	SET_STRING_ELT(R_clust_stats_names, 7, mkChar("sum_dists"));
	SET_STRING_ELT(R_clust_stats_names, 8, mkChar("min_dist"));
	SET_STRING_ELT(R_clust_stats_names, 9, mkChar("max_dist"));
	SET_STRING_ELT(R_clust_stats_names, 10, mkChar("cl_avg_min_dist"));
	SET_STRING_ELT(R_clust_stats_names, 11, mkChar("cl_avg_max_dist"));
	SET_STRING_ELT(R_clust_stats_names, 12, mkChar("cl_avg_dist_weighted"));
	SET_STRING_ELT(R_clust_stats_names, 13, mkChar("cl_avg_dist_unweighted"));
	setAttrib(R_clust_stats, R_NamesSymbol, R_clust_stats_names);

	UNPROTECT(2);
	return R_clust_stats;
}