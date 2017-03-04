/* =============================================================================
 * scclust for R -- R wrapper for the scclust library
 * https://github.com/fsavje/scclust-R
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

#include "utilities.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <R.h>
#include <Rinternals.h>
#include <scclust.h>
#include "error.h"
#include "internal.h"


// =============================================================================
// External function implementations
// =============================================================================

SEXP Rscc_check_clustering(const SEXP R_clustering,
                           const SEXP R_size_constraint,
                           const SEXP R_type_labels,
                           const SEXP R_type_constraints)
{
	if (!isInteger(R_clustering)) {
		iRscc_error("`R_clustering` is not a valid clustering object.");
	}
	if (!isInteger(getAttrib(R_clustering, install("cluster_count")))) {
		iRscc_error("`R_clustering` is not a valid clustering object.");
	}
	if (!isInteger(R_size_constraint)) {
		iRscc_error("`R_size_constraint` must be integer.");
	}
	if (isNull(R_type_labels)) {
		if (!isNull(R_type_constraints)) {
			iRscc_error("`R_type_constraints` must be NULL when no types are supplied.");
		}
	} else {
		if (!isInteger(R_type_labels)) {
			iRscc_error("`R_type_labels` must be factor, integer or NULL.");
		}
		if (!isInteger(R_type_constraints)) {
			iRscc_error("`R_type_constraints` must be integer.");
		}
	}

	const uintmax_t num_data_points = (uintmax_t) xlength(R_clustering);
	const uintmax_t num_clusters = (uintmax_t) asInteger(getAttrib(R_clustering, install("cluster_count")));
	if (num_clusters == 0) {
		iRscc_error("`R_clustering` is empty.");
	}

	const uint32_t size_constraint = (uint32_t) asInteger(R_size_constraint);
   uintmax_t num_types = 0;
   uint32_t* type_constraints = NULL;
   size_t len_type_labels = 0;
   const int* type_labels = NULL;

	if (isInteger(R_type_labels) && isInteger(R_type_constraints)) {
		num_types = (uintmax_t) xlength(R_type_constraints);
		len_type_labels = (size_t) xlength(R_type_labels);
		if (len_type_labels != num_data_points) {
			iRscc_error("`R_type_labels` does not match `R_clustering`.");
		}
		if (num_types >= 2) {
			type_constraints = (uint32_t*) R_alloc(num_types, sizeof(uint32_t)); // Automatically freed by R on return
			if (type_constraints == NULL) iRscc_error("Could not allocate memory.");
			const int* const tmp_type_constraints = INTEGER(R_type_constraints);
			for (size_t i = 0; i < num_types; ++i) {
				if (tmp_type_constraints[i] < 0) {
				  iRscc_error("Negative type size constraint.");
				}
				type_constraints[i] = (uint32_t) tmp_type_constraints[i];
			}
			type_labels = INTEGER(R_type_labels);
		}
	}

	scc_ErrorCode ec;
	scc_Clustering* clustering;
	if ((ec = scc_init_existing_clustering(num_data_points,
	                                       num_clusters,
	                                       INTEGER(R_clustering),
	                                       false,
	                                       &clustering)) != SCC_ER_OK) {
		iRscc_scc_error();
	}

	bool is_OK = false;
	if ((ec = scc_check_clustering(clustering,
	                               size_constraint,
	                               num_types,
	                               type_constraints,
	                               len_type_labels,
	                               type_labels,
	                               &is_OK)) != SCC_ER_OK) {
		scc_free_clustering(&clustering);
		iRscc_scc_error();
	}

	scc_free_clustering(&clustering);

	return ScalarLogical((int) is_OK);
}


SEXP Rscc_get_scclust_stats(const SEXP R_clustering,
                            const SEXP R_distances)
{
	Rscc_set_dist_functions();

	if (!isInteger(R_clustering)) {
		iRscc_error("`R_clustering` is not a valid clustering object.");
	}
	if (!isInteger(getAttrib(R_clustering, install("cluster_count")))) {
		iRscc_error("`R_clustering` is not a valid clustering object.");
	}
	if (!idist_check_distance_object(R_distances)) {
		iRscc_error("`R_distances` is not a valid distance object.");
	}

	const uintmax_t num_data_points = (uintmax_t) idist_num_data_points(R_distances);
	const uintmax_t num_clusters = (uintmax_t) asInteger(getAttrib(R_clustering, install("cluster_count")));

	if (((uintmax_t) xlength(R_clustering)) != num_data_points) {
		iRscc_error("`R_distances` does not match `R_clustering`.");
	}
	if (num_clusters == 0) {
		iRscc_error("`R_clustering` is empty.");
	}

	scc_ErrorCode ec;
	scc_Clustering* clustering;
	if ((ec = scc_init_existing_clustering(num_data_points,
	                                       num_clusters,
	                                       INTEGER(R_clustering),
	                                       false,
	                                       &clustering)) != SCC_ER_OK) {
		iRscc_scc_error();
	}

	scc_ClusteringStats clust_stats;
	if ((ec = scc_get_clustering_stats(clustering,
	                                   Rscc_get_distances_pointer(R_distances),
	                                   &clust_stats)) != SCC_ER_OK) {
		scc_free_clustering(&clustering);
		iRscc_scc_error();
	}

	scc_free_clustering(&clustering);

	if (clust_stats.num_data_points > INT_MAX) iRscc_error("Too many data points.");
	if (clust_stats.num_assigned > INT_MAX) iRscc_error("Too many assigned data points.");
	if (clust_stats.num_clusters > INT_MAX) iRscc_error("Too many clusters.");
	if (clust_stats.num_populated_clusters > INT_MAX) iRscc_error("Too many populated clusters.");
	if (clust_stats.min_cluster_size > INT_MAX) iRscc_error("Too large clusters.");
	if (clust_stats.max_cluster_size > INT_MAX) iRscc_error("Too large clusters.");

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
