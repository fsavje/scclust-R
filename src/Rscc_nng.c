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

#include "Rscc_nng.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <R.h>
#include <Rinternals.h>
#include <scclust.h>
#include <scc_data_obj.h>
#include "Rscc_misc.h"


// ==============================================================================
// Internal function prototypes
// ==============================================================================

scc_SeedMethod iRsccwrap_parse_seed_method(SEXP R_seed_method);

scc_UnassignedMethod iRsccwrap_parse_unassigned_method(SEXP R_unassigned_method);


// ==============================================================================
// External function implementations
// ============================================================================== 

SEXP Rsccwrap_nng_clustering(const SEXP R_distance_object,
                             const SEXP R_size_constraint,
                             const SEXP R_seed_method,
                             const SEXP R_main_unassigned_method,
                             const SEXP R_main_radius,
                             const SEXP R_main_data_points,
                             const SEXP R_secondary_unassigned_method,
                             const SEXP R_secondary_radius)
{
	if (!isMatrix(R_distance_object) || !isReal(R_distance_object)) error("Invalid distance object.");
	if (!isInteger(R_size_constraint)) error("`R_size_constraint` must be integer.");
	if (!isString(R_seed_method)) error("`R_seed_method` must be string.");
	if (!isString(R_main_unassigned_method)) error("`R_main_unassigned_method` must be string.");
	if (!isNull(R_main_radius) && !isReal(R_main_radius)) error("`R_main_radius` must be real.");
	if (!isNull(R_main_data_points) && !isLogical(R_main_data_points)) error("`R_target_types` must be logical.");
	if (!isString(R_secondary_unassigned_method)) error("`R_secondary_unassigned_method` must be string.");
	if (!isNull(R_secondary_radius) && !isReal(R_secondary_radius)) error("`R_secondary_radius` must be real.");

	const uintmax_t num_data_points = (uintmax_t) INTEGER(getAttrib(R_distance_object, R_DimSymbol))[1];
	const uintmax_t num_dimensions = (uintmax_t) INTEGER(getAttrib(R_distance_object, R_DimSymbol))[0];
	const uint32_t size_constraint = (uint32_t) asInteger(R_size_constraint);
	const scc_SeedMethod seed_method = iRsccwrap_parse_seed_method(R_seed_method);
	const scc_UnassignedMethod main_unassigned_method = iRsccwrap_parse_unassigned_method(R_main_unassigned_method);
	const scc_UnassignedMethod secondary_unassigned_method = iRsccwrap_parse_unassigned_method(R_secondary_unassigned_method);
	
	size_t len_main_data_points = 0;
	bool* main_data_points = NULL;
	if (isLogical(R_main_data_points)) {
		len_main_data_points = (size_t) xlength(R_main_data_points);
		if (len_main_data_points < num_data_points) {
			error("Invalid `R_main_data_points`.");
		}
		bool* const main_data_points = (bool*) R_alloc(len_main_data_points, sizeof(bool)); // Automatically freed by R on return
		if (main_data_points == NULL) error("Could not allocate memory.");
		const int* const tmp_main_data_points = LOGICAL(R_main_data_points);
		for (size_t i = 0; i < len_main_data_points; ++i) {
			main_data_points[i] = (tmp_main_data_points[i] == 1);
		}
	}

	bool main_radius_constraint = false;
	double main_radius = 0.0;
	if (isReal(R_main_radius)) {
		main_radius_constraint = true;
		main_radius = asReal(R_main_radius);
	}

	bool secondary_radius_constraint = false;
	double secondary_radius = 0.0;
	if (isReal(R_secondary_radius)) {
		secondary_radius_constraint = true;
		secondary_radius = asReal(R_secondary_radius);
	}


	scc_ErrorCode ec;
	scc_DataSetObject* data_set_object;
	if ((ec = scc_get_data_set_object(num_data_points,
	                                  num_dimensions,
	                                  (size_t) xlength(R_distance_object),
	                                  REAL(R_distance_object),
	                                  false,
	                                  &data_set_object)) != SCC_ER_OK) {
		print_error_and_return();
	}

	SEXP R_cluster_labels = PROTECT(allocVector(INTSXP, (R_xlen_t) num_data_points));
	scc_Clustering* clustering;
	if ((ec = scc_init_empty_clustering(num_data_points,
	                                    INTEGER(R_cluster_labels),
	                                    &clustering)) != SCC_ER_OK) {
		scc_free_data_set_object(&data_set_object);
		UNPROTECT(1);
		print_error_and_return();
	}

	if ((ec = scc_nng_clustering(clustering,
	                             data_set_object,
	                             size_constraint,
	                             seed_method,
	                             main_unassigned_method,
	                             main_radius_constraint,
	                             main_radius,
	                             len_main_data_points,
	                             main_data_points,
	                             secondary_unassigned_method,
	                             secondary_radius_constraint,
	                             secondary_radius)) != SCC_ER_OK) {
		scc_free_clustering(&clustering);
		scc_free_data_set_object(&data_set_object);
		UNPROTECT(1);
		print_error_and_return();
	}

	scc_free_data_set_object(&data_set_object);

	uintmax_t num_clusters = 0;
	if ((ec = scc_get_clustering_info(clustering,
	                                  NULL,
	                                  &num_clusters)) != SCC_ER_OK) {
		scc_free_clustering(&clustering);
		UNPROTECT(1);
		print_error_and_return();
	}

	scc_free_clustering(&clustering);

	if (num_clusters > INT_MAX) error("Too many clusters.");
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


// ==============================================================================
// Internal function implementations 
// ==============================================================================

scc_SeedMethod iRsccwrap_parse_seed_method(const SEXP R_seed_method)
{
	if (!isString(R_seed_method)) error("`R_seed_method` must be string.");

	scc_SeedMethod seed_method;
	const char* seed_method_string = CHAR(asChar(R_seed_method));
	if (strcmp(seed_method_string, "lexical") == 0) {
		seed_method = SCC_SM_LEXICAL;
	} else if (strcmp(seed_method_string, "inwards_order") == 0) {
		seed_method = SCC_SM_INWARDS_ORDER;
	} else if (strcmp(seed_method_string, "inwards_updating") == 0) {
		seed_method = SCC_SM_INWARDS_UPDATING;
	} else if (strcmp(seed_method_string, "exclusion_order") == 0) {
		seed_method = SCC_SM_EXCLUSION_ORDER;
	} else if (strcmp(seed_method_string, "exclusion_updating") == 0) {
		seed_method = SCC_SM_EXCLUSION_UPDATING;
	} else {
		error("Not a valid seed method.");
	}

	return seed_method;
}


scc_UnassignedMethod iRsccwrap_parse_unassigned_method(const SEXP R_unassigned_method)
{
	if (!isString(R_unassigned_method)) error("`R_unassigned_method` must be string.");

	scc_UnassignedMethod unassigned_method;
	const char* unassigned_method_string = CHAR(asChar(R_unassigned_method));
	if (strcmp(unassigned_method_string, "ignore") == 0) {
		unassigned_method = SCC_UM_IGNORE;
	} else if (strcmp(unassigned_method_string, "by_nng") == 0) {
		unassigned_method = SCC_UM_ASSIGN_BY_NNG;
	} else if (strcmp(unassigned_method_string, "closest_assigned") == 0) {
		unassigned_method = SCC_UM_CLOSEST_ASSIGNED;
	} else if (strcmp(unassigned_method_string, "closest_seed") == 0) {
		unassigned_method = SCC_UM_CLOSEST_SEED;
	} else if (strcmp(unassigned_method_string, "estimated_radius_closest_seed") == 0) {
		unassigned_method = SCC_UM_CLOSEST_SEED_EST_RADIUS;
	} else {
		error("Not a valid main unassigned method.");
	}

	return unassigned_method;
}
