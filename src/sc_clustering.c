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

#include "sc_clustering.h"
#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <R.h>
#include <Rinternals.h>
#include <scclust.h>
#include "error.h"
#include "internal.h"


// =============================================================================
// Internal function prototypes
// =============================================================================

static scc_SeedMethod iRscc_parse_seed_method(SEXP R_seed_method);

static scc_UnassignedMethod iRscc_parse_unassigned_method(SEXP R_unassigned_method);


// =============================================================================
// External function implementations
// =============================================================================

SEXP Rscc_sc_clustering(const SEXP R_distances,
                        const SEXP R_size_constraint,
                        const SEXP R_type_labels,
                        const SEXP R_type_constraints,
                        const SEXP R_seed_method,
                        const SEXP R_primary_data_points,
                        const SEXP R_primary_unassigned_method,
                        const SEXP R_secondary_unassigned_method,
                        const SEXP R_seed_radius,
                        const SEXP R_primary_radius,
                        const SEXP R_secondary_radius,
                        const SEXP R_batch_size)
{
	Rscc_set_dist_functions();

	if (!idist_check_distance_object(R_distances)) {
		iRscc_error("`R_distances` is not a valid distance object.");
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
	if (!isString(R_seed_method)) {
		iRscc_error("`R_seed_method` must be string.");
	}
	if (!isNull(R_primary_data_points) && !isInteger(R_primary_data_points)) {
		iRscc_error("`R_primary_data_points` must be NULL or integer.");
	}
	if (!isString(R_primary_unassigned_method)) {
		iRscc_error("`R_primary_unassigned_method` must be string.");
	}
	if (!isString(R_secondary_unassigned_method)) {
		iRscc_error("`R_secondary_unassigned_method` must be string.");
	}
	if (!isNull(R_seed_radius) && !isReal(R_seed_radius)) {
		iRscc_error("`R_seed_radius` must be NULL or double.");
	}
	if (!isNull(R_primary_radius) && !isString(R_primary_radius) && !isReal(R_primary_radius)) {
		iRscc_error("`R_primary_radius` must be NULL, string or double.");
	}
	if (!isNull(R_secondary_radius) && !isString(R_secondary_radius) && !isReal(R_secondary_radius)) {
		iRscc_error("`R_secondary_radius` must be NULL, string or double.");
	}
	if (!isNull(R_batch_size) && !isInteger(R_batch_size)) {
		iRscc_error("`R_batch_size` must be NULL or integer.");
	}

	const uintmax_t num_data_points = (uintmax_t) idist_num_data_points(R_distances);

	scc_ClusterOptions options = scc_default_cluster_options;

	options.size_constraint = (uint32_t) asInteger(R_size_constraint);
	options.seed_method = iRscc_parse_seed_method(R_seed_method);
	options.primary_unassigned_method = iRscc_parse_unassigned_method(R_primary_unassigned_method);
	options.secondary_unassigned_method = iRscc_parse_unassigned_method(R_secondary_unassigned_method);

	if (isInteger(R_type_labels) && isInteger(R_type_constraints)) {
		const uintmax_t num_types = (uintmax_t) xlength(R_type_constraints);
		const size_t len_type_labels = (size_t) xlength(R_type_labels);
		if (len_type_labels != num_data_points) {
			iRscc_error("`R_type_labels` does not match `R_distances`.");
		}
		if (num_types >= 2) {
			uint32_t* const type_constraints = (uint32_t*) R_alloc(num_types, sizeof(uint32_t)); // Automatically freed by R on return
			if (type_constraints == NULL) iRscc_error("Could not allocate memory.");
			const int* const tmp_type_constraints = INTEGER(R_type_constraints);
			for (size_t i = 0; i < num_types; ++i) {
				if (tmp_type_constraints[i] < 0) {
				  iRscc_error("Negative type size constraint.");
				}
				type_constraints[i] = (uint32_t) tmp_type_constraints[i];
			}

			options.num_types = num_types;
			options.type_constraints = type_constraints;
			options.len_type_labels = len_type_labels;
			options.type_labels = INTEGER(R_type_labels);
		}
	}

	if (isInteger(R_primary_data_points)) {
		options.len_primary_data_points = (size_t) xlength(R_primary_data_points);
		options.primary_data_points = INTEGER(R_primary_data_points);
	}

	if (isReal(R_seed_radius)) {
		options.seed_radius = SCC_RM_USE_SUPPLIED;
		options.seed_supplied_radius = asReal(R_seed_radius);
	} else if (isNull(R_seed_radius)) {
		options.seed_radius = SCC_RM_NO_RADIUS;
	} else if (isString(R_seed_radius)) {
		if (strcmp(CHAR(asChar(R_seed_radius)), "no_radius") == 0) {
			options.seed_radius = SCC_RM_NO_RADIUS;
		} else {
			iRscc_error("Not a valid radius method.");
		}
	}

	if (isReal(R_primary_radius)) {
		options.primary_radius = SCC_RM_USE_SUPPLIED;
		options.primary_supplied_radius = asReal(R_primary_radius);
	} else if (isNull(R_primary_radius)) {
		options.primary_radius = SCC_RM_NO_RADIUS;
	} else if (isString(R_primary_radius)) {
		if (strcmp(CHAR(asChar(R_primary_radius)), "no_radius") == 0) {
			options.primary_radius = SCC_RM_NO_RADIUS;
		} else if (strcmp(CHAR(asChar(R_primary_radius)), "seed_radius") == 0) {
			options.primary_radius = SCC_RM_USE_SEED_RADIUS;
		} else if (strcmp(CHAR(asChar(R_primary_radius)), "estimated_radius") == 0) {
			options.primary_radius = SCC_RM_USE_ESTIMATED;
		} else {
			iRscc_error("Not a valid radius method.");
		}
	}

	if (isReal(R_secondary_radius)) {
		options.secondary_radius = SCC_RM_USE_SUPPLIED;
		options.secondary_supplied_radius = asReal(R_secondary_radius);
	} else if (isNull(R_secondary_radius)) {
		options.secondary_radius = SCC_RM_NO_RADIUS;
	} else if (isString(R_secondary_radius)) {
		if (strcmp(CHAR(asChar(R_secondary_radius)), "no_radius") == 0) {
			options.secondary_radius = SCC_RM_NO_RADIUS;
		} else if (strcmp(CHAR(asChar(R_secondary_radius)), "seed_radius") == 0) {
			options.secondary_radius = SCC_RM_USE_SEED_RADIUS;
		} else if (strcmp(CHAR(asChar(R_secondary_radius)), "estimated_radius") == 0) {
			options.secondary_radius = SCC_RM_USE_ESTIMATED;
		} else {
			iRscc_error("Not a valid radius method.");
		}
	}

	if (isInteger(R_batch_size)) {
		options.batch_size = (uint32_t) asInteger(R_batch_size);
	}

	scc_ErrorCode ec;
	SEXP R_cluster_labels = PROTECT(allocVector(INTSXP, (R_xlen_t) num_data_points));
	scc_Clustering* clustering;
	if ((ec = scc_init_empty_clustering(num_data_points,
	                                    INTEGER(R_cluster_labels),
	                                    &clustering)) != SCC_ER_OK) {
		iRscc_scc_error();
	}

	if ((ec = scc_make_clustering(Rscc_get_distances_pointer(R_distances),
	                              clustering,
	                              &options)) != SCC_ER_OK) {
		scc_free_clustering(&clustering);
		iRscc_scc_error();
	}

	uintmax_t num_clusters = 0;
	if ((ec = scc_get_clustering_info(clustering,
	                                  NULL,
	                                  &num_clusters)) != SCC_ER_OK) {
		scc_free_clustering(&clustering);
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


// =============================================================================
// Internal function implementations
// =============================================================================

static scc_SeedMethod iRscc_parse_seed_method(const SEXP R_seed_method)
{
	if (!isString(R_seed_method)) iRscc_error("`R_seed_method` must be string.");

	const char* seed_method_string = CHAR(asChar(R_seed_method));
	if (strcmp(seed_method_string, "lexical") == 0) {
		return SCC_SM_LEXICAL;
	} else if (strcmp(seed_method_string, "batches") == 0) {
		return SCC_SM_BATCHES;
	} else if (strcmp(seed_method_string, "inwards_order") == 0) {
		return SCC_SM_INWARDS_ORDER;
	} else if (strcmp(seed_method_string, "inwards_updating") == 0) {
		return SCC_SM_INWARDS_UPDATING;
	} else if (strcmp(seed_method_string, "inwards_alt_updating") == 0) {
		return SCC_SM_INWARDS_ALT_UPDATING;
	} else if (strcmp(seed_method_string, "exclusion_order") == 0) {
		return SCC_SM_EXCLUSION_ORDER;
	} else if (strcmp(seed_method_string, "exclusion_updating") == 0) {
		return SCC_SM_EXCLUSION_UPDATING;
	} else {
		iRscc_error("Not a valid seed method.");
	}

	return 999; // Unreachable, but needed to silence compiler warning
}


static scc_UnassignedMethod iRscc_parse_unassigned_method(const SEXP R_unassigned_method)
{
	if (!isString(R_unassigned_method)) iRscc_error("`R_unassigned_method` must be string.");

	const char* unassigned_method_string = CHAR(asChar(R_unassigned_method));
	if (strcmp(unassigned_method_string, "ignore") == 0) {
		return SCC_UM_IGNORE;
	} else if (strcmp(unassigned_method_string, "any_neighbor") == 0) {
		return SCC_UM_ANY_NEIGHBOR;
	} else if (strcmp(unassigned_method_string, "closest_assigned") == 0) {
		return SCC_UM_CLOSEST_ASSIGNED;
	} else if (strcmp(unassigned_method_string, "closest_seed") == 0) {
		return SCC_UM_CLOSEST_SEED;
	} else {
		iRscc_error("Not a valid unassigned method.");
	}

	return 999; // Unreachable, but needed to silence compiler warning
}
