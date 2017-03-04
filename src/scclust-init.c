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

#include <R_ext/Rdynload.h>
#include "hierarchical.h"
#include "sc_clustering.h"
#include "utilities.h"


static const R_CallMethodDef callMethods[] = {
	{"Rscc_hierarchical_clustering",  (DL_FUNC) &Rscc_hierarchical_clustering,  4},
	{"Rscc_sc_clustering",            (DL_FUNC) &Rscc_sc_clustering,           12},
	{"Rscc_check_scclust",            (DL_FUNC) &Rscc_check_scclust,            4},
	{"Rscc_get_scclust_stats",        (DL_FUNC) &Rscc_get_scclust_stats,        2},
	{NULL,                            NULL,                                     0}
};


void R_init_scclust(DllInfo *info) {
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
}
