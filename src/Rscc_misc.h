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

#ifndef RSCCWRAP_MISC_HG
#define RSCCWRAP_MISC_HG

#include <R.h>
#include <Rinternals.h>
#include <scclust.h>

static inline void print_error_and_return(void) {
	char error_buffer[255];
	scc_get_latest_error(255, error_buffer);
	error(error_buffer);
}

#endif // ifndef RSCCWRAP_MISC_HG
