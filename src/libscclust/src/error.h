/* =============================================================================
 * scclust -- A C library for size-constrained clustering
 * https://github.com/fsavje/scclust
 *
 * Copyright (C) 2015-2017  Fredrik Savje -- http://fredriksavje.com
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
 * ========================================================================== */

#ifndef SCC_ERROR_HG
#define SCC_ERROR_HG

#include "../include/scclust.h"


// =============================================================================
// Macros
// =============================================================================

#define iscc_make_error(ec) iscc_make_error__(ec, NULL, __FILE__, __LINE__)

#define iscc_make_error_msg(ec, msg) iscc_make_error__(ec, msg, __FILE__, __LINE__)

#define iscc_no_error() (SCC_ER_OK)


// =============================================================================
// Function prototypes
// =============================================================================

scc_ErrorCode iscc_make_error__(scc_ErrorCode ec,
                                const char* msg,
                                const char* file,
                                int line);


void iscc_reset_error(void);


#endif // ifndef SCC_ERROR_HG
