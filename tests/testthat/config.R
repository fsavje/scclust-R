# ==============================================================================
# Rscclust -- R wrapper for the scclust library
# https://github.com/fsavje/Rscclust
#
# Copyright (C) 2016  Fredrik Savje -- http://fredriksavje.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/
# ==============================================================================

if (identical(Sys.getenv("TRAVIS"), "true")) {
  run_slow_tests <- TRUE
} else if (identical(Sys.getenv("APPVEYOR"), "True")) {
  run_slow_tests <- TRUE
} else {
  run_slow_tests <- FALSE
}

compiled_with_stable_nng <- FALSE

compiled_with_stable_findseed <- FALSE
