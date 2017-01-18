# ==============================================================================
# scclust -- A C library for size constrained clustering
# https://github.com/fsavje/scclust
#
# Copyright (C) 2015-2016  Fredrik Savje -- http://fredriksavje.com
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library. If not, see http://www.gnu.org/licenses/
# ==============================================================================

CFLAGS = -std=c99 -O2 -pedantic -Wall -Wextra -Wconversion -Wfloat-equal -Werror
XTRA_FLAGS = -DNDEBUG

LIBDIR = lib
DOCSDIR = doc

OBJECTS = \
	data_set.o \
	digraph_core.o \
	digraph_operations.o \
	dist_search_imp.o \
	error.o \
	hierarchical_clustering.o \
	nng_batch_clustering.o \
	nng_clustering.o \
	nng_core.o \
	nng_findseeds.o \
	scclust_spi.o \
	scclust.o

.PHONY: all clean docs library

all: library

library: $(LIBDIR)/libscclust.a

docs: $(DOCSDIR)
	doxygen DoxyAPI

clean:
	$(RM) -R $(DOCSDIR) $(LIBDIR) src/*.o

$(LIBDIR)/libscclust.a: $(addprefix src/,$(OBJECTS))
	$(AR) rcs $(LIBDIR)/libscclust.a $^

$(LIBDIR)/libscclust.a: | $(LIBDIR)

%.o: %.c
	$(CC) -c $(CFLAGS) $(XTRA_FLAGS) $< -o $@

$(DOCSDIR):
	mkdir -p $(DOCSDIR)

$(LIBDIR):
	mkdir -p $(LIBDIR)
