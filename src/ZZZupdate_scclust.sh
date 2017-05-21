#!/bin/bash

rm -rf libscclust
mkdir libscclust
wget https://github.com/fsavje/scclust/archive/master.zip
unzip master.zip
cd libscclust
../scclust-master/configure --with-clabel=int --with-clabel-na=min --with-typelabel=int --with-pointindex=int
cd ..


rm libscclust/Makefile
cat <<EOF > libscclust/Makefile
include \$(MAKECONF)

# Use 64-bit arc ref: -DSCC_ARC64
# Use stable findseed: -DSCC_STABLE_FINDSEED
# Use stable NNG: -DSCC_STABLE_NNG
XTRA_FLAGS = -DNDEBUG

LIBOUT = lib/libscclust.a

LIBOBJS = \\
	src/data_set.o \\
	src/digraph_core.o \\
	src/digraph_operations.o \\
	src/dist_search_imp.o \\
	src/error.o \\
	src/hierarchical_clustering.o \\
	src/nng_batch_clustering.o \\
	src/nng_clustering.o \\
	src/nng_core.o \\
	src/nng_findseeds.o \\
	src/scclust_spi.o \\
	src/scclust.o \\
	src/utilities.o

.PHONY: all clean

all: \$(LIBOUT)

\$(LIBOUT): \$(LIBOBJS)
	mkdir -p lib
	\$(AR) -rcs \$(LIBOUT) \$^

%.o: %.c
	\$(CC) -c \$(ALL_CPPFLAGS) \$(ALL_CFLAGS) \$(XTRA_FLAGS) \$< -o \$@

clean:
	\$(RM) -rf lib src/*.o
EOF

rm -r scclust-master master.zip
rm -r libscclust/examples libscclust/DoxyAPI libscclust/doxyrefs.bib libscclust/README.md
