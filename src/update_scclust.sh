#!/bin/bash

rm -r -f libscc/*
mkdir -p libscc/exlib/libANN
wget https://github.com/fsavje/scclust/archive/master.zip
unzip master.zip
cp -r scclust-master/include scclust-master/src libscc
cp -r scclust-master/exlib/libANN/include scclust-master/exlib/libANN/src libscc/exlib/libANN
rm -r scclust-master master.zip

sed -i '40,80s/typedef uint32_t scc_Clabel;/typedef int scc_Clabel;/' libscc/include/scclust.h
sed -i '40,80s/static const scc_Clabel SCC_CLABEL_MAX = UINT32_MAX;/static const scc_Clabel SCC_CLABEL_MAX = INT_MAX;/' libscc/include/scclust.h
sed -i '40,80s/static const scc_Clabel SCC_CLABEL_NA = UINT32_MAX;/static const scc_Clabel SCC_CLABEL_NA = INT_MIN;/' libscc/include/scclust.h
sed -i '40,80s/typedef uint_fast16_t scc_TypeLabel;/typedef int scc_TypeLabel;/' libscc/include/scclust.h


cat <<EOF > libscc/Makefile
include \$(MAKECONF)

LIBOUT = libscc.a
SCC_OBJ = \$(addprefix src/,digraph_core.o digraph_operations.o dist_nnsearch_ANN.o dist_search.o error.o \\
                           hierarchical_clustering.o nng_clustering.o nng_batch_clustering.o \\
                           nng_core.o nng_findseeds.o scc_data_set.o scclust.o) \\
          \$(addprefix exlib/libANN/src/,ANN.o brute.o kd_tree.o kd_util.o kd_split.o kd_dump.o \\
                                        kd_search.o kd_pr_search.o kd_fix_rad_search.o bd_tree.o \\
                                        bd_search.o bd_pr_search.o bd_fix_rad_search.o perf.o)

# Use 64-bit arc ref: -DSCC_ARC64
# Use BD-tree: -DSCC_ANN_BDTREE
# Use stable findseed: -DSCC_STABLE_FINDSEED
# Use stable NNG: -DSCC_STABLE_NNG
XTRA_FLAGS = -DNDEBUG -DSCC_DPID_INT

.PHONY: clean library

library: \$(LIBOUT)

\$(LIBOUT): \$(SCC_OBJ)
	\$(AR) -rcs \$(LIBOUT) \$^

src/%.o: src/%.c
	\$(CC) -c \$(ALL_CPPFLAGS) \$(ALL_CFLAGS) \$(XTRA_FLAGS) \$< -o \$@

src/%.o: src/%.cpp
	\$(CXX) -c \$(ALL_CPPFLAGS) \$(ALL_CXXFLAGS) \$(XTRA_FLAGS) \$< -o \$@

exlib/libANN/src/%.o: exlib/libANN/src/%.cpp
	\$(CXX) -c \$(ALL_CPPFLAGS) \$(ALL_CXXFLAGS) -DNDEBUG -Iexlib/libANN/include -Wno-unused-const-variable \$< -o \$@

clean:
	\$(RM) -f \$(LIBOUT) src/*.o exlib/libANN/src/*.o
EOF
