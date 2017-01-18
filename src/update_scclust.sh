#!/bin/bash

rm -rf libscclust
mkdir libscclust
rm -rf annwrapper
mkdir annwrapper
wget https://github.com/fsavje/scclust/archive/master.zip
unzip master.zip
cd libscclust
../scclust-master/configure --with-clabel=int --with-clabel-na=min --with-typelabel=int --with-pointindex=int
cd ..


rm libscclust/Makefile
cat <<EOF > libscclust/Makefile
include \$(MAKECONF)

LIBOUT = lib/libscclust.a

LIBOBJ = \\
	data_set.o \\
	digraph_core.o \\
	digraph_operations.o \\
	dist_search_imp.o \\
	error.o \\
	hierarchical_clustering.o \\
	nng_batch_clustering.o \\
	nng_clustering.o \\
	nng_core.o \\
	nng_findseeds.o \\
	scclust_spi.o \\
	scclust.o

.PHONY : all clean

all : \$(LIBOUT)

\$(LIBOUT) : \$(addprefix src/,\$(LIBOBJ))
	mkdir -p lib
	\$(AR) -rcs \$(LIBOUT) \$^

%.o: %.c
	\$(CC) -c \$(ALL_CPPFLAGS) \$(ALL_CFLAGS) -DNDEBUG \$< -o \$@

clean :
	\$(RM) -rf lib src/*.o
EOF


mv libscclust/examples/ann/ann_wrapper.cpp libscclust/examples/ann/ann_wrapper.h annwrapper/
cat <<EOF > annwrapper/Makefile
include \$(MAKECONF)

.PHONY : all clean

all : ann_wrapper.o

ann_wrapper.o: ann_wrapper.cpp ../libscclust ../libann
	\$(CXX) -c \$(ALL_CPPFLAGS) \$(ALL_CXXFLAGS) -DNDEBUG -I../libscclust -I../libann \$< -o \$@

clean :
	\$(RM) -f ann_wrapper.o
EOF

rm -r scclust-master master.zip
rm -r libscclust/examples libscclust/DoxyAPI libscclust/doxyrefs.bib libscclust/README.md
