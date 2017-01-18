#!/bin/sh

rm -rf libann
mkdir libann
wget https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.zip
unzip -q ann_1.1.2.zip
rm ann_1.1.2/src/Makefile
cp -R ann_1.1.2/include ann_1.1.2/src ann_1.1.2/Copyright.txt ann_1.1.2/License.txt libann
rm -rf ann_1.1.2 ann_1.1.2.zip


cat <<EOF > libann/Makefile
include \$(MAKECONF)

LIBOUT = lib/libann.a

LIBOBJ = \\
	ANN.o \\
	brute.o \\
	kd_tree.o \\
	kd_util.o \\
	kd_split.o \\
	kd_dump.o \\
	kd_search.o \\
	kd_pr_search.o \\
	kd_fix_rad_search.o \\
	bd_tree.o \\
	bd_search.o \\
	bd_pr_search.o \\
	bd_fix_rad_search.o \\
	perf.o

.PHONY : all clean

all : \$(LIBOUT)

\$(LIBOUT) : \$(addprefix src/,\$(LIBOBJ))
	mkdir -p lib
	\$(AR) -rcs \$(LIBOUT) \$^

%.o : %.cpp
	\$(CXX) -c \$(ALL_CPPFLAGS) \$(ALL_CXXFLAGS) -DNDEBUG -Iinclude -Wno-unused-const-variable \$< -o \$@

clean :
	\$(RM) -rf lib src/*.o
EOF
