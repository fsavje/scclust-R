# Use 64-bit arc ref: -DSCC_ARC64
# Use stable NNG: -DSCC_STABLE_NNG
# Use stable findseed: -DSCC_STABLE_FINDSEED
XTRA_FLAGS =

LIBOBJS = \
	src/data_set.o \
	src/digraph_core.o \
	src/digraph_operations.o \
	src/dist_search_imp.o \
	src/error.o \
	src/hierarchical_clustering.o \
	src/nng_batch_clustering.o \
	src/nng_clustering.o \
	src/nng_core.o \
	src/nng_findseeds.o \
	src/scclust_spi.o \
	src/scclust.o \
	src/utilities.o

libscclust.a: $(LIBOBJS)
	$(R_AR) -rcs libscclust.a $^

%.o: %.c
	$(R_CC) $(R_CPPFLAGS) $(R_CFLAGS) $(XTRA_FLAGS) -c $< -o $@

clean:
	$(R_RM) libscclust.a $(LIBOBJS)

.PHONY: clean
