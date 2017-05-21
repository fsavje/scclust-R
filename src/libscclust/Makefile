include $(MAKECONF)

# Use 64-bit arc ref: -DSCC_ARC64
# Use stable findseed: -DSCC_STABLE_FINDSEED
# Use stable NNG: -DSCC_STABLE_NNG
XTRA_FLAGS = -DNDEBUG

LIBOUT = lib/libscclust.a

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

.PHONY: all clean

all: $(LIBOUT)

$(LIBOUT): $(LIBOBJS)
	mkdir -p lib
	$(AR) -rcs $(LIBOUT) $^

%.o: %.c
	$(CC) -c $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(XTRA_FLAGS) $< -o $@

clean:
	$(RM) -rf lib src/*.o
