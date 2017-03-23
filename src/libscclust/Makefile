include $(MAKECONF)

# Use 64-bit arc ref: -DSCC_ARC64
# Use stable findseed: -DSCC_STABLE_FINDSEED
# Use stable NNG: -DSCC_STABLE_NNG
XTRA_FLAGS = -DNDEBUG

LIBOUT = lib/libscclust.a

LIBOBJ = \
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
	scclust.o \
	utilities.o

.PHONY : all clean

all : $(LIBOUT)

$(LIBOUT) : $(addprefix src/,$(LIBOBJ))
	mkdir -p lib
	$(AR) -rcs $(LIBOUT) $^

%.o: %.c
	$(CC) -c $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(XTRA_FLAGS) $< -o $@

clean :
	$(RM) -rf lib src/*.o
