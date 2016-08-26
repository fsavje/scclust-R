include $(MAKECONF)

LIBOUT = libscc.a
OBJECTS = $(addprefix src/,digraph_core.o digraph_operations.o dist_nnsearch_ANN.o dist_search.o error.o \
                           greedy_bottom_clustering.o greedy_top_clustering.o nng_clustering.o \
                           nng_batch_clustering.o nng_core.o nng_findseeds.o scc_data_obj.o scclust.o) \
          $(addprefix exlib/libANN/src/,ANN.o brute.o kd_tree.o kd_util.o kd_split.o kd_dump.o \
                                        kd_search.o kd_pr_search.o kd_fix_rad_search.o bd_tree.o \
                                        bd_search.o bd_pr_search.o bd_fix_rad_search.o perf.o)

# Use 64-bit arc ref "-DSCC_ARC64"
# Use BD-tree "-DSCC_ANN_BDTREE"
XTRA_FLAGS = -DNDEBUG -DSCC_DPID_INT

.PHONY: clean library

library: $(LIBOUT)

$(LIBOUT): $(OBJECTS)
	$(AR) -rcs $(LIBOUT) $^

src/%.o: src/%.c
	$(CC) -c $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(XTRA_FLAGS) $< -o $@

src/%.o: src/%.cpp
	$(CXX) -c $(ALL_CPPFLAGS) $(ALL_CXXFLAGS) $(XTRA_FLAGS) $< -o $@

exlib/libANN/src/%.o: exlib/libANN/src/%.cpp
	$(CXX) -c $(ALL_CPPFLAGS) $(ALL_CXXFLAGS) -DNDEBUG -Iexlib/libANN/include -Wno-unused-const-variable $< -o $@

clean:
	$(RM) -f $(LIBOUT) src/*.o exlib/libANN/src/*.o
