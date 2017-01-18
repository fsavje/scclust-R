include $(MAKECONF)

LIBOUT = lib/libann.a

OBJECTS = \
	ANN.o \
	brute.o \
	kd_tree.o \
	kd_util.o \
	kd_split.o \
	kd_dump.o \
	kd_search.o \
	kd_pr_search.o \
	kd_fix_rad_search.o \
	bd_tree.o \
	bd_search.o \
	bd_pr_search.o \
	bd_fix_rad_search.o \
	perf.o

.PHONY : all clean

all : $(LIBOUT)

$(LIBOUT) : $(addprefix src/,$(OBJECTS))
	mkdir -p lib
	$(AR) -rcs $(LIBOUT) $^

src/%.o : src/%.cpp
	$(CXX) -c $(ALL_CPPFLAGS) $(ALL_CXXFLAGS) -DNDEBUG -Iinclude -Wno-unused-const-variable $< -o $@

clean :
	$(RM) -rf lib src/*.o
