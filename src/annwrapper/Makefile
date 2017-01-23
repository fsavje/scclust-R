include $(MAKECONF)

.PHONY : all clean

all : ann_wrapper.o

ann_wrapper.o: ann_wrapper.cpp ../libscclust ../libann
	$(CXX) -c $(ALL_CPPFLAGS) $(ALL_CXXFLAGS) -DNDEBUG -I../libscclust -I../libann $< -o $@

clean :
	$(RM) -f ann_wrapper.o
