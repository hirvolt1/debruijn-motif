# to compile, run `make' or `make gvmotif'
#
# gvmotif has support for rendering graphs using graphviz.
# (requires libgraphviz-dev package)

CXXFLAGS = -O3 -std=c++0x -Wall
CXX = g++
GVCFLAGS = `pkg-config libgvc --cflags` -DGVIZ
GVLDFLAGS = `pkg-config libgvc --libs`
EXECUTABLE = motif
GVEXECUTABLE = gvmotif

SRCS = main.cpp similarity.cpp sw_graph.cpp sim_mat.cpp trav.cpp util.cpp gviz.cpp work_for_node_pair.cpp similarity_finder.cpp
OBJS := $(patsubst %.cpp, %.o, $(SRCS))
GVOBJS := $(patsubst %.cpp, gv%.o, $(SRCS))

all: $(SRCS) $(EXECUTABLE)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@
$(EXECUTABLE): $(OBJS)
	$(CXX) $^ -o $@

gv%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(GVCFLAGS) $< -c -o $@
$(GVEXECUTABLE): $(GVOBJS)
	$(CXX) $^ -o $@ $(GVLDFLAGS)

clean:
	rm -f *.o $(EXECUTABLE) $(GVEXECUTABLE)
