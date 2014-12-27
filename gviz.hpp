#ifndef GVIZ_HPP_INCLUDED
#define GVIZ_HPP_INCLUDED
#ifdef GVIZ

#include "sw_graph.hpp"

void render_graph(char **seqs, int seq_count, sw_graph<> *G, char *out);
std::string get_node_name(sw_node *node);

#endif // #ifdef GVIZ
#endif
