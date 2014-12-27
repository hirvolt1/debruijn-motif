#ifndef DEBRUIJN_PMOTIF_DO_WORK_FOR_NODE_PAIR_H_
#define DEBRUIJN_PMOTIF_DO_WORK_FOR_NODE_PAIR_H_

#include <string>
#include "sw_graph.hpp"
#include "to_string.hpp"


void do_work_for_node_pair(sw_node &node1, sw_node &node2, sim_type sim_value);
void do_work_for_node_pair_by_pointers(sw_node *node1_p, sw_node *node2_p,
                                       sim_type sim_value);

inline std::string get_identifying_string(sw_node* node) {
  return to_string(node->seq).append(":").append(to_string(node->subword).substr(0,g_subword_length));
}


#endif  // DEBRUIJN_PMOTIF_DO_WORK_FOR_NODE_PAIR_H_
