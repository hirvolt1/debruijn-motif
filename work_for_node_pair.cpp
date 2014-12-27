#include <iostream>
#include <string>
#include "work_for_node_pair.hpp"
#include "to_string.hpp"

//#define db_weight_adjust
//#define db_adj_ic_adjust
//#define db_verbose_high

int get_vector_map_size(map<int, ic_edge_cont> &tmp_ic_edges) {
  if (tmp_ic_edges.size() == 0) {
    return 0;
  }
  int size_sum = 0;
  for (map<int, ic_edge_cont>::iterator it = tmp_ic_edges.begin();
       it != tmp_ic_edges.end();
       it++) {
    size_sum += it->second.size();
  }
  return size_sum;
}

void cout_vector_map(map<int, sim_type> &max_sims) {
  if (max_sims.size() == 0) {
    std::cout << "map empty" << std::endl;
    return;
  }

  std::cout << "Map:" << std::endl;
  for (map<int, sim_type>::iterator it = max_sims.begin();
       it != max_sims.end();
       it++) {
    std::cout << it->first << ":" << it->second << std::endl;
  }
  return;
}




void cout_vector_map(map<int, ic_edge_cont> &tmp_ic_edges) {
  if (tmp_ic_edges.size() == 0) {
    std::cout << "map empty" << std::endl;
    return;
  }

  std::cout << "Map:" << std::endl;
  for (map<int, ic_edge_cont>::iterator it = tmp_ic_edges.begin();
       it != tmp_ic_edges.end();
       it++) {
    std::cout << it->first << ":" << it->second.size();
    for (ic_edge_cont::iterator it2 = it->second.begin();
         it2 != it->second.end();
         it2++) {
      std::cout << " " << get_identifying_string(*it2);
    }
    std::cout << std::endl;
  }
  return;
}


inline sim_type similarity_real(const char *a, const char *b) {
  sim_type score;
  sim_type sim = 0;
  for(unsigned int i = 0; i < g_subword_length; i++) {
    score = g_scoring_matrix[a[i]-'A'][b[i]-'A'];  // implicit type conversion
    sim += score;
  }
  return sim;
}

void do_work_for_node_pair_by_pointers(sw_node *node1_p, sw_node *node2_p,
                                       sim_type sim_value) {

  bool add_ic_edge_for1 = false;
  //bool add_ic_edge_for2 = false;
  bool remove_old_from1 = false;
  //bool remove_old_from2 = false;
  bool first_for_this_seq1 = false;
  //bool first_for_this_seq2 = false;
  //NOTE this is already checked before calling this function, heh
  //if(node1_p->seq != node2_p->seq) { // Not in the same sequence



  sim_value = similarity_real(node1_p->subword, node2_p->subword);
#ifdef db_weight_adjust
  std::cout << "n1.w: " << node1_p->weight << ", n2.w: " << node2_p->weight
      << ", sim: " << sim_value << std::endl;
#endif
  node1_p->weight += sim_value;
  //node2_p->weight += sim_value;
#ifdef db_weight_adjust
  std::cout << "n1.w: " << node1_p->weight << ", n2.w: " << node2_p->weight
      << ", sim: " << sim_value << std::endl;
#endif


#ifdef db_verbose_high
  std::string n1w = to_string(node1_p->subword).substr(0,4);
  std::string n2w = to_string(node2_p->subword).substr(0,4);

  std::cout << "pre : n1: " << node1_p->seq << ":" << n1w << ";" 
      << get_vector_map_size(node1_p->tmp_ic_edges) << ";" << node1_p->gen_mul
      << std::endl << "pre : n2: " << node2_p->seq << ":" << n2w << ";"
      << get_vector_map_size(node2_p->tmp_ic_edges) << ";" << node2_p->gen_mul
      << std::endl;

  std::cout << "n1 map... ";
  cout_vector_map(node1_p->tmp_ic_edges);
  std::cout << "... maxes:";
  cout_vector_map(node1_p->max_sims);
  std::cout << "n2 map... ";
  cout_vector_map(node2_p->tmp_ic_edges);
  std::cout << "... maxes:";
  cout_vector_map(node2_p->max_sims);
  std::cout << "sim: " << sim_value << std::endl;

#endif
  // Update max_sims for node1 if needed
  if ( node1_p->max_sims.find(node2_p->seq) == node1_p->max_sims.end() ) {
    node1_p->max_sims[node2_p->seq] = sim_value;
    first_for_this_seq1 = true;
    add_ic_edge_for1 = true;
  } else if ( node1_p->max_sims[node2_p->seq] < sim_value) {
    node1_p->max_sims[node2_p->seq] = sim_value;
    add_ic_edge_for1 = true;
    remove_old_from1 = true;
  } else if ( node1_p->max_sims[node2_p->seq] == sim_value) { //more than 1 
    add_ic_edge_for1 = true;
  }

  // Update max_sims for node2 if needed
  /*if ( node2_p->max_sims.find(node1_p->seq) == node2_p->max_sims.end() ) {
    node2_p->max_sims[node1_p->seq] = sim_value;
    first_for_this_seq2 = true;
    add_ic_edge_for2 = true;
    } else if ( node2_p->max_sims[node1_p->seq] < sim_value) {
    node2_p->max_sims[node1_p->seq] = sim_value;
    add_ic_edge_for2 = true;
    remove_old_from2 = true;
    }*/

  // Update tmp_ic_edges for node1 if needed
  if (add_ic_edge_for1) {
#ifdef db_adj_ic_adjust
    std::cout << "attempting pushing: " << get_identifying_string(node2_p);
    std::cout << ", before pushing map state is: ";
    cout_vector_map(node1_p->tmp_ic_edges);
#endif
    if (remove_old_from1) {  // ok, not really necessary
      (node1_p->tmp_ic_edges[node2_p->seq]).clear();
    }
    if (node1_p->tmp_ic_edges.find(node2_p->seq) 
        == node1_p->tmp_ic_edges.end()) {
      ic_edge_cont iec_tmp;
      iec_tmp.push_back(node2_p);
      node1_p->tmp_ic_edges[node2_p->seq] = iec_tmp;
    } else {
      (node1_p->tmp_ic_edges[node2_p->seq]).push_back(node2_p);
    }
#ifdef db_adj_ic_adjust
    std::cout << "attempted pushing: " << get_identifying_string(node2_p);
    std::cout << ", resulting map: ";
    cout_vector_map(node1_p->tmp_ic_edges);
#endif

  }

  // Update tmp_ic_edges for node2 if needed
  /*if (add_ic_edge_for2) {
    if (remove_old_from2) {  // ok, not really necessary
    (node2_p->tmp_ic_edges[node1_p->seq]).clear();
    }
    if (node2_p->tmp_ic_edges.find(node1_p->seq) 
    == node2_p->tmp_ic_edges.end()) {
    ic_edge_cont iec_tmp;
    iec_tmp.push_back(node1_p);
    node2_p->tmp_ic_edges[node1_p->seq] = iec_tmp;
    } else {
    (node2_p->tmp_ic_edges[node1_p->seq]).push_back(node1_p);
    }
    }*/

  // increment general multiplicity if first edges assigned for these seqs
  if (first_for_this_seq1) {
    node1_p->gen_mul++;
  }

  /*if (first_for_this_seq2) {
    node2_p->gen_mul++;
    }*/
#ifdef db_adj_ic_adjust
  int size_sum = 0;
  for (map<int, ic_edge_cont>::iterator it = node2_p->tmp_ic_edges.begin();
       it != node2_p->tmp_ic_edges.end();
       it++) {
    //for (ice_edge_cont::iterator it2 = it->second.begin();
    //    it2 != it->second.end();
    //    it2++) {

    size_sum += it->second.size();
  }
  std::cout << "n2 tmp ics: " << size_sum << ", gm: "<< node2_p->gen_mul 
      << std::endl;

#endif
#ifdef db_verbose_high
  //std::string n1w(node1_p->subword, node1_p->seq + 4);
  //std::string n2w(node2_p->subword, node2_p->seq + 4);
  std::cout << "post : n1: " << node1_p->seq << ":" << n1w << ";" 
      << get_vector_map_size(node1_p->tmp_ic_edges) << ";" << node1_p->gen_mul
      << std::endl << "post : n2: " << node2_p->seq << ":" << n2w << ";"
      << get_vector_map_size(node2_p->tmp_ic_edges) << ";" << node2_p->gen_mul
      << std::endl;

  std::cout << "n1 map... ";
  cout_vector_map(node1_p->tmp_ic_edges);
  std::cout << "... maxes:";
  cout_vector_map(node1_p->max_sims);
  std::cout << "n2 map... ";
  cout_vector_map(node2_p->tmp_ic_edges);
  std::cout << "... maxes:";
  cout_vector_map(node2_p->max_sims);


#endif

  //}
}

//void adjust_weights(sw_node *node1_p, sim_type to_deduce) {
//  node1_p->weight -= sim_value;
//}

