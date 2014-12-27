#include <iostream>
#include <algorithm>
#include <string>
#include <map>

#include "sw_graph.hpp"
#include "sim_mat.hpp"
#include "util.hpp"
#include "similarity_finder.hpp"
#include "work_for_node_pair.hpp"

//#define SIM_TYPE_MIN LONG_MIN
#define depth2

//#define db_cout_similar_nodes

//#define db_cout_both_similarity_finder_processes

//#define db_weight


#ifdef db_cout_both_similarity_finder_processes
int g_traditional_pair_counts = 0;
#endif


using std::string;

void 
init_map_amino_acids(map<char,int> &background_map)
{
  string amino_acids = "ACDEFGHIKLMNPQRSTVWYZ";
  size_t amino_acids_size = amino_acids.size();
  for (size_t i = 0; i < amino_acids_size; i++)
  {
    background_map[amino_acids[i]] = 0;
  }
}

/*vector<long double> 
  map_to_vector_amino_acids(map<char,int> &background_map)
  {
  string amino_acids = "ACDEFGHIKLMNPQRSTVWYZ";
  size_t amino_acids_size = amino_acids.size();
  vector<long double> out(amino_acids_size);
  for (int i = 0; i < amino_acids_size; i++)
  {
  out[i] = background_map[amino_acids[i]];
  background_map[amino_acids[i]] = 0;
  }
  }*/

sw_graph<node_cont>* build_debruijn_graph(char **seqs, int seq_count,
                                          SIM_mat *SIM, map<char,int> *background) {
  //map<char,int> background_map;
  init_map_amino_acids(*background);

  sw_graph<node_cont> *Gret;
  sw_graph<seq_node_cont> *G;
  int seq_idx;
  unsigned int seq_len;
  size_t i;
  sw_node *node, *prev;
  sw_graph<seq_node_cont>::iterator it;

  printf("Building de Bruijn graph ... ");

  // loop through the sequences
  Gret = new sw_graph<node_cont>();
  G = new sw_graph<seq_node_cont>();
  char *seq;
  for(seq_idx = 0; seq_idx < seq_count; seq_idx++) {
    print_progress(seq_idx, seq_count);

    // check that seq_len >= g_subword_length
    seq = seqs[seq_idx];
    seq_len = strlen(seq);
    if(seq_len < g_subword_length) continue;

    // build de bruijn graph for this sequence
    G->nodes.clear();
    node = new sw_node;
    prev = NULL;
    for(i = 0; i+g_subword_length <= seq_len; i++) {
      if(!validate_character(seq[i])) {
        fprintf(stderr, "\nInvalid amino acid: `%c'\n", seq[i]);
        fprintf(stderr, "Sequence: %d & Position: %lu\n", seq_idx+1, i);
        delete node; delete G;
        destroy_graph(Gret);
        return NULL;
      }
      (*background)[seq[i]]++;
      node->set_subword(&seq[i]);
      node->seq = seq_idx;

      // have we seen this subword in this sequence?
      it = G->find(node);
      if(it != G->nodes.end()) {
        (*it)->mul++;

        // add an edge if it isn't there already.
        // O(n), n is the size of the alphabet (small?)
        if(std::find(prev->adj_right.begin(),
                     prev->adj_right.end(),
                     *it) == prev->adj_right.end()) {
          prev->add_edge(*it);
        }
        prev = *it;
        continue;
      }

      // add this node to the graph of the current sequence
      G->add_node(node);
      if(prev) prev->add_edge(node);

      // SIM matrix
      if(SIM) node->SIM_row = SIM->add_row(node->subword);

      // add node to the final graph.
      Gret->add_node(node);

      prev = node;
      node = new sw_node;
    }
    for(i = seq_len - g_subword_length; i < seq_len; i++) {
      if(!validate_character(seq[i])) {
        fprintf(stderr, "\nInvalid amino acid: `%c'\n", seq[i]);
        fprintf(stderr, "Sequence: %d & Position: %lu\n", seq_idx+1, i);
        return NULL;
      }
      (*background)[seq[i]]++;
    }

    delete node;
  }
  print_progress(1, 1, true); // 100 %\n
  delete G;
  return Gret;
}

void destroy_graph(sw_graph<node_cont> *G) {
  node_cont::iterator it;
  for(it = G->nodes.begin(); it != G->nodes.end(); ++it)
    delete *it;
  delete G;
  return;
}

/**
 * old version of add_ic_edges that uses "symmetry" to speed up the computation.
 * does not parallelize nicely.
 bool add_ic_edges(sw_graph<node_cont>* G, SIM_mat *SIM) {
 if(SIM) return add_ic_edges_SIM(G, SIM);

 printf("Calculating node weights and adding IC edges ... ");
 node_cont::iterator it, it2;
 sw_node *node;
 sim_type sim;
 unsigned long i = 0, j = G->nodes.size(); // for print_progress

// O(n^2)
for(it2 = G->nodes.begin(); it2 != G->nodes.end(); ++it2) {
print_progress(1.0*i*(i+1), 1.0*j*(j+1)); i++;
node = *it2;

// take advantage of the symmetry to improve execution time.
for(it = G->nodes.begin(); it != it2; ++it) { // it != it2
sim = node->similarity(*it);
if(sim < g_theta) continue;

node->weight += sim;
(*it)->weight += sim;

if(node->seq != (*it)->seq) {
try_add_ic_edge(node, *it, sim);
try_add_ic_edge(*it, node, sim);

// node->gen_mul++ if we are interested in the number of
// similar nodes instead of sequences
node->try_inc_gen_mul(*it);
(*it)->try_inc_gen_mul(node);
// ^- doesn't parallelize nicely (see try_inc_gen_mul definition
// in sw_graph.hpp)
}
}
}
print_progress(1, 1, true);

printf("Calculating neighbourhood general multplicities ...");
calculate_ngen_muls(G);
printf("\n");
return true;
}
inline bool try_add_ic_edge(sw_node *to, int sim) {
if(sim >= this->max_sim) {
if(sim > this->max_sim) {
this->max_sim = sim;
adj_ic.clear(); // O(n) !
}
adj_ic.push_back(to);
} else return false;
return true;
}
inline bool try_inc_gen_mul(sw_node *to) {
return try_inc_gen_mul(to->seq);
}
inline bool try_inc_gen_mul(int seq) {
if(highest_gm_seq < seq) {
highest_gm_seq = seq;
gen_mul++;
return true;
}
return false;
}
*/

bool add_ic_edges(sw_graph<node_cont>* G, int seq_count, SIM_mat *SIM) {
  if(SIM) return add_ic_edges_SIM(G, seq_count, SIM);

  printf("Calculating node weights and adding IC edges ... ");
  unsigned long i = 0, j = G->nodes.size(); // for print_progress
  for(node_cont::iterator it = G->nodes.begin();
      it != G->nodes.end(); ++it) {
    print_progress(i++, j);
    add_node_ic_edges(G, seq_count, *it);
  }
  print_progress(1, 1, true);

  printf("Calculating neighbourhood general multplicities ...");
  calculate_ngen_muls(G);
  printf("\n");
#ifdef db_cout_both_similarity_finder_processes
  std::cout << "add_ic_edges trad, similar pairs: "
      << g_traditional_pair_counts << std::endl;
#endif
  return true;
}

bool add_node_ic_edges(sw_graph<node_cont> *G, int seq_count, sw_node *node) {
  // keep track of maximal similarity nodes for each graph/sequence.
  // this could be optimized to store only one sim_type and ic_edge_cont.
  // (nodes in G are sorted by sequence indexes in ascending order) 
  std::vector<ic_edge_cont> adj_ics(seq_count);
  std::vector<sim_type> max_sims(seq_count, SIM_TYPE_MIN);
  int highest_gm_seq = -1;
  int similar_count = 0;

  sim_type sim;
  node_cont::iterator it;
  for(it = G->nodes.begin(); it != G->nodes.end(); ++it) {
    sim = node->similarity(*it);
    //if(sim < g_theta) continue; //This makes no sense at all? -kk
    if(sim == SIM_TYPE_MIN) continue; //This ought to be correct
    if(node->seq != (*it)->seq) { // Not in the same sequence
#ifdef db_weight
      std::cout << "nw: " << node->weight << ", sim: " << sim << std::endl;
#endif
      node->weight += sim;
#ifdef db_weight
      std::cout << "new nw:" <<  node->weight << std::endl; 
#endif
      similar_count++;
#ifdef db_cout_similar_nodes
      std::cout << "a: " << node->seq << ":" << node->subword << ", b: " << (*it)->seq << ":" << (*it)->subword << ", s: " << sim << std::endl;
#endif
      // add ic edge if sim >= max_sim
      sim_type &max_sim = max_sims[(*it)->seq];
      ic_edge_cont &adj_ic = adj_ics[(*it)->seq];
      //sim_type &max_sim = max_sims[0]; // old version
      //ic_edge_cont &adj_ic = adj_ics[0];
      if(sim >= max_sim) {
        if(sim > max_sim) {
          max_sim = sim;
          adj_ic.clear(); // O(n)!

          //adj_ic.push_back(*it);  // ONLY THE TOP
        }
        adj_ic.push_back(*it);  // ALL WITH SAME VALUE
      }

      // increment general multiplicity
      if(highest_gm_seq < (*it)->seq) {
        highest_gm_seq = (*it)->seq;
        node->gen_mul++;
      }
    }
  }

  // add inter-component edges to the node.
  // this can be done very fast if ic_edge_cont == std::list
  // is std::copy slow?

  node->adj_ic.clear();
  for(std::vector<ic_edge_cont>::iterator it = adj_ics.begin();
      it != adj_ics.end(); ++it) {
    std::copy((*it).begin(), (*it).end(), std::back_inserter(node->adj_ic));

  }
  //std::cout << similar_count << std::endl;
#ifdef db_cout_both_similarity_finder_processes
  g_traditional_pair_counts += similar_count;
#endif
  return true;
}


// it on yks, node on toinen... olkoon node nyt meidan node1. it node2.
/*void do_work_for_node_pair(sw_node &node1, sw_node &node2, 
  sim_type sim_value) {
  node1.weight += sim_value;
  node2.weight += sim_value;
  bool add_ic_edge_for1 = false;
  bool add_ic_edge_for2 = false;
  bool remove_old_from1 = false;
  bool remove_old_from2 = false;
  bool first_for_this_seq1 = false;
  bool first_for_this_seq2 = false;

  if(node1.seq != node2.seq) { // Not in the same sequence

// Update max_sims for node1 if needed
if (node1.max_sims == NULL) {
node1.max_sims = new map<int, sim_type>;
(*node1.max_sims)[node2.seq] = sim_value;
add_ic_edge_for1 = true;
first_for_this_seq1 = true;
} else {
if ( node1.max_sims->find(node2.seq) == node1.max_sims->end() ) {
(*node1.max_sims)[node2.seq] = sim_value;
first_for_this_seq1 = true;
} else if ( (*node1.max_sims)[node2.seq] < sim_value) {
(*node1.max_sims)[node2.seq] = sim_value;
add_ic_edge_for1 = true;
remove_old_from1 = true;
}
}

// Update max_sims for node2 if needed
if (node2.max_sims == NULL) {
node2.max_sims = new map<int, sim_type>;
(*node2.max_sims)[node1.seq] = sim_value;
add_ic_edge_for2 = true;
first_for_this_seq2 = true;
} else {
if ( node2.max_sims->find(node1.seq) == node2.max_sims->end() ) {
(*node2.max_sims)[node1.seq] = sim_value;
first_for_this_seq2 = true;
} else if ( (*node2.max_sims)[node1.seq] < sim_value) {
(*node2.max_sims)[node1.seq] = sim_value;
add_ic_edge_for2 = true;
remove_old_from2 = true;
}
}

// Update tmp_ic_edges for node1 if needed
if (add_ic_edge_for1) {
if (node1.tmp_ic_edges == NULL) {
node1.tmp_ic_edges = new map<int, ic_edge_cont>;
(*node1.tmp_ic_edges)[node2.seq].push_back(&node2);
} else {
if (remove_old_from1) {  // ok, not really necessary
(*node1.tmp_ic_edges)[node2.seq].clear();
}
(*node1.tmp_ic_edges)[node2.seq].push_back(&node2);
}
}

// Update tmp_ic_edges for node2 if needed
if (add_ic_edge_for2) {
if (node2.tmp_ic_edges == NULL) {
node2.tmp_ic_edges = new map<int, ic_edge_cont>;
(*node2.tmp_ic_edges)[node1.seq].push_back(&node1);
} else {
if (remove_old_from2) {  // ok, not really necessary
(*node2.tmp_ic_edges)[node1.seq].clear();
}
(*node2.tmp_ic_edges)[node1.seq].push_back(&node1);
}
}

// increment general multiplicity if first edges assigned for these seqs
if (first_for_this_seq1) {
  node1.gen_mul++;
}

if (first_for_this_seq2) {
  node2.gen_mul++;
}
}
}*/
// TODO;
// * poistettava nuo mapit ja pointterit lopuksi niistä joilla ne on
// * heitettävä adj_ic:hen nuo mapissa olevat ja poistettava neki mapit
// * tupla map-accesseja voi siivota

bool add_ic_edges_treelike(sw_graph<node_cont> *G, int seq_count) {
  // use SimilarityFinder to do the right things.
  // change two last parameters to optimize
  SimilarityFinder sf(g_max_distance, g_max_distance*4.0, 2.0);
  std::cout << "ProcessNew: " << sf.ProcessNew(G->nodes) << std::endl;
  std::cout << "Finished treelike" << std::endl;
  return true;
}


// builds the subword graph using SIM matrix
bool add_ic_edges_SIM(sw_graph<node_cont> *G, int seq_count, SIM_mat *SIM) {
  unsigned long sw_count = SIM->get_sw_count();
  printf("Distinct subwords: %lu\n", sw_count);

  // calculate size for the SIM matrix
  printf("The SIM matrix size will be: ");
  float SIM_sz;
  SIM_sz = sw_count * sizeof(sim_type*);
  SIM_sz += (((sw_count-1)*sw_count) / 2.0) * sizeof(sim_type);

  // print in KB, MB or GB
  SIM_sz /= 1024.0; // in KB
  if(SIM_sz < 1024.0)
    printf("%.2f KB\n", SIM_sz);
  else if(SIM_sz < 1024.0*1024.0)
    printf("%.2f MB\n", SIM_sz/1024.0);
  else printf("%.2f GB\n", SIM_sz/(1024.0*1024.0));

  // allocate and fill the SIM matrix
  printf("Computing the SIM matrix ... ");
  SIM->allocate();
  SIM->compute(true);

  // add inter-component edges and calculate weights.
  // duplicate code.. this should be merged with add_node_ic_edges
  printf("Calculating node weights and adding IC edges ... ");
  std::vector<ic_edge_cont> adj_ics(seq_count);
  std::vector<sim_type> max_sims(seq_count);
  int highest_gm_seq;
  node_cont::iterator it, it2;
  sw_node *node;
  sim_type sim;
  unsigned long i = 0, j = G->nodes.size(); // for print_progress
  for(it2 = G->nodes.begin(); it2 != G->nodes.end(); ++it2) {
    print_progress(i++, j);
    node = *it2;

    highest_gm_seq = -1;
    std::fill(max_sims.begin(), max_sims.end(), SIM_TYPE_MIN);
    for(it = G->nodes.begin(); it != G->nodes.end(); ++it) {
      sim = SIM->get_sim(node->SIM_row, (*it)->SIM_row);
      //if(sim < g_theta) continue; 
      if(sim == SIM_TYPE_MIN) continue; //kk: fixed to new sim

      node->weight += sim;
      if(node->seq != (*it)->seq) {
        // add ic edge if sim >= max_sim
        sim_type &max_sim = max_sims[(*it)->seq];
        ic_edge_cont &adj_ic = adj_ics[(*it)->seq];
        //sim_type &max_sim = max_sims[0]; // old..
        //ic_edge_cont &adj_ic = adj_ics[0];
        if(sim >= max_sim) {
          if(sim > max_sim) {
            max_sim = sim;
            adj_ic.clear(); // O(n)!
          }
          adj_ic.push_back(*it);
        }

        if(highest_gm_seq < (*it)->seq) {
          highest_gm_seq = (*it)->seq;
          node->gen_mul++;
        }
      }
    }

    // add inter-component edges to the node.
    for(std::vector<ic_edge_cont>::iterator it = adj_ics.begin();
        it != adj_ics.end(); ++it) {
      std::copy((*it).begin(), (*it).end(), std::back_inserter(node->adj_ic));
      (*it).clear();
    }
  }
  print_progress(1, 1, true);

  printf("Calculating neighbourhood general multiplicities ...");
  calculate_ngen_muls(G);
  printf("\n");

  //convert formed map into a nice vector.


  return G;
}

void calculate_ngen_muls(sw_graph<node_cont> *G, bool print) {
  node_cont::iterator it2;
  edge_cont::iterator it, it3;
  std::set<sw_node*> seen, seen2;
  sw_node *node, *node2;
  float sum, sum2;
  int i, j; // print_progress

  if(print) {
    i = 0;
    j = G->nodes.size();
  } else j = 0; // silence warning
  for(it2 = G->nodes.begin(); it2 != G->nodes.end(); ++it2) {
    node = *it2;
    seen.clear();
    seen2.clear();
    sum = sum2 = 0;
    if(print) print_progress(i++, j);

    for(it = node->adj_left.begin(); it != node->adj_left.end(); ++it) {
      if(seen.insert(*it).second) {
        sum += (*it)->gen_mul;
        node2 = *it;
        for(it3 = node2->adj_left.begin(); it3 != node2->adj_left.end(); ++it3) {
          if(seen2.insert(*it3).second) {
            sum2 += (*it3)->gen_mul;
          }
        }

      }
    }
    for(it = node->adj_right.begin(); it != node->adj_right.end(); ++it) {
      if(seen.insert(*it).second) {
        sum += (*it)->gen_mul; 
        node2 = *it;
        for(it3 = node2->adj_right.begin(); it3 != node2->adj_right.end(); ++it3) {
          if(seen2.insert(*it3).second) {  
            sum2 += (*it3)->gen_mul;
          }
        }

      }
    }

    node->ngen_mul = 4.0/7.0 * node->gen_mul + 2.0/7.0 * (sum / seen.size()) + 1.0/7.0 * (sum2 / seen2.size());
  }
  if(print) print_progress(1, 1, true);
  return;
}
