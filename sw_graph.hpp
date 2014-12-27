#ifndef SW_GRAPH_HPP_INCLUDED
#define SW_GRAPH_HPP_INCLUDED

#include <string>
#include <string.h>
#include <map>

#include "common.hpp"
#include "similarity.hpp"

using std::map;

// default (stable!) container for nodes.
typedef vec_cont node_cont;

// used during building the graph to check whether a node has been already
// encountered. should support fast inserts and lookups
typedef set_cont seq_node_cont;
//typedef hash_cont seq_node_cont;
#define HASHTABLE_BUCKETS 10000

// edge container
typedef vec_cont edge_cont;
typedef vec_cont ic_edge_cont;

// a subword node.
struct sw_node {
    // edges
    edge_cont adj_right;        // edges originating from this node
    edge_cont adj_left;         // list of nodes pointing to this node
    ic_edge_cont adj_ic;        // adjacent inter-component nodes
    map<int, sim_type> max_sims; // Null pointer = no value assigned
    map<int, ic_edge_cont> tmp_ic_edges; // Null pointer = no value assigned

    char *subword;              // pointer to the k-mer
    int seq;                    // seq idx
    unsigned int mul;           // multiplicity
    sim_type weight;            // sum of similarities >= theta
    unsigned int gen_mul;       // general multiplicity
    float ngen_mul;             // neighbourhood general multiplicity
    unsigned int SIM_row;       // row/index in the SIM matrix
    unsigned int visits;        // m'(v_i). counts visits during graph traversal

    // constructor
    sw_node() {
        mul = gen_mul = 1;
    }
    sw_node(char *sw) { sw_node(); set_subword(sw); }

    // set subword
    inline void set_subword(char *sw) {
        this->subword = sw;
        this->weight = similarityfn(sw, sw);
    }

    // returns similarity
    inline sim_type similarity(sw_node *n) {
        return similarityfn(this->subword, n->subword);
    }

    // method for adding edges
    inline bool add_edge(sw_node *to) {
        this->adj_right.push_back(to);
        to->adj_left.push_back(this);
        return true;
    }
};

// sw_node* compare classes
struct sw_node_cmp {
    bool operator() (const sw_node *lhs, const sw_node *rhs) const {
        return strncmp(lhs->subword, rhs->subword, g_subword_length) < 0;
    }
};
struct sw_node_eq {
    bool operator ()(const sw_node *lhs, const sw_node *rhs) const {
        return strncmp(lhs->subword, rhs->subword, g_subword_length) == 0;
    }
};
struct sw_node_hash { // TODO: a better hash function?
    std::size_t operator ()(const sw_node *node) const {
        return std::unordered_set<std::string>::hasher()
            (std::string(node->subword, g_subword_length));
    }
};

// subword graph.
// templates make changing the type of the node container "easy"..
template<class Container = node_cont>
struct sw_graph {
    Container nodes;

    inline sw_node* add_node(sw_node *n) {
        nodes.insert(nodes.end(), n); // remove the first param?
        return n;
    }

    inline typename Container::iterator find(sw_node *n) {
        return this->nodes.find(n);
    }

    inline typename Container::iterator equal_range(sw_node *n) {
        return this->nodes.equal_range(n);
    }
};

// template specializations
template<>
struct sw_graph<set_cont> {
    typedef set_cont::iterator iterator;
    set_cont nodes;
    inline sw_node* add_node(sw_node *n) {
        nodes.insert(n);
        return n;
    }
    inline iterator find(sw_node *n) { return nodes.find(n); }
};

template<>
struct sw_graph<hash_cont> {
    typedef hash_cont::iterator iterator;
    hash_cont nodes;
    sw_graph(int buckets=HASHTABLE_BUCKETS) {
        nodes.rehash(buckets); // set the number of buckets. STL's default: 10
    }
    inline sw_node* add_node(sw_node *n) {
        nodes.insert(n);
        return n;
    }
    inline iterator find(sw_node *n) { return nodes.find(n); }
};

template<>
struct sw_graph<vec_cont> {
    typedef vec_cont::iterator iterator;
    vec_cont nodes;
    inline sw_node* add_node(sw_node *n) {
        nodes.insert(nodes.end(), n);
        return n;
    }
    inline vec_cont::iterator find(sw_node *n) {
        vec_cont::iterator it;
        for(it = this->nodes.begin(); it != nodes.end(); it++) {
            if(!strncmp((*it)->subword, n->subword, g_subword_length))
                break;
        }
        return it;
    }
};

template<>
struct sw_graph<list_cont> {
    list_cont nodes;
    inline sw_node* add_node(sw_node *n) {
        nodes.insert(nodes.end(), n);
        return n;
    }
    inline list_cont::iterator find(sw_node *n) {
        list_cont::iterator it;
        for(it = this->nodes.begin(); it != nodes.end(); it++) {
            if(!strncmp((*it)->subword, n->subword, g_subword_length))
                break;
        }
        return it;
    }
};

// functions

// builds the de Bruijn graphs. does not calculate node weights, general
// multiplicities or add inter-component edges.
//
// the returned nodes are in ascending order in terms of seq indexes.
// other code relies on this behavior.
// Calculates background distribution / model.
sw_graph<node_cont>* build_debruijn_graph(char **seqs, int seq_count,
        SIM_mat *SIM=NULL, map<char,int> *background=NULL);

// destroys a graph returned by build_debruijn_graph
void destroy_graph(sw_graph<node_cont> *G);

// creates inter-component edges and calculates node weights and general
// multiplicities. G is a graph returned by build_debruijn_graphs
bool add_ic_edges(sw_graph<node_cont> *G, int seq_count, SIM_mat *SIM=NULL);

// adds inter-component edges to one node. also calculates gen_muls, weight etc
bool add_node_ic_edges(sw_graph<node_cont> *G, int seq_count, sw_node *node);

// same as add_ic_edges but uses SIM matrix
bool add_ic_edges_SIM(sw_graph<node_cont> *G, int seq_count, SIM_mat *SIM);

// same as add_ic_edges, but uses a similarity tree to speed up process
bool add_ic_edges_treelike(sw_graph<node_cont> *G, int seq_count);//, 
                           //SIM_mat *SIM=NULL);
void do_work_for_node_pair(sw_node &node1, sw_node &node2, 
                           sim_type sim_value);


// calculates neighbourhood multiplicities for all nodes in G
void calculate_ngen_muls(sw_graph<node_cont> *G, bool print=false);

#endif
