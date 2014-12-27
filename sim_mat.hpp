#ifndef SW_GRAPH_SIM_HPP_INCLUDED
#define SW_GRAPH_SIM_HPP_INCLUDED

#include <string.h>
#include "sw_graph.hpp"

// matrix SIM.
// subwords are added with method add_row, and after that the actual matrix
// can be allocated using allocate() and filled with compute().
struct SIM_mat {
    // sim_type limits the size of g_subword_length
    //typedef unsigned char sim_type; // defined in similarity.hpp
    typedef std::map<char*, int, sw_cmp>::iterator iterator;

private:
    sim_type **mat; // symmetric. the other half is not stored

    // maps a subword to index/row in mat
    std::map<char*, int, sw_cmp> sw_map;

public:
    SIM_mat() { mat = NULL; }
    ~SIM_mat() { free(); }

    // returns row/index of the inserted subword.
    // if the matrix already had the subword, then the existing row is returned
    inline int add_row(char *sw) {
        int i = sw_map.size();
        std::pair<iterator,bool> it;
        it = sw_map.insert(std::pair<char*, int>(sw, sw_map.size()));
        return (it.second) ? i : (*it.first).second;
    }

    // returns row of the subword in the matrix.
    // return -1 if the subword is not in the matrix
    inline int get_row(char *sw) {
        iterator it = sw_map.find(sw);
        return (it != sw_map.end()) ? it->second : -1;
    }

    // allocates the matrix.
    // add subwords with add_row() before allocating.
    void allocate();
    void free();

    // fills the matrix after it has been allocated
    bool compute(bool print=false);

    inline unsigned int get_sw_count() { return sw_map.size(); } // O(1)
    inline sim_type** get_matrix() { return mat; }
    
    inline sim_type get_sim(unsigned int a, unsigned int b) {
        return (a > b) ? mat[a][b] : mat[b][a];
    }
    inline sim_type get_sim(char *a, char *b) {
        return get_sim(sw_map.find(a)->second, sw_map.find(b)->second);
    }
};

#endif
