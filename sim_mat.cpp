#include <iostream>
#include <algorithm>
#include "similarity.hpp"
#include "sim_mat.hpp"
#include "util.hpp"

void SIM_mat::allocate() {
    int dim = sw_map.size();
    if(mat) this->free();
    mat = new sim_type*[dim];
    for(int i = 0; i < dim; i++) {
        mat[i] = new sim_type[i+1];
    }
    return;
}

void SIM_mat::free() {
    if(mat) {
        int i, j;
        for(i = 0, j = sw_map.size(); i < j; i++)
            delete [] mat[i];
        delete [] mat;
        mat = NULL;
    }
    return;
}

// computes the SIM matrix O(n^2).
bool SIM_mat::compute(bool print) {
    SIM_mat::iterator it, it2;
    unsigned int i, j, c, sw_count;
    if(!mat) allocate();

    sw_count = get_sw_count();
    for(c = 0, it = sw_map.begin(); it != sw_map.end(); c++, ++it) {
        if(print) print_progress(c*(c+1), sw_count*(sw_count+1));

        it2 = sw_map.begin();
        while(1) {
            i = it->second; j = it2->second;
            if(i < j) {
                mat[j][i] = similarityfn(it->first, it2->first);
            } else {
                mat[i][j] = similarityfn(it->first, it2->first);
            }

            if(it == it2) break;
            ++it2;
        }
    }
    print_progress(1, 1, true);

    return true;
}
