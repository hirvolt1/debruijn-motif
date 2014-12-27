#include <math.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "util.hpp"

void print_progress(float cur, float max, bool newline) {
    static int last_progress = -1;
    static bool last_nl;
    static int last_chars = 0;
    int progress, i;
    char bsbuf[16];

    progress = round(100.0 * cur/max);
    if(progress != last_progress || newline != last_nl) {

        // printf '\b' * last_chars (remove old text)
        for(i = 0; i < last_chars; bsbuf[i++] = '\b');
        bsbuf[i] = '\0';
        printf("%s", bsbuf);

        if(!newline) {
            last_chars = printf("%d %%", progress);
            fflush(stdout);
        } else {
            printf("%d %%\n", progress);
            last_chars = 0;
        }

        last_progress = progress;
        last_nl = newline;
    }
    return;
}

void print_stats(char **seqs, unsigned int seq_count, sw_graph<> *G) {
    unsigned int i, j, len;
    unsigned long total_len, distinct, nodes, edges, ic_edges;
    unsigned long sw_count;
    std::set<char*, sw_cmp> sws;

    // calculate the number of subwords
    total_len = distinct =  0;
    for(i = 0; i < seq_count; i++) {
        len = strlen(seqs[i]);
        for(j = 0; j+g_subword_length <= len; j++) {
            if(sws.insert(&seqs[i][j]).second) {
                distinct++;
            }
        }
        total_len += len;
    }
    sws.clear();
    if(total_len < seq_count*(g_subword_length-1)) {
        sw_count = 0;
    } else sw_count = total_len - seq_count*(g_subword_length-1);

    // count nodes & edges
    node_cont::iterator it;
    std::vector<edge_cont>::iterator it2;
    sw_node *node;
    nodes = edges = ic_edges = 0;
    for(it = G->nodes.begin(); it != G->nodes.end(); it++) {
        node = *it;

        nodes++;
        edges += node->adj_right.size();
        ic_edges += node->adj_ic.size();
    }

    printf("Subwords:\t\t%lu\n"
           "Distinct subwords:\t%lu\n"
           "Nodes:\t\t\t%lu\n"
           "Edges:\t\t\t%lu\n"
           "Inter-component edges:\t%lu\n"
           "Average out-degree:\t%.2f\n",
           sw_count, distinct, nodes,
           edges+ic_edges, ic_edges,
           (1.0f*edges)/nodes);
    return;
}

// trivial algorithm..
std::vector<size_t>* search_mismatch(const char *text, const char *pat,
        unsigned int mismatches) {
    std::vector<size_t> *ret = new std::vector<size_t>;
    size_t m, n, i, j;

    // calculate lengths
    m = strlen(pat);
    n = strlen(text);
    
    // search
    size_t mism;
    for(j = 0; j < n-m+1; j++) {
        mism = 0;
        for(i = 0; i < m; i++) {
            if(pat[i] != text[i+j]) {
                mism++;
                if(mism > mismatches)
                    break;
            }
        }

        if(mism <= mismatches) {
            ret->push_back(j);
        }
    }

    return ret;
}

std::vector<size_t>* search_errors(const char *text, const char *pat,
        unsigned int errors) {
    std::vector<size_t> *ret = new std::vector<size_t>;
    size_t m, n, i, j;
    unsigned int *Dcol, *delta;

    // calculate lengths
    m = strlen(pat);
    n = strlen(text);

    // allocate Dcol and delta arrays.
    // delta keeps track of the number of insertions and deletions and is used
    // to determine starting positions of the occurences
    Dcol = new unsigned int[m+1]; // indexing starts from 1
    delta = new unsigned int[m+1];
    for(i = 0; i < m+1; i++) {
        Dcol[i] = delta[i] = i;
    }

    // search
    unsigned int top, prev, prev_delta, tmp, tmp_delta, diag;
    top = (errors+1 < m) ? errors+1 : m;
    Dcol[0] = delta[0] = 0;
    for(j = 0; j < n; j++) {

        // prev = D[i-1][j-1]
        prev = prev_delta = 0;
        for(i = 1; i <= top; i++) {
            tmp = Dcol[i];
            tmp_delta = delta[i];
            
            // D[i][j] = min(D[i-1][j-1] + (if P[i] = S[j] then 1 else 0),
            //               D[i-1][j] + 1, D[i][j-1] + 1)
            diag = (pat[i-1] == text[j]) ? prev : prev + 1;
            if(diag < Dcol[i]+1 && diag < Dcol[i-1]+1) {
                Dcol[i] = diag;
                delta[i] = prev_delta;
            } else {
                if(Dcol[i] < Dcol[i-1]) {
                    Dcol[i] = Dcol[i] + 1;
                    delta[i] = delta[i] - 1;
                } else {
                    Dcol[i] = Dcol[i-1] + 1;
                    delta[i] = delta[i-1] + 1;
                }
            }

            prev = tmp;
            prev_delta = tmp_delta;
        }

        while(Dcol[top] > errors) top--;

        if(top == m) {
            // match
            ret->push_back(j-m + delta[top]+1);
        } else top++;
    }

    delete [] delta;
    delete [] Dcol;
    return ret;
}

unsigned int best_ham_dist(const char *text, const char *pat) {
    unsigned int i, j,  m, n, tmp, best, ham;
    const char *str;

    m = strlen(pat); n = strlen(text);
    if(strlen(text) < strlen(pat)) {
        // swap text & pat
        tmp = m; m = n; n = tmp;
        str = text; text = pat; pat = str;
    }

    best = UINT_MAX;
    for(i = 0; i <= n-m && best; i++) {
        ham = m;
        for(j = 0; j < m; j++) {
            if(text[i+j] == pat[j]) {
                ham--;
            }
        }

        if(ham < best) best = ham;
    }
    return best;
}


