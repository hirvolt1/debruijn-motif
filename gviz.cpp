#ifdef GVIZ
#include <stdio.h>
#include <sstream>
#include <gvc.h>
#include "gviz.hpp"

// the code is horrible :)

// disable (const char*) -> (char*) conversion warning
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wwrite-strings"

// the names must be unique !
std::string get_node_name(sw_node *node) {
    std::string ret;
    std::stringstream stm;
    stm.str("");
    ret.assign(node->subword, g_subword_length);
    stm << ret << "\n"
        << "(" << node->seq << "," << node->weight << ","
        << node->mul << "," << node->gen_mul << ")";
    ret = stm.str();
    return ret;
}

void render_graph(char **seqs, int seq_count, sw_graph<> *G, char *out) {
    GVC_t *gvc;
    std::string oout, tmp;
    Agraph_t *g, **clusters, *c;
    int i;
    std::stringstream sstm;
    std::map<sw_node*, Agnode_t*> agnodes; // maps sw_nodes to ag_nodes
    Agnode_t *n;

    // setup context and parameters
    gvc = gvContext();
    oout = "-o"; oout += out; oout+= ".png";
    char *args[] = {
        "dot",
        "-Tpng",
        //"-v",
        (char*)oout.c_str()
    };
    gvParseArgs(gvc, sizeof(args)/sizeof(char*), args);

    // create directed graph
    g = agopen("g", AGDIGRAPH); // ADIGRAPHSTRICT?
    agsafeset(g, "compound", "true", "");
    agsafeset(g, "labelloc", "top", "");

    sstm.str("");
    sstm << "k=" << g_subword_length << ", t=" << g_theta << "\n"
        << "(seq, weight, multiplicity, generalized multiplicity)";
    agsafeset(g, "label", (char*)sstm.str().c_str(), "");

    // subgraphs (clusters)
    clusters = new Agraph_t*[seq_count];
    for(i = 0; i < seq_count; i++) {
        sstm.str("");
        sstm << "cluster_" << i;
        clusters[i] = agsubg(g, (char*)sstm.str().c_str());
        agsafeset(clusters[i], "label", seqs[i], "");
        //agsafeset(clusters[i], "pencolor", "white", "");
    }

    // first pass: add nodes
    node_cont::iterator it;
    edge_cont::iterator it2, it3;
    sw_node *node;
    for(it = G->nodes.begin(); it != G->nodes.end(); it++) {
        node = *it;

        // create agnode
        tmp = get_node_name(node);
        agnodes[node] = agnode(clusters[node->seq], (char*)tmp.c_str());
    }
    
    // second pass: add edges
    for(it = G->nodes.begin(); it != G->nodes.end(); it++) {
        node = *it;
        c = clusters[node->seq];
        n = agnodes[node];

        for(it2 = node->adj_right.begin(); it2 != node->adj_right.end(); it2++)
            agedge(c, n, agnodes[*it2]);
        for(it2 = node->adj_ic.begin(); it2 != node->adj_ic.end(); it2++)
            agsafeset(agedge(g, n, agnodes[*it2]), "color", "red", "");
    }

    // render
    gvLayoutJobs(gvc, g);
    gvRenderJobs(gvc, g);

    // free resources
    gvFreeLayout(gvc, g);
    delete [] clusters;
    agclose(g);
    gvFreeContext(gvc);
    return;
}
#pragma GCC diagnostic pop
#endif
