#ifndef TRAV_HPP_INCLUDED
#define TRAV_HPP_INCLUDED

#include <deque>
#include <list>
#include "sw_graph.hpp"

typedef std::pair<sw_node*, sw_node*> node_pair;

// list is faster than vector here (O(1) delete)
struct seq_path;
typedef std::list<seq_path> path_cont;

// counters are defined as <char, unsigned int> pairs.
// TODO: create counter-class
typedef std::map<char, unsigned int> counter_map;
typedef std::pair<char, unsigned int> counter;
typedef std::vector<counter> counter_cont;

enum direction { LEFT = 0, RIGHT };

// represents a path in one sequence
struct seq_path {
    union {
        sw_node *node[2];
        struct {
            sw_node *left; // node[LEFT]
            sw_node *right; // node[RIGHT]
        };
    };
    unsigned int errors;

    seq_path(sw_node *node) {
        left = right = node;
        errors = 0;
    }
    inline unsigned int get_seq_idx() { return left->seq; } // or right->seq

    // moves left/right nodes according to the counters.
    // the counters have to be sorted in desceding order.
    // returns false if the move fails or is impossible.
    bool move(counter_cont *counters, direction dir);
};

// position weight matrix vector
struct PWM_vec {
    // again, std::map is a little heavy structure for our small alphabet.
    std::map<char, unsigned int> freqs;
    //std::vector<unsigned int> freqs;

    PWM_vec() {}
    ~PWM_vec() {}

    unsigned int& operator[] (const char key) { return freqs[key]; }
    char get_char();

    //draw_logo etc
};

// find_consensus() returns an instance of this class
struct consensus_info {
    path_cont paths;

    // position weight matix
    std::deque<PWM_vec> PWM;
    std::string consensus; // initialize_consensus_word()
    std::string consensus_top_aas;

    // max node's general multiplicity and weight before adjustments
    sw_node *max_node;
    sim_type max_node_w;
    unsigned int max_node_gm;
    float max_node_ngm;

    // call initialize_score() before using these variables
    long double score;
    long double credability, credability_normalized;
    long double log_odds, log_odds_normalized;
    long double log_p_background_; //the following 3 are just for printing
    long double log_support_;
    long double log_number_of_seqs_;
    unsigned int err_total;
    unsigned int path_count;
    long double weighted_log_odds;
    long double multiplied;
    long double gm, weight, gm_normalized, weight_normalized;
    unsigned int rank_weighted_log_odds, rank_score, rank_log_odds, 
                 rank_credability, rank_multiplied, rank_gm, rank_weight,
                 rank_gm_then_weight;

    consensus_info() {};
    ~consensus_info() {};

    // initializes current consensus word from PWM
    void initialize_consensus_word();

    // calculates score and related info
    void initialize_score(map<char,int> &background, unsigned int seq_count);

    void print();
    void printPWM();

    inline long int add_score(long int add) {
        this->score += add;
        return this->score;
    }
    inline long int add_score(consensus_info *cinfo) {
        this->score += cinfo->score;
        this->path_count += cinfo->path_count;
        this->err_total += cinfo->err_total;
        return this->score;
    }
};

// sorts consensus_infos by score.
// call consensus_info::initialize_score before using this class
struct pcinfo_score_cmp_desc {
    bool operator()(const consensus_info *a, const consensus_info *b) const {
        return a->score > b->score;
    }
};

struct pcinfo_normalized_weight_and_gm_cmp_desc {
    bool operator()(const consensus_info *a, const consensus_info *b) const 
    {
        return ((a->gm_normalized + a->weight_normalized) > 
            (b->gm_normalized + b->weight_normalized));
    }
};

struct pcinfo_normalized_weight_gm_log_odds_credability_cmp_desc {
    bool operator()(const consensus_info *a, const consensus_info *b) const 
    {
        return ((a->gm_normalized + a->weight_normalized + 
              a->log_odds_normalized + a->credability_normalized) > 
            (b->gm_normalized + b->weight_normalized + 
             b->log_odds_normalized + b->credability_normalized));
    }
};

struct pcinfo_normalized_weight_gm_weighted_log_odds_cmp_desc {
    bool operator()(const consensus_info *a, const consensus_info *b) const 
    {
        return ((a->gm_normalized + a->weight_normalized + 
              a->weighted_log_odds * 2.0) > 
            (b->gm_normalized + b->weight_normalized + 
              b->weighted_log_odds * 2.0));
    }
};


struct pcinfo_weighted_log_odds_cmp_desc {
    bool operator()(const consensus_info *a, const consensus_info *b) const {
        return a->weighted_log_odds > b->weighted_log_odds;
    }
};

struct pcinfo_log_odds_cmp_desc {
    bool operator()(const consensus_info *a, const consensus_info *b) const {
        return a->log_odds > b->log_odds;
    }
};

struct pcinfo_credability_cmp_desc {
    bool operator()(const consensus_info *a, const consensus_info *b) const {
        return a->credability > b->credability;
    }
};

struct pcinfo_multiplied_cmp_desc {
    bool operator()(const consensus_info *a, const consensus_info *b) const {
        return a->multiplied > b->multiplied;
    }
};

struct pcinfo_weight_cmp_desc {
    bool operator()(const consensus_info *a, const consensus_info *b) const {
        return a->weight > b->weight;
    }
};

struct pcinfo_gm_cmp_desc {
    bool operator()(const consensus_info *a, const consensus_info *b) const {
        return a->gm > b->gm;
    }
};

struct pcinfo_gm_then_weight_desc {
    bool operator()(const consensus_info *a, const consensus_info *b) const {
        if(a->gm == b->gm)
        {
          return a->weight > b->weight;
        }
        else
        {
          return a->gm > b->gm;
        }
    }
};


bool is_a_gap(counter_cont &counts_of_position);


void AssignRanks(std::list<consensus_info*> &cinfos);

// initializes counters and returns them in descending order.
// returns the number of counters
int init_counters(counter_cont *C, path_cont *paths, direction dir);

struct counter_cmp_desc {
    bool operator() (const counter &lhs, const counter &rhs) const {
        return lhs.second > rhs.second;
    }
};

// repeatedly moves the left or right nodes of given paths according to dir.
// the function may erase items from paths if a path gets stuck or the error
// threshold is met (if g_error_threshold is non-zero). the function stops
// when the highest counter value is less than
// g_stop_threshold * cinfo.paths.size()
//
// returns the number of successful moves.
int move_paths(consensus_info &cinfo, bool adjust=false, SIM_mat *SIM=NULL);

// adjusts node weights
void adjust_weights(path_cont &paths, direction dir, SIM_mat *SIM=NULL);

// returns consensus word and related info (see consensus_info)
consensus_info* find_consensus(sw_graph<node_cont> *G, unsigned int seq_count,
        bool adjust=false, sw_node *start=NULL, SIM_mat *SIM=NULL);

// print occurrences
void print_occ(char *seqs[], path_cont &paths, const char *consensus);

void FilterOut(std::list<consensus_info*> &cinfos);


#endif
