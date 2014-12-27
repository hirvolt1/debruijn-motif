// declares global variables and commonly used containers
#ifndef COMMON_HPP_INCLUDED
#define COMMON_HPP_INCLUDED

#include <list>
#include <set>
#include <map>
#include <unordered_set>
#include <vector>
#include <string>
#include <string.h>

// global variables, defined in main.cpp
extern unsigned int g_subword_length;
extern float g_theta;
extern float g_tau;
extern float g_stop_threshold;
extern unsigned int g_gaps;
extern unsigned int g_motif_length;
extern unsigned int g_error_threshold;
extern int g_verb;
extern bool g_SS_tree_flag;

// debug printf
#define printfv1(M, ...) if(g_verb >= 1) { printf(M, ##__VA_ARGS__); }
#define printfv2(M, ...) if(g_verb >= 2) { printf(M, ##__VA_ARGS__); }

// forward declarations
struct sw_node;
struct SIM_mat;
struct sw_node_cmp;
struct sw_node_eq;
struct sw_node_hash;

// binary search tree
typedef std::multiset<sw_node*, sw_node_cmp> set_cont;

// dynamic array
typedef std::vector<sw_node*> vec_cont;

// linked list
typedef std::list<sw_node*> list_cont;

// hash table (requires -std=c++0x or -std=c++11)
typedef std::unordered_multiset<sw_node*, sw_node_hash, sw_node_eq> hash_cont;
#define HASHTABLE_BUCKETS 10000

// compare class for subwords
struct sw_cmp {
    bool operator() (const char *lhs, const char *rhs) const {
        return strncmp(lhs, rhs, g_subword_length) < 0;
    }
};
struct sw_eq {
    bool operator ()(const char *lhs, const char *rhs) const {
        return strncmp(lhs, rhs, g_subword_length) == 0;
    }
};
struct sw_hash { // TODO: a better hash function?
    std::size_t operator ()(const char *subword) const {
        return std::unordered_set<std::string>::hasher()
            (std::string(subword, g_subword_length));
    }
};

#endif
