#ifndef UTIL_HPP_INCLUDED
#define UTIL_HPP_INCLUDED

#include "sw_graph.hpp"

void print_progress(float cur, float max, bool newline=false);
void print_stats(char **seqs, unsigned int seq_count, sw_graph<> *G);

// returns approximate occurrences of pat in text with at most k mismatches.
// trivial algorithm (SLOW)
std::vector<size_t>* search_mismatch(const char *text, const char *pat,
        unsigned int mismatches);

// returns approximate occurrences of pat in text with at most k errors.
// not currently used for anything.
std::vector<size_t>* search_errors(const char *text, const char *pat,
        unsigned int errors);

// finds the best hamming distance between pat and all |pat| length substrings
// of text.
unsigned int best_ham_dist(const char *text, const char *pat);

template <class T>
void UpdateMax(T &compared, T &maxval)
{
  if(compared > maxval)
  {
    maxval = compared;
  }
}

template <class T>
void UpdateMin(T &compared, T &minval)
{
  if(compared < minval)
  {
    minval = compared;
  }
}

#endif
