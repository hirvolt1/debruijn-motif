#ifndef SIMILARITY_HPP_INCLUDED
#define SIMILARITY_HPP_INCLUDED

#include "limits.h"
#include "common.hpp"
#include "float.h"

// similarity type
typedef long double sim_type;
#define SIM_TYPE_MIN FLT_MIN
#define SIM_TYPE_MAX FLT_MAX

//#define SIM_TYPE_MIN LONG_MIN
// TODO: define separate weight_type...

// ~2.64 KB per matrix
extern int g_scoring_matrix[26][26];
extern bool g_matches_in_scoring_matrix[26][26];
extern sim_type g_distance_matrix[26][26];
extern const int g_blosum50[26][26];
extern const int g_blosum62[26][26];
extern const int g_pam30[26][26];
//extern float g_mismatch_avg;
//extern float g_match_avg; 
extern sim_type g_min_similarity;// = (g_theta * g_match_avg + (1-g_theta) * g_mismatch_avg) * (float)g_subword_length;
extern sim_type g_max_distance;

void set_scoring_matrix(const int matrix[26][26]);

//requires set_scoring_matrix first -> derived from that.
void set_distance_matrix();

//checks that triangle inequality holds for the created g_distance_matrix
bool check_distance_matrix_triangle_inequality();

// similarity function for subwords.
inline sim_type similarityfn(const char *a, const char *b) {
  sim_type score;
  sim_type sim = 0;
  for(unsigned int i = 0; i < g_subword_length; i++) {
    score = g_scoring_matrix[a[i]-'A'][b[i]-'A'];  // implicit type conversion
    sim += score;
  }
  if(sim < g_min_similarity) {
    return SIM_TYPE_MIN;
  }
  return sim;
}

inline sim_type distancefn(const char *a, const char* b) {
  sim_type score;
  sim_type dist = 0;
  for(unsigned int i = 0; i < g_subword_length; i++) {
    score = g_distance_matrix[a[i]-'A'][b[i]-'A'];
    dist += score;
  }
  if(dist > g_max_distance) {
    return SIM_TYPE_MAX;
  }
  return dist;
}

inline sim_type similarityfn(std::string &a, std::string &b) {
  sim_type score;
  sim_type sim = 0;
  for(unsigned int i = 0; i < g_subword_length; i++) {
    score = g_scoring_matrix[a[i]-'A'][b[i]-'A'];  // implicit type conversion
    sim += score;
  }
  if(sim < g_min_similarity) {
    return SIM_TYPE_MIN;
  }
  return sim;
}

inline sim_type distancefn(std::string &a, std::string &b) {
  sim_type score;
  sim_type dist = 0;
  for(unsigned int i = 0; i < g_subword_length; i++) {
    score = g_distance_matrix[a[i]-'A'][b[i]-'A'];
    dist += score;
  }
  if(dist > g_max_distance) {
    return SIM_TYPE_MAX;
  }
  return dist;
}


bool validate_character(char c);

float find_best_distance_threshold(float original_g_theta);

#endif

