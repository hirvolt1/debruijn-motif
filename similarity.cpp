#include "similarity.hpp"
#include <iostream>
#include <list>


#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define db_announce_thresholds

using std::cout;
using std::endl;
using std::list;


int g_scoring_matrix[26][26] = {{0}};
bool g_matches_in_scoring_matrix[26][26] = {{false}};
sim_type g_distance_matrix[26][26] = {{0.0f}};

sim_type g_min_similarity = 0.0;
sim_type g_max_distance = 0.0;

//float g_mismatch_avg = 0.0;
//float g_match_avg = 0.0;
// blosum, columns and rows in order: ABCDEFGHIJKLMNOPQRSTUVWXYZ
//                                     |       |    |     |  | |
// blosum / amino-acid skips BJOUX =  01       9    14    20 2325  

//TODO is duplicate of "of inequality can be wrong)


bool is_real_aa(int i) {
  if(i > -1  && i < 27) {
    if(i != 1 && i != 9 && i != 14 && i != 20 && i != 23 && i != 25) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

bool are_real_aas(int i, int j) {
  if(is_real_aa(i) && is_real_aa(j)) {
    return true;
  } else {
    return false;
  }
}

bool are_real_aas(int i, int j, int k) {
  if(is_real_aa(i) && is_real_aa(j) && is_real_aa(k)) {
     return true;
  } else {
    return false;
  }
}   

bool are_real_aas(int i, int j, int k, int l) {
  if(is_real_aa(i) && is_real_aa(j) && is_real_aa(k) && is_real_aa(l)) {
     return true;
  } else {
    return false;
  }
}   


const int g_blosum50[26][26] = {
  { 5, -2, -1, -2, -1, -3,  0, -2, -1, -2, -1, -2, -1, -1,  0, -1, -1, -2,  1,  0,  0,  0, -3, -1, -2, -1},
  {-2,  6, -3,  6,  1, -4, -1,  0, -4, -4,  0, -4, -3,  5,  0, -2,  0, -1,  0,  0,  0, -3, -5, -1, -3,  1},
  {-1, -3, 13, -4, -3, -2, -3, -3, -2, -2, -3, -2, -2, -2,  0, -4, -3, -4, -1, -1,  0, -1, -5, -1, -3, -3},
  {-2,  6, -4,  8,  2, -5, -1, -1, -4, -4, -1, -4, -4,  2,  0, -1,  0, -2,  0, -1,  0, -4, -5, -1, -3,  1},
  {-1,  1, -3,  2,  6, -3, -3,  0, -4, -3,  1, -3, -2,  0,  0, -1,  2,  0, -1, -1,  0, -3, -3, -1, -2,  5},
  {-3, -4, -2, -5, -3,  8, -4, -1,  0,  1, -4,  1,  0, -4,  0, -4, -4, -3, -3, -2,  0, -1,  1, -1,  4, -4},
  { 0, -1, -3, -1, -3, -4,  8, -2, -4, -4, -2, -4, -3,  0,  0, -2, -2, -3,  0, -2,  0, -4, -3, -1, -3, -2},
  {-2,  0, -3, -1,  0, -1, -2, 10, -4, -3,  0, -3, -1,  1,  0, -2,  1,  0, -1, -2,  0, -4, -3, -1,  2,  0},
  {-1, -4, -2, -4, -4,  0, -4, -4,  5,  4, -3,  2,  2, -3,  0, -3, -3, -4, -3, -1,  0,  4, -3, -1, -1, -3},
  {-2, -4, -2, -4, -3,  1, -4, -3,  4,  4, -3,  4,  2, -4,  0, -3, -3, -3, -3, -1,  0,  2, -2, -1, -1, -3},
  {-1,  0, -3, -1,  1, -4, -2,  0, -3, -3,  6, -3, -2,  0,  0, -1,  2,  3,  0, -1,  0, -3, -3, -1, -2,  1},
  {-2, -4, -2, -4, -3,  1, -4, -3,  2,  4, -3,  5,  3, -4,  0, -4, -2, -3, -3, -1,  0,  1, -2, -1, -1, -3},
  {-1, -3, -2, -4, -2,  0, -3, -1,  2,  2, -2,  3,  7, -2,  0, -3,  0, -2, -2, -1,  0,  1, -1, -1,  0, -1},
  {-1,  5, -2,  2,  0, -4,  0,  1, -3, -4,  0, -4, -2,  7,  0, -2,  0, -1,  1,  0,  0, -3, -4, -1, -2,  0},
  { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  {-1, -2, -4, -1, -1, -4, -2, -2, -3, -3, -1, -4, -3, -2,  0, 10, -1, -3, -1, -1,  0, -3, -4, -1, -3, -1},
  {-1,  0, -3,  0,  2, -4, -2,  1, -3, -3,  2, -2,  0,  0,  0, -1,  7,  1,  0, -1,  0, -3, -1, -1, -1,  4},
  {-2, -1, -4, -2,  0, -3, -3,  0, -4, -3,  3, -3, -2, -1,  0, -3,  1,  7, -1, -1,  0, -3, -3, -1, -1,  0},
  { 1,  0, -1,  0, -1, -3,  0, -1, -3, -3,  0, -3, -2,  1,  0, -1,  0, -1,  5,  2,  0, -2, -4, -1, -2,  0},
  { 0,  0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1,  2,  5,  0,  0, -3, -1, -2, -1},
  { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  { 0, -3, -1, -4, -3, -1, -4, -4,  4,  2, -3,  1,  1, -3,  0, -3, -3, -3, -2,  0,  0,  5, -3, -1, -1, -3},
  {-3, -5, -5, -5, -3,  1, -3, -3, -3, -2, -3, -2, -1, -4,  0, -4, -1, -3, -4, -3,  0, -3, 15, -1,  2, -2},
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0, -1, -1, -1, -1, -1,  0, -1, -1, -1, -1, -1},
  {-2, -3, -3, -3, -2,  4, -3,  2, -1, -1, -2, -1,  0, -2,  0, -3, -1, -1, -2, -2,  0, -1,  2, -1,  8, -2},
  {-1,  1, -3,  1,  5, -4, -2,  0, -3, -3,  1, -3, -1,  0,  0, -1,  4,  0,  0, -1,  0, -3, -2, -1, -2,  5},
};

// blosum, columns and rows in order: ABCDEFGHIJKLMNOPQRSTUVWXYZ
const int g_blosum62[26][26] = {
  { 4, -2,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -1, -2,  0, -1, -1, -1,  1,  0,  0,  0, -3, -1, -2, -1},
  {-2,  4, -3,  4,  1, -3, -1,  0, -3, -3,  0, -4, -3,  4,  0, -2,  0, -1,  0, -1,  0, -3, -4, -1, -3,  0},
  { 0, -3,  9, -3, -4, -2, -3, -3, -1, -1, -3, -1, -1, -3,  0, -3, -3, -3, -1, -1,  0, -1, -2, -1, -2, -3},
  {-2,  4, -3,  6,  2, -3, -1, -1, -3, -3, -1, -4, -3,  1,  0, -1,  0, -2,  0, -1,  0, -3, -4, -1, -3,  1},
  {-1,  1, -4,  2,  5, -3, -2,  0, -3, -3,  1, -3, -2,  0,  0, -1,  2,  0,  0, -1,  0, -2, -3, -1, -2,  4},
  {-2, -3, -2, -3, -3,  6, -3, -1,  0,  0, -3,  0,  0, -3,  0, -4, -3, -3, -2, -2,  0, -1,  1, -1,  3, -3},
  { 0, -1, -3, -1, -2, -3,  6, -2, -4, -4, -2, -4, -3,  0,  0, -2, -2, -2,  0, -2,  0, -3, -2, -1, -3, -2},
  {-2,  0, -3, -1,  0, -1, -2,  8, -3, -3, -1, -3, -2,  1,  0, -2,  0,  0, -1, -2,  0, -3, -2, -1,  2,  0},
  {-1, -3, -1, -3, -3,  0, -4, -3,  4,  3, -3,  2,  1, -3,  0, -3, -3, -3, -2, -1,  0,  3, -3, -1, -1, -3},
  {-1, -3, -1, -3, -3,  0, -4, -3,  3,  3, -3,  3,  2, -3,  0, -3, -2, -2, -2, -1,  0,  2, -2, -1, -1, -3},
  {-1,  0, -3, -1,  1, -3, -2, -1, -3, -3,  5, -2, -1,  0,  0, -1,  1,  2,  0, -1,  0, -2, -3, -1, -2,  1},
  {-1, -4, -1, -4, -3,  0, -4, -3,  2,  3, -2,  4,  2, -3,  0, -3, -2, -2, -2, -1,  0,  1, -2, -1, -1, -3},
  {-1, -3, -1, -3, -2,  0, -3, -2,  1,  2, -1,  2,  5, -2,  0, -2,  0, -1, -1, -1,  0,  1, -1, -1, -1, -1},
  {-2,  4, -3,  1,  0, -3,  0,  1, -3, -3,  0, -3, -2,  6,  0, -2,  0,  0,  1,  0,  0, -3, -4, -1, -2,  0},
  { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  {-1, -2, -3, -1, -1, -4, -2, -2, -3, -3, -1, -3, -2, -2,  0,  7, -1, -2, -1, -1,  0, -2, -4, -1, -3, -1},
  {-1,  0, -3,  0,  2, -3, -2,  0, -3, -2,  1, -2,  0,  0,  0, -1,  5,  1,  0, -1,  0, -2, -2, -1, -1,  4},
  {-1, -1, -3, -2,  0, -3, -2,  0, -3, -2,  2, -2, -1,  0,  0, -2,  1,  5, -1, -1,  0, -3, -3, -1, -2,  0},
  { 1,  0, -1,  0,  0, -2,  0, -1, -2, -2,  0, -2, -1,  1,  0, -1,  0, -1,  4,  1,  0, -2, -3, -1, -2,  0},
  { 0, -1, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1,  1,  5,  0,  0, -2, -1, -2, -1},
  { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  { 0, -3, -1, -3, -2, -1, -3, -3,  3,  2, -2,  1,  1, -3,  0, -2, -2, -3, -2,  0,  0,  4, -3, -1, -1, -2},
  {-3, -4, -2, -4, -3,  1, -2, -2, -3, -2, -3, -2, -1, -4,  0, -4, -2, -3, -3, -2,  0, -3, 11, -1,  2, -2},
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0, -1, -1, -1, -1, -1,  0, -1, -1, -1, -1, -1},
  {-2, -3, -2, -3, -2,  3, -3,  2, -1, -1, -2, -1, -1, -2,  0, -3, -1, -2, -2, -2,  0, -1,  2, -1,  7, -2},
  {-1,  0, -3,  1,  4, -3, -2,  0, -3, -3,  1, -3, -1,  0,  0, -1,  4,  0,  0, -1,  0, -2, -2, -1, -2,  4}
};

const int g_pam30[26][26] = {
  {  6,  -3,  -6,  -3,  -2,  -8,  -2,  -7,  -5,  -6,  -7,  -6,  -5,  -4,   0,  -2,  -4,  -7,   0,  -1,   0,  -2, -13,  -1,  -8,  -3},
  { -3,   6, -12,   6,   1, -10,  -3,  -1,  -6,  -8,  -2,  -9, -10,   6,   0,  -7,  -3,  -7,  -1,  -3,   0,  -8, -10,  -1,  -6,   0},
  { -6, -12,  10, -14, -14, -13,  -9,  -7,  -6,  -9, -14, -15, -13, -11,   0,  -8, -14,  -8,  -3,  -8,   0,  -6, -15,  -1,  -4, -14},
  { -3,   6, -14,   8,   2, -15,  -3,  -4,  -7, -10,  -4, -12, -11,   2,   0,  -8,  -2, -10,  -4,  -5,   0,  -8, -15,  -1, -11,   1},
  { -2,   1, -14,   2,   8, -14,  -4,  -5,  -5,  -7,  -4,  -9,  -7,  -2,   0,  -5,   1,  -9,  -4,  -6,   0,  -6, -17,  -1,  -8,   6},
  { -8, -10, -13, -15, -14,   9,  -9,  -6,  -2,  -2, -14,  -3,  -4,  -9,   0, -10, -13,  -9,  -6,  -9,   0,  -8,  -4,  -1,   2, -13},
  { -2,  -3,  -9,  -3,  -4,  -9,   6,  -9, -11, -10,  -7, -10,  -8,  -3,   0,  -6,  -7,  -9,  -2,  -6,   0,  -5, -15,  -1, -14,  -5},
  { -7,  -1,  -7,  -4,  -5,  -6,  -9,   9,  -9,  -7,  -6,  -6, -10,   0,   0,  -4,   1,  -2,  -6,  -7,   0,  -6,  -7,  -1,  -3,  -1},
  { -5,  -6,  -6,  -7,  -5,  -2, -11,  -9,   8,   5,  -6,  -1,  -1,  -5,   0,  -8,  -8,  -5,  -7,  -2,   0,   2, -14,  -1,  -6,  -6},
  { -6,  -8,  -9, -10,  -7,  -2, -10,  -7,   5,   6,  -7,   6,   0,  -6,   0,  -7,  -5,  -7,  -8,  -5,   0,   0,  -7,  -1,  -7,  -6},
  { -7,  -2, -14,  -4,  -4, -14,  -7,  -6,  -6,  -7,   7,  -8,  -2,  -1,   0,  -6,  -3,   0,  -4,  -3,   0,  -9, -12,  -1,  -9,  -4},
  { -6,  -9, -15, -12,  -9,  -3, -10,  -6,  -1,   6,  -8,   7,   1,  -7,   0,  -7,  -5,  -8,  -8,  -7,   0,  -2,  -6,  -1,  -7,  -7},
  { -5, -10, -13, -11,  -7,  -4,  -8, -10,  -1,   0,  -2,   1,  11,  -9,   0,  -8,  -4,  -4,  -5,  -4,   0,  -1, -13,  -1, -11,  -5},
  { -4,   6, -11,   2,  -2,  -9,  -3,   0,  -5,  -6,  -1,  -7,  -9,   8,   0,  -6,  -3,  -6,   0,  -2,   0,  -8,  -8,  -1,  -4,  -3},
  {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
  { -2,  -7,  -8,  -8,  -5, -10,  -6,  -4,  -8,  -7,  -6,  -7,  -8,  -6,   0,   8,  -3,  -4,  -2,  -4,   0,  -6, -14,  -1, -13,  -4},
  { -4,  -3, -14,  -2,   1, -13,  -7,   1,  -8,  -5,  -3,  -5,  -4,  -3,   0,  -3,   8,  -2,  -5,  -5,   0,  -7, -13,  -1, -12,   6},
  { -7,  -7,  -8, -10,  -9,  -9,  -9,  -2,  -5,  -7,   0,  -8,  -4,  -6,   0,  -4,  -2,   8,  -3,  -6,   0,  -8,  -2,  -1, -10,  -4},
  {  0,  -1,  -3,  -4,  -4,  -6,  -2,  -6,  -7,  -8,  -4,  -8,  -5,   0,   0,  -2,  -5,  -3,   6,   0,   0,  -6,  -5,  -1,  -7,  -5},
  { -1,  -3,  -8,  -5,  -6,  -9,  -6,  -7,  -2,  -5,  -3,  -7,  -4,  -2,   0,  -4,  -5,  -6,   0,   7,   0,  -3, -13,  -1,  -6,  -6},
  {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
  { -2,  -8,  -6,  -8,  -6,  -8,  -5,  -6,   2,   0,  -9,  -2,  -1,  -8,   0,  -6,  -7,  -8,  -6,  -3,   0,   7, -15,  -1,  -7,  -6},
  {-13, -10, -15, -15, -17,  -4, -15,  -7, -14,  -7, -12,  -6, -13,  -8,   0, -14, -13,  -2,  -5, -13,   0, -15,  13,  -1,  -5, -14},
  { -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1,  -1,  -1},
  { -8,  -6,  -4, -11,  -8,   2, -14,  -3,  -6,  -7,  -9,  -7, -11,  -4,   0, -13, -12, -10,  -7,  -6,   0,  -7,  -5,  -1,  10,  -9},
  { -3,   0, -14,   1,   6, -13,  -5,  -1,  -6,  -6,  -4,  -7,  -5,  -3,   0,  -4,   6,  -4,  -5,  -6,   0,  -6, -14,  -1,  -9,   6},
};

list<int>::iterator ListIteratorAt(list<int> &l, int p) {
  list<int>::iterator it = l.begin();
  for (int i = 0; i < p; i++) {
    it++;
  }
  return it;
}

int GetMedian(list<int> &matrix_values, int n_values) {
  int median = 0;
  matrix_values.sort();
  if (n_values == 0) {
    std::cerr << "ERROR: MATRIX DOES NOT CONTAIN POSITIVE VALUES" << endl;
  } else if (n_values == 1) {
    median = *(matrix_values.begin());
  } else if (n_values % 2 == 0 && n_values) {  // is even
    median = *(ListIteratorAt(matrix_values, (n_values - 2)/2));
  } else {  // is odd
    median = *(ListIteratorAt(matrix_values, (n_values - 1)/2));
  }
  return median;
}

void set_scoring_matrix(const int matrix[26][26]) {
  int i, j, n_values;
  float mismatch_avg = 0.0;
  float match_avg = 0.0;
  float match_count = 0.0;
  float mismatch_count = 0.0;
  float median = 0.0;
  list<int> matrix_values;
  //find median
  n_values = 0;
  for(i = 0; i < 26; i++) for(j = 0; j < 26; j++) {
    if(is_real_aa(i) && is_real_aa(j)) {
      if(matrix[i][j] > 0) {
        n_values++;
        matrix_values.push_back(matrix[i][j]);
      }
    }
  }
  median = (float)GetMedian(matrix_values, n_values);

  for(i = 0; i < 26; i++) for(j = 0; j < 26; j++) {
    g_scoring_matrix[i][j] = matrix[i][j];
    if(are_real_aas(i,j)) {
      if (matrix[i][j] > median) {
        match_avg += matrix[i][j];
        match_count += 1.0;
        g_matches_in_scoring_matrix[i][j] = true;
      } else {
        mismatch_avg += matrix[i][j];
        mismatch_count += 1.0;
        g_matches_in_scoring_matrix[i][j] = false;
      }
    }
  }
  match_avg = match_avg / match_count;
  mismatch_avg = mismatch_avg / mismatch_count;
  g_min_similarity = (g_theta * match_avg + (1-g_theta) * mismatch_avg) * (float)g_subword_length;
#ifdef db_announce_thresholds
  std::cout << "g_min_similarity: " << g_min_similarity << std::endl;
#endif
  return;
}


void set_distance_matrix() {
  // use the similarity matrix to correctly set values to distance matrix
  // distance with self is set to zero
  // d(a,b) = d(b,a) = kahden avg. = tarpeeks?
  // eli d(x,y) = 0.5 * (s(x,x) - s(x,y) + s(y,y) - s(y,x))
  //

  int i, j;
  float mismatch_avg = 0.0;
  float match_avg = 0.0;
  float match_count = 0.0;
  float mismatch_count = 0.0;
  list<int> matrix_values;

  //trivial speedup by going [i = 0..26 and j = i..26], not needed
  for(i = 0; i < 26; i++) for(j = 0; j < 26; j++) {
    g_distance_matrix[i][j] = static_cast<sim_type>(g_scoring_matrix[i][i] 
        + g_scoring_matrix[j][j] - g_scoring_matrix[i][j] 
        - g_scoring_matrix[j][i]) * 0.5;
    //g_distance_matrix[i][j] = static_cast<sim_type>(g_scoring_matrix[i][i]
    //                                                - g_scoring_matrix[i][j]);
    if(are_real_aas(i,j)) {
      if (g_matches_in_scoring_matrix[i][j]) {
        match_avg += g_distance_matrix[i][j];
        match_count += 1.0;
      } else {
        mismatch_avg += g_distance_matrix[i][j];
        mismatch_count += 1.0;
      }
    }
  }
  match_avg = match_avg / match_count;
  mismatch_avg = mismatch_avg / mismatch_count;
  g_max_distance = (g_theta * match_avg + (1-g_theta) * mismatch_avg) 
      * static_cast<float>(g_subword_length);
#ifdef db_announce_thresholds
  std::cout << "g_max_distance: " << g_max_distance << ", ma_avg: "
      << match_avg << ", mi_avg: " << mismatch_avg << std::endl;
#endif
  return;
}

bool check_distance_matrix_triangle_inequality() {
  // this will be O(n^3) but then again we only have a 26 by 26 matrix, not that
  // many operations after all.
  int i, j, k;
  bool first_error = true;

  for(i = 0; i < 26; i++) for(j = 0; j < 26; j++) for (k = 0; k < 26; k++) {
    if (g_distance_matrix[i][j] > g_distance_matrix[i][k] 
        + g_distance_matrix[k][j]) {
      if (are_real_aas(i,j,k)) {
        if (first_error) {
          std::cerr << "triangle inequality does not hold for:" << std::endl;
          first_error = false;
        }
        std::cerr << "i: " << i << ", j: " << j << ", k: " << k << ", d[i][j]: " 
            << g_distance_matrix[i][j] << ", d[i][k]: " << g_distance_matrix[i][k]
            << ", d[k][j]: " << g_distance_matrix[k][j] << std::endl;
      }
    }
  }

  if (!first_error) {
    return false;
  }
  return true;
}


float find_best_distance_threshold(float original_g_theta) {

  srand(time(NULL));
  int randten;

  std::string str_tmp = "";
  int i, j, k, l;
  std::vector<std::string> subwords;
  k = 0;
  l = 0;
  for (i = 0; i < 26; i++) for(j = 0; j < 26; j++) for (k = 0; k < 26; k++) for (l = 0; l < 26; l++) {
    if (are_real_aas(i,j,k,l)) {
      randten = rand() % 10;
      if (randten == 0) {
        str_tmp.push_back(static_cast<char>('A' + i));
        str_tmp.push_back(static_cast<char>('A' + j));
        str_tmp.push_back(static_cast<char>('A' + k));
        str_tmp.push_back(static_cast<char>('A' + l));
    
        subwords.push_back(str_tmp);
        str_tmp = "";
      }
    }
  }
  g_subword_length = 4;

  float best_theta = 0.0f;
  sim_type smallest_value = 0.835;
  sim_type increment = 0.0015;
  sim_type highest_value = 0.845;
  int matching_results, best_matching_results;
  best_matching_results = 0;
  bool similarity_matches, distance_matches;
  int n_subwords = subwords.size();

  for (sim_type theta = smallest_value; theta <= highest_value;
       theta += increment) {
    matching_results = 0;
    g_theta = original_g_theta;
    set_scoring_matrix(g_blosum62);
    g_theta = theta;
    set_distance_matrix();
    std::cout << "testing with" << theta << std::endl;
    int counter = 0;
    for (std::vector<std::string>::iterator it = subwords.begin();
         it != subwords.end(); it++) {
      counter++;
      if (counter % 1000 == 0) {
        std::cout << "Prog: " 
            << static_cast<double>(counter)/static_cast<double>(n_subwords)
            << std::endl;
      }
      for (std::vector<std::string>::iterator it2 = subwords.begin();
           it2 != subwords.end(); it2++) {
        if (distancefn(*it, *it2) == SIM_TYPE_MAX) {
          distance_matches = false; 
        } else {
          distance_matches = true;
        }
        if (similarityfn(*it, *it2) == SIM_TYPE_MIN) {
          similarity_matches = false;
        } else {
          similarity_matches = true;
        }
        if (similarity_matches == distance_matches) {
          matching_results++;
        }
      }
    }
    std::cout << "theta: " << theta << ", " << matching_results << " / "
        << n_subwords * n_subwords << std::endl;
    if (matching_results > best_matching_results) {
      best_matching_results = matching_results;
      best_theta = theta;
    }
  }
  std::cout << "best theta was " << best_theta << ", matching results: "
      << best_matching_results << " / " << subwords.size() * subwords.size()
      << std::endl;
  return best_theta;
}

bool validate_character(char c) {
  //static const char valid[] = { 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, // A...J
  //                              1, 1, 1, 1, 0, 1, 1, 1, 1, 1, // K...T
  //                              0, 1, 1, 1, 1, 0 };           // U...Z
  if(c < 'A' || c > 'Z') return false;
  return true; //valid[c-'A'];
}
