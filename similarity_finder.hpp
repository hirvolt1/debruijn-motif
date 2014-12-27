#ifndef DEBRUIJN_PMOTIF_SIMILARITY_FINDER_H_
#define DEBRUIJN_PMOTIF_SIMILARITY_FINDER_H_

#include <iostream>
#include <vector>
#include <string>
#include "style_macros.h"
#include "similarity_finder_typedefs.hpp"
#include "similarity_tree.hpp"

class SimilarityFinder {
 public:
  SimilarityFinder(double distance_threshold,
                   double branching_initial_distance,
                   double branching_divider) {
    distance_threshold_ = distance_threshold;
    branching_initial_distance_ = branching_initial_distance;
    branching_divider_ = branching_divider;
    steps_taken_by_preprocess_ = 0;
    similar_pairs_found_ = 0;
    depth_of_best_parent_candidate_ = 0;
    distance_of_best_parent_candidate_ = branching_initial_distance_;
    work_is_adjust_weights_ = false;
  }
  ~SimilarityFinder() {}

  // public methods
  int AddNodes(SFInputContainerType &nodes);
  int Process(SFInputContainerType &nodes);
  int ProcessNew(SFInputContainerType &nodes);
  double get_steps_taken_by_preprocess() {
    return steps_taken_by_preprocess_;
  }
  double get_similar_pairs_found() {
    return similar_pairs_found_;
  }
  void set_work_is_adjust_weights() {
    work_is_adjust_weights_ = true;
  }


  // public data members
private:

  // private methods
  // 
  void RecursivelyFindSimiliarNodesAndAdd(SimilarityNode &parent,
      double required_similarity_at_this_depth,
      int this_depth, bool can_add_under);

  void FindSimiliarNodesAndAdd(SimilarityNode &root);

  
  // private data members
  int steps_taken_by_preprocess_;
  int similar_pairs_found_;
  int depth_of_best_parent_candidate_;
  double distance_of_best_parent_candidate_;
  double branching_initial_distance_;
  double branching_divider_;
  double distance_threshold_;
  bool add_nodes_;
  bool do_work_;
  bool work_is_adjust_weights_;
  SFInputCellType *node_ptr_;// = NULL;
  SimilarityNode *current_content_to_be_assigned_under_;// = NULL;
  //int n_nodes_;
  // ensure correct use
  DISALLOW_COPY_AND_ASSIGN(SimilarityFinder);

};
#endif  // DEBRUIJN_PMOTIF_SIMILARITY_FINDER_H_
