#include "similarity_finder.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "work_for_node_pair.hpp"
#include "similarity.hpp"

//#define db_verbose_CalculateDistanceWrapper

//#define db_verbose_CalculateDistance
//#define db_1
//#define db_cout_pairs_found
//#define db_verbose_recursive
//#define db_weight
//#define db_report_tree_size

inline sim_type similarity_real(const char *a, const char *b) {
  sim_type score;
  sim_type sim = 0;
  for(unsigned int i = 0; i < g_subword_length; i++) {
    score = g_scoring_matrix[a[i]-'A'][b[i]-'A'];  // implicit type conversion
    sim += score;
  }
  return sim;
}


inline int SingleSymbolDistance(char c1, char c2) {
  return abs(static_cast<int>(c1) - static_cast<int>(c2));  
}

// Will be inherited from another method. Used initially for developing, this
// works only with strings
double CalculateDistance(std::string &node1, std::string &node2) {
  std::string &smaller = node1;
  std::string &larger = node2;  // rename
  // NOTE here might be happening something reaaally nasty. Consider reference
  // to be just an alias: node2 = node1; node1 = node2; //(=node1)
  //Flip if needed
  if (node1.size() > node2.size()) {
    larger = node1;  // by reference... ?
    smaller = node2;
  }  // else they are ok

  int distance = 0;
  for (size_t i = 0; i < smaller.size(); i++) {
    distance += SingleSymbolDistance(smaller[i], larger[i]);
  }

#ifdef db_verbose_CalculateDistance
  std::cout << "n1: " << node1 << " , n2: " << node2 << ", d: " 
      << static_cast<double>(distance) << std::endl;
#endif

  return static_cast<double>(distance);
}

sim_type CalculateDistanceWrapper(sw_node* node1, sw_node* node2) {
  const char *a = node1->subword;
  const char *b = node2->subword;
  sim_type score;
  sim_type dist = 0;
  for(unsigned int i = 0; i < g_subword_length; i++) {
    score = g_distance_matrix[a[i]-'A'][b[i]-'A'];
    dist += score;
  }
  //if(dist > g_max_distance) {
  //  return SIM_TYPE_MAX;
  //}
#ifdef db_verbose_CalculateDistanceWrapper
  if (dist <= g_max_distance) {
    std::cout << "a: " << node1->seq << ":" << a << ", b: " << node2->seq << ":"
        << b << ", d: " << dist << std::endl;
  }
#endif

  return dist;
}

int SimilarityFinder::Process(SFInputContainerType &nodes) {
  steps_taken_by_preprocess_ = 0;
  similar_pairs_found_ = 0;
  for (SFInputContainerType::iterator node_it = nodes.begin(); 
       node_it != nodes.end(); node_it++) {
    for (SFInputContainerType::iterator node_it2 = node_it + 1; 
         node_it2 != nodes.end(); node_it2++) {
      steps_taken_by_preprocess_++;
      sim_type distance = CalculateDistanceWrapper(*node_it, *node_it2);
      if (distance <= distance_threshold_) {
        if((*node_it)->seq != (*node_it2)->seq) {

          similar_pairs_found_++;
          //#ifdef db_cout_pairs_found
          //std::cout << "Old, pair: " << *node_it << ", " << *node_it2
          //    << std::endl;
          //#endif
        }
        // iterator points to a sw_node-pointer , we want the sw_nodes to get 
        // ... a reference to them :). So pretty.
        //do_work_for_node_pair(*(*(node_it)), *(*(node_it2)), distance);
      }
    }
  }
  return similar_pairs_found_;
}

void SimilarityFinder::RecursivelyFindSimiliarNodesAndAdd(
    SimilarityNode &parent, 
    double required_similarity_at_this_depth,
    int this_depth, bool can_add_under) {

#ifdef db_verbose_recursive
  std::cout << "recursive, parent: " << parent.content_ << std::endl;
#endif

  bool can_add_under_child;
  sim_type new_to_child_distance, new_to_grandkids_min_distance; 
  //sim_type new_to_grandkids_max_distance;
  required_similarity_at_this_depth = required_similarity_at_this_depth
      / branching_divider_;
  this_depth++;

  for (std::list<SimilarityNode>::iterator child_it = parent.children_.begin();
       child_it != parent.children_.end(); child_it++) {
    new_to_child_distance 
        = CalculateDistanceWrapper(child_it->content_, *node_ptr_);
    steps_taken_by_preprocess_++;
    //new_to_grandkids_max_distance = new_to_child_distance 
    //    + required_similarity_at_this_depth;
    new_to_grandkids_min_distance = new_to_child_distance
        - required_similarity_at_this_depth;
    can_add_under_child = can_add_under;

    // TODO if it is possible to do with a dynamic structure: if the child and
    // all its children are fully under/inside the distance_threshold_ range
    // from node "new" they should be all assigned to be similar with the node.
    // This can be done by (hopefully) accessing left (and right) leaves of
    // node in constant time and processing these leaves.

    //if node could be assigned under this child

#ifdef db_1
    std::cout << "new_to_child_distance: " << new_to_child_distance
        << ", required_similarity_at_this_depth: " 
        << required_similarity_at_this_depth 
        << ", this depth " << this_depth 
        << ", depth_of_best_parent_candidate_: " 
        << depth_of_best_parent_candidate_
        << ", new_to_child_distance: "
        << new_to_child_distance
        << ", distance_of_best_parent_candidate_: " 
        << distance_of_best_parent_candidate_ << std::endl;
#endif
    if(add_nodes_) {
      if (new_to_child_distance <= required_similarity_at_this_depth 
          && can_add_under_child) {
        if (this_depth > depth_of_best_parent_candidate_) {
          current_content_to_be_assigned_under_ = &(*child_it);
#ifdef db_1
          std::cout << "going to assign under: " << child_it->content_ 
              << std::endl;
#endif
          distance_of_best_parent_candidate_ = new_to_child_distance;
          depth_of_best_parent_candidate_ = this_depth;
        } else if (this_depth == depth_of_best_parent_candidate_ &&
                   new_to_child_distance < distance_of_best_parent_candidate_) {
          current_content_to_be_assigned_under_ = &(*child_it);
          distance_of_best_parent_candidate_ = new_to_child_distance;
          depth_of_best_parent_candidate_ = this_depth;
        }
      } else {
        can_add_under_child = false;
      }
    }

    // can this child have similar (matching) children?

#ifdef db_verbose_recursive
    std::cout << "child: " << child_it->content_ << " to grandkids min dist: "
        << new_to_grandkids_min_distance << ", threshold: " 
        << distance_threshold_ << std::endl;
#endif
    if (new_to_grandkids_min_distance <= distance_threshold_) {
      // grandkids can contain a match

      // is the child itself a match?
      if (new_to_child_distance <= distance_threshold_ && do_work_) {
        // child is a match
        if(work_is_adjust_weights_) {

        } else
        if(child_it->content_->seq != (*node_ptr_)->seq) {
          if(work_is_adjust_weights_) {
            sim_type sim_value = similarity_real((*node_ptr_)->subword,
                                                 child_it->content_->subword);
            if (sim_value > (*node_ptr_)->weight) {
              (*node_ptr_)->weight = 0;
            } else {
              (*node_ptr_)->weight -= sim_value;
            }

          } else {
            do_work_for_node_pair_by_pointers(*node_ptr_, child_it->content_,
                                              new_to_child_distance);
          }
          //do_work_for_node_pair(*node_ref, *(child_it->content_),
          //                      new_to_child_distance);

#ifdef db_weight
          std::cout << "in recu, node_ref->weight: " << (*node_ptr_)->weight
              << ", child_it->content_->weight: " 
              << child_it->content_->weight << std::endl;
          std::cout << "ics, nref: " << (*node_ptr_)->adj_ic.size() << std::endl;
#endif
          similar_pairs_found_++;
#ifdef db_cout_pairs_found
          std::cout << "New, pair: " << child_it->content_ << ", " << node_ref 
              << std::endl;
#endif
        }
      }
      // proceed to children if any
      if (!child_it->children_.empty()) {
        RecursivelyFindSimiliarNodesAndAdd(*child_it,
                                           required_similarity_at_this_depth,
                                           this_depth, can_add_under_child);
      }
    }


  }
  // TODO: keep track of the parent where to add the new node.
  // TODO: max distance can be used to realize that everything from here down 
  // matches
}


void SimilarityFinder::FindSimiliarNodesAndAdd(SimilarityNode &root) {
  current_content_to_be_assigned_under_ = &root;
  //SFInputCellType &node_ref = *node_it; //node_ref points to type sw_node*

  double required_similarity_at_this_depth = branching_initial_distance_;
  distance_of_best_parent_candidate_ = branching_initial_distance_;
  depth_of_best_parent_candidate_ = 0;
#ifdef db_1
  std::cout << "find and add: " << *node_it << std::endl;
#endif

  RecursivelyFindSimiliarNodesAndAdd(root, required_similarity_at_this_depth,
                                     0, true);
  // add new node with matching content to the correct spot
  if (add_nodes_) {
    SimilarityNode new_node(*node_ptr_);

#ifdef db_1
    std::cout << new_node.content_ << " under " << 
        current_content_to_be_assigned_under->content_ << std::endl;
#endif

    current_content_to_be_assigned_under_->children_.push_back(new_node);
  }
}

int SimilarityFinder::ProcessNew(SFInputContainerType &nodes) {
  //#ifdef db_report_tree_size
  //n_nodes_ = 0;
  //#endif

  steps_taken_by_preprocess_ = 0;
  similar_pairs_found_ = 0;

  // initialize tree structure

  SFInputCellType root_content = new sw_node;  
  // TODO needs modification if this is not a string
  SimilarityNode root(root_content);

  add_nodes_ = true;
  do_work_ = false;
  for (SFInputContainerType::iterator node_it = nodes.begin(); 
       node_it != nodes.end(); node_it++) {
    node_ptr_ = &(*node_it);
    FindSimiliarNodesAndAdd(root);
#ifdef db_weight
    std::cout << "PN (*node_it)->weight: " << (*node_it)->weight << std::endl;
    std::cout << "ics: " << (*node_it)->adj_ic.size() << std::endl;
#endif
  }
#ifdef db_report_tree_size
  std::cerr << "Tree contained " << nodes.size() << " nodes. " << std::endl;
  std::cerr << "Should total " << nodes.size() * (sizeof(SFInputCellType)
                                                  + 2*sizeof(SFInputCellType*));
#endif

  add_nodes_ = false;
  do_work_ = true;

  for (SFInputContainerType::iterator node_it = nodes.begin(); 
       node_it != nodes.end(); node_it++) {
    node_ptr_ = &(*node_it);

    FindSimiliarNodesAndAdd(root);

    sw_node *current_node_p;
    current_node_p = *node_it;

    for(map<int, ic_edge_cont>::iterator it = current_node_p->tmp_ic_edges.begin();
        it != current_node_p->tmp_ic_edges.end();
        it++) {
      for(ic_edge_cont::iterator it2 = it->second.begin();
          it2 != it->second.end();
          it2++) {
        current_node_p->adj_ic.push_back((*it2));
      }
#ifdef db_weight
      std::cout << "tmp: " << std::endl;
      for(ic_edge_cont::iterator it2 = it->second.begin();
          it2 != it->second.end();
          it2++) {
        std::cout << get_identifying_string(*it2) << " ";
      }
      std::cout << std::endl;
#endif
      it->second.clear();
    }
    current_node_p->max_sims.clear();
    current_node_p->tmp_ic_edges.clear();
#ifdef db_weight
    std::cout << "adj_ic: " << std::endl;
    for(ic_edge_cont::iterator it = current_node_p->adj_ic.begin();
          it != current_node_p->adj_ic.end();
          it++) {
      std::cout << get_identifying_string(*it) << " ";
    }
    std::cout << std::endl;
#endif
  }  

  delete root_content; root_content = NULL;
  return similar_pairs_found_;
}
// TODO: some double vs. sim_type fiddling would possibly be cool.
