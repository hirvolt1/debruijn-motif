#ifndef DEBRUIJN_PMOTIF_SIMILARITY_TREE_H_
#define DEBRUIJN_PMOTIF_SIMILARITY_TREE_H_

#include <list>

#include "style_macros.h"
#include "similarity_finder_typedefs.hpp"

class SimilarityNode {
 public:
  SimilarityNode(SFInputCellType &content) : content_(content) {}
  ~SimilarityNode() {}

  // public methods
  
  // 1 parameter constructor and copy constructor needed for list.push_back()
  SimilarityNode(const SimilarityNode &n) : content_(n.content_) {
    // content_ = n.content_;
    children_ = n.children_;
  }

  void operator=(const SimilarityNode &n) {
    content_ = n.content_;
    children_ = n.children_;
  }

  // public data members
  SFInputCellType &content_;
  std::list<SimilarityNode> children_;
 private:

  // ensure correct use
 

  // empty constructor not allowed...
  SimilarityNode();
};  

#endif  // DEBRUIJN_PMOTIF_SIMILARITY_TREE_H_

