/**
 * @file subgraph.hpp
 * Implimentation file for viewing a subgraph from our Graph
 */


#include <fstream>
#include <iterator>
#include <iostream>
#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"

#include "Graph.hpp"

/** An iterator that skips over elements of another iterator based on whether
 * those elements satisfy a predicate.
 *
 * Given an iterator range [@a first, @a last) and a predicate @a pred,
 * this iterator models a filtered range such that all i with
 * @a first <= i < @a last and @a pred(*i) appear in order of the original range.
 */
template <typename Pred, typename It>
class filter_iterator : private equality_comparable<filter_iterator<Pred,It>>
{
 public:
  // Get all of the iterator traits and make them our own
  using value_type        = typename std::iterator_traits<It>::value_type;
  using pointer           = typename std::iterator_traits<It>::pointer;
  using reference         = typename std::iterator_traits<It>::reference;
  using difference_type   = typename std::iterator_traits<It>::difference_type;
  using iterator_category = typename std::input_iterator_tag;

  // Constructor
  filter_iterator(const Pred& p, const It& first, const It& last)
      : p_(p), it_(first), end_(last) {

    // HW1 #4: YOUR CODE HERE
    //std::cout<< "out of while" << std::endl;
    while (it_ !=end_ and (!p_(*it_)))
    {
      //std::cout << "in while loop" << std::endl;

          ++it_;
      std::cout << (*it_).index() << std::endl;

    //std::cout << "after ++i" << std::endl;
   }
    //std::cout << "after while" << std::endl;
  }

  // HW1 #4: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
   value_type operator*() const{
     
      return *it_;
   }
   filter_iterator& operator++(){
     do{
       ++(this->it_);
     }
     while(!this->p_(*this->it_) and this->it_!=this->end_);

     return *this;
   }
   bool operator==(const filter_iterator &filit) const{
      return this->it_ == filit.it_;
   }
 

 private:
  Pred p_;
  It it_;
  It end_;
};

/** Helper function for constructing filter_iterators. This deduces the type of
 * the predicate function and the iterator so the user doesn't have to write it.
 * This also allows the use of lambda functions as predicates.
 *
 * Usage:
 * // Construct an iterator that filters odd values out and keeps even values.
 * std::vector<int> a = ...;
 * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
 */
template <typename Pred, typename Iter>
filter_iterator<Pred,Iter> make_filtered(const Iter& it, const Iter& end,
                                         const Pred& p) {
  return filter_iterator<Pred,Iter>(p, it, end);
}

// HW1 #4: YOUR CODE HERE
// Specify and write an interesting predicate on the nodes.
// Explain what your predicate is intended to do and test it.
// If you'd like you may create new nodes and tets files.

/** Test predicate for HW1 #4 */
struct OddPredicate{
  template <typename NODE>
  bool operator()(const NODE &n)
  {
    return n.index() %2 ==0;
  }
};


struct SlicePredicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
    return n.position().x < 0;
  }
};



 