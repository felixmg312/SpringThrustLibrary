#pragma once
/** @file SpaceSearcher.hpp
 * @brief Define the SpaceSearcher class for making efficient spatial searches.
 */

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"
#include "MortonCoder.hpp"
#include <thrust/tuple.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

/** @class SpaceSearcher
 * @brief Class for making spatial searches, which uses the MortonCoder
 *        class as a backend.
 *
 * Given a range of data items and a mapping between these
 * data items and Points, the SpaceSearcher class can be used to quickly
 * iterate over data items which are contained (or almost conatined) inside
 * any given BoundingBox.
 *
 * See "space_search_test.cpp" for a usage example.
 */

// using GraphType = Graph<NodeData, EdgeData>;
// using Node = typename GraphType::node_type;
// GraphType graph;

// #ifndef k
// using GraphType = Graph<NodeData,EdgeData>;
// #define k
// #endif

template <typename T, int L = 7>
class SpaceSearcher
{
private:
  // Implementation types

  /** The number of levels in the MortonCoder. This controls the "accuracy" of
   * the searching iterators (10 -- high accuracy, 1 -- low accuracy), but
   * can also impose more expensive searching.
   */
  static constexpr int NumLevels = L;
  /** Type of MortonCoder. */
  using MortonCoderType = MortonCoder<NumLevels>;
  /** Type of the Morton codes. */
  using code_type = typename MortonCoderType::code_type;
  /** Helper struct: (code_type,T) pair */
  struct morton_pair;

public:
  ////////////////////////////////////
  // TYPE DEFINITIONS AND CONSTANTS //
  ////////////////////////////////////

  /** The type of values that are stored and returned in the spacial search */
  using value_type = T;

  /** Type of iterators, which iterate over items inside a BoundingBox. */
  struct NeighborhoodIterator;

  /** Synonym for NeighborhoodIterator */
  using iterator = NeighborhoodIterator;
  using const_iterator = NeighborhoodIterator;

public:
  /////////////////
  // CONSTRUCTOR //
  /////////////////

  /** @brief SpaceSearcher Constructor.
   *
   * For a range of data items of type @a T given by [@a first, @a last)
   * and a function object @a t2p that maps between data items and @a Points, we
   * arrange the data along a space filling curve which allows all of the data
   * contained withing a given bounding box to be iterated over in less than
   * linear time.
   *
   * @param[in] bb      The "parent" bounding box which this SpaceSearcher
   *                      functions within. All data and queries should
   *                      be contained within.
   * @param[in] t_begin Iterator to first data item of type @a T.
   * @param[in] t_end   Iterator to one past the last data item.
   * @param[in] t2p     A functor that maps data items to @a Points.
   *                      Provides an interface equivalent to
   *                        Point t2p(const T& t) const
   *
   * @pre For all i in [@a first,@a last), @a bb.contains(@a t2p(*i)).
   */

  // template<typename Node>
  // struct T2Point{
  //   __host__ __device__
  //   Point operator()(Node n){
  //     return n.position();
  //   }

  // };

  template <typename TIter, typename T2Point>
  SpaceSearcher(const Box3D &bb, TIter first, TIter last, T2Point t2p) : SpaceSearcher(bb, first, last, thrust::make_transform_iterator<T2Point, TIter>(first, t2p), thrust::make_transform_iterator<T2Point, TIter>(last, t2p))
  {
  }

  /** @brief SpaceSearcher Constructor.
   *
   * For a range of data items of type @a T given by [@a tfirst, @a tlast)
   * and a corresponding range of @a Points given by [@a pfirst, @a plast),
   * we arrange the data along a space filling curve which allows all of the
   * data contained withing a given bounding box to be iterated over in less
   * than linear time.
   *
   * @param[in] bb      The "parent" bounding box which this SpaceSearcher
   *                      functions within. All data and queries should
   *                      be contained within.
   * @param[in] tfirst  Iterator to first data item of type T.
   * @param[in] tlast   Iterator to one past the last data item.
   * @param[in] pfirst  Iterator to first Point corresponding to the position
   *                      of the first data item, *tfirst.
   * @param[in] tlast   Iterator to one past the last @a Point.
   *
   * @pre std::distance(tfirst,tlast) == std::distance(pfirst,plast).
   * @pre For all i in [@a pfirst,@a plast), bb.contains(*i).
   */

  struct Space_Helper
  {
    __host__ __device__
    Space_Helper(MortonCoderType mc_) : mc{mc_} {}
    code_type operator()(Point p)
    {
      return mc.code(p);
    }

  public:
    MortonCoderType mc;
  };

  template <typename TIter, typename PointIter>
  SpaceSearcher(const Box3D &bb, TIter tfirst, TIter tlast, PointIter pfirst, PointIter plast) : mc_(bb)
  {
    /// Without parallel
    // (void) plast;
    // z_data_ = std::vector<morton_pair>();
    // for (TIter it= tfirst; it != tlast; ++it)
    // {
    //   z_data_.push_back(thrust::tuple<code_type, T>(mc_.code(*pfirst), *it));
    //   std::cout<< z_data_.size() << std::endl;
    //   ++pfirst;
    // }

    // With parallel
    Space_Helper sp = Space_Helper(mc_);
    thrust::transform_iterator<Space_Helper, PointIter> code_first(pfirst, sp);
    thrust::transform_iterator<Space_Helper, PointIter> code_last(plast, sp);
    typedef thrust::transform_iterator<Space_Helper, PointIter> CodeIter;
    typedef thrust::tuple<CodeIter, TIter> IteratorTuple;
    typedef thrust::zip_iterator<IteratorTuple> ZipIterator;
    ZipIterator iter_begin(thrust::make_tuple(code_first, tfirst));
    ZipIterator iter_end(thrust::make_tuple(code_last, tlast));
    z_data_ = std::vector<morton_pair>(iter_begin, iter_end);
    thrust::sort(thrust::omp::par, z_data_.begin(), z_data_.end(), [](morton_pair a, morton_pair b)
                 { return a.code_ < b.code_; });
    assert(std::is_sorted(z_data_.begin(), z_data_.end(), [](morton_pair a, morton_pair b)
                          { return a.code_ < b.code_; }));
  }

  ///////////////
  // Accessors //
  ///////////////

  /** The bounding box this SpaceSearcher functions within. */
  Box3D
  bounding_box() const
  {
    return mc_.bounding_box();
  }

  //////////////
  // Iterator //
  //////////////

  /** @class SpaceSearcher::NeighborhoodIterator
   * @brief NeighborhoodIterator class for data items. A forward iterator.
   *
   * Iterates over data items of type @a T contained
   * within epsilon of a given bounding box.
   */
  struct NeighborhoodIterator
  {
    using value_type = T;
    using pointer = T *;
    using reference = T &;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;

    // Default constructor
    NeighborhoodIterator() = default;

    // Iterator operators
    const value_type &operator*() const
    {
      return (*i_).value_;
    }
    NeighborhoodIterator &operator++()
    {
      //    std::cout << (*i_).value_.position()<< std::endl;
      ++i_;
      fix();
      return *this;
    }
    bool operator==(const NeighborhoodIterator &other) const
    {
      return i_ == other.i_;
    }
    bool operator!=(const NeighborhoodIterator &other) const
    {
      return !(*this == other);
    }

  private:
    friend SpaceSearcher;
    using MortonIter = typename std::vector<morton_pair>::const_iterator;
    // RI: i_ == end_ || MortonCoderType::is_in_box(*i_, min_, max_)
    MortonIter i_, end_;
    code_type min_, max_;
    NeighborhoodIterator(MortonIter i, MortonIter end,
                         code_type min, code_type max)
        : i_(i), end_(end), min_(min), max_(max)
    {
      fix();
    }
    // @post RI
    void fix()
    {
      while (i_ < end_)
      {
        code_type c = MortonCoderType::advance_to_box(*i_, min_, max_);
        if (c == *i_)
          break;
        i_ = std::lower_bound(i_, end_, c);
      }
    }
  };

  /** Iterator to the first item contained
   *   within some epsilon of a bounding box.
   * @param bb The bounding box to iterate over.
   * @pre bounding_box.contains(bb)
   */
  const_iterator begin(const Box3D &bb) const
  {
    assert(bounding_box().contains(bb));
    code_type morton_min = mc_.code(bb.min());
    code_type morton_max = mc_.code(bb.max());
    auto mit_end = std::lower_bound(z_data_.begin(), z_data_.end(), morton_max);
    return NeighborhoodIterator(z_data_.begin(), mit_end, morton_min, morton_max);
  }

  /** Iterator to one-past-the-last item contained
   *   within some epsilon of a bounding box.
   * @param bb The bounding box to iterate over.
   * @pre bounding_box.contains(bb)
   */
  const_iterator end(const Box3D &bb) const
  {
    assert(bounding_box().contains(bb));
    code_type morton_min = mc_.code(bb.min());
    code_type morton_max = mc_.code(bb.max());
    auto mit_end = std::lower_bound(z_data_.begin(), z_data_.end(), morton_max);
    return NeighborhoodIterator(mit_end, mit_end, morton_min, morton_max);
  }

private:
  // MortonCoder instance associated with this SpaceSearcher.
  MortonCoderType mc_;

  // A (code_type,value_type) pair that can be used as a MortonCode
  struct morton_pair
  {
    code_type code_;
    value_type value_;
    // Cast operator to treat a morton_pair as a code_type in std::algorithms
    operator const code_type &() const { return code_; }
    // HW4: YOUR CODE HERE

    morton_pair(thrust::tuple<code_type, value_type> pairs)
    {
      code_ = thrust::get<0>(pairs);
      value_ = thrust::get<1>(pairs);
    }
  };

  // Pairs of Morton codes and data items of type T.
  // RI: std::is_sorted(z_data_.begin(), z_data_.end())
  std::vector<morton_pair> z_data_;
};
