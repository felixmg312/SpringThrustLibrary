/**
 * @file GraphSymmetricMatrix.hpp
 * Implimentation file for treating the Graph as a MTL Matrix
 */
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <fstream>
#include "CME212/Color.hpp"
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"
#include <iostream>
#include "Graph.hpp"

struct NodeData
{
  NodeData() : val{0}, isBoundary(false) {}
  double val;
  bool isBoundary;
};

// HW3: YOUR CODE HERE
// Define a GraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!
using GraphType = Graph<NodeData, char>; //<  DUMMY Placeholder
using NodeType = typename GraphType::node_type;
class GraphSymmetricMatrix
{

public:
  GraphSymmetricMatrix(GraphType *graph) : graph_(graph)
  {
  }
  unsigned int get_size() const
  {
    return graph_->size();
  }
  template <typename VectorIn, typename VectorOut, typename Assign>
  void mult(const VectorIn &v, VectorOut &w, Assign) const
  {
    for (auto it = graph_->node_begin(); it != graph_->node_end(); ++it)
    {
      NodeType node = *it;
      double temp = 0;
      // If i is on the boundary
      if (node.value().isBoundary)
      {
        temp += v[node.index()];
      }
      else
      {
        temp -= node.degree() * v[node.index()];
      }
      //iterate through the adjacent list
      for (auto inc = node.edge_begin(); inc != node.edge_end(); ++inc)
      {

        if (!node.value().isBoundary and !(*inc).node2().value().isBoundary)
        {
          temp += v[(*inc).node2().index()];
        }
      }
      Assign::apply(w[node.index()], temp);
    }
  }

  template <typename Vector>
  mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector>
  operator*(const Vector &v) const
  {
    return {*this, v};
  }

private:
  GraphType *graph_;
};

inline std::size_t size(const GraphSymmetricMatrix A)
{
  return A.get_size() * A.get_size();
}
inline std::size_t num_rows(const GraphSymmetricMatrix A)
{
  return A.get_size();
}
inline std::size_t num_cols(const GraphSymmetricMatrix A)
{
  return A.get_size();
}

namespace mtl
{
  namespace ashape
  {
    template <>
    struct ashape_aux<GraphSymmetricMatrix>
    {
      typedef nonscal type;
    };
  }

  template <>
  struct Collection<GraphSymmetricMatrix>
  {
    typedef double value_type;
    typedef unsigned size_type;
  };
};

/** Remove all the nodes in graph @a g whose posiiton is within Box3D @a bb.
 * @param[in,out] g  The Graph to remove nodes from
 * @param[in]    bb  The BoundingBox, all nodes inside this box will be removed
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 */
void remove_box(GraphType &g, const Box3D &bb)
{
  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
    NodeType node = *it;
    if (bb.contains(node.position()))
    {
      g.remove_node(node);
    }
  }
  return;
}

// HW3: YOUR CODE HERE
// Define NodeColor and NodePosition functors

// Color Functor
struct ColorFn
{
  double max = 0;

  CME212::Color operator()(const NodeType node)
  {
    if (node.value().val <= -0.75)
    {
      return CME212::Color::make_heat(0.1);
    }
    else if (node.value().val <= -0.5)
    {
      return CME212::Color::make_heat(0.2);
    }
    else if (node.value().val <= 0)
    {
      return CME212::Color::make_heat(0.4);
    }
    else if (node.value().val <= 0.4)
    {
      return CME212::Color::make_heat(0.6);
    }
    else if (node.value().val <= 0.75)
    {
      return CME212::Color::make_heat(0.75);
    }

    else if (node.value().val <= 0.9)
    {
      return CME212::Color::make_heat(0.9);
    }
    return CME212::Color::make_heat(1);
  }
};

// Position Functor
struct NodePosition
{
  double xi;
  double yi;
  double ui;
  Point operator()(const NodeType node)
  {
    xi = node.position().x;
    yi = node.position().y;
    ui = node.value().val;
    Point p = Point(xi, yi, ui);
    return p;
  }
};

template <typename Real, typename ColorFn, typename NodePosition>
class visual_iteration : public itl::cyclic_iteration<Real>
{

  typedef itl::cyclic_iteration<Real> super;
  typedef visual_iteration self;

public:
  template <class Vector>
  visual_iteration(GraphType *graph_, CME212::SFML_Viewer *viewer_, mtl::dense_vector<double> *x_, const Vector &r0, int max_iter_, Real tol_, Real atol_ = Real(0), int cycle_ = 100)
      : super(r0, max_iter_, tol_, atol_, cycle_), graph(graph_), viewer(viewer_), x(x_) {}


  //run this in each iteration of finished
  void update_viewer()
  {
    viewer->clear();
    auto node_map = viewer->empty_node_map(*graph);
    viewer->add_nodes(graph->node_begin(), graph->node_end(), ColorFn(), NodePosition(), node_map);
    viewer->add_edges(graph->edge_begin(), graph->edge_end(), node_map);
    for (auto it = graph->node_begin(); it != graph->node_end(); ++it)
    {
      (*it).value().val = (*x)[(*it).index()];
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(250));
  }

  bool finished()
  {
    update_viewer();

    return super::finished();
  }

  template <typename T>
  bool finished(const T &r)
  {
    bool ret = super::finished(r);
    update_viewer();
    return ret;
  }

protected:
  GraphType *graph;
  CME212::SFML_Viewer *viewer;
  mtl::dense_vector<double> *x;
};
