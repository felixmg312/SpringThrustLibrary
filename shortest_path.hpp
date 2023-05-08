/**
 * @file shortest_path.cpp
 * Implimentation file for using our templated Graph to determine shortest paths.
 */

#include <vector>
#include <fstream>
#include <queue>
#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"

#include "Graph.hpp"

// Define our types
using GraphType = Graph<int,int>;
using NodeType = typename GraphType::node_type;
using NodeIter = typename GraphType::node_iterator;

/** Find the node with the minimum euclidean distance to a point.
 * @param g  The graph of nodes to search.
 * @param point  The point to use as the query.
 * @return An iterator to the node of @a g with the minimun Eucliean
 *           distance to @a point.
 *           graph.node_end() if graph.num_nodes() == 0.
 *
 * @post For all i, 0 <= i < graph.num_nodes(),
 *          norm(point - *result) <= norm(point - g.node(i).position())
 */

struct compare_functor
{
  const Point *point;
  GraphType *graph;
  compare_functor(const Point *point_) : point{point_} {}
  bool operator()(NodeType node1, NodeType node2) const
  {
    Point point1 = node1.position();
    Point point2 = node2.position();
    Point pt1_sub = point1 - *point;
    Point pt2_sub = point2 - *point;
    if (norm(pt1_sub) < norm(pt2_sub))
    {
      return true;
    }
    else
    {
      return false;
    }
  }
};

NodeIter nearest_node(const GraphType &g, const Point &point)
{

  // HW1 #3: YOUR CODE HERE
  // Quiet compiler warning
  return std::min_element(g.node_begin(), g.node_end(), compare_functor(&point));
}

/** Update a graph with the shortest path lengths from a root node.
 * @param[in,out] g     Input graph
 * @param[in,out] root  Root node to start the search.
 * @return The maximum path length found.
 *
 * @post root.value() == 0
 * @post Graph has modified node values indicating the minimum path length
 *           to the root.
 * @post Graph nodes that are unreachable from the root have value() == -1.
 *
 * This sets all nodes' value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(GraphType &g, NodeType &root)
{
  // HW1 #3: YOUR CODE HERE
  //(void)g, (void)root; // Quiet compiler warnings

  std::vector<bool> visited(g.size(), false);
  visited[root.index()] = true;
  std::queue<NodeType> que;
  que.push(root);
  int max = 0;
  //Using BFS to find the shortest path
  while (!que.empty())
  {
    NodeType node = que.front();
    que.pop();
    int path_length = node.value();
    if (path_length > max)
    {
      max = path_length;
    }
    GraphType::IncidentIterator it; 
    for (it = node.edge_begin(); it != node.edge_end(); ++it)
    {
      //std::cout << (*it).node1().index() << std::endl;
      //std::cout << (*it).node2().index() << std::endl;

      if (!visited[(*it).node2().index()])
      {
        que.push((*it).node2());
        visited[(*it).node2().index()] = true;
        (*it).node2().value()=path_length+1;
      }
    }
  }

  for (unsigned int i = 0; i < visited.size(); i++)
  {
    if (!visited[i])
    {
      g.node(i).value() = -1;
    }
  }
  return max;
}
