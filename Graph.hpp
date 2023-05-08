#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <set>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/iterator/counting_iterator.h"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename v>
void swap(v &val1, v &val2)
{
  v temp = val1;
  val1 = val2;
  val2 = temp;
}

void swap_vector(std::vector<unsigned> &vec, unsigned index1, unsigned index2)
{
  unsigned temp = vec[index1];
  vec[index1] = vec[index2];
  vec[index2] = temp;
}
template <typename V, typename E>
class Graph
{
private:
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  unsigned node_size_;
  unsigned edge_size_;
  struct proxyNode;
  struct proxyEdge;
  std::vector<proxyNode> Nodes;
  std::vector<proxyEdge> Edges;
  std::map<std::vector<unsigned>, unsigned> e2ID;
  std::vector<unsigned> i2u_edge;
  std::vector<unsigned> i2u_node;

public:
  //
  // PUBLIC TYPE DEFINITIONS
  //
  // typedef V node_value_type;
  // typedef E edge_value_type;

  /** Type of this graph. */
  using graph_type = Graph;
  using node_value_type = V;
  using edge_value_type = E;

  /** Predeclaration of Node type. */
  // class Node;
  class Node;

  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  struct NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
  {
    node_size_ = 0;
    edge_size_ = 0;
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node>
  {
  public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */

    Node()
    {
      this->graph_ = nullptr;
      this->id = 0;
    }

    Point &position()
    {
      return graph_->Nodes[id].pt;
    }
    /** Return this node's position. */
    const Point &position() const
    {
      return graph_->Nodes[id].pt;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const
    {

      return graph_->Nodes[id].uid;
    }

    node_value_type &value()
    {
      return graph_->Nodes[id].val;
    }
    const node_value_type &value() const
    {
      return graph_->Nodes[id].val;
    }
    size_type degree() const
    {
      return graph_->Nodes[id].neighbors.size();
    }
    incident_iterator edge_begin() const
    {
      return incident_iterator(this->graph_, this->id, 0);
    }
    incident_iterator edge_end() const
    {
      return incident_iterator(this->graph_, this->id, this->degree());
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const
    {
      // HW0: YOUR CODE HERE
      return this->graph_ == n.graph_ and this->id == n.id;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node &n) const
    {
      if (this->graph_ != n.graph_)
      {
        return this->graph_ < n.graph_;
      }
      return this->id < n.id;
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph *graph_;
    size_type id;
    Node(const Graph *graph, size_type id) : graph_{const_cast<Graph *>(graph)}, id{id} {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const
  {
    return i2u_node.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const
  {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point &position, const node_value_type &v = node_value_type())
  {
    std::vector<size_type> neighbors = std::vector<size_type>();
    node_size_ += 1;
    proxyNode proxynode = proxyNode(Nodes.size(), position, v, neighbors);
    Nodes.push_back(proxynode);
    i2u_node.push_back(Nodes.size() - 1);
    return Node(this, Nodes.size() - 1);
  }

  /** Remove a node in the graph, returning a boolean value, 0 if failed remove and 1 if successfull removal.
   * @param[in] node The node that we will want to remove
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() <= old num_edges() if there are edges attached to the node, then we will also remove the edges
   * Complexity: O(1) amortized operations.
   */
  size_type remove_node(const Node &node)
  {
    if (!has_node(node))
    {
      return 0;
    }
    node_size_ -= 1;
    std::vector<size_type> junk; // this is a vector to store the node id that's gonna get removed
    size_type neighborsize = Nodes[node.id].neighbors.size();
    Node a = Node(this, node.id);
    for (size_type i = 0; i < neighborsize; ++i)
    {

      Node b = Node(this, Nodes[node.id].neighbors[i]);
      junk.push_back(b.id);
    }
    for (size_type i = 0; i < junk.size(); i++)
    {

      Node b = Node(this, junk[i]);

      this->remove_edge(a, b);
    }
    Nodes[node.id].isvalid = false;
    unsigned last_index = Nodes[i2u_node[i2u_node.size() - 1]].index;
    swap_vector(i2u_node, Nodes[node.id].uid, i2u_node.size() - 1);
    swap(Nodes[node.id].uid, Nodes[last_index].uid);
    i2u_node.pop_back();
    return 1;
  }

  /** Remove a node in the graph, returning a boolean value, 0 if failed remove and 1 if successfull removal.
   * @param[in] node The iterator of the node that we will want to remove
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() <= old num_edges() if there are edges attached to the node, then we will also remove the edges
   * Complexity: O(1) amortized operations.
   */
  node_iterator remove_node(node_iterator n_it)
  {
    Node n = *n_it;
    if (remove_node(n))
    {
      return n_it;
    }
    else
    {
      return this->node_begin();
    }
  }
  /* Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const
  {
    return Nodes[n.id].isvalid and n.id >= 0 and n.id <= Nodes.size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  // construct a node with internal index
  // return node with internal id
  Node node(size_type i) const
  {

    if (i >= i2u_node.size())
    {
      return Node();
    }
    else
    {
      return Node(this, i2u_node[i]);
    }
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Node>
  {
  public:
    /** Construct an invalid Edge. */
    Edge()
    {
      this->graph_ = nullptr;
      this->id = 0;
    }

    const edge_value_type &value() const
    {
      return graph_->Edges[id].val;
    }

    edge_value_type &value()
    {
      return graph_->Edges[id].val;
    }

    double length() const
    {
      return norm(node1().position() - node2().position());
    }

    /** Return a node of this Edge */
    Node node1() const
    {
      size_type node_id = graph_->Edges[this->id].node1;

      return Node(graph_, node_id);
    }

    /** Return the other node of this Edge */
    Node node2() const
    {
      size_type node_id = graph_->Edges[id].node2;
      return Node(graph_, node_id);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const
    {

      return this->id == e.id and this->graph_ == e.graph_;
    }

    /** Test whether this edge and @a e are not equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator!=(const Edge &e) const
    {

      return this->id != e.id or this->graph_ != e.graph_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const
    {
      if (this->graph_ != e.graph_)
      {
        return this->graph_ < e.graph_;
      }
      else
      {
        return this->id < e.id;
      }
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph *graph_;
    size_type id;
    Edge(const Graph *graph, size_type id) : graph_{const_cast<Graph *const>(graph)}, id{id} {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const
  {
    return i2u_edge.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const
  {
    return Edge(this, i2u_edge[i]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const
  {
    if (!has_node(a) or !has_node(b))
    {
      return false;
    }

    std::vector<size_type> node_set = {a.id, b.id};

    auto search = e2ID.find(node_set);

    if (search != e2ID.end() and Edges[search->second].isvalid)
    {
      return true;
    }
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node &a, const Node &b, const edge_value_type &v = edge_value_type())
  {

    size_type converted_id1 = a.id;
    size_type converted_id2 = b.id;

    std::vector<size_type> key1 = {converted_id1, converted_id2};
    std::vector<size_type> key2 = {converted_id2, converted_id1};

    if (has_edge(a, b))
    {
      size_type id = e2ID.find(key1)->second;
      if (Edges[id].node1 != converted_id1)
      {
        swap(Edges[id].node1, Edges[id].node2);
      }
      return Edge(this, id);
    }

    else
    {
      this->edge_size_ += 1;
      proxyEdge proxyedge = proxyEdge(converted_id1, converted_id2, Edges.size(), v);
      e2ID[key1] = Edges.size();
      e2ID[key2] = Edges.size();
      Edges.push_back(proxyedge);

      Nodes[converted_id1].neighbors.push_back(converted_id2);
      Nodes[converted_id2].neighbors.push_back(converted_id1);
      i2u_edge.push_back(Edges.size() - 1);
      return Edge(this, Edges.size() - 1);
    }
  }
  /** Remove an edge of the graph, return 0 if failed remove and return 1 if successfully reemoved
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return 1 if successfully removed or 0 if failed remove
   * @post has_edge(@a a, @a b) == true
   * @post new num_edges() == old num_edges()-1.
   *       will also remove the adjacent list of the corresponding @a a, @a b nodes
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type remove_edge(const Node &a, const Node &b)
  {

    if (!has_edge(a, b))
    {
      return 0;
    }

    std::vector<size_type> key1 = {a.id, b.id};
    std::vector<size_type> key2 = {b.id, a.id};
    auto id1 = e2ID.find(key1)->second;
    unsigned last_index = Edges[i2u_edge[i2u_edge.size() - 1]].index;
    Edges[id1].isvalid = false;
    swap_vector(i2u_edge, Edges[id1].uid, i2u_edge.size() - 1);
    swap(Edges[id1].uid, Edges[last_index].uid);
    i2u_edge.pop_back();
    e2ID.erase(key1);
    e2ID.erase(key2);
    edge_size_ -= 1;

    for (size_type i = 0; i < this->Nodes[a.id].neighbors.size(); i++)
    {
      if (b.id == this->Nodes[a.id].neighbors[i])
      {
        auto temp = this->Nodes[a.id].neighbors.begin() + i;
        this->Nodes[a.id].neighbors.erase(temp);
      }
    }
    for (size_type i = 0; i < this->Nodes[b.index()].neighbors.size(); i++)
    {
      if (a.id == this->Nodes[b.id].neighbors[i])
      {
        auto temp = this->Nodes[b.id].neighbors.begin() + i;
        this->Nodes[b.id].neighbors.erase(temp);
      }
    }

    return 1;
  }

  /** Remove an edge of the graph, return the iterator if the removal is sucessfull, else return the edge begin iterator
   * @pre edge e is the iterator of a valid edge
   * @return 1 if successfully removed or 0 if failed remove
   * @post has_edge(@a a, @a b) == true
   * @post new num_edges() == old num_edges()-1.
   *       will also remove the corresponding neighbors of the adjacent list of the corresponding @a a, @a b nodes
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type remove_edge(const Edge &e)
  {
    return remove_edge(e.node1(), e.node2());
  }

  edge_iterator remove_edge(edge_iterator e_it)
  {
    Edge e = *e_it;
    if (remove_edge(e))
    {
      return e_it;
    }
    else
    {
      return this->edge_begin();
    }
  }
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    Nodes.clear();
    Edges.clear();
    e2ID.clear();
    i2u_edge.clear();
    i2u_node.clear();
    node_size_ = 0;
    edge_size_ = 0;
  }

  
  // Node Iterator
  

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  // class NodeIterator
  // {
  // public:
  //   // These type definitions let us use STL's iterator_traits.
  //   using value_type = Node;                             // Element type
  //   using pointer = Node *;                              // Pointers to elements
  //   using reference = Node &;                            // Reference to elements
  //   using difference_type = std::ptrdiff_t;              // Signed difference
  //   using iterator_category = std::forward_iterator_tag; // Weak Category, Proxy

  //   /** Construct an invalid NodeIterator. */
  //   NodeIterator()
  //   {
  //   }

  //   Node operator*() const
  //   {

  //     return Node(graph_, graph_->i2u_node[id]);
  //   }
  //   NodeIterator &operator++()
  //   {
  //     if (id == graph_->size())
  //     {
  //       return *this;
  //     }
  //     else
  //     {
  //       this->id++;
  //       return *this;
  //     }
  //   }
  //   bool operator==(const NodeIterator &node) const
  //   {
  //     return this->id == node.id;
  //   }

  //   bool operator!=(const NodeIterator &node) const
  //   {
  //     return this->id != node.id;
  //   }

  // private:
  //   friend class Graph;
  //   Graph *graph_;
  //   size_type id;

  //   NodeIterator(const Graph *graph, size_type id_) : graph_{const_cast<Graph *>(graph)}, id{id_} {}
  // };

  struct id2node
  {
    __host__ __device__
    //overloaded operator
    Node operator()(size_type id) const{
      return graph_->node(id);
    }
  
    // constructor
    Graph *graph_;
    id2node(const graph_type *graph) : graph_{const_cast<Graph *>(graph)} {}
  };
  
 
  struct NodeIterator : thrust::transform_iterator<id2node, thrust::counting_iterator<size_type>,Node>{
    using super_t = thrust::transform_iterator<id2node, thrust::counting_iterator<size_type>,Node>;
    NodeIterator(){} //Default Invalid counstructor
    NodeIterator(const super_t& ti): super_t{ti}{}
    private:
      friend class Graph;
      NodeIterator(const graph_type *g, size_type id_) :super_t(thrust::counting_iterator<size_type>(id_), id2node(g)) {}
  };
  // using node_iterator=NodeIterator;

  node_iterator node_begin() const
  {
    return NodeIterator(this, 0);
  }
  node_iterator node_end() const
  {
    return NodeIterator(this, this->size());
  }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    IncidentIterator()
    {
      this->graph_ = nullptr;
      this->id = 0;
    }

    Edge operator*() const
    {
      Node a = Node(this->graph_, this->Node_id);
      Node b = Node(this->graph_, graph_->Nodes[a.id].neighbors[this->id]);
      std::vector<unsigned int> key = {a.id, b.id};
      size_type edge_id = this->graph_->e2ID.find(key)->second;

      if (graph_->Edges[edge_id].node1 != a.id)
      {
        swap(graph_->Edges[edge_id].node1, graph_->Edges[edge_id].node2);
      }

      return Edge(this->graph_, edge_id);
    }

    IncidentIterator &operator++()
    {
      this->id++;
      return *this;
    }
    bool operator==(const IncidentIterator &inc) const
    {
      return this->Node_id == inc.Node_id and this->id == inc.id;
    }

    bool operator!=(const IncidentIterator &inc) const
    {
      return this->Node_id != inc.Node_id or this->id != inc.id;
    }

  private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph *graph_;
    size_type Node_id;
    size_type id;
    IncidentIterator(const Graph *graph, int Node_id_, size_type id_) : graph_(const_cast<Graph *>(graph)), Node_id(Node_id_), id(id_) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator()
    {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const
    {
      return Edge(this->graph_, this->graph_->i2u_edge[id]);
    }
    EdgeIterator &operator++()
    {
      this->id++;
      return *this;
    }
    bool operator==(const EdgeIterator &edgeit) const
    {
      return this->id == edgeit.id and this->graph_ == edgeit.graph_;
    }

    bool operator!=(const EdgeIterator &edgeit) const
    {
      return this->id != edgeit.id or this->graph_ != edgeit.graph_;
    }

  private:
    friend class Graph;
    Graph *graph_;
    size_type id;
    EdgeIterator(const Graph *graph, size_type id_) : graph_{const_cast<Graph *>(graph)}, id{id_} {}
  };

  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const
  {
    return EdgeIterator(this, 0);
  }
  edge_iterator edge_end() const
  {
    return EdgeIterator(this, this->edge_size_);
  }

private:
  struct proxyNode
  {
    size_type index;
    Point pt;
    node_value_type val;
    std::vector<size_type> neighbors;
    size_type uid;
    bool isvalid;

    proxyNode(size_type index, Point pt, node_value_type val, std::vector<size_type> neighbors)
    {
      this->index = index;
      this->pt = pt;
      this->val = val;
      this->neighbors = neighbors;
      this->uid = index;
      this->isvalid = true;
    }
  };

  struct proxyEdge
  {
    size_type index;
    size_type node1;
    size_type node2;
    edge_value_type val;
    size_type uid;
    bool isvalid;
    proxyEdge(size_type node1, size_type node2, size_type index, edge_value_type val)
    {
      this->node1 = node1;
      this->node2 = node2;
      this->index = index;
      this->val = val;
      this->uid = index;
      this->isvalid = true;
    }
  };
};

#endif // CME212_GRAPH_HPP