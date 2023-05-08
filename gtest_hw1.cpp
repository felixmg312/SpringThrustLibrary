#include <gtest/gtest.h>
#include <fstream>

#include "CME212/Util.hpp"
#include "Graph.hpp"
#include "shortest_path.hpp"
#include "subgraph.hpp"


class GraphPointFixture : public ::testing::Test {
 protected:
   //Define types
  using GraphType = Graph<int,int>;
  using NodeType  = typename GraphType::node_type;
  using NodeIter = typename GraphType::node_iterator;
  using EdgeType  = typename GraphType::edge_type;

  //Set up Graph and Points
  GraphType graph;
  std::vector<Point> points;
  virtual void SetUp() {
    for(int i = 0; i < 10; i++)
      points.push_back(Point(i));
  }
  
};

// Test adding node with default value
TEST_F(GraphPointFixture, DefaultNodeVal){
  graph.add_node(points[0]);
  EXPECT_EQ( graph.node(0).value(), 0 ) << "add_node does not intalize node vale with a default 0 value";
}

// Test degree function
TEST_F(GraphPointFixture, Degree){
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  graph.add_edge(n0, n1);

  EXPECT_EQ(n2.degree(),0)  << "n2 degree is 0";
  EXPECT_EQ(n1.degree(), 1) << "n1 degree is 1";
}

// Test node iterator
TEST_F(GraphPointFixture, NodeIter){
  graph.add_node(points[0]);
  graph.add_node(points[1]);
  graph.add_node(points[2]);
  
  int iter = 0;
  for(auto ni = graph.node_begin(); ni != graph.node_end(); ++ni){
    ++iter;
  }
  EXPECT_EQ(iter, 3) << " error in node iteration " ;
}
TEST_F(GraphPointFixture, EdgeIter)
{
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);

  graph.add_edge(n0, n1);
  graph.add_edge(n2, n3);
  int iter = 0;
  for (auto ni = graph.edge_begin(); ni != graph.edge_end(); ++ni)
  {
    ++iter;
  }
  EXPECT_EQ(iter, 2) << " error in node iteration ";
  EXPECT_TRUE(graph.node_end() == graph.node_end()) << " error in node_end ";
}
//Test nearest node
TEST_F(GraphPointFixture, ShortestPath){
  graph.add_node(points[0]);
  graph.add_node(points[1]);
  graph.add_node(points[2]);
  
  NodeIter nearest = nearest_node(graph, Point(0));
  EXPECT_EQ( *nearest, graph.node(0)) << " error finding nearest node " ;
}

TEST_F(GraphPointFixture, NodeIter2)
{
  graph.add_node(points[0]);
  graph.add_node(points[1]);
  graph.add_node(points[2]);
  graph.add_node(points[3]);
  graph.add_node(points[4]);
  graph.add_node(points[5]);

  int iter = 0;
  for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni)
  {
    ++iter;
  }
  EXPECT_EQ(iter, 6) << " error in node iteration ";
}

TEST_F(GraphPointFixture, Remove_edge)
{
  
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);

  graph.add_edge(n0,n1);
  graph.add_edge(n0,n2);
  // for(auto it=n0.edge_begin(); it!= n0.edge_end(); ++it){
  //   std::cout<< (*it).node2().index() << std::endl;
  // } 
  graph.remove_node(n0);
  // graph.add_edge(n0, n1);
  // graph.remove_edge(n0, n1);
 
  EXPECT_EQ(graph.num_edges(), 0);
}

TEST_F(GraphPointFixture, Remove_node)
{
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);


  graph.remove_node(n0);
  graph.remove_node(n1);
  graph.remove_node(n2);
  graph.remove_node(n3);

  EXPECT_EQ(graph.size(), 0);
}

TEST_F(GraphPointFixture, Remove_node1)
{
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  graph.add_edge(n0, n1);
  graph.add_edge(n0, n2);
  graph.add_edge(n1,n2);
  graph.remove_node(n0);

  EXPECT_EQ(graph.num_edges(), 1);
}


TEST_F(GraphPointFixture, Remove_node2)
{
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);

  graph.add_edge(n0, n1);
  graph.add_edge(n0, n2);
  graph.add_edge(n0, n3);
  graph.add_edge(n1,n3);
  // for (auto it = graph.node_begin(); it != graph.node_end(); ++it)
  // {
  //   std::cout << (*it).index() << std::endl;
    
  // }
  // for (auto it =n0.edge_begin(); it != n0.edge_end(); ++it)
  // {
  //   std::cout << (*it).node2().index() << std::endl;
  // }
  graph.remove_node(n3);

  // for(auto it = graph.edge_begin(); it!=graph.edge_end(); ++it ){
  //   std::cout<< (*it).node1().index()<<std::endl;
  //   std::cout << (*it).node2().index() << std::endl;
  // }

  unsigned count_edges=0;
  for (unsigned k = 0; k < graph.num_nodes(); ++k)
  {
    for (unsigned j = k + 1; j < graph.num_nodes(); ++j)
    {
      if (graph.has_edge(graph.node(k), graph.node(j)))
      {
        ++count_edges;
         // std::cout << graph.node(k).index() << "next" << graph.node(j).index() << std::endl;
        // std::cout<< " in here " << count_edges<<std::endl;
      }
    }
  }
  EXPECT_EQ(graph.num_edges(), count_edges);
}





TEST_F(GraphPointFixture, Remove_node3)
{
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);
  NodeType n4 = graph.add_node(points[4]);
  NodeType n5 = graph.add_node(points[5]);
  NodeType n6 = graph.add_node(points[6]);
  NodeType n7= graph.add_node(points[7]);

  
  graph.add_edge(n4, n3);
  graph.add_edge(n0, n3);
  graph.add_edge(n1, n3);
  graph.add_edge(n1, n4);

  graph.add_edge(n2, n3);
  graph.add_edge(n0, n1);
  graph.add_edge(n0, n2);
  graph.add_edge(n5, n4);
  graph.add_edge(n6, n7);
  
  //-----------------------Debug-------------------------------
  // for (auto it =n0.edge_begin(); it != n0.edge_end(); ++it)
  // {
  //   std::cout << (*it).node2().index() << std::endl;
  // }
  // graph.remove_node(graph.node(0));
  // for (auto it = n1.edge_begin(); it != n1.edge_end(); ++it)
  // {
  //   std::cout << (*it).node2().index() << std::endl;
  // }
  // graph.remove_node(graph.node(1));

  // for (auto it = n1.edge_begin(); it != n1.edge_end(); ++it)
  // {
  //   std::cout << (*it).node2().index() << std::endl;
  // }
  //-----------------------Debug-------------------------------

  graph.remove_node(graph.node(1));

  graph.remove_node(graph.node(2));
  graph.remove_node(graph.node(3));
  graph.remove_node(graph.node(1));
  graph.remove_node(graph.node(1));

  unsigned count_edges = 0;
  for (unsigned k = 0; k < graph.num_nodes(); ++k)
  {
    for (unsigned j = k + 1; j < graph.num_nodes(); ++j)
    {

      if (graph.has_edge(graph.node(k), graph.node(j)))
      {
        ++count_edges;
      }
    }
  }


 EXPECT_EQ(graph.num_edges(), count_edges);
}
