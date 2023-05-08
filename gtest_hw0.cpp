#include <gtest/gtest.h>
#include <fstream>
#include "CME212/Util.hpp"
#include "Graph.hpp"
#include<iostream>



class GraphPointFixture : public ::testing::Test {
 protected:
   //Define types
   using GraphType = Graph<int, int>;
   using NodeType = typename GraphType::node_type;
   using NodeIter = typename GraphType::node_iterator;
   using EdgeType = typename GraphType::edge_type;

   //Set up Graph and Points
   GraphType graph;
   std::vector<Point> points;
   virtual void SetUp()
   {
     for (int i = 0; i < 10; i++)
       points.push_back(Point(i));
  }
  
};

// Test has_node function
TEST_F(GraphPointFixture, HasNode){
  GraphType::node_type n0 = graph.add_node(points[0]);
  EXPECT_TRUE( graph.has_node(n0) ) << "has_node did not find n0";
}

// Test num nodes/size functions
TEST_F(GraphPointFixture, Size){

  EXPECT_EQ(graph.num_nodes(),graph.size()) << "num_nodes and size are different"  ;
  EXPECT_EQ(graph.size(), 0) << "starting size is not 0" ;

  graph.add_node(points[0]);
  graph.add_node(points[1]);

  EXPECT_EQ(graph.num_nodes(),graph.size()) << "num_nodes and size are different";
  EXPECT_EQ(graph.size(), 2) << "size is incorrect";
}

// Test edge function

TEST_F(GraphPointFixture, Edge){

  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  graph.add_node(points[2]);
  EdgeType e0 = graph.add_edge(n0, n1);
  
  EXPECT_EQ(e0, graph.edge(0)) << "error in edge retreval"  ;
}

// Verify only one of e0 < e1 or e1 < e0 is true
TEST_F(GraphPointFixture, Tricotomy){
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  
  EdgeType e0 = graph.add_edge(n0, n1);
  EdgeType e1 = graph.add_edge(n1, n2);

  EXPECT_TRUE( (e0 < e1) ^ (e1 < e0) ) << "error in edge comparison";
}

// Add existing edge doesn't create another edge
TEST_F(GraphPointFixture, AddEdge){
  graph.clear();
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);
  NodeType n4 = graph.add_node(points[4]);
  NodeType n5 = graph.add_node(points[5]);
  NodeType n6 = graph.add_node(points[6]);
  NodeType n7 = graph.add_node(points[7]);
  graph.add_edge(n0, n1);
  graph.add_edge(n2, n3);
  graph.add_edge(n4,n5);
  graph.add_edge(n6, n7);//edge3
  graph.add_edge(n7, n6);//edge4

  //graph.add_edge(n1, n0);
  //std::cout<< graph.num_edges() <<std::endl;
  EXPECT_EQ(graph.num_edges(), 4) << " adding size doesn't match" ;
  //EXPECT_TRUE(graph.add_edge(n6,n7).node1()==n6) << "(a,b)==(b,a)";
}
TEST_F(GraphPointFixture, Node_12)
{
  graph.clear();
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);
  NodeType n4 = graph.add_node(points[4]);
  NodeType n5 = graph.add_node(points[5]);
  NodeType n6 = graph.add_node(points[6]);
  NodeType n7 = graph.add_node(points[7]);
  graph.add_edge(n0, n1);
  graph.add_edge(n2, n3);
  graph.add_edge(n4, n5);
  graph.add_edge(n6, n7); //edge3
  graph.add_edge(n7, n6); //edge4

  //graph.add_edge(n1, n0);
 // std::cout << graph.num_edges() << std::endl;
  EXPECT_TRUE(n1==n1) << "(a,b)==(b,a)";
  //std::cout <<"index is "<< graph.add_edge(n6, n7).node1().index()<<std::endl;
  EXPECT_TRUE(graph.add_edge(n7, n6).node2()==n6)<< "(a,b)==(b,a)";
  EXPECT_TRUE(graph.add_edge(n1, n0).node1() == n1) << "(a,b)==(b,a)";
}
