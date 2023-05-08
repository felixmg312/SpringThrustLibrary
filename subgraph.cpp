/**
 * @file subgraph.cpp
 * Test script for viewing a subgraph from our Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include "subgraph.hpp"

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Define our types
  using GraphType = Graph<int, int>;
  using NodeType = typename GraphType::node_type;
  using NodeIter = typename GraphType::node_iterator;

  // Construct a Graph
  GraphType graph;
  std::vector<NodeType> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SFML_Viewer
  CME212::SFML_Viewer viewer;

  // HW1 #4: YOUR CODE HERE
  // Use the filter_iterator to plot an induced subgraph.
  GraphType::NodeIterator it;
  SlicePredicate slp= SlicePredicate();
  auto node_map = viewer.empty_node_map(graph);
  filter_iterator<SlicePredicate, GraphType::NodeIterator> filt1 = make_filtered(graph.node_begin(), graph.node_end(), slp);
  filter_iterator<SlicePredicate, GraphType::NodeIterator> filt2 = make_filtered(graph.node_end(), graph.node_end(), slp);
  viewer.add_nodes(filt1,filt2,node_map);
  viewer.add_edges(graph.edge_begin(),graph.edge_end(),node_map);
  // Center the view and enter the event loop for interactivity

  //OddPredicate odp;
  //filter_iterator<OddPredicate, GraphType::NodeIterator> myfilt1 = make_filtered(graph.node_begin(), graph.node_end(), odp);
  //filter_iterator<OddPredicate, GraphType::NodeIterator> myfilt2 = make_filtered(graph.node_end(), graph.node_end(), odp);
  
  //GraphType::EdgeIterator its=graph.edge_begin();
  //for(;its!=graph.edge_end();++its){
   // std::cout<< (*its).node1().index() << std::endl;
  //}

  //auto node_map = viewer.empty_node_map(graph);
  //viewer.add_nodes(myfilt1, myfilt2, node_map);
  //viewer.add_edges(graph.edge_begin(),graph.edge_begin(),node_map);

  viewer.center_view();
  viewer.event_loop();

  return 0;
}
