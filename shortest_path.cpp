#include "shortest_path.hpp"
/**
 * @file shortest_path.cpp
 * Test script for using our templated Graph to determine shortest paths.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

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

  // HW1 #3: YOUR CODE HERE
  // Use nearest_node and shortest_path_lengths to set the node values
  // Construct a Color functor and view with the SFML_Viewer


  Point originpt = Point(-1, 0, 1);
  NodeIter node_it = nearest_node(graph, originpt);
  NodeType root_node=*node_it;

  int max_len= shortest_path_lengths(graph, root_node);

  struct ColorFn
  {
    int max_len;
    ColorFn(int max_len_):max_len{max_len_}{}
    CME212::Color operator()(const NodeType node)
    {
      const int node_value= node.value();
      //std::cout<< "node value" << node_value << std::endl;
      //std::cout<<"max_len" << max_len << std::endl;
     // std::cout << "divided" << node_value / max_len<< std::endl;
      
      double divided = 1.0- (double)node_value /  (double)max_len;
      //std::cout << "divided" << divided << std::endl;

      return CME212::Color::make_heat(divided);
    }
  };
  ColorFn colrfn(max_len);
  std::map<NodeType, unsigned> node_map= viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), ColorFn(max_len),node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  // Center the view and enter the event loop for interactivity
  viewer.center_view();
  viewer.event_loop();

  return 0;
}
