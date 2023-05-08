/**
 * @file poisson.cpp
 * Test script for treating for using the GraphSymmetricMatrix class
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SFML_Viewer to visualize the solution.
 */

#include "CME212/SFML_Viewer.hpp"
#include "GraphSymmetricMatrix.hpp"
#include <math.h>
// HW3: YOUR CODE HERE
// Define visual_iteration that inherits from cyclic_iteration

int main(int argc, char **argv)
{
  // Check arguments
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Define an empty Graph
  GraphType graph;

  {
    // Create a nodes_file from the first input argument
    std::ifstream nodes_file(argv[1]);
    // Interpret each line of the nodes_file as a 3D Point and add to the Graph
    std::vector<NodeType> node_vec;
    Point p;
    while (CME212::getline_parsed(nodes_file, p))
      node_vec.push_back(graph.add_node(2 * p - Point(1, 1, 0)));

    // Create a tets_file from the second input argument
    std::ifstream tets_file(argv[2]);
    // Interpret each line of the tets_file as four ints which refer to nodes
    std::array<int, 4> t;
    while (CME212::getline_parsed(tets_file, t))
    {
      graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
      graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
      graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
      graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
    }
  }

  // Get the edge length, should be the same for each edge
  auto it = graph.edge_begin();
  assert(it != graph.edge_end());
  double h = norm((*it).node1().position() - (*it).node2().position());

  // Make holes in our Graph
  remove_box(graph, Box3D(Point(-0.8 + h, -0.8 + h, -1), Point(-0.4 - h, -0.4 - h, 1)));
  remove_box(graph, Box3D(Point(0.4 + h, -0.8 + h, -1), Point(0.8 - h, -0.4 - h, 1)));
  remove_box(graph, Box3D(Point(-0.8 + h, 0.4 + h, -1), Point(-0.4 - h, 0.8 - h, 1)));
  remove_box(graph, Box3D(Point(0.4 + h, 0.4 + h, -1), Point(0.8 - h, 0.8 - h, 1)));
  remove_box(graph, Box3D(Point(-0.6 + h, -0.2 + h, -1), Point(0.6 - h, 0.2 - h, 1)));

  // HW3: YOUR CODE HERE
  // Define b using the graph, f, and g.
  // Construct the GraphSymmetricMatrix A using the graph
  // Solve Au = b using MTL.

  struct Boundary
  {
    double operator()(NodeType node)
    {
      CME212::BoundingBox<Point> bb = CME212::BoundingBox<Point>(Point(-0.6, -0.2, -1), Point(0.6, 0.2, 1));

      if (norm_inf(node.position()) == 1)
      {
        return 0;
      }
      else if (norm_inf(node.position() - Point(0.6, 0.6, 0)) < 0.2 || norm_inf(node.position() - Point(0.6, -0.6, 0)) < 0.2 || norm_inf(node.position() - Point(-0.6, 0.6, 0)) < 0.2 || norm_inf(node.position() - Point(-0.6, -0.6, 0)) < 0.2)
      {
        return -0.2;
      }
      else if (bb.contains(node.position()))
      {
        return 1;
      }
      else
      {
        return -100; // indicated that it is not in the boundary condition
      }
    }
  };

  struct Force
  {
    double operator()(NodeType node)
    {
      return 5 * cos(norm_1(node.position()));
    }
  };

  // set up isBoundary
  Boundary g;
  Force f;
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it)
  {
    NodeType node = *it;
    if (g(node) == -100)
    {
      node.value().isBoundary = false;
    }
    else
    {
      node.value().isBoundary = true;
    }
  }
  // set up b:
  mtl::dense_vector<double> b(graph.size());
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it)
  {
    NodeType node = *it;
    if (node.value().isBoundary)
    {
      b[node.index()] = g(node);
    }
    else
    {
      double temp = 0;
      for (auto it = node.edge_begin(); it != node.edge_end(); ++it)
      {
        if ((*it).node2().value().isBoundary)
        {
          temp += g((*it).node2());
        }
      }
      b[node.index()] += h * h * f(node) - temp;
    }
  }

  // Conjugate gradient solver
  const int size_ = graph.size();
  GraphSymmetricMatrix A = GraphSymmetricMatrix(&graph);

  //Initializing the viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();
  mtl::dense_vector<double> x(size_, 0.0);
  assert(x != b);
  itl::pc::identity<GraphSymmetricMatrix> P(A);

  // itl::cyclic_iteration<double> iter(b, 5000, 1.e-10, 1.e-10,50);
  // itl::cg(A, x, b, P, iter);

  //declaring my visual iterator
  visual_iteration<double, ColorFn, NodePosition> iter2(&graph, &viewer, &x, b, 5000, 1.e-10, 1.e-10, 1);


 // using thread to animate the cg solver
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&]()
  { itl::cg(A, x, b, P, iter2); });
  viewer.event_loop();
  interrupt_sim_thread = true;
  sim_thread.join();
}

