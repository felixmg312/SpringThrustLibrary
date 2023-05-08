/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include "CME212/SFML_Viewer.hpp"
#include "mass_spring.hpp"

int main(int argc, char **argv)
{
  // Check arguments
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct an empty graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  std::vector<typename GraphType::node_type> nodes;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int, 4> t;
  while (CME212::getline_parsed(tets_file, t))
  {
    graph.add_edge(nodes[t[0]], nodes[t[1]]);
    graph.add_edge(nodes[t[0]], nodes[t[2]]);
#if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it)
  {
    Node node = *it;
    node.value().vel = Point(0);
    node.value().mass = 1.0 / double(graph.size());
    node.value().initial_pos = node.position();
  }
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it)
  {
    Edge e = *it;
    e.value().K = 100;
    e.value().L = norm(e.node1().position() - e.node2().position());
  }

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // We want viewer interaction and the simulation at the same time
  // Viewer is thread-safe, so launch the simulation in a child thread
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&]()
                                {
    // Begin the mass-spring simulation
    double dt = 0.0005;
    double t_start = 0;
    double t_end = 1.0;
    

    for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt)
    {

      // viewer.clear();
      // node_map.clear();
      //Different test for my symp_euler_step ->> not sure if i need to remove it
      //symp_euler_step(graph, t, dt, Problem1Force());
      //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce()));
      //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce(0.01)), makeCombinedConstraints(SphereConstraint(), PlaneConstraint()));
      //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce(0.01)), SphereConstraint());

      symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce()),makeCombinedConstraints(SelfCollisionConstraint(),SphereConstraint()));
      // symp_euler_step(graph, t, dt, make_combined_force(DampingForce(0.01), GravityForce(), MassSpringForce()), makeCombinedConstraints(SelfCollisionConstraint(), PlaneConstraint(), PinConstraint()));
      // symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce()), PinConstraint());

      // Update viewer with nodes' new positions
      viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
      viewer.set_label(t);

      // These lines slow down the animation for small graphs, like grid0_*.
      // Feel free to remove them or tweak the constants.
      if (graph.size() < 100)
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    } }); // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}