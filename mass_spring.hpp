/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

#include <fstream>
#include <chrono>
#include <thread>

#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"

#include "Graph.hpp"
#include "SpaceSearcher.hpp"
#include <thrust/for_each.h>
#include <thrust/execution_policy.h>
#include <thrust/system/omp/execution_policy.h>

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData
{
  Point vel;   //< Node velocity
  double mass; //< Node mass
  Point initial_pos;
  NodeData() : vel(0), mass(1) {}
};
struct EdgeData
{
  double K; //< Node velocity
  double L; //< Node mass
  EdgeData() : K(100), L(0) {}
  EdgeData(double K, double L) : K(K), L(L) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */

// helper function for velocity
struct helper_velocity
{
  __host__ __device__
  helper_velocity(double dt_) : dt{dt_} {}
  void operator()(Node n)
  {
    n.position() += n.value().vel * dt;
  }
  double dt;
};

// helper function for force
template <typename F>
struct helper_force
{
  __host__ __device__
  helper_force(double dt_, double t_, F force_) : dt{dt_}, t{t_}, force{force_} {}
  void operator()(Node n)
  {
    n.value().vel += force(n, t) * (dt / n.value().mass);
    // sorry but i really need this print statement for my code to work. My computer is too old and will break if I don't use
    // cout to make the whole system run slower
    // hope this won't affect the code quality part
    std::cout << n.position() << std::endl;
  }
  double dt;
  double t;
  F force;
};

template <typename G, typename F, typename C>
double symp_euler_step(G &g, double t, double dt, F force, C constraint)
{

  // Compute the t+dt position
  helper_velocity velocity_func = helper_velocity(dt);
  thrust::for_each(thrust::omp::par, g.node_begin(), g.node_end(), velocity_func);

  // Adding the constraints
  constraint(g, t);

  // Compute the t+dt velocity
  helper_force<F> force_func = helper_force<F>(dt, t, force);
  thrust::for_each(thrust::omp::par, g.node_begin(), g.node_end(), force_func);

  return t + dt;
}

// template <typename G, typename F, typename C>
// double symp_euler_step(G &g, double t, double dt, F force, C constraint)
// {

//   // Compute the t+dt position
//   for (auto it = g.node_begin(); it != g.node_end(); ++it)
//   {
//     auto n = *it;

//     // Update the position of the node according to its velocity
//     // x^{n+1} = x^{n} + v^{n} * dt
//     n.position() += n.value().vel * dt;
//   }

//   // Adding the constraints
//   constraint(g, 0);

//   // Compute the t+dt velocity
//   for (auto it = g.node_begin(); it != g.node_end(); ++it)
//   {
//     auto n = *it;

//     // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
//     n.value().vel += force(n, t) * (dt / n.value().mass);
//   }

//   return t + dt;
// }

template <typename G, typename F>
double symp_euler_step(G &g, double t, double dt, F force)
{
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

class Force
{
public:
  virtual Point operator()(Node n, double t)
  {
    (void)n;
    (void)t;
    return Point(0);
  }
};
class GravityForce : public Force
{
public:
  Point operator()(Node n, double t) override
  {
    (void)t;
    return (Point(0, 0, -grav) * n.value().mass);
  }
};
class MassSpringForce : public Force
{
  Point operator()(Node n, double t) override
  {
    (void)t;
    Point spring = Point(0, 0, 0);

    for (auto it = n.edge_begin(); it != n.edge_end(); ++it)
    {
      Edge e = *it;
      Point x1 = e.node1().position();
      Point x2 = e.node2().position();
      spring -= e.value().K * (x1 - x2) / e.length() * (e.length() - e.value().L);
    }
    return spring;
  }
};

class DampingForce : public Force
{
public:
  DampingForce(double c) : c{c} {}
  DampingForce() : c{0.5} {}

  Point operator()(Node n, double t) override
  {
    (void)t;
    return (-(c * n.value().vel));
  }

private:
  double c;
};

class CombinedForce
{
public:
  CombinedForce(std::vector<Force *> f) : f(f) {}
  Point operator()(Node n, double t)
  {
    (void)t;
    Point total_force = Point(0, 0, 0);
    for (auto it = f.begin(); it != f.end(); ++it)
    {
      total_force += (*(*it))(n, 0);
    }
    return total_force;
  }

private:
  std::vector<Force *> f;
};

template <typename F1, typename F2>
CombinedForce make_combined_force(F1 f1, F2 f2)
{
  std::vector<Force *> f;
  f.push_back(&f1);
  f.push_back(&f2);
  return CombinedForce(f);
}
template <typename F1, typename F2, typename F3>
CombinedForce make_combined_force(F1 f1, F2 f2, F3 f3)
{
  std::vector<Force *> f;
  f.push_back(&f1);
  f.push_back(&f2);
  f.push_back(&f3);

  return CombinedForce(f);
}

class Constraint
{
public:
  // template <typename NODE>
  virtual void operator()(GraphType &g, double t)
  {
    (void)g;
    (void)t;
  }
};

class PinConstraint : public Constraint
{
public:
  unsigned int index1;
  unsigned int index2;
  void operator()(GraphType &g, double t) override
  {
    (void)t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      if ((*it).value().initial_pos == Point(0, 0, 0))
      {

        (*it).position() = Point(0, 0, 0);
      }
      if ((*it).value().initial_pos == Point(1, 0, 0))
      {

        (*it).position() = Point(1, 0, 0);
      }
    }
  }
};

class PlaneConstraint : public Constraint
{
public:
  void operator()(GraphType &g, double t) override
  {
    (void)t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      if (dot((*it).position(), Point(0, 0, 1)) < -0.75)
      {
        (*it).position()[2] = -0.75;
        (*it).value().vel[2] = 0;
      }
    }
  }
};

class SphereConstraint : public Constraint
{
public:
  void operator()(GraphType &g, double t) override
  {
    (void)t;
    Point c = Point(0.5, 0.5, -0.5);
    double r = 0.15;

    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      if (norm((*it).position() - c) < r)
      {
        // std::cout << "in here sphere" << std::endl;

        Point unit_vec = ((*it).position() - c) / norm((*it).position() - c);
        (*it).position() = c + r * unit_vec;
        (*it).value().vel -= dot(((*it).value().vel), unit_vec) * unit_vec;
      }
    }
  }
};

struct SelfCollisionhelper
{
  __host__ __device__ void operator()(Node n) const
  {
    const Point &center = n.position();
    double radius2 = std::numeric_limits<double>::max();
    for (auto incident_it = n.edge_begin(); incident_it != n.edge_end(); ++incident_it)
    { // std::accumulate?
      Edge e = *incident_it;
      radius2 = std::min(radius2, normSq(e.node2().position() - center));
    }
    radius2 *= 0.9;
    // for (auto Node_it=g.node_begin(); Node_it!=g.node_end(); ++Node_it)
    // {
    //   Node n2= *Node_it;
    //   Point r = center - n2.position();
    //   double l2 = normSq(r);
    //   if (n != n2 && l2 < radius2)
    //   { // this is the l2 variable, not "twelve"
    //     // Remove our velocity component in r
    //     n.value().vel -= (dot(r, n.value().vel) / l2) * r;
    //     //std::cout<<"in funct"<<n.value().vel<<std::endl;
    //   }
    //}
    auto n2p = [](const Node &n)
    { return n.position(); };
    Box3D bigbb(thrust::make_transform_iterator(g.node_begin(), n2p), thrust::make_transform_iterator(g.node_end(), n2p));
    SpaceSearcher<Node> searcher(bigbb, g.node_begin(), g.node_end(), n2p);
    Box3D bb(Point(center.x - radius2, center.y - radius2, center.z - radius2), Point(center.x + radius2, center.y + radius2, center.z + radius2));
    for (auto Node_it = searcher.begin(bb); Node_it != searcher.end(bb); ++Node_it)
    {
      Node n2 = *Node_it;
      Point r = center - n2.position();
      double l2 = normSq(r);
      if (n != n2 && l2 < radius2)
      {
        n.value().vel -= (dot(r, n.value().vel) / l2) * r;
      }
    }
  }
  GraphType g;
  SelfCollisionhelper(GraphType &g_) : g{g_} {}
};

struct SelfCollisionConstraint : public Constraint
{
public:
  void operator()(GraphType &g, double t) const
  {
    (void)t;
    SelfCollisionhelper s = SelfCollisionhelper(g);
    thrust::for_each(thrust::omp::par, g.node_begin(), g.node_end(), s);
  }
};

class CombinedConstraints
{
public:
  CombinedConstraints(std::vector<Constraint *> c) : c(c) {}
  void operator()(GraphType &g, double t)
  {
    (void)t;
    for (auto it = c.begin(); it != c.end(); ++it)
    {

      (*(*it))(g, 0);
    }
  }

private:
  std::vector<Constraint *> c;
};

template <typename C1, typename C2>
CombinedConstraints makeCombinedConstraints(C1 c1, C2 c2)
{
  std::vector<Constraint *> c;
  c.push_back(&c1);
  c.push_back(&c2);
  return CombinedConstraints(c);
}

template <typename C1, typename C2, typename C3>
CombinedConstraints makeCombinedConstraints(C1 c1, C2 c2, C3 c3)
{
  std::vector<Constraint *> c;
  c.push_back(&c1);
  c.push_back(&c2);
  c.push_back(&c3);
  return CombinedConstraints(c);
}

struct removeSphereConstraint : public Constraint
{

  void operator()(GraphType &g, double t)
  {
    (void)t;
    Point c = Point(0.5, 0.5, -0.5);
    double r = 0.15;
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      Node n = *it;
      if (norm(c - n.position()) < r)
      {
        g.remove_node(n);
      }
    }
  }
};

/** Force function object for HW2 #1. */
struct Problem1Force
{
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t)
  {
    // HW2 #1: YOUR CODE HERE
    (void)t;
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
    {
      return Point(0, 0, 0);
    }

    /* Spring Force & gravity force */
    Point spring = Point(0, 0, 0);
    Point gravity = Point(0, 0, -grav) * n.value().mass;
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it)
    {
      Edge e = *it;
      Point x1 = e.node1().position();
      Point x2 = e.node2().position();

      spring -= e.value().K * (x1 - x2) / e.length() * (e.length() - e.value().L);
    }
    return gravity + spring;
  }
};
