#include <iostream>
#include <fstream>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/rmat_graph_generator.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/random_layout.hpp>

using namespace boost;

struct position_t {
  typedef vertex_property_tag kind;
};


int main(int argc, char* argv[])
{
  int verts = 10;
  if(argc > 1){
    verts = atoi(argv[1]);
  }
  minstd_rand gen;

  typedef property<position_t, circle_topology<>::point_type> PositionProperty;
  typedef adjacency_list<vecS, vecS, directedS, PositionProperty> Graph;
  typedef rmat_iterator<minstd_rand, Graph> RMATGen;

  Graph g(RMATGen(gen, verts, 8*verts, 0.57, 0.19, 0.19, 0.05), RMATGen(), verts);

  property_map<Graph, position_t>::type position = get(position_t(), g);

  circle_topology<>* my_h = new circle_topology<>(gen);

  random_graph_layout(g, position, *my_h);
 
  //dynamic_properties dp;
  //dp.property("pos", get(position, g));
  //write_graphviz_dp(std::cout, g, dp);
  write_graphviz(std::cout, g);  
 return 0;
}
