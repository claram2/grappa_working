#include <iostream>
#include <fstream>
#include <vector>
#include <boost/graph/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/rmat_graph_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/graph/graphviz.hpp>
/*
namespace boost {
  enum vertex_position_t { vertex_position };
};*/

using namespace boost;

int main(int argc, char* argv[])
{
  int verts = 10;
  if(argc > 1){
    verts = atoi(argv[1]);
  }

  minstd_rand gen;

  typedef heart_topology<> topology;

  typedef topology::point_type my_p_type;
  //std::vector<my_p_type> my_positions(verts);
  

  heart_topology<> my_t();
  //topology my_t = new topology();


typedef adjacency_list<> Graph;
//typedef adjacency_list<vecS, vecS, directedS, property<vertex_position_t, my_p_type> > Graph;
typedef rmat_iterator<minstd_rand, Graph> RMATGen;

// Create graph with 100 nodes and 400 edges
  Graph g(RMATGen(gen, verts, 4*verts, 0.57, 0.19, 0.19, 0.05), RMATGen(), verts);

write_graphviz(std::cout, g);

/*






  random_graph_layout(g, get(vertex_position, g), my_t);

 typedef property_map<Graph, vertex_index_t>::type IndexMap;
 IndexMap index = get(vertex_index, g);

 std::ofstream outfile("rmat_graph.gv");



 outfile << "digraph my_rmat_graph {" << std::endl;
 graph_traits<Graph>::edge_iterator e_i, e_end;
 for(tie(e_i, e_end) = edges(g); e_i != e_end; ++e_i){
  outfile << index[source(*e_i, g)] << " -> " << index[target(*e_i, g)] << ";" << std::endl;
 }
 outfile << "}" << std::endl;

 outfile.close();

*/
  return 0;
}
