#include <iostream>
#include <fstream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/rmat_graph_generator.hpp>
#include <boost/random/linear_congruential.hpp>

typedef boost::adjacency_list<> Graph;
typedef boost::rmat_iterator<boost::minstd_rand, Graph> RMATGen;

int main(int argc, char* argv[])
{
  int verts = 10;
  if(argc > 1){
    verts = atoi(argv[1]);
  }
  boost::minstd_rand gen;
  // Create graph with 100 nodes and 400 edges
  Graph g(RMATGen(gen, verts, 4*verts, 0.57, 0.19, 0.19, 0.05), RMATGen(), verts);


 typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;
 IndexMap index = get(boost::vertex_index, g);

 std::ofstream outfile("rmat_graph.gv");


 outfile << "digraph my_rmat_graph {" << std::endl;
 boost::graph_traits<Graph>::edge_iterator e_i, e_end;
 for(tie(e_i, e_end) = edges(g); e_i != e_end; ++e_i){
  outfile << index[source(*e_i, g)] << " -> " << index[target(*e_i, g)] << ";" << std::endl;
 }
 outfile << "}" << std::endl;

 outfile.close();


  return 0;
}
