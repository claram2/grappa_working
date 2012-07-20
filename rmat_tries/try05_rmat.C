#include <iostream>
#include <fstream>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/rmat_graph_generator.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/random_layout.hpp>

using namespace boost;

int main(int argc, char* argv[])
{
  int verts = 10;
  if(argc > 1){
    verts = atoi(argv[1]);
  }
  minstd_rand gen;

typedef adjacency_list<> Graph;
  typedef rmat_iterator<minstd_rand, Graph> RMATGen;

  Graph g(RMATGen(gen, verts, 8*verts, 0.57, 0.19, 0.19, 0.05), RMATGen(), verts);

  write_graphviz(std::cout, g);  
 return 0;
}
