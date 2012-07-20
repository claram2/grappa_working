#include <iostream>
#include <fstream>
#include <iomanip>
#include "mpi.h"
#include <boost/graph/use_mpi.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/graph/distributed/graphviz.hpp>
#include "topology.hpp"
#include <boost/graph/random_layout.hpp>
#include <boost/graph/point_traits.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>

using namespace boost;
using namespace graph;
using namespace distributed;

struct position_t {
  typedef vertex_property_tag kind;
};

int main(int argc, char* argv[])
{
MPI_Init(&argc, &argv);
int display_size = 20;

  minstd_rand gen;

  typedef square_topology<> topology;
  typedef topology::point_type Point;
  typedef property<position_t, Point, property <vertex_name_t, std::string > > VertexProperty;

  typedef adjacency_list<listS, distributedS<mpi_process_group, vecS>, undirectedS, VertexProperty> Graph;
  Graph g(0);
  dynamic_properties dp;

  property_map<Graph, vertex_name_t>::type name = get(vertex_name, g);
  dp.property("node_id",name);

//  if(process_id(g.process_group()) == 0) {
  std::ifstream infile("ncvxqp9.gv");

  bool status = read_graphviz(infile, g, dp, "node_id");
  std::cout << "read was successfull: " << status << std::endl;
  //}

  graph_traits<Graph>::edge_iterator e_i, e_end, next;

  graph_traits<Graph>::vertex_iterator v_i, v_end, v_next;

  topology* my_h = new topology(gen);
  typedef property_map<Graph, position_t>::type PositionMap;
PositionMap position = get(position_t(), g);

  std::cout << "right before first synchronization" << std::endl;
  synchronize(g.process_group());
  std::cout << "right after first synchronization" << std::endl;
  random_graph_layout(g, position, *my_h);
  synchronize(g.process_group());
  std::cout << "randomly layed out" << std::endl;
  fruchterman_reingold_force_directed_layout(g, position, *my_h, cooling(linear_cooling<double>(10)));
  synchronize(g.process_group());
  std::cout << "force directed layed out" << std::endl;
  if(process_id(g.process_group()) == 0){
 std::ofstream outfile("ncvxqp9_new.gv");

 outfile << "strict graph my_rmat_graph {" << std::endl;
 outfile << "size=\"" << display_size << "!\"" << std::endl;
 outfile << "graph [bgcolor = \"black\", fontcolor=\"white\"]" << std::endl; // size=\"6,6\" ratio=fill]" << std::endl;
 outfile << "node [style=invis width=0, height=0 fixedsize=true label=\"\"]" << std::endl;


 double x,y;

 for(boost::tie(v_i, v_end) = vertices(g); v_i != v_end; ++v_i){
  Point p = position[*v_i];
  x = p[0];
  y = p[1];

  outfile << name[*v_i] << " [pos = \""; 
  outfile << (x*display_size) << ",";
  outfile << (y*display_size);
  outfile << "!\"]" << std::endl;
}

 std::cout << "vertices printed" << std::endl;

double distance;
 for(tie(e_i, e_end) = edges(g); e_i != e_end; ++e_i){
  distance = sqrt(((position[source(*e_i, g)])[0] - (position[target(*e_i, g)])[0])*((position[source(*e_i, g)])[0] - (position[target(*e_i, g)])[0]) + ((position[source(*e_i, g)])[1] - (position[target(*e_i, g)])[1])*((position[source(*e_i, g)])[1] - (position[target(*e_i, g)])[1]));
  outfile << name[source(*e_i, g)] << " -- " << name[target(*e_i, g)] << " [color = \"#";
  outfile << std::setfill ('0') << std::setw(2) << std::hex << (int)(((int)(distance*7*100))*0.001*255); 
  outfile << std::setfill ('0') << std::setw(2) << std::hex << (int)((((int)(distance*7*100000))% 1000)*0.001*255);
  outfile << std::setfill ('0') << std::setw(2) << std::hex << (int)((((int)(distance*7*100000000))% 1000)*0.001*255);
  
  outfile << "\"]" << std::dec << std::endl;
 }
 outfile << "}" << std::endl;

 outfile.close();
}

 MPI_Finalize();
 return 0;
}
