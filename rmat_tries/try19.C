#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/rmat_graph_generator.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/gursoy_atun_layout.hpp>
#include <boost/graph/point_traits.hpp>
#include <boost/graph/graphviz.hpp>

using namespace boost;
using namespace graph;

struct position_t {
  typedef vertex_property_tag kind;
};

int main(int argc, char* argv[])
{

int display_size = 20;

  minstd_rand gen;

  typedef square_topology<> topology;
  typedef topology::point_type Point;
  typedef property<position_t, Point, property <vertex_name_t, std::string > > VertexProperty;

  typedef adjacency_list<listS, vecS, undirectedS, VertexProperty> Graph;
  Graph g(0);
  dynamic_properties dp;
  typedef property_map<Graph, position_t>::type PositionMap;
  PositionMap position = get(position_t(), g);

  property_map<Graph, vertex_name_t>::type name = get(vertex_name, g);
  dp.property("node_id",name);

  std::ifstream infile("caida.gv");

  bool status = read_graphviz(infile, g, dp, "node_id");
  std::cout << "read was successfull: " << status << std::endl;

  graph_traits<Graph>::edge_iterator e_i, e_end, next;

  graph_traits<Graph>::vertex_iterator v_i, v_end, v_next;
  typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;
  IndexMap index = get(boost::vertex_index, g);

  int k;


  topology* my_h = new topology(gen);

  random_graph_layout(g, position, *my_h);
  std::cout << "randomly layed out" << std::endl;
  fruchterman_reingold_force_directed_layout(g, position, *my_h);
  std::cout << "force directed layed out" << std::endl;
 std::ofstream outfile("caida_new.gv");

 outfile << "strict graph my_rmat_graph {" << std::endl;
 outfile << "size=\"" << display_size << "!\"" << std::endl;
 outfile << "graph [bgcolor = \"black\", fontcolor=\"white\"]" << std::endl; // size=\"6,6\" ratio=fill]" << std::endl;
 outfile << "node [style=invis width=0, height=0 fixedsize=true label=\"\"]" << std::endl;


 double x,y;
 k = 0;

 for(tie(v_i, v_end) = vertices(g); v_i != v_end; ++v_i){
  x = (position(*v_i))[0];
  y = (position(*v_i))[1];
  outfile << k << " [pos = \""; 
  outfile << (x*display_size/2) << ",";
  outfile << (y*display_size/2);
  outfile << "!\"]" << std::endl;
  ++k;
 }

 std::cout << "vertices printed" << std::endl;

double distance;
 for(tie(e_i, e_end) = edges(g); e_i != e_end; ++e_i){
  distance = sqrt(((position[source(*e_i, g)])[0] - (position[target(*e_i, g)])[0])*((position[source(*e_i, g)])[0] - (position[target(*e_i, g)])[0]) + ((position[source(*e_i, g)])[1] - (position[target(*e_i, g)])[1])*((position[source(*e_i, g)])[1] - (position[target(*e_i, g)])[1]));
  outfile << name[source(*e_i, g)] << " -- " << name[target(*e_i, g)] << " [color = \"#";
  outfile << std::setfill ('0') << std::setw(2) << std::hex << (int)(((int)(distance*16*100))*0.001*255); 
  outfile << std::setfill ('0') << std::setw(2) << std::hex << (int)((((int)(distance*16*100000))% 1000)*0.001*255);
  outfile << std::setfill ('0') << std::setw(2) << std::hex << (int)((((int)(distance*16*100000000))% 1000)*0.001*255);
  
  outfile << "\"]" << std::dec << std::endl;
 }
 outfile << "}" << std::endl;

 outfile.close();



 return 0;
}
