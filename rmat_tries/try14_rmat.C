#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/rmat_graph_generator.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/gursoy_atun_layout.hpp>
#include <boost/graph/point_traits.hpp>

using namespace boost;
using namespace graph;

struct position_t {
  typedef vertex_property_tag kind;
};


int main(int argc, char* argv[])
{
  int verts = 10;
  //Dim width(2);
  //Dim height(2);
  //float width = 2.0;
  //float height = 2.0;
  if(argc > 1){
    verts = atoi(argv[1]);
  }
  minstd_rand gen;

  typedef square_topology<> topology;
  typedef topology::point_type Point;
  typedef property<position_t, Point> PositionProperty;
  typedef adjacency_list<vecS, vecS, undirectedS, PositionProperty> Graph;
  typedef rmat_iterator<minstd_rand, Graph> RMATGen;

  Graph g(RMATGen(gen, verts, 8*verts, 0.15, 0.22, 0.30, 0.33), RMATGen(), verts);

  typedef property_map<Graph, position_t>::type PositionMap;
  PositionMap position = get(position_t(), g);

  topology* my_h = new topology(gen);

  typedef property_traits<PositionMap>::value_type Point2;
  //typedef point_traits<Point2>::component_type Dim;
//  point_traits<Point2>::component_type* width = new point_traits<Point2>::component_type(2);
  //Dim* width = new Dim(2);  

  random_graph_layout(g, position, *my_h);
  //fruchterman_reingold_force_directed_layout(g, position);//, 2.0, 2.0, cooling(linear_cooling<double>(100)));//, *width, height);
  gursoy_atun_layout(g, *my_h, position);
  //dynamic_properties dp;
  //dp.property("pos", get(position, g));
  //write_graphviz_dp(std::cout, g, dp);

 std::ofstream outfile("rmat_graph.gv");


 outfile << "graph my_rmat_graph {" << std::endl;
 outfile << "size=\"6!\"" << std::endl;
 outfile << "graph [bgcolor = \"black\", fontcolor=\"white\"]" << std::endl; // size=\"6,6\" ratio=fill]" << std::endl;
 outfile << "node [style=invis width=0, height=0 fixedsize=true label=\"\"]" << std::endl;
 graph_traits<Graph>::edge_iterator e_i, e_end;
 
 typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;
 IndexMap index = get(boost::vertex_index, g);

 for(int i = 0; i < verts; ++i){
  outfile << index[i] << " [pos = \"" << (position[i])[0]<< "," << (position[i])[1] << "\"]" << std::endl;
 }

double distance;
 for(tie(e_i, e_end) = edges(g); e_i != e_end; ++e_i){
  distance = sqrt(((position[source(*e_i, g)])[0] - (position[target(*e_i, g)])[0])*((position[source(*e_i, g)])[0] - (position[target(*e_i, g)])[0]) + ((position[source(*e_i, g)])[1] - (position[target(*e_i, g)])[1])*((position[source(*e_i, g)])[1] - (position[target(*e_i, g)])[1]));
  outfile << index[source(*e_i, g)] << " -- " << index[target(*e_i, g)] << " [color = \"#";
  outfile << std::setfill ('0') << std::setw(2) << std::hex << (int)(((int)(distance*16*100))*0.001*255); 
  outfile << std::setfill ('0') << std::setw(2) << std::hex << (int)((((int)(distance*16*100000))% 1000)*0.001*255);
  outfile << std::setfill ('0') << std::setw(2) << std::hex << (int)((((int)(distance*16*100000000))% 1000)*0.001*255);
  
  outfile << "\"]" << std::dec << std::endl;
 }
 outfile << "}" << std::endl;

 outfile.close();



 return 0;
}
