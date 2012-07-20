#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/rmat_graph_generator.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/point_traits.hpp>

using namespace boost;
using namespace graph;

struct position_t {
  typedef vertex_property_tag kind;
};


int main(int argc, char* argv[])
{

  int verts = 10;
  std::string file_out;
  int display_size = 20;
 
  if(argc > 1){
    file_out = argv[1];
  }else{
    std::cout << "Program requires a filename to print to as a command line argument" << std::endl;
    return 1;
  }
  if(argc > 2){
    verts = atoi(argv[2]);
  }else{
    std::cout << "Using default of 10 vertices because there was no command line argument" << std::endl;
  }
  int edg = 20*verts;
  if(argc > 3){
    edg = atoi(argv[3]);
  }else{
    std::cout << "Using default of 20 edges per vertex since there was no command line argument" << std::endl;
  }
  minstd_rand gen;

  typedef circle_topology<> topology;
  typedef topology::point_type Point;
  typedef property<position_t, Point> PositionProperty;
  typedef adjacency_list<listS, vecS, undirectedS, PositionProperty> Graph;
  typedef rmat_iterator<minstd_rand, Graph> RMATGen;

  Graph g(RMATGen(gen, verts, edg, 0.15, 0.22, 0.30, 0.33), RMATGen(), verts);



  graph_traits<Graph>::edge_iterator e_i, e_end, next;
  tie(e_i, e_end) = edges(g);
  for(next = e_i; e_i != e_end; e_i = next){
    ++next;
    
    if(source(*e_i, g) == target(*e_i, g)){
      remove_edge(*e_i, g);
    }
  }
  graph_traits<Graph>::vertex_iterator v_i, v_end, v_next;
  typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;
  IndexMap index = get(boost::vertex_index, g);
  graph_traits<Graph>::adjacency_iterator n_i, n_end;
/*
  int k = 0;
  tie(v_i, v_end) = vertices(g);
  while(v_i != v_end){
   tie(n_i, n_end) = adjacent_vertices(*v_i, g);
   if(n_i == n_end){
     remove_vertex(*v_i, g);
     tie(v_i, v_end) = vertices(g);
     v_i += k;
   }else{
     ++k;
     ++v_i;
   }
  }*/

 typedef property_map<Graph, position_t>::type PositionMap;
  PositionMap position = get(position_t(), g);

  topology* my_h = new topology(gen);

  random_graph_layout(g, position, *my_h);
  fruchterman_reingold_force_directed_layout(g, position, *my_h);

 std::ofstream outfile(file_out.c_str());

 outfile << "strict graph my_rmat_graph {" << std::endl;
 outfile << "size=\"" << display_size << "!\"" << std::endl;
 outfile << "graph [bgcolor = \"black\", fontcolor=\"white\"]" << std::endl; // size=\"6,6\" ratio=fill]" << std::endl;
 outfile << "node [style=invis width=0, height=0 fixedsize=true label=\"\"]" << std::endl;


 double x,y;

 for(tie(v_i, v_end) = vertices(g); v_i != v_end; ++v_i){
  x = (position(*v_i))[0];
  y = (position(*v_i))[1];
  outfile << index[*v_i] << " [pos = \""; 
  outfile << (x*display_size) << ",";
  outfile << (y*display_size);
  outfile << "!\"]" << std::endl;
 }

double distance;
 for(tie(e_i, e_end) = edges(g); e_i != e_end; ++e_i){
  distance = sqrt(((position[source(*e_i, g)])[0] - (position[target(*e_i, g)])[0])*((position[source(*e_i, g)])[0] - (position[target(*e_i, g)])[0]) + ((position[source(*e_i, g)])[1] - (position[target(*e_i, g)])[1])*((position[source(*e_i, g)])[1] - (position[target(*e_i, g)])[1]));
  outfile << index[source(*e_i, g)] << " -- " << index[target(*e_i, g)] << " [color = \"#";
  outfile << std::setfill ('0') << std::setw(2) << std::hex << (int)(((int)(distance*3.5*100))*0.001*255); 
  outfile << std::setfill ('0') << std::setw(2) << std::hex << (int)((((int)(distance*3.5*100000))% 1000)*0.001*255);
  outfile << std::setfill ('0') << std::setw(2) << std::hex << (int)((((int)(distance*3.5*100000000))% 1000)*0.001*255);
  
  outfile << "\"]" << std::dec << std::endl;
 }
 outfile << "}" << std::endl;

 outfile.close();



 return 0;
}
