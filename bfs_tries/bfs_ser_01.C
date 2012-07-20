#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/point_traits.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>

using namespace boost;
using namespace graph;













struct finish_time_t {
  typedef vertex_property_tag kind;
};

template <typename TimeMap> class bfs_time_visitor:public default_bfs_visitor {
  public:
    bfs_time_visitor(TimeMap tmap):m_tmap(tmap) {}
    template <typename Vertex, typename Graph>
    void finish_vertex(Vertex u, const Graph & g) const {
      put(m_tmap, u, boost::posix_time::microsec_clock::local_time());
    }
    TimeMap m_tmap;
};








int main(int argc, char* argv[])
{
std::string file;
  std::string file_in;
  if(argc > 1){
    file = argv[1];
    file_in = file + ".gv";
  }else{
std::cout << std::endl;
    std::cout << "Possible command line arguments are as follows:" << std::endl;
    std::cout << "Filename (file of name Filename.gv assumed to be in the same folder)		required" << std::endl;
std::cout << std::endl;
    return 1;
  }

  std::cout << "graph is " << file_in << std::endl;

//typedef property<position_t, Point, property <vertex_name_t, std::string > > VertexProperty;
  typedef property<finish_time_t, boost::posix_time::ptime, property <vertex_name_t, std::string > > VertexProperty;
//  typedef property<vertex_name_t, std::string> VertexProperty;

  typedef adjacency_list<listS, vecS, undirectedS, VertexProperty> Graph;
  Graph g(0);
  dynamic_properties dp;
  typedef property_map<Graph, finish_time_t>::type FTimeMap;
  FTimeMap f_time = get(finish_time_t(), g);
  property_map<Graph, vertex_name_t>::type name = get(vertex_name, g);
  dp.property("node_id",name);

  std::ifstream infile(file_in.c_str());

  bool status = read_graphviz(infile, g, dp, "node_id");
  std::cout << "read was successfull: " << status << std::endl;

//  graph_traits<Graph>::edge_iterator e_i, e_end, next;

  graph_traits<Graph>::vertex_iterator v_i, v_end, v_next;
  //typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;
  //IndexMap index = get(boost::vertex_index, g);
  tie(v_i, v_end) = vertices(g);
  bfs_time_visitor<FTimeMap> vis(f_time);
  breadth_first_search(g, *v_i, visitor(vis));

  for(tie(v_i, v_end) = vertices(g); v_i != v_end; ++v_i){
    std::cout << "Vertex " << name[*v_i] << " finished by " << f_time[*v_i] << std::endl;
  }









/*
std::string file_out = file + "_new.gv";
 std::ofstream outfile(file_out.c_str());

 outfile << "strict graph my_rmat_graph {" << std::endl;
 outfile << "size=\"" << display_size << "!\"" << std::endl;
 outfile << "graph [bgcolor = \"black\", fontcolor=\"white\"]" << std::endl; // size=\"6,6\" ratio=fill]" << std::endl;
 outfile << "node [style=invis width=0, height=0 fixedsize=true label=\"\"]" << std::endl;


 double x,y;

 for(tie(v_i, v_end) = vertices(g); v_i != v_end; ++v_i){
  x = (position(*v_i))[0];
  y = (position(*v_i))[1];
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
*/


 return 0;
}
