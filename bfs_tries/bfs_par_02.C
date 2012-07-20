#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <vector>
#include "mpi.h"
#include <boost/graph/use_mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/graph/distributed/concepts.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
//#include "posix_time.hpp"
#include <boost/random/linear_congruential.hpp>
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/distributed/graphviz.hpp>
#include <boost/graph/point_traits.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/graph/distributed/concepts.hpp>
#include <boost/property_map/parallel/global_index_map.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

//#include <boost/serialization/serialization.hpp>
//#include "posix_time.hpp"

using namespace boost;
using namespace graph;
using namespace distributed;

/*
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
*//*
struct vertex_proc_num_t {
  typedef vertex_property_tag kind;
};
*/






int main(int argc, char* argv[])
{

  mpi::environment env(argc, argv);
  mpi::communicator world;
  int my_proc_id = process_id(world);

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
  std::ifstream infile(file_in.c_str());
  char_separator<char> sep;
  std::string sep_also;
  int verts = 0;
  int lines = 0;
  std::string temp;
  infile >> temp;
  typedef adjacency_list<listS, distributedS<mpi_process_group, vecS>, undirectedS> Graph;
  if(temp.c_str() == "graph"){
//    typedef adjacency_list<listS, distributedS<mpi_process_group, vecS>, undirectedS> Graph;
    sep = char_separator<char>("--");
    sep_also = "--";
  }else{
//    typedef adjacency_list<listS, distributedS<mpi_process_group, vecS>, directedS> Graph;
    sep = char_separator<char>("->");
    sep_also = "->";
  }
  getline(infile, temp);
  //std::cout << "Remains of first line are " << temp << std::endl;
  while(getline(infile, temp)){
    ++lines;
    tokenizer<char_separator<char> > tokens(temp, sep);
    BOOST_FOREACH(const std::string& t, tokens){
      verts = std::max(verts, atoi(t.c_str()));
    }
  }

  verts += 1;
  std::cout << "Number of vertices = " << verts << " from process " << my_proc_id << std::endl;
  std::cout << "Number of lines = " << lines << " from process " << my_proc_id << std::endl;


/*
//typedef property<position_t, Point, property <vertex_name_t, std::string > > VertexProperty;
//  typedef property<finish_time_t, boost::posix_time::ptime, property <vertex_name_t, std::string, property<vertex_proc_num_t, int> > > VertexProperty;
//  typedef property<vertex_name_t, std::string> VertexProperty;
  typedef property<vertex_name_t, std::string, property<vertex_proc_num_t, int> > VertexProperty;

typedef adjacency_list<listS, distributedS<mpi_process_group, vecS>, undirectedS, VertexProperty> Graph;
  Graph g(0);
  dynamic_properties dp;
property_map<Graph, vertex_name_t>::type name = get(vertex_name, g);
  dp.property("node_id",name);
  if(my_proc_id == 0){
  std::ifstream infile(file_in.c_str());

  bool status = read_graphviz(infile, g, dp, "node_id");
  std::cout << "read was successfull: " << status << std::endl;
  }*/

  Graph g(verts);
  std::ifstream infile2(file_in.c_str());
  getline(infile2, temp);
//  std::vector<std::string> this_edge;
  int source, target;
  if(my_proc_id == 0){
  for(int i = 1; i < lines; ++i){
    getline(infile2, temp);
//    split(this_edge, temp, is_any_of(sep_also.c_str()));
//    if(temp.c_str() != "}\n"){
   // std::cout << "Line is " << temp << std::endl;
    tokenizer<char_separator<char> > tokens(temp, sep);
    tokenizer<char_separator<char> >::iterator beg=tokens.begin();
    source = atoi((*beg).c_str());
   // if(beg != tokens.end()){
    ++beg;
    target = atoi((*beg).c_str());
//    std::cout << "Adding edge " << source << "," << target << std::endl;// << " or " << this_edge[0] << "," << this_edge[1] << std::endl;
    add_edge(vertex(source, g), vertex(target, g), g);
    //add_edge(vertex(atoi((this_edge[0]).c_str()), g), vertex((this_edge[1]).c_str()), g), g);
   // }
  }
  std::cout << "Added all edges" << std::endl;
  }

  world.barrier();
  synchronize(g.process_group());

 // graph_traits<Graph>::vertex_descriptor try_1 = boost::vertex(0,g);
//  graph_traits<Graph>::vertex_descriptor try = vertex(0,g);
  //std::cout << "Worked before the redistribute" << std::endl;
/*  property_map<Graph, vertex_proc_num_t>::type proc_map = get(vertex_proc_num_t(), g);
  graph_traits<Graph>::vertex_iterator v_i, v_end;
  for(tie(v_i, v_end) = vertices(g); v_i != v_end; ++v_i){
    std::string vert_name = name[*v_i];
    proc_map[*v_i] = atoi(vert_name.c_str()) % num_processes(g.process_group());
  }
  g.redistribute(proc_map);*/
/*
  typedef property_map<Graph, finish_time_t>::type FTimeMap;
  FTimeMap f_time = get(finish_time_t(), g);
*/
/*
  graph_traits<Graph>::edge_iterator e_i, e_end, next;

//  graph_traits<Graph>::vertex_iterator v_i, v_end, v_next;
  //typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;
  //IndexMap index = get(boost::vertex_index, g);
//  tie(v_i, v_end) = vertices(g);
  //bfs_time_visitor<FTimeMap> vis(f_time);


  typedef property_map<Graph, boost::vertex_index_t>::type IndexMap;
  IndexMap index = get(boost::vertex_index, g);
  typedef property_map<Graph, vertex_global_t>::type GlobalMap;
  boost::parallel::global_index_map<IndexMap, GlobalMap> global_index(g.process_group(), num_vertices(g), get(vertex_index, g), get(vertex_global, g));




  world.barrier();
  synchronize(g.process_group());
  world.barrier();




  std::ostringstream local_out;
  std::string file_out = file + "_new.gv";
  std::ofstream outfile(file_out.c_str());
  int k;



  for(tie(e_i, e_end) = edges(g); e_i != e_end; ++e_i){
//    p = get(position, source(*e_i, g));
  //  q = get(position, target(*e_i, g));
    //distance = sqrt((p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]));//sqrt(((position[source(*e_i, g)])[0] - (position[target(*e_i, g)])[0])*((position[source(*e_i, g)])[0] - (position[target(*e_i, g)])[0]) + ((position[source(*e_i, g)])[1] - (position[target(*e_i, g)])[1])*((position[source(*e_i, g)])[1] - (position[target(*e_i, g)])[1]));
    //if( distance != 0.0){
      local_out << get(global_index, source(*e_i, g)) << " -- " << get(global_index, target(*e_i, g)) << std::endl;//" [color = \"#";
  */  /*  local_out << std::setfill ('0') << std::setw(2) << std::hex << (int)(((int)(distance*3.5*100))*0.001*255);
      local_out << std::setfill ('0') << std::setw(2) << std::hex << (int)((((int)(distance*3.5*100000))% 1000)*0.001*255);
      local_out << std::setfill ('0') << std::setw(2) << std::hex << (int)((((int)(distance*3.5*100000000))% 1000)*0.001*255);
      local_out << "\"]" << std::dec << std::endl;
    }*//*
  }

  if(my_proc_id == 0){
    outfile << local_out.str();
    synchronize(g.process_group());
    for(k = 1; k < num_processes(g.process_group()); ++k){
      int len;
      receive(g.process_group(), k, 0, len);
      char* data = new char [len+1];
      data[len] = 0;
      receive(g.process_group(), k, 1, data, len);
      outfile << data;
      delete [] data;
    }
    outfile << "}" << std::endl;
  }else{
    std::string local_result = local_out.str();
    const char* data = local_result.c_str();
    int len = local_result.length();
    send(g.process_group(), 0, 0, len);
    send(g.process_group(), 0, 1, data, len);
    synchronize(g.process_group());
  }








*/
  world.barrier();
  synchronize(g.process_group());
 // world.barrier();


  std::cout << "About to search" << std::endl;
  graph_traits<Graph>::vertex_descriptor start = boost::vertex(0, g);
  std::cout << "Gotten past making start" << std::endl;
//  std::cout << "Pointer is " << start << std::endl;
  breadth_first_search(g, start, visitor(bfs_visitor<null_visitor>()));
/*
  for(tie(v_i, v_end) = vertices(g); v_i != v_end; ++v_i){
    std::cout << "Vertex " << name[*v_i] << " finished by " << f_time[*v_i] << std::endl;
  }

*/
  std::cout << "Breadth first search completed" << std::endl;






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
