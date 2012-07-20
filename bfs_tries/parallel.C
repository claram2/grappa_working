#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include "mpi.h"
#include <boost/graph/use_mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/graph/distributed/concepts.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include "topology.hpp"
#include <boost/graph/distributed/graphviz.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/graph/point_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/graph/distributed/concepts.hpp>
#include <boost/property_map/parallel/global_index_map.hpp>

using namespace boost;
using namespace graph;
using namespace distributed;

struct position_t {
  typedef vertex_property_tag kind;
};

struct vertex_proc_num_t {
  typedef vertex_property_tag kind;
};

int main(int argc, char* argv[]){
  mpi::environment env(argc, argv);
  mpi::communicator world;
  int my_proc_id = process_id(world);
  int iter = 10;
  std::string file;
  std::string file_in;
  if(argc > 1){
    file = argv[1];
    file_in = file + ".gv";
  }else{
    if(my_proc_id == 0){
      std::cout << std::endl;
      std::cout << "Possible command line arguments are as follows:" << std::endl;
      std::cout << "Filename (file of name Filename.gv assumed to be in the same folder)		required" << std::endl;
      std::cout << "Number of inches that the .gv file should be scaled to when it is processed	optional	default of 20" << std::endl;
      std::cout << std::endl;
    }
    return 1;
  }
  int display_size = 20;
  if(argc > 2){
    display_size = atoi(argv[2]);
    if(display_size < 1){
      if(my_proc_id == 0){
        std::cout << "Display size must be positive" << std::endl;
      }
    }
  }else{
    if(my_proc_id == 0){
      std::cout << "Program can take an integer inch display size to scale to, using default of 20 inches" << std::endl;
    }
  }
  if(my_proc_id == 0){
    std::cout << "Reading file " << file_in << std::endl;
  }
  world.barrier();

  typedef circle_topology<> topology;
  typedef topology::point_type Point;
  typedef property<position_t, Point, property<vertex_name_t, std::string, property< vertex_proc_num_t, int> > > VertexProperty;
  
  typedef adjacency_list<listS, distributedS<mpi_process_group, vecS>, undirectedS, VertexProperty> Graph;

  Graph g(0);
  dynamic_properties dp;

  property_map<Graph, vertex_name_t>::type name = get(vertex_name, g);
  dp.property("node_id",name);
  
  if(my_proc_id == 0){
    std::ifstream infile(file_in.c_str());
    read_graphviz(infile, g, dp, "node_id");
  }

  world.barrier();
  synchronize(g.process_group());  

  typedef property_map<Graph, boost::vertex_index_t>::type IndexMap;
  IndexMap index_1 = get(boost::vertex_index, g);
  typedef property_map<Graph, vertex_global_t>::type GlobalMap;
  boost::parallel::global_index_map<IndexMap, GlobalMap> global_index_1(g.process_group(), num_vertices(g), get(vertex_index, g), get(vertex_global, g));

  property_map<Graph, vertex_proc_num_t>::type proc_map = get(vertex_proc_num_t(), g);
  graph_traits<Graph>::vertex_iterator v_i, v_end;
  for(tie(v_i, v_end) = vertices(g); v_i != v_end; ++v_i){
    proc_map[*v_i] = get(global_index_1, *v_i) % num_processes(g.process_group());
  }
  g.redistribute(proc_map);
  graph_traits<Graph>::edge_iterator e_i, e_end, next;

  typedef property_map<Graph, position_t>::type PositionMap;
  PositionMap position = get(position_t(), g);

  IndexMap index = get(boost::vertex_index, g);
boost::parallel::global_index_map<IndexMap, GlobalMap> global_index(g.process_group(), num_vertices(g), get(vertex_index, g), get(vertex_global, g));
  minstd_rand gen2((unsigned int) (my_proc_id+1));
  topology* my_h = new topology(gen2);

  world.barrier();
  synchronize(g.process_group());
  random_graph_layout(g, position, *my_h);
  world.barrier();
  synchronize(g.process_group());
  
  std::ostringstream local_out;
  std::string file_out = file + "_new.gv";
  std::ofstream outfile(file_out.c_str());
  Point p;
  Point q;
  int k;
  double distance;
  synchronize(g.process_group());
  double x,y;

  for(boost::tie(v_i, v_end) = vertices(g); v_i != v_end; ++v_i){
    p = get(position, *v_i);//position[*v_i];
    x = p[0];
    y = p[1];
    local_out << get(global_index, *v_i) << " [pos = \""; 
    local_out << (x*display_size) << ",";
    local_out << (y*display_size);
    local_out << "!\"]" << std::endl;
  }

  if(my_proc_id == 0){
    outfile << "strict graph my_rmat_graph {" << std::endl;
    outfile << "size=\"" << display_size << "!\"" << std::endl;
    outfile << "graph [bgcolor = \"black\", fontcolor=\"white\"]" << std::endl; // size=\"6,6\" ratio=fill]" << std::endl;
    outfile << "node [style=invis width=0, height=0 fixedsize=true label=\"\"]" << std::endl;
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
  }else{
    std::string local_result = local_out.str();
    const char* data = local_result.c_str();
    int len = local_result.length();
    send(g.process_group(), 0, 0, len);
    send(g.process_group(), 0, 1, data, len);
    synchronize(g.process_group());
  }

  synchronize(g.process_group());
  synchronize(g.process_group());
  synchronize(g.process_group());
  local_out.str("");

  for(tie(e_i, e_end) = edges(g); e_i != e_end; ++e_i){
    request(position, target(*e_i, g));
  }
  world.barrier();
  synchronize(g.process_group());
  synchronize(g.process_group());

  for(tie(e_i, e_end) = edges(g); e_i != e_end; ++e_i){
    p = get(position, source(*e_i, g));
    q = get(position, target(*e_i, g));
    distance = sqrt((p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]));//sqrt(((position[source(*e_i, g)])[0] - (position[target(*e_i, g)])[0])*((position[source(*e_i, g)])[0] - (position[target(*e_i, g)])[0]) + ((position[source(*e_i, g)])[1] - (position[target(*e_i, g)])[1])*((position[source(*e_i, g)])[1] - (position[target(*e_i, g)])[1]));
    if( distance != 0.0){
      local_out << get(global_index, source(*e_i, g)) << " -- " << get(global_index, target(*e_i, g)) << " [color = \"#";
      local_out << std::setfill ('0') << std::setw(2) << std::hex << (int)(((int)(distance*3.5*100))*0.001*255); 
      local_out << std::setfill ('0') << std::setw(2) << std::hex << (int)((((int)(distance*3.5*100000))% 1000)*0.001*255);
      local_out << std::setfill ('0') << std::setw(2) << std::hex << (int)((((int)(distance*3.5*100000000))% 1000)*0.001*255);
      local_out << "\"]" << std::dec << std::endl;
    }
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

  return 0;
}

