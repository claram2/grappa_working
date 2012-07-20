#include <iostream>
#include <fstream>
#include <iomanip>
#include "mpi.h"
#include <boost/graph/use_mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/graph/rmat_graph_generator.hpp>
//#include <boost/graph/distributed/graphviz.hpp>
#include "topology.hpp"
#include <boost/graph/random_layout.hpp>
#include <boost/graph/point_traits.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
//#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/distributed/concepts.hpp>
#include <boost/property_map/parallel/global_index_map.hpp>

using namespace boost;
using namespace graph;
using namespace distributed;

struct position_t {
  typedef vertex_property_tag kind;
};
/*
struct global_position_t {
  typedef vertex_property_tag kind;
};
*/
int main(int argc, char* argv[]){
  mpi::environment env(argc, argv);
  mpi::communicator world;
//  MPI_Init(&argc, &argv);
  int display_size = 20;
  int verts = 40;
  int edg = 20*verts;
  minstd_rand gen;

  typedef square_topology<> topology;
  typedef topology::point_type Point;
//  typedef property<position_t, Point, property <vertex_name_t, std::string > > VertexProperty;
  typedef property<position_t, Point> VertexProperty;
  typedef adjacency_list<listS, distributedS<mpi_process_group, vecS>, undirectedS, VertexProperty> Graph;

  typedef rmat_iterator<minstd_rand, Graph> RMATGen;

  Graph g(RMATGen(gen, verts, edg, 0.15, 0.22, 0.30, 0.33), RMATGen(), verts);
  world.barrier();
  synchronize(g.process_group());  
  std::cout << "Hello from process " << process_id(g.process_group()) << std::endl;

  graph_traits<Graph>::edge_iterator e_i, e_end, next;
  tie(e_i, e_end) = edges(g);
  for(next = e_i; e_i != e_end; e_i = next){
    ++next;
    
    if(source(*e_i, g) == target(*e_i, g)){
      remove_edge(*e_i, g);
    }
  }
  graph_traits<Graph>::vertex_iterator v_i, v_end, v_next;



  graph_traits<Graph>::adjacency_iterator n_i, n_end;

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
  }

  typedef property_map<Graph, boost::vertex_index_t>::type IndexMap;
//  IndexMap index = get(boost::vertex_index, g);
  typedef property_map<Graph, vertex_global_t>::type GlobalMap;
  boost::parallel::global_index_map<IndexMap, GlobalMap> global_index(g.process_group(), num_vertices(g), get(vertex_index, g), get(vertex_global, g));

  topology* my_h = new topology(gen);
  typedef property_map<Graph, position_t>::type PositionMap;
  PositionMap position = get(position_t(), g);
//  typedef property_map<Graph, position_t*>::type PositionMap;
//  PositionMap position = get(&position_t, g);
  //typedef property_map<Graph, global_position_t>::type GlobalPositionMap;
//  typedef property_map<Graph, Point*>::type global_pos = get(&position_t(), g);
  
//  std::cout << "right before first synchronization" << std::endl;
  world.barrier();
  synchronize(g.process_group());
//  std::cout << "right after first synchronization" << std::endl;
  random_graph_layout(g, position, *my_h);
  world.barrier();
  synchronize(g.process_group());
  std::cout << "randomly layed out, process " << process_id(g.process_group()) << std::endl;

  Point origin;
  origin[0] = 0;
  origin[1] = 0;

  Point extent;
  extent[0] = 1;
  extent[1] = 1;






//  fruchterman_reingold_force_directed_layout(g, position, *my_h);//, origin, extent, attractive_force(square_distance_attractive_force()), repulsive_force(square_distance_repulsive_force()), force_pairs(all_force_pairs()), cooling(linear_cooling<double>(10)));
  world.barrier();
  synchronize(g.process_group());
//  std::cout << "force directed layed out, process " << process_id(g.process_group()) << std::endl;
//  property_map<Graph, Point position_t::*>::type 
  std::ofstream outfile("rmat_par.gv");
  Point p;
  Point q;
  double distance;
  if(process_id(g.process_group()) == 1){
//    std::ofstream outfile("rmat_par.gv");

    outfile << "strict graph my_rmat_graph {" << std::endl;
    outfile << "size=\"" << display_size << "!\"" << std::endl;
    outfile << "graph [bgcolor = \"black\", fontcolor=\"white\"]" << std::endl; // size=\"6,6\" ratio=fill]" << std::endl;
    outfile << "node [style=invis width=0, height=0 fixedsize=true label=\"\"]" << std::endl;


    double x,y;
  //  Point p;
    //Point q;
    for(boost::tie(v_i, v_end) = vertices(g); v_i != v_end; ++v_i){
      p = get(position, *v_i);//position[*v_i];
      x = p[0];
      y = p[1];

      outfile << get(global_index, *v_i) << " [pos = \""; 
      outfile << (x*display_size) << ",";
      outfile << (y*display_size);
      outfile << "!\"]" << std::endl;
    }

    std::cout << "vertices printed" << std::endl;

//    double distance;
    for(tie(e_i, e_end) = edges(g); e_i != e_end; ++e_i){
      request(position, target(*e_i, g));
    }
    std::cout << "Got through first loop" << std::endl;
  }
    world.barrier();
    synchronize(g.process_group());
  if(process_id(g.process_group()) == 1){
    for(tie(e_i, e_end) = edges(g); e_i != e_end; ++e_i){
//      std::cout << "Got into loop" << std::endl;
      p = get(position, source(*e_i, g));
  //    std::cout << "Got p" << std::endl;
    //  std::cout << "index of q is " << get(global_index, target(*e_i, g)) << std::endl;
      q = get(position, target(*e_i, g));
      //std::cout << "Got p and q" << std::endl;*/
      distance = sqrt((p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]));//sqrt(((position[source(*e_i, g)])[0] - (position[target(*e_i, g)])[0])*((position[source(*e_i, g)])[0] - (position[target(*e_i, g)])[0]) + ((position[source(*e_i, g)])[1] - (position[target(*e_i, g)])[1])*((position[source(*e_i, g)])[1] - (position[target(*e_i, g)])[1]));
      outfile << get(global_index, source(*e_i, g)) << " -- " << get(global_index, target(*e_i, g)) << " [color = \"#";
      outfile << std::setfill ('0') << std::setw(2) << std::hex << (int)(((int)(distance*7*100))*0.001*255); 
      outfile << std::setfill ('0') << std::setw(2) << std::hex << (int)((((int)(distance*7*100000))% 1000)*0.001*255);
      outfile << std::setfill ('0') << std::setw(2) << std::hex << (int)((((int)(distance*7*100000000))% 1000)*0.001*255);
  
      outfile << "\"]" << std::dec << std::endl;
    }
    outfile << "}" << std::endl;

    outfile.close();
  }
  std::cout << "Almost goodbye from process " << process_id(g.process_group()) << std::endl;
  synchronize(g.process_group());
  std::cout << "Goodbye from process " << process_id(g.process_group()) << std::endl;
 // MPI_Finalize();
  return 0;
}
