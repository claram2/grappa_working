#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include "mpi.h"
#include <boost/graph/use_mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/graph/distributed/concepts.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/random/linear_congruential.hpp>
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

using namespace boost;
using namespace graph;
using namespace distributed;

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
    sep = char_separator<char>("--");
  }else{
    sep = char_separator<char>("->");
  }
  getline(infile, temp);
  while(getline(infile, temp)){
    ++lines;
    tokenizer<char_separator<char> > tokens(temp, sep);
    BOOST_FOREACH(const std::string& t, tokens){
      verts = std::max(verts, atoi(t.c_str()));
    }
  }

  verts += 1;

  Graph g(verts);
  std::ifstream infile2(file_in.c_str());
  getline(infile2, temp);
  int source, target;
  if(my_proc_id == 0){
    for(int i = 1; i < lines; ++i){
      getline(infile2, temp);
      tokenizer<char_separator<char> > tokens(temp, sep);
      tokenizer<char_separator<char> >::iterator beg=tokens.begin();
      source = atoi((*beg).c_str());
      ++beg;
      target = atoi((*beg).c_str());
      add_edge(vertex(source, g), vertex(target, g), g);
    }
  }

  synchronize(g.process_group());
  world.barrier();
  //get time
  breadth_first_search(g, start, visitor(bfs_visitor<null_visitor>()));
  world.barrier();
  //get time

  return 0;
}
