#include <iostream>
#include <utility>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace boost;

int main(int argc, char* argv[]){
 int num_verts = 5;
 if(argc > 1){
  num_verts = atoi(argv[1]);
 }
// std::cout << num_verts << std::endl;
 typedef adjacency_list<vecS,vecS,directedS> Graph;
 Graph g(num_verts);
 typedef std::pair<int, int> Edge;

 int num_edg = 0;

 for(int i = 0; i < num_verts - 2; ++i){
  add_edge(i, i+2, g);
  ++num_edg;
  if(i % 5 == 0){
   add_edge(i, i+1, g);
   ++num_edg;
  }
 }

 typedef property_map<Graph, vertex_index_t>::type IndexMap;
 IndexMap index = get(vertex_index, g);

 std::cout << "my edges are: ";
 graph_traits<Graph>::edge_iterator e_i, e_end;
 for(tie(e_i, e_end) = edges(g); e_i != e_end; ++e_i){
  std::cout << "(" << index[source(*e_i, g)] << "," << index[target(*e_i, g)] << ") ";
 }
 std::cout << std::endl;

}
