 /**
 * Copyright (c) 2009 Carnegie Mellon University.
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */
#ifndef GRAPHLAB_GRAPH_JSON_PARSER_HPP
#define GRAPHLAB_GRAPH_JSON_PARSER_HPP

#include <string>
#include <sstream>
#include <iostream>
#include <vector>

#if defined(__clang) || (defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)))
#pragma GCC diagnostic push
#endif

#pragma GCC diagnostic ignored "-Wreorder"
#include <libjson/libjson.h>

#if defined(__clang) || (defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)))
#pragma GCC diagnostic pop
#endif

#include <boost/functional.hpp>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/util/hdfs.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/graph/ingress/distributed_identity_ingress.hpp>

namespace graphlab {

  namespace builtin_parsers {
    template <typename EdgeData>
    bool empty_edge_parser(EdgeData& e, const std::string& line) {
      return true;
    }

    template <typename VertexData>
    bool empty_vertex_parser(VertexData& v, const std::string& line) {
      return true;
    }

  }



template<typename VertexData, typename EdgeData>
class distributed_graph;


template <typename VertexData, typename EdgeData>
class json_parser {

  public:
    typedef distributed_graph<VertexData, EdgeData> graph_type;
    typedef EdgeData edge_data_type;
    typedef VertexData vertex_data_type;

    typedef typename graph_type::vertex_id_type vertex_id_type;
    typedef typename graph_type::lvid_type lvid_type;

    typedef boost::function<bool(edge_data_type&, const std::string&)> edge_parser_type;
    typedef boost::function<bool(vertex_data_type&, const std::string&)> vertex_parser_type;
    typedef boost::function<bool(graph_type&, const std::string&)> line_parser_type;


  public:
  json_parser (graph_type& graph, const std::string& prefix, bool gzip=false,
      edge_parser_type edge_parser=builtin_parsers::empty_edge_parser<EdgeData>,
      vertex_parser_type vertex_parser=builtin_parsers::empty_vertex_parser<VertexData>) :
    graph(graph), prefix(prefix), gzip(gzip), edge_parser(edge_parser), vertex_parser(vertex_parser) {
    }


  bool load() {
    line_parser_type graph_structure_parser = parse_graph_structure_from_json;
    line_parser_type vid2lvid_parser = parse_vid2lvid_from_json;
    line_parser_type edata_parser = boost::bind(parse_edatalist_from_json, _1, _2, edge_parser);
    line_parser_type vrecord_parser = boost::bind(parse_vrecord_from_json, _1, _2, vertex_parser);

    bool success = parse_by_line(graphfilename(), graph_structure_parser);
    success = success & parse_by_line(vid2lvidfilename(), vid2lvid_parser); 
    success = success & parse_by_line(edatafilename(), edata_parser); 
    success = success & parse_by_line(vrecordfilename(), vrecord_parser);

    if (!success) {
      logstream(LOG_FATAL) << "Fail parsing graph json" << std::endl;
      return false;
    }

    graph.local_graph.finalized = true;

    
    ASSERT_GE(graph.local_graph.num_vertices(), graph.local_graph.gstore.num_vertices);
    ASSERT_EQ(graph.vid2lvid.size(), graph.local_graph.num_vertices());
    ASSERT_EQ(graph.lvid2record.size(), graph.local_graph.num_vertices());



    logstream(LOG_INFO) << "Finished loading graph" << graph.procid() 
      << "\n\tnverts: " << graph.num_local_own_vertices() 
      << "\n\tnreplicas: " << graph.local_graph.num_vertices() 
      << "\n\tnedges: " << graph.local_graph.num_edges() 
      << std::endl;


    if (graph.ingress_ptr == NULL) {
      graph.ingress_ptr = new distributed_identity_ingress<VertexData, EdgeData>(graph.rpc.dc(), graph);
    }

    graph.ingress_ptr->exchange_global_info();
    delete graph.ingress_ptr;


    graph.finalized = true;
    return success;
  }

  bool parse_by_line (const std::string& srcfilename, line_parser_type line_parser) {
    std::string fname;
    //check for "/" ending in directory"
    if(!boost::ends_with(prefix,"/"))
        fname = prefix + "/" + srcfilename;
    else 
        fname = prefix + srcfilename;

    logstream(LOG_INFO) << "Load graph json from " << fname << std::endl;

    boost::iostreams::filtering_stream<boost::iostreams::input> fin;
    // loading from hdfs
    if (boost::starts_with(prefix, "hdfs://")) {
      graphlab::hdfs& hdfs = hdfs::get_hdfs();
      graphlab::hdfs::fstream in_file(hdfs, fname);
      if (!in_file.good()) {
        logstream(LOG_FATAL) << "Fail to open file " << fname << std::endl;
        return false;
      }

      if (gzip) fin.push(boost::iostreams::gzip_decompressor());
      fin.push(in_file);

      if (!fin.good()) {
        logstream(LOG_FATAL) << "Fail to read from stream " << fname << std::endl;
        return false;
      }

      load_from_stream(fname, fin, line_parser);

      if (gzip) fin.pop();
      fin.pop();
    } else { // loading from disk
      std::ifstream in_file(fname.c_str(),
          std::ios_base::in | std::ios_base::binary);

      if (!in_file.good()) {
        logstream(LOG_FATAL) << "Fail to open file " << fname << std::endl;
        return false;
      }

      if (gzip) fin.push(boost::iostreams::gzip_decompressor());
      fin.push(in_file);

      if (!fin.good()) {
        logstream(LOG_FATAL) << "Fail to read from stream " << fname << std::endl;
        return false;
      }

      load_from_stream(fname, fin, line_parser);

      if (gzip) fin.pop();
      fin.pop();
    }


    return true;
  }
    /**
       \internal
       This internal function is used to load a single line from an input stream
     */
    template<typename Fstream>
    bool load_from_stream(std::string filename, Fstream& fin,
                          line_parser_type& line_parser) {
      size_t linecount = 0;
      timer ti; ti.start();
      while(fin.good() && !fin.eof()) {
        std::string line;
        std::getline(fin, line);
        if(line.empty()) continue;
        if(fin.fail()) break;
        const bool success = line_parser(graph, line);
        if (!success) {
          logstream(LOG_WARNING)
            << "Error parsing line " << linecount << " in "
            << filename << ": " << std::endl
            << "\t\"" << line << "\"" << std::endl;
          return false;
        }
        ++linecount;
        if (ti.current_time() > 5.0) {
          logstream(LOG_INFO) << linecount << " Lines read" << std::endl;
          ti.start();
        }
      }
      return true;
    } // end of load from stream



  /* Parse the graph structure from json */
  static bool parse_graph_structure_from_json (graph_type& graph, const std::string& str) {
    JSONNode n = libjson::parse(str);
    JSONNode::const_iterator i = n.begin();
    typedef typename graph_type::local_graph_type local_graph_type;
    local_graph_type& local_graph = graph.get_local_graph();
    while(i != n.end()) {
      if (i->name() == "numEdges") {
        local_graph.gstore.num_edges = i->as_int();
      } else if (i->name() == "numVertices") {
        local_graph.gstore.num_vertices= i->as_int();
      } else if (i->name() == "csr") {
        // parse rowIndex -> graph.local_graph.gstore.csr_source
        // parse colIndex -> graph.local_graph.gstore.csr_target
        JSONNode csr = *i;
        JSONNode::const_iterator j = csr.begin();
        while (j != csr.end()) {
          if (j->name() == "rowIndex") {
              parse_vid_array (local_graph.gstore.CSR_src, *j);
          } else if (j->name() == "colIndex") {
              parse_vid_array (local_graph.gstore.CSR_dst, *j);
          } else {
              logstream(LOG_ERROR) << "Error parsing json into graph. Unknown json node name:"
                << "CSR:" << j->name() << std::endl;
          }
          ++j;
        }
      } else if (i->name() == "csc") {
        // parse rowIndex -> graph.local_graph.gstore.csc_target
        // parse colIndex -> graph.local_graph.gstore.csc_source
        JSONNode csc = *i;
        JSONNode::const_iterator j = csc.begin();

        while (j != csc.end()) {
          if (j->name() == "rowIndex") {
              parse_vid_array (local_graph.gstore.CSC_dst, *j);
          } else if (j->name() == "colIndex") {
              parse_vid_array (local_graph.gstore.CSC_src, *j);
          } else {
              logstream(LOG_ERROR) << "Error parsing json into graph. Unknown json node name:"
                << "CSC:"<<j->name() << std::endl;
          }
          ++j;
        }
      } else if (i->name() == "c2rMap") {
        parse_vid_array (local_graph.gstore.c2r_map, *i);
      } else {
        logstream(LOG_ERROR) << "Error parsing json into graph. Unknown json node name:" <<
          i->name() << std::endl;
      }
      ++i;
    } // end while


    ASSERT_EQ(local_graph.gstore.num_edges, local_graph.gstore.c2r_map.size());
    ASSERT_EQ(local_graph.gstore.num_edges, local_graph.gstore.CSR_dst.size());
    ASSERT_EQ(local_graph.gstore.num_edges, local_graph.gstore.CSC_src.size());

    graph.lvid2record.reserve(local_graph.gstore.num_vertices);
    graph.lvid2record.resize(local_graph.gstore.num_vertices);
    local_graph.reserve(local_graph.gstore.num_vertices);
    // local_graph.finalized = true;
    return true;
  }

  /* Parse the vid2lvid map from json */
  static bool parse_vid2lvid_from_json (graph_type& graph, const std::string& str) {
    JSONNode n = libjson::parse(str);
    JSONNode::const_iterator i = n.begin();
    typedef typename graph_type::local_graph_type local_graph_type;
    while(i != n.end()) {
      if (i->name() == "vid2lvid") {
        JSONNode::const_iterator j = i->begin();
        typename graph_type::vid2lvid_map_type & map = graph.vid2lvid;
        map.clear();
        while(j != i->end()) {
          graph.vid2lvid[boost::lexical_cast<vertex_id_type>(j->name())] = (boost::lexical_cast<lvid_type>)(j->as_int());
          ++j;
        }
      } else {
        // report error
        return false;
      }
      ++i;
    }
    return true;
  }

  /* Parse the edata list from json */
  static bool parse_edatalist_from_json (graph_type& graph, const std::string& str,
      edge_parser_type edge_parser) {
    JSONNode n = libjson::parse(str);
    JSONNode::const_iterator i = n.begin();
    typedef typename graph_type::local_graph_type local_graph_type;
    local_graph_type& local_graph = graph.get_local_graph();
    while(i != n.end()) {
       if (i->name() == "edataList") {
        // parse  edatalist -> graph.local_graph.gstore.edata
        JSONNode edatanode= *i;
        JSONNode::const_iterator j = edatanode.begin();
        std::vector<edge_data_type>& edatalist = local_graph.gstore.edge_data_list;
        edatalist.clear();
        edatalist.reserve(local_graph.gstore.num_edges);
        edge_data_type e;
        while (j != edatanode.end()) {
          edge_parser(e, j->as_string());
          edatalist.push_back(e);
          ++j;
        }
       } else {
         return false;
         // report error
       }
      ++i;
    }
    return true;
  }


  /* Parse the vertex record list from json */
  static bool parse_vrecord_from_json (graph_type& graph, const std::string& str,
      vertex_parser_type vertex_parser) {

    typedef typename graph_type::local_graph_type local_graph_type;
    local_graph_type& local_graph = graph.get_local_graph();

    JSONNode n = libjson::parse(str);
    JSONNode::const_iterator i = n.begin();

    vertex_data_type vdata;
    typename graph_type::vertex_record vrecord;

    while (i != n.end()) {
      if (i->name() == "mirrors") {
        JSONNode::const_iterator j = (*i).begin();
        while (j != (*i).end()) {
          int mirror = j->as_int();
          vrecord._mirrors.set_bit((procid_t)mirror);
          ++j;
        }
      } else if (i->name() == "inEdges") {
        vrecord.num_in_edges = i->as_int();
      } else if (i->name() == "outEdges") {
        vrecord.num_out_edges = i->as_int();
      } else if (i->name() == "gvid") {
        // Check unsafe
        vrecord.gvid = boost::lexical_cast<vertex_id_type>(i->as_int());
      } else if (i->name() == "owner") {
        vrecord.owner = (procid_t)i->as_int();
      } else if (i->name() == "VertexData") {
        if (!(i->type() == JSON_NULL))
          vertex_parser(vdata, i->as_string());
      } else {
        logstream(LOG_ERROR) << "Error parsing json into vrecord. Unknown json node name:" <<
          i->name() << std::endl;
      }
      ++i;
    }

      if (graph.vid2lvid.find(vrecord.gvid) == graph.vid2lvid.end()) {
        // Check if this a singlton node
        // ignore for now
        logstream(LOG_WARNING) << "Singleton node detected: gvid = " << vrecord.gvid << ". Ignored" << std::endl;
      } else {
        lvid_type lvid = graph.vid2lvid[vrecord.gvid];
        graph.lvid2record[lvid] = vrecord;
        local_graph.add_vertex(lvid, vdata);
        if (vrecord.owner == graph.procid()) ++graph.local_own_nverts;
      }

    return true;
  }


  /*  Helper function starts here  */
  private:
  /* Parse an json integer array and copy to a int vector */
  static bool parse_vid_array (std::vector<vertex_id_type>& to, const JSONNode& n) {
    if (n.type() != JSON_ARRAY) return false;
    to.clear();
    JSONNode::const_iterator i = n.begin();
    while (i != n.end()) {
      to.push_back(i->as_int());
      ++i;
    }
    return true;
  }


  std::string zeropadding(const std::string& s, int width) {
      ASSERT_LE(s.length(), width);
      std::ostringstream ss;
      ss << std::setw(width) << std::setfill('0') << s;
      return ss.str();
  } 


  const std::string graphfilename() {
    procid_t pid = graph.procid();
    std::string suffix =  gzip ? ".gz" : "";
    return "graph/graph"+tostr(pid)+"-r-"+zeropadding(tostr(pid), 5)+suffix;
    // return "graph/graph"+tostr(pid)+"-r-00000"+suffix;
  }

  const std::string vid2lvidfilename() {
    procid_t pid = graph.procid();
    std::string suffix =  gzip ? ".gz" : "";
    return "graph/vid2lvid"+tostr(pid)+"-r-"+zeropadding(tostr(pid), 5)+suffix;
    // return "graph/vid2lvid"+tostr(pid)+"-r-00000"+suffix;
  }

  const std::string edatafilename() {
    procid_t pid = graph.procid();
    std::string suffix =  gzip ? ".gz" : "";
    return "graph/edata"+tostr(pid)+"-r-"+zeropadding(tostr(pid),5)+suffix;
    // return "graph/edata"+tostr(pid)+"-r-00000"+suffix;
  }

  const std::string vrecordfilename() {
    procid_t pid = graph.procid();
    std::string suffix =  gzip ? ".gz" : "";
    return "vrecord/vdata"+tostr(pid)+"-r-"+zeropadding(tostr(pid),5)+suffix;
    // return "vrecord/vdata"+tostr(pid)+"-r-00000"+ suffix;
  }

  private:
    graph_type& graph;
    std::string prefix;
    bool gzip;
    edge_parser_type edge_parser;
    vertex_parser_type vertex_parser;
}; // json_parser


} // namespace graphlab
#endif
