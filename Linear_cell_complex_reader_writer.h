/*
 Copyright (c) 2011-2014 Ken Arroyo Ohori
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

#ifndef LCC_READER_WRITER_H
#define LCC_READER_WRITER_H

#include <fstream>
#include <sstream>
#include <time.h>

#include <gdal/ogrsf_frmts.h>

#include <CGAL/Linear_cell_complex.h>

#include "Centroid_computer.h"
#include "Point_fitter.h"

template <class LCC, unsigned int dimension = LCC::dimension>
struct Obj_writer {
  static void write_obj(LCC &lcc, const char *path) {
    std::cerr << "Error: Cannot write non 2D LCC!" << std::endl;
  }
};

template <class LCC>
struct Obj_writer<LCC, 2> {
  
  static void write_obj(LCC &lcc, const char *path) {
    CGAL_assertion(LCC::dimension >= 2 && LCC::dimension <= 3);
    
    srand(time(NULL));
    
    std::cout << "Writing file " << path << "..." << std::endl;
    std::ofstream output_file_stream;
    output_file_stream.open(path);
    if (!output_file_stream.is_open()) {
      std::cerr << "Could not create file." << std::endl;
      return;
    }
    
    std::string material_path(path);
    material_path.replace(material_path.end()-3, material_path.end(), "mtl");
    std::cout << "Writing file " << material_path << "..." << std::endl;
    std::ofstream material_file_stream;
    material_file_stream.open(material_path);
    if (!material_file_stream.is_open()) {
      std::cerr << "Could not create file." << std::endl;
      return;
    }
    
    output_file_stream << "##" << std::endl;
    output_file_stream << "## lcc-qb" << std::endl;
    output_file_stream << "## " << path << std::endl;
    output_file_stream << "##" << std::endl;
    output_file_stream << std::endl;
    
    material_file_stream << "##" << std::endl;
    material_file_stream << "## lcc-qb" << std::endl;
    material_file_stream << "## " << material_path << std::endl;
    material_file_stream << "##" << std::endl;
    material_file_stream << std::endl;
    
    output_file_stream << "mtllib " << material_path.substr(material_path.rfind("/")+1) << std::endl;
    output_file_stream << std::endl;
    
    std::map<typename LCC::Vertex_attribute_const_handle, unsigned int> vertex_map;
    unsigned int material_counter = 0;
    
    typename LCC::Vertex_attribute_const_range &vertices = lcc.vertex_attributes();
    Point_fitter<LCC, typename LCC::Point, typename LCC::FT> pf;
    for (typename LCC::Vertex_attribute_const_range::const_iterator current_vertex = vertices.begin(); current_vertex != vertices.end(); ++current_vertex) {
      pf.add_point(current_vertex->point());
    }
    
    for (typename LCC::Vertex_attribute_const_range::const_iterator current_vertex = vertices.begin(); current_vertex != vertices.end(); ++current_vertex) {
      typename LCC::Point new_point = pf.get_point(current_vertex->point());
      output_file_stream << "v " << new_point.x() << " " << new_point.y() << " 0" << std::endl;
      vertex_map[current_vertex] = vertex_map.size()+1;
    }
    
    output_file_stream << std::endl;
    output_file_stream << "g lcc" << std::endl;
    
    typename LCC::template Attribute_const_range<2>::type &facets = lcc.template attributes<2>();
    for (typename LCC::template Attribute_const_range<2>::type::const_iterator current_facet = facets.begin(); current_facet != facets.end(); ++current_facet) {
      
      material_file_stream << "newmtl mat" << material_counter << std::endl;
      material_file_stream << "Ka 0.2 0.2 0.2" << std::endl;
      material_file_stream << "Kd " << float(rand())/RAND_MAX << " " << float(rand())/RAND_MAX << " " << float(rand())/RAND_MAX << std::endl;
      output_file_stream << "usemtl mat" << material_counter << std::endl;
      ++material_counter;
      
      typename LCC::Dart_const_handle dart_in_facet = current_facet->dart();
      typename LCC::template Dart_of_cell_const_range<2, 2> darts_in_facet = lcc.template darts_of_cell<2, 2>(dart_in_facet);
      output_file_stream << "f";
      for (typename LCC::template Dart_of_cell_const_range<2, 2>::const_iterator current_dart = darts_in_facet.begin(); current_dart != darts_in_facet.end(); ++current_dart) {
        output_file_stream << " " << vertex_map[current_dart->template attribute<0>()];
      } output_file_stream << std::endl;
    }
    
    output_file_stream << std::endl;
    material_file_stream << std::endl;
    
    output_file_stream.close();
    material_file_stream.close();
  }
  
  static void write_darts_as_obj(LCC &lcc, const char *path) {
    CGAL_assertion(LCC::dimension >= 2 && LCC::dimension <= 3);
    
    srand(time(NULL));
    
    std::cout << "Writing file " << path << "..." << std::endl;
    std::ofstream output_file_stream;
    output_file_stream.open(path);
    if (!output_file_stream.is_open()) {
      std::cerr << "Could not create file." << std::endl;
      return;
    }
    
    std::string material_path(path);
    material_path.replace(material_path.end()-3, material_path.end(), "mtl");
    std::cout << "Writing file " << material_path << "..." << std::endl;
    std::ofstream material_file_stream;
    material_file_stream.open(material_path);
    if (!material_file_stream.is_open()) {
      std::cerr << "Could not create file." << std::endl;
      return;
    }
    
    output_file_stream << "##" << std::endl;
    output_file_stream << "## lcc-tools" << std::endl;
    output_file_stream << "## " << path << std::endl;
    output_file_stream << "##" << std::endl;
    output_file_stream << std::endl;
    
    material_file_stream << "##" << std::endl;
    material_file_stream << "## lcc-tools" << std::endl;
    material_file_stream << "## " << material_path << std::endl;
    material_file_stream << "##" << std::endl;
    material_file_stream << std::endl;
    
    output_file_stream << "mtllib " << material_path.substr(material_path.rfind("/")+1) << std::endl;
    output_file_stream << std::endl;
    
    std::map<typename LCC::Vertex_attribute_const_handle, unsigned int> vertex_map;
    //    std::map<typename LCC::template Attribute_const_handle<1>::type, unsigned int> edge_map;
    std::map<typename LCC::template Attribute_const_handle<2>::type, unsigned int> facet_map;
    unsigned int material_counter = 0;
    typename LCC::Point p;
    
    typename LCC::Vertex_attribute_const_range &vertices = lcc.vertex_attributes();
    Point_fitter<LCC, typename LCC::Point, typename LCC::FT> pf;
    for (typename LCC::Vertex_attribute_const_range::const_iterator current_vertex = vertices.begin(); current_vertex != vertices.end(); ++current_vertex) {
      pf.add_point(current_vertex->point());
    }
    
    for (typename LCC::Vertex_attribute_const_range::const_iterator current_vertex = vertices.begin(); current_vertex != vertices.end(); ++current_vertex) {
      typename LCC::Point new_point = pf.get_point(current_vertex->point());
      output_file_stream << "v " << new_point.x() << " " << new_point.y() << " 0" << std::endl;
      vertex_map[current_vertex] = vertex_map.size()+1;
    }
    
    //    typename LCC::template Attribute_const_range<1>::type &edges = lcc.template attributes<1>();
    //    for (typename LCC::template Attribute_const_range<1>::type::const_iterator current_edge = edges.begin(); current_edge != edges.end(); ++current_edge) {
    //      p = lcc.template barycenter<1>(current_edge->dart());
    //      output_file_stream << "v " << p.x() << " " << p.y() << " 0" << std::endl;
    //      edge_map[current_edge] = vertex_map.size()+edge_map.size()+1;
    //    }
    
    typename LCC::template Attribute_const_range<2>::type &facets = lcc.template attributes<2>();
    for (typename LCC::template Attribute_const_range<2>::type::const_iterator current_facet = facets.begin(); current_facet != facets.end(); ++current_facet) {
      p = pf.get_point(lcc.template barycenter<2>(current_facet->dart()));
      output_file_stream << "v " << p.x() << " " << p.y() << " 0" << std::endl;
      facet_map[current_facet] = vertex_map.size()+facet_map.size()+1;
    }
    
    output_file_stream << std::endl;
    output_file_stream << "g lcc" << std::endl;
    
    int mark = lcc.get_new_mark();
    CGAL_precondition(lcc.is_whole_map_unmarked(mark));
    for (typename LCC::Dart_const_range::const_iterator current_dart = lcc.darts().begin(); current_dart != lcc.darts().end(); ++current_dart) {
      if (lcc.is_marked(current_dart, mark)) continue;
      
      material_file_stream << "newmtl mat" << material_counter << std::endl;
      material_file_stream << "Ka 0.2 0.2 0.2" << std::endl;
      material_file_stream << "Kd " << float(rand())/RAND_MAX << " " << float(rand())/RAND_MAX << " " << float(rand())/RAND_MAX << std::endl;
      output_file_stream << "usemtl mat" << material_counter << std::endl;
      ++material_counter;
      
      for (typename LCC::template Dart_of_cell_const_range<2, 2>::const_iterator current_dart_in_cell = lcc.template darts_of_cell<2, 2>(current_dart).begin();
           current_dart_in_cell != lcc.template darts_of_cell<2, 2>(current_dart).end();
           ++current_dart_in_cell) {
        lcc.mark(current_dart_in_cell, mark);
        output_file_stream << "f " << vertex_map[current_dart_in_cell->template attribute<0>()] << " " << vertex_map[current_dart_in_cell->beta(1)->template attribute<0>()] << " " << facet_map[current_dart_in_cell->template attribute<2>()] << std::endl;
      }
    } lcc.negate_mark(mark);
    CGAL_postcondition(lcc.is_whole_map_unmarked(mark));
    
    output_file_stream << std::endl;
    material_file_stream << std::endl;
    
    output_file_stream.close();
    material_file_stream.close();
  }
};

template <class LCC>
struct Obj_writer<LCC, 3> {
  
  static void write_obj(LCC &lcc, const char *path) {
    CGAL_assertion(LCC::dimension == 3);
    
    srand(time(NULL));
    
    std::cout << "Writing file " << path << "..." << std::endl;
    std::ofstream output_file_stream;
    output_file_stream.open(path);
    if (!output_file_stream.is_open()) {
      std::cerr << "Could not create file." << std::endl;
      return;
    }
    
    std::string material_path(path);
    material_path.replace(material_path.end()-3, material_path.end(), "mtl");
    std::cout << "Writing file " << material_path << "..." << std::endl;
    std::ofstream material_file_stream;
    material_file_stream.open(material_path);
    if (!material_file_stream.is_open()) {
      std::cerr << "Could not create file." << std::endl;
      return;
    }
    
    output_file_stream << "##" << std::endl;
    output_file_stream << "## lcc-tools" << std::endl;
    output_file_stream << "## " << path << std::endl;
    output_file_stream << "##" << std::endl;
    output_file_stream << std::endl;
    
    material_file_stream << "##" << std::endl;
    material_file_stream << "## lcc-tools" << std::endl;
    material_file_stream << "## " << material_path << std::endl;
    material_file_stream << "##" << std::endl;
    material_file_stream << std::endl;
    
    material_file_stream << "newmtl black" << std::endl;
    material_file_stream << "Ka 0.2 0.2 0.2" << std::endl;
    material_file_stream << "Kd 0 0 0" << std::endl;
    
    output_file_stream << "mtllib " << material_path.substr(material_path.rfind("/")+1) << std::endl;
    output_file_stream << std::endl;
    
    std::map<typename LCC::Vertex_attribute_const_handle, unsigned int> vertex_map;
    std::map<typename LCC::template Attribute_const_handle<3>::type, unsigned int> volume_map;
    
    typename LCC::Vertex_attribute_const_range &vertices = lcc.vertex_attributes();
    Point_fitter<LCC, typename LCC::Point, typename LCC::FT> pf;
    for (typename LCC::Vertex_attribute_const_range::const_iterator current_vertex = vertices.begin(); current_vertex != vertices.end(); ++current_vertex) {
      pf.add_point(current_vertex->point());
    }
    
    for (typename LCC::Vertex_attribute_const_range::const_iterator current_vertex = vertices.begin(); current_vertex != vertices.end(); ++current_vertex) {
      typename LCC::Point new_point = pf.get_point(current_vertex->point());
      output_file_stream << "v " << new_point.x() << " " << new_point.y() << " " << new_point.z() <<  std::endl;
      vertex_map[current_vertex] = vertex_map.size()+1;
    }
    
    output_file_stream << std::endl;
    output_file_stream << "g lcc" << std::endl;
    
    typename LCC::template Attribute_const_range<2>::type &facets = lcc.template attributes<2>();
    for (typename LCC::template Attribute_const_range<2>::type::const_iterator current_facet = facets.begin(); current_facet != facets.end(); ++current_facet) {
      typename LCC::Dart_const_handle dart_in_facet = current_facet->dart();
      
      if (!lcc.is_free(dart_in_facet, 3)) output_file_stream << "usemtl black" << std::endl;
      else if (volume_map.count(lcc.template attribute<3>(dart_in_facet)) > 0) output_file_stream << "usemtl mat" << volume_map[lcc.template attribute<3>(dart_in_facet)] << std::endl;
      else {
        volume_map[lcc.template attribute<3>(dart_in_facet)] = volume_map.size();
        material_file_stream << "newmtl mat" << volume_map[lcc.template attribute<3>(dart_in_facet)] << std::endl;
        material_file_stream << "Ka 0.2 0.2 0.2" << std::endl;
        material_file_stream << "Kd " << float(rand())/RAND_MAX << " " << float(rand())/RAND_MAX << " " << float(rand())/RAND_MAX << std::endl;
        output_file_stream << "usemtl mat" << volume_map[lcc.template attribute<3>(dart_in_facet)] << std::endl;
      }
      
      typename LCC::template Dart_of_cell_const_range<2, 2> darts_in_facet = lcc.template darts_of_cell<2, 2>(dart_in_facet);
      output_file_stream << "f";
      for (typename LCC::template Dart_of_cell_const_range<2, 2>::const_iterator current_dart = darts_in_facet.begin(); current_dart != darts_in_facet.end(); ++current_dart) {
        output_file_stream << " " << vertex_map[current_dart->template attribute<0>()];
      } output_file_stream << std::endl;
    }
    
    output_file_stream << std::endl;
    material_file_stream << std::endl;
    
    output_file_stream.close();
    material_file_stream.close();
  }
  
  static void write_darts_as_obj(LCC &lcc, const char *path) {
    CGAL_assertion(LCC::dimension == 3);
    
    srand(time(NULL));
    
    std::cout << "Writing file " << path << "..." << std::endl;
    std::ofstream output_file_stream;
    output_file_stream.open(path);
    if (!output_file_stream.is_open()) {
      std::cerr << "Could not create file." << std::endl;
      return;
    }
    
    std::string material_path(path);
    material_path.replace(material_path.end()-3, material_path.end(), "mtl");
    std::cout << "Writing file " << material_path << "..." << std::endl;
    std::ofstream material_file_stream;
    material_file_stream.open(material_path);
    if (!material_file_stream.is_open()) {
      std::cerr << "Could not create file." << std::endl;
      return;
    }
    
    output_file_stream << "##" << std::endl;
    output_file_stream << "## lcc-tools" << std::endl;
    output_file_stream << "## " << path << std::endl;
    output_file_stream << "##" << std::endl;
    output_file_stream << std::endl;
    
    material_file_stream << "##" << std::endl;
    material_file_stream << "## lcc-tools" << std::endl;
    material_file_stream << "## " << material_path << std::endl;
    material_file_stream << "##" << std::endl;
    material_file_stream << std::endl;
    
    material_file_stream << "newmtl black" << std::endl;
    material_file_stream << "Ka 0.2 0.2 0.2" << std::endl;
    material_file_stream << "Kd 0 0 0" << std::endl;
    
    std::map<typename LCC::Vertex_attribute_const_handle, unsigned int> vertex_map;
    //    std::map<typename LCC::template Attribute_const_handle<1>::type, unsigned int> edge_map;
    std::map<typename LCC::template Attribute_const_handle<2>::type, unsigned int> facet_map;
    std::map<typename LCC::template Attribute_const_handle<3>::type, unsigned int> volume_map;
    typename LCC::Point p;
    
    typename LCC::Vertex_attribute_const_range &vertices = lcc.vertex_attributes();
    Point_fitter<LCC, typename LCC::Point, typename LCC::FT> pf;
    for (typename LCC::Vertex_attribute_const_range::const_iterator current_vertex = vertices.begin(); current_vertex != vertices.end(); ++current_vertex) {
      pf.add_point(current_vertex->point());
    }
    
    for (typename LCC::Vertex_attribute_const_range::const_iterator current_vertex = vertices.begin(); current_vertex != vertices.end(); ++current_vertex) {
      typename LCC::Point new_point = pf.get_point(current_vertex->point());
      output_file_stream << "v " << new_point.x() << " " << new_point.y() << " " << new_point.z() <<  std::endl;
      vertex_map[current_vertex] = vertex_map.size()+1;
    }
    
    //    typename LCC::template Attribute_const_range<1>::type &edges = lcc.template attributes<1>();
    //    for (typename LCC::template Attribute_const_range<1>::type::const_iterator current_edge = edges.begin(); current_edge != edges.end(); ++current_edge) {
    //      p = lcc.template barycenter<1>(current_edge->dart());
    //      output_file_stream << "v " << p.x() << " " << p.y() << " 0" << std::endl;
    //      edge_map[current_edge] = vertex_map.size()+edge_map.size()+1;
    //    }
    
    typename LCC::template Attribute_const_range<2>::type &facets = lcc.template attributes<2>();
    for (typename LCC::template Attribute_const_range<2>::type::const_iterator current_facet = facets.begin(); current_facet != facets.end(); ++current_facet) {
      //      p = lcc.template barycenter<2>(current_facet->dart());
      p = pf.get_point(Centroid_computer<LCC, typename LCC::Point, typename LCC::FT>::template centroid_of_cell<2>(lcc, current_facet->dart()));
      output_file_stream << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
      facet_map[current_facet] = vertex_map.size()+facet_map.size()+1;
    }
    
    output_file_stream << std::endl;
    output_file_stream << "g lcc" << std::endl;
    
    for (typename LCC::Dart_const_range::const_iterator current_dart = lcc.darts().begin(); current_dart != lcc.darts().end(); ++current_dart) {
      if (lcc.is_free(current_dart, 1)) continue;
      if (current_dart->template attribute<0>() == LCC::null_handle) continue;
      if (current_dart->beta(1)->template attribute<0>() == LCC::null_handle) continue;
      if (current_dart->template attribute<2>() == LCC::null_handle) continue;
      
      if (!lcc.is_free(current_dart, 3)) output_file_stream << "usemtl black" << std::endl;
      else if (volume_map.count(lcc.template attribute<3>(current_dart)) > 0) output_file_stream << "usemtl mat" << volume_map[lcc.template attribute<3>(current_dart)] << std::endl;
      else {
        volume_map[lcc.template attribute<3>(current_dart)] = volume_map.size();
        material_file_stream << "newmtl mat" << volume_map[lcc.template attribute<3>(current_dart)] << std::endl;
        material_file_stream << "Ka 0.2 0.2 0.2" << std::endl;
        material_file_stream << "Kd " << float(rand())/RAND_MAX << " " << float(rand())/RAND_MAX << " " << float(rand())/RAND_MAX << std::endl;
        output_file_stream << "usemtl mat" << volume_map[lcc.template attribute<3>(current_dart)] << std::endl;
      }
      
      output_file_stream << "f " << vertex_map[current_dart->template attribute<0>()] << " " << vertex_map[current_dart->beta(1)->template attribute<0>()] << " " << facet_map[current_dart->template attribute<2>()] << std::endl;
    }
    
    output_file_stream << std::endl;
    material_file_stream << std::endl;
    
    output_file_stream.close();
    material_file_stream.close();
  }
};

template <class LCC_>
class Linear_cell_complex_reader_writer {
public:
  typedef LCC_ LCC;
  typedef typename LCC::Point Point;
  typedef typename LCC::Traits::FT FT;
  
  Linear_cell_complex_reader_writer() {
    
  }
  
  static void write_obj(const char *path, LCC &lcc) {
    Obj_writer<LCC>::write_obj(lcc, path);
  }
  
  static void write_darts_as_obj(const char *path, LCC &lcc) {
    Obj_writer<LCC>::write_darts_as_obj(lcc, path);
  }
};

#endif
