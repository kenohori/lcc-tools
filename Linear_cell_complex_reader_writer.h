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

#include "Linear_cell_complex_incremental_builder.h"

template <class LCCRW, unsigned int dimension = 0, class Attribute_type = typename LCCRW::LCC::Vertex_attribute>
struct Cell_complex_parser;

template <class LCCRW, unsigned int dimension, class Attribute_type>
struct Cell_complex_parser {
  static void parse(LCCRW &lccrw, std::ifstream &input_file_stream, std::vector<std::map<long, typename LCCRW::LCC::Dart_handle> > &cells) {
    std::string current_line, current_word;
    unsigned long number_of_attributes, current_index;
    std::istringstream line_stream;
    unsigned long current_number, remaining_cells_in_this_dimension;
    bool reading_header = true;
    
    while (!input_file_stream.eof()) {
      getline(input_file_stream, current_line);
      std::cout << "Line: \'" << current_line << "\'" << std::endl;
      
      // Empty line
      if (current_line.size() == 0) {
        continue;
      }
      
      // Just a newline sequence
      if (current_line[0] == '\r' || current_line[0] == '\n') {
        continue;
      }
      
      line_stream = std::istringstream(current_line + "\n");
      
      if (reading_header) {
        line_stream >> current_word;
        current_number = string_to_any(current_word, current_number);
        remaining_cells_in_this_dimension = current_number;
        
        line_stream >> current_word;
        current_number = string_to_any(current_word, current_number);
        number_of_attributes = current_number;
        
        std::cout << remaining_cells_in_this_dimension << " " << dimension << "-cells with " << number_of_attributes << " attributes." << std::endl;
        reading_header = false;
        cells.push_back(std::map<long, typename LCCRW::LCC::Dart_handle>());
        Index_builder<typename LCCRW::LCCB, dimension>::build_for_dimension(lccrw.builder);
      }
      
      else {
        line_stream >> current_word;
        current_index = string_to_any(current_word, current_number);
        
        std::list<typename LCCRW::LCC::Dart_handle> faces_in_cell;
        while (line_stream.peek() != '\n') {
          line_stream >> current_word;
          current_number = string_to_any(current_word, current_number);
          faces_in_cell.push_back(cells[dimension-1][current_number]);
        }
        
        typename LCCRW::LCC::template Attribute_handle<dimension>::type cell_attribute = lccrw.lcc.template create_attribute<dimension>(-1);
        std::cout << "Created " << dimension << "-cell<" << &*cell_attribute << "> (" << std::endl;
        cells[dimension][current_index] = lccrw.builder.template get_cell<dimension>(cell_attribute, faces_in_cell.begin(), faces_in_cell.end()).first;
        cells[dimension][current_index]->template attribute<dimension>()->info() = current_index;
        std::cout << dimension << "-cell<" << &*cells[dimension][current_index]->template attribute<dimension>() << ">[" << cells[dimension][current_index]->template attribute<dimension>()->info() << "]" << std::endl;
        --remaining_cells_in_this_dimension;
      }
      
      if (remaining_cells_in_this_dimension == 0) Cell_complex_parser<LCCRW, dimension+1, typename LCCRW::LCC::template Attribute_type<dimension+1>::type>::parse(lccrw, input_file_stream, cells);
    }
  }
};

template <class LCCRW, unsigned int dimension>
struct Cell_complex_parser<LCCRW, dimension, CGAL::Void> {
  static void parse(LCCRW &lccrw, std::ifstream &input_file_stream, std::vector<std::map<long, typename LCCRW::LCC::Dart_handle> > &cells) {
    std::string current_line;
    
    while (!input_file_stream.eof()) {
      getline(input_file_stream, current_line);
      std::cout << "Line: \'" << current_line << "\' (unprocessed)" << std::endl;
    }
  }
};

template <class LCCRW>
struct Cell_complex_parser<LCCRW, 0, typename LCCRW::LCC::Vertex_attribute> {
  static void parse(LCCRW &lccrw, std::ifstream &input_file_stream, std::vector<std::map<long, typename LCCRW::LCC::Dart_handle> > &cells) {
    std::string current_line, current_word;
    unsigned long number_of_attributes, current_index;
    std::istringstream line_stream;
    unsigned long dimensions, current_number, remaining_cells_in_this_dimension;
    bool reading_header = true;
    typename LCCRW::FT current_value;
    
    while (!input_file_stream.eof()) {
      getline(input_file_stream, current_line);
      std::cout << "Line: \'" << current_line << "\'" << std::endl;
      
      // Empty line
      if (current_line.size() == 0) {
        continue;
      }
      
      // Just a newline sequence
      if (current_line[0] == '\r' || current_line[0] == '\n') {
        continue;
      }
      
      line_stream = std::istringstream(current_line + "\n");
      
      if (reading_header) {
        line_stream >> current_word;
        current_number = string_to_any(current_word, current_number);
        remaining_cells_in_this_dimension = current_number;
        
        line_stream >> current_word;
        dimensions = string_to_any(current_word, dimensions);
        std::cout << "Cell complex has " << dimensions << " dimensions." << std::endl;
        
        line_stream >> current_word;
        current_number = string_to_any(current_word, current_number);
        number_of_attributes = current_number;
        
        std::cout << remaining_cells_in_this_dimension << " vertices with " << number_of_attributes << " attributes." << std::endl;
        reading_header = false;
        cells.push_back(std::map<long, typename LCCRW::LCC::Dart_handle>());
        Index_builder<typename LCCRW::LCCB, 0>::build_for_dimension(lccrw.builder);
      }
      
      else {
        line_stream >> current_word;
        current_index = string_to_any(current_word, current_index);
        
        Point_creator<typename LCCRW::Point, typename LCCRW::FT> creator;
        for (unsigned int current_coordinate = 0; current_coordinate < LCCRW::dimension && line_stream.peek() != '\n'; ++current_coordinate) {
          line_stream >> current_word;
          current_value = string_to_any(current_word, current_value);
          creator.push_back(current_value);
        }
        
        typename LCCRW::LCC::Vertex_attribute_handle vertex_attribute = lccrw.lcc.template create_attribute<0>(creator.create_point(), -1);
        //std::cout << "Created vertex<" << &*vertex_attribute << "> (" << vertex_attribute->point() << ")" << std::endl;
        cells[0][current_index] = lccrw.builder.get_vertex(vertex_attribute);
        cells[0][current_index]->template attribute<0>()->info() = current_index;
        std::cout << "\tvertex<" << &*cells[0][current_index]->template attribute<0>() << ">[" << cells[0][current_index]->template attribute<0>()->info() << "] (" << cells[0][current_index]->template attribute<0>()->point() << ")" << std::endl;
        --remaining_cells_in_this_dimension;
      }
      
      if (remaining_cells_in_this_dimension == 0) Cell_complex_parser<LCCRW, 1, typename LCCRW::LCC::template Attribute_type<1>::type>::parse(lccrw, input_file_stream, cells);
    }
  }
};

template <class LCCRW>
struct Cell_complex_parser<LCCRW, 1, typename LCCRW::LCC::template Attribute_type<1>::type> {
  static void parse(LCCRW &lccrw, std::ifstream &input_file_stream, std::vector<std::map<long, typename LCCRW::LCC::Dart_handle> > &cells) {
    std::string current_line, current_word;
    unsigned long number_of_attributes, current_index;
    std::istringstream line_stream;
    unsigned long current_number, remaining_cells_in_this_dimension;
    bool reading_header = true;
    
    while (!input_file_stream.eof()) {
      getline(input_file_stream, current_line);
      std::cout << "Line: \'" << current_line << "\'" << std::endl;
      
      // Empty line
      if (current_line.size() == 0) {
        continue;
      }
      
      // Just a newline sequence
      if (current_line[0] == '\r' || current_line[0] == '\n') {
        continue;
      }
      
      line_stream = std::istringstream(current_line + "\n");
      
      if (reading_header) {
        line_stream >> current_word;
        current_number = string_to_any(current_word, current_number);
        remaining_cells_in_this_dimension = current_number;
        
        line_stream >> current_word;
        current_number = string_to_any(current_word, current_number);
        number_of_attributes = current_number;
        
        std::cout << remaining_cells_in_this_dimension << " edges with " << number_of_attributes << " attributes." << std::endl;
        reading_header = false;
        cells.push_back(std::map<long, typename LCCRW::LCC::Dart_handle>());
        Index_builder<typename LCCRW::LCCB, 1>::build_for_dimension(lccrw.builder);
      }
      
      else {
        line_stream >> current_word;
        current_index = string_to_any(current_word, current_number);
        
        line_stream >> current_word;
        current_number = string_to_any(current_word, current_number);
        typename LCCRW::LCC::Dart_handle start_vertex = cells[0][current_number];
        
        line_stream >> current_word;
        current_number = string_to_any(current_word, current_number);
        typename LCCRW::LCC::Dart_handle end_vertex = cells[0][current_number];
        
        typename LCCRW::LCC::template Attribute_handle<1>::type edge_attribute = lccrw.lcc.template create_attribute<1>(-1);
        //std::cout << "Created edge<" << &*edge_attribute << "> (" << std::endl;
        cells[1][current_index] = lccrw.builder.get_edge(edge_attribute, start_vertex, end_vertex).first;
        cells[1][current_index]->template attribute<1>()->info() = current_index;
        std::cout << "\tedge<" << &*cells[1][current_index]->template attribute<1>() << ">[" << cells[1][current_index]->template attribute<1>()->info() << "]" << std::endl;
        --remaining_cells_in_this_dimension;
      }
      
      if (remaining_cells_in_this_dimension == 0) Cell_complex_parser<LCCRW, 2, typename LCCRW::LCC::template Attribute_type<2>::type>::parse(lccrw, input_file_stream, cells);
    }
  }
};

template <class LCCRW>
struct Cell_complex_parser<LCCRW, 2, typename LCCRW::LCC::template Attribute_type<2>::type> {
  static void parse(LCCRW &lccrw, std::ifstream &input_file_stream, std::vector<std::map<long, typename LCCRW::LCC::Dart_handle> > &cells) {
    std::string current_line, current_word;
    unsigned long number_of_attributes, current_index;
    std::istringstream line_stream;
    unsigned long current_number, remaining_cells_in_this_dimension;
    bool reading_header = true;
    
    while (!input_file_stream.eof()) {
      getline(input_file_stream, current_line);
      std::cout << "Line: \'" << current_line << "\'" << std::endl;
      
      // Empty line
      if (current_line.size() == 0) {
        continue;
      }
      
      // Just a newline sequence
      if (current_line[0] == '\r' || current_line[0] == '\n') {
        continue;
      }
      
      line_stream = std::istringstream(current_line + "\n");
      
      if (reading_header) {
        line_stream >> current_word;
        current_number = string_to_any(current_word, current_number);
        remaining_cells_in_this_dimension = current_number;
        
        line_stream >> current_word;
        current_number = string_to_any(current_word, current_number);
        number_of_attributes = current_number;
        
        std::cout << remaining_cells_in_this_dimension << " facets with " << number_of_attributes << " attributes." << std::endl;
        reading_header = false;
        cells.push_back(std::map<long, typename LCCRW::LCC::Dart_handle>());
        Index_builder<typename LCCRW::LCCB, 2>::build_for_dimension(lccrw.builder);
      }
      
      else {
        
        // If there are no edges, assume the indices refer to vertices
        if (cells[1].size() == 0) {
          line_stream >> current_word;
          current_index = string_to_any(current_word, current_index);
          
          std::list<typename LCCRW::LCC::Dart_handle> vertices_in_facet;
          while (line_stream.peek() != '\n') {
            line_stream >> current_word;
            current_number = string_to_any(current_word, current_number);
            vertices_in_facet.push_back(cells[0][current_number]);
          }
          
          typename LCCRW::LCC::template Attribute_handle<2>::type facet_attribute = lccrw.lcc.template create_attribute<2>(-1);
          //std::cout << "Created facet<" << &*facet_attribute << "> (" << std::endl;
          cells[2][current_index] = lccrw.builder.get_facet_from_vertices(facet_attribute, vertices_in_facet.begin(), vertices_in_facet.end()).first;
          cells[2][current_index]->template attribute<2>()->info() = current_index;
          std::cout << "\tfacet<" << &*cells[2][current_index]->template attribute<2>() << ">[" << cells[2][current_index]->template attribute<2>()->info() << "]" << std::endl;
        }
        
        // TODO: Build from edges
        else {
          
        }
        
        --remaining_cells_in_this_dimension;
      }
      
      if (remaining_cells_in_this_dimension == 0) Cell_complex_parser<LCCRW, 3, typename LCCRW::LCC::template Attribute_type<3>::type>::parse(lccrw, input_file_stream, cells);
    }
  }
};

template <class LCC, unsigned int dimension = LCC::dimension>
struct Ply_writer {
  static void write_ply(LCC &lcc, const char *path) {
    std::cerr << "Error: Cannot write non 2D/3D LCC!" << std::endl;
  }
};

template <class LCC>
struct Ply_writer<LCC, 2> {
  
  static void write_ply(LCC &lcc, const char *path) {
    CGAL_assertion(LCC::dimension >= 2 && LCC::dimension <= 3);
    
    std::cout << "Writing file " << path << "..." << std::endl;
    std::ofstream output_file_stream;
    output_file_stream.open(path);
    if (!output_file_stream.is_open()) {
      std::cerr << "Could not create file." << std::endl;
      return;
    }
    
    output_file_stream << "ply" << std::endl;
    output_file_stream << "format ascii 1.0" << std::endl;
    output_file_stream << "comment created by ken" << std::endl;
    
    typename LCC::Vertex_attribute_const_range &vertices = lcc.vertex_attributes();
    output_file_stream << "element vertex " << vertices.size() << std::endl;
    output_file_stream << "property float x" << std::endl;
    output_file_stream << "property float y" << std::endl;
    output_file_stream << "property float z" << std::endl;
    
    typename LCC::template Attribute_const_range<1>::type &edges = lcc.template attributes<1>();
    output_file_stream << "element edge " << edges.size() << std::endl;
    output_file_stream << "property int vertex1" << std::endl;
    output_file_stream << "property int vertex2" << std::endl;
    
    typename LCC::template Attribute_const_range<2>::type &facets = lcc.template attributes<2>();
    output_file_stream << "element face " << facets.size() << std::endl;
    output_file_stream << "property list uchar int vertex_index" << std::endl;
    
    output_file_stream << "end_header" << std::endl;
    
    std::map<typename LCC::Vertex_attribute_const_handle, unsigned int> vertex_map;
    for (typename LCC::Vertex_attribute_const_range::const_iterator current_vertex = vertices.begin(); current_vertex != vertices.end(); ++current_vertex) {
      output_file_stream << current_vertex->point().x() << " " << current_vertex->point().y() << " 0" << std::endl;
      vertex_map[current_vertex] = vertex_map.size();
    }
    
    for (typename LCC::template Attribute_const_range<1>::type::const_iterator current_edge = edges.begin(); current_edge != edges.end(); ++current_edge) {
      typename LCC::Dart_const_handle dart_in_edge = current_edge->dart();
      output_file_stream << vertex_map[dart_in_edge->template attribute<0>()] << " " << vertex_map[dart_in_edge->beta(1)->template attribute<0>()] << std::endl;
    }
    
    for (typename LCC::template Attribute_const_range<2>::type::const_iterator current_facet = facets.begin(); current_facet != facets.end(); ++current_facet) {
      typename LCC::Dart_const_handle dart_in_facet = current_facet->dart();
      typename LCC::template Dart_of_cell_const_range<2, 2> darts_in_facet = lcc.template darts_of_cell<2, 2>(dart_in_facet);
      output_file_stream << darts_in_facet.size();
      for (typename LCC::template Dart_of_cell_const_range<2, 2>::const_iterator current_dart = darts_in_facet.begin(); current_dart != darts_in_facet.end(); ++current_dart) {
        output_file_stream << " " << vertex_map[current_dart->template attribute<0>()];
      } output_file_stream << std::endl;
    }
    
    output_file_stream.close();
  }
};

template <class LCC>
struct Ply_writer<LCC, 3> {
  
  static void write_ply(LCC &lcc, const char *path) {
    CGAL_assertion(LCC::dimension >= 2 && LCC::dimension <= 3);
    
    std::cout << "Writing file " << path << "..." << std::endl;
    std::ofstream output_file_stream;
    output_file_stream.open(path);
    if (!output_file_stream.is_open()) {
      std::cerr << "Could not create file." << std::endl;
      return;
    }
    
    output_file_stream << "ply" << std::endl;
    output_file_stream << "format ascii 1.0" << std::endl;
    output_file_stream << "comment created by ken" << std::endl;
    
    typename LCC::Vertex_attribute_const_range &vertices = lcc.vertex_attributes();
    output_file_stream << "element vertex " << vertices.size() << std::endl;
    output_file_stream << "property float x" << std::endl;
    output_file_stream << "property float y" << std::endl;
    output_file_stream << "property float z" << std::endl;
    
    typename LCC::template Attribute_const_range<1>::type &edges = lcc.template attributes<1>();
    output_file_stream << "element edge " << edges.size() << std::endl;
    output_file_stream << "property int vertex1" << std::endl;
    output_file_stream << "property int vertex2" << std::endl;
    
    typename LCC::template Attribute_const_range<2>::type &facets = lcc.template attributes<2>();
    output_file_stream << "element face " << facets.size() << std::endl;
    output_file_stream << "property list uchar int vertex_index" << std::endl;
    
    output_file_stream << "end_header" << std::endl;
    
    std::map<typename LCC::Vertex_attribute_const_handle, unsigned int> vertex_map;
    for (typename LCC::Vertex_attribute_const_range::const_iterator current_vertex = vertices.begin(); current_vertex != vertices.end(); ++current_vertex) {
      output_file_stream << current_vertex->point().x() << " " << current_vertex->point().y() << " " << current_vertex->point().z() << std::endl;
      vertex_map[current_vertex] = vertex_map.size();
    }
    
    for (typename LCC::template Attribute_const_range<1>::type::const_iterator current_edge = edges.begin(); current_edge != edges.end(); ++current_edge) {
      typename LCC::Dart_const_handle dart_in_edge = current_edge->dart();
      output_file_stream << vertex_map[dart_in_edge->template attribute<0>()] << " " << vertex_map[dart_in_edge->beta(1)->template attribute<0>()] << std::endl;
    }
    
    for (typename LCC::template Attribute_const_range<2>::type::const_iterator current_facet = facets.begin(); current_facet != facets.end(); ++current_facet) {
      typename LCC::Dart_const_handle dart_in_facet = current_facet->dart();
      typename LCC::template Dart_of_cell_const_range<2, 2> darts_in_facet = lcc.template darts_of_cell<2, 2>(dart_in_facet);
      output_file_stream << darts_in_facet.size();
      for (typename LCC::template Dart_of_cell_const_range<2, 2>::const_iterator current_dart = darts_in_facet.begin(); current_dart != darts_in_facet.end(); ++current_dart) {
        output_file_stream << " " << vertex_map[current_dart->template attribute<0>()];
      } output_file_stream << std::endl;
    }
    
    output_file_stream.close();
  }
};

template <class LCCRW>
struct OGR_reader {
  static void load(LCCRW &lccrw, const char *path) {
    std::cout << "Reading " << path << "..." << std::endl;
    
    OGRRegisterAll();
    OGRDataSource *data_source = OGRSFDriverRegistrar::Open(path, false);
    if (data_source == NULL) {
      std::cerr << "Could not open file." << std::endl;
      return;
    }
    
    std::cout << "\tType: " << data_source->GetDriver()->GetName() << std::endl;
    std::cout << "\tLayers: " << data_source->GetLayerCount() << std::endl;
    
    Index_builder<typename LCCRW::LCCB, 0>::build_for_dimension(lccrw.builder);
    Index_builder<typename LCCRW::LCCB, 2>::build_for_dimension(lccrw.builder);
    std::list<typename LCCRW::LCC::Dart_handle> facets;
    
    for (int current_layer = 0; current_layer < data_source->GetLayerCount(); ++current_layer) {
      OGRLayer *layer = data_source->GetLayer(current_layer);
      layer->ResetReading();
      std::cout << "\tReading layer " << current_layer << " (" << layer->GetFeatureCount(true) << " features)..." << std::endl;
      
      OGRFeature *feature;
      while ((feature = layer->GetNextFeature()) != NULL) {
        switch(feature->GetGeometryRef()->getGeometryType()) {
            
          case wkbPolygon: {
            OGRPolygon *geometry = static_cast<OGRPolygon *>(feature->GetGeometryRef());
            
            // Outer ring
            std::list<typename LCCRW::LCC::Dart_handle> vertices;
            for (int current_point = 0; current_point < geometry->getExteriorRing()->getNumPoints(); ++current_point) {
              Point_creator<typename LCCRW::LCC::Point, typename LCCRW::LCC::FT> creator;
              creator.push_back(geometry->getExteriorRing()->getX(current_point));
              creator.push_back(geometry->getExteriorRing()->getY(current_point));
              typename LCCRW::LCC::Vertex_attribute_handle vertex_attribute = lccrw.lcc.template create_attribute<0>(creator.create_point(), lccrw.lcc.template attributes<0>().size());
              vertices.push_back(lccrw.builder.get_vertex(vertex_attribute));
            }
            
            // TODO: Inner rings
            
            if (vertices.front() == vertices.back()) vertices.pop_back();
            typename LCCRW::LCC::template Attribute_handle<2>::type facet_attribute = lccrw.lcc.template create_attribute<2>(feature->GetFID());
            facets.push_back(lccrw.builder.get_facet_from_vertices(facet_attribute, vertices.begin(), vertices.end(), false).first);
            
            break;
          }
            
          default:
            std::cerr << "Feature type not supported." << std::endl;
            break;
        }
        
        OGRFeature::DestroyFeature(feature);
      }
    }
    
    lccrw.builder.template link_cells<2>(facets.begin(), facets.end());
    //CGAL_postcondition(lccrw.builder.lcc.is_valid());
    
    OGRDataSource::DestroyDataSource(data_source);
  }
};

template <class LCC, unsigned int dimension = 0, class Attribute_type = typename LCC::Vertex_attribute>
struct Map_attributes_writer;

template <class LCC, unsigned int dimension, class Attribute_type>
struct Map_attributes_writer {
  static void write_attributes_as_gmap(LCC &lcc, std::ofstream &output_file_stream, std::map<typename LCC::Dart_const_handle, unsigned int> &dart_handle_to_dart_index) {
    //std::map<typename LCCRW::LCC::template Attribute_const_handle<dimension>, unsigned int> attribute_handle_to_attribute_index;
    //unsigned int current_attribute_index = 0;
    
    output_file_stream << "===" << dimension << "-cells===" << std::endl;
    output_file_stream << "address id dart" << std::endl;
    
    for (typename LCC::template Attribute_const_range<dimension>::type::const_iterator current_attribute = lcc.template attributes<dimension>().begin();
         current_attribute != lcc.template attributes<dimension>().end();
         ++current_attribute) {
      output_file_stream << &*current_attribute << " " << current_attribute->info() << " " << dart_handle_to_dart_index[current_attribute->dart()] << std::endl;
    }
    
    Map_attributes_writer<LCC, dimension+1, typename LCC::template Attribute_type<dimension+1>::type>::write_attributes_as_gmap(lcc, output_file_stream, dart_handle_to_dart_index);
  }
  
  static void write_attributes_as_cmap(LCC &lcc, std::ofstream &output_file_stream) {
    
    output_file_stream << "===" << dimension << "-cells===" << std::endl;
    output_file_stream << "address id dart" << std::endl;
    
    for (typename LCC::template Attribute_const_range<dimension>::type::const_iterator current_attribute = lcc.template attributes<dimension>().begin();
         current_attribute != lcc.template attributes<dimension>().end();
         ++current_attribute) {
      output_file_stream << &*current_attribute << " " << current_attribute->info() << " " << current_attribute->dart()->id << std::endl;
    }
    
    Map_attributes_writer<LCC, dimension+1, typename LCC::template Attribute_type<dimension+1>::type>::write_attributes_as_cmap(lcc, output_file_stream);
  }
};

template <class LCC>
struct Map_attributes_writer<LCC, 0, typename LCC::Vertex_attribute> {
  static void write_attributes_as_gmap(LCC &lcc, std::ofstream &output_file_stream, std::map<typename LCC::Dart_const_handle, unsigned int> &dart_handle_to_dart_index) {
    //std::map<typename LCCRW::LCC::Vertex_attribute, unsigned int> attribute_handle_to_attribute_index;
    //unsigned int current_attribute_index = 0;
    
    output_file_stream << "===0-cells===" << std::endl;
    output_file_stream << "address id coords dart" << std::endl;
    
    for (typename LCC::Vertex_attribute_const_range::const_iterator current_attribute = lcc.vertex_attributes().begin();
         current_attribute != lcc.vertex_attributes().end();
         ++current_attribute) {
      //output_file_stream << current_attribute_index << " " << current_attribute->point() << std::endl;
      output_file_stream << &*current_attribute << " " << current_attribute->info() << " " << current_attribute->point() << " " << dart_handle_to_dart_index[current_attribute->dart()] << std::endl;
      //attribute_handle_to_attribute_index[current_attribute] = current_attribute_index;
      //++current_attribute_index;
    }
    
    Map_attributes_writer<LCC, 1, typename LCC::template Attribute_type<1>::type>::write_attributes_as_gmap(lcc, output_file_stream, dart_handle_to_dart_index);
  }
  
  static void write_attributes_as_cmap(LCC &lcc, std::ofstream &output_file_stream) {
    
    output_file_stream << "===0-cells===" << std::endl;
    output_file_stream << "address id coords dart" << std::endl;
    
    for (typename LCC::Vertex_attribute_const_range::const_iterator current_attribute = lcc.vertex_attributes().begin();
         current_attribute != lcc.vertex_attributes().end();
         ++current_attribute) {
      output_file_stream << &*current_attribute << " " << current_attribute->info() << " " << current_attribute->point() << " " << current_attribute->dart()->id << std::endl;
    }
    
    Map_attributes_writer<LCC, 1, typename LCC::template Attribute_type<1>::type>::write_attributes_as_cmap(lcc, output_file_stream);
  }
};

template <class LCC, unsigned int dimension>
struct Map_attributes_writer<LCC, dimension, CGAL::Void> {
  static void write_attributes_as_gmap(LCC &lcc, std::ofstream &output_file_stream, std::map<typename LCC::Dart_const_handle, unsigned int> &dart_handle_to_dart_index) {
    
  }
  
  static void write_attributes_as_cmap(LCC &lcc, std::ofstream &output_file_stream) {
    
  }
};

template <class LCC, unsigned int dimension = LCC::dimension>
struct Map_dart_link_writer;

template <class LCC, unsigned int dimension>
struct Map_dart_link_writer {
  static void write_as_gmap_origin(LCC &lcc, std::map<typename LCC::Dart_const_handle, unsigned int> &dart_handle_to_dart_index, typename LCC::Dart_const_handle &dart, std::ofstream &output_file_stream) {
    Map_dart_link_writer<LCC, dimension-1>::write_as_gmap_origin(lcc, dart_handle_to_dart_index, dart, output_file_stream);
    if (!lcc.is_free(dart, dimension)) output_file_stream << " " << dart_handle_to_dart_index[dart->beta(dimension)]+1;
    else output_file_stream << " -1";
  }
  
  static void write_as_gmap_destination(LCC &lcc, std::map<typename LCC::Dart_const_handle, unsigned int> &dart_handle_to_dart_index, typename LCC::Dart_const_handle &dart, std::ofstream &output_file_stream) {
    Map_dart_link_writer<LCC, dimension-1>::write_as_gmap_destination(lcc, dart_handle_to_dart_index, dart, output_file_stream);
    if (!lcc.is_free(dart, dimension)) output_file_stream << " " << dart_handle_to_dart_index[dart->beta(dimension)];
    else output_file_stream << " -1";
  }
  
  static void write_as_cmap(LCC &lcc, typename LCC::Dart_const_handle &dart, std::ofstream &output_file_stream) {
    Map_dart_link_writer<LCC, dimension-1>::write_as_cmap(lcc, dart, output_file_stream);
    if (!lcc.is_free(dart, dimension)) output_file_stream << " " << dart->beta(dimension)->id;
    else output_file_stream << " -1";
  }
};

template <class LCC>
struct Map_dart_link_writer<LCC, 0> {
  static void write_as_gmap_origin(LCC &lcc, std::map<typename LCC::Dart_const_handle, unsigned int> &dart_handle_to_dart_index, typename LCC::Dart_const_handle &dart, std::ofstream &output_file_stream) {
    output_file_stream << dart_handle_to_dart_index[dart] << " " << dart_handle_to_dart_index[dart]+1;
  }
  
  static void write_as_gmap_destination(LCC &lcc, std::map<typename LCC::Dart_const_handle, unsigned int> &dart_handle_to_dart_index, typename LCC::Dart_const_handle &dart, std::ofstream &output_file_stream) {
    output_file_stream << dart_handle_to_dart_index[dart]+1 << " " << dart_handle_to_dart_index[dart];
  }
  
  static void write_as_cmap(LCC &lcc, typename LCC::Dart_const_handle &dart, std::ofstream &output_file_stream) {
    output_file_stream << dart->id;
    if (!lcc.is_free(dart, 0)) output_file_stream << " " << dart->beta(0)->id;
    else output_file_stream << " -1";
  }
};

template <class LCC, unsigned int dimension = LCC::dimension>
struct Map_dart_attribute_writer;

template <class LCC, unsigned int dimension>
struct Map_dart_attribute_writer {
  static void write_as_gmap_origin(std::map<typename LCC::Dart_const_handle, unsigned int> &dart_handle_to_dart_index, typename LCC::Dart_const_handle &dart, std::ofstream &output_file_stream) {
    Map_dart_attribute_writer<LCC, dimension-1>::write_as_gmap_origin(dart_handle_to_dart_index, dart, output_file_stream);
    output_file_stream << " " << &*dart->template attribute<dimension>();
  }
  
  static void write_as_gmap_destination(std::map<typename LCC::Dart_const_handle, unsigned int> &dart_handle_to_dart_index, typename LCC::Dart_const_handle &dart, std::ofstream &output_file_stream) {
    Map_dart_attribute_writer<LCC, dimension-1>::write_as_gmap_destination(dart_handle_to_dart_index, dart, output_file_stream);
    output_file_stream << " " << &*dart->template attribute<dimension>();
  }
  
  static void write_as_cmap(typename LCC::Dart_const_handle &dart, std::ofstream &output_file_stream) {
    Map_dart_attribute_writer<LCC, dimension-1>::write_as_cmap(dart, output_file_stream);
    output_file_stream << " " << &*dart->template attribute<dimension>();
  }
};

template <class LCC>
struct Map_dart_attribute_writer<LCC, 0> {
  static void write_as_gmap_origin(std::map<typename LCC::Dart_const_handle, unsigned int> &dart_handle_to_dart_index, typename LCC::Dart_const_handle &dart, std::ofstream &output_file_stream) {
    output_file_stream << " " << &*dart->template attribute<0>();
  }
  
  static void write_as_gmap_destination(std::map<typename LCC::Dart_const_handle, unsigned int> &dart_handle_to_dart_index, typename LCC::Dart_const_handle &dart, std::ofstream &output_file_stream) {
    output_file_stream << " " << &*dart->beta(1)->template attribute<0>();
  }
  
  static void write_as_cmap(typename LCC::Dart_const_handle &dart, std::ofstream &output_file_stream) {
    output_file_stream << " " << &*dart->template attribute<0>();
  }
};

template <class LCC>
struct Map_darts_writer {
  static void write_darts_as_gmap(LCC &lcc, std::ofstream &output_file_stream, std::map<typename LCC::Dart_const_handle, unsigned int> &dart_handle_to_dart_index) {
    output_file_stream << "===Darts===" << std::endl;
    
    for (typename LCC::Dart_const_range::const_iterator current_dart = lcc.darts().begin();
         current_dart != lcc.darts().end();
         ++current_dart) {
      Map_dart_link_writer<LCC>::write_as_gmap_origin(lcc, dart_handle_to_dart_index, current_dart, output_file_stream);
      Map_dart_attribute_writer<LCC>::write_as_gmap_origin(dart_handle_to_dart_index, current_dart, output_file_stream);
      output_file_stream << std::endl;
      Map_dart_link_writer<LCC>::write_as_gmap_destination(lcc, dart_handle_to_dart_index, current_dart, output_file_stream);
      Map_dart_attribute_writer<LCC>::write_as_gmap_destination(dart_handle_to_dart_index, current_dart, output_file_stream);
      output_file_stream << std::endl;
    }
  }
  
  static void write_darts_as_cmap(LCC &lcc, std::ofstream &output_file_stream) {
    output_file_stream << "===Darts===" << std::endl;
    
    for (typename LCC::Dart_const_range::const_iterator current_dart = lcc.darts().begin();
         current_dart != lcc.darts().end();
         ++current_dart) {
      Map_dart_link_writer<LCC>::write_as_cmap(lcc, current_dart, output_file_stream);
      Map_dart_attribute_writer<LCC>::write_as_cmap(current_dart, output_file_stream);
      output_file_stream << std::endl;
    }
  }
};

template <class LCC>
struct Map_writer {
  static void write_as_gmap(LCC &lcc, const char *path) {
    //CGAL_precondition(lcc.is_valid());
    std::cout << "Writing file " << path << "..." << std::endl;
    std::ofstream output_file_stream;
    output_file_stream.open(path);
    if (!output_file_stream.is_open()) {
      std::cerr << "Could not create file." << std::endl;
      return;
    }
    
    std::map<typename LCC::Dart_const_handle, unsigned int> dart_handle_to_dart_index;
    unsigned int current_dart_index = 0;
    
    for (typename LCC::Dart_const_range::const_iterator current_dart = lcc.darts().begin();
         current_dart != lcc.darts().end();
         ++current_dart) {
      if (dart_handle_to_dart_index.count(current_dart) == 0) {
        dart_handle_to_dart_index[current_dart] = current_dart_index;
        current_dart_index += 2;
      }
    }
    
    output_file_stream << "===LCC info===" << std::endl;
    output_file_stream << "dimension: " << LCC::dimension << std::endl;
    output_file_stream << "embedding dimension: " << LCC::ambient_dimension << std::endl;
    
    Map_attributes_writer<LCC, 0>::write_attributes_as_gmap(lcc, output_file_stream, dart_handle_to_dart_index);
    Map_darts_writer<LCC>::write_darts_as_gmap(lcc, output_file_stream, dart_handle_to_dart_index);
  }
  
  static void write_as_cmap(LCC &lcc, const char *path) {
    //CGAL_precondition(lcc.is_valid());
    std::cout << "Writing file " << path << "..." << std::endl;
    std::ofstream output_file_stream;
    output_file_stream.open(path);
    if (!output_file_stream.is_open()) {
      std::cerr << "Could not create file." << std::endl;
      return;
    }
    
    output_file_stream << "===LCC info===" << std::endl;
    output_file_stream << "Dimension=" << LCC::dimension << ", Ambient dimension=" << LCC::ambient_dimension << std::endl;
    lcc.display_characteristics(output_file_stream);
    output_file_stream << std::endl;
    
    Map_attributes_writer<LCC, 0>::write_attributes_as_cmap(lcc, output_file_stream);
    Map_darts_writer<LCC>::write_darts_as_cmap(lcc, output_file_stream);
  }
};

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

template <unsigned int d>
class Linear_cell_complex_reader_writer {
public:
  typedef typename Linear_cell_complex_with_ids<d>::type LCC;
  typedef typename LCC::Point Point;
  typedef typename LCC::Traits::FT FT;
  typedef Linear_cell_complex_incremental_builder<LCC> LCCB;
  typedef Linear_cell_complex_reader_writer<d> Self;
  static const unsigned int dimension = d;
  
  LCC lcc;
  LCCB builder;
  
  Linear_cell_complex_reader_writer() : builder(lcc) {
    
  }
  
  void open_cell_complex(const char *path) {
    std::cout << "Opening file " << path << "..." << std::endl;
    std::ifstream input_file_stream;
    input_file_stream.open(path);
    if (!input_file_stream.is_open()) {
      std::cerr << "Could not open file." << std::endl;
      return;
    }
    
    std::vector<std::map<long, typename LCC::Dart_handle> > cells;
    Cell_complex_parser<Self>::parse(*this, input_file_stream, cells);
    std::list<typename LCC::Dart_handle> highest_dimensional_cells;
    for (typename std::map<long, typename LCC::Dart_handle>::iterator current_cell = cells[LCC::dimension].begin();
         current_cell != cells[LCC::dimension].end();
         ++current_cell)
      highest_dimensional_cells.push_back(current_cell->second);
    builder.template link_cells<LCC::dimension>(highest_dimensional_cells.begin(), highest_dimensional_cells.end());
    input_file_stream.close();
  }
  
  void print_cell_complex() {
    std::cout << "========== Cell complex ==========" << std::endl;
    Cell_complex_printer<LCC>::print(lcc);
  }
  
  void print_facets() {
    std::cout << "========== Facets ==========" << std::endl;
    typename LCC::template Attribute_range<2>::type &facets = lcc.template attributes<2>();
    for (typename LCC::template Attribute_range<2>::type::iterator current_facet = facets.begin(); current_facet != facets.end(); ++current_facet) {
      std::cout << "facet<" << &*current_facet << ">[" << current_facet->info() << "]:" << std::endl;
      typename LCC::Dart_handle first_dart = current_facet->dart();
      typename LCC::Dart_handle current_dart = first_dart;
      do {
        std::cout << "\tdart<" << &*current_dart << ">[" << current_dart->id << "] (" << current_dart->template attribute<0>()->point() << ")" << std::endl;
        current_dart = current_dart->beta(1);
      } while (current_dart != first_dart);
    }
  }
  
  void write_ply(const char *path) {
    Ply_writer<LCC>::write_ply(lcc, path);
  }
  
  void load_from_ogr(const char *path) {
    OGR_reader<Self>::load(*this, path);
  }
  
  void write_gmap(const char *path) {
    Map_writer<LCC>::write_as_gmap(lcc, path);
  }
  
  static void write_gmap(const char *path, LCC &l) {
    Map_writer<LCC>::write_as_gmap(l, path);
  }
  
  void write_cmap(const char *path) {
    Map_writer<LCC>::write_as_cmap(lcc, path);
  }
  
  static void write_cmap(const char *path, LCC &l) {
    Map_writer<LCC>::write_as_cmap(l, path);
  }
  
  void write_obj(const char *path) {
    Obj_writer<LCC>::write_obj(lcc, path);
  }
  
  static void write_obj(const char *path, LCC &l) {
    Obj_writer<LCC>::write_obj(l, path);
  }
  
  void write_darts_as_obj(const char *path) {
    Obj_writer<LCC>::write_darts_as_obj(lcc, path);
  }
  
  static void write_darts_as_obj(const char *path, LCC &l) {
    Obj_writer<LCC>::write_darts_as_obj(l, path);
  }
};

#endif
