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

#include <fstream>

#include <gdal/ogrsf_frmts.h>

#include "Linear_cell_complex_builder.h"
#include "Linear_cell_complex_builder_with_index.h"
#include "Linear_cell_complex_printer.h"
#include "Linear_cell_complex_reader_writer.h"

void load_Martijns_format(const char *path) {
  typedef Linear_cell_complex_with_ids<3> LCC;
  //  typedef CGAL::Linear_cell_complex<3> LCC;
  LCC lcc;
  LCC lcci;
  Linear_cell_complex_builder<LCC> builder(lcc);
  Linear_cell_complex_builder_with_index<LCC> builder_with_index(lcci);
  
  std::ifstream input_file_stream(path);
  std::string current_line;
  
  std::map<int, std::list<OGRGeometry *> > polyhedra;
  
  // Parse
  while (!input_file_stream.eof()) {
    getline(input_file_stream, current_line);
    //    std::cout << "Line: \'" << current_line << "\'" << std::endl;
    
    // Empty line
    if (current_line.size() == 0) {
      continue;
    }
    
    // Just a newline sequence
    if (current_line[0] == '\r' || current_line[0] == '\n') {
      continue;
    }
    
    // Comment
    if (current_line[0] == '#') {
      continue;
    }
    
    std::size_t separator = current_line.find(" ; ");
    if (separator == current_line.npos) {
      continue;
    }
    
    int polyhedron_id = std::stoi(current_line.substr(0, separator));
    std::string polygon_wkt = current_line.substr(separator+3);
    char *polygon_cstr = new char[polygon_wkt.length()+1];
    std::strcpy(polygon_cstr, polygon_wkt.c_str());
    //    std::cout << "\tID   : \"" << polyhedron_id << "\"" << std::endl;
    //    std::cout << "\tGeom : \"" << polygon_wkt << "\"" << std::endl;
    
    OGRGeometry *polygon;
    char *polygon_cstr_pointer = polygon_cstr;
    OGRErr error = OGRGeometryFactory::createFromWkt(&polygon_cstr_pointer, NULL, &polygon);
    if (error == OGRERR_NONE) {
      polyhedra[polyhedron_id].push_back(polygon);
      //      std::cout << "\tSuccessfully created OGRGeometry<" << polygon << ">" << std::endl;
    }
    delete[] polygon_cstr;
  }
  
  // Create lcc
  std::list<LCC::Dart_handle> volumes, volumes_with_index;
  std::clock_t start_time;
  for (std::map<int, std::list<OGRGeometry *> >::iterator current_polyhedron = polyhedra.begin(); current_polyhedron != polyhedra.end(); ++current_polyhedron) {
    std::cout << "Polyhedron " << current_polyhedron->first << " has " << current_polyhedron->second.size() << " faces." << std::endl;
    std::list<LCC::Dart_handle> faces, faces_with_index;
    std::clock_t get_vertex_total_time = 0, get_face_total_time = 0, validate_face_total_time = 0, validate_volume_total_time = 0, get_cell_total_time;
    std::clock_t get_vertex_with_index_total_time = 0, get_face_with_index_total_time = 0, validate_face_with_index_total_time = 0, validate_volume_with_index_total_time = 0, get_cell_with_index_total_time;
    
    // Each face
    for (std::list<OGRGeometry *>::iterator current_polygon = current_polyhedron->second.begin(); current_polygon != current_polyhedron->second.end(); ++current_polygon) {
      if ((*current_polygon)->getGeometryType() != wkbPolygon25D) {
        std::cerr << "\tIncorrect geometry type" << std::endl;
        continue;
      } OGRPolygon *polygon = static_cast<OGRPolygon *>(*current_polygon);
      //      std::cout << "\tPolygon with " << polygon->getExteriorRing()->getNumPoints()-1 << " vertices" << std::endl;
      
      // Vertices
      start_time = clock();
      std::vector<LCC::Dart_handle> vertices;
      vertices.reserve(polygon->getExteriorRing()->getNumPoints());
      for (int current_vertex = 0; current_vertex < polygon->getExteriorRing()->getNumPoints(); ++current_vertex) {
        LCC::Point vertex(polygon->getExteriorRing()->getX(current_vertex),
                          polygon->getExteriorRing()->getY(current_vertex),
                          polygon->getExteriorRing()->getZ(current_vertex));
        vertices.push_back(builder.get_vertex(vertex));
      } get_vertex_total_time += clock() - start_time;
      
      // Vertices with index
      start_time = clock();
      std::vector<LCC::Dart_handle> vertices_with_index;
      vertices_with_index.reserve(polygon->getExteriorRing()->getNumPoints());
      for (int current_vertex = 0; current_vertex < polygon->getExteriorRing()->getNumPoints(); ++current_vertex) {
        LCC::Point vertex(polygon->getExteriorRing()->getX(current_vertex),
                          polygon->getExteriorRing()->getY(current_vertex),
                          polygon->getExteriorRing()->getZ(current_vertex));
        vertices_with_index.push_back(builder_with_index.get_vertex(vertex));
      } get_vertex_with_index_total_time += clock() - start_time;
      
      // Face
      start_time = clock();
      faces.push_back(builder.get_face(vertices.begin(), vertices.end()));
      get_face_total_time += clock() - start_time;
      
      // Face with index
      start_time = clock();
      faces_with_index.push_back(builder_with_index.get_face(vertices_with_index.begin(), vertices_with_index.end()));
      get_face_with_index_total_time += clock() - start_time;
      
      // Validate that faces created using both methods are equal
      CGAL_assertion(Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry<2>(lcc, faces.back(), lcci, faces_with_index.back()));
      
      // Validate face creation
      start_time = clock();
      LCC::Dart_handle start_dart = faces.back();
      while (lcc.point(start_dart) != lcc.point(vertices[0])) start_dart = start_dart->beta(1);
      LCC::Dart_handle current_dart = start_dart;
      bool same_orientation;
      if (lcc.point(vertices[1]) == lcc.point(current_dart->beta(1))) same_orientation = true;
      else same_orientation = false;
      unsigned int vertices_in_face = 0;
      do {
        CGAL_postcondition(lcc.point(current_dart) == lcc.point(vertices[vertices_in_face]));
        CGAL_postcondition(current_dart->beta(0) != current_dart);
        CGAL_postcondition(current_dart->beta(1) != current_dart);
        CGAL_postcondition(current_dart->beta(0)->beta(1) == current_dart);
        CGAL_postcondition(current_dart->beta(1)->beta(0) == current_dart);
        if (same_orientation) current_dart = current_dart->beta(1);
        else current_dart = current_dart->beta(0);
        ++vertices_in_face;
      } while (current_dart != start_dart);
      CGAL_postcondition(vertices_in_face == polygon->getExteriorRing()->getNumPoints()-1);
      vertices_in_face = 0;
      for (LCC::template Dart_of_cell_range<2, 2>::iterator current_dart_in_cell = lcc.template darts_of_cell<2, 2>(start_dart).begin(); current_dart_in_cell != lcc.template darts_of_cell<2, 2>(start_dart).end(); ++current_dart_in_cell) {
        if (same_orientation) CGAL_postcondition(lcc.point(current_dart_in_cell) == lcc.point(vertices[vertices_in_face]));
        else CGAL_postcondition(lcc.point(current_dart_in_cell) == lcc.point(vertices[vertices.size()-vertices_in_face-1]));
        ++vertices_in_face;
      } CGAL_postcondition(vertices_in_face == polygon->getExteriorRing()->getNumPoints()-1);
      validate_face_total_time += clock() - start_time;
      std::cout << ".";
      
      // Validate face creation with index
      start_time = clock();
      LCC::Dart_handle start_dart_with_index = faces_with_index.back();
      while (lcci.point(start_dart_with_index) != lcci.point(vertices_with_index[0])) start_dart_with_index = start_dart_with_index->beta(1);
      LCC::Dart_handle current_dart_with_index = start_dart_with_index;
      bool same_orientation_with_index;
      if (lcci.point(vertices_with_index[1]) == lcci.point(current_dart_with_index->beta(1))) same_orientation_with_index = true;
      else same_orientation_with_index = false;
      unsigned int vertices_in_face_with_index = 0;
      do {
        CGAL_postcondition(lcci.point(current_dart_with_index) == lcci.point(vertices_with_index[vertices_in_face_with_index]));
        CGAL_postcondition(current_dart_with_index->beta(0) != current_dart_with_index);
        CGAL_postcondition(current_dart_with_index->beta(1) != current_dart_with_index);
        CGAL_postcondition(current_dart_with_index->beta(0)->beta(1) == current_dart_with_index);
        CGAL_postcondition(current_dart_with_index->beta(1)->beta(0) == current_dart_with_index);
        if (same_orientation_with_index) current_dart_with_index = current_dart_with_index->beta(1);
        else current_dart_with_index = current_dart_with_index->beta(0);
        ++vertices_in_face_with_index;
      } while (current_dart_with_index != start_dart_with_index);
      CGAL_postcondition(vertices_in_face_with_index == polygon->getExteriorRing()->getNumPoints()-1);
      vertices_in_face_with_index = 0;
      for (LCC::template Dart_of_cell_range<2, 2>::iterator current_dart_in_cell = lcci.template darts_of_cell<2, 2>(start_dart_with_index).begin(); current_dart_in_cell != lcci.template darts_of_cell<2, 2>(start_dart_with_index).end(); ++current_dart_in_cell) {
        if (same_orientation_with_index) CGAL_postcondition(lcci.point(current_dart_in_cell) == lcci.point(vertices_with_index[vertices_in_face_with_index]));
        else CGAL_postcondition(lcci.point(current_dart_in_cell) == lcci.point(vertices_with_index[vertices_with_index.size()-vertices_in_face_with_index-1]));
        ++vertices_in_face_with_index;
      } CGAL_postcondition(vertices_in_face_with_index == polygon->getExteriorRing()->getNumPoints()-1);
      validate_face_with_index_total_time += clock() - start_time;
      std::cout << ",";
    }
    
    std::cout << std::endl;
    
    // Validate that the faces form a closed and quasi-manifold volume
    start_time = clock();
    for (typename std::list<LCC::Dart_handle>::iterator dart_in_current_facet_1 = faces.begin(); dart_in_current_facet_1 != faces.end(); ++dart_in_current_facet_1) {
      for (typename LCC::template Dart_of_cell_range<2, 2>::iterator current_ridge_1 = lcc.template darts_of_cell<2, 2>(*dart_in_current_facet_1).begin(); current_ridge_1 != lcc.template darts_of_cell<2, 2>(*dart_in_current_facet_1).end(); ++current_ridge_1) {
        int matches = 0;
        
        for (typename std::list<LCC::Dart_handle>::iterator dart_in_current_facet_2 = faces.begin(); dart_in_current_facet_2 != faces.end(); ++dart_in_current_facet_2) {
          if (dart_in_current_facet_1 == dart_in_current_facet_2) continue;
          for (typename LCC::template Dart_of_cell_range<2, 2>::iterator current_ridge_2 = lcc.template darts_of_cell<2, 2>(*dart_in_current_facet_2).begin(); current_ridge_2 != lcc.template darts_of_cell<2, 2>(*dart_in_current_facet_2).end(); ++current_ridge_2) {
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<1>(lcc, current_ridge_1, current_ridge_2)) ++matches;
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<1>(lcc, current_ridge_1, current_ridge_2)) ++matches;
          }
        }
        
        CGAL_precondition(matches == 1);
      }
    } validate_volume_total_time += clock() - start_time;
    
    // Validate that the faces form a closed and quasi-manifold volume with index
    start_time = clock();
    for (typename std::list<LCC::Dart_handle>::iterator dart_in_current_facet_1 = faces_with_index.begin(); dart_in_current_facet_1 != faces_with_index.end(); ++dart_in_current_facet_1) {
      for (typename LCC::template Dart_of_cell_range<2, 2>::iterator current_ridge_1 = lcci.template darts_of_cell<2, 2>(*dart_in_current_facet_1).begin(); current_ridge_1 != lcci.template darts_of_cell<2, 2>(*dart_in_current_facet_1).end(); ++current_ridge_1) {
        int matches = 0;
        
        for (typename std::list<LCC::Dart_handle>::iterator dart_in_current_facet_2 = faces_with_index.begin(); dart_in_current_facet_2 != faces_with_index.end(); ++dart_in_current_facet_2) {
          if (dart_in_current_facet_1 == dart_in_current_facet_2) continue;
          for (typename LCC::template Dart_of_cell_range<2, 2>::iterator current_ridge_2 = lcci.template darts_of_cell<2, 2>(*dart_in_current_facet_2).begin(); current_ridge_2 != lcci.template darts_of_cell<2, 2>(*dart_in_current_facet_2).end(); ++current_ridge_2) {
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<1>(lcci, current_ridge_1, current_ridge_2)) ++matches;
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<1>(lcci, current_ridge_1, current_ridge_2)) ++matches;
          }
        }
        
        CGAL_precondition(matches == 1);
      }
    } validate_volume_with_index_total_time += clock() - start_time;
    
    // Volume
    start_time = clock();
    volumes.push_back(builder.template get_cell<3>(faces.begin(), faces.end()));
    get_cell_total_time = clock() - start_time;
    
    // Volume with index
    start_time = clock();
    volumes_with_index.push_back(builder_with_index.template get_cell<3>(faces_with_index.begin(), faces_with_index.end()));
    get_cell_with_index_total_time = clock() - start_time;
    
    // Validate that volumes created using both methods are equal
    CGAL_assertion(Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry<3>(lcc, volumes.back(), lcci, volumes_with_index.back()));
    
    // Times
    std::cout << "\tget_vertex():      " << get_vertex_total_time / double(CLOCKS_PER_SEC) << " s." << std::endl;
//    std::cout << "\tget_vertex():      " << get_vertex_with_index_total_time / double(CLOCKS_PER_SEC) << " s." << std::endl;
    std::cout << "\tget_vertex(i):     " << get_vertex_with_index_total_time / double(CLOCKS_PER_SEC) << " s. (" << double(get_vertex_total_time)/double(get_vertex_with_index_total_time) << "x)" << std::endl;
    
    std::cout << "\tget_face():        " << get_face_total_time / double(CLOCKS_PER_SEC) << " s." << std::endl;
//    std::cout << "\tget_face():        " << get_face_with_index_total_time / double(CLOCKS_PER_SEC) << " s." << std::endl;
    std::cout << "\tget_face(i):       " << get_face_with_index_total_time / double(CLOCKS_PER_SEC) << " s. (" << double(get_face_total_time)/double(get_face_with_index_total_time) << "x)" << std::endl;
    
//    std::cout << "\tvalidate faces:    " << validate_face_total_time / double(CLOCKS_PER_SEC) << " s." << std::endl;
//    std::cout << "\tvalidate faces:    " << validate_face_with_index_total_time / double(CLOCKS_PER_SEC) << " s." << std::endl;
//    std::cout << "\tvalidate faces i:  " << validate_face_with_index_total_time / double(CLOCKS_PER_SEC) << " s. (" << double(validate_face_total_time)/double(validate_face_with_index_total_time) << "x)" << std::endl;
    
//    std::cout << "\tvalidate volume:   " << validate_volume_total_time / double(CLOCKS_PER_SEC) << " s." << std::endl;
//    std::cout << "\tvalidate volume:   " << validate_volume_with_index_total_time / double(CLOCKS_PER_SEC) << " s." << std::endl;
//    std::cout << "\tvalidate volume i: " << validate_volume_with_index_total_time / double(CLOCKS_PER_SEC) << " s. (" << double(validate_volume_total_time)/double(validate_volume_with_index_total_time) << "x)" << std::endl;
    
    std::cout << "\tget_cell():        " << get_cell_total_time / double(CLOCKS_PER_SEC) << " s." << std::endl;
//    std::cout << "\tget_cell():        " << get_cell_with_index_total_time / double(CLOCKS_PER_SEC) << " s." << std::endl;
    std::cout << "\tget_cell(i):       " << get_cell_with_index_total_time / double(CLOCKS_PER_SEC) << " s. (" << double(get_cell_total_time)/double(get_cell_with_index_total_time) << "x)" << std::endl;
  }
  
  start_time = clock();
  builder.template link_cells<3>(volumes.begin(), volumes.end());
  std::cout << "\tw/o index:  " << (clock()-start_time) / double(CLOCKS_PER_SEC) << " s." << std::endl;
  CGAL_postcondition(lcc.is_valid());
  start_time = clock();
  builder_with_index.template link_cells<3>(volumes_with_index.begin(), volumes_with_index.end());
  std::cout << "\twith index: " << (clock()-start_time) / double(CLOCKS_PER_SEC) << " s." << std::endl;
  CGAL_postcondition(lcci.is_valid());
  
  // Validate that the output is equivalent
  std::cout << "Comparing built LCCs..." << std::endl;
  for (LCC::Dart_range::iterator current_dart = lcci.darts().begin(); current_dart != lcci.darts().end(); ++current_dart) {
    if (Linear_cell_complex_cell_comparator<LCC>::have_same_geometry_and_opposite_orientation_from<4>(lcc, lcc.darts().begin(), lcci, current_dart)) {
      std::cout << "\tMatched with opposite orientation!" << std::endl;
    } if (Linear_cell_complex_cell_comparator<LCC>::have_same_geometry_and_orientation_from<4>(lcc, lcc.darts().begin(), lcci, current_dart)) {
      std::cout << "\tMatched with same orientation!" << std::endl;
    }
  }
  
  // Output
  std::cout << "Output:" << std::endl;
  lcc.display_characteristics(std::cout);
  std::cout << std::endl;
  lcci.display_characteristics(std::cout);
  std::cout << std::endl;
  Linear_cell_complex_reader_writer<LCC>::write_obj("/Users/ken/Desktop/polyhedra.obj", lcc);
//  Linear_cell_complex_printer<LCC>::print(lcc);
}

int main(int argc, const char *argv[]) {
//  load_Martijns_format("/Users/ken/Doctorado/data/varioscale/utm.txt");
  load_Martijns_format("/Users/ken/Doctorado/data/varioscale/utm2.txt");
//  load_Martijns_format("/Users/ken/Doctorado/data/varioscale/atkis.txt");
  return 0;
}

