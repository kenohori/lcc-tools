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

#ifndef LCC_VALIDATOR_H
#define LCC_VALIDATOR_H

#include "Linear_cell_complex_with_ids.h"

template <class LCC_>
class Linear_cell_complex_validator {
public:
  typedef LCC_ LCC;
  typedef typename LCC::Dart_handle Dart_handle;
  typedef typename LCC::Point Point;
  
  static bool are_dart_ids_valid(LCC &lcc) {
    return true;
  }
  
  template <unsigned int dimension>
  static bool are_attribute_ids_valid(LCC &lcc) {
    return true;
  }
  
  static bool are_vertices_unique(LCC &lcc) {
    std::set<typename LCC::Point> vertices;
    for (typename LCC::Vertex_attribute_range::iterator current_vertex = lcc.vertex_attributes().begin(); current_vertex != lcc.vertex_attributes().end(); ++current_vertex) {
      if (vertices.count(current_vertex->point())) return false;
      else vertices.insert(current_vertex->point());
    } return true;
  }
  
  static void validate_map_betas(LCC &lcc) {
    for (typename LCC::Dart_range::iterator current_dart = lcc.darts().begin(); current_dart != lcc.darts().end(); ++current_dart) {
      if (!lcc.is_free(current_dart, 0)) CGAL_assertion(current_dart->beta(0)->beta(1) == current_dart);
      if (!lcc.is_free(current_dart, 1)) CGAL_assertion(current_dart->beta(1)->beta(0) == current_dart);
      for (unsigned int dimension = 2; dimension <= LCC::dimension; ++dimension) {
        if (!lcc.is_free(current_dart, dimension)) CGAL_assertion(current_dart->beta(dimension)->beta(dimension) == current_dart);
      }
    }
  }
};

template <unsigned int d>
class Linear_cell_complex_validator<Linear_cell_complex_with_ids<d> > {
public:
  typedef Linear_cell_complex_with_ids<d> LCC;
  typedef typename LCC::Dart_handle Dart_handle;
  typedef typename LCC::Point Point;
  
  static bool are_dart_ids_valid(LCC &lcc) {
    std::set<int> dart_ids;
    for (typename LCC::Dart_range::iterator current_dart = lcc.darts().begin(); current_dart != lcc.darts().end(); ++current_dart) {
      if (current_dart->id == -1) continue;
      if (dart_ids.count(current_dart->id)) return false;
      dart_ids.insert(current_dart->id);
    } if (*dart_ids.begin() != 0) return false;
    if (*dart_ids.rbegin() != dart_ids.size()-1) return false;
    return true;
  }
  
  template <unsigned int dimension>
  static bool are_attribute_ids_valid(LCC &lcc) {
    std::set<int> attribute_ids;
    for (typename LCC::template Attribute_range<dimension>::type::iterator current_attribute = lcc.template attributes<dimension>().begin(); current_attribute != lcc.template attributes<dimension>().end(); ++current_attribute) {
      if (current_attribute->info() == -1) continue;
      if (attribute_ids.count(current_attribute->info())) {
        return false;
      } attribute_ids.insert(current_attribute->info());
    } if (*attribute_ids.begin() != 0) return false;
    if (*attribute_ids.rbegin() != attribute_ids.size()-1) return false;
    return true;
  }
  
  static bool are_vertices_unique(LCC &lcc) {
    std::set<typename LCC::Point> vertices;
    for (typename LCC::Vertex_attribute_range::iterator current_vertex = lcc.vertex_attributes().begin(); current_vertex != lcc.vertex_attributes().end(); ++current_vertex) {
      if (vertices.count(current_vertex->point())) return false;
      else vertices.insert(current_vertex->point());
    } return true;
  }
  
  static void validate_map_betas(LCC &lcc) {
    for (typename LCC::Dart_range::iterator current_dart = lcc.darts().begin(); current_dart != lcc.darts().end(); ++current_dart) {
      if (!lcc.is_free(current_dart, 0)) CGAL_assertion(current_dart->beta(0)->beta(1) == current_dart);
      if (!lcc.is_free(current_dart, 1)) CGAL_assertion(current_dart->beta(1)->beta(0) == current_dart);
      for (unsigned int dimension = 2; dimension <= LCC::dimension; ++dimension) {
        if (!lcc.is_free(current_dart, dimension)) CGAL_assertion(current_dart->beta(dimension)->beta(dimension) == current_dart);
      }
    }
  }
  
  template <unsigned int dimension>
  static void validate_map_attributes(LCC &lcc) {
    for (typename LCC::template Attribute_range<dimension>::type::iterator current_attribute = lcc.template attributes<dimension>().begin(); current_attribute != lcc.template attributes<dimension>().end(); ++current_attribute) {
      CGAL_assertion(lcc.template attribute<dimension>(current_attribute->dart()) == current_attribute);
    }
  }
};

#endif
