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

#ifndef LCC_SMALLEST_VERTEX_INDEX_H
#define LCC_SMALLEST_VERTEX_INDEX_H

#include <CGAL/Linear_cell_complex.h>

template <class LCC>
struct Linear_cell_complex_smallest_vertex_index_for_vertices {
public:
  std::map<typename LCC::Point, typename LCC::Dart_handle> contents;
  
  void validate(LCC &lcc) {
    
    // Check that the vertices are really at this point
    for (typename std::map<typename LCC::Point, typename LCC::Dart_handle>::iterator current_vertex = contents.begin(); current_vertex != contents.end(); ++current_vertex) {
      CGAL_assertion(current_vertex->first == lcc.point(current_vertex->second));
    }
  }
  
  void update_index_after_reversing_orientation(LCC &lcc) {
    
    // Smart method
    for (typename std::map<typename LCC::Point, typename LCC::Dart_handle>::iterator current_vertex = contents.begin(); current_vertex != contents.end(); ++current_vertex) {
      if (current_vertex->first != lcc.point(current_vertex->second)) {
        current_vertex->second = current_vertex->second->beta(1);
        CGAL_postcondition(current_vertex->first == lcc.point(current_vertex->second));
      }
    }
  }
};

template <class LCC>
struct Linear_cell_complex_smallest_vertex_index_per_dimension {
public:
  std::map<typename LCC::Point, std::list<typename LCC::Dart_handle> > contents;
  
  void validate(LCC &lcc) {
    for (typename std::map<typename LCC::Point, std::list<typename LCC::Dart_handle> >::iterator current_vertex = contents.begin(); current_vertex != contents.end(); ++current_vertex) {
      
      // Check that the cells are really at this point
      for (typename std::list<typename LCC::Dart_handle>::iterator current_cell = current_vertex->second.begin(); current_cell != current_vertex->second.end(); ++current_cell) {
        CGAL_assertion(current_vertex->first == lcc.point(*current_cell));
      }
    }
  }
  
  void update_index_after_reversing_orientation(LCC &lcc) {
    for (typename std::map<typename LCC::Point, std::list<typename LCC::Dart_handle> >::iterator current_vertex = contents.begin(); current_vertex != contents.end(); ++current_vertex) {
      for (typename std::list<typename LCC::Dart_handle>::iterator current_cell = current_vertex->second.begin(); current_cell != current_vertex->second.end(); ++current_cell) {
        if (current_vertex->first != lcc.point(*current_cell)) {
          (*current_cell) = (*current_cell)->beta(1);
          CGAL_postcondition(current_vertex->first == lcc.point(*current_cell));
        }
      }
    }
  }
};

template <class LCC, class Index, unsigned int dimension = LCC::dimension>
struct Do_all {
  static void validate(LCC &lcc, Index &index) {
    Do_all<LCC, Index, dimension-1>::validate(lcc, index);
    std::get<dimension>(index).validate(lcc);
  }
  
  static void update_index_after_reversing_orientation(LCC &lcc, Index &index) {
    Do_all<LCC, Index, dimension-1>::update_index_after_reversing_orientation(lcc, index);
    std::get<dimension>(index).update_index_after_reversing_orientation(lcc);
  }
};

template <class LCC, class Index>
struct Do_all<LCC, Index, 0> {
  static void validate(LCC &lcc, Index &index) {
    std::get<0>(index).validate(lcc);
  }
  
  static void update_index_after_reversing_orientation(LCC &lcc, Index &index) {
    std::get<0>(index).update_index_after_reversing_orientation(lcc);
  }
};

template <class LCC>
class Linear_cell_complex_smallest_vertex_index {
protected:
  template <unsigned int dimensions_to_add, class Result = CGAL::cpp11::tuple<> >
  struct Linear_cell_complex_smallest_vertex_index_up_to_dimension;
  
  template <class ... Result>
  struct Linear_cell_complex_smallest_vertex_index_up_to_dimension<0, CGAL::cpp11::tuple<Result ...> > {
    typedef CGAL::cpp11::tuple<Linear_cell_complex_smallest_vertex_index_for_vertices<LCC>, Result ...> tuple;
  };
  
  template <unsigned int dimensions_to_add, class ... Result>
  struct Linear_cell_complex_smallest_vertex_index_up_to_dimension<dimensions_to_add, CGAL::cpp11::tuple<Result ...> > {
    typedef typename Linear_cell_complex_smallest_vertex_index_up_to_dimension<dimensions_to_add-1, CGAL::cpp11::tuple<Linear_cell_complex_smallest_vertex_index_per_dimension<LCC>, Result ...> >::tuple tuple;
  };
  
  typedef typename Linear_cell_complex_smallest_vertex_index_up_to_dimension<LCC::dimension>::tuple Index;
  Index index;
  
public:
  typedef std::list<typename LCC::Dart_handle> &Results;
  typedef typename std::list<typename LCC::Dart_handle>::iterator Result_iterator;
  
  std::list<typename LCC::Dart_handle> no_cells_found;
  
  void insert_vertex(LCC &lcc, typename LCC::Dart_handle dart) {
    CGAL_precondition(std::get<0>(index).contents.count(lcc.point(dart)) == 0);
    std::get<0>(index).contents[lcc.point(dart)] = dart;
  }
  
  template <unsigned int dimension>
  void insert_cell(LCC &lcc, typename LCC::Dart_handle dart) {
    CGAL_precondition(dimension >= 1);
    
    if (dimension == 1) {
      CGAL_precondition(!lcc.is_free(dart, 1));
      if (lcc.point(dart) < lcc.point(dart->beta(1))) std::get<1>(index).contents[lcc.point(dart)].push_back(dart);
      else std::get<1>(index).contents[lcc.point(dart->beta(1))].push_back(dart);
    }
    
    else {
      typename LCC::Dart_handle smallest_vertex = this->template find_smallest_vertex<dimension>(lcc, dart);
      std::get<dimension>(index).contents[lcc.point(smallest_vertex)].push_back(smallest_vertex);
      CGAL_postcondition(lcc.template attribute<dimension>(dart) == lcc.template attribute<dimension>(smallest_vertex));
    }
  }
  
  typename LCC::Dart_handle find_vertex(const typename LCC::Point &p) {
    if (std::get<0>(index).contents.count(p) != 0) return std::get<0>(index).contents[p];
    else return LCC::null_handle;
  }
  
  template <unsigned int dimension>
  Results find_cells(const typename LCC::Point &p) {
    CGAL_precondition(dimension >= 1);
    if (std::get<dimension>(index).contents.count(p) != 0) return std::get<dimension>(index).contents[p];
    else return no_cells_found;
  }
  
  void update_index_after_reversing_orientation(LCC &lcc) {
    Do_all<LCC, Index>::update_index_after_reversing_orientation(lcc, index);
  }
  
  template <unsigned int dimension>
  void validate_index(LCC &lcc) {
    std::get<dimension>(index).validate(lcc);
  }
  
  template <unsigned int dimension>
  static typename LCC::Dart_handle find_smallest_vertex(LCC &lcc, typename LCC::Dart_handle dh) {
    CGAL_precondition(dimension >= 2);
    typename LCC::Dart_handle smallest_vertex = dh;
    for (typename LCC::template Dart_of_cell_range<dimension, dimension>::iterator current_dart = lcc.template darts_of_cell<dimension, dimension>(dh).begin(); current_dart != lcc.template darts_of_cell<dimension, dimension>(dh).end(); ++current_dart) {
      if (lcc.point(current_dart) < lcc.point(smallest_vertex)) {
        smallest_vertex = current_dart;
      }
    } CGAL_postcondition(lcc.template attribute<dimension>(dh) == lcc.template attribute<dimension>(smallest_vertex));
    return smallest_vertex;
  }
};

#endif
