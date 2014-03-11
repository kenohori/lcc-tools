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

#ifndef LCC_INCREMENTAL_BUILDER_H
#define LCC_INCREMENTAL_BUILDER_H

#include <CGAL/Linear_cell_complex.h>
#include "Linear_cell_complex_utilities.h"

template <class LCC_>
class Linear_cell_complex_incremental_builder;

// Used to have a well defined order on Point_d. Already present in CGAL for Point_2 and Point_3.
template <class Traits>
bool operator<(const CGAL::Point_d<Traits> &p1, const CGAL::Point_d<Traits> &p2) {
  for (typename CGAL::Point_d<Traits>::Cartesian_const_iterator current_coordinate_p1 = p1.cartesian_begin(), current_coordinate_p2 = p2.cartesian_begin();
       current_coordinate_p1 != p1.cartesian_end() && current_coordinate_p2 != p2.cartesian_end();
       ++current_coordinate_p1, ++current_coordinate_p2) {
    if (*current_coordinate_p1 < *current_coordinate_p2) return true;
    if (*current_coordinate_p1 > *current_coordinate_p2) return false;
  } return false;
}

// Index using the smallest vertex of a cell, could be exchanged for something else
template <class LCC>
struct Index {
public:
  
  // TODO: What to do for attributes? Should we compare them as well?
  struct Compare_vertex_attribute_handle {
    bool operator() (const typename LCC::Vertex_attribute_handle &vah1, const typename LCC::Vertex_attribute_handle &vah2) const {
      return vah1->point() < vah2->point();
    }
  };
  
  void clear() {
    contents.clear();
  }
  
  std::size_t count(const typename LCC::Vertex_attribute_handle &vah) {
    return contents.count(vah);
  }
  
  // Warning: this returns the total number of unique indices (keys), NOT the number of indexed features!
  std::size_t size() {
    return contents.size();
  }
  
  std::size_t entries() {
    std::size_t total_entries = 0;
    for (typename std::map<typename LCC::Vertex_attribute_handle, std::list<typename LCC::Dart_handle> >::iterator current_vertex = contents.begin(); current_vertex != contents.end(); ++current_vertex) {
      total_entries += current_vertex->second.size();
    } return total_entries;
  }
  
  std::list<typename LCC::Dart_handle> &operator[](const typename LCC::Vertex_attribute_handle &vah) {
    return contents[vah];
  }
  
  typename std::map<typename LCC::Vertex_attribute_handle, std::list<typename LCC::Dart_handle> >::iterator begin() {
    return contents.begin();
  }
  
  typename std::map<typename LCC::Vertex_attribute_handle, std::list<typename LCC::Dart_handle> >::iterator end() {
    return contents.end();
  }
  
  typename std::map<typename LCC::Vertex_attribute_handle, std::list<typename LCC::Dart_handle> >::const_iterator begin() const {
    return contents.begin();
  }
  
  typename std::map<typename LCC::Vertex_attribute_handle, std::list<typename LCC::Dart_handle> >::const_iterator end() const {
    return contents.end();
  }
  
  template <class LCC_>
  friend std::ostream &operator<<(std::ostream &output, const Index<LCC_> &index);
  
private:
  std::map<typename LCC::Vertex_attribute_handle, std::list<typename LCC::Dart_handle>, Compare_vertex_attribute_handle> contents;
};

// Build a smallest vertex index in this dimension, O(n) on the number of darts
template <class LCC_builder, unsigned int dimension>
struct Index_builder {
  
  // Build an index for all the dimension-cells in the Linear_cell_complex lcc
  static void build_for_dimension(LCC_builder &lcc_builder) {
    //std::cout << "\tbuild_for_dimension<LCC, " << dimension << ">(...)" << std::endl;
    
    if (lcc_builder.indices.count(dimension) > 0) {
      //std::cout << "\t" << dimension << "-index already exists." << std::endl;
      return;
    }
    
    std::cout << "\tbuilding " << dimension << "-index..." << std::endl;
    lcc_builder.indices[dimension] = Index<typename LCC_builder::LCC>();
    int mark = lcc_builder.lcc.get_new_mark();
    CGAL_precondition(lcc_builder.lcc.is_whole_map_unmarked(mark));
    
    // Look for all possible dimension-cells
    typename LCC_builder::LCC::Dart_range &darts = lcc_builder.lcc.darts();
    for (typename LCC_builder::LCC::Dart_range::iterator current_dart = darts.begin(), last_dart = darts.end(); current_dart != last_dart; ++current_dart) {
      if (lcc_builder.lcc.is_marked(current_dart, mark)) continue;
      typename LCC_builder::LCC::template Dart_of_cell_range<dimension, dimension> darts_in_cell = lcc_builder.lcc.template darts_of_cell<dimension, dimension>(current_dart);
      typename LCC_builder::LCC::Dart_handle smallest_vertex_in_cell = darts_in_cell.begin();
      bool closed_cell = true;
      
      // For every dart in the dimension-cell
      for (typename LCC_builder::LCC::template Dart_of_cell_range<dimension, dimension>::iterator current_dart_in_cell = darts_in_cell.begin(); current_dart_in_cell != darts_in_cell.end(); ++current_dart_in_cell) {
        lcc_builder.lcc.mark(current_dart_in_cell, mark);
        if (!closed_cell) continue;
        
        // Check if the cell is not open (topologically closed) at lower dimensions
        for (int current_dimension = 0; current_dimension < dimension; ++current_dimension) {
          if (lcc_builder.lcc.is_free(current_dart_in_cell, current_dimension)) {
            closed_cell = false;
            break;
          }
        }
        
        // Check if the point at this dart is the smallest in the dimension-cell
        if (current_dart_in_cell->template attribute<0>()->point() < smallest_vertex_in_cell->template attribute<0>()->point())
          smallest_vertex_in_cell = current_dart_in_cell;
      }
      
      // Add to index
      if (closed_cell) lcc_builder.indices[dimension][smallest_vertex_in_cell->template attribute<0>()].push_back(smallest_vertex_in_cell);
    }
    
    lcc_builder.lcc.negate_mark(mark);
    CGAL_postcondition(lcc_builder.lcc.is_whole_map_unmarked(mark));
    lcc_builder.lcc.free_mark(mark);
    CGAL_postcondition(lcc_builder.indices.count(dimension) > 0);
    //std::cout << "\tIndex has " << lcc_builder.indices[dimension].entries() << " entries." << std::endl;
  }
  
  template <class Index, class InputIterator>
  static void build_for_cell_boundaries(LCC_builder &lcc_builder, Index &index, InputIterator first, InputIterator last) {
    //std::cout << "\tbuild_for_cell_boundaries<LCC, " << dimension << ">(...)" << std::endl;
    CGAL_precondition(dimension > 0);
    index.clear();
    int mark = lcc_builder.lcc.get_new_mark();
    CGAL_precondition(lcc_builder.lcc.is_whole_map_unmarked(mark));
    for (InputIterator current = first; current != last; ++current) {
      for (typename LCC_builder::LCC::template Dart_of_cell_range<dimension, dimension>::iterator current_dart_in_cell = lcc_builder.lcc.template darts_of_cell<dimension, dimension>(*current).begin(),
           last_dart_in_cell = lcc_builder.lcc.template darts_of_cell<dimension, dimension>(*current).end(); current_dart_in_cell != last_dart_in_cell; ++current_dart_in_cell) {
        if (lcc_builder.lcc.is_marked(current_dart_in_cell, mark)) continue;
        if (dimension == 2) {
          typename LCC_builder::LCC::Dart_handle smallest_vertex = lcc_builder.template find_smallest_vertex_in_cell<1>(current_dart_in_cell);
          index[smallest_vertex->template attribute<0>()].push_back(current_dart_in_cell);
          lcc_builder.lcc.mark(current_dart_in_cell, mark);
        } else {
          typename LCC_builder::LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart_in_boundary = lcc_builder.lcc.template darts_of_cell<dimension-1, dimension-1>(current_dart_in_cell).begin();
          typename LCC_builder::LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator last_dart_in_boundary = lcc_builder.lcc.template darts_of_cell<dimension-1, dimension-1>(current_dart_in_cell).end();
          typename LCC_builder::LCC::Dart_handle smallest_vertex = current_dart_in_boundary;
          lcc_builder.lcc.mark(current_dart_in_boundary, mark);
          ++current_dart_in_boundary;
          
          while (current_dart_in_boundary != last_dart_in_boundary) {
            if (current_dart_in_boundary->template attribute<0>()->point() < smallest_vertex->template attribute<0>()->point())
              smallest_vertex = current_dart_in_boundary;
            lcc_builder.lcc.mark(current_dart_in_boundary, mark);
            ++current_dart_in_boundary;
          } index[smallest_vertex->template attribute<0>()].push_back(smallest_vertex);
        }
      } for (typename LCC_builder::LCC::template Dart_of_cell_range<dimension, dimension>::iterator current_dart_in_cell = lcc_builder.lcc.template darts_of_cell<dimension, dimension>(*current).begin(),
             last_dart_in_cell = lcc_builder.lcc.template darts_of_cell<dimension, dimension>(*current).end(); current_dart_in_cell != last_dart_in_cell; ++current_dart_in_cell) {
        lcc_builder.lcc.unmark(current_dart_in_cell, mark);
      }
    } CGAL_postcondition(lcc_builder.lcc.is_whole_map_unmarked(mark));
    //std::cout << "\tRidge index has " << index.entries() << " entries." << std::endl;
  }
};

// Special case for dimension = 1. At this point, the orbit of an edge probably only includes one of its end points, so we need to use beta(1) too.
template <class LCC_builder>
struct Index_builder<LCC_builder, 1> {
  static void build_for_dimension(LCC_builder &lcc_builder) {
    //std::cout << "\tbuild_for_dimension<LCC, 1 (sp)>(...)" << std::endl;
    if (lcc_builder.indices.count(1) > 0) {
      std::cout << "\t1-index already exists." << std::endl;
      return;
    }
    
    //std::cout << "\tbuilding 1-index..." << std::endl;
    lcc_builder.indices[1] = Index<typename LCC_builder::LCC>();
    for (typename LCC_builder::LCC::Dart_range::iterator current_dart = lcc_builder.lcc.darts().begin(), last_dart = lcc_builder.lcc.darts().end(); current_dart != last_dart; ++current_dart) {
      if (lcc_builder.lcc.is_free(current_dart, 1)) continue;
      typename LCC_builder::LCC::Dart_handle smallest_vertex = lcc_builder.template find_smallest_vertex_in_cell<1>(current_dart);
      lcc_builder.indices[1][smallest_vertex->template attribute<0>()].push_back(current_dart);
    } CGAL_postcondition(lcc_builder.indices.count(1) > 0);
  }
  
  template <class Index, class InputIterator>
  static void build_for_cell_boundaries(LCC_builder &lcc_builder, Index &index, InputIterator first, InputIterator last) {
    //std::cout << "\tbuild_for_cell_boundaries<LCC, 1 (sp)>(...)" << std::endl;
    for (InputIterator current = first; current != last; ++current) {
      index[(*current)->template attribute<0>()].push_back(*current);
      if (lcc_builder.lcc.is_free(*current, 1)) index[(*current)->beta(1)->template attribute<0>()].push_back(*current);
    }
  }
};

template <class LCC>
std::ostream &operator<<(std::ostream &output, const Index<LCC> &index) {
  for (typename std::map<typename LCC::Vertex_attribute_handle, std::list<typename LCC::Dart_handle> >::const_iterator current_vertex = index.begin(); current_vertex != index.end(); ++current_vertex) {
    output << "\tVertex: " << current_vertex->first->point() << std::endl;
    for (typename std::list<typename LCC::Dart_handle>::const_iterator current_dart = current_vertex->second.begin(); current_dart != current_vertex->second.end(); ++current_dart) {
      output << "\t\tDart<" << &**current_dart << ">: " << (*current_dart)->template attribute<0>()->point() << std::endl;
      //output << "\t\tDart<" << &**current_dart << ">: " << (*current_dart)->template attribute<0>()->point() << " -> " << (*current_dart)->beta(1)->template attribute<0>()->point() << std::endl;
    }
  } return output;
}

template <class LCC_>
class Linear_cell_complex_incremental_builder {
public:
  typedef LCC_ LCC;
  typedef typename LCC::Dart_handle Dart_handle;
  typedef typename LCC::Dart Dart;
  typedef typename LCC::Point Point;
  typedef typename LCC::Helper Helper;
  typedef Linear_cell_complex_incremental_builder<LCC> Self;
  
protected:
  // Copies the attributes of dart_from to dart_to, changing the orientation, except for the dimension-dimensional one, if it exists.
  template <class LCC, unsigned int copy_dimension>
  struct Copy_flipped_attributes {
    template <unsigned int dimension>
    static void run(LCC *lcc, Dart_handle dart_from, Dart_handle dart_to) {
      if (dimension == 0) lcc->template set_attribute<0>(dart_to, dart_from->beta(1)->template attribute<0>());
      else if (dimension != copy_dimension+1 && dart_from->template attribute<dimension>() != NULL) lcc->template set_attribute<dimension>(dart_to, dart_from->template attribute<dimension>());
    }
  };
  
public:
  LCC &lcc;
  int number_of_darts, implicit_edges_created;
  std::map<unsigned int, Index<LCC> > indices;
  
  Linear_cell_complex_incremental_builder(LCC &alcc) : lcc(alcc) {
    number_of_darts = 0;
    implicit_edges_created = 0;
  }
  
  template <unsigned int dimension>
  void print_index() {
    if (indices.count(dimension) > 0) std::cout << indices[dimension] << std::endl;
  }
  
  // Make a copy of a cell with reversed orientation and link it to the original one along the beta of the current_dimension
  template <unsigned int dimension>
  Dart_handle clone_and_link_flipped_cell(Dart_handle &dart_in_cell) {
    //std::cout << "clone_and_link_flipped_cell<" << dimension << ">(...)" << std::endl;
    CGAL_precondition(dimension > 1);
    CGAL_precondition(LCC::dimension > dimension);
    
    // Create a copy of each dart
    std::map<Dart_handle, Dart_handle> dart_map;
    for (typename LCC::template Dart_of_cell_range<dimension, dimension>::iterator current_dart = lcc.template darts_of_cell<dimension, dimension>(dart_in_cell).begin(),
         last_dart = lcc.template darts_of_cell<dimension, dimension>(dart_in_cell).end(); current_dart != last_dart; ++current_dart) {
      dart_map[current_dart] = lcc.create_dart(number_of_darts);
      ++number_of_darts;
    }
    
    // Link them according to the combinatorial structure of the cell of dart_in_cell, and link the two cells together along the dimension-th dimension
    for (typename LCC::template Dart_of_cell_range<dimension, dimension>::iterator current_dart = lcc.template darts_of_cell<dimension, dimension>(dart_in_cell).begin(),
         last_dart = lcc.template darts_of_cell<dimension, dimension>(dart_in_cell).end(); current_dart != last_dart; ++current_dart) {
      lcc.basic_link_beta_1(dart_map[current_dart], dart_map[current_dart->beta(0)]);
      for (unsigned int current_dimension = 2; current_dimension < dimension; ++current_dimension) {
        lcc.basic_link_beta_for_involution(dart_map[current_dart], dart_map[current_dart->beta(current_dimension)], current_dimension);
      } //lcc.basic_link_beta_for_involution(current_dart, dart_map[current_dart], dimension+1); // TODO: Should we do this? It's more efficient here, but maybe not. It breaks the incremental construction concept.
      Helper::template Foreach_enabled_attributes<Copy_flipped_attributes<LCC, dimension> >::run(&lcc, current_dart, dart_map[current_dart]);
    } return dart_map[dart_in_cell];
  }
  
  // Find the smallest vertex in a cell, so as to build a smallest vertex index
  template <unsigned int dimension>
  Dart_handle find_smallest_vertex_in_cell(Dart_handle &dart_in_cell) {
    //std::cout << "find_smallest_vertex_in_cell<" << dimension << ">(...)" << std::endl;
    CGAL_precondition(dimension <= LCC::dimension);
    if (dimension > 1) {
      typename LCC::template Dart_of_cell_range<dimension, dimension>::iterator current_dart = lcc.template darts_of_cell<dimension, dimension>(dart_in_cell).begin();
      typename LCC::template Dart_of_cell_range<dimension, dimension>::iterator last_dart = lcc.template darts_of_cell<dimension, dimension>(dart_in_cell).end();
      Dart_handle smallest_vertex = current_dart;
      ++current_dart;
      while (current_dart != last_dart) {
        if (current_dart->template attribute<0>()->point() < smallest_vertex->template attribute<0>()->point())
          smallest_vertex = current_dart;
        ++current_dart;
      } return smallest_vertex;
    } else if (dimension == 1) {
      if (lcc.is_free(dart_in_cell, 1)) return dart_in_cell;
      if (dart_in_cell->template attribute<0>()->point() < dart_in_cell->beta(1)->template attribute<0>()->point())
        return dart_in_cell;
      else return dart_in_cell->beta(1);
    } else return dart_in_cell;
  }
  
  // Check if two cells are equal (combinatorially and geometrically)
  template <unsigned int dimension>
  bool are_orbits_equal(Dart_handle dh1, Dart_handle dh2) {
    /*std::cout << "are_orbits_equal?" << std::endl;
     std::cout << "\t[(" << dh1->template attribute<0>()->point() << ")->(" << dh1->beta(1)->template attribute<0>()->point() << ")]" << std::endl;
     std::cout << "\t[(" << dh2->template attribute<0>()->point() << ")->(" << dh2->beta(1)->template attribute<0>()->point() << ")]" << std::endl;*/
    CGAL_precondition(dimension <= LCC::dimension);
    if (dh1 == dh2) return true;
    
    if (dimension == 0) return dh1->template attribute<0>() == dh2->template attribute<0>();
    if (dimension == 1) {
      if (dh1->template attribute<0>() != dh2->template attribute<0>()) return false;
      if (dh1->beta(1)->template attribute<0>() != dh2->beta(1)->template attribute<0>()) return false;
      return true;
    }
    
    bool match = true;
    
    int mark1 = lcc.get_new_mark();
    int mark2 = lcc.get_new_mark();
    CGAL_precondition(lcc.is_whole_map_unmarked(mark1));
    CGAL_precondition(lcc.is_whole_map_unmarked(mark2));
    std::stack<Dart_handle> to_check1, to_check2, to_unmark1, to_unmark2;
    
    to_check1.push(dh1);
    to_check2.push(dh2);
    lcc.mark(dh1, mark1);
    lcc.mark(dh2, mark2);
    to_unmark1.push(dh1);
    to_unmark2.push(dh2);
    
    Dart_handle dart1, dart2;
    std::map<Dart_handle, Dart_handle> bijection;
    
    while (match && !to_check1.empty()) {
      dart1 = to_check1.top();
      to_check1.pop();
      dart2 = to_check2.top();
      to_check2.pop();
      
      // Check vertices
      if (dart1->template attribute<0>() != dart2->template attribute<0>()) { // use attribtue directly
        match = false;
        break;
      }
      
      bijection[dart1] = dart2;
      
      for (unsigned int current_dimension = 1; current_dimension < dimension; ++current_dimension) {
        // Check betas
        if (lcc.is_free(dart1, current_dimension)) {
          if (!lcc.is_free(dart2, current_dimension)) {
            match = false;
            break;
          } continue;
        } if (lcc.is_free(dart2, current_dimension)) {
          match = false;
          break;
        }
        
        // Check neighbours
        if (lcc.is_marked(dart1->beta(current_dimension), mark1)) {
          if (!lcc.is_marked(dart2->beta(current_dimension), mark2)) {
            match = false;
            break;
          } if (bijection.count(dart1->beta(current_dimension)) && bijection.count(dart2->beta(current_dimension))) {
            if (bijection[dart1->beta(current_dimension)] != dart2->beta(current_dimension)) {
              match = false;
              break;
            }
          } continue;
        } if (lcc.is_marked(dart2->beta(current_dimension), mark2)) {
          match = false;
          break;
        }
        
        to_check1.push(dart1->beta(current_dimension));
        to_check2.push(dart2->beta(current_dimension));
        lcc.mark(dart1->beta(current_dimension), mark1);
        lcc.mark(dart2->beta(current_dimension), mark2);
        to_unmark1.push(dart1->beta(current_dimension));
        to_unmark2.push(dart2->beta(current_dimension));
      }
    }
    
    // Unmark darts
    while (!to_unmark1.empty()) {
      lcc.unmark(to_unmark1.top(), mark1);
      to_unmark1.pop();
      lcc.unmark(to_unmark2.top(), mark2);
      to_unmark2.pop();
    }
    
    // Free marks
    CGAL_postcondition(lcc.is_whole_map_unmarked(mark1));
    CGAL_postcondition(lcc.is_whole_map_unmarked(mark2));
    lcc.free_mark(mark1);
    lcc.free_mark(mark2);
    
    return match;
  }
  
  // Check if two cells are equal but with opposite orientations (combinatorially and geometrically)
  template <unsigned int dimension>
  bool are_orbits_equal_but_oppositely_oriented(Dart_handle dh1, Dart_handle dh2) {
    /*std::cout << "are_orbits_equal_but_oppositely_oriented?" << std::endl;
     std::cout << "\t[(" << dh1->template attribute<0>()->point() << ")->(" << dh1->beta(1)->template attribute<0>()->point() << ")]" << std::endl;
     std::cout << "\t[(" << dh2->beta(1)->template attribute<0>()->point() << ")<-(" << dh2->template attribute<0>()->point() << ")]" << std::endl;*/
    CGAL_precondition(dimension <= LCC::dimension);
    
    if (dh1 == dh2) return false;
    
    if (dimension == 0) return dh1->template attribute<0>() == dh2->template attribute<0>();
    if (dimension == 1) {
      if (dh1->template attribute<0>() != dh2->beta(1)->template attribute<0>()) return false;
      if (dh1->beta(1)->template attribute<0>() != dh2->template attribute<0>()) return false;
      return true;
    }
    
    bool match = true;
    
    int mark1 = lcc.get_new_mark();
    int mark2 = lcc.get_new_mark();
    CGAL_precondition(lcc.is_whole_map_unmarked(mark1));
    CGAL_precondition(lcc.is_whole_map_unmarked(mark2));
    std::stack<Dart_handle> to_check1, to_check2, to_unmark1, to_unmark2;
    
    to_check1.push(dh1);
    to_check2.push(dh2);
    lcc.mark(dh1, mark1);
    lcc.mark(dh2, mark2);
    to_unmark1.push(dh1);
    to_unmark2.push(dh2);
    
    Dart_handle dart1, dart2;
    std::map<Dart_handle, Dart_handle> bijection;
    
    while (match && !to_check1.empty()) {
      dart1 = to_check1.top();
      to_check1.pop();
      dart2 = to_check2.top();
      to_check2.pop();
      
      // Check vertices
      if (!lcc.is_free(dart1, 0) && !lcc.is_free(dart2, 1)) {
        if (dart1->template attribute<0>() != dart2->beta(1)->template attribute<0>()) {
          match = false;
          break;
        }
      }
      
      bijection[dart1] = dart2;
      
      for (unsigned int current_dimension = 1; current_dimension < dimension; ++current_dimension) {
        // Check betas
        if (lcc.is_free(dart1, current_dimension)) {
          if (!lcc.is_free(dart2, CGAL_BETAINV(current_dimension))) {
            match = false;
            break;
          } continue;
        } if (lcc.is_free(dart2, CGAL_BETAINV(current_dimension))) {
          match = false;
          break;
        }
        
        // Check neighbours
        if (lcc.is_marked(dart1->beta(current_dimension), mark1)) {
          if (!lcc.is_marked(dart2->beta(CGAL_BETAINV(current_dimension)), mark2)) {
            match = false;
            break;
          } if (bijection.count(dart1->beta(current_dimension)) && bijection.count(dart2->beta(CGAL_BETAINV(current_dimension)))) {
            if (bijection[dart1->beta(current_dimension)] != dart2->beta(CGAL_BETAINV(current_dimension))) {
              match = false;
              break;
            }
          } continue;
        } if (lcc.is_marked(dart2->beta(CGAL_BETAINV(current_dimension)), mark2)) {
          match = false;
          break;
        }
        
        to_check1.push(dart1->beta(current_dimension));
        to_check2.push(dart2->beta(CGAL_BETAINV(current_dimension)));
        lcc.mark(dart1->beta(current_dimension), mark1);
        lcc.mark(dart2->beta(CGAL_BETAINV(current_dimension)), mark2);
        to_unmark1.push(dart1->beta(current_dimension));
        to_unmark2.push(dart2->beta(CGAL_BETAINV(current_dimension)));
      }
    }
    
    // Unmark darts
    while (!to_unmark1.empty()) {
      lcc.unmark(to_unmark1.top(), mark1);
      to_unmark1.pop();
      lcc.unmark(to_unmark2.top(), mark2);
      to_unmark2.pop();
    }
    
    // Free marks
    CGAL_postcondition(lcc.is_whole_map_unmarked(mark1));
    CGAL_postcondition(lcc.is_whole_map_unmarked(mark2));
    lcc.free_mark(mark1);
    lcc.free_mark(mark2);
    
    return match;
  }
  
  // After reversing the orientation of a cell, it is necessary to update the indices
  
  void update_index_after_flipping(unsigned int index_to_check_dimension) {
    //std::cout << "update_index_after_flipping(" << index_to_check_dimension << ")" << std::endl;
    if (indices.count(index_to_check_dimension) > 0) {
      for (typename std::map<typename LCC::Vertex_attribute_handle, std::list<typename LCC::Dart_handle> >::iterator current_vertex = indices[index_to_check_dimension].begin(); current_vertex != indices[index_to_check_dimension].end(); ++current_vertex) {
        for (typename std::list<typename LCC::Dart_handle>::iterator current_dart = current_vertex->second.begin(); current_dart != current_vertex->second.end(); ++current_dart) {
          if ((*current_dart)->template attribute<0>() != current_vertex->first) {
            if ((*current_dart)->beta(1)->template attribute<0>() == current_vertex->first) {
              //std::cout << "\tFlipped " << index_to_check_dimension << "-cell." << std::endl;
              *current_dart = (*current_dart)->beta(1);
            } else {
              std::cout << "Error: invalid index!" << std::endl;
            }
          }
        }
      }
    }
  }
  
  void update_index_after_flipping(Index<LCC> &index_to_check) {
    //std::cout << "update_index_after_flipping(...)" << std::endl;
    for (typename std::map<typename LCC::Vertex_attribute_handle, std::list<typename LCC::Dart_handle> >::iterator current_vertex = index_to_check.begin(); current_vertex != index_to_check.end(); ++current_vertex) {
      for (typename std::list<typename LCC::Dart_handle>::iterator current_dart = current_vertex->second.begin(); current_dart != current_vertex->second.end(); ++current_dart) {
        if ((*current_dart)->template attribute<0>() != current_vertex->first) {
          if ((*current_dart)->beta(1)->template attribute<0>() == current_vertex->first) {
            *current_dart = (*current_dart)->beta(1);
          } else {
            std::cout << "Error: invalid index!" << std::endl;
          }
        }
      }
    }
  }
  
  // -------------------------------------------------------------------------------------------------------------------
  // FIND FUNCTIONS
  //
  // Used to check if a certain cell exists beforehand, so as to not create duplicate cells
  // For 0-, 1- and 2-cells, there is a natural order of their boundaries, so this operation makes sense.
  // For higher dimensional cells, it is simpler to create a cell first and then use the isomorphism test to see if they are equal
  // An index is used if there is one
  //
  // The Dart_handle returned is for a Dart in the cell if it exists, null_dart_handle otherwise
  // A bool indicates whether the cell has the same orientation as the input
  // -------------------------------------------------------------------------------------------------------------------
  
  // Find if a vertex (given by a 0-attribute) already exists
  Dart_handle find_vertex(typename LCC::Vertex_attribute_handle vertex_handle, bool free = false) {
    
    if (indices.count(0) > 0) {
      //std::cout << "\tusing 0-index." << std::endl;
      std::list<Dart_handle> existing_vertices = indices[0][vertex_handle];
      for (typename std::list<Dart_handle>::iterator current_dart = existing_vertices.begin(); current_dart != existing_vertices.end(); ++current_dart) {
        if (!free || (lcc.is_free(*current_dart, 0) && lcc.is_free(*current_dart, 1))) {
          return *current_dart;
        }
      }
    }
    
    else {
      std::cout << "\tWarning: no 0-index available!" << std::endl;
      for (typename LCC::Dart_range::iterator current_dart = lcc.darts().begin(), last_dart = lcc.darts().end(); current_dart != last_dart; ++current_dart) {
        if (current_dart->template attribute<0>()->point() == vertex_handle->point()) {
          if (!free || (lcc.is_free(current_dart, 0) && lcc.is_free(current_dart, 1))) {
            return current_dart;
          }
        }
      }
    } return LCC::null_handle;
  }
  
  // Find if an edge already exists, as defined by the two end vertices specified by darts
  std::pair<Dart_handle, bool> find_edge(Dart_handle start, Dart_handle end, bool free = false) {
    
    if (indices.count(1) > 0) {
      if (start->template attribute<0>()->point() < end->template attribute<0>()->point()) {
        std::list<Dart_handle> possible_matches = indices[1][start->template attribute<0>()];
        for (typename std::list<Dart_handle>::iterator current_match = possible_matches.begin(); current_match != possible_matches.end(); ++current_match) {
          if (!lcc.is_free(*current_match, 1)) {
            if ((*current_match)->beta(1)->template attribute<0>()->point() == end->template attribute<0>()->point()) {
              if (!free || (lcc.is_free(*current_match, 0) && lcc.is_free((*current_match)->beta(1), 1))) {
                return std::pair<Dart_handle, bool>(*current_match, true);
              }
            }
          } else {
            if ((*current_match)->beta(0)->template attribute<0>()->point() == end->template attribute<0>()->point()) {
              if (!free || (lcc.is_free(*current_match, 1) && lcc.is_free((*current_match)->beta(0), 0))) {
                return std::pair<Dart_handle, bool>(*current_match, false);
              }
            }
          }
        }
      }
      
      else {
        std::list<Dart_handle> possible_matches = indices[1][end->template attribute<0>()];
        for (typename std::list<Dart_handle>::iterator current_match = possible_matches.begin(); current_match != possible_matches.end(); ++current_match) {
          if (!lcc.is_free(*current_match, 1)) {
            if ((*current_match)->beta(1)->template attribute<0>()->point() == start->template attribute<0>()->point()) {
              if (!free || (lcc.is_free(*current_match, 0) && lcc.is_free((*current_match)->beta(1), 1))) {
                return std::pair<Dart_handle, bool>((*current_match)->beta(1), false);
              }
            }
          } else {
            if ((*current_match)->beta(0)->template attribute<0>()->point() == start->template attribute<0>()->point()) {
              if (!free || (lcc.is_free(*current_match, 1) && lcc.is_free((*current_match)->beta(0), 0))) {
                return std::pair<Dart_handle, bool>((*current_match)->beta(0), true);
              }
            }
          }
        }
      }
    }
    
    else {
      for (typename LCC::Dart_range::iterator current_dart =  lcc.darts().begin(), last_dart =  lcc.darts().end(); current_dart != last_dart; ++current_dart) {
        if (current_dart->template attribute<0>()->point() == start->template attribute<0>()->point()) {
          if (!lcc.is_free(current_dart, 1)) {
            if (current_dart->beta(1)->template attribute<0>()->point() == end->template attribute<0>()->point()) {
              if (!free || (lcc.is_free(current_dart, 0) && lcc.is_free(current_dart->beta(1), 1))) {
                return std::pair<Dart_handle, bool>(current_dart, true);
              }
            }
          }
        } if (current_dart->template attribute<0>()->point() == end->template attribute<0>()->point()) {
          if (!lcc.is_free(current_dart, 1)) {
            if (current_dart->beta(1)->template attribute<0>()->point() == start->template attribute<0>()->point()) {
              if (!free || (lcc.is_free(current_dart, 0) && lcc.is_free(current_dart->beta(1), 1))) {
                return std::pair<Dart_handle, bool>(current_dart->beta(1), false);
              }
            }
          }
        }
      }
    } return std::pair<Dart_handle, bool>(LCC::null_handle, true);
  }
  
  // Find if a facet already exists, defined as a loop of darts given by the iterators. Both orientations are checked.
  template <class InputIterator>
  std::pair<Dart_handle, bool> find_facet_from_vertices(InputIterator first, InputIterator last, bool free = false) {
    
    if (indices.count(2) > 0) {
      //std::cout << "\tusing 2-index." << std::endl;
      
      // Find the smallest vertex among the given ones
      InputIterator current = first;
      Dart_handle smallest = *current;
      ++current;
      while (current != last) {
        if ((*current)->template attribute<0>()->point() < smallest->template attribute<0>()->point())
          smallest = *current;
        ++current;
      } bool found;
      
      // Search among the possible matches in the index
      std::list<Dart_handle> possible_matches = indices[2][smallest->template attribute<0>()];
      for (typename std::list<Dart_handle>::iterator current_match = possible_matches.begin(); current_match != possible_matches.end(); ++current_match) {
        Dart_handle first_dart = *current_match;
        found = true;
        
        // See if the first vertex is among the matched facet, if so we start the comparison there
        while (first_dart->template attribute<0>()->point() != (*first)->template attribute<0>()->point()) {
          first_dart = first_dart->beta(1);
          if (first_dart == *current_match) {
            found = false;
            break;
          }
        } if (!found) continue;
        
        // Comparison with equal orientation
        Dart_handle current_dart = first_dart;
        current = first;
        do {
          if ((*current)->template attribute<0>()->point() != current_dart->template attribute<0>()->point()) {
            found = false;
            break;
          } ++current;
          current_dart = current_dart->beta(1);
        } while (current != first);
        if (found) {
          if (!free || lcc.is_free(current_dart, 2)) {
            return std::pair<Dart_handle, bool>(current_dart, true);
          }
        }
        
        // Comparison with opposite orientation
        current_dart = first_dart;
        current = first;
        found = true;
        do {
          if ((*current)->template attribute<0>()->point() != current_dart->template attribute<0>()->point()) {
            found = false;
            break;
          } ++current;
          current_dart = current_dart->beta(0);
        } while (current != first);
        if (found) {
          if (!free || lcc.is_free(current_dart, 2)) {
            return std::pair<Dart_handle, bool>(current_dart, false);
          }
        }
      }
    }
    
    else {
      std::cout << "\tWarning: no 2-index available!" << std::endl;
      for (typename LCC::template One_dart_per_cell_range<2>::iterator current_facet = lcc.template one_dart_per_cell<2>().begin(), last_facet = lcc.template one_dart_per_cell<2>().end(); current_facet != last_facet; ++current_facet) {
        
        // Turn and try to find a dart that matches the first dart of the given facet
        Dart_handle first_dart = current_facet;
        bool found = false;
        do {
          if (first_dart->template attribute<0>()->point() == (*first)->template attribute<0>()->point()) {
            found = true;
            break;
          } if (lcc.is_free(first_dart, 1)) break;
          first_dart = first_dart->beta(1);
        } while (first_dart != current_facet);
        if (!found) continue;
        
        // Try matching it with equal orientation
        Dart_handle current_dart = first_dart;
        InputIterator current = first;
        do {
          if ((*current)->template attribute<0>()->point() != current_dart->template attribute<0>()->point()) {
            found = false;
            break;
          } ++current;
          if (lcc.is_free(current_dart, 1)) {
            found = false;
            break;
          } current_dart = current_dart->beta(1);
        } while (current != first);
        if (found) {
          if (!free || lcc.is_free(current_dart, 2)) {
            return std::pair<Dart_handle, bool>(current_dart, true);
          }
        }
        
        // Try matching it with opposite orientation
        current_dart = first_dart;
        current = first;
        found = true;
        do {
          if ((*current)->template attribute<0>()->point() != current_dart->template attribute<0>()->point()) {
            found = false;
            break;
          } ++current;
          if (lcc.is_free(current_dart, 0)) {
            found = false;
            break;
          } current_dart = current_dart->beta(0);
        } while (current != first);
        if (found) {
          if (!free || lcc.is_free(current_dart, 2)) {
            return std::pair<Dart_handle, bool>(current_dart, false);
          }
        }
      }
    } return std::pair<Dart_handle, bool>(LCC::null_handle, true);
  }
  
  // TODO: Find if a facet already exists, defined as a loop of darts given by the iterators. Both orientations are checked.
  template <class InputIterator>
  std::pair<Dart_handle, bool> find_facet_from_edges(InputIterator first, InputIterator last, bool free = false) {
    
  }
  
  // NOTE: Other find functions are not useful since there is no natural order of the i-1 cells around an i-cell for i > 2
  
  // -------------------------------------------------------------------------------------------------------------------
  // LINK FUNCTIONS
  //
  // Used to link a set of i-cells along their common boundaries
  // The orientation of some of the i-cells is reversed if necessary
  //
  // The Dart_handle returned is for a Dart in the linked cells
  // -------------------------------------------------------------------------------------------------------------------
  
  // TODO: link_edges
  template <class InputIterator>
  Dart_handle link_edges(InputIterator first, InputIterator last) {
    
  }
  
  template <unsigned int dimension, class InputIterator>
  Dart_handle link_cells(InputIterator first, InputIterator last) {
    CGAL_precondition(dimension > 1);
    //std::cout << "link_cells<" << dimension << ">(...)" << std::endl;
    
    // Build a ridge index to make the process fast
    Index<LCC> ridge_index;
    Index_builder<Self, dimension>::build_for_cell_boundaries(*this, ridge_index, first, last);
    
    // Go face by face
    InputIterator current = first;
    while (current != last) {
      //std::cout << "\tProcessing " << dimension << "-cell..." << std::endl;
      std::list<Dart_handle> boundary_cells;
      
      // Get one dart per ridge in the face
      for (typename LCC::template One_dart_per_incident_cell_range<dimension-1, dimension>::iterator current_boundary_cell = lcc.template one_dart_per_incident_cell<dimension-1, dimension>(*current).begin(),
           last_boundary_cell = lcc.template one_dart_per_incident_cell<dimension-1, dimension>(*current).end(); current_boundary_cell != last_boundary_cell; ++current_boundary_cell)
        boundary_cells.push_back(current_boundary_cell);
      //std::cout << "\t\t" << boundary_cells.size() << " " << dimension-1 << "D ridges" << std::endl;
      
      // Go ridge by ridge (for a single face)
      for (typename std::list<Dart_handle>::iterator current_boundary_cell = boundary_cells.begin(); current_boundary_cell != boundary_cells.end(); ++current_boundary_cell) {
        //std::cout << "\t\tTesting " << dimension-1 << "-cell in the boundary of the " << dimension << "-cell..." << std::endl;
        
        // The ridge may have already been matched
        // TODO: Make sure that ridges are detached before?
        if (!lcc.is_free(*current_boundary_cell, dimension)) {
          //std::cout << "\t\t\tAlready matched" << std::endl;
          continue;
        }
        
        // Build a smallest vertex index to match the current ridge
        Dart_handle smallest_vertex = find_smallest_vertex_in_cell<dimension-1>(*current_boundary_cell);
        std::list<Dart_handle> possible_matches = ridge_index[smallest_vertex->template attribute<0>()];
        //std::cout << "\t\t\t" << possible_matches.size() << " possible matches found in index" << std::endl;
        
        // Possible results of the ridge matching process
        // TODO: Enforce manifoldness here by checking if there are multiple matches. It goes in place of the "already used" optimisation.
        for (typename std::list<Dart_handle>::iterator current_match = possible_matches.begin(); current_match != possible_matches.end() && lcc.is_free(*current_boundary_cell, dimension); ++current_match) { // current_match = dimension-2 cells in index
          
          // NOTE: For dimension = 2 the smallest vertex cannot be used since ridges are edges
          if (dimension == 2) {
            if (*current_match == *current_boundary_cell) continue; // Same dart
            if (!lcc.is_free(*current_match, dimension)) continue; // Already used (optimisation)
            if (are_orbits_equal<dimension-1>(*current_match, *current_boundary_cell)) {
              //std::cout << "\t\t\tMatch! Flipping..." << std::endl;
              lcc.reverse_orientation_connected_component(*current_boundary_cell);
              CGAL_assertion(lcc.is_free(*current_boundary_cell, dimension));
              lcc.template link_beta<2>(*current_match, *current_boundary_cell, false);
              break;
            } else if (are_orbits_equal_but_oppositely_oriented<1>(*current_match, *current_boundary_cell)) {
              //std::cout << "\t\t\tMatch opposite!" << std::endl;
              CGAL_assertion(lcc.is_free(*current_boundary_cell, dimension));
              lcc.template link_beta<2>(*current_match, *current_boundary_cell, false);
              break;
            }
          }
          
          // Optimisation for dimension = 3, since ridges are faces and therefore have a single dart at the smallest vertex
          else if (dimension == 3) {
            if (*current_match == smallest_vertex) continue; // Same dart
            if (!lcc.is_free(*current_match, dimension)) continue; // Already used (optimisation)
            if (are_orbits_equal<dimension-1>(*current_match, smallest_vertex)) {
              //std::cout << "\t\t\tMatch! Flipping..." << std::endl;
              lcc.reverse_orientation_connected_component(smallest_vertex);
              update_index_after_flipping(ridge_index);
              CGAL_assertion(lcc.is_free(*current_match, dimension));
              CGAL_assertion(are_orbits_equal_but_oppositely_oriented<dimension-1>(*current_match, smallest_vertex));
              lcc.template sew<dimension>(*current_match, smallest_vertex);
              break;
            } else if (are_orbits_equal_but_oppositely_oriented<dimension-1>((*current_match)->beta(0), smallest_vertex)) {
              //std::cout << "\t\t\tMatch opposite!" << std::endl;
              CGAL_assertion(lcc.is_free((*current_match)->beta(0), dimension));
              lcc.template sew<dimension>((*current_match)->beta(0), smallest_vertex);
              break;
            }
          }
          
          // General case. Looping around the smallest vertex is required
          else {
            for (typename LCC::template Dart_of_cell_range<0>::iterator current_dart = lcc.template darts_of_cell<0>(*current_match).begin(),
                 last_dart = lcc.template darts_of_cell<0>(*current_match).end(); current_dart != last_dart; ++current_dart) {
              // current_dart = darts around the smallest vertex (any of them is a valid starting point for orbit comparisons)
              // TODO: This iteration can be improved by following only a single beta.
              if (current_dart->template attribute<0>() != (*current_match)->template attribute<0>()) continue; // no incoming darts, just outgoing
              if (current_dart == smallest_vertex) continue; // Same dart
              if (!lcc.is_free(current_dart, dimension)) continue;  // Already used (optimisation)
              if (are_orbits_equal<dimension-1>(current_dart, smallest_vertex)) {
                //std::cout << "\t\t\tMatch! Flipping..." << std::endl;
                lcc.reverse_orientation_connected_component(smallest_vertex);
                update_index_after_flipping(ridge_index);
                CGAL_assertion(lcc.is_free(current_dart, dimension));
                CGAL_assertion(are_orbits_equal_but_oppositely_oriented<dimension-1>(current_dart, smallest_vertex));
                lcc.template sew<dimension>(current_dart, smallest_vertex);
                break;
              } else if (are_orbits_equal_but_oppositely_oriented<dimension-1>(current_dart->beta(0), smallest_vertex)) {
                //std::cout << "\t\t\tMatch opposite!" << std::endl;
                CGAL_assertion(lcc.is_free(current_dart->beta(0), dimension));
                lcc.template sew<dimension>(current_dart->beta(0), smallest_vertex);
                break;
              }
            }
          }
        }
        
        // Check if the ridge remains unmatched
        if (lcc.is_free(*current_boundary_cell, dimension)) {
          // TODO: Meaningful reporting? Not necessary when the faces don't form a closed cell
          //std::cout << "\t\t\tWarning: Unmatched ridge!" << std::endl;
        }
      } ++current;
    }
    
    // Create unique 1-attributes for each edge
    if (dimension == 2) {
      int mark = lcc.get_new_mark();
      CGAL_precondition(lcc.is_whole_map_unmarked(mark));
      current = first;
      
      // Check every input facet
      while (current != last) {
        for (typename LCC::template Dart_of_cell_range<2, 2>::iterator current_dart_in_facet = lcc.template darts_of_cell<2, 2>(*current).begin(),
             last_dart_in_facet = lcc.template darts_of_cell<2, 2>(*current).end();
             current_dart_in_facet != last_dart_in_facet;
             ++current_dart_in_facet) {
          
          // An unmarked dart means an unprocessed edge
          if (lcc.is_marked(current_dart_in_facet, mark)) continue;
          typename LCC::template Attribute_handle<1>::type edge_attribute = lcc.template create_attribute<1>(implicit_edges_created);
          ++implicit_edges_created;
          lcc.template set_attribute<1>(current_dart_in_facet, edge_attribute);
          for (typename LCC::template Dart_of_cell_range<1>::iterator current_dart_in_edge = lcc.template darts_of_cell<1>(current_dart_in_facet).begin(),
               last_dart_in_edge = lcc.template darts_of_cell<1>(current_dart_in_facet).end();
               current_dart_in_edge != last_dart_in_edge;
               ++current_dart_in_edge) {
            lcc.mark(current_dart_in_edge, mark);
          }
        }
        
        ++current;
      }
      
      // Unmark darts
      current = first;
      while (current != last) {
        for (typename LCC::template Dart_of_cell_range<2, 2>::iterator current_dart_in_facet = lcc.template darts_of_cell<2, 2>(*current).begin(),
             last_dart_in_facet = lcc.template darts_of_cell<2, 2>(*current).end();
             current_dart_in_facet != last_dart_in_facet;
             ++current_dart_in_facet) {
          for (typename LCC::template Dart_of_cell_range<1>::iterator current_dart_in_edge = lcc.template darts_of_cell<1>(current_dart_in_facet).begin(),
               last_dart_in_edge = lcc.template darts_of_cell<1>(current_dart_in_facet).end();
               current_dart_in_edge != last_dart_in_edge;
               ++current_dart_in_edge) {
            lcc.unmark(current_dart_in_edge, mark);
          }
        } ++current;
      }
      
      CGAL_postcondition(lcc.is_whole_map_unmarked(mark));
      lcc.free_mark(mark);
    }
    
    return *first;
  }
  
  // -------------------------------------------------------------------------------------------------------------------
  // CREATE FUNCTIONS
  //
  // Used to create a new i-cell from the given i-1 cells on its boundary (or a 0-attribute for 0-cells)
  // The orientation of some of the i-1 cells is reversed if necessary
  // The created cell is added to the i-index (if it exists)
  // A dart in the newly created cell is returned
  //
  // The Dart_handle returned is for a Dart in the newly created cell
  // -------------------------------------------------------------------------------------------------------------------
  
  Dart_handle create_dart_for_vertex_and_index(typename LCC::Vertex_attribute_handle vertex_handle) {
    Dart_handle new_dart = lcc.create_dart(number_of_darts);
    ++number_of_darts;
    lcc.set_vertex_attribute(new_dart, vertex_handle);
    if (indices.count(0) > 0) indices[0][vertex_handle].push_back(new_dart);
    return new_dart;
  }
  
  Dart_handle create_edge_and_index(typename LCC::template Attribute_handle<1>::type edge_handle, Dart_handle start, Dart_handle end) {
    lcc.template link_beta<1>(start, end);
    if (indices.count(1) > 0) {
      if (start->template attribute<0>()->point() < end->template attribute<0>()->point()) indices[1][start->template attribute<0>()].push_back(start);
      else indices[1][start->template attribute<0>()].push_back(end);
    } if (edge_handle != NULL) lcc.template set_attribute<1>(start, edge_handle);
    return start;
  }
  
  template <class InputIterator>
  Dart_handle create_facet_from_vertices_and_index(typename LCC::template Attribute_handle<2>::type facet_handle, InputIterator first, InputIterator last) {
    InputIterator current = first;
    InputIterator next = current;
    ++next;
    Dart_handle smallest = *current;
    while (next != last) {
      lcc.template link_beta<1>(*current, *next);
      //lcc.template set_attribute<1>(*current, lcc.template create_attribute<1>()); Don't create attributes here to avoid creating duplicates for (yet) unliked edges!!!
      if ((*current)->template attribute<0>()->point() < smallest->template attribute<0>()->point()) smallest = *current;
      ++current;
      ++next;
    } if ((*current)->template attribute<0>()->point() < smallest->template attribute<0>()->point()) smallest = *current;
    lcc.template link_beta<1>(*current, *first);
    //lcc.template set_attribute<1>(*current, lcc.template create_attribute<1>()); // Same as above!!!
    if (facet_handle != NULL) lcc.template set_attribute<2>(*current, facet_handle);
    if (indices.count(2) > 0) indices[2][smallest->template attribute<0>()].push_back(smallest);
    return *current;
  }
  
  // TODO: create_facet_from_edges_and_index
  template <class InputIterator>
  Dart_handle create_facet_from_edges_and_index(typename LCC::template Attribute_handle<2>::type facet_handle, InputIterator first, InputIterator last) {
    
  }
  
  template <unsigned int dimension, class InputIterator>
  Dart_handle create_cell_and_index(typename LCC::template Attribute_handle<dimension>::type cell_handle, InputIterator first, InputIterator last) {
    CGAL_precondition(dimension > 2);
    //std::cout << "create_cell_and_index<" << dimension << ">(...)" << std::endl;
    
    link_cells<dimension-1>(first, last);
    
    if (cell_handle != NULL) lcc.template set_attribute<dimension>(*first, cell_handle);
    if (indices.count(dimension) > 0) {
      Dart_handle smallest_vertex = find_smallest_vertex_in_cell<dimension>(*first);
      indices[dimension][smallest_vertex->template attribute<0>()].push_back(smallest_vertex);
    }
    
    return *first;
  }
  
  // -------------------------------------------------------------------------------------------------------------------
  // GET FUNCTIONS
  //
  // Used to get an i-cell from the given i-1 cells on its boundary (or a 0-attribute for 0-cells)
  // If there is no such existing i-cell, it is created and returned
  // If it already exists, the existing one is returned
  // Existing i-1 cells are used whenever possible
  //
  // The Dart_handle returned is for a Dart in the newly created cell
  // The bool indicates if the cell returned has been newly created
  // -------------------------------------------------------------------------------------------------------------------
  
  Dart_handle get_vertex(typename LCC::Vertex_attribute_handle vertex_handle, bool check_for_existing_cell = true) {
    std::cout << "get_vertex(" << vertex_handle->point() << ")" << std::endl;
    Dart_handle existing_vertex = LCC::null_handle;
    if (check_for_existing_cell) existing_vertex = find_vertex(vertex_handle);
    if (existing_vertex == LCC::null_handle) {
      Dart_handle new_dart = create_dart_for_vertex_and_index(vertex_handle);
      std::cout << "\tcreated new dart<" << &*new_dart << ">[" << new_dart->id << "] for vertex<" << &*vertex_handle << ">[" << vertex_handle->info() << "] (" << vertex_handle->point() << ")" << std::endl;
      return new_dart;
    } else {
      std::cout << "\tdeleted given vertex<" << &*vertex_handle << ">[" << vertex_handle->info() << "] and returning existing dart<" << &*existing_vertex << ">[" << existing_vertex->id << "] for vertex<" << &*existing_vertex->template attribute<0>() << ">[" << existing_vertex->template attribute<0>()->info() << "] (" << existing_vertex->template attribute<0>()->point() << ")" << std::endl;
      lcc.erase_vertex_attribute(vertex_handle);
      return existing_vertex;
    }
  }
  
  std::pair<Dart_handle, bool> get_edge(typename LCC::template Attribute_handle<1>::type edge_handle, Dart_handle start, Dart_handle end, bool check_for_existing_cell = true, bool attempt_to_reuse_cells = true) {
    std::cout << "get_edge(...)" << std::endl;
    if (check_for_existing_cell) {
      std::pair<Dart_handle, bool> existing_edge = find_edge(start, end, false);
      if (existing_edge.first != LCC::null_handle) {
        std::cout << "\treturning existing edge..." << std::endl;
        return existing_edge;
      }
    }
    
    Dart_handle free_start, free_end;
    if (lcc.is_free(start, 0) && lcc.is_free(start, 1)) {
      std::cout << "\tusing given dart<" << &*start << "> for start vertex<" << &*start->template attribute<0>() << ">..." << std::endl;
      free_start = start;
    } else if (attempt_to_reuse_cells) free_start = find_vertex(start->template attribute<0>(), true);
    if (free_start == LCC::null_handle) {
      free_start = lcc.create_dart(number_of_darts);
      ++number_of_darts;
      lcc.set_vertex_attribute(free_start, start->template attribute<0>());
      std::cout << "\tcreated new dart<" << &*free_start << "> for start vertex<" << &*start->template attribute<0>() << ">..." << std::endl;
    }
    
    if (lcc.is_free(end, 0) && lcc.is_free(end, 1)) {
      std::cout << "\tusing given dart<" << &*end << "> for end vertex<" << &*end->template attribute<0>() << ">..." << std::endl;
      free_end = end;
    } else if (attempt_to_reuse_cells) free_end = find_vertex(end->template attribute<0>(), true);
    if (free_end == LCC::null_handle) {
      free_end = lcc.create_dart(number_of_darts);
      ++number_of_darts;
      lcc.set_vertex_attribute(free_end, end->template attribute<0>());
      std::cout << "\tcreated new dart<" << &*free_end << "> for end vertex<" << &*end->template attribute<0>() << ">..." << std::endl;
    }
    
    Dart_handle new_edge = create_edge_and_index(edge_handle, free_start, free_end);
    std::cout << "\tcreated structure at dart<" << &*new_edge << "> for edge<" << &*new_edge->template attribute<1>() << ">" << std::endl;
    return std::pair<Dart_handle, bool>(new_edge, true);
  }
  
  template <class InputIterator>
  std::pair<Dart_handle, bool> get_facet_from_vertices(typename LCC::template Attribute_handle<2>::type facet_handle, InputIterator first, InputIterator last, bool check_for_existing_cell = true, bool attempt_to_reuse_cells = true) {
    std::cout << "get_facet_from_vertices(...)" << std::endl;
    if (check_for_existing_cell) {
      std::pair<Dart_handle, bool> existing_facet = find_facet_from_vertices(first, last, false);
      if (existing_facet.first != LCC::null_handle) {
        std::cout << "\treturning existing facet at dart<" << &*existing_facet.first << ">..." << std::endl;
        return existing_facet;
      }
    }
    
    InputIterator current = first;
    InputIterator next = current;
    ++next;
    std::list<Dart_handle> free_vertices;
    while (current != last) {
      if (next == last) next = first;
      if (lcc.is_free(*current, 0) && lcc.is_free(*current, 1)) {
        std::cout << "\tusing given dart<" << &**current << ">[" << (*current)->id << "] at vertex<" << &*(*current)->template attribute<0>() << ">[" << (*current)->template attribute<0>()->info() << "] (" << (*current)->template attribute<0>()->point() << ")" << std::endl;
        free_vertices.push_back(*current);
        ++current;
        ++next;
        continue;
      } if (!attempt_to_reuse_cells) {
        free_vertices.push_back(lcc.create_dart(number_of_darts));
        ++number_of_darts;
        lcc.set_vertex_attribute(free_vertices.back(), (*current)->template attribute<0>());
        std::cout << "\tcreated new dart<" << &*free_vertices.back() << ">[" << free_vertices.back()->id << "] at vertex<" << &*(*current)->template attribute<0>() << ">[" << (*current)->template attribute<0>()->info() << "] (" << (*current)->template attribute<0>()->point() << ")" << std::endl;
        //CGAL_assertion(lcc.is_valid()); NOTE: It might not be valid since it might not be manifold!
        ++current;
        ++next;
        continue;
      } std::pair<Dart_handle, bool> free_edge = find_edge(*current, *next, true);
      if (free_edge.first != LCC::null_handle) {
        std::cout << "\treusing existing edge..." << std::endl;
        if (free_edge.second) {
          free_vertices.push_back(free_edge.first);
          free_vertices.push_back(free_edge.first->beta(1));
          lcc.template unlink_beta<1>(free_edge.first);
          ++current;
          ++next;
          if (next == last) break;
        } else {
          free_vertices.push_back(free_edge.first->beta(1));
          free_vertices.push_back(free_edge.first);
          lcc.template unlink_beta<1>(free_edge.first);
          ++current;
          ++next;
          if (next == last) break;
        } ++current;
        ++next;
        continue;
      }
      Dart_handle free_vertex = find_vertex((*current)->template attribute<0>(), true);
      if (free_vertex != LCC::null_handle) {
        free_vertices.push_back(free_vertex);
        std::cout << "\reusing free dart<" << &*free_vertices.back() << ">[" << free_vertices.back()->id << "] at vertex<" << &*(*current)->template attribute<0>() << ">[" << (*current)->template attribute<0>()->info() << "] (" << (*current)->template attribute<0>()->point() << ")" << std::endl;
      } else {
        free_vertices.push_back(lcc.create_dart(number_of_darts));
        ++number_of_darts;
        lcc.set_vertex_attribute(free_vertices.back(), (*current)->template attribute<0>());
        std::cout << "\tcreated new dart<" << &*free_vertices.back() << ">[" << free_vertices.back()->id << "] at vertex<" << &*(*current)->template attribute<0>() << ">[" << (*current)->template attribute<0>()->info() << "] (" << (*current)->template attribute<0>()->point() << ")" << std::endl;
      } ++current;
      ++next;
    }
    
    Dart_handle new_facet = create_facet_from_vertices_and_index(facet_handle, free_vertices.begin(), free_vertices.end());
    std::cout << "\tcreated structure at dart<" << &*new_facet << ">[" << new_facet->id << "] for facet<" << &*facet_handle << ">[" << facet_handle->info() << "]" << std::endl;
    return std::pair<Dart_handle, bool>(new_facet, true);
  }
  
  template <class InputIterator>
  std::pair<Dart_handle, bool> get_facet_from_edges(typename LCC::template Attribute_handle<2>::type facet_handle, InputIterator first, InputIterator last, bool check_for_existing_cell = true, bool attempt_to_reuse_cells = true) {
    
  }
  
  template <unsigned int dimension, class InputIterator>
  std::pair<Dart_handle, bool> get_cell(typename LCC::template Attribute_handle<dimension>::type cell_handle, InputIterator first, InputIterator last, bool check_if_it_exists = true, bool attempt_to_reuse_cells = true) {
    std::cout << "get_cell<" << dimension << ">(...)" << std::endl;
    CGAL_precondition(dimension > 2);
    
    // Reuse given or existing cells, if possible
    std::list<Dart_handle> free_cells;
    for (InputIterator current = first; current != last; ++current) {
      if (lcc.is_free(*current, dimension-1)) {
        std::cout << "\tusing given dart<" << &**current << ">[" << (*current)->id << "] at " << dimension-1 << "-cell<" << &*(*current)->template attribute<2>() << ">[" << (*current)->template attribute<2>()->info() << "]" << std::endl;
        free_cells.push_back(*current);
        continue;
      } else if (!attempt_to_reuse_cells) {
        Dart_handle free_cell = clone_and_link_flipped_cell<dimension-1>(*current);
        std::cout << "\tcreated new dart<" << &*free_cell << ">[" << free_cell->id << "] at " << dimension-1 << "-cell<" << &*free_cell->template attribute<dimension-1>() << ">[" << free_cell->template attribute<dimension-1>()->info() << "]" << std::endl;
        free_cells.push_back(free_cell);
        continue;
      }
      
      // TODO: reusing existing cells
      
      Dart_handle free_cell = clone_and_link_flipped_cell<dimension-1>(*current);
      std::cout << "\tcreated new dart<" << &*free_cell << ">[" << free_cell->id << "] at " << dimension-1 << "-cell<" << &*free_cell->template attribute<dimension-1>() << ">[" << free_cell->template attribute<dimension-1>()->info() << "]" << std::endl;
      free_cells.push_back(free_cell);
      continue;
    }
    
    Dart_handle new_cell = create_cell_and_index<dimension>(cell_handle, free_cells.begin(), free_cells.end());
    std::cout << "\tcreated structure at dart<" << &*new_cell << ">[" << new_cell->id << "] for " << dimension << "-cell<" << &*new_cell->template attribute<dimension>() << ">[" << new_cell->template attribute<dimension>()->info() << "]" << std::endl;
    return std::pair<Dart_handle, bool>(new_cell, true);
  }
};

#endif
