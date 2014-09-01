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

#ifndef LCC_ORBIT_COMPARATOR_H
#define LCC_ORBIT_COMPARATOR_H

#include <CGAL/Linear_cell_complex.h>

template <class LCC1, class LCC2 = LCC1>
class Linear_cell_complex_cell_comparator {
public:
  template <unsigned int dimension>
  static bool have_same_geometry(LCC1 &lcc, typename LCC1::Dart_handle dh1, typename LCC2::Dart_handle dh2) {
    return Linear_cell_complex_cell_comparator::template have_same_geometry<dimension>(lcc, dh1, lcc, dh2);
  }
  
  template <unsigned int dimension>
  static bool have_same_geometry(LCC1 &lcc1, typename LCC1::Dart_handle dh1, LCC2 &lcc2, typename LCC2::Dart_handle dh2) {
    for (typename LCC1::template Dart_of_cell_range<dimension, dimension>::iterator current_dart1 = lcc1.template darts_of_cell<dimension, dimension>(dh1).begin(); current_dart1 != lcc1.template darts_of_cell<dimension, dimension>(dh1).end(); ++current_dart1) {
      if (Linear_cell_complex_cell_comparator::template have_same_geometry_and_orientation_from<dimension>(lcc1, current_dart1, lcc2, dh2)) return true;
      if (Linear_cell_complex_cell_comparator::template have_same_geometry_and_opposite_orientation_from<dimension>(lcc1, current_dart1, lcc2, dh2)) return true;
    } return false;
  }
  
  template <unsigned int dimension>
  static bool have_same_geometry_and_orientation_from(LCC1 &lcc, typename LCC1::Dart_handle dh1, typename LCC2::Dart_handle dh2) {
    return Linear_cell_complex_cell_comparator::template have_same_geometry_and_orientation_from<dimension>(lcc, dh1, lcc, dh2);
  }
  
  template <unsigned int dimension>
  static bool have_same_geometry_and_orientation_from(LCC1 &lcc1, typename LCC1::Dart_handle dh1, LCC2 &lcc2, typename LCC2::Dart_handle dh2) {
    CGAL_precondition(dimension <= LCC1::dimension+1);
    CGAL_precondition(dimension <= LCC2::dimension+1);
    if (dh1 == dh2) return true;
    
    if (dimension == 0) return lcc1.point(dh1) == lcc2.point(dh2);
    if (dimension == 1) {
      if (lcc1.point(dh1) != lcc2.point(dh2)) return false;
      if (lcc1.point(dh1->beta(1)) != lcc2.point(dh2->beta(1))) return false;
      return true;
    }
    
    bool match = true;
    
    int mark1 = lcc1.get_new_mark();
    int mark2 = lcc2.get_new_mark();
    CGAL_precondition(lcc1.is_whole_map_unmarked(mark1));
    CGAL_precondition(lcc2.is_whole_map_unmarked(mark2));
    std::stack<typename LCC1::Dart_handle> to_check1, to_unmark1;
    std::stack<typename LCC2::Dart_handle> to_check2, to_unmark2;
    
    to_check1.push(dh1);
    to_check2.push(dh2);
    lcc1.mark(dh1, mark1);
    lcc2.mark(dh2, mark2);
    to_unmark1.push(dh1);
    to_unmark2.push(dh2);
    
    typename LCC1::Dart_handle dart1;
    typename LCC2::Dart_handle dart2;
    std::map<typename LCC1::Dart_handle, typename LCC2::Dart_handle> bijection; // TODO: Slow
    
    while (match && !to_check1.empty()) {
      dart1 = to_check1.top();
      to_check1.pop();
      dart2 = to_check2.top();
      to_check2.pop();
      
      // Check vertices
      if (lcc1.point(dart1) != lcc2.point(dart2)) {
        match = false;
        break;
      }
      
      bijection[dart1] = dart2;
      
      for (unsigned int current_dimension = 1; current_dimension < dimension; ++current_dimension) {
        // Check betas
        if (lcc1.is_free(dart1, current_dimension)) {
          if (!lcc2.is_free(dart2, current_dimension)) {
            match = false;
            break;
          } continue;
        } if (lcc2.is_free(dart2, current_dimension)) {
          match = false;
          break;
        }
        
        // Check neighbours
        if (lcc1.is_marked(dart1->beta(current_dimension), mark1)) {
          if (!lcc2.is_marked(dart2->beta(current_dimension), mark2)) {
            match = false;
            break;
          } if (bijection.count(dart1->beta(current_dimension)) && bijection.count(dart2->beta(current_dimension))) {
            if (bijection[dart1->beta(current_dimension)] != dart2->beta(current_dimension)) {
              match = false;
              break;
            }
          } continue;
        } if (lcc2.is_marked(dart2->beta(current_dimension), mark2)) {
          match = false;
          break;
        }
        
        to_check1.push(dart1->beta(current_dimension));
        to_check2.push(dart2->beta(current_dimension));
        lcc1.mark(dart1->beta(current_dimension), mark1);
        lcc2.mark(dart2->beta(current_dimension), mark2);
        to_unmark1.push(dart1->beta(current_dimension));
        to_unmark2.push(dart2->beta(current_dimension));
      }
    }
    
    // Unmark darts
    CGAL_assertion(to_unmark1.size() == to_unmark2.size());
    while (!to_unmark1.empty()) {
      CGAL_assertion(lcc1.is_marked(to_unmark1.top(), mark1));
      lcc1.unmark(to_unmark1.top(), mark1);
      to_unmark1.pop();
      CGAL_assertion(lcc2.is_marked(to_unmark2.top(), mark2));
      lcc2.unmark(to_unmark2.top(), mark2);
      to_unmark2.pop();
    }
    
    // Free marks
    CGAL_postcondition(lcc1.is_whole_map_unmarked(mark1));
    CGAL_postcondition(lcc2.is_whole_map_unmarked(mark2));
    lcc1.free_mark(mark1);
    lcc2.free_mark(mark2);
    
    return match;
  }
  
  template <unsigned int dimension>
  static bool have_same_geometry_and_opposite_orientation_from(LCC1 &lcc, typename LCC1::Dart_handle dh1, typename LCC2::Dart_handle dh2) {
    return Linear_cell_complex_cell_comparator::template have_same_geometry_and_opposite_orientation_from<dimension>(lcc, dh1, lcc, dh2);
  }
  
  template <unsigned int dimension>
  static bool have_same_geometry_and_opposite_orientation_from(LCC1 &lcc1, typename LCC1::Dart_handle dh1, LCC2 &lcc2, typename LCC2::Dart_handle dh2) {
    CGAL_precondition(dimension <= LCC1::dimension+1);
    CGAL_precondition(dimension <= LCC2::dimension+1);
    if (dh1 == dh2) return false;
    
    if (dimension == 0) return lcc1.point(dh1) == lcc2.point(dh2);
    if (dimension == 1) {
      if (lcc1.point(dh1) != lcc2.point(dh2->beta(1))) return false;
      if (lcc1.point(dh1->beta(1)) != lcc2.point(dh2)) return false;
      return true;
    }
    
    bool match = true;
    
    int mark1 = lcc1.get_new_mark();
    int mark2 = lcc2.get_new_mark();
    CGAL_precondition(lcc1.is_whole_map_unmarked(mark1));
    CGAL_precondition(lcc2.is_whole_map_unmarked(mark2));
    std::stack<typename LCC1::Dart_handle> to_check1, to_unmark1;
    std::stack<typename LCC2::Dart_handle> to_check2, to_unmark2;
    
    to_check1.push(dh1);
    to_check2.push(dh2);
    lcc1.mark(dh1, mark1);
    lcc2.mark(dh2, mark2);
    to_unmark1.push(dh1);
    to_unmark2.push(dh2);
    
    typename LCC1::Dart_handle dart1;
    typename LCC2::Dart_handle dart2;
    std::map<typename LCC1::Dart_handle, typename LCC2::Dart_handle> bijection; // TODO: Slow
    
    while (match && !to_check1.empty()) {
      dart1 = to_check1.top();
      to_check1.pop();
      dart2 = to_check2.top();
      to_check2.pop();
      
      // Check vertices
      if (!lcc1.is_free(dart1, 0) && !lcc2.is_free(dart2, 1)) {
        if (lcc1.point(dart1) != lcc2.point(dart2->beta(1))) {
          match = false;
          break;
        }
      }
      
      bijection[dart1] = dart2;
      
      for (unsigned int current_dimension = 1; current_dimension < dimension; ++current_dimension) {
        // Check betas
        if (lcc1.is_free(dart1, current_dimension)) {
          if (!lcc2.is_free(dart2, CGAL_BETAINV(current_dimension))) {
            match = false;
            break;
          } continue;
        } if (lcc2.is_free(dart2, CGAL_BETAINV(current_dimension))) {
          match = false;
          break;
        }
        
        // Check neighbours
        if (lcc1.is_marked(dart1->beta(current_dimension), mark1)) {
          if (!lcc2.is_marked(dart2->beta(CGAL_BETAINV(current_dimension)), mark2)) {
            match = false;
            break;
          } if (bijection.count(dart1->beta(current_dimension)) && bijection.count(dart2->beta(CGAL_BETAINV(current_dimension)))) {
            if (bijection[dart1->beta(current_dimension)] != dart2->beta(CGAL_BETAINV(current_dimension))) {
              match = false;
              break;
            }
          } continue;
        } if (lcc2.is_marked(dart2->beta(CGAL_BETAINV(current_dimension)), mark2)) {
          match = false;
          break;
        }
        
        to_check1.push(dart1->beta(current_dimension));
        to_check2.push(dart2->beta(CGAL_BETAINV(current_dimension)));
        lcc1.mark(dart1->beta(current_dimension), mark1);
        lcc2.mark(dart2->beta(CGAL_BETAINV(current_dimension)), mark2);
        to_unmark1.push(dart1->beta(current_dimension));
        to_unmark2.push(dart2->beta(CGAL_BETAINV(current_dimension)));
      }
    }
    
    // Unmark darts
    CGAL_assertion(to_unmark1.size() == to_unmark2.size());
    while (!to_unmark1.empty()) {
      lcc1.unmark(to_unmark1.top(), mark1);
      to_unmark1.pop();
      lcc2.unmark(to_unmark2.top(), mark2);
      to_unmark2.pop();
    }
    
    // Free marks
    CGAL_postcondition(lcc1.is_whole_map_unmarked(mark1));
    CGAL_postcondition(lcc2.is_whole_map_unmarked(mark2));
    lcc1.free_mark(mark1);
    lcc2.free_mark(mark2);
    
    return match;
  }
};

#endif