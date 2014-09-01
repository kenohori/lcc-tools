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

#ifndef LCC_CLONER_H
#define LCC_CLONER_H

#include <map>

#include <CGAL/Linear_cell_complex.h>

template <class LCC, unsigned int copy_dimension>
struct Copy_flipped_attributes {
  template <unsigned int dimension>
  static void run(LCC *lcc, typename LCC::Dart_handle dart_from, typename LCC::Dart_handle dart_to) {
    if (dimension == 0) lcc->template set_attribute<0>(dart_to, dart_from->beta(1)->template attribute<0>());
    else if (dimension != copy_dimension+1 && dart_from->template attribute<dimension>() != NULL) lcc->template set_attribute<dimension>(dart_to, dart_from->template attribute<dimension>());
  }
};

template <class LCC_>
class Linear_cell_complex_cell_cloner {
public:
  typedef LCC_ LCC;
  typedef typename LCC::Dart_handle Dart_handle;
  typedef typename LCC::Helper Helper;
  
  template <unsigned int dimension>
  static Dart_handle clone_cell_with_reversed_orientation(LCC &lcc, Dart_handle &dart_in_cell) {
    CGAL_precondition(dimension > 1);
    CGAL_precondition(LCC::dimension >= dimension);
    
    // Create a copy of each dart
    std::map<Dart_handle, Dart_handle> dart_map;
    for (typename LCC::template Dart_of_cell_range<dimension, dimension>::iterator current_dart = lcc.template darts_of_cell<dimension, dimension>(dart_in_cell).begin(); current_dart != lcc.template darts_of_cell<dimension, dimension>(dart_in_cell).end(); ++current_dart) {
      dart_map[current_dart] = lcc.create_temporary_dart();
    }
    
    // Link them according to the combinatorial structure of the cell of dart_in_cell
    for (typename LCC::template Dart_of_cell_range<dimension, dimension>::iterator current_dart = lcc.template darts_of_cell<dimension, dimension>(dart_in_cell).begin(); current_dart != lcc.template darts_of_cell<dimension, dimension>(dart_in_cell).end(); ++current_dart) {
      lcc.basic_link_beta_1(dart_map[current_dart], dart_map[current_dart->beta(0)]);
      for (unsigned int current_dimension = 2; current_dimension < dimension; ++current_dimension) {
        lcc.basic_link_beta_for_involution(dart_map[current_dart], dart_map[current_dart->beta(current_dimension)], current_dimension);
      } Helper::template Foreach_enabled_attributes<Copy_flipped_attributes<LCC, dimension> >::run(&lcc, current_dart, dart_map[current_dart]);
    } return dart_map[dart_in_cell];
  }
};

#endif
