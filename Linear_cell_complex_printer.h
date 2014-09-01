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

#ifndef LCC_PRINTER_H
#define LCC_PRINTER_H

#include <CGAL/Linear_cell_complex.h>

#include "Linear_cell_complex_with_ids.h"

template <class LCC, unsigned int dimension = LCC::dimension>
struct Linear_cell_complex_attribute_printer;

template <class LCC, unsigned int dimension>
struct Linear_cell_complex_attribute_printer {
  static void print(LCC &lcc) {
    Linear_cell_complex_attribute_printer<LCC, dimension-1>::print(lcc);
    typename LCC::template Attribute_const_range<dimension>::type &cells = lcc.template attributes<dimension>();
    std::cout << dimension << "-attributes (" << cells.size() << "):" << std::endl;
    for (typename LCC::template Attribute_const_range<dimension>::type::const_iterator current_cell = cells.begin(); current_cell != cells.end(); ++current_cell) {
      std::cout << "\t" << dimension << "-attribute<" << &*current_cell << ">[" << current_cell->info() << "] -> ";
      if (&*current_cell->dart() != 0x00) std::cout << "dart<" << &*current_cell->dart() << ">[" << current_cell->dart()->id << "]" << std::endl;
      else std::cout << "null" << std::endl;
    }
  }
};

template <class LCC>
struct Linear_cell_complex_attribute_printer<LCC, 0> {
  static void print(LCC &lcc) {
    typename LCC::Vertex_attribute_const_range &vertices = lcc.vertex_attributes();
    std::cout << "0-attributes (" << vertices.size() << "):" << std::endl;
    for (typename LCC::Vertex_attribute_const_range::const_iterator current_vertex = vertices.begin(); current_vertex != vertices.end(); ++current_vertex) {
      std::cout << "\t0-attribute<" << &*current_vertex << ">[" << current_vertex->info() << "] (" << current_vertex->point() << ") -> ";
      if (&*current_vertex->dart() != 0x00) std::cout << "dart<" << &*current_vertex->dart() << ">[" << current_vertex->dart()->id << "]" << std::endl;
      else std::cout << "null" << std::endl;
    }
  }
};

template <class LCC, unsigned int dimension = LCC::dimension>
struct Linear_cell_complex_dart_attribute_printer;

template <class LCC, unsigned int dimension>
struct Linear_cell_complex_dart_attribute_printer {
  static void print(typename LCC::Dart_const_handle &dart) {
    Linear_cell_complex_dart_attribute_printer<LCC, dimension-1>::print(dart);
    if (dart->template attribute<dimension>() != NULL) std::cout << "\t\t" << dimension << "-attribute<" << &*dart->template attribute<dimension>() << ">[" << dart->template attribute<dimension>()->info() << "]" << std::endl;
    else std::cout << "\t\t" << dimension << "-attribute<NULL>" << std::endl;
  }
};

template <class LCC>
struct Linear_cell_complex_dart_attribute_printer<LCC, 0> {
  static void print(typename LCC::Dart_const_handle &dart) {
    std::cout << "\t\t0-attribute<" << &*dart->template attribute<0>() << ">[" << dart->template attribute<0>()->info() << "] (" << dart->template attribute<0>()->point() << ")" << std::endl;
  }
};

template <class LCC, unsigned int dimension = LCC::dimension>
struct Linear_cell_complex_dart_beta_printer;

template <class LCC, unsigned int dimension>
struct Linear_cell_complex_dart_beta_printer {
  static void print(typename LCC::Dart_const_handle &dart) {
    Linear_cell_complex_dart_beta_printer<LCC, dimension-1>::print(dart);
    if (dart->beta(dimension) != LCC::null_handle) std::cout << "\t\tbeta[" << dimension << "]: dart<" << &*dart->beta(dimension) << ">[" << dart->beta(dimension)->id << "]" << std::endl;
    else std::cout << "\t\tbeta[" << dimension << "]: null" << std::endl;
  }
};

template <class LCC>
struct Linear_cell_complex_dart_beta_printer<LCC, 0> {
  static void print(typename LCC::Dart_const_handle &dart) {
    if (dart->beta(0) != LCC::null_handle) std::cout << "\t\tbeta[0]: dart<" << &*dart->beta(0) << ">[" << dart->beta(0)->id << "]" << std::endl;
    else std::cout << "\t\tbeta[0]: null" << std::endl;
  }
};

template <class LCC>
class Linear_cell_complex_printer {
public:
  static void print(LCC &lcc) {
    Linear_cell_complex_attribute_printer<LCC>::print(lcc);
    
    typename LCC::Dart_const_range &darts = lcc.darts();
    std::cout << "Darts (" << darts.size() << "):" << std::endl;
    for (typename LCC::Dart_const_range::const_iterator current_dart = darts.begin(); current_dart != darts.end(); ++current_dart) {
      std::cout << "\tdart<" << &*current_dart << ">[" << current_dart->id << "]" << std::endl;
      Linear_cell_complex_dart_attribute_printer<LCC>::print(current_dart);
      Linear_cell_complex_dart_beta_printer<LCC>::print(current_dart);
    }
  }
};

#endif
