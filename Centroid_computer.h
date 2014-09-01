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

#ifndef CENTROID_COMPUTER_H
#define CENTROID_COMPUTER_H

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

template <class LCC, class Point, class FT>
struct Centroid_computer {
  template <unsigned int dimension>
  static Point centroid_of_cell(LCC &lcc, typename LCC::Dart_const_handle dart) {
    std::cerr << "Error: Cannot compute centroid of unknown type." << std::endl;
    return Point();
  }
};

template <class LCC, class Kernel>
struct Centroid_computer<LCC, typename CGAL::Point_2<Kernel>, typename Kernel::FT> {
  template <unsigned int dimension>
  static typename LCC::Point centroid_of_cell(LCC &lcc, typename LCC::Dart_const_handle dart) {
    typename Kernel::FT sum[2];
    unsigned int number_of_points = 0;
    for (unsigned int current_dimension = 0; current_dimension < 2; ++current_dimension) sum[current_dimension] = 0;
    for (typename LCC::template Dart_of_cell_const_range<dimension, dimension>::const_iterator current_dart = lcc.template darts_of_cell<dimension, dimension>(dart).begin(); current_dart != lcc.template darts_of_cell<dimension, dimension>(dart).end(); ++current_dart) {
      if (current_dart->template attribute<0>() == LCC::null_handle) continue;
      for (unsigned int current_dimension = 0; current_dimension < 2; ++current_dimension) {
        sum[current_dimension] += lcc.point(current_dart)[current_dimension];
      } ++number_of_points;
    } for (unsigned int current_dimension = 0; current_dimension < 2; ++current_dimension) {
      sum[current_dimension] /= number_of_points;
    } return CGAL::Point_2<Kernel>(sum[0], sum[1]);
  }
};

template <class LCC, class Kernel>
struct Centroid_computer<LCC, typename CGAL::Point_3<Kernel>, typename Kernel::FT> {
  template <unsigned int dimension>
  static typename LCC::Point centroid_of_cell(LCC &lcc, typename LCC::Dart_const_handle dart) {
    typename Kernel::FT sum[3];
    unsigned int number_of_points = 0;
    for (unsigned int current_dimension = 0; current_dimension < 3; ++current_dimension) sum[current_dimension] = 0;
    for (typename LCC::template Dart_of_cell_const_range<dimension, dimension>::const_iterator current_dart = lcc.template darts_of_cell<dimension, dimension>(dart).begin(); current_dart != lcc.template darts_of_cell<dimension, dimension>(dart).end(); ++current_dart) {
      if (current_dart->template attribute<0>() == LCC::null_handle) continue;
      for (unsigned int current_dimension = 0; current_dimension < 3; ++current_dimension) {
        sum[current_dimension] += lcc.point(current_dart)[current_dimension];
      } ++number_of_points;
    } for (unsigned int current_dimension = 0; current_dimension < 3; ++current_dimension) {
      sum[current_dimension] /= number_of_points;
    } //std::cout << "Centroid of " << dimension << "-cell<" << &*dart->template attribute<dimension>() << "> at " << CGAL::Point_3<Kernel>(sum[0], sum[1], sum[2]) << std::endl;
    return CGAL::Point_3<Kernel>(sum[0], sum[1], sum[2]);
  }
};

template <class LCC, class Kernel>
struct Centroid_computer<LCC, typename CGAL::Point_d<Kernel>, typename Kernel::FT> {
  template <unsigned int dimension>
  static typename LCC::Point centroid_of_cell(LCC &lcc, typename LCC::Dart_const_handle dart) {
    typename Kernel::FT sum[LCC::d2];
    unsigned int number_of_points = 0;
    for (unsigned int current_dimension = 0; current_dimension < LCC::FT2; ++current_dimension) sum[current_dimension] = 0;
    for (typename LCC::template Dart_of_cell_const_range<dimension, dimension>::const_iterator current_dart = lcc.template darts_of_cell<dimension, dimension>(dart).begin(); current_dart != lcc.template darts_of_cell<dimension, dimension>(dart).end(); ++current_dart) {
      for (unsigned int current_dimension = 0; current_dimension < LCC::FT2; ++current_dimension) {
        sum[current_dimension] += lcc.point(current_dart)[current_dimension];
      } ++number_of_points;
    } for (unsigned int current_dimension = 0; current_dimension < LCC::FT2; ++current_dimension) {
      sum[current_dimension] /= number_of_points;
    } return CGAL::Point_d<Kernel>(LCC::d2, sum, sum+LCC::d2);
  }
};

#endif
