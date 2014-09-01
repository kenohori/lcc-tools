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

#ifndef POINT_FITTER_H
#define POINT_FITTER_H

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

template <class LCC, class Point, class FT>
struct Point_fitter {
  Point_fitter() {
    std::cerr << "Error: Cannot fit point of unknown type." << std::endl;
  }
  
  void add_point(const Point &p) {
    std::cerr << "Error: Cannot fit point of unknown type." << std::endl;
  }
  
  Point get_point(const Point &p) {
    std::cerr << "Error: Cannot fit point of unknown type." << std::endl;
    return Point();
  }
};

template <class LCC, class Kernel>
struct Point_fitter<LCC, typename CGAL::Point_2<Kernel>, typename Kernel::FT> {
  typename CGAL::Point_2<Kernel> min, max;
  unsigned int num_points = 0;
  
  Point_fitter() {
    
  }
  
  void add_point(const typename CGAL::Point_2<Kernel> &p) {
    if (num_points == 0) {
      min = max = p;
      ++num_points;
    } else {
      if (p < min) min = p;
      if (p > max) max = p;
      ++num_points;
    }
  }
  
  typename CGAL::Point_2<Kernel> get_point(const typename CGAL::Point_2<Kernel> &p) {
    typename CGAL::Vector_2<Kernel> ranges = max - min;
    typename Kernel::FT scale_factor;
    if (ranges.y() > ranges.x()) scale_factor = 10.0/ranges[1];
    else scale_factor = scale_factor = 10.0/ranges[0];
    typename CGAL::Vector_2<Kernel> translation_to_centre = CGAL::ORIGIN - CGAL::midpoint(min, max);
    return CGAL::ORIGIN+((p+translation_to_centre)-CGAL::ORIGIN)*scale_factor;
  }
};

template <class LCC, class Kernel>
struct Point_fitter<LCC, typename CGAL::Point_3<Kernel>, typename Kernel::FT> {
  typename CGAL::Point_3<Kernel> min, max;
  unsigned int num_points = 0;
  
  Point_fitter() {
    
  }
  
  void add_point(const typename CGAL::Point_3<Kernel> &p) {
    if (num_points == 0) {
      min = max = p;
      ++num_points;
    } else {
      if (p < min) min = p;
      if (p > max) max = p;
      ++num_points;
    }
  }
  
  typename CGAL::Point_3<Kernel> get_point(const typename CGAL::Point_3<Kernel> &p) {
    typename CGAL::Vector_3<Kernel> ranges = max - min;
    typename Kernel::FT scale_factor;
    if (ranges.z() > ranges.x() && ranges.z() > ranges.y()) scale_factor = 10.0/ranges[2];
    else if (ranges.y() > ranges.x()) scale_factor = 10.0/ranges[1];
    else scale_factor = scale_factor = 10.0/ranges[0];
    typename CGAL::Vector_3<Kernel> translation_to_centre = CGAL::ORIGIN - CGAL::midpoint(min, max);
    return CGAL::ORIGIN+((p+translation_to_centre)-CGAL::ORIGIN)*scale_factor;
  }
};

#endif
