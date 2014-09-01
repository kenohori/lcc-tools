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

#ifndef LCC_WITH_IDS_H
#define LCC_WITH_IDS_H

#include <set>

#include <CGAL/Linear_cell_complex.h>

template <int d, class Refs>
struct Dart_with_id : public CGAL::Dart<d, Refs> {
public:
  typedef CGAL::Dart<d, Refs> Dart;
  typedef typename Refs::size_type size_type;
  static const size_type NB_MARKS = Refs::NB_MARKS;
  int id;
  
  Dart_with_id() : Dart() {
    id = -1;
  }
  
  Dart_with_id(int id) : Dart() {
    this->id = id;
  }
  
  Dart_with_id(const Dart& adart) : Dart(adart) {
    id = -1;
  }
};

template <unsigned int d>
struct Linear_cell_complex_items_with_id {
  template <class LCC>
  struct Dart_wrapper {
    typedef CGAL::Cell_attribute_with_point<LCC, int> Point_attribute_with_id;
    typedef CGAL::Cell_attribute<LCC, int> Attribute_with_id;
    
    template <unsigned int attributes_to_add, class Result = CGAL::cpp11::tuple<> >
    struct Linear_cell_complex_items_with_id_attributes;
    
    template <class ... Result>
    struct Linear_cell_complex_items_with_id_attributes<0, CGAL::cpp11::tuple<Result ...> > {
      typedef CGAL::cpp11::tuple<Point_attribute_with_id, Result ...> tuple;
    };
    
    template <unsigned int attributes_to_add, class ... Result>
    struct Linear_cell_complex_items_with_id_attributes<attributes_to_add, CGAL::cpp11::tuple<Result ...> > {
      typedef typename Linear_cell_complex_items_with_id_attributes<attributes_to_add-1, CGAL::cpp11::tuple<Attribute_with_id, Result ...> >::tuple tuple;
    };
    
    typedef Dart_with_id<d, LCC> Dart;
    typedef typename Linear_cell_complex_items_with_id_attributes<d>::tuple Attributes;
  };
};

template <unsigned int d>
class Linear_cell_complex_with_ids : public CGAL::Linear_cell_complex<d, d, CGAL::Linear_cell_complex_traits<d>, Linear_cell_complex_items_with_id<d> > {
public:
  typedef CGAL::Linear_cell_complex<d, d, CGAL::Linear_cell_complex_traits<d>, Linear_cell_complex_items_with_id<d> > LCC;
  typedef typename LCC::Dart_handle Dart_handle;
  typedef typename LCC::Point Point;
  typedef typename LCC::Vertex_attribute_handle Vertex_attribute_handle;
  template<int i>
  struct Attribute_handle : public LCC::template Attribute_handle<i> {};
  
  Vertex_attribute_handle create_vertex_attribute(const Point &apoint) {
    //    std::cout << "Creating vertex[" << this->vertex_attributes().size() << "]" << std::endl;
    return LCC::create_vertex_attribute(apoint, this->vertex_attributes().size());
  }
  
  template <unsigned int i>
  typename Attribute_handle<i>::type create_attribute() {
    //    std::cout << "Creating " << i << "-attribute[" << this->template attributes<i>().size() << "]" << std::endl;
    return LCC::template create_attribute<i>(this->template attributes<i>().size());
  }
  
  Dart_handle create_dart() {
    Dart_handle new_dart = LCC::create_dart();
    new_dart->id = this->darts().size()-1;
    return new_dart;
  }
  
  Dart_handle create_dart(const Point &apoint) {
    Vertex_attribute_handle vah = create_vertex_attribute(apoint);
    Dart_handle new_dart = LCC::create_dart(vah);
    new_dart->id = this->darts().size()-1;
    return new_dart;
  }
  
  Dart_handle create_dart(Vertex_attribute_handle vh) {
    Dart_handle new_dart = LCC::create_dart(vh);
    new_dart->id = this->darts().size()-1;
    return new_dart;
  }
  
  Dart_handle create_temporary_dart() {
    Dart_handle new_dart = LCC::create_dart();
    CGAL_postcondition(new_dart->id == -1);
    return new_dart;
  }
  
  Dart_handle create_temporary_dart(Vertex_attribute_handle vh) {
    Dart_handle dart_of_vertex = vh->dart();
    CGAL_precondition(dart_of_vertex != LCC::null_handle);
    Dart_handle new_dart = LCC::create_dart(vh);
    vh->set_dart(dart_of_vertex);
    CGAL_postcondition(new_dart->id == -1);
    CGAL_postcondition(this->vertex_attribute(new_dart) == vh);
    CGAL_postcondition(vh->dart() != new_dart);
    return new_dart;
  }
  
  int make_temporary_darts_permanent() {
    int changed_darts = 0;
    int highest_id = -1;
    for (typename LCC::Dart_range::iterator current_dart = this->darts().begin(); current_dart != this->darts().end(); ++current_dart) {
      if (current_dart->id > highest_id) highest_id = current_dart->id;
    } for (typename LCC::Dart_range::iterator current_dart = this->darts().begin(); current_dart != this->darts().end(); ++current_dart) {
      if (current_dart->id == -1) {
        current_dart->id = highest_id+1;
        ++changed_darts;
        ++highest_id;
      }
    } return changed_darts;
  }
};

#endif