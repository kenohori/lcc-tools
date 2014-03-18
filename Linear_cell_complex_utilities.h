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

#ifndef LCC_UTILITIES_H
#define LCC_UTILITIES_H

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
struct Linear_cell_complex_items_with_label {
  template <class LCC>
  struct Dart_wrapper {
    typedef CGAL::Cell_attribute_with_point<LCC, std::string> Point_attribute_with_label;
    typedef CGAL::Cell_attribute<LCC, std::string> Attribute_with_label;
    
    template <unsigned int attributes_to_add, class Result = CGAL::cpp11::tuple<> >
    struct Linear_cell_complex_items_with_label_attributes;
    
    template <class ... Result>
    struct Linear_cell_complex_items_with_label_attributes<0, CGAL::cpp11::tuple<Result ...> > {
      typedef CGAL::cpp11::tuple<Point_attribute_with_label, Result ...> tuple;
    };
    
    template <unsigned int attributes_to_add, class ... Result>
    struct Linear_cell_complex_items_with_label_attributes<attributes_to_add, CGAL::cpp11::tuple<Result ...> > {
      typedef typename Linear_cell_complex_items_with_label_attributes<attributes_to_add-1, CGAL::cpp11::tuple<Attribute_with_label, Result ...> >::tuple tuple;
    };
    
    typedef Dart_with_id<d, LCC> Dart;
    typedef typename Linear_cell_complex_items_with_label_attributes<d>::tuple Attributes;
  };
};

template <unsigned int d>
struct Linear_cell_complex_with_ids {
public:
  typedef CGAL::Linear_cell_complex<d, d, CGAL::Linear_cell_complex_traits<d>, Linear_cell_complex_items_with_id<d> > type;
};

template <unsigned int d>
struct Linear_cell_complex_with_labelled_embeddings_and_ids {
public:
  typedef CGAL::Linear_cell_complex<d, d, CGAL::Linear_cell_complex_traits<d>, Linear_cell_complex_items_with_label<d> > type;
};

template <class Point, class FT>
struct Point_creator {
public:
  void push_back(FT ft);
  void pop_back();
  Point create_point();
};

template <class Kernel>
struct Point_creator<CGAL::Point_2<Kernel>, typename Kernel::FT> {
private:
  std::vector<typename Kernel::FT> coordinates;
public:
  void push_back(typename Kernel::FT ft) {
    coordinates.push_back(ft);
  }
  
  void pop_back() {
    coordinates.pop_back();
  }
  
  CGAL::Point_2<Kernel> create_point() {
    return CGAL::Point_2<Kernel>(coordinates[0], coordinates[1]);
  }
};

template <class Kernel>
struct Point_creator<CGAL::Point_3<Kernel>, typename Kernel::FT> {
private:
  std::vector<typename Kernel::FT> coordinates;
public:
  void push_back(typename Kernel::FT ft) {
    coordinates.push_back(ft);
  }
  
  void pop_back() {
    coordinates.pop_back();
  }
  
  CGAL::Point_3<Kernel> create_point() {
    return CGAL::Point_3<Kernel>(coordinates[0], coordinates[1], coordinates[2]);
  }
};

template <class Kernel>
struct Point_creator<CGAL::Point_d<Kernel>, typename Kernel::FT> {
private:
  std::vector<typename Kernel::FT> coordinates;
public:
  void push_back(typename Kernel::FT ft) {
    coordinates.push_back(ft);
  }
  
  void pop_back() {
    coordinates.pop_back();
  }
  
  CGAL::Point_d<Kernel> create_point() {
    //std::cout << "Creating Point_d(" << coordinates.size() << ")" << std::endl;
    return CGAL::Point_d<Kernel>(coordinates.size(), coordinates.begin(), coordinates.end());
  }
};

template <typename T>
T string_to_any(const std::string &text, T value) {
  std::stringstream ss;
  for (std::string::const_iterator i = text.begin(); i != text.end(); ++i)
    if (isdigit(*i) || *i=='e' || *i=='-' || *i=='+' || *i=='.') ss << *i;
  T result;
  return ss >> result ? result : value;
}

template <class LCC>
std::ostream &operator<<(std::ostream &output, const CGAL::Dart<LCC::dimension, LCC> &dart) {
  output << "Dart<" << &dart << "> {" << std::endl;
  output << "\tvertex: " << dart.template attribute<0>()->point() << std::endl;
  for (unsigned int current_beta = 0; current_beta < LCC::dimension+1; ++current_beta) {
    output << "\tbeta[" << current_beta << "]: " << &*dart.beta(current_beta) << std::endl;
  } output << "}";
  return output;
}

template <class LCC, unsigned int dimension = LCC::dimension>
struct Cell_complex_attribute_printer;

template <class LCC, unsigned int dimension>
struct Cell_complex_attribute_printer {
  static void print(LCC &lcc) {
    Cell_complex_attribute_printer<LCC, dimension-1>::print(lcc);
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
struct Cell_complex_attribute_printer<LCC, 0> {
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
struct Cell_complex_dart_attribute_printer;

template <class LCC, unsigned int dimension>
struct Cell_complex_dart_attribute_printer {
  static void print(typename LCC::Dart_const_handle &dart) {
    Cell_complex_dart_attribute_printer<LCC, dimension-1>::print(dart);
    if (dart->template attribute<dimension>() != NULL) std::cout << "\t\t" << dimension << "-attribute<" << &*dart->template attribute<dimension>() << ">[" << dart->template attribute<dimension>()->info() << "]" << std::endl;
    else std::cout << "\t\t" << dimension << "-attribute<NULL>" << std::endl;
  }
};

template <class LCC>
struct Cell_complex_dart_attribute_printer<LCC, 0> {
  static void print(typename LCC::Dart_const_handle &dart) {
    std::cout << "\t\t0-attribute<" << &*dart->template attribute<0>() << ">[" << dart->template attribute<0>()->info() << "] (" << dart->template attribute<0>()->point() << ")" << std::endl;
  }
};

template <class LCC, unsigned int dimension = LCC::dimension>
struct Cell_complex_dart_beta_printer;

template <class LCC, unsigned int dimension>
struct Cell_complex_dart_beta_printer {
  static void print(typename LCC::Dart_const_handle &dart) {
    Cell_complex_dart_beta_printer<LCC, dimension-1>::print(dart);
    if (dart->beta(dimension) != LCC::null_handle) std::cout << "\t\tbeta[" << dimension << "]: dart<" << &*dart->beta(dimension) << ">[" << dart->beta(dimension)->id << "]" << std::endl;
    else std::cout << "\t\tbeta[" << dimension << "]: null" << std::endl;
  }
};

template <class LCC>
struct Cell_complex_dart_beta_printer<LCC, 0> {
  static void print(typename LCC::Dart_const_handle &dart) {
    if (dart->beta(0) != LCC::null_handle) std::cout << "\t\tbeta[0]: dart<" << &*dart->beta(0) << ">[" << dart->beta(0)->id << "]" << std::endl;
    else std::cout << "\t\tbeta[0]: null" << std::endl;
  }
};

template <class LCC>
struct Cell_complex_printer {
  static void print(LCC &lcc) {
    Cell_complex_attribute_printer<LCC>::print(lcc);
    
    typename LCC::Dart_const_range &darts = lcc.darts();
    std::cout << "Darts (" << darts.size() << "):" << std::endl;
    for (typename LCC::Dart_const_range::const_iterator current_dart = darts.begin(); current_dart != darts.end(); ++current_dart) {
      std::cout << "\tdart<" << &*current_dart << ">[" << current_dart->id << "]" << std::endl;
      Cell_complex_dart_attribute_printer<LCC>::print(current_dart);
      Cell_complex_dart_beta_printer<LCC>::print(current_dart);
    }
  }
};

// The preconditions avoids linking already used darts and enforces orientability of the map
template <class LCC, class Dart_handle, unsigned int dimension>
struct Gmap_dart_linker {
  static void basic_link_alpha(LCC &lcc, Dart_handle &dart1, Dart_handle &dart2) {
    CGAL_precondition(dart1.cmap_dart() != dart2.cmap_dart());
    CGAL_precondition(lcc.is_free(dart1.cmap_dart(), dimension) ||
                      dart1.cmap_dart()->beta(dimension) == dart2.cmap_dart());
    CGAL_precondition(lcc.is_free(dart2.cmap_dart(), dimension) ||
                      dart2.cmap_dart()->beta(dimension) == dart1.cmap_dart());
    CGAL_precondition(dart1.is_destination() != dart2.is_destination());
    lcc.template basic_link_beta<dimension>(dart1.cmap_dart(), dart2.cmap_dart());
  }
};

template <class LCC, class Dart_handle>
struct Gmap_dart_linker<LCC, Dart_handle, 0> {
  static void basic_link_alpha(LCC &lcc, Dart_handle &dart1, Dart_handle &dart2) {
    CGAL_precondition(dart1.is_destination() != dart2.is_destination());
    CGAL_precondition(dart1.cmap_dart() == dart2.cmap_dart());
  }
};

template <class LCC, class Dart_handle>
struct Gmap_dart_linker<LCC, Dart_handle, 1> {
  static void basic_link_alpha(LCC &lcc, Dart_handle &dart1, Dart_handle &dart2) {
    CGAL_precondition(dart1.cmap_dart() != dart2.cmap_dart());
    if (dart1.is_destination()) {
      CGAL_precondition(lcc.is_free(dart1.cmap_dart(), 1) ||
                        dart1.cmap_dart()->beta(1) == dart2.cmap_dart());
      CGAL_precondition(dart1.is_destination() != dart2.is_destination());
      lcc.template basic_link_beta_1(dart1.cmap_dart(), dart2.cmap_dart());
    } else {
      CGAL_precondition(lcc.is_free(dart2.cmap_dart(), 1) ||
                        dart2.cmap_dart()->beta(1) == dart1.cmap_dart());
      CGAL_precondition(dart1.is_destination() != dart2.is_destination());
      lcc.template basic_link_beta_1(dart2.cmap_dart(), dart1.cmap_dart());
    }
  }
};

template <class LCC>
struct Gmap_faker {
public:
  struct Dart_handle {
  private:
    typename LCC::Dart_handle dart;
    bool destination;
    
  public:
    struct Dart_handle_comparator {
      bool operator()(const Dart_handle &dh1, const Dart_handle &dh2) const {
        if (dh1.cmap_dart() < dh2.cmap_dart()) return true;
        if (dh1.cmap_dart() > dh2.cmap_dart()) return false;
        if (dh1.is_destination() == dh2.is_destination()) return false;
        if (!dh1.is_destination()) return true;
        return false;
      }
    };
    
    // TODO: support passing things by reference for efficiency?
//    Dart_handle() {
//      
//    }
//    
//    Dart_handle(Dart_handle &dh) {
//      dart = dh.cmap_dart();
//      destination = dh.is_destination();
//    }
//    
//    Dart_handle &operator=(Dart_handle &dh) {
//      if (this == &dh) return *this;
//      dart = dh.cmap_dart();
//      destination = dh.is_destination();
//      return *this;
//    }
    
    Dart_handle(typename LCC::Dart_handle da, bool de) {
      CGAL_precondition(da != LCC::null_handle);
      dart = da;
      destination = de;
    }
    
    typename LCC::Dart_handle cmap_dart() const {
//      std::cout << "Dart<" << &*dart << ">" << std::endl;
      return dart;
    }
    
    bool is_destination() const {
      return destination;
    }
  };
  
  static Dart_handle dart_from(typename LCC::Dart_handle dart) {
    Dart_handle new_dart(dart, false);
    return new_dart;
  }
  
  static Dart_handle dart_to(typename LCC::Dart_handle dart) {
    Dart_handle new_dart(dart, true);
    return new_dart;
  }
  
  static bool is_free(LCC &lcc, Dart_handle dart, unsigned int i) {
    if (i == 0) return false; // Warning! This is actually unknown
    if (i == 1 && !dart.is_destination()) return lcc.is_free(dart.cmap_dart(), 0);
    return lcc.is_free(dart.cmap_dart(), i);
  }
  
  static Dart_handle alpha(Dart_handle dart, unsigned int i) {
    if (i == 0) {
      Dart_handle new_dart(dart.cmap_dart(), !dart.is_destination());
      return new_dart;
    } if (i == 1 && !dart.is_destination()) {
      Dart_handle new_dart(dart.cmap_dart()->beta(0), !dart.is_destination());
      return new_dart;
    } Dart_handle new_dart(dart.cmap_dart()->beta(i), !dart.is_destination());
    return new_dart;
  }
  
  template <unsigned int i>
  static void basic_link_alpha(LCC &lcc, Dart_handle &dart1, Dart_handle &dart2) {
    Gmap_dart_linker<LCC, Dart_handle, i>::basic_link_alpha(lcc, dart1, dart2);
  }
  
  template <unsigned int i>
  static typename LCC::template Attribute_const_handle<i>::type attribute(Dart_handle dart) {
    if (i == 0 && dart.is_destination()) return dart.cmap_dart()->beta(1)->template attribute<i>();
    return dart.cmap_dart()->template attribute<i>();
  }
  
  template <unsigned int i>
  static void set_attribute(LCC &lcc, Dart_handle &dart, typename LCC::template Attribute_handle<i>::type &attribute) {
    if (i >= 1 || !dart.is_destination()) {
      lcc.template set_attribute<i>(dart.cmap_dart(), attribute);
    }
    
    else if (i == 0) {
      if (!lcc.is_free(dart.cmap_dart(), 1)) {
        lcc.template set_attribute<i>(dart.cmap_dart()->beta(1), attribute);
      }
    }
  }
  
  // NOTE: Hard-coded for the extrusion case, cannot be generalised easily!
  // Check: If this does not work, the number of darts in a layer will be wrong!
  template <unsigned int i>
  static Dart_handle create_dart_linked_by_alpha(LCC &lcc, Dart_handle dart, bool create = true) {
//    std::cout << "create_dart_linked_by_alpha<" << i << ">" << std::endl;
//    std::cout << i;
    
    // For involutions it's the same behaviour
    if (i > 1) {
      
      // If it has not already been created
      if (is_free(lcc, alpha(dart, 0), i)) {
        std::cout << "Dd";
        Dart_handle new_dart = Dart_handle(lcc.create_dart(), !dart.is_destination());
        lcc.template basic_link_beta<i>(dart.cmap_dart(), new_dart.cmap_dart());
        CGAL_postcondition(dart.is_destination() != new_dart.is_destination());
        return new_dart;
      }
      
      // It it has already been created
      else {
        Dart_handle old_dart = Dart_handle(alpha(alpha(dart, 0), i).cmap_dart(), !alpha(alpha(dart, 0), i).is_destination());
        CGAL_postcondition(dart.is_destination() != old_dart.is_destination());
        return old_dart;
      }
    }
    
    // For an oriented edge
    else if (i == 1) {
      
      if (create) {
        std::cout << "D";
        Dart_handle new_dart = Dart_handle(lcc.create_dart(), !dart.is_destination());
        if (dart.is_destination()) {
          lcc.template basic_link_beta_1(dart.cmap_dart(), new_dart.cmap_dart());
        } else {
          lcc.template basic_link_beta_1(new_dart.cmap_dart(), dart.cmap_dart());
        } CGAL_postcondition(dart.is_destination() != new_dart.is_destination());
        return new_dart;
      }
      
      // It it has already been created
      else {
        std::cout << "d";
        Dart_handle old_dart = alpha(dart, 0);
        while (!is_free(lcc, old_dart, 1)) {
          old_dart = alpha(alpha(old_dart, 1), 0);
        }

        if (dart.is_destination()) {
          lcc.template basic_link_beta_1(dart.cmap_dart(), old_dart.cmap_dart());
        } else {
          lcc.template basic_link_beta_1(old_dart.cmap_dart(), dart.cmap_dart());
        } CGAL_postcondition(dart.is_destination() != old_dart.is_destination());
        return old_dart;
      }
    }
    
    // For alpha_0
    else {
      Dart_handle old_dart = Dart_handle(dart.cmap_dart(), !dart.is_destination());
      std::cout << "d";
      CGAL_postcondition(dart.is_destination() != old_dart.is_destination());
      return old_dart;
    }
  }
};

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
