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

#include "Linear_cell_complex_reader_writer.h"
#include "Linear_cell_complex_extrusion_ranges.h"

#ifndef LCC_EXTRUDER_H
#define LCC_EXTRUDER_H

template <class LCCE, unsigned int dimension>
struct Extruded_embeddings_of_dimension {
public:
  typedef typename LCCE::Lower_dimensional_cell_complex::template Attribute_const_handle<dimension>::type Unextruded_attribute_const_handle;
  typedef typename LCCE::Higher_dimensional_cell_complex::template Attribute_handle<dimension>::type Extruded_lower_dimensional_attribute_handle;
  typedef typename LCCE::Higher_dimensional_cell_complex::template Attribute_handle<dimension+1>::type Extruded_higher_dimensional_attribute_handle;
  
  struct Extruded_embedding {
    Extruded_lower_dimensional_attribute_handle base, top;
    Extruded_higher_dimensional_attribute_handle ex;
  };
  
  typedef std::map<Unextruded_attribute_const_handle, Extruded_embedding> Embeddings;
  
  Embeddings embeddings;
};

template <class LCCE, unsigned int dimension = LCCE::Lower_dimensional_cell_complex::dimension, class Result = CGAL::cpp11::tuple<> >
struct Extruded_embeddings_of_dimension_and_lower;

template <class LCCE, class ... Result>
struct Extruded_embeddings_of_dimension_and_lower<LCCE, 0, CGAL::cpp11::tuple<Result ...> > {
  typedef CGAL::cpp11::tuple<Extruded_embeddings_of_dimension<LCCE, 0>, Result ...> type;
};

template <class LCCE, unsigned int dimension, class ... Result>
struct Extruded_embeddings_of_dimension_and_lower<LCCE, dimension, CGAL::cpp11::tuple<Result ...> > {
  typedef typename Extruded_embeddings_of_dimension_and_lower<LCCE, dimension-1, CGAL::cpp11::tuple<Extruded_embeddings_of_dimension<LCCE, dimension>, Result ...> >::type type;
};

template <class LCCE, unsigned int dimension = LCCE::Lower_dimensional_cell_complex::dimension>
struct Extrude_embeddings_of_dimension_or_lower;

template <class LCCE, unsigned int dimension>
struct Extrude_embeddings_of_dimension_or_lower {
  static void extrude(typename LCCE::Extruded_embeddings &ee, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Higher_dimensional_cell_complex::FT min, typename LCCE::Higher_dimensional_cell_complex::FT max) {
    Extrude_embeddings_of_dimension_or_lower<LCCE, dimension-1>::extrude(ee, ldcc, hdcc, min, max);
    typename LCCE::Lower_dimensional_cell_complex::template Attribute_const_range<dimension>::type &ld_cells = ldcc.template attributes<dimension>();
    for (typename LCCE::Lower_dimensional_cell_complex::template Attribute_const_range<dimension>::type::const_iterator current_cell = ld_cells.begin();
         current_cell != ld_cells.end();
         ++current_cell) {
      typename Extruded_embeddings_of_dimension<LCCE, dimension>::Extruded_embedding &e = std::get<dimension>(ee).embeddings[current_cell];
      e.base = hdcc.template create_attribute<dimension>(current_cell->info());
      e.ex = hdcc.template create_attribute<dimension+1>(current_cell->info()+100);
      e.top = hdcc.template create_attribute<dimension>(current_cell->info()+10);
    }
  }
};

template <class LCCE>
struct Extrude_embeddings_of_dimension_or_lower<LCCE, 0> {
  static void extrude(typename LCCE::Extruded_embeddings &ee, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Higher_dimensional_cell_complex::FT min, typename LCCE::Higher_dimensional_cell_complex::FT max) {
    typename LCCE::Lower_dimensional_cell_complex::Vertex_attribute_const_range &ld_vertices = ldcc.vertex_attributes();
    for (typename LCCE::Lower_dimensional_cell_complex::Vertex_attribute_const_range::const_iterator current_vertex = ld_vertices.begin();
         current_vertex != ld_vertices.end();
         ++current_vertex) {
      typename Extruded_embeddings_of_dimension<LCCE, 0>::Extruded_embedding &e = std::get<0>(ee).embeddings[current_vertex];
      Point_creator<typename LCCE::Higher_dimensional_cell_complex::Point, typename LCCE::Higher_dimensional_cell_complex::FT> creator;
      for (typename LCCE::Lower_dimensional_cell_complex::Point::Cartesian_const_iterator current_coordinate = current_vertex->point().cartesian_begin();
           current_coordinate != current_vertex->point().cartesian_end();
           ++current_coordinate) {
        typename LCCE::Lower_dimensional_cell_complex::FT coordinate = *current_coordinate;
        creator.push_back(coordinate);
      } creator.push_back(min);
      e.base = hdcc.template create_attribute<0>(creator.create_point(), current_vertex->info());
      e.ex = hdcc.template create_attribute<1>(current_vertex->info()+100);
      creator.pop_back();
      creator.push_back(max);
      e.top = hdcc.template create_attribute<0>(creator.create_point(), current_vertex->info()+10);
    }
  }
};

template <class LCCE, unsigned int from_dimension, unsigned int to_dimension>
struct Link_betas_of_dart_according_to_pattern {
  static void link(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle current_dart) {
    if (!ldcc.is_free(current_dart, from_dimension)) hdcc.template basic_link_beta<from_dimension>(layers[current_dart].back(), layers[current_dart->beta(from_dimension)].back());
    Link_betas_of_dart_according_to_pattern<LCCE, from_dimension+1, to_dimension>::link(layers, ldcc, hdcc, current_dart);
  }
};

// Do nothing
template <class LCCE, unsigned int to_dimension>
struct Link_betas_of_dart_according_to_pattern<LCCE, to_dimension, to_dimension> {
  static void link(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle current_dart) {
    
  }
};

template <class LCCE, unsigned int from_dimension, unsigned int to_dimension>
struct Link_betas_of_dart_according_to_lower_pattern {
  static void link(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle current_dart) {
    if (!ldcc.is_free(current_dart, from_dimension-1)) hdcc.template basic_link_beta<from_dimension>(layers[current_dart].back(), layers[current_dart->beta(from_dimension-1)].back());
    Link_betas_of_dart_according_to_lower_pattern<LCCE, from_dimension+1, to_dimension>::link(layers, ldcc, hdcc, current_dart);
  }
};

template <class LCCE, unsigned int to_dimension>
struct Link_betas_of_dart_according_to_lower_pattern<LCCE, to_dimension, to_dimension> {
  static void link(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle current_dart) {
    // Do nothing
  }
};

template <class LCCE, unsigned int from_dimension, unsigned int to_dimension, unsigned int extruding_dimension>
struct Set_extruded_embeddings_base {
  static void set(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle current_dart, typename LCCE::Extruded_embeddings &ee) {
    if (from_dimension <= extruding_dimension) {
      hdcc.template set_attribute<from_dimension>(layers[current_dart].back(), std::get<from_dimension>(ee).embeddings[current_dart->template attribute<from_dimension>()].base);
      std::get<from_dimension>(ee).embeddings[current_dart->template attribute<from_dimension>()].base->set_dart(layers[current_dart].back());
    } else {
      hdcc.template set_attribute<from_dimension>(layers[current_dart].back(), std::get<from_dimension-1>(ee).embeddings[current_dart->template attribute<from_dimension-1>()].ex);
      std::get<from_dimension-1>(ee).embeddings[current_dart->template attribute<from_dimension-1>()].ex->set_dart(layers[current_dart].back());
    } Set_extruded_embeddings_base<LCCE, from_dimension+1, to_dimension, extruding_dimension>::set(layers, hdcc, current_dart, ee);
  }
};

template <class LCCE, unsigned int from_dimension, unsigned int extruding_dimension>
struct Set_extruded_embeddings_base<LCCE, from_dimension, from_dimension, extruding_dimension> {
  static void set(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle current_dart, typename LCCE::Extruded_embeddings &ee) {
    hdcc.template set_attribute<from_dimension>(layers[current_dart].back(), std::get<from_dimension-1>(ee).embeddings[current_dart->template attribute<from_dimension-1>()].ex);
    std::get<from_dimension-1>(ee).embeddings[current_dart->template attribute<from_dimension-1>()].ex->set_dart(layers[current_dart].back());
  }
};

template <class LCCE, unsigned int from_dimension, unsigned int to_dimension, unsigned int extruding_dimension>
struct Set_extruded_embeddings_top {
  static void set(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle current_dart, typename LCCE::Extruded_embeddings &ee) {
    if (from_dimension <= extruding_dimension) {
      hdcc.template set_attribute<from_dimension>(layers[current_dart].back(), std::get<from_dimension>(ee).embeddings[current_dart->template attribute<from_dimension>()].top);
      std::get<from_dimension>(ee).embeddings[current_dart->template attribute<from_dimension>()].top->set_dart(layers[current_dart].back());
    } else {
      hdcc.template set_attribute<from_dimension>(layers[current_dart].back(), std::get<from_dimension-1>(ee).embeddings[current_dart->template attribute<from_dimension-1>()].ex);
      std::get<from_dimension-1>(ee).embeddings[current_dart->template attribute<from_dimension-1>()].ex->set_dart(layers[current_dart].back());
    } Set_extruded_embeddings_top<LCCE, from_dimension+1, to_dimension, extruding_dimension>::set(layers, hdcc, current_dart, ee);
  }
};

template <class LCCE, unsigned int from_dimension, unsigned int extruding_dimension>
struct Set_extruded_embeddings_top<LCCE, from_dimension, from_dimension, extruding_dimension> {
  static void set(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle current_dart, typename LCCE::Extruded_embeddings &ee) {
    hdcc.template set_attribute<from_dimension>(layers[current_dart].back(), std::get<from_dimension-1>(ee).embeddings[current_dart->template attribute<from_dimension-1>()].ex);
    std::get<from_dimension-1>(ee).embeddings[current_dart->template attribute<from_dimension-1>()].ex->set_dart(layers[current_dart].back());
  }
};

template <class LCCE, unsigned int dimension, unsigned int dimension_to>
struct Extrude_darts_top {
  static void create_darts_layer(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extruded_embeddings &ee) {
    std::cout << "Extruding top at dimension " << dimension << "..." << std::endl;
    
    // Make a new darts layer
    typename LCCE::Lower_dimensional_cell_complex::Dart_const_range &ld_darts = ldcc.darts();
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_const_range::const_iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      layers[current_dart].push_back(hdcc.create_dart(hdcc.darts().size()));
    }
    
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_const_range::const_iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Link darts
      if (dimension % 2 == 1) hdcc.basic_link_beta_1(layers[current_dart].back(), layers[current_dart->beta(0)].back());
      else hdcc.basic_link_beta_1(layers[current_dart].back(), layers[current_dart->beta(1)].back());
      Link_betas_of_dart_according_to_pattern<LCCE, 2, dimension>::link(layers, ldcc, hdcc, current_dart);
      hdcc.template basic_link_beta<dimension>(layers[current_dart].back(), layers[current_dart][layers[current_dart].size()-2]);
      Link_betas_of_dart_according_to_lower_pattern<LCCE, dimension+2, LCCE::Lower_dimensional_cell_complex::dimension+2>::link(layers, ldcc, hdcc, current_dart);
      
      // Set correct embeddings
      if (dimension % 2 == 1) {
        hdcc.template set_attribute<0>(layers[current_dart].back(), std::get<0>(ee).embeddings[current_dart->beta(1)->template attribute<0>()].top);
        std::get<0>(ee).embeddings[current_dart->beta(1)->template attribute<0>()].top->set_dart(layers[current_dart].back());
      } else {
        hdcc.template set_attribute<0>(layers[current_dart].back(), std::get<0>(ee).embeddings[current_dart->template attribute<0>()].top);
        std::get<0>(ee).embeddings[current_dart->template attribute<0>()].top->set_dart(layers[current_dart].back());
      } Set_extruded_embeddings_top<LCCE, 1, LCCE::Higher_dimensional_cell_complex::dimension, dimension>::set(layers, hdcc, current_dart, ee);
    }
    
    Extrude_darts_top<LCCE, dimension+1, dimension_to>::create_darts_layer(layers, ldcc, hdcc, ee);
  }
};

template <class LCCE, unsigned int dimension_to>
struct Extrude_darts_top<LCCE, 0, dimension_to> {
  static void create_darts_layer(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extruded_embeddings &ee) {
    std::cout << "Extruding top at dimension 0..." << std::endl;
    
    // Make a new darts layer
    typename LCCE::Lower_dimensional_cell_complex::Dart_const_range &ld_darts = ldcc.darts();
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_const_range::const_iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      layers[current_dart].push_back(hdcc.create_dart(hdcc.darts().size()));
    }
    
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_const_range::const_iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Link darts
      hdcc.basic_link_beta_1(layers[current_dart].back(), layers[current_dart][layers[current_dart].size()-3]);
      hdcc.template basic_link_beta<2>(layers[current_dart].back(), layers[current_dart->beta(0)][layers[current_dart].size()-2]);
      Link_betas_of_dart_according_to_lower_pattern<LCCE, 3, LCCE::Lower_dimensional_cell_complex::dimension+2>::link(layers, ldcc, hdcc, current_dart);
      
      // Set correct embeddings
      hdcc.template set_attribute<0>(layers[current_dart].back(), std::get<0>(ee).embeddings[current_dart->template attribute<0>()].top);
      std::get<0>(ee).embeddings[current_dart->template attribute<0>()].top->set_dart(layers[current_dart].back());
      Set_extruded_embeddings_top<LCCE, 1, LCCE::Higher_dimensional_cell_complex::dimension, 0>::set(layers, hdcc, current_dart, ee);
    }
    
    Extrude_darts_top<LCCE, 1, dimension_to>::create_darts_layer(layers, ldcc, hdcc, ee);
  }
};

template <class LCCE, unsigned int dimension_to>
struct Extrude_darts_top<LCCE, 1, dimension_to> {
  static void create_darts_layer(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extruded_embeddings &ee) {
    std::cout << "Extruding top at dimension 1..." << std::endl;
    
    // Make a new darts layer
    typename LCCE::Lower_dimensional_cell_complex::Dart_const_range &ld_darts = ldcc.darts();
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_const_range::const_iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      layers[current_dart].push_back(hdcc.create_dart(hdcc.darts().size()));
    }
    
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_const_range::const_iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Link darts
      hdcc.basic_link_beta_1(layers[current_dart][layers[current_dart].size()-3], layers[current_dart].back());
      hdcc.basic_link_beta_1(layers[current_dart].back(), layers[current_dart][layers[current_dart].size()-2]);
      Link_betas_of_dart_according_to_lower_pattern<LCCE, 3, LCCE::Lower_dimensional_cell_complex::dimension+2>::link(layers, ldcc, hdcc, current_dart);
      
      // Set correct embeddings
      hdcc.template set_attribute<0>(layers[current_dart].back(), std::get<0>(ee).embeddings[current_dart->beta(1)->template attribute<0>()].top);
      std::get<0>(ee).embeddings[current_dart->beta(1)->template attribute<0>()].top->set_dart(layers[current_dart].back());
      Set_extruded_embeddings_top<LCCE, 1, LCCE::Higher_dimensional_cell_complex::dimension, 1>::set(layers, hdcc, current_dart, ee);
    }
    
    Extrude_darts_top<LCCE, 2, dimension_to>::create_darts_layer(layers, ldcc, hdcc, ee);
  }
};

template <class LCCE, unsigned int dimension_to>
struct Extrude_darts_top<LCCE, dimension_to, dimension_to> {
  static void create_darts_layer(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extruded_embeddings &ee) {
    std::cout << "Extruding top at dimension " << dimension_to << "..." << std::endl;
    
    // Make a new darts layer
    typename LCCE::Lower_dimensional_cell_complex::Dart_const_range &ld_darts = ldcc.darts();
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_const_range::const_iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      layers[current_dart].push_back(hdcc.create_dart(hdcc.darts().size()));
    }
    
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_const_range::const_iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Link darts
      if (dimension_to % 2 == 1) hdcc.basic_link_beta_1(layers[current_dart].back(), layers[current_dart->beta(0)].back());
      else hdcc.basic_link_beta_1(layers[current_dart].back(), layers[current_dart->beta(1)].back());
      Link_betas_of_dart_according_to_pattern<LCCE, 2, dimension_to>::link(layers, ldcc, hdcc, current_dart);
      hdcc.template basic_link_beta<dimension_to>(layers[current_dart].back(), layers[current_dart][layers[current_dart].size()-2]);
      Link_betas_of_dart_according_to_lower_pattern<LCCE, dimension_to+2, LCCE::Lower_dimensional_cell_complex::dimension+2>::link(layers, ldcc, hdcc, current_dart);
      
      // Set correct embeddings
      if (dimension_to % 2 == 1) {
        hdcc.template set_attribute<0>(layers[current_dart].back(), std::get<0>(ee).embeddings[current_dart->beta(1)->template attribute<0>()].top);
        std::get<0>(ee).embeddings[current_dart->beta(1)->template attribute<0>()].top->set_dart(layers[current_dart].back());
      } else {
        hdcc.template set_attribute<0>(layers[current_dart].back(), std::get<0>(ee).embeddings[current_dart->template attribute<0>()].top);
        std::get<0>(ee).embeddings[current_dart->template attribute<0>()].top->set_dart(layers[current_dart].back());
      } Set_extruded_embeddings_top<LCCE, 1, LCCE::Higher_dimensional_cell_complex::dimension, dimension_to>::set(layers, hdcc, current_dart, ee);
    }
  }
};

template <class LCCE, unsigned int dimension>
struct Extrude_darts_base {
  static void create_darts_layer(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extruded_embeddings &ee) {
    std::cout << "Extruding base at dimension " << dimension << "..." << std::endl;
    
    // Make a new darts layer
    typename LCCE::Lower_dimensional_cell_complex::Dart_const_range &ld_darts = ldcc.darts();
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_const_range::const_iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      layers[current_dart].push_back(hdcc.create_dart(hdcc.darts().size()));
    }
    
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_const_range::const_iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Link darts
      if (dimension % 2 == 0) hdcc.basic_link_beta_1(layers[current_dart].back(), layers[current_dart->beta(0)].back());
      else hdcc.basic_link_beta_1(layers[current_dart].back(), layers[current_dart->beta(1)].back());
      Link_betas_of_dart_according_to_pattern<LCCE, 2, dimension>::link(layers, ldcc, hdcc, current_dart);
      if (dimension < LCCE::Lower_dimensional_cell_complex::dimension) hdcc.template basic_link_beta<dimension+1>(layers[current_dart].back(), layers[current_dart][layers[current_dart].size()-2]);
      Link_betas_of_dart_according_to_lower_pattern<LCCE, dimension+2, LCCE::Lower_dimensional_cell_complex::dimension+2>::link(layers, ldcc, hdcc, current_dart);
      
      // Set correct embeddings
      if (dimension % 2 == 0) {
        hdcc.template set_attribute<0>(layers[current_dart].back(), std::get<0>(ee).embeddings[current_dart->beta(1)->template attribute<0>()].base);
        std::get<0>(ee).embeddings[current_dart->beta(1)->template attribute<0>()].base->set_dart(layers[current_dart].back());
      } else {
        hdcc.template set_attribute<0>(layers[current_dart].back(), std::get<0>(ee).embeddings[current_dart->template attribute<0>()].base);
        std::get<0>(ee).embeddings[current_dart->template attribute<0>()].base->set_dart(layers[current_dart].back());
      } Set_extruded_embeddings_base<LCCE, 1, LCCE::Higher_dimensional_cell_complex::dimension, dimension>::set(layers, hdcc, current_dart, ee);
    }
    
    Extrude_darts_base<LCCE, dimension-1>::create_darts_layer(layers, ldcc, hdcc, ee);
  }
};

template <class LCCE>
struct Extrude_darts_base<LCCE, 1> {
  static void create_darts_layer(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extruded_embeddings &ee) {
    std::cout << "Extruding base at dimension 1..." << std::endl;
    
    // Make a new darts layer
    typename LCCE::Lower_dimensional_cell_complex::Dart_const_range &ld_darts = ldcc.darts();
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_const_range::const_iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      layers[current_dart].push_back(hdcc.create_dart(hdcc.darts().size()));
    }
    
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_const_range::const_iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Link darts
      hdcc.template basic_link_beta<2>(layers[current_dart].back(), layers[current_dart][layers[current_dart].size()-2]);
      Link_betas_of_dart_according_to_lower_pattern<LCCE, 3, LCCE::Lower_dimensional_cell_complex::dimension+2>::link(layers, ldcc, hdcc, current_dart);
      
      // Set correct embeddings
      hdcc.template set_attribute<0>(layers[current_dart].back(), std::get<0>(ee).embeddings[current_dart->template attribute<0>()].base);
      std::get<0>(ee).embeddings[current_dart->template attribute<0>()].base->set_dart(layers[current_dart].back());
      Set_extruded_embeddings_base<LCCE, 1, LCCE::Higher_dimensional_cell_complex::dimension, 1>::set(layers, hdcc, current_dart, ee);
    }
    
    Extrude_darts_base<LCCE, 0>::create_darts_layer(layers, ldcc, hdcc, ee);
  }
};

template <class LCCE>
struct Extrude_darts_base<LCCE, 0> {
  static void create_darts_layer(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extruded_embeddings &ee) {
    std::cout << "Extruding base at dimension 0..." << std::endl;
    
    // Make a new darts layer
    typename LCCE::Lower_dimensional_cell_complex::Dart_const_range &ld_darts = ldcc.darts();
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_const_range::const_iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      layers[current_dart].push_back(hdcc.create_dart(hdcc.darts().size()));
    }
    
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_const_range::const_iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Link darts
      hdcc.basic_link_beta_1(layers[current_dart][layers[current_dart].size()-2], layers[current_dart].back());
      Link_betas_of_dart_according_to_lower_pattern<LCCE, 3, LCCE::Lower_dimensional_cell_complex::dimension+2>::link(layers, ldcc, hdcc, current_dart);
      
      // Set correct embeddings
      hdcc.template set_attribute<0>(layers[current_dart].back(), std::get<0>(ee).embeddings[current_dart->beta(1)->template attribute<0>()].base);
      std::get<0>(ee).embeddings[current_dart->beta(1)->template attribute<0>()].base->set_dart(layers[current_dart].back());
      hdcc.template set_attribute<1>(layers[current_dart].back(), std::get<0>(ee).embeddings[current_dart->beta(1)->template attribute<0>()].ex);
      std::get<0>(ee).embeddings[current_dart->beta(1)->template attribute<0>()].ex->set_dart(layers[current_dart].back());
      Set_extruded_embeddings_base<LCCE, 2, LCCE::Higher_dimensional_cell_complex::dimension, 0>::set(layers, hdcc, current_dart, ee);
    }
    
    Extrude_darts_top<LCCE, 0, LCCE::Lower_dimensional_cell_complex::dimension>::create_darts_layer(layers, ldcc, hdcc, ee);
  }
};

template <class LCCE, unsigned int dimension = LCCE::Lower_dimensional_cell_complex::dimension>
struct Extrude_darts {
  static void extrude(std::map<typename LCCE::Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename LCCE::Higher_dimensional_cell_complex::Dart_handle> > &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extruded_embeddings &ee) {
    Extrude_darts_base<LCCE, dimension>::create_darts_layer(layers, ldcc, hdcc, ee);
  }
};



template <unsigned int unextruded_dimension>
class Linear_cell_complex_extruder {
public:
  typedef typename Linear_cell_complex_with_ids<unextruded_dimension>::type Lower_dimensional_cell_complex;
  typedef typename Linear_cell_complex_with_ids<unextruded_dimension+1>::type Higher_dimensional_cell_complex;
  typedef Linear_cell_complex_extruder<unextruded_dimension> Self;
  typedef typename Lower_dimensional_cell_complex::FT FT;
  
  typedef typename Extruded_embeddings_of_dimension_and_lower<Self>::type Extruded_embeddings;
  
  Higher_dimensional_cell_complex extrude(Lower_dimensional_cell_complex &ldcc, FT min, FT max) {
    CGAL_precondition(Lower_dimensional_cell_complex::dimension > 1);
    Higher_dimensional_cell_complex hdcc;
    
    // Extrude the embeddings
    Extruded_embeddings extruded_embeddings;
    Extrude_embeddings_of_dimension_or_lower<Self>::extrude(extruded_embeddings, ldcc, hdcc, min, max);
    
    // Extrude the darts
    std::map<typename Lower_dimensional_cell_complex::Dart_const_handle, std::vector<typename Higher_dimensional_cell_complex::Dart_handle> > layers;
    typename Lower_dimensional_cell_complex::Dart_const_range &ld_darts = ldcc.darts();
    for (typename Lower_dimensional_cell_complex::Dart_const_range::const_iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      layers[current_dart] = std::vector<typename Higher_dimensional_cell_complex::Dart_handle>();
    } Extrude_darts<Self>::extrude(layers, ldcc, hdcc, extruded_embeddings);
    
    return hdcc;
  }
};

#endif
