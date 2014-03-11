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

#include "Linear_cell_complex_extruder.h"

#ifndef LCC_EXTRUDER_WITH_RANGE_H
#define LCC_EXTRUDER_WITH_RANGE_H

template <class ERTPD, unsigned int dimension = ERTPD::LCC::dimension>
struct Pass_extrusion_ranges_to_boundary;

// TODO: Add a special functions for 1-cells?
// Maybe not needed since the dual (any beta_2 to beta_n) of a dart will add the opposite vertex of an edge
template <class ERTPD, unsigned int dimension>
struct Pass_extrusion_ranges_to_boundary {
  static void pass(typename ERTPD::LCC &lcc, ERTPD &ertpd) {
    
    // Check every dimension-cell
    for (typename Extrusion_ranges_map_of_dimension<typename ERTPD::LCC, dimension>::type::const_iterator current_cell_range = std::get<dimension>(ertpd.ranges).ranges_map.begin();
         current_cell_range != std::get<dimension>(ertpd.ranges).ranges_map.end();
         ++current_cell_range) {
      
      // Find all the (dimension-1)-cells on its boundary
      std::set<typename ERTPD::LCC::template Attribute_const_handle<dimension-1>::type> boundary_cells;
      typename ERTPD::LCC::Dart_const_handle dart_in_cell = current_cell_range->first->dart();
      typename ERTPD::LCC::template Dart_of_cell_const_range<dimension> cell_darts_range = lcc.template darts_of_cell<dimension>(dart_in_cell);
      for (typename ERTPD::LCC::template Dart_of_cell_const_range<dimension>::const_iterator current_dart = cell_darts_range.begin();
           current_dart != cell_darts_range.end();
           ++current_dart) {
        boundary_cells.insert(current_dart->template attribute<dimension-1>());
      } //std::cout << dimension << "-cell<" << &*current_cell_range->first << "> has " << boundary_cells.size() << " " << dimension-1 << "-cells in its boundary" << std::endl;
      
      // Pass the ranges of to its boundary cells
      for (typename std::set<typename ERTPD::LCC::template Attribute_const_handle<dimension-1>::type>::const_iterator current_boundary_cell = boundary_cells.begin();
           current_boundary_cell != boundary_cells.end();
           ++current_boundary_cell) {
        std::get<dimension-1>(ertpd.ranges).ranges_map[*current_boundary_cell].add_ranges(current_cell_range->second);
      }
    }
    
    Pass_extrusion_ranges_to_boundary<ERTPD, dimension-1>::pass(lcc, ertpd);
  }
};

template <class ERTPD>
struct Pass_extrusion_ranges_to_boundary<ERTPD, 0> {
  static void pass(typename ERTPD::LCC &lcc, ERTPD &ertpd) {
    // Do nothing, these have no boundary
  }
};

template <class LCCE, unsigned int dimension>
struct Extruded_embeddings_with_range_of_dimension {
public:
  typedef typename LCCE::Lower_dimensional_cell_complex::template Attribute_const_handle<dimension>::type Unextruded_attribute_const_handle;
  typedef typename LCCE::Higher_dimensional_cell_complex::template Attribute_handle<dimension>::type Extruded_lower_dimensional_attribute_handle;
  typedef typename LCCE::Higher_dimensional_cell_complex::template Attribute_handle<dimension+1>::type Extruded_higher_dimensional_attribute_handle;
  
  struct Extruded_embedding {
  public:
    std::map<typename LCCE::Lower_dimensional_cell_complex::FT, Extruded_lower_dimensional_attribute_handle> nonex;
    std::map<typename LCCE::Lower_dimensional_cell_complex::FT, Extruded_higher_dimensional_attribute_handle> ex; // use the lower end of the range
    
    // TODO: More complex checking? Optimised search?
    Extruded_higher_dimensional_attribute_handle &find_ex(typename LCCE::Lower_dimensional_cell_complex::FT &min,
                                                          typename LCCE::Lower_dimensional_cell_complex::FT &max) {
      if (ex.count(min) > 0) return ex[min];
      else {
        for (typename std::map<typename LCCE::Lower_dimensional_cell_complex::FT, Extruded_higher_dimensional_attribute_handle>::const_reverse_iterator current_attribute = ex.rbegin();
             current_attribute != ex.rend();
             ++current_attribute) {
          if (current_attribute->first < min) return ex[current_attribute->first];
        } return ex[ex.begin()->first];
      }
    }
  };
  
  typedef std::map<Unextruded_attribute_const_handle, Extruded_embedding> Embeddings;
  
  Embeddings embeddings;
};

template <class LCCE, unsigned int dimension = LCCE::Lower_dimensional_cell_complex::dimension, class Result = CGAL::cpp11::tuple<> >
struct Extruded_embeddings_with_range_of_dimension_and_lower;

template <class LCCE, class ... Result>
struct Extruded_embeddings_with_range_of_dimension_and_lower<LCCE, 0, CGAL::cpp11::tuple<Result ...> > {
  typedef CGAL::cpp11::tuple<Extruded_embeddings_with_range_of_dimension<LCCE, 0>, Result ...> type;
};

template <class LCCE, unsigned int dimension, class ... Result>
struct Extruded_embeddings_with_range_of_dimension_and_lower<LCCE, dimension, CGAL::cpp11::tuple<Result ...> > {
  typedef typename Extruded_embeddings_with_range_of_dimension_and_lower<LCCE, dimension-1, CGAL::cpp11::tuple<Extruded_embeddings_with_range_of_dimension<LCCE, dimension>, Result ...> >::type type;
};

template <class LCCE, unsigned int dimension = LCCE::Lower_dimensional_cell_complex::dimension>
struct Extrude_embeddings_with_range_of_dimension_or_lower;

template <class LCCE, unsigned int dimension>
struct Extrude_embeddings_with_range_of_dimension_or_lower {
  static void extrude(typename LCCE::Extruded_embeddings_with_range &ee, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extrusion_ranges &er) {
    Extrude_embeddings_with_range_of_dimension_or_lower<LCCE, dimension-1>::extrude(ee, ldcc, hdcc, er);
    typename LCCE::Lower_dimensional_cell_complex::template Attribute_const_range<dimension>::type &ld_cells = ldcc.template attributes<dimension>();
    for (typename LCCE::Lower_dimensional_cell_complex::template Attribute_const_range<dimension>::type::const_iterator current_cell = ld_cells.begin();
         current_cell != ld_cells.end();
         ++current_cell) {
      typename Extruded_embeddings_with_range_of_dimension<LCCE, dimension>::Extruded_embedding &e = std::get<dimension>(ee).embeddings[current_cell];
      Extrusion_ranges<typename LCCE::Lower_dimensional_cell_complex> &r = std::get<dimension>(er.ranges).ranges_map[current_cell];
      for (typename Extrusion_ranges<typename LCCE::Lower_dimensional_cell_complex>::Range_set::const_iterator current_range = r.begin();
           current_range != r.end();
           ++current_range) {
        for (typename Extrusion_range<typename LCCE::Lower_dimensional_cell_complex>::Range::const_iterator current_element = current_range->begin();
             current_element != current_range->end();
             ++current_element) {
          e.nonex[*current_element] = hdcc.template create_attribute<dimension>(current_cell->info());
          typename Extrusion_range<typename LCCE::Lower_dimensional_cell_complex>::Range::const_iterator next_element = current_element;
          ++next_element;
          if (next_element != current_range->end()) e.ex[*current_element] = hdcc.template create_attribute<dimension+1>(current_cell->info());
        }
      }
    }
  }
};

template <class LCCE>
struct Extrude_embeddings_with_range_of_dimension_or_lower<LCCE, 0> {
  static void extrude(typename LCCE::Extruded_embeddings_with_range &ee, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extrusion_ranges &er) {
    typename LCCE::Lower_dimensional_cell_complex::Vertex_attribute_const_range &ld_vertices = ldcc.vertex_attributes();
    for (typename LCCE::Lower_dimensional_cell_complex::Vertex_attribute_const_range::const_iterator current_vertex = ld_vertices.begin();
         current_vertex != ld_vertices.end();
         ++current_vertex) {
      typename Extruded_embeddings_with_range_of_dimension<LCCE, 0>::Extruded_embedding &e = std::get<0>(ee).embeddings[current_vertex];
      Extrusion_ranges<typename LCCE::Lower_dimensional_cell_complex> &r = std::get<0>(er.ranges).ranges_map[current_vertex];
      for (typename Extrusion_ranges<typename LCCE::Lower_dimensional_cell_complex>::Range_set::const_iterator current_range = r.begin();
           current_range != r.end();
           ++current_range) {
        for (typename Extrusion_range<typename LCCE::Lower_dimensional_cell_complex>::Range::const_iterator current_element = current_range->begin();
             current_element != current_range->end();
             ++current_element) {
          Point_creator<typename LCCE::Higher_dimensional_cell_complex::Point, typename LCCE::Higher_dimensional_cell_complex::FT> creator;
          for (typename LCCE::Lower_dimensional_cell_complex::Point::Cartesian_const_iterator current_coordinate = current_vertex->point().cartesian_begin();
               current_coordinate != current_vertex->point().cartesian_end();
               ++current_coordinate) {
            typename LCCE::Lower_dimensional_cell_complex::FT coordinate = *current_coordinate;
            creator.push_back(coordinate);
          } creator.push_back(*current_element);
          e.nonex[*current_element] = hdcc.template create_attribute<0>(creator.create_point(), current_vertex->info());
          typename Extrusion_range<typename LCCE::Lower_dimensional_cell_complex>::Range::const_iterator next_element = current_element;
          ++next_element;
          if (next_element != current_range->end()) e.ex[*current_element] = hdcc.template create_attribute<1>(current_vertex->info());
        }
      }
    }
  }
};

template <class LCCE, unsigned int from_dimension, unsigned int to_dimension, unsigned int extruding_dimension>
struct Set_extruded_embeddings_base_with_range {
  static void set(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::LDGF::Dart_handle current_dart, typename LCCE::Extruded_embeddings_with_range &ee, typename LCCE::FT current_range_min, typename LCCE::FT current_range_max) {
    if (from_dimension <= extruding_dimension) {
      LCCE::HDGF::template set_attribute<from_dimension>(hdcc, layers[current_dart].back(), std::get<from_dimension>(ee).embeddings[LCCE::LDGF::template attribute<from_dimension>(current_dart)].nonex[current_range_min]);
      std::get<from_dimension>(ee).embeddings[LCCE::LDGF::template attribute<from_dimension>(current_dart)].nonex[current_range_min]->set_dart(layers[current_dart].back().cmap_dart());
    } else {
      LCCE::HDGF::template set_attribute<from_dimension>(hdcc, layers[current_dart].back(), std::get<from_dimension-1>(ee).embeddings[LCCE::LDGF::template attribute<from_dimension-1>(current_dart)].find_ex(current_range_min, current_range_max));
      if (from_dimension == 1 && current_dart.is_destination()) {
        if (!hdcc.is_free(layers[current_dart].back().cmap_dart(), 1)) std::get<from_dimension-1>(ee).embeddings[LCCE::LDGF::template attribute<from_dimension-1>(current_dart)].find_ex(current_range_min, current_range_max)->set_dart(layers[current_dart].back().cmap_dart()->beta(1));
      } else std::get<from_dimension-1>(ee).embeddings[LCCE::LDGF::template attribute<from_dimension-1>(current_dart)].find_ex(current_range_min, current_range_max)->set_dart(layers[current_dart].back().cmap_dart());
    } Set_extruded_embeddings_base_with_range<LCCE, from_dimension+1, to_dimension, extruding_dimension>::set(layers, hdcc, current_dart, ee, current_range_min, current_range_max);
  }
};

template <class LCCE, unsigned int to_dimension, unsigned int extruding_dimension>
struct Set_extruded_embeddings_base_with_range<LCCE, 0, to_dimension, extruding_dimension> {
  static void set(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::LDGF::Dart_handle current_dart, typename LCCE::Extruded_embeddings_with_range &ee, typename LCCE::FT current_range_min, typename LCCE::FT current_range_max) {
    LCCE::HDGF::template set_attribute<0>(hdcc, layers[current_dart].back(), std::get<0>(ee).embeddings[LCCE::LDGF::template attribute<0>(current_dart)].nonex[current_range_min]);
    if (!layers[current_dart].back().is_destination()) std::get<0>(ee).embeddings[LCCE::LDGF::template attribute<0>(current_dart)].nonex[current_range_min]->set_dart(layers[current_dart].back().cmap_dart());
    else if (!hdcc.is_free(layers[current_dart].back().cmap_dart(), 1)) std::get<0>(ee).embeddings[LCCE::LDGF::template attribute<0>(current_dart)].nonex[current_range_min]->set_dart(layers[current_dart].back().cmap_dart()->beta(1));
    Set_extruded_embeddings_base_with_range<LCCE, 1, to_dimension, extruding_dimension>::set(layers, hdcc, current_dart, ee, current_range_min, current_range_max);
  }
};

template <class LCCE, unsigned int from_dimension, unsigned int extruding_dimension>
struct Set_extruded_embeddings_base_with_range<LCCE, from_dimension, from_dimension, extruding_dimension> {
  static void set(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::LDGF::Dart_handle current_dart, typename LCCE::Extruded_embeddings_with_range &ee, typename LCCE::FT current_range_min, typename LCCE::FT current_range_max) {
    LCCE::HDGF::template set_attribute<from_dimension>(hdcc, layers[current_dart].back(), std::get<from_dimension-1>(ee).embeddings[LCCE::LDGF::template attribute<from_dimension-1>(current_dart)].find_ex(current_range_min, current_range_max));
    std::get<from_dimension-1>(ee).embeddings[LCCE::LDGF::template attribute<from_dimension-1>(current_dart)].find_ex(current_range_min, current_range_max)->set_dart(layers[current_dart].back().cmap_dart());
  }
};

template <class LCCE, unsigned int from_dimension, unsigned int to_dimension, unsigned int extruding_dimension>
struct Set_extruded_embeddings_top_with_range {
  static void set(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::LDGF::Dart_handle current_dart, typename LCCE::Extruded_embeddings_with_range &ee, typename LCCE::FT current_range_min, typename LCCE::FT current_range_max) {
    if (from_dimension <= extruding_dimension) {
      LCCE::HDGF::template set_attribute<from_dimension>(hdcc, layers[current_dart].back(), std::get<from_dimension>(ee).embeddings[LCCE::LDGF::template attribute<from_dimension>(current_dart)].nonex[current_range_max]);
      std::get<from_dimension>(ee).embeddings[LCCE::LDGF::template attribute<from_dimension>(current_dart)].nonex[current_range_max]->set_dart(layers[current_dart].back().cmap_dart());
    } else {
      LCCE::HDGF::template set_attribute<from_dimension>(hdcc, layers[current_dart].back(), std::get<from_dimension-1>(ee).embeddings[LCCE::LDGF::template attribute<from_dimension-1>(current_dart)].find_ex(current_range_min, current_range_max));
      if (from_dimension == 1 && current_dart.is_destination()) {
        if (!hdcc.is_free(layers[current_dart].back().cmap_dart(), 1)) std::get<from_dimension-1>(ee).embeddings[LCCE::LDGF::template attribute<from_dimension-1>(current_dart)].find_ex(current_range_min, current_range_max)->set_dart(layers[current_dart].back().cmap_dart()->beta(1));
      } else std::get<from_dimension-1>(ee).embeddings[LCCE::LDGF::template attribute<from_dimension-1>(current_dart)].find_ex(current_range_min, current_range_max)->set_dart(layers[current_dart].back().cmap_dart());
    } Set_extruded_embeddings_top_with_range<LCCE, from_dimension+1, to_dimension, extruding_dimension>::set(layers, hdcc, current_dart, ee, current_range_min, current_range_max);
  }
};

template <class LCCE, unsigned int to_dimension, unsigned int extruding_dimension>
struct Set_extruded_embeddings_top_with_range<LCCE, 0, to_dimension, extruding_dimension> {
  static void set(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::LDGF::Dart_handle current_dart, typename LCCE::Extruded_embeddings_with_range &ee, typename LCCE::FT current_range_min, typename LCCE::FT current_range_max) {
    LCCE::HDGF::template set_attribute<0>(hdcc, layers[current_dart].back(), std::get<0>(ee).embeddings[LCCE::LDGF::template attribute<0>(current_dart)].nonex[current_range_max]);
    if (!layers[current_dart].back().is_destination()) std::get<0>(ee).embeddings[LCCE::LDGF::template attribute<0>(current_dart)].nonex[current_range_max]->set_dart(layers[current_dart].back().cmap_dart());
    else if (!hdcc.is_free(layers[current_dart].back().cmap_dart(), 1)) {
      std::get<0>(ee).embeddings[LCCE::LDGF::template attribute<0>(current_dart)].nonex[current_range_max]->set_dart(layers[current_dart].back().cmap_dart()->beta(1));
    } Set_extruded_embeddings_top_with_range<LCCE, 1, to_dimension, extruding_dimension>::set(layers, hdcc, current_dart, ee, current_range_min, current_range_max);
  }
};

template <class LCCE, unsigned int from_dimension, unsigned int extruding_dimension>
struct Set_extruded_embeddings_top_with_range<LCCE, from_dimension, from_dimension, extruding_dimension> {
  static void set(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::LDGF::Dart_handle current_dart, typename LCCE::Extruded_embeddings_with_range &ee, typename LCCE::FT current_range_min, typename LCCE::FT current_range_max) {
    LCCE::HDGF::template set_attribute<from_dimension>(hdcc, layers[current_dart].back(), std::get<from_dimension-1>(ee).embeddings[LCCE::LDGF::template attribute<from_dimension-1>(current_dart)].find_ex(current_range_min, current_range_max));
    std::get<from_dimension-1>(ee).embeddings[LCCE::LDGF::template attribute<from_dimension-1>(current_dart)].find_ex(current_range_min, current_range_max)->set_dart(layers[current_dart].back().cmap_dart());
  }
};

template <class LCCE, unsigned int from_dimension, unsigned int to_dimension>
struct Link_alphas_of_dart_according_to_pattern {
  static void link(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::LDGF::Dart_handle current_dart) {
    if (!LCCE::LDGF::is_free(ldcc, current_dart, from_dimension)) LCCE::HDGF::template basic_link_alpha<from_dimension>(hdcc, layers[current_dart].back(), layers[LCCE::LDGF::alpha(current_dart, from_dimension)].back());
    Link_alphas_of_dart_according_to_pattern<LCCE, from_dimension+1, to_dimension>::link(layers, ldcc, hdcc, current_dart);
  }
};

// Do nothing
template <class LCCE, unsigned int to_dimension>
struct Link_alphas_of_dart_according_to_pattern<LCCE, to_dimension, to_dimension> {
  static void link(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::LDGF::Dart_handle current_dart) {
    
  }
};

template <class LCCE, unsigned int from_dimension, unsigned int to_dimension>
struct Link_alphas_of_dart_according_to_lower_pattern {
  static void link(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::LDGF::Dart_handle current_dart) {
    if (!LCCE::LDGF::is_free(ldcc, current_dart, from_dimension-1)) LCCE::HDGF::template basic_link_alpha<from_dimension>(hdcc, layers[current_dart].back(), layers[LCCE::LDGF::alpha(current_dart, from_dimension-1)].back());
    Link_alphas_of_dart_according_to_lower_pattern<LCCE, from_dimension+1, to_dimension>::link(layers, ldcc, hdcc, current_dart);
  }
};

template <class LCCE, unsigned int to_dimension>
struct Link_alphas_of_dart_according_to_lower_pattern<LCCE, to_dimension, to_dimension> {
  static void link(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::LDGF::Dart_handle current_dart) {
    // Do nothing
  }
};

template <class LCCE, unsigned int dimension, unsigned int dimension_to>
struct Extrude_darts_top_with_range {
  static void create_darts_layer(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extruded_embeddings_with_range &ee, typename LCCE::FT current_range_min, typename LCCE::FT current_range_max, typename LCCE::Extrusion_ranges &er) {
    std::cout << "\tExtruding top (gmaps) at dimension " << dimension << "..." << std::endl;
    
    // Make a new darts layer connected by dimension+1 to the previous one
    typename LCCE::Lower_dimensional_cell_complex::Dart_range &ld_darts = ldcc.darts();
    std::cout << "\t\t[";
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_range::iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Check if the extrusion range of this dart is in the current range
      if (std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
        
        // Check if the minimum of the current range is in the extrusion range
        if (std::get<dimension>(er.ranges).ranges_map[current_dart->template attribute<dimension>()].is_in(current_range_max)) {
          typename LCCE::LDGF::Dart_handle current_dart_from = LCCE::LDGF::dart_from(current_dart);
          layers[current_dart_from].push_back(LCCE::HDGF::template create_dart_linked_by_alpha<dimension>(hdcc, layers[current_dart_from].back(), true));
          typename LCCE::LDGF::Dart_handle current_dart_to = LCCE::LDGF::dart_to(current_dart);
          layers[current_dart_to].push_back(LCCE::HDGF::template create_dart_linked_by_alpha<dimension>(hdcc, layers[current_dart_to].back(), false));
        } else {
          std::cout << "--";
        }
      } else {
        std::cout << "  ";
      }
    } std::cout << "]" << std::endl;
    
    // Set links and embeddings for the layer
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_range::iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Check if the extrusion range of this dart is in the current range
      if (std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
        
        // Check if the minimum of the current range is in the extrusion range
        if (std::get<dimension>(er.ranges).ranges_map[current_dart->template attribute<dimension>()].is_in(current_range_max)) {
          
          // Link darts within this layer
          // The orientation has already been defined by linking it to the previous layer when creating the darts
          Link_alphas_of_dart_according_to_pattern<LCCE, 1, dimension>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_from(current_dart));
          Link_alphas_of_dart_according_to_pattern<LCCE, 1, dimension>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_to(current_dart));
          if (!ldcc.is_free(current_dart, LCCE::Lower_dimensional_cell_complex::dimension) &&
              std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->beta(LCCE::Lower_dimensional_cell_complex::dimension)->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
            Link_alphas_of_dart_according_to_lower_pattern<LCCE, dimension+2, LCCE::Lower_dimensional_cell_complex::dimension+2>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_from(current_dart));
            Link_alphas_of_dart_according_to_lower_pattern<LCCE, dimension+2, LCCE::Lower_dimensional_cell_complex::dimension+2>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_to(current_dart));
          } else {
            Link_alphas_of_dart_according_to_lower_pattern<LCCE, dimension+2, LCCE::Lower_dimensional_cell_complex::dimension+1>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_from(current_dart));
            Link_alphas_of_dart_according_to_lower_pattern<LCCE, dimension+2, LCCE::Lower_dimensional_cell_complex::dimension+1>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_to(current_dart));
          }
          
          
          // Set correct embeddings
          Set_extruded_embeddings_top_with_range<LCCE, 0, LCCE::Higher_dimensional_cell_complex::dimension, dimension>::set(layers, hdcc, LCCE::LDGF::dart_from(current_dart), ee, current_range_min, current_range_max);
          Set_extruded_embeddings_top_with_range<LCCE, 0, LCCE::Higher_dimensional_cell_complex::dimension, dimension>::set(layers, hdcc, LCCE::LDGF::dart_to(current_dart), ee, current_range_min, current_range_max);
        }
      }
    }
    
    Extrude_darts_top_with_range<LCCE, dimension+1, LCCE::Lower_dimensional_cell_complex::dimension>::create_darts_layer(layers, ldcc, hdcc, ee, current_range_min, current_range_max, er);
  }
};

template <class LCCE, unsigned int dimension_to>
struct Extrude_darts_top_with_range<LCCE, 0, dimension_to> {
  static void create_darts_layer(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extruded_embeddings_with_range &ee, typename LCCE::FT current_range_min, typename LCCE::FT current_range_max, typename LCCE::Extrusion_ranges &er) {
    std::cout << "\tExtruding top (gmaps) at dimension 0..." << std::endl;
    
    // Make a new darts layer connected by dimension+1 to the previous one
    typename LCCE::Lower_dimensional_cell_complex::Dart_range &ld_darts = ldcc.darts();
    std::cout << "\t\t[";
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_range::iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Check if the extrusion range of this dart is in the current range
      if (std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
        
        // Check if the maximum of the current range is in the extrusion range
        if (std::get<0>(er.ranges).ranges_map[LCCE::LDGF::template attribute<0>(LCCE::LDGF::dart_from(current_dart))].is_in(current_range_max)) {
          typename LCCE::LDGF::Dart_handle current_dart_from = LCCE::LDGF::dart_from(current_dart);
          layers[current_dart_from].push_back(LCCE::HDGF::template create_dart_linked_by_alpha<0>(hdcc, layers[current_dart_from].back()));
        } else {
          std::cout << "-";
        }
        
        if (std::get<0>(er.ranges).ranges_map[LCCE::LDGF::template attribute<0>(LCCE::LDGF::dart_to(current_dart))].is_in(current_range_max)) {
          typename LCCE::LDGF::Dart_handle current_dart_to = LCCE::LDGF::dart_to(current_dart);
          layers[current_dart_to].push_back(LCCE::HDGF::template create_dart_linked_by_alpha<0>(hdcc, layers[current_dart_to].back()));
        } else {
          std::cout << "-";
        }
      } else {
        std::cout << "  ";
      }
    } std::cout << "]" << std::endl;
    
    // Set links and embeddings for the layer
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_range::iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Check if the extrusion range of this dart is in the current range
      if (std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
        
        // Check if the minimum of the current range is in the extrusion range
        if (std::get<0>(er.ranges).ranges_map[LCCE::LDGF::template attribute<0>(LCCE::LDGF::dart_from(current_dart))].is_in(current_range_max)) {
          
          // Linking darts within this layer is not needed
          
          // Set correct embeddings
          Set_extruded_embeddings_top_with_range<LCCE, 0, LCCE::Higher_dimensional_cell_complex::dimension, 0>::set(layers, hdcc, LCCE::LDGF::dart_from(current_dart), ee, current_range_min, current_range_max);
        }
        
        if (std::get<0>(er.ranges).ranges_map[LCCE::LDGF::template attribute<0>(LCCE::LDGF::dart_to(current_dart))].is_in(current_range_max)) {
          
          // Linking darts within this layer is not needed
          
          // Set correct embeddings
          Set_extruded_embeddings_top_with_range<LCCE, 0, LCCE::Higher_dimensional_cell_complex::dimension, 0>::set(layers, hdcc, LCCE::LDGF::dart_to(current_dart), ee, current_range_min, current_range_max);
        }
      }
    }
    
    Extrude_darts_top_with_range<LCCE, 1, LCCE::Lower_dimensional_cell_complex::dimension>::create_darts_layer(layers, ldcc, hdcc, ee, current_range_min, current_range_max, er);
  }
};

template <class LCCE, unsigned int dimension>
struct Extrude_darts_top_with_range<LCCE, dimension, dimension> {
  static void create_darts_layer(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extruded_embeddings_with_range &ee, typename LCCE::FT current_range_min, typename LCCE::FT current_range_max, typename LCCE::Extrusion_ranges &er) {
    std::cout << "\tExtruding top (gmaps) at dimension " << dimension << "..." << std::endl;
    
    // Make a new darts layer connected by dimension+1 to the previous one
    typename LCCE::Lower_dimensional_cell_complex::Dart_range &ld_darts = ldcc.darts();
    std::cout << "\t\t[";
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_range::iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Check if the extrusion range of this dart is in the current range
      if (std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
        
        // Check if the maximum of the current range is in the extrusion range
        if (std::get<dimension>(er.ranges).ranges_map[current_dart->template attribute<dimension>()].is_in(current_range_max)) {
          typename LCCE::LDGF::Dart_handle current_dart_from = LCCE::LDGF::dart_from(current_dart);
          layers[current_dart_from].push_back(LCCE::HDGF::template create_dart_linked_by_alpha<dimension>(hdcc, layers[current_dart_from].back()));
          typename LCCE::LDGF::Dart_handle current_dart_to = LCCE::LDGF::dart_to(current_dart);
          layers[current_dart_to].push_back(LCCE::HDGF::template create_dart_linked_by_alpha<dimension>(hdcc, layers[current_dart_to].back()));
        } else {
          std::cout << "--";
        }
      } else {
        std::cout << "  ";
      }
    } std::cout << "]" << std::endl;
    
    // Set links and embeddings for the layer
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_range::iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Check if the extrusion range of this dart is in the current range
      if (std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
        
        // Check if the maximum of the current range is in the extrusion range
        if (std::get<dimension>(er.ranges).ranges_map[current_dart->template attribute<dimension>()].is_in(current_range_max)) {
          
          // Link darts within this layer
          // The orientation has already been defined by linking it to the previous layer when creating the darts
          Link_alphas_of_dart_according_to_pattern<LCCE, 1, dimension>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_from(current_dart));
          Link_alphas_of_dart_according_to_pattern<LCCE, 1, dimension>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_to(current_dart));
          
          // Set correct embeddings
          Set_extruded_embeddings_top_with_range<LCCE, 0, LCCE::Higher_dimensional_cell_complex::dimension, dimension>::set(layers, hdcc, LCCE::LDGF::dart_from(current_dart), ee, current_range_min, current_range_max);
          Set_extruded_embeddings_top_with_range<LCCE, 0, LCCE::Higher_dimensional_cell_complex::dimension, dimension>::set(layers, hdcc, LCCE::LDGF::dart_to(current_dart), ee, current_range_min, current_range_max);
        }
      }
    }
    
  }
};

template <class LCCE, unsigned int dimension = LCCE::Lower_dimensional_cell_complex::dimension>
struct Extrude_darts_base_with_range {
  static void create_darts_layer(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extruded_embeddings_with_range &ee, typename LCCE::FT current_range_min, typename LCCE::FT current_range_max, typename LCCE::Extrusion_ranges &er) {
    std::cout << "\tExtruding base (gmaps) at dimension " << dimension << "..." << std::endl;
    
    // Make a new darts layer connected by dimension+1 to the previous one
    typename LCCE::Lower_dimensional_cell_complex::Dart_range &ld_darts = ldcc.darts();
    std::cout << "\t\t[";
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_range::iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Check if the extrusion range of this dart is in the current range
      if (std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
        
        // Check if the minimum of the current range is in the extrusion range
        if (std::get<dimension>(er.ranges).ranges_map[current_dart->template attribute<dimension>()].is_in(current_range_min)) {
          typename LCCE::LDGF::Dart_handle current_dart_from = LCCE::LDGF::dart_from(current_dart);
          layers[current_dart_from].push_back(LCCE::HDGF::template create_dart_linked_by_alpha<dimension+1>(hdcc, layers[current_dart_from].back()));
          typename LCCE::LDGF::Dart_handle current_dart_to = LCCE::LDGF::dart_to(current_dart);
          layers[current_dart_to].push_back(LCCE::HDGF::template create_dart_linked_by_alpha<dimension+1>(hdcc, layers[current_dart_to].back()));
        } else {
          std::cout << "--";
        }
      } else {
        std::cout << "  ";
      }
    } std::cout << "]" << std::endl;
    
    // Set links and embeddings for the layer
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_range::iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Check if the extrusion range of this dart is in the current range
      if (std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
        
        // Check if the minimum of the current range is in the extrusion range
        if (std::get<dimension>(er.ranges).ranges_map[current_dart->template attribute<dimension>()].is_in(current_range_min)) {
          
          // Link darts within this layer
          // The orientation has already been defined by linking it to the previous layer when creating the darts
          Link_alphas_of_dart_according_to_pattern<LCCE, 1, dimension>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_from(current_dart));
          Link_alphas_of_dart_according_to_pattern<LCCE, 1, dimension>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_to(current_dart));
          if (!ldcc.is_free(current_dart, LCCE::Lower_dimensional_cell_complex::dimension) &&
              std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->beta(LCCE::Lower_dimensional_cell_complex::dimension)->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
            Link_alphas_of_dart_according_to_lower_pattern<LCCE, dimension+2, LCCE::Lower_dimensional_cell_complex::dimension+2>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_from(current_dart));
            Link_alphas_of_dart_according_to_lower_pattern<LCCE, dimension+2, LCCE::Lower_dimensional_cell_complex::dimension+2>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_to(current_dart));
          } else {
            Link_alphas_of_dart_according_to_lower_pattern<LCCE, dimension+2, LCCE::Lower_dimensional_cell_complex::dimension+1>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_from(current_dart));
            Link_alphas_of_dart_according_to_lower_pattern<LCCE, dimension+2, LCCE::Lower_dimensional_cell_complex::dimension+1>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_to(current_dart));
          }
          
          // Set correct embeddings
          Set_extruded_embeddings_base_with_range<LCCE, 0, LCCE::Higher_dimensional_cell_complex::dimension, dimension>::set(layers, hdcc, LCCE::LDGF::dart_from(current_dart), ee, current_range_min, current_range_max);
          Set_extruded_embeddings_base_with_range<LCCE, 0, LCCE::Higher_dimensional_cell_complex::dimension, dimension>::set(layers, hdcc, LCCE::LDGF::dart_to(current_dart), ee, current_range_min, current_range_max);
        }
      }
    }
    
    Extrude_darts_base_with_range<LCCE, dimension-1>::create_darts_layer(layers, ldcc, hdcc, ee, current_range_min, current_range_max, er);
  }
};

template <class LCCE>
struct Extrude_darts_base_with_range<LCCE, 0> {
  static void create_darts_layer(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extruded_embeddings_with_range &ee, typename LCCE::FT current_range_min, typename LCCE::FT current_range_max, typename LCCE::Extrusion_ranges &er) {
    std::cout << "\tExtruding base (gmaps) at dimension 0..." << std::endl;
    
    // Make a new darts layer connected by dimension+1 to the previous one
    typename LCCE::Lower_dimensional_cell_complex::Dart_range &ld_darts = ldcc.darts();
    std::cout << "\t\t[";
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_range::iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Check if the extrusion range of this dart is in the current range
      if (std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
        
        // Check if the minimum of the current range is in the extrusion range
        if (std::get<0>(er.ranges).ranges_map[LCCE::LDGF::template attribute<0>(LCCE::LDGF::dart_from(current_dart))].is_in(current_range_min)) {
          typename LCCE::LDGF::Dart_handle current_dart_from = LCCE::LDGF::dart_from(current_dart);
          layers[current_dart_from].push_back(LCCE::HDGF::template create_dart_linked_by_alpha<1>(hdcc, layers[current_dart_from].back(), true));
        } else {
          std::cout << "-";
        }
        
        if (std::get<0>(er.ranges).ranges_map[LCCE::LDGF::template attribute<0>(LCCE::LDGF::dart_to(current_dart))].is_in(current_range_min)) {
          typename LCCE::LDGF::Dart_handle current_dart_to = LCCE::LDGF::dart_to(current_dart);
          layers[current_dart_to].push_back(LCCE::HDGF::template create_dart_linked_by_alpha<1>(hdcc, layers[current_dart_to].back(), true));
        } else {
          std::cout << "-";
        }
      } else {
        std::cout << "  ";
      }
    } std::cout << "]" << std::endl;
    
    // Set links and embeddings for the layer
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_range::iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Check if the extrusion range of this dart is in the current range
      if (std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
        
        // Check if the minimum of the current range is in the extrusion range
        if (std::get<0>(er.ranges).ranges_map[LCCE::LDGF::template attribute<0>(LCCE::LDGF::dart_from(current_dart))].is_in(current_range_min)) {
          
          // Link darts within this layer
          // The orientation has already been defined by linking it to the previous layer when creating the darts
          if (!ldcc.is_free(current_dart, LCCE::Lower_dimensional_cell_complex::dimension) &&
              std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->beta(LCCE::Lower_dimensional_cell_complex::dimension)->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
            Link_alphas_of_dart_according_to_lower_pattern<LCCE, 2, LCCE::Lower_dimensional_cell_complex::dimension+2>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_from(current_dart));
          } else {
            Link_alphas_of_dart_according_to_lower_pattern<LCCE, 2, LCCE::Lower_dimensional_cell_complex::dimension+1>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_to(current_dart));
          }
          
          // Set correct embeddings
          Set_extruded_embeddings_base_with_range<LCCE, 0, LCCE::Higher_dimensional_cell_complex::dimension, 0>::set(layers, hdcc, LCCE::LDGF::dart_from(current_dart), ee, current_range_min, current_range_max);
        }
        
        if (std::get<0>(er.ranges).ranges_map[LCCE::LDGF::template attribute<0>(LCCE::LDGF::dart_to(current_dart))].is_in(current_range_min)) {
          
          // Link darts within this layer
          // The orientation has already been defined by linking it to the previous layer when creating the darts
          if (!ldcc.is_free(current_dart, LCCE::Lower_dimensional_cell_complex::dimension) &&
              std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->beta(LCCE::Lower_dimensional_cell_complex::dimension)->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
            Link_alphas_of_dart_according_to_lower_pattern<LCCE, 2, LCCE::Lower_dimensional_cell_complex::dimension+2>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_to(current_dart));
          } else {
            Link_alphas_of_dart_according_to_lower_pattern<LCCE, 2, LCCE::Lower_dimensional_cell_complex::dimension+1>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_to(current_dart));
          }
          
          // Set correct embeddings
          Set_extruded_embeddings_base_with_range<LCCE, 0, LCCE::Higher_dimensional_cell_complex::dimension, 0>::set(layers, hdcc, LCCE::LDGF::dart_to(current_dart), ee, current_range_min, current_range_max);
        }
      }
    }
    
    Extrude_darts_top_with_range<LCCE, 0, LCCE::Lower_dimensional_cell_complex::dimension>::create_darts_layer(layers, ldcc, hdcc, ee, current_range_min, current_range_max, er);
  }
};

// Use combinatorial maps directly only for this layer
template <class LCCE, unsigned int dimension = LCCE::Lower_dimensional_cell_complex::dimension>
struct Extrude_darts_first_layer_with_range {
  static void create_darts_layer(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extruded_embeddings_with_range &ee, typename LCCE::FT current_range_min, typename LCCE::FT current_range_max, typename LCCE::Extrusion_ranges &er) {
    std::cout << "\tExtruding first layer (cmaps) at dimension " << dimension << "..." << std::endl;
    
    // Make a new darts layer
    typename LCCE::Lower_dimensional_cell_complex::Dart_range &ld_darts = ldcc.darts();
    std::cout << "\t\t[";
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_range::iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Check if the extrusion range of this dart is in the current range
      if (std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
        
        // Check if the minimum of the current range is in the extrusion range
        if (std::get<dimension>(er.ranges).ranges_map[current_dart->template attribute<dimension>()].is_in(current_range_min)) {
          typename LCCE::Higher_dimensional_cell_complex::Dart_handle new_dart = hdcc.create_dart(hdcc.darts().size());
          if (dimension % 2 == 0) {
            layers[LCCE::LDGF::dart_from(current_dart)].push_back(LCCE::HDGF::dart_to(new_dart));
            layers[LCCE::LDGF::dart_to(current_dart)].push_back(LCCE::HDGF::dart_from(new_dart));
          } else {
            layers[LCCE::LDGF::dart_from(current_dart)].push_back(LCCE::HDGF::dart_from(new_dart));
            layers[LCCE::LDGF::dart_to(current_dart)].push_back(LCCE::HDGF::dart_to(new_dart));
          } std::cout << "Dd";
        } else {
          std::cout << "--";
        }
      } else {
        std::cout << "  ";
      }
    } std::cout << "]" << std::endl;
    
    // Set links and embeddings for the layer
    for (typename LCCE::Lower_dimensional_cell_complex::Dart_range::iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      
      // Check if the extrusion range of this dart is in the current range
      if (std::get<LCCE::Lower_dimensional_cell_complex::dimension>(er.ranges).ranges_map[current_dart->template attribute<LCCE::Lower_dimensional_cell_complex::dimension>()].covers(current_range_min, current_range_max)) {
        
        // Check if the minimum of the current range is in the extrusion range
        if (std::get<dimension>(er.ranges).ranges_map[current_dart->template attribute<dimension>()].is_in(current_range_min)) {
          
          // Link darts within this layer
          // Note that this defines the orientation, so cmaps are used here!
          if (dimension % 2 == 0) hdcc.basic_link_beta_1(layers[LCCE::LDGF::dart_from(current_dart)].back().cmap_dart(), layers[LCCE::LDGF::dart_from(current_dart->beta(0))].back().cmap_dart());
          else hdcc.basic_link_beta_1(layers[LCCE::LDGF::dart_from(current_dart)].back().cmap_dart(), layers[LCCE::LDGF::dart_from(current_dart->beta(1))].back().cmap_dart());
          Link_alphas_of_dart_according_to_pattern<LCCE, 2, dimension>::link(layers, ldcc, hdcc, LCCE::LDGF::dart_from(current_dart));
          
          // Connect to previous layer if there is one just adjacent to it using beta_dimension+1
          // Note that we only need to use the cmap darts here
          if (std::get<dimension>(er.ranges).ranges_map[current_dart->template attribute<dimension>()].is_surrounded(current_range_min)) {
            LCCE::HDGF::template basic_link_alpha<dimension+1>(hdcc, layers[LCCE::LDGF::dart_from(current_dart)].back(), layers[LCCE::LDGF::dart_from(current_dart)][layers[LCCE::LDGF::dart_from(current_dart)].size()-2]);
          }
          
          // Set correct embeddings
          if (dimension % 2 == 0) {
            hdcc.template set_attribute<0>(layers[LCCE::LDGF::dart_from(current_dart)].back().cmap_dart(), std::get<0>(ee).embeddings[current_dart->beta(1)->template attribute<0>()].nonex[current_range_min]);
            std::get<0>(ee).embeddings[current_dart->beta(1)->template attribute<0>()].nonex[current_range_min]->set_dart(layers[LCCE::LDGF::dart_from(current_dart)].back().cmap_dart());
          } else {
            hdcc.template set_attribute<0>(layers[LCCE::LDGF::dart_from(current_dart)].back().cmap_dart(), std::get<0>(ee).embeddings[current_dart->template attribute<0>()].nonex[current_range_min]);
            std::get<0>(ee).embeddings[current_dart->template attribute<0>()].nonex[current_range_min]->set_dart(layers[LCCE::LDGF::dart_from(current_dart)].back().cmap_dart());
          } Set_extruded_embeddings_base_with_range<LCCE, 1, LCCE::Higher_dimensional_cell_complex::dimension, dimension>::set(layers, hdcc, LCCE::LDGF::dart_from(current_dart), ee, current_range_min, current_range_max);
        }
      }
    }
    
    Extrude_darts_base_with_range<LCCE, dimension-1>::create_darts_layer(layers, ldcc, hdcc, ee, current_range_min, current_range_max, er);
  }
};

template <class LCCE, unsigned int dimension = LCCE::Lower_dimensional_cell_complex::dimension>
struct Extrude_darts_with_range {
  static void extrude(std::map<typename LCCE::LDGF::Dart_handle, std::vector<typename LCCE::HDGF::Dart_handle>, typename LCCE::LDGF::Dart_handle::Dart_handle_comparator> &layers, typename LCCE::Lower_dimensional_cell_complex &ldcc, typename LCCE::Higher_dimensional_cell_complex &hdcc, typename LCCE::Extruded_embeddings_with_range &ee, std::set<typename LCCE::FT> &all_ranges, typename LCCE::Extrusion_ranges &er) {
    
    // Extrude range by range
    for (typename std::set<typename LCCE::FT>::const_iterator current_value = all_ranges.begin(); current_value != all_ranges.end(); ++current_value) {
      typename std::set<typename LCCE::FT>::const_iterator next_value = current_value;
      ++next_value;
      if (next_value == all_ranges.end()) break;
      std::cout << "Extruding at [" << *current_value << "," << *next_value << "]" << std::endl;
      Extrude_darts_first_layer_with_range<LCCE>::create_darts_layer(layers, ldcc, hdcc, ee, *current_value, *next_value, er);
      //      break;
    }
  }
};

template <unsigned int unextruded_dimension>
class Linear_cell_complex_extruder_with_range {
public:
  typedef typename Linear_cell_complex_with_ids<unextruded_dimension>::type Lower_dimensional_cell_complex;
  typedef typename Linear_cell_complex_with_ids<unextruded_dimension+1>::type Higher_dimensional_cell_complex;
  typedef Linear_cell_complex_extruder_with_range<unextruded_dimension> Self;
  typedef typename Lower_dimensional_cell_complex::FT FT;
  typedef Gmap_faker<Higher_dimensional_cell_complex> HDGF;
  typedef Gmap_faker<Lower_dimensional_cell_complex> LDGF;
  
  typedef typename Extruded_embeddings_with_range_of_dimension_and_lower<Self>::type Extruded_embeddings_with_range;
  typedef Extrusion_ranges_tuple_per_dimension<Lower_dimensional_cell_complex> Extrusion_ranges;
  
private:
  Extrusion_ranges extrusion_ranges;
  std::set<typename Lower_dimensional_cell_complex::FT> all_ranges;
  Higher_dimensional_cell_complex hdcc;
  
public:
  void load_ranges_file(const char *ranges_path, Lower_dimensional_cell_complex &ldcc) {
    Extrusion_ranges_reader<Lower_dimensional_cell_complex>::load_from_ranges_file(ranges_path, ldcc, extrusion_ranges, all_ranges);
  }
  
  void load_ranges_from_ogr(const char *ranges_path, Lower_dimensional_cell_complex &ldcc, const char *min_field_name, const char *max_field_name) {
    Extrusion_ranges_reader<Lower_dimensional_cell_complex>::load_from_ogr(ranges_path, ldcc, extrusion_ranges, all_ranges, min_field_name, max_field_name);
  }
  
  Higher_dimensional_cell_complex extrude_using_ranges(Lower_dimensional_cell_complex &ldcc) {
    CGAL_precondition(Lower_dimensional_cell_complex::dimension > 1);
    
    // Spread the ranges to the lower dimensional cells (adjacency)
    Pass_extrusion_ranges_to_boundary<Extrusion_ranges>::pass(ldcc, extrusion_ranges);
    Print_all_extrusion_ranges<Extrusion_ranges>::print(extrusion_ranges);
    
    // Extrude the embeddings
    Extruded_embeddings_with_range extruded_embeddings;
    Extrude_embeddings_with_range_of_dimension_or_lower<Self>::extrude(extruded_embeddings, ldcc, hdcc, extrusion_ranges);
    
    // Extrude the darts
    std::map<typename LDGF::Dart_handle, std::vector<typename HDGF::Dart_handle>, typename LDGF::Dart_handle::Dart_handle_comparator> layers;
    typename Lower_dimensional_cell_complex::Dart_range &ld_darts = ldcc.darts();
    for (typename Lower_dimensional_cell_complex::Dart_range::iterator current_dart = ld_darts.begin();
         current_dart != ld_darts.end();
         ++current_dart) {
      typename LDGF::Dart_handle d = LDGF::dart_from(current_dart);
      layers[d] = std::vector<typename HDGF::Dart_handle>();
      d = LDGF::dart_to(current_dart);
      layers[d] = std::vector<typename HDGF::Dart_handle>();
    } Extrude_darts_with_range<Self>::extrude(layers, ldcc, hdcc, extruded_embeddings, all_ranges, extrusion_ranges);
    
    return hdcc;
  }
};

#endif