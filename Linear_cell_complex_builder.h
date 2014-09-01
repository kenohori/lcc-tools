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

#ifndef LCC_BUILDER_H
#define LCC_BUILDER_H

#include "Linear_cell_complex_cell_comparator.h"
#include "Linear_cell_complex_cell_cloner.h"
#include "Linear_cell_complex_validator.h"

template <class LCC>
class Call_if_lcc_with_ids {
public:
  static typename LCC::Dart_handle create_temporary_dart(LCC &lcc, typename LCC::Vertex_attribute_handle vah) {
    return lcc.create_dart(vah);
  }
  
  static int make_temporary_darts_permanent(LCC &lcc) {
    return -1;
  }
};

template <unsigned int d>
class Call_if_lcc_with_ids<Linear_cell_complex_with_ids<d> > {
  typedef Linear_cell_complex_with_ids<d> LCC;
public:
  static typename LCC::Dart_handle create_temporary_dart(LCC &lcc, typename LCC::Vertex_attribute_handle vah) {
    return lcc.create_temporary_dart(vah);
  }
  
  static int make_temporary_darts_permanent(LCC &lcc) {
    return lcc.make_temporary_darts_permanent();
  }
};

template <class LCC_>
class Linear_cell_complex_builder {
public:
  typedef LCC_ LCC;
  typedef typename LCC::Dart_handle Dart_handle;
  typedef typename LCC::Point Point;
  
protected:
  LCC &lcc;
  
public:
  Linear_cell_complex_builder(LCC &alcc) : lcc(alcc) {

  }
  
  Dart_handle get_vertex(const Point &p) {
    
    // Brute-force search
    for (typename LCC::Dart_range::iterator current_dart = lcc.darts().begin(), last_dart = lcc.darts().end(); current_dart != last_dart; ++current_dart) {
      if (lcc.point(current_dart) == p) {
        return current_dart;
      }
    }
    
    // Not found, create a new one
    Dart_handle new_dart = lcc.create_dart(p);
    
    // Validation
//    CGAL_postcondition(Linear_cell_complex_validator<LCC>::are_dart_ids_valid(lcc));
//    CGAL_postcondition(Linear_cell_complex_validator<LCC>::are_vertices_unique(lcc));
//    CGAL_postcondition(Linear_cell_complex_validator<LCC>::template are_attribute_ids_valid<0>(lcc));
    
    return new_dart;
  }
  
  template <class InputIterator>
  Dart_handle get_face(InputIterator begin, InputIterator end) {
    
    // Validate input
    InputIterator last = end;
    --last;
    CGAL_precondition(*begin == *last); // check that it's a closed face
    
    // Use free darts, make new (temporary) darts for those that aren't free
    std::list<Dart_handle> free_vertices, temporary_vertices;
    for (InputIterator current = begin; current != last; ++current) {
      if (lcc.is_free(*current, 0) && lcc.is_free(*current, 1)) free_vertices.push_back(*current);
      else {
        free_vertices.push_back(Call_if_lcc_with_ids<LCC>::create_temporary_dart(lcc, (*current)->template attribute<0>()));
        temporary_vertices.push_back(free_vertices.back());
      }
    }
    
    // Create temporary structure
    typename std::list<Dart_handle>::iterator current_vertex = free_vertices.begin();
    typename std::list<Dart_handle>::iterator next_vertex = current_vertex;
    ++next_vertex;
    while (next_vertex != free_vertices.end()) {
      lcc.link_beta_1(*current_vertex, *next_vertex);
      ++current_vertex;
      ++next_vertex;
    } lcc.link_beta_1(free_vertices.back(), free_vertices.front());
    
    // Brute-force search
    Dart_handle found_face_brute_force = LCC::null_handle;
    for (typename LCC::Dart_range::iterator current_dart = lcc.darts().begin(), last_dart = lcc.darts().end(); current_dart != last_dart; ++current_dart) {
      if (&*current_dart == &*free_vertices.front()) continue; // TODO: Compare by memory address for now
      if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<2>(lcc, current_dart, free_vertices.front()) || Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<2>(lcc, current_dart, free_vertices.front())) {
        CGAL_precondition(free_vertices.size() == temporary_vertices.size()); // TODO: If always true (which I think it is), this can be used as an optimisation
        for (current_vertex = temporary_vertices.begin(); current_vertex != temporary_vertices.end(); ++current_vertex) {
          lcc.erase_dart(*current_vertex);
        } found_face_brute_force = current_dart;
        break;
      }
    }
    
    // Found, return existing
    if (found_face_brute_force != LCC::null_handle) {
      
      // Validation
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::are_dart_ids_valid(lcc));
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::template are_attribute_ids_valid<2>(lcc));
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<0>(lcc);
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<2>(lcc);
      
      return found_face_brute_force;
    }
    
    // Not found, create a new one
    else {
      typename LCC::template Attribute_handle<2>::type facet_attribute = lcc.template create_attribute<2>();
      lcc.template set_attribute<2>(free_vertices.front(), facet_attribute);
      int changed_darts = Call_if_lcc_with_ids<LCC>::make_temporary_darts_permanent(lcc);
      
      // Validation
      CGAL_postcondition(changed_darts == -1 || changed_darts == temporary_vertices.size());
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::are_dart_ids_valid(lcc));
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::template are_attribute_ids_valid<2>(lcc));
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<0>(lcc);
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<2>(lcc);
      
      return free_vertices.front();
    }
  }
  
  template <unsigned int dimension, class InputIterator>
  Dart_handle get_cell(InputIterator begin, InputIterator end) {
    CGAL_precondition(dimension > 2);
    unsigned int number_of_facets = 0;
    for (InputIterator current = begin; current != end; ++current) {
      ++number_of_facets;
    }
//    std::cout << "get_cell<" << dimension << ">(" << number_of_facets << " " << dimension-1 << "-cells)" << std::endl;
    
    // Validate that all the facets are closed dimension-1 cells
    for (InputIterator current = begin; current != end; ++current) {
      for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart = lcc.template darts_of_cell<dimension-1, dimension-1>(*current).begin(); current_dart != lcc.template darts_of_cell<dimension-1, dimension-1>(*current).end(); ++current_dart) {
        for (unsigned int b = 0; b < dimension-1; ++b) {
          CGAL_precondition(!lcc.is_free(current_dart, b));
        }
      }
    }
    
    // Use free facets, make new (temporary) facets for those that aren't free
    std::list<Dart_handle> free_facets, temporary_facets;
    for (InputIterator current = begin; current != end; ++current) {
      bool all_free = true, all_not_free = true;
      for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart = lcc.template darts_of_cell<dimension-1, dimension-1>(*current).begin(); current_dart != lcc.template darts_of_cell<dimension-1, dimension-1>(*current).end(); ++current_dart) {
        if (lcc.is_free(current_dart, dimension-1)) all_not_free = false;
        else all_free = false;
      } CGAL_precondition(all_free != all_not_free);
      if (all_free) {
//        std::cout << ".";
        free_facets.push_back(*current);
      } else {
//        std::cout << "c";
        free_facets.push_back(Linear_cell_complex_cell_cloner<LCC>::template clone_cell_with_reversed_orientation<dimension-1>(lcc, *current));
        temporary_facets.push_back(free_facets.back());
      }
    }
//    std::cout << std::endl;
    
    // Validate that all the facets are beta(dimension-1)-free
    for (typename std::list<Dart_handle>::iterator dart_in_current_facet = free_facets.begin(); dart_in_current_facet != free_facets.end(); ++dart_in_current_facet) {
      for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart_in_facet = lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_facet).begin(); current_dart_in_facet != lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_facet).end(); ++current_dart_in_facet) {
        CGAL_precondition(lcc.is_free(current_dart_in_facet, dimension-1));
      }
    }
    
    // Validate that this will create a closed (topologically open) quasi-manifold cell
    for (typename std::list<Dart_handle>::iterator dart_in_current_facet_1 = free_facets.begin(); dart_in_current_facet_1 != free_facets.end(); ++dart_in_current_facet_1) {
      for (typename LCC::template One_dart_per_incident_cell_range<dimension-2, dimension-1, dimension-1>::iterator dart_in_current_ridge_1 = lcc.template one_dart_per_incident_cell<dimension-2, dimension-1, dimension-1>(*dart_in_current_facet_1).begin(); dart_in_current_ridge_1 != lcc.template one_dart_per_incident_cell<dimension-2, dimension-1, dimension-1>(*dart_in_current_facet_1).end(); ++dart_in_current_ridge_1) {
        int matches = 0;
        
        for (typename std::list<Dart_handle>::iterator dart_in_current_facet_2 = free_facets.begin(); dart_in_current_facet_2 != free_facets.end(); ++dart_in_current_facet_2) {
          if (&*dart_in_current_facet_1 == &*dart_in_current_facet_2) continue; // TODO: Compare by memory address for now
          for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart_in_facet_2 = lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_facet_2).begin(); current_dart_in_facet_2 != lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_facet_2).end(); ++current_dart_in_facet_2) {
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<dimension-2>(lcc,dart_in_current_ridge_1, current_dart_in_facet_2)) ++matches;
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<dimension-2>(lcc, dart_in_current_ridge_1, current_dart_in_facet_2)) ++matches;
          }
        }
        
        CGAL_precondition(matches == 1);
      }
    }
    
    // Create temporary structure
    for (typename std::list<Dart_handle>::iterator dart_in_current_facet_1 = free_facets.begin(); dart_in_current_facet_1 != free_facets.end(); ++dart_in_current_facet_1) {
      std::list<Dart_handle> ridges_in_current_facet_1;
      for (typename LCC::template One_dart_per_incident_cell_range<dimension-2, dimension-1, dimension-1>::iterator dart_in_current_ridge_1 = lcc.template one_dart_per_incident_cell<dimension-2, dimension-1, dimension-1>(*dart_in_current_facet_1).begin(); dart_in_current_ridge_1 != lcc.template one_dart_per_incident_cell<dimension-2, dimension-1, dimension-1>(*dart_in_current_facet_1).end(); ++dart_in_current_ridge_1) {
        ridges_in_current_facet_1.push_back(dart_in_current_ridge_1);
      } for (typename std::list<Dart_handle>::iterator dart_in_current_ridge_1 = ridges_in_current_facet_1.begin(); dart_in_current_ridge_1 != ridges_in_current_facet_1.end(); ++dart_in_current_ridge_1) {
        bool matched_ridge = !lcc.is_free(*dart_in_current_ridge_1, dimension-1);
        
        for (typename std::list<Dart_handle>::iterator dart_in_current_facet_2 = dart_in_current_facet_1; !matched_ridge && dart_in_current_facet_2 != free_facets.end(); ++dart_in_current_facet_2) {
          if (&*dart_in_current_facet_1 == &*dart_in_current_facet_2) continue; // TODO: Compare by memory address for now
          for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart_in_facet_2 = lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_facet_2).begin(), last_dart_in_facet_2 = lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_facet_2).end(); !matched_ridge && current_dart_in_facet_2 != last_dart_in_facet_2; ++current_dart_in_facet_2) {
            
            // Match with same orientation
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<dimension-2>(lcc, *dart_in_current_ridge_1, current_dart_in_facet_2)) {
              matched_ridge = true;
              CGAL_assertion(lcc.vertex_attribute(*dart_in_current_ridge_1) == lcc.vertex_attribute(current_dart_in_facet_2));
              lcc.reverse_orientation_connected_component(current_dart_in_facet_2);
              CGAL_assertion(lcc.vertex_attribute(*dart_in_current_ridge_1) == lcc.vertex_attribute(current_dart_in_facet_2->beta(1)));
              CGAL_assertion(lcc.vertex_attribute((*dart_in_current_ridge_1)->beta(1)) == lcc.vertex_attribute(current_dart_in_facet_2));
              CGAL_assertion(lcc.template is_sewable<dimension-1>(*dart_in_current_ridge_1, current_dart_in_facet_2));
              lcc.template sew<dimension-1>(*dart_in_current_ridge_1, current_dart_in_facet_2);
            }
            
            // Match with opposite orientation
            else if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<dimension-2>(lcc, *dart_in_current_ridge_1, current_dart_in_facet_2)) {
              matched_ridge = true;
              CGAL_assertion(lcc.vertex_attribute(*dart_in_current_ridge_1) == lcc.vertex_attribute(current_dart_in_facet_2->beta(1)));
              CGAL_assertion(lcc.vertex_attribute((*dart_in_current_ridge_1)->beta(1)) == lcc.vertex_attribute(current_dart_in_facet_2));
              CGAL_assertion(lcc.template is_sewable<dimension-1>(*dart_in_current_ridge_1, current_dart_in_facet_2));
              lcc.template sew<dimension-1>(*dart_in_current_ridge_1, current_dart_in_facet_2);
            }
            
          }
        }
        
      }
    }
    
    // Validate the created temporary structure
    for (typename std::list<Dart_handle>::iterator dart_in_current_facet = free_facets.begin(); dart_in_current_facet != free_facets.end(); ++dart_in_current_facet) {
      for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart_in_facet = lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_facet).begin(); current_dart_in_facet != lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_facet).end(); ++current_dart_in_facet) {
        CGAL_assertion(!lcc.is_free(current_dart_in_facet, dimension-1));
        CGAL_assertion(current_dart_in_facet->beta(dimension-1) != current_dart_in_facet);
        CGAL_assertion(current_dart_in_facet->beta(dimension-1)->beta(dimension-1) == current_dart_in_facet);
      }
    }
    
    // Brute-force search
    Dart_handle found_face_brute_force = LCC::null_handle;
    for (typename LCC::Dart_range::iterator current_dart = lcc.darts().begin(), last_dart = lcc.darts().end(); current_dart != last_dart; ++current_dart) {
      if (&*current_dart == &*free_facets.front()) continue; // TODO: Compare by memory address for now
      if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<dimension>(lcc, current_dart, free_facets.front()) || Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<dimension>(lcc, current_dart, free_facets.front())) {
        CGAL_assertion(free_facets.size() == temporary_facets.size()); // TODO: If always true (which I think it is), this can be used as an optimisation
        for (typename std::list<Dart_handle>::iterator current_facet = temporary_facets.begin(); current_facet != temporary_facets.end(); ++current_facet) {
          std::list<Dart_handle> darts_to_erase;
          for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart_in_facet = lcc.template darts_of_cell<dimension-1, dimension-1>(*current_facet).begin(); current_dart_in_facet != lcc.template darts_of_cell<dimension-1, dimension-1>(*current_facet).end(); ++current_dart_in_facet) {
            darts_to_erase.push_back(current_dart_in_facet);
          } for (typename std::list<Dart_handle>::iterator current_dart_to_erase = darts_to_erase.begin(); current_dart_to_erase != darts_to_erase.end(); ++current_dart_to_erase) {
            lcc.erase_dart(*current_dart_to_erase);
          }
        } found_face_brute_force = current_dart;
      }
    }
    
    // Found, return existing
    if (found_face_brute_force != LCC::null_handle) {
      
      // Validation
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::are_dart_ids_valid(lcc));
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::template are_attribute_ids_valid<dimension>(lcc));
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<0>(lcc);
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<dimension-1>(lcc);
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<dimension>(lcc);
      
      return found_face_brute_force;
    }
  
    // Not found, create a new one
    // TODO: When to create the 1-attributes?
    else {
      typename LCC::template Attribute_handle<dimension>::type cell_attribute = lcc.template create_attribute<dimension>();
      lcc.template set_attribute<dimension>(free_facets.front(), cell_attribute);
      unsigned int temporary_darts = 0;
      for (typename std::list<Dart_handle>::iterator current_facet = temporary_facets.begin(); current_facet != temporary_facets.end(); ++current_facet) {
        for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart_in_facet = lcc.template darts_of_cell<dimension-1, dimension-1>(*current_facet).begin(); current_dart_in_facet != lcc.template darts_of_cell<dimension-1, dimension-1>(*current_facet).end(); ++current_dart_in_facet) {
          ++temporary_darts;
        }
      } int changed_darts = Call_if_lcc_with_ids<LCC>::make_temporary_darts_permanent(lcc);
      
      // Validation
      CGAL_postcondition(changed_darts == -1 || changed_darts == temporary_darts);
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::are_dart_ids_valid(lcc));
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::template are_attribute_ids_valid<dimension>(lcc));
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<0>(lcc);
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<dimension-1>(lcc);
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<dimension>(lcc);
      
      return free_facets.front();
    }
  }
  
  template <unsigned int dimension, class InputIterator>
  Dart_handle link_cells(InputIterator begin, InputIterator end) {
    CGAL_precondition(dimension >= 2);
    unsigned int number_of_facets = 0;
    for (InputIterator current = begin; current != end; ++current) {
      ++number_of_facets;
    } std::cout << "link_cells<" << dimension << ">(" << number_of_facets << " " << dimension << "-cells)" << std::endl;
    
    // Validate that all the facets are closed dimension-cells, unique and are all beta(dimension)-free
    int mark = this->lcc.get_new_mark();
    CGAL_precondition(this->lcc.is_whole_map_unmarked(mark));
    for (InputIterator current = begin; current != end; ++current) {
      for (typename LCC::template Dart_of_cell_range<dimension, dimension>::iterator current_dart = this->lcc.template darts_of_cell<dimension, dimension>(*current).begin(); current_dart != this->lcc.template darts_of_cell<dimension, dimension>(*current).end(); ++current_dart) {
        CGAL_precondition(!this->lcc.is_marked(current_dart, mark));
        this->lcc.mark(current_dart, mark);
        for (unsigned int b = 0; b < dimension; ++b) {
          CGAL_precondition(!this->lcc.is_free(current_dart, b));
        } CGAL_precondition(this->lcc.is_free(current_dart, dimension));
      }
    } this->lcc.negate_mark(mark);
    CGAL_postcondition(this->lcc.is_whole_map_unmarked(mark));
    this->lcc.free_mark(mark);
    
    // Validate that this will create a quasi-manifold cell complex
    for (InputIterator dart_in_current_facet_1 = begin; dart_in_current_facet_1 != end; ++dart_in_current_facet_1) {
      for (typename LCC::template One_dart_per_incident_cell_range<dimension-1, dimension, dimension>::iterator dart_in_current_ridge_1 = lcc.template one_dart_per_incident_cell<dimension-1, dimension, dimension>(*dart_in_current_facet_1).begin(); dart_in_current_ridge_1 != lcc.template one_dart_per_incident_cell<dimension-1, dimension, dimension>(*dart_in_current_facet_1).end(); ++dart_in_current_ridge_1) {
        int matches = 0;
        
        for (InputIterator dart_in_current_facet_2 = begin; dart_in_current_facet_2 != end; ++dart_in_current_facet_2) {
          if (&*dart_in_current_facet_1 == &*dart_in_current_facet_2) continue; // TODO: Compare by memory address for now
          for (typename LCC::template Dart_of_cell_range<dimension, dimension>::iterator current_dart_in_facet_2 = lcc.template darts_of_cell<dimension, dimension>(*dart_in_current_facet_2).begin(); current_dart_in_facet_2 != lcc.template darts_of_cell<dimension, dimension>(*dart_in_current_facet_2).end(); ++current_dart_in_facet_2) {
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<dimension-1>(lcc,dart_in_current_ridge_1, current_dart_in_facet_2)) ++matches;
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<dimension-1>(lcc, dart_in_current_ridge_1, current_dart_in_facet_2)) ++matches;
          }
        }
        
        CGAL_precondition(matches <= 1);
      }
    }
    
    // Create structure
    for (InputIterator dart_in_current_facet_1 = begin; dart_in_current_facet_1 != end; ++dart_in_current_facet_1) {
      std::list<Dart_handle> ridges_in_current_facet_1;
      for (typename LCC::template One_dart_per_incident_cell_range<dimension-1, dimension, dimension>::iterator dart_in_current_ridge_1 = lcc.template one_dart_per_incident_cell<dimension-1, dimension, dimension>(*dart_in_current_facet_1).begin(); dart_in_current_ridge_1 != lcc.template one_dart_per_incident_cell<dimension-1, dimension, dimension>(*dart_in_current_facet_1).end(); ++dart_in_current_ridge_1) {
        ridges_in_current_facet_1.push_back(dart_in_current_ridge_1);
      } for (typename std::list<Dart_handle>::iterator dart_in_current_ridge_1 = ridges_in_current_facet_1.begin(); dart_in_current_ridge_1 != ridges_in_current_facet_1.end(); ++dart_in_current_ridge_1) {
        bool matched_ridge = !lcc.is_free(*dart_in_current_ridge_1, dimension);
        
        for (InputIterator dart_in_current_facet_2 = dart_in_current_facet_1; !matched_ridge && dart_in_current_facet_2 != end; ++dart_in_current_facet_2) {
          if (&*dart_in_current_facet_1 == &*dart_in_current_facet_2) continue; // TODO: Compare by memory address for now
          for (typename LCC::template Dart_of_cell_range<dimension, dimension>::iterator current_dart_in_facet_2 = lcc.template darts_of_cell<dimension, dimension>(*dart_in_current_facet_2).begin(), last_dart_in_facet_2 = lcc.template darts_of_cell<dimension, dimension>(*dart_in_current_facet_2).end(); !matched_ridge && current_dart_in_facet_2 != last_dart_in_facet_2; ++current_dart_in_facet_2) {
            
            // Match with same orientation
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<dimension-1>(lcc, *dart_in_current_ridge_1, current_dart_in_facet_2)) {
              matched_ridge = true;
              CGAL_assertion(lcc.vertex_attribute(*dart_in_current_ridge_1) == lcc.vertex_attribute(current_dart_in_facet_2));
              lcc.reverse_orientation_connected_component(current_dart_in_facet_2);
              CGAL_assertion(lcc.vertex_attribute(*dart_in_current_ridge_1) == lcc.vertex_attribute(current_dart_in_facet_2->beta(1)));
              CGAL_assertion(lcc.vertex_attribute((*dart_in_current_ridge_1)->beta(1)) == lcc.vertex_attribute(current_dart_in_facet_2));
              CGAL_assertion(lcc.template is_sewable<dimension>(*dart_in_current_ridge_1, current_dart_in_facet_2));
              lcc.template sew<dimension>(*dart_in_current_ridge_1, current_dart_in_facet_2);
            }
            
            // Match with opposite orientation
            else if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<dimension-1>(lcc, *dart_in_current_ridge_1, current_dart_in_facet_2)) {
              matched_ridge = true;
              CGAL_assertion(lcc.vertex_attribute(*dart_in_current_ridge_1) == lcc.vertex_attribute(current_dart_in_facet_2->beta(1)));
              CGAL_assertion(lcc.vertex_attribute((*dart_in_current_ridge_1)->beta(1)) == lcc.vertex_attribute(current_dart_in_facet_2));
              CGAL_assertion(lcc.template is_sewable<dimension>(*dart_in_current_ridge_1, current_dart_in_facet_2));
              lcc.template sew<dimension>(*dart_in_current_ridge_1, current_dart_in_facet_2);
            }
            
          }
        }
        
      }
    }
    
    // Validate the created structure
    for (InputIterator dart_in_current_facet = begin; dart_in_current_facet != end; ++dart_in_current_facet) {
      for (typename LCC::template Dart_of_cell_range<dimension, dimension>::iterator current_dart_in_facet = lcc.template darts_of_cell<dimension, dimension>(*dart_in_current_facet).begin(); current_dart_in_facet != lcc.template darts_of_cell<dimension, dimension>(*dart_in_current_facet).end(); ++current_dart_in_facet) {
        if (lcc.is_free(current_dart_in_facet, dimension)) continue;
        CGAL_assertion(current_dart_in_facet->beta(dimension) != current_dart_in_facet);
        CGAL_assertion(current_dart_in_facet->beta(dimension)->beta(dimension) == current_dart_in_facet);
      }
    } CGAL_postcondition(Linear_cell_complex_validator<LCC>::are_dart_ids_valid(lcc));
    CGAL_postcondition(Linear_cell_complex_validator<LCC>::template are_attribute_ids_valid<dimension>(lcc));
    Linear_cell_complex_validator<LCC>::template validate_map_attributes<0>(lcc);
    Linear_cell_complex_validator<LCC>::template validate_map_attributes<dimension-1>(lcc);
    Linear_cell_complex_validator<LCC>::template validate_map_attributes<dimension>(lcc);
    
    return *begin;
  }
};

#endif