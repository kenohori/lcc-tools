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

#ifndef LCC_BUILDER_WITH_INDEX_H
#define LCC_BUILDER_WITH_INDEX_H

#include "Linear_cell_complex_builder.h"
#include "Linear_cell_complex_smallest_vertex_index.h"

template <class LCC_>
class Linear_cell_complex_builder_with_index : public Linear_cell_complex_builder<LCC_> {
public:
  typedef Linear_cell_complex_builder<LCC_> Base;
  typedef LCC_ LCC;
  typedef typename LCC::Dart_handle Dart_handle;
  typedef typename LCC::Point Point;
  typedef Linear_cell_complex_smallest_vertex_index<LCC> Index;
  
protected:
  Index index;
  
public:
  Linear_cell_complex_builder_with_index(LCC &alcc) : Linear_cell_complex_builder<LCC>(alcc) {

  }
  
  Dart_handle get_vertex(const Point &p) {
    
    // Validate index
//    index.template validate_index<0>(this->lcc);
    
    // Index search
    Dart_handle result = index.find_vertex(p);
    if (result != LCC::null_handle) return result;
    
    // Not found, create a new one
    Dart_handle new_dart = this->lcc.create_dart(p);
    index.insert_vertex(this->lcc, new_dart);
    
    // Validation
//    CGAL_postcondition(Linear_cell_complex_validator<LCC>::are_dart_ids_valid(this->lcc));
//    CGAL_postcondition(Linear_cell_complex_validator<LCC>::are_vertices_unique(this->lcc));
//    CGAL_postcondition(Linear_cell_complex_validator<LCC>::template are_attribute_ids_valid<0>(this->lcc));
    
    return new_dart;
  }
  
  template <class InputIterator>
  Dart_handle get_face(InputIterator begin, InputIterator end) {
    
    // Validate input
    InputIterator last = end;
    --last;
    CGAL_precondition(*begin == *last); // check that it's a closed face
    
    // Validate index
//    index.template validate_index<2>(this->lcc);
    
    // Use free darts, make new (temporary) darts for those that aren't free
    std::list<Dart_handle> free_vertices, temporary_vertices;
    for (InputIterator current = begin; current != last; ++current) {
      if (this->lcc.is_free(*current, 0) && this->lcc.is_free(*current, 1)) free_vertices.push_back(*current);
      else {
        free_vertices.push_back(Call_if_lcc_with_ids<LCC>::create_temporary_dart(this->lcc, (*current)->template attribute<0>()));
        temporary_vertices.push_back(free_vertices.back());
      }
    }
    
    // Create temporary structure
    typename std::list<Dart_handle>::iterator current_vertex = free_vertices.begin();
    typename std::list<Dart_handle>::iterator next_vertex = current_vertex;
    ++next_vertex;
    while (next_vertex != free_vertices.end()) {
      this->lcc.basic_link_beta_1(*current_vertex, *next_vertex);
      ++current_vertex;
      ++next_vertex;
    } this->lcc.basic_link_beta_1(free_vertices.back(), free_vertices.front());
    
    // Index search
    Dart_handle found_face_with_index = LCC::null_handle;
    Dart_handle smallest_vertex = Index::template find_smallest_vertex<2>(this->lcc, free_vertices.front());
    typename Index::Results &results = index.template find_cells<2>(this->lcc.point(smallest_vertex));
    for (typename Index::Result_iterator current_face = results.begin(); current_face != results.end(); ++current_face) {
      CGAL_assertion(this->lcc.point(*current_face) == this->lcc.point(smallest_vertex));
      if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<2>(this->lcc, *current_face, smallest_vertex) || Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<2>(this->lcc, (*current_face)->beta(0), smallest_vertex)) {
        CGAL_precondition(free_vertices.size() == temporary_vertices.size());
        for (current_vertex = temporary_vertices.begin(); current_vertex != temporary_vertices.end(); ++current_vertex) {
          this->lcc.erase_dart(*current_vertex);
        } found_face_with_index = *current_face;
        break;
      }
    }
    
    // Brute force search
//    Dart_handle found_face_brute_force = LCC::null_handle;
//    for (typename LCC::Dart_range::iterator current_dart = this->lcc.darts().begin(), last_dart = this->lcc.darts().end(); current_dart != last_dart; ++current_dart) {
//      if (&*current_dart == &*free_vertices.front()) continue; // TODO: Compare by memory address for now
//      if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation<2>(this->lcc, current_dart, free_vertices.front()) || Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation<2>(this->lcc, current_dart, free_vertices.front())) {
//        CGAL_precondition(free_vertices.size() == temporary_vertices.size()); // TODO: If always true (which I think it is), this can be used as an optimisation
//        for (current_vertex = temporary_vertices.begin(); current_vertex != temporary_vertices.end(); ++current_vertex) {
//          this->lcc.erase_dart(*current_vertex);
//        } found_face_brute_force = current_dart;
//        break;
//      }
//    }
    
    // Validate that brute-force and index searches give equivalent results
//    if (found_face_brute_force == LCC::null_handle) {
//      CGAL_assertion(found_face_with_index == LCC::null_handle);
//    } else {
//      CGAL_assertion(found_face_with_index != LCC::null_handle);
//      bool equivalent = false;
//      for (typename LCC::template Dart_of_cell_range<2, 2>::iterator current_dart = this->lcc.template darts_of_cell<2, 2>(found_face_brute_force).begin(); current_dart != this->lcc.template darts_of_cell<2, 2>(found_face_brute_force).end(); ++current_dart) {
//        if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation<2>(this->lcc, found_face_with_index, current_dart) || Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation<2>(this->lcc, found_face_with_index, current_dart)) {
//          equivalent = true;
//          break;
//        }
//      } CGAL_assertion(equivalent);
//    }

    // Found, return existing
    if (found_face_with_index != LCC::null_handle) {
      
      // Validation
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::are_dart_ids_valid(this->lcc));
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::template are_attribute_ids_valid<2>(this->lcc));
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<0>(this->lcc);
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<2>(this->lcc);
      
      return found_face_with_index;
    }
    
    // Not found, create a new one
    else {
      typename LCC::template Attribute_handle<2>::type facet_attribute = this->lcc.template create_attribute<2>();
      this->lcc.template set_attribute<2>(free_vertices.front(), facet_attribute);
      index.template insert_cell<2>(this->lcc, free_vertices.front());
      int changed_darts = Call_if_lcc_with_ids<LCC>::make_temporary_darts_permanent(this->lcc);
    
      // Validation
      CGAL_postcondition(changed_darts == -1 || changed_darts == temporary_vertices.size());
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::are_dart_ids_valid(this->lcc));
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::template are_attribute_ids_valid<2>(this->lcc));
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<0>(this->lcc);
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<2>(this->lcc);
    
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
    
    // Validate index
//    index.template validate_index<dimension>(this->lcc);
    
    // Validate that all the facets are closed dimension-1 cells
    for (InputIterator current = begin; current != end; ++current) {
      for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart = this->lcc.template darts_of_cell<dimension-1, dimension-1>(*current).begin(); current_dart != this->lcc.template darts_of_cell<dimension-1, dimension-1>(*current).end(); ++current_dart) {
        for (unsigned int b = 0; b < dimension-1; ++b) {
          CGAL_precondition(!this->lcc.is_free(current_dart, b));
        }
      }
    }
    
    // Use free facets, make new (temporary) facets for those that aren't free
    std::list<Dart_handle> free_facets, temporary_facets;
    for (InputIterator current = begin; current != end; ++current) {
      bool all_free = true, all_not_free = true;
      for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart = this->lcc.template darts_of_cell<dimension-1, dimension-1>(*current).begin(); current_dart != this->lcc.template darts_of_cell<dimension-1, dimension-1>(*current).end(); ++current_dart) {
        if (this->lcc.is_free(current_dart, dimension-1)) all_not_free = false;
        else all_free = false;
      } CGAL_precondition(all_free != all_not_free);
      if (all_free) {
        //        std::cout << ".";
        free_facets.push_back(*current);
      } else {
        //        std::cout << "c";
        free_facets.push_back(Linear_cell_complex_cell_cloner<LCC>::template clone_cell_with_reversed_orientation<dimension-1>(this->lcc, *current));
        temporary_facets.push_back(free_facets.back());
      }
    }
    //    std::cout << std::endl;
    
    // Validate that all the facets and beta(dimension-1)-free
    for (typename std::list<Dart_handle>::iterator dart_in_current_facet = free_facets.begin(); dart_in_current_facet != free_facets.end(); ++dart_in_current_facet) {
      for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart_in_facet = this->lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_facet).begin(); current_dart_in_facet != this->lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_facet).end(); ++current_dart_in_facet) {
        CGAL_precondition(this->lcc.is_free(current_dart_in_facet, dimension-1));
      }
    }
    
    // Build an index on the ridges
    Index ridge_index;
    for (typename std::list<Dart_handle>::iterator dart_in_current_facet = free_facets.begin(); dart_in_current_facet != free_facets.end(); ++dart_in_current_facet) {
      for (typename LCC::template One_dart_per_incident_cell_range<dimension-2, dimension-1, dimension-1>::iterator dart_in_current_ridge = this->lcc.template one_dart_per_incident_cell<dimension-2, dimension-1, dimension-1>(*dart_in_current_facet).begin(); dart_in_current_ridge != this->lcc.template one_dart_per_incident_cell<dimension-2, dimension-1, dimension-1>(*dart_in_current_facet).end(); ++dart_in_current_ridge) {
        ridge_index.template insert_cell<dimension-2>(this->lcc, dart_in_current_ridge);
      }
    }
    
    // Validate that this will create a closed (topologically open) quasi-manifold cell
    for (typename std::list<Dart_handle>::iterator dart_in_current_facet_1 = free_facets.begin(); dart_in_current_facet_1 != free_facets.end(); ++dart_in_current_facet_1) {
      for (typename LCC::template One_dart_per_incident_cell_range<dimension-2, dimension-1, dimension-1>::iterator dart_in_current_ridge_1 = this->lcc.template one_dart_per_incident_cell<dimension-2, dimension-1, dimension-1>(*dart_in_current_facet_1).begin(); dart_in_current_ridge_1 != this->lcc.template one_dart_per_incident_cell<dimension-2, dimension-1, dimension-1>(*dart_in_current_facet_1).end(); ++dart_in_current_ridge_1) {
        int matches = 0;
        
        typename LCC::Point smallest_vertex;
        if (dimension == 3) {
          if (this->lcc.point(dart_in_current_ridge_1) < this->lcc.point(dart_in_current_ridge_1->beta(1))) smallest_vertex = this->lcc.point(dart_in_current_ridge_1);
          else smallest_vertex = this->lcc.point(dart_in_current_ridge_1->beta(1));
        } else {
          smallest_vertex = this->lcc.point(Index::template find_smallest_vertex<dimension-2>(this->lcc, dart_in_current_ridge_1));
        } typename Index::Results &results = ridge_index.template find_cells<dimension-2>(smallest_vertex);
        
        for (typename Index::Result_iterator dart_in_current_ridge_2 = results.begin(); dart_in_current_ridge_2 != results.end(); ++dart_in_current_ridge_2) {
          if (&*dart_in_current_ridge_1 == &**dart_in_current_ridge_2) continue; // TODO: Compare by memory address for now
          if (dimension == 3) {
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<1>(this->lcc, dart_in_current_ridge_1, *dart_in_current_ridge_2)) ++matches;
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<1>(this->lcc, dart_in_current_ridge_1, *dart_in_current_ridge_2)) ++matches;
          } else {
            for (typename LCC::template Dart_of_cell_range<dimension-2, dimension-2>::iterator current_dart_in_ridge_2 = this->lcc.template darts_of_cell<dimension-2, dimension-2>(*dart_in_current_ridge_2).begin(); current_dart_in_ridge_2 != this->lcc.template darts_of_cell<dimension-2, dimension-2>(*dart_in_current_ridge_2).end(); ++current_dart_in_ridge_2) {
              CGAL_assertion(this->lcc.template attribute<dimension-2>(*dart_in_current_ridge_2) == this->lcc.template attribute<dimension-2>(current_dart_in_ridge_2));
              if (&*dart_in_current_ridge_1 == &*current_dart_in_ridge_2) continue; // TODO: Compare by memory address for now
              if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<dimension-2>(this->lcc, dart_in_current_ridge_1, current_dart_in_ridge_2)) ++matches;
              if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<dimension-2>(this->lcc, dart_in_current_ridge_1, current_dart_in_ridge_2)) ++matches;
            }
          }
        }
        
        CGAL_precondition(matches == 1);
      }
    }
    
    // Create temporary structure
    for (typename std::list<Dart_handle>::iterator dart_in_current_facet_1 = free_facets.begin(); dart_in_current_facet_1 != free_facets.end(); ++dart_in_current_facet_1) {
      std::list<Dart_handle> ridges_in_current_facet_1;
      for (typename LCC::template One_dart_per_incident_cell_range<dimension-2, dimension-1, dimension-1>::iterator dart_in_current_ridge_1 = this->lcc.template one_dart_per_incident_cell<dimension-2, dimension-1, dimension-1>(*dart_in_current_facet_1).begin(); dart_in_current_ridge_1 != this->lcc.template one_dart_per_incident_cell<dimension-2, dimension-1, dimension-1>(*dart_in_current_facet_1).end(); ++dart_in_current_ridge_1) {
        ridges_in_current_facet_1.push_back(dart_in_current_ridge_1);
      } for (typename std::list<Dart_handle>::iterator dart_in_current_ridge_1 = ridges_in_current_facet_1.begin(); dart_in_current_ridge_1 != ridges_in_current_facet_1.end(); ++dart_in_current_ridge_1) {
        bool matched_ridge = !this->lcc.is_free(*dart_in_current_ridge_1, dimension-1);
        
        typename LCC::Point smallest_vertex;
        if (dimension == 3) {
          if (this->lcc.point(*dart_in_current_ridge_1) < this->lcc.point((*dart_in_current_ridge_1)->beta(1))) smallest_vertex = this->lcc.point(*dart_in_current_ridge_1);
          else smallest_vertex = this->lcc.point((*dart_in_current_ridge_1)->beta(1));
        } else {
          smallest_vertex = this->lcc.point(Index::template find_smallest_vertex<dimension-2>(this->lcc, *dart_in_current_ridge_1));
        } typename Index::Results &results = ridge_index.template find_cells<dimension-2>(smallest_vertex);
        
        for (typename Index::Result_iterator dart_in_current_ridge_2 = results.begin(); !matched_ridge && dart_in_current_ridge_2 != results.end(); ++dart_in_current_ridge_2) {
          if (&**dart_in_current_ridge_1 == &**dart_in_current_ridge_2) continue; // TODO: Compare by memory address for now
          
          if (dimension == 3) {
            
            // Match with same orientation
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<1>(this->lcc, *dart_in_current_ridge_1, *dart_in_current_ridge_2)) {
              matched_ridge = true;
              CGAL_assertion(this->lcc.vertex_attribute(*dart_in_current_ridge_1) == this->lcc.vertex_attribute(*dart_in_current_ridge_2));
              this->lcc.reverse_orientation_connected_component(*dart_in_current_ridge_2);
              index.update_index_after_reversing_orientation(this->lcc);
              ridge_index.update_index_after_reversing_orientation(this->lcc);
              CGAL_assertion(this->lcc.vertex_attribute(*dart_in_current_ridge_1) == this->lcc.vertex_attribute((*dart_in_current_ridge_2)->beta(1)));
              CGAL_assertion(this->lcc.vertex_attribute((*dart_in_current_ridge_1)->beta(1)) == this->lcc.vertex_attribute(*dart_in_current_ridge_2));
              this->lcc.template link_beta<2>(*dart_in_current_ridge_1, *dart_in_current_ridge_2);
            }
            
            // Match with opposite orientation
            else if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<1>(this->lcc, *dart_in_current_ridge_1, *dart_in_current_ridge_2)) {
              matched_ridge = true;
              CGAL_assertion(this->lcc.vertex_attribute(*dart_in_current_ridge_1) == this->lcc.vertex_attribute((*dart_in_current_ridge_2)->beta(1)));
              CGAL_assertion(this->lcc.vertex_attribute((*dart_in_current_ridge_1)->beta(1)) == this->lcc.vertex_attribute(*dart_in_current_ridge_2));
              this->lcc.template link_beta<2>(*dart_in_current_ridge_1, *dart_in_current_ridge_2);
            }
          }
          
          // General case for dimension > 3
          else {
            for (typename LCC::template Dart_of_cell_range<dimension-2, dimension-2>::iterator current_dart_in_ridge_2 = this->lcc.template darts_of_cell<dimension-2, dimension-2>(*dart_in_current_ridge_2).begin(); !matched_ridge && current_dart_in_ridge_2 != this->lcc.template darts_of_cell<dimension-2, dimension-2>(*dart_in_current_ridge_2).end(); ++current_dart_in_ridge_2) {
              CGAL_assertion(this->lcc.template attribute<dimension-2>(*dart_in_current_ridge_2) == this->lcc.template attribute<dimension-2>(current_dart_in_ridge_2));
              if (&**dart_in_current_ridge_1 == &*current_dart_in_ridge_2) continue;

              // Match with same orientation
              if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<dimension-2>(this->lcc, *dart_in_current_ridge_1, current_dart_in_ridge_2)) {
                matched_ridge = true;
                CGAL_assertion(this->lcc.vertex_attribute(*dart_in_current_ridge_1) == this->lcc.vertex_attribute(current_dart_in_ridge_2));
                this->lcc.reverse_orientation_connected_component(current_dart_in_ridge_2);
                index.update_index_after_reversing_orientation(this->lcc);
                ridge_index.update_index_after_reversing_orientation(this->lcc);
                CGAL_assertion(this->lcc.vertex_attribute(*dart_in_current_ridge_1) == this->lcc.vertex_attribute(current_dart_in_ridge_2->beta(1)));
                CGAL_assertion(this->lcc.vertex_attribute((*dart_in_current_ridge_1)->beta(1)) == this->lcc.vertex_attribute(current_dart_in_ridge_2));
                CGAL_assertion(this->lcc.template is_sewable<dimension-1>(*dart_in_current_ridge_1, current_dart_in_ridge_2));
                this->lcc.template sew<dimension-1>(*dart_in_current_ridge_1, current_dart_in_ridge_2);
              }
              
              // Match with opposite orientation
              else if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<dimension-2>(this->lcc, *dart_in_current_ridge_1, current_dart_in_ridge_2)) {
                matched_ridge = true;
                CGAL_assertion(this->lcc.vertex_attribute(*dart_in_current_ridge_1) == this->lcc.vertex_attribute(current_dart_in_ridge_2->beta(1)));
                CGAL_assertion(this->lcc.vertex_attribute((*dart_in_current_ridge_1)->beta(1)) == this->lcc.vertex_attribute(current_dart_in_ridge_2));
                CGAL_assertion(this->lcc.template is_sewable<dimension-1>(*dart_in_current_ridge_1, current_dart_in_ridge_2));
                this->lcc.template sew<dimension-1>(*dart_in_current_ridge_1, current_dart_in_ridge_2);
              }
            }
          }
        }
        
      }
    }
    
    // Validate the created temporary structure
    for (typename std::list<Dart_handle>::iterator dart_in_current_facet = free_facets.begin(); dart_in_current_facet != free_facets.end(); ++dart_in_current_facet) {
      for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart_in_facet = this->lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_facet).begin(); current_dart_in_facet != this->lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_facet).end(); ++current_dart_in_facet) {
        CGAL_assertion(!this->lcc.is_free(current_dart_in_facet, dimension-1));
        CGAL_assertion(current_dart_in_facet->beta(dimension-1) != current_dart_in_facet);
        CGAL_assertion(current_dart_in_facet->beta(dimension-1)->beta(dimension-1) == current_dart_in_facet);
      }
    }
    
    // Index search
    Dart_handle found_cell_with_index = LCC::null_handle;
    Dart_handle smallest_vertex = Index::template find_smallest_vertex<dimension>(this->lcc, free_facets.front());
    typename Index::Results &results = index.template find_cells<dimension>(this->lcc.point(smallest_vertex));
    for (typename Index::Result_iterator current_cell = results.begin(); found_cell_with_index == LCC::null_handle && current_cell != results.end(); ++current_cell) {
      CGAL_assertion(this->lcc.point(*current_cell) == this->lcc.point(smallest_vertex));
      for (typename LCC::template Dart_of_cell_range<0, dimension>::iterator current_dart = this->lcc.template darts_of_cell<0, dimension>(free_facets.front()).begin(); current_dart != this->lcc.template darts_of_cell<0, dimension>(free_facets.front()).end(); ++current_dart) {
        if (&*current_dart == &*free_facets.front()) continue; // TODO: Compare by memory address for now
        if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<dimension>(this->lcc, current_dart, free_facets.front()) || Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<dimension>(this->lcc, current_dart, free_facets.front())) {
          CGAL_assertion(free_facets.size() == temporary_facets.size()); // TODO: If always true (which I think it is), this can be used as an optimisation
          for (typename std::list<Dart_handle>::iterator current_facet = temporary_facets.begin(); current_facet != temporary_facets.end(); ++current_facet) {
            std::list<Dart_handle> darts_to_erase;
            for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart_in_facet = this->lcc.template darts_of_cell<dimension-1, dimension-1>(*current_facet).begin(); current_dart_in_facet != this->lcc.template darts_of_cell<dimension-1, dimension-1>(*current_facet).end(); ++current_dart_in_facet) {
              darts_to_erase.push_back(current_dart_in_facet);
            } for (typename std::list<Dart_handle>::iterator current_dart_to_erase = darts_to_erase.begin(); current_dart_to_erase != darts_to_erase.end(); ++current_dart_to_erase) {
              this->lcc.erase_dart(*current_dart_to_erase);
            }
          } found_cell_with_index = current_dart;
          break;
        }
      }
    }
    
    // Brute-force search
//    Dart_handle found_cell_brute_force = LCC::null_handle;
//    for (typename LCC::Dart_range::iterator current_dart = this->lcc.darts().begin(), last_dart = this->lcc.darts().end(); current_dart != last_dart; ++current_dart) {
//      if (&*current_dart == &*free_facets.front()) continue; // TODO: Compare by memory address for now
//      if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation<dimension>(this->lcc, current_dart, free_facets.front()) || Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation<dimension>(this->lcc, current_dart, free_facets.front())) {
//        CGAL_assertion(free_facets.size() == temporary_facets.size()); // TODO: If always true (which I think it is), this can be used as an optimisation
//        for (typename std::list<Dart_handle>::iterator current_facet = temporary_facets.begin(); current_facet != temporary_facets.end(); ++current_facet) {
//          std::list<Dart_handle> darts_to_erase;
//          for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart_in_facet = this->lcc.template darts_of_cell<dimension-1, dimension-1>(*current_facet).begin(); current_dart_in_facet != this->lcc.template darts_of_cell<dimension-1, dimension-1>(*current_facet).end(); ++current_dart_in_facet) {
//            darts_to_erase.push_back(current_dart_in_facet);
//          } for (typename std::list<Dart_handle>::iterator current_dart_to_erase = darts_to_erase.begin(); current_dart_to_erase != darts_to_erase.end(); ++current_dart_to_erase) {
//            this->lcc.erase_dart(*current_dart_to_erase);
//          }
//        } found_cell_brute_force = current_dart;
//        break;
//      }
//    }
    
    // Validate that brute-force and index searches give equivalent results
//    if (found_cell_brute_force == LCC::null_handle) {
//      CGAL_assertion(found_cell_with_index == LCC::null_handle);
//    } else {
//      CGAL_assertion(found_cell_with_index != LCC::null_handle);
//      bool equivalent = false;
//      for (typename LCC::template Dart_of_cell_range<dimension, dimension>::iterator current_dart = this->lcc.template darts_of_cell<dimension, dimension>(found_cell_brute_force).begin(); current_dart != this->lcc.template darts_of_cell<dimension, dimension>(found_cell_brute_force).end(); ++current_dart) {
//        if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation<dimension>(this->lcc, found_cell_with_index, current_dart) || Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation<dimension>(this->lcc, found_cell_with_index, current_dart)) {
//          equivalent = true;
//          break;
//        }
//      } CGAL_assertion(equivalent);
//    }
    
    // Found, return existing
    if (found_cell_with_index != LCC::null_handle) {
      
      // Validation
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::are_dart_ids_valid(this->lcc));
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::template are_attribute_ids_valid<dimension>(this->lcc));
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<0>(this->lcc);
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<dimension-1>(lcc);
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<dimension>(this->lcc);
      
      return found_cell_with_index;
    }
    
    // Not found, create a new one
    // TODO: When to create the 1-attributes?
    else {
      typename LCC::template Attribute_handle<dimension>::type cell_attribute = this->lcc.template create_attribute<dimension>();
      this->lcc.template set_attribute<dimension>(free_facets.front(), cell_attribute);
      index.template insert_cell<dimension>(this->lcc, free_facets.front());
      unsigned int temporary_darts = 0;
      for (typename std::list<Dart_handle>::iterator current_facet = temporary_facets.begin(); current_facet != temporary_facets.end(); ++current_facet) {
        for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart_in_facet = this->lcc.template darts_of_cell<dimension-1, dimension-1>(*current_facet).begin(); current_dart_in_facet != this->lcc.template darts_of_cell<dimension-1, dimension-1>(*current_facet).end(); ++current_dart_in_facet) {
          ++temporary_darts;
        }
      } int changed_darts = Call_if_lcc_with_ids<LCC>::make_temporary_darts_permanent(this->lcc);
      
      // Validation
      CGAL_postcondition(changed_darts == -1 || changed_darts == temporary_darts);
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::are_dart_ids_valid(this->lcc));
//      CGAL_postcondition(Linear_cell_complex_validator<LCC>::template are_attribute_ids_valid<dimension>(this->lcc));
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<0>(this->lcc);
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<dimension-1>(lcc);
//      Linear_cell_complex_validator<LCC>::template validate_map_attributes<dimension>(this->lcc);
      
      return free_facets.front();
    }
  }
  
  template <unsigned int dimension, class InputIterator>
  Dart_handle link_cells(InputIterator begin, InputIterator end) {
    // TODO: Doesn't use index yet
    CGAL_precondition(dimension >= 2);
    unsigned int number_of_facets = 0;
    for (InputIterator current = begin; current != end; ++current) {
      ++number_of_facets;
    } std::cout << "link_cells<" << dimension << ">(" << number_of_facets << " " << dimension << "-cells)" << std::endl;
    
    // Validate index
    index.template validate_index<dimension>(this->lcc);
    
    // Validate that all the facets are closed dimension-cells, unique, contain everything and are all beta(dimension)-free
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
    
    // Build an index on the ridges
    Index ridge_index;
    for (InputIterator dart_in_current_facet = begin; dart_in_current_facet != end; ++dart_in_current_facet) {
      for (typename LCC::template One_dart_per_incident_cell_range<dimension-1, dimension, dimension>::iterator dart_in_current_ridge = this->lcc.template one_dart_per_incident_cell<dimension-1, dimension, dimension>(*dart_in_current_facet).begin(); dart_in_current_ridge != this->lcc.template one_dart_per_incident_cell<dimension-1, dimension, dimension>(*dart_in_current_facet).end(); ++dart_in_current_ridge) {
        ridge_index.template insert_cell<dimension-1>(this->lcc, dart_in_current_ridge);
      }
    }
    
    // Validate ridge index
    ridge_index.template validate_index<dimension-1>(this->lcc);
    mark = this->lcc.get_new_mark();
    CGAL_precondition(this->lcc.is_whole_map_unmarked(mark));
    for (InputIterator dart_in_current_facet = begin; dart_in_current_facet != end; ++dart_in_current_facet) {
      for (typename LCC::template One_dart_per_incident_cell_range<dimension-1, dimension, dimension>::iterator dart_in_current_ridge = this->lcc.template one_dart_per_incident_cell<dimension-1, dimension, dimension>(*dart_in_current_facet).begin(); dart_in_current_ridge != this->lcc.template one_dart_per_incident_cell<dimension-1, dimension, dimension>(*dart_in_current_facet).end(); ++dart_in_current_ridge) {
        for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart = this->lcc.template darts_of_cell<dimension-1, dimension-1>(dart_in_current_ridge).begin(); current_dart != this->lcc.template darts_of_cell<dimension-1, dimension-1>(dart_in_current_ridge).end(); ++current_dart) {
          CGAL_precondition(!this->lcc.is_marked(current_dart, mark));
          this->lcc.mark(current_dart, mark);
        }
      }
    } this->lcc.negate_mark(mark);
    CGAL_postcondition(this->lcc.is_whole_map_unmarked(mark));
    this->lcc.free_mark(mark);
    
    // Validate that this will create a quasi-manifold cell complex
    for (InputIterator dart_in_current_facet_1 = begin; dart_in_current_facet_1 != end; ++dart_in_current_facet_1) {
      for (typename LCC::template One_dart_per_incident_cell_range<dimension-1, dimension, dimension>::iterator dart_in_current_ridge_1 = this->lcc.template one_dart_per_incident_cell<dimension-1, dimension, dimension>(*dart_in_current_facet_1).begin(); dart_in_current_ridge_1 != this->lcc.template one_dart_per_incident_cell<dimension-1, dimension, dimension>(*dart_in_current_facet_1).end(); ++dart_in_current_ridge_1) {
        int matches = 0;
        
        typename LCC::Point smallest_vertex;
        if (dimension == 2) {
          if (this->lcc.point(dart_in_current_ridge_1) < this->lcc.point(dart_in_current_ridge_1->beta(1))) smallest_vertex = this->lcc.point(dart_in_current_ridge_1);
          else smallest_vertex = this->lcc.point(dart_in_current_ridge_1->beta(1));
        } else {
          smallest_vertex = this->lcc.point(Index::template find_smallest_vertex<dimension-1>(this->lcc, dart_in_current_ridge_1));
        } typename Index::Results &results = ridge_index.template find_cells<dimension-1>(smallest_vertex);
        
        for (typename Index::Result_iterator dart_in_current_ridge_2 = results.begin(); dart_in_current_ridge_2 != results.end(); ++dart_in_current_ridge_2) {
          if (&*dart_in_current_ridge_1 == &**dart_in_current_ridge_2) continue; // TODO: Compare by memory address for now
          if (dimension == 2) {
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<1>(this->lcc,dart_in_current_ridge_1, *dart_in_current_ridge_2)) ++matches;
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<1>(this->lcc, dart_in_current_ridge_1, *dart_in_current_ridge_2)) ++matches;
          } else {
            for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart_in_ridge_2 = this->lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_ridge_2).begin(); current_dart_in_ridge_2 != this->lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_ridge_2).end(); ++current_dart_in_ridge_2) {
              CGAL_assertion(this->lcc.template attribute<dimension-1>(*dart_in_current_ridge_2) == this->lcc.template attribute<dimension-1>(current_dart_in_ridge_2));
              if (&*dart_in_current_ridge_1 == &*current_dart_in_ridge_2) continue; // TODO: Compare by memory address for now
              if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<dimension-1>(this->lcc, dart_in_current_ridge_1, current_dart_in_ridge_2)) ++matches;
              if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<dimension-1>(this->lcc, dart_in_current_ridge_1, current_dart_in_ridge_2)) ++matches;
            }
          }
        }
        
        CGAL_precondition(matches <= 1);
      }
    }
    
    // Create structure
    for (InputIterator dart_in_current_facet_1 = begin; dart_in_current_facet_1 != end; ++dart_in_current_facet_1) {
      std::list<Dart_handle> ridges_in_current_facet_1;
      for (typename LCC::template One_dart_per_incident_cell_range<dimension-1, dimension, dimension>::iterator dart_in_current_ridge_1 = this->lcc.template one_dart_per_incident_cell<dimension-1, dimension, dimension>(*dart_in_current_facet_1).begin(); dart_in_current_ridge_1 != this->lcc.template one_dart_per_incident_cell<dimension-1, dimension, dimension>(*dart_in_current_facet_1).end(); ++dart_in_current_ridge_1) {
        ridges_in_current_facet_1.push_back(dart_in_current_ridge_1);
      } for (typename std::list<Dart_handle>::iterator dart_in_current_ridge_1 = ridges_in_current_facet_1.begin(); dart_in_current_ridge_1 != ridges_in_current_facet_1.end(); ++dart_in_current_ridge_1) {
        bool matched_ridge = !this->lcc.is_free(*dart_in_current_ridge_1, dimension);
        
        typename LCC::Point smallest_vertex;
        if (dimension == 2) {
          if (this->lcc.point(*dart_in_current_ridge_1) < this->lcc.point((*dart_in_current_ridge_1)->beta(1))) smallest_vertex = this->lcc.point(*dart_in_current_ridge_1);
          else smallest_vertex = this->lcc.point((*dart_in_current_ridge_1)->beta(1));
        } else {
          smallest_vertex = this->lcc.point(Index::template find_smallest_vertex<dimension-1>(this->lcc, *dart_in_current_ridge_1));
        } typename Index::Results &results = ridge_index.template find_cells<dimension-1>(smallest_vertex);
        
        for (typename Index::Result_iterator dart_in_current_ridge_2 = results.begin(); !matched_ridge && dart_in_current_ridge_2 != results.end(); ++dart_in_current_ridge_2) {
          if (&**dart_in_current_ridge_1 == &**dart_in_current_ridge_2) continue; // TODO: Compare by memory address for now
          
          if (dimension == 2) {
            
            // Match with same orientation
            if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<1>(this->lcc, *dart_in_current_ridge_1, *dart_in_current_ridge_2)) {
              matched_ridge = true;
              CGAL_assertion(this->lcc.vertex_attribute(*dart_in_current_ridge_1) == this->lcc.vertex_attribute(*dart_in_current_ridge_2));
              this->lcc.reverse_orientation_connected_component(*dart_in_current_ridge_2);
              index.update_index_after_reversing_orientation(this->lcc);
              ridge_index.update_index_after_reversing_orientation(this->lcc);
              CGAL_assertion(this->lcc.vertex_attribute(*dart_in_current_ridge_1) == this->lcc.vertex_attribute((*dart_in_current_ridge_2)->beta(1)));
              CGAL_assertion(this->lcc.vertex_attribute((*dart_in_current_ridge_1)->beta(1)) == this->lcc.vertex_attribute(*dart_in_current_ridge_2));
              this->lcc.template link_beta<2>(*dart_in_current_ridge_1, *dart_in_current_ridge_2);
            }
            
            // Match with opposite orientation
            else if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<1>(this->lcc, *dart_in_current_ridge_1, *dart_in_current_ridge_2)) {
              matched_ridge = true;
              CGAL_assertion(this->lcc.vertex_attribute(*dart_in_current_ridge_1) == this->lcc.vertex_attribute((*dart_in_current_ridge_2)->beta(1)));
              CGAL_assertion(this->lcc.vertex_attribute((*dart_in_current_ridge_1)->beta(1)) == this->lcc.vertex_attribute(*dart_in_current_ridge_2));
              this->lcc.template link_beta<2>(*dart_in_current_ridge_1, *dart_in_current_ridge_2);
            }
          }
          
          // General case for dimension > 2
          else {
            for (typename LCC::template Dart_of_cell_range<dimension-1, dimension-1>::iterator current_dart_in_ridge_2 = this->lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_ridge_2).begin(); !matched_ridge && current_dart_in_ridge_2 != this->lcc.template darts_of_cell<dimension-1, dimension-1>(*dart_in_current_ridge_2).end(); ++current_dart_in_ridge_2) {
              CGAL_assertion(this->lcc.template attribute<dimension-1>(*dart_in_current_ridge_2) == this->lcc.template attribute<dimension-1>(current_dart_in_ridge_2));
              if (&**dart_in_current_ridge_1 == &*current_dart_in_ridge_2) continue; // TODO: Compare by memory address for now
              
              // Match with same orientation
              if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_orientation_from<dimension-1>(this->lcc, *dart_in_current_ridge_1, current_dart_in_ridge_2)) {
                matched_ridge = true;
                CGAL_assertion(this->lcc.vertex_attribute(*dart_in_current_ridge_1) == this->lcc.vertex_attribute(current_dart_in_ridge_2));
                this->lcc.reverse_orientation_connected_component(current_dart_in_ridge_2);
                index.update_index_after_reversing_orientation(this->lcc);
                ridge_index.update_index_after_reversing_orientation(this->lcc);
                CGAL_assertion(this->lcc.vertex_attribute(*dart_in_current_ridge_1) == this->lcc.vertex_attribute(current_dart_in_ridge_2->beta(1)));
                CGAL_assertion(this->lcc.vertex_attribute((*dart_in_current_ridge_1)->beta(1)) == this->lcc.vertex_attribute(current_dart_in_ridge_2));
                CGAL_assertion(this->lcc.template is_sewable<dimension>(*dart_in_current_ridge_1, current_dart_in_ridge_2));
                this->lcc.template sew<dimension>(*dart_in_current_ridge_1, current_dart_in_ridge_2);
              }
              
              // Match with opposite orientation
              else if (Linear_cell_complex_cell_comparator<LCC>::template have_same_geometry_and_opposite_orientation_from<dimension-1>(this->lcc, *dart_in_current_ridge_1, current_dart_in_ridge_2)) {
                matched_ridge = true;
                CGAL_assertion(this->lcc.vertex_attribute(*dart_in_current_ridge_1) == this->lcc.vertex_attribute(current_dart_in_ridge_2->beta(1)));
                CGAL_assertion(this->lcc.vertex_attribute((*dart_in_current_ridge_1)->beta(1)) == this->lcc.vertex_attribute(current_dart_in_ridge_2));
                CGAL_assertion(this->lcc.template is_sewable<dimension>(*dart_in_current_ridge_1, current_dart_in_ridge_2));
                this->lcc.template sew<dimension>(*dart_in_current_ridge_1, current_dart_in_ridge_2);
              }
            }
          }
        }
        
      }
    }
    
    // Validate the created temporary structure
    for (InputIterator dart_in_current_facet = begin; dart_in_current_facet != end; ++dart_in_current_facet) {
      for (typename LCC::template Dart_of_cell_range<dimension, dimension>::iterator current_dart_in_facet = this->lcc.template darts_of_cell<dimension, dimension>(*dart_in_current_facet).begin(); current_dart_in_facet != this->lcc.template darts_of_cell<dimension, dimension>(*dart_in_current_facet).end(); ++current_dart_in_facet) {
        if (this->lcc.is_free(current_dart_in_facet, dimension)) continue;
        CGAL_assertion(current_dart_in_facet->beta(dimension) != current_dart_in_facet);
        CGAL_assertion(current_dart_in_facet->beta(dimension)->beta(dimension) == current_dart_in_facet);
      }
    }
    
    return *begin;
  }
};

#endif
