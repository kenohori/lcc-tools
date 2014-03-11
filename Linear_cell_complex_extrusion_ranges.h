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

#ifndef LCC_EXTRUSION_RANGES_H
#define LCC_EXTRUSION_RANGES_H

#include <set>

template <class LCC>
struct Print_extrusion_range;

// Abstracts a range, possibly including intermediate key points
template <class LCC>
struct Extrusion_range {
private:
  struct Extrusion_range_comparator {
    bool operator()(const typename LCC::FT &ft1, const typename LCC::FT &ft2) const {
      return ft1 < ft2;
    }
  };
  
public:
  typedef std::set<typename LCC::FT, Extrusion_range_comparator> Range;
  
  template <class LCC_>
  friend struct Print_extrusion_range;
  
private:
  Range range;
  
public:
  Extrusion_range() {
    //    std::cout << "Extrusion_range() -> Extrusion_range<" << this << ">";
    //    Print_extrusion_range<LCC>::print(*this);
    //    std::cout << std::endl;
  }
  
  Extrusion_range(const Extrusion_range &er) {
    //    std::cout << "Extrusion_range(Extrusion_range) -> Extrusion_range<" << this << ">" << std::endl;
    for (typename Range::iterator current_element = er.begin(); current_element != er.end(); ++current_element) {
      range.insert(*current_element);
    } /*std::cout << "\t->";
       Print_extrusion_range<LCC>::print(*this);
       std::cout << std::endl;*/
  }
  
  Extrusion_range(typename LCC::FT &min, typename LCC::FT &max) {
    //    std::cout << "Extrusion_range(" << min << ", " << max << ") -> Extrusion_range<" << this << ">" << std::endl;
    range.insert(min);
    range.insert(max);
  }
  
  Extrusion_range &operator=(const Extrusion_range &er) {
    //    std::cout << "Extrusion_range = Extrusion_range" << std::endl;
    if (this == &er) return *this;
    for (typename Range::iterator current_element = er.begin(); current_element != er.end(); ++current_element) {
      range.insert(*current_element);
    } return *this;
  }
  
  typename LCC::FT min() const {
    return *range.begin();
  }
  
  typename LCC::FT max() const {
    return *range.rbegin();
  }
  
  typename Range::const_iterator begin() const {
    return range.begin();
  }
  
  typename Range::const_iterator end() const {
    return range.end();
  }
  
  void insert_elements(const Extrusion_range &er) {
    //    std::cout << "Extrusion_range<" << this << ">";
    //    Print_extrusion_range<LCC>::print(*this);
    //    std::cout << ".insert_elements(Extrusion_range<" << &er << ">";
    //    Print_extrusion_range<LCC>::print(er);
    //    std::cout << ")" << std::endl;
    range.insert(er.begin(), er.end());
    //    std::cout << "\t->";
    //    Print_extrusion_range<LCC>::print(*this);
    //    std::cout << std::endl;
  }
  
  bool overlaps_or_touches(const Extrusion_range &other) const {
    if (other.min() >= min() && other.min() <= max()) return true;
    if (other.max() >= min() && other.max() <= max()) return true;
    if (other.min() <= min() && other.max() >= max()) return true;
    return false;
  }
  
  bool covers(typename LCC::FT &rmin, typename LCC::FT &rmax) const {
    if (rmin >= min() && rmax <= max()) return true;
    return false;
  }
  
  bool is_in(typename LCC::FT &value) const {
    if (range.count(value) > 0) return true;
    return false;
  }
  
  bool is_in(typename LCC::FT &rmin, typename LCC::FT &rmax) const {
    for (typename Range::const_iterator current_value = begin(); current_value != end(); ++current_value) {
      if (*current_value == rmin) {
        typename Range::const_iterator next_value = current_value;
        ++next_value;
        if (next_value != end()) {
          if (*next_value == rmax) return true;
        }
      }
    } return false;
  }
  
  bool is_surrounded(typename LCC::FT &value) const {
    if (range.count(value) > 0 && min() != value && max() != value) return true;
    return false;
  }
};

template <class LCC>
struct Print_extrusion_range {
  static void print(const Extrusion_range<LCC> &er) {
    std::cout << "[";
    for (typename Extrusion_range<LCC>::Range::const_iterator current_element = er.range.begin(); current_element != er.range.end(); ++current_element) {
      std::cout << *current_element << " ";
    } std::cout << "]";
  }
};

template <class LCC>
struct Print_extrusion_ranges;

// Abstracts a set of ranges (belonging to a single object)
template <class LCC>
struct Extrusion_ranges {
public:
  typedef Extrusion_ranges<LCC> Self;
  typedef Extrusion_range<LCC> Range;
  
private:
  // Caution: this only compared the first element since sets are assumed to be disjoint!
  struct Extrusion_ranges_comparator {
    bool operator()(const Range &cell1, const Range &cell2) {
      return cell1.min() < cell2.min();
    }
  };
  
public:
  typedef std::set<Range, Extrusion_ranges_comparator> Range_set;
  
  template <class LCC_>
  friend struct Print_extrusion_ranges;
  
protected:
  Range_set ranges;
  
public:
  Extrusion_ranges() {
    //    std::cout << "Extrusion_ranges() -> Extrusion_ranges<" << this << ">";
    //    Print_extrusion_ranges<LCC>::print(*this);
    //    std::cout << std::endl;
  }
  
  void add_range(typename LCC::FT min, typename LCC::FT max) {
    //    std::cout << "Extrusion_ranges<" << this << ">.add_range(" << min << ", " << max << ")" << std::endl;
    Range new_range(min, max);
    add_range(new_range);
  }
  
  void add_range(const Range &r) {
    //    std::cout << "Extrusion_ranges<" << this << ">.add_range(const Extrusion_range<" << &r << ">)" << std::endl;
    Range new_range(r);
    add_range(new_range);
  }
  
  void add_range(Range &r) {
    //    std::cout << "Extrusion_ranges<" << this << ">.add_range(Extrusion_range<" << &r << ">";
    //    Print_extrusion_range<LCC>::print(r);
    //    std::cout << ")" << std::endl;
    
    // Check if the range overlaps or touches any other
    std::list<typename Range_set::iterator> overlapping_ranges;
    typename Range_set::iterator current_range = ranges.begin();
    while (current_range != ranges.end()) {
      if (current_range->overlaps_or_touches(r)) overlapping_ranges.push_back(current_range);
      ++current_range;
    }
    
    // If there are overlaps, merge them into new_range and remove them
    if (!overlapping_ranges.empty()) {
      //      std::cout << "\toverlapping ranges!" << std::endl;
      for (typename std::list<typename Range_set::iterator>::const_iterator current_overlapping_range = overlapping_ranges.begin();
           current_overlapping_range != overlapping_ranges.end();
           ++current_overlapping_range) {
        r.insert_elements(**current_overlapping_range);
      } /*std::cout << "\tranges before: ";
         Print_extrusion_ranges<LCC>::print(*this);
         std::cout << std::endl;*/
      overlapping_ranges.push_back(overlapping_ranges.back()); // Copy the last iterator (otherwise the code breaks when there's only one)
      ++overlapping_ranges.back(); //Advance this iterator to be one position after the last overlapping range
      ranges.erase(overlapping_ranges.front(), overlapping_ranges.back());
      //      std::cout << "\tranges after: ";
      //      Print_extrusion_ranges<LCC>::print(*this);
      //      std::cout << std::endl;
    }
    
    //    std::cout << "\t r = ";
    //    Print_extrusion_range<LCC>::print(r);
    //    std::cout << std::endl;
    
    ranges.insert(r);   // Note that this makes a copy of the range!
    //    std::cout << "\t->";
    //    Print_extrusion_ranges<LCC>::print(*this);
    //    std::cout << std::endl;
  }
  
  void add_ranges(Extrusion_ranges<LCC> &er) {
    //    std::cout << "Extrusion_ranges<" << this << ">.add_ranges(Extrusion_ranges)" << std::endl;
    for (typename Range_set::iterator current_range = er.ranges.begin(); current_range != er.ranges.end(); ++current_range) {
      add_range(*current_range);
    }
  }
  
  void add_ranges(const Extrusion_ranges<LCC> &er) {
    //    std::cout << "Extrusion_ranges<" << this << ">.add_ranges(const Extrusion_ranges)" << std::endl;
    for (typename Range_set::const_iterator current_range = er.ranges.begin(); current_range != er.ranges.end(); ++current_range) {
      add_range(*current_range);
    }
  }
  
  typename LCC::FT min() const {
    return ranges.begin()->min();
  }
  
  typename LCC::FT max() const {
    return ranges.rbegin()->max();
  }
  
  typename Range_set::const_iterator begin() const {
    return ranges.begin();
  }
  
  typename Range_set::const_iterator end() const {
    return ranges.end();
  }
  
  bool covers(typename LCC::FT &rmin, typename LCC::FT &rmax) const {
    for (typename Range_set::const_iterator current_range = begin(); current_range != end(); ++current_range) {
      if (current_range->covers(rmin, rmax)) return true;
    } return false;
  }
  
  bool is_in(typename LCC::FT &value) const {
    for (typename Range_set::const_iterator current_range = begin(); current_range != end(); ++current_range) {
      if (current_range->is_in(value)) return true;
    } return false;
  }
  
  bool is_in(typename LCC::FT &rmin, typename LCC::FT &rmax) const {
    for (typename Range_set::const_iterator current_range = begin(); current_range != end(); ++current_range) {
      if (current_range->is_in(rmin, rmax)) return true;
    } return false;
  }
  
  bool is_surrounded(typename LCC::FT &value) const {
    for (typename Range_set::const_iterator current_range = begin(); current_range != end(); ++current_range) {
      if (current_range->is_surrounded(value)) return true;
    } return false;
  }
};

template <class LCC>
struct Print_extrusion_ranges {
  static void print(const Extrusion_ranges<LCC> &er) {
    std::cout << "[";
    for (typename Extrusion_ranges<LCC>::Range_set::const_iterator current_range = er.ranges.begin(); current_range != er.ranges.end(); ++current_range) {
      Print_extrusion_range<LCC>::print(*current_range);
      std::cout << " ";
    } std::cout << "]";
  }
};

// Abstracts a map of cell->ranges for a particular dimension
template <class LCC, unsigned int dimension>
struct Extrusion_ranges_map_of_dimension {
public:
  typedef std::map<typename LCC::template Attribute_const_handle<dimension>::type, Extrusion_ranges<LCC> > type;
  
  type ranges_map;
};

// Abstracts a tuple of maps of cell->ranges, each element containing the map of a particular dimension
template <class LCC, unsigned int dimension = LCC::dimension, class Result = CGAL::cpp11::tuple<> >
struct Extrusion_ranges_tuple_per_dimension_up_to;

template <class LCC, class ... Result>
struct Extrusion_ranges_tuple_per_dimension_up_to<LCC, 0, CGAL::cpp11::tuple<Result ...> > {
  typedef CGAL::cpp11::tuple<Extrusion_ranges_map_of_dimension<LCC, 0>, Result ...> type;
};

template <class LCC, unsigned int dimension, class ... Result>
struct Extrusion_ranges_tuple_per_dimension_up_to<LCC, dimension, CGAL::cpp11::tuple<Result ...> > {
  typedef typename Extrusion_ranges_tuple_per_dimension_up_to<LCC, dimension-1, CGAL::cpp11::tuple<Extrusion_ranges_map_of_dimension<LCC, dimension>, Result ...> >::type type;
};

// Abstracts a tuple of maps of cell->ranges, each element containing the map of a particular dimension
template <class LCC_>
struct Extrusion_ranges_tuple_per_dimension {
public:
  typedef LCC_ LCC;
  typedef typename Extrusion_ranges_tuple_per_dimension_up_to<LCC>::type type;
  
  type ranges;
  
  Extrusion_ranges_tuple_per_dimension() {
    std::cout << "Extrusion_ranges_tuple_per_dimension(dim = " << LCC::dimension << ")" << std::endl;
  }
};

template <class ERTPD, unsigned int dimension>
struct Print_extrusion_ranges_of_dimension {
  static void print(ERTPD &ertpd) {
    std::cout << dimension << "-cells (" << std::get<dimension>(ertpd.ranges).ranges_map.size() << "):" << std::endl;
    for (typename Extrusion_ranges_map_of_dimension<typename ERTPD::LCC, dimension>::type::const_iterator current_cell_range = std::get<dimension>(ertpd.ranges).ranges_map.begin();
         current_cell_range != std::get<dimension>(ertpd.ranges).ranges_map.end();
         ++current_cell_range) {
      std::cout << "\t" << dimension << "-cell<" << &*current_cell_range->first << "> -> ";
      Print_extrusion_ranges<typename ERTPD::LCC>::print(current_cell_range->second);
      std::cout << std::endl;
    }
  }
};

template <class ERTPD, unsigned int dimension>
struct Print_extrusion_ranges_of_dimension_and_lower {
  static void print(ERTPD &ertpd) {
    Print_extrusion_ranges_of_dimension_and_lower<ERTPD, dimension-1>::print(ertpd);
    Print_extrusion_ranges_of_dimension<ERTPD, dimension>::print(ertpd);
  }
};

template <class ERTPD>
struct Print_extrusion_ranges_of_dimension_and_lower<ERTPD, 0> {
  static void print(ERTPD &ertpd) {
    Print_extrusion_ranges_of_dimension<ERTPD, 0>::print(ertpd);
  }
};

template <class ERTPD>
struct Print_all_extrusion_ranges {
  static void print(ERTPD &ertpd) {
    Print_extrusion_ranges_of_dimension_and_lower<ERTPD, ERTPD::LCC::dimension>::print(ertpd);
  }
};

template <class LCC>
struct Extrusion_ranges_reader {
public:
  typedef Extrusion_ranges_tuple_per_dimension<LCC> Extrusion_ranges;
  
  static void load_from_ranges_file(const char *path, LCC &lcc, Extrusion_ranges &extrusion_ranges, std::set<typename LCC::FT> &all_ranges) {
    
    // Make a map with the ids of the cells of the highest dimension
    std::map<int, typename LCC::template Attribute_const_handle<LCC::dimension>::type> cells_map;
    typename LCC::template Attribute_const_range<LCC::dimension>::type &cells = lcc.template attributes<LCC::dimension>();
    for (typename LCC::template Attribute_const_range<LCC::dimension>::type::const_iterator current_cell = cells.begin();
         current_cell != cells.end();
         ++current_cell) {
      cells_map[current_cell->info()] = current_cell;
      //std::cout << "mapping cell<" << &*current_cell << ">[" << current_cell->info() << "]" << std::endl;
    }
    
    
    std::cout << "Opening file " << path << "..." << std::endl;
    std::ifstream input_file_stream;
    input_file_stream.open(path);
    if (!input_file_stream.is_open()) {
      std::cerr << "Could not open file." << std::endl;
      return;
    }
    
    std::string current_line, current_word;
    int cell_id;
    typename LCC::FT cell_min, cell_max;
    std::istringstream line_stream;
    
    while (!input_file_stream.eof()) {
      getline(input_file_stream, current_line);
      std::cout << "Line: \'" << current_line << "\'" << std::endl;
      
      // Empty line
      if (current_line.size() == 0) {
        continue;
      }
      
      // Just a newline sequence
      if (current_line[0] == '\r' || current_line[0] == '\n') {
        continue;
      }
      
      line_stream = std::istringstream(current_line + "\n");
      
      line_stream >> current_word;
      cell_id = string_to_any(current_word, cell_id);
      line_stream >> current_word;
      cell_min = string_to_any(current_word, cell_min);
      line_stream >> current_word;
      cell_max = string_to_any(current_word, cell_max);
      
      if (cells_map.count(cell_id) <= 0) {
        std::cerr << "\tNo matching cell found for id=" << cell_id << std::endl;
        continue;
      }
      
      // Assign the ranges to the highest dimensional cells
      std::cout << "\t" << cell_id << " -> cell<" << &*cells_map[cell_id] << ">" << std::endl;
      std::get<LCC::dimension>(extrusion_ranges.ranges).ranges_map[cells_map[cell_id]].add_range(cell_min, cell_max);
      all_ranges.insert(cell_min);
      all_ranges.insert(cell_max);
    }
  }
  
  static void load_from_ogr(const char *path, LCC &lcc, Extrusion_ranges &extrusion_ranges, std::set<typename LCC::FT> &all_ranges, const char *min_field_name, const char *max_field_name) {
    
    // Make a map with the ids of the cells of the highest dimension
    std::map<int, typename LCC::template Attribute_const_handle<LCC::dimension>::type> cells_map;
    typename LCC::template Attribute_const_range<LCC::dimension>::type &cells = lcc.template attributes<LCC::dimension>();
    for (typename LCC::template Attribute_const_range<LCC::dimension>::type::const_iterator current_cell = cells.begin();
         current_cell != cells.end();
         ++current_cell) {
      cells_map[current_cell->info()] = current_cell;
      //std::cout << "mapping cell<" << &*current_cell << ">[" << current_cell->info() << "]" << std::endl;
    }
    
    std::cout << "Opening file " << path << "..." << std::endl;
    
    OGRRegisterAll();
    OGRDataSource *data_source = OGRSFDriverRegistrar::Open(path, false);
    if (data_source == NULL) {
      std::cerr << "Could not open file." << std::endl;
      return;
    }
    
    std::cout << "\tType: " << data_source->GetDriver()->GetName() << std::endl;
    std::cout << "\tLayers: " << data_source->GetLayerCount() << std::endl;
    
    for (int current_layer = 0; current_layer < data_source->GetLayerCount(); ++current_layer) {
      OGRLayer *layer = data_source->GetLayer(current_layer);
      layer->ResetReading();
      std::cout << "\tReading layer " << current_layer << " (" << layer->GetFeatureCount(true) << " features)..." << std::endl;
      
      OGRFeature *feature;
      while ((feature = layer->GetNextFeature()) != NULL) {
        switch(feature->GetGeometryRef()->getGeometryType()) {
            
          case wkbPolygon: {
            
            if (max_field_name == NULL) continue;
            
            int rangemin_index = -1;
            if (min_field_name != NULL) rangemin_index = feature->GetFieldIndex(min_field_name);
            int rangemax_index = feature->GetFieldIndex(max_field_name);
            if (rangemax_index == -1) continue;
            
            long cell_id = feature->GetFID();
            typename LCC::FT cell_min = 0, cell_max;
            if (rangemin_index != -1) cell_min = feature->GetFieldAsDouble(rangemin_index);
            cell_max = feature->GetFieldAsDouble(rangemax_index);
            std::cout << "\t" << cell_id << " -> cell<" << &*cells_map[cell_id] << ">" << std::endl;
            std::get<LCC::dimension>(extrusion_ranges.ranges).ranges_map[cells_map[cell_id]].add_range(cell_min, cell_max);
            all_ranges.insert(cell_min);
            all_ranges.insert(cell_max);
            
            break;
          }
            
          default:
            std::cerr << "Feature type not supported." << std::endl;
            break;
        }
        
        OGRFeature::DestroyFeature(feature);
      }
    }
    
    OGRDataSource::DestroyDataSource(data_source);
  }
};

#endif
