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

#include "Linear_cell_complex_extruder_with_range.h"

int main(int argc, const char * argv[]) {
  Linear_cell_complex_reader_writer<2> lccrw;
  lccrw.load_from_ogr("tud_merged.shp");
  
  Linear_cell_complex_extruder_with_range<2> lcce2;
  lcce2.load_ranges_from_ogr("tud_merged.shp", lccrw.lcc, NULL, "ex_height");
  Linear_cell_complex_extruder<2>::Higher_dimensional_cell_complex extruded = lcce2.extrude_using_ranges(lccrw.lcc);
  
  Linear_cell_complex_extruder_with_range<3> lcce3;
  lcce3.load_ranges_from_ogr("tud_merged.shp", extruded, "ex_built", "ex_destroy");
  Linear_cell_complex_extruder<3>::Higher_dimensional_cell_complex reextruded = lcce3.extrude_using_ranges(extruded);
  Cell_complex_printer<Linear_cell_complex_extruder<3>::Higher_dimensional_cell_complex>::print(reextruded);
  
  return 0;
}

