/**************************************************************************
***
*** Copyright (c) 2026 IonQ, Inc.
*** Copyright (c) 1995-2000 Regents of the University of California,
***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
*** Copyright (c) 2000-2007 Regents of the University of Michigan,
***               Saurabh N. Adya, Jarrod A. Roy, David A. Papa and
***               Igor L. Markov
***
***  Contact author(s): abk@cs.ucsd.edu, imarkov@umich.edu
***  Original Affiliation:   UCLA, Computer Science Department,
***                          Los Angeles, CA 90095-1596 USA
***
***  Permission is hereby granted, free of charge, to any person obtaining
***  a copy of this software and associated documentation files (the
***  "Software"), to deal in the Software without restriction, including
***  without limitation
***  the rights to use, copy, modify, merge, publish, distribute, sublicense,
***  and/or sell copies of the Software, and to permit persons to whom the
***  Software is furnished to do so, subject to the following conditions:
***
***  The above copyright notice and this permission notice shall be included
***  in all copies or substantial portions of the Software.
***
*** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
*** EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
*** OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
*** IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
*** CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
*** OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
*** THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***
***
***************************************************************************/

#ifndef _ISPD06DENSITYMAP_H_
#define _ISPD06DENSITYMAP_H_

// The purpose of this utility is to provide a C++ alternative to the perl
// script associated with the ISPD 2006 placement competition.  The script
// computes the "density" of a design by dividing the core area into a regular
// grid and computing the average overflow of the grid cells.

// Documentation for the perl script:
// #  Density target checking perl script: check_density_target.pl
//     * Usage: check_density_target.pl node-file solution-pl-file scl-file
//     density-target
//     * example: %check_density_target.pl adaptec1.nodes adaptec1-yourplacer.pl
//     adaptec1.scl 0.5
//     * Simply superimpose a bin grid (120 x 120, i.e., 10 circuit row heights)
//     over placement solution to calculate bin utilization. The bin utilization
//     is defined as the movable cell usage divided by the free space (i.e., bin
//     area - fixed usage in bin). The script reports the number of violation
//     bins, the average overflow etc.
//     * Pay extra attention to "Scaled Overflow per bin" value in the last
//     line. This is the one we're going to use for score.

// More description of its function from email:
// Probably you want to ask what is the "scaled overflow per bin" value.

// scaled_overflow_per_bin = some_constant * total_overflow_amount / [
// total_movable_cell_area / (single_bin_area * density_target)]

//"total_movable_cell_area / (single_bin_area * density_target)" gives us the
//approximate number of bins to pack all the movable cells, and
//"scaled_overflow_per_bin" is simply the total overflow divided by this number.
//The reason we're using "total_movable_cell_area / (single_bin_area *
//density_target)" instead of total number of bins is to make
//"scaled_overflow_per_bin" independent of design density, size of placement
//area, size of bin grid etc.

// here is some of the code:
//$bin_musage = the total area of movables
//$target_density = the required density
//$bin_free_space = this is the space available for movable objects,
//                   which is the area of the rectangle minus the area of fixed
//                   objects
//$total_overflow_amount += ($bin_musage - $target_density * $bin_free_space);

// printf "\tOverflow per bin: %f\tTotal overflow amount: %f\n",
// $total_overflow_amount/$total_bin_count, $total_overflow_amount;
//$scaled_overflow_per_bin =
//($total_overflow_amount*$XUNIT*$YUNIT*$target_density)/($total_movable_area *
//400); printf "\tScaled Overflow per bin: %f\n", $scaled_overflow_per_bin *
// $scaled_overflow_per_bin;

#include <string>

#include "Geoms/ibbox.h"
#include "RBPlace/RBPlacement.h"
#include "baseGeneric2DResourceMap.h"

class ISPD06DensityMap
    : public Generic2DResourceMap<double, DoubleDistT, DoubleZeroT> {
 public:
  ISPD06DensityMap(const RBPlacement& rbplace, bool makeMacrosFixed);

  void compute(void);

  void setTargetDensity(double targetDensity);

  // inspectors
  Palette getPalette(void) const;
  double setTargetDensity(void) const { return _targetDensity; }

  void plotXPM(const std::string& baseFileName) const;
  void plotXPM(const std::string& baseFileName, double maxResourceVal) const;

  void printInfo(std::ostream& os) const;
  void printHistogramInfo(std::ostream& os) const;

  double getScaledOverflow(void) const;
  double getTotalOverfill(void) const;
  double getMaxOverfill(void) const;
  unsigned getNumOverfullTiles(void) const;

  void removeCellUsage(unsigned cellId);
  void addCellUsage(unsigned cellId);
  double estimateOverfullnessChangeFromCellMove(unsigned cellId,
                                                const Point& newPos) const;

  BBox getGridCell(const Point& p) const;

 private:  // functions
  // takes a BBox box from the layout region and returns the IBBox that contains
  // box in tile coordinates
  IBBox getTileCoords(const BBox& box) const;

 private:  // data
  const RBPlacement& _rbplace;
  bool _makeMacrosFixed;
  double _targetDensity;

  double _total_overflow_amount;
  unsigned _numTiles;
  double _total_movable_area;
};

#endif
