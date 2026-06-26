/**************************************************************************
***
*** Copyright (c) 2026 IonQ, Inc.
*** Copyright (c) 1995-2000 Regents of the University of California,
***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
*** Copyright (c) 2000-2012 Regents of the University of Michigan,
***               Saurabh N. Adya, Jarrod A. Roy, David A. Papa and
***               Igor L. Markov
***
***  Contact author(s): abk@cs.ucsd.edu, igor.markov1@gmail.com
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

// Placement transform driver.
// It builds a grid placement, rotates it, and writes the original and
// transformed placements so the transform helpers can be checked.

// Created by Igor Markov, July 1997, VLSICAD ABKGroup UCLA/CS
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <cmath>
#include <fstream>

#include "placement.h"
#include "transf.h"

using std::cout;
using std::endl;
using std::ofstream;

static inline double normalizedCoord(double value) {
  return (std::fabs(value) < 1e-12) ? 0.0 : value;
}

int main() {
  Timer tm1;
  Placement pl(Placement::_Grid1, 50);
  RectRotation rectRotate(Pi / 3, pl.getBBox());
  tm1.stop();
  cout << "Creating : " << tm1 << endl;
  tm1.start();
  ofstream ofs("orig.gen.pla");
  ofs << pl;
  tm1.stop();
  cout << "Saving   : " << tm1 << endl;
  tm1.start();
  rectRotate(pl);
  tm1.stop();
  cout << "Rotating : " << tm1 << endl;
  tm1.start();
  ofstream ofs1("rotd.gen.pla");
  ofs1 << pl;
  tm1.stop();
  cout << "Saving   : " << tm1 << endl;
  BBox bbox = pl.getBBox();
  bbox.xMin = normalizedCoord(bbox.xMin);
  bbox.yMin = normalizedCoord(bbox.yMin);
  bbox.xMax = normalizedCoord(bbox.xMax);
  bbox.yMax = normalizedCoord(bbox.yMax);
  cout << "Resulting BBox : " << bbox << endl;
  return 0;
}
