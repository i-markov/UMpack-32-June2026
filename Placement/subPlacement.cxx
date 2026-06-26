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

//  Created:  Igor Markov, VLSI CAD ABKGROUP UCLA, June 15, 1997.
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <ABKCommon/abkassert.h>
#include <ABKCommon/infolines.h>
#include <Geoms/plGeom.h>
#include <Placement/placementUtils.h>
#include <Placement/subPlacement.h>

#include <iostream>

using std::endl;
using std::ostream;

namespace {

struct SubPlacementPointAccessor {
  const SubPlacement& pl;
  SubPlacementPointAccessor(const SubPlacement& p) : pl(p) {}
  const Point& operator()(unsigned i) const { return pl[i]; }
};

}  // namespace

Point SubPlacement::getCenterOfGravity() const {
  return placement_utils::centerOfGravity(getSize(),
                                          SubPlacementPointAccessor(*this));
}

Point SubPlacement::getCenterOfGravity(double* weights) const {
  return placement_utils::weightedCenterOfGravity(
      getSize(), SubPlacementPointAccessor(*this), weights);
}

double SubPlacement::getPolygonArea() const {
  return placement_utils::polygonArea(getSize(),
                                      SubPlacementPointAccessor(*this));
}

bool SubPlacement::isInsidePolygon(const Point& p) const {
  return placement_utils::isInsidePolygon(getSize(),
                                          SubPlacementPointAccessor(*this), p);
}

ostream& operator<<(ostream& out, const SubPlacement& arg) {
  TimeStamp tm;
  out << tm;
  out << " Total points in sub placement : " << arg.getSize() << endl;
  for (unsigned i = 0; i < arg.getSize(); i++) out << arg.points(i) << endl;
  return out;
}

#ifndef AUTHMAR
#define AUTHMAR
static const char AuthorMarkov[50] =
    "062897, Igor Markov, VLSI CAD ABKGROUP UCLA";
#endif

double operator-(const SubPlacement& arg1, const SubPlacement& arg2)
/* returns the Linf distance between placements
   which you can compare then to epsilon     */
{
  abkfatal(arg1.getSize() == arg2.getSize(),
           " Comparing placements of different size");
  struct LeftAccessor {
    const SubPlacement& pl;
    LeftAccessor(const SubPlacement& p) : pl(p) {}
    const Point& operator()(unsigned i) const { return pl[i]; }
  };
  return placement_utils::linfDistance(arg1.getSize(), LeftAccessor(arg1),
                                       LeftAccessor(arg2));
}
