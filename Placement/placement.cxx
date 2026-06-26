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

//  971206 ilm  Placement::_Grid2 now constructs a regular grid size*size
//                filling the given BBox
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <ABKCommon/abkassert.h>
#include <ABKCommon/abkio.h>
#include <ABKCommon/infolines.h>
#include <Geoms/plGeom.h>
#include <Placement/placement.h>
#include <Placement/placementUtils.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <strings.h>

using std::endl;
using std::istream;
using std::ostream;
using std::string;

namespace {

struct PlacementPointAccessor {
  const Placement& pl;
  PlacementPointAccessor(const Placement& p) : pl(p) {}
  const Point& operator()(unsigned i) const { return pl[i]; }
};

}  // namespace

Placement::Placement(unsigned n, Point pt) : nPts(n) {
  abkfatal(n > 0, "Empty placement");
  points = new Point[n];
  for (unsigned k = 0; k < n; k++) points[k] = pt;
}

Point Placement::getCenterOfGravity() const {
  return placement_utils::centerOfGravity(nPts, PlacementPointAccessor(*this));
}

Point Placement::getCenterOfGravity(double* weights) const {
  return placement_utils::weightedCenterOfGravity(
      nPts, PlacementPointAccessor(*this), weights);
}

double Placement::getPolygonArea() const {
  return placement_utils::polygonArea(nPts, PlacementPointAccessor(*this));
}

bool Placement::isInsidePolygon(const Point& p) const {
  return placement_utils::isInsidePolygon(nPts, PlacementPointAccessor(*this),
                                          p);
}

Placement::Placement(const Mapping& pullBackMap, const Placement& from)
    : nPts(pullBackMap.getSourceSize())
// pullback
{
  if (nPts == 0) {
    points = NULL;
    return;
  }
  points = new Point[nPts];
  for (unsigned k = 0; k < nPts; k++) points[k] = from[pullBackMap[k]];
}

Placement::Placement(Placement::Examples ex, unsigned size, BBox bb) {
  unsigned i, j, nXFixed, nYFixed, nXandYFixed, num = 0;

  // setting sizes;
  switch (ex) {
    case _Random:
      nPts = size;
      if (size / 10 > 20)
        nYFixed = 20;
      else
        nYFixed = size / 10;
      nXandYFixed = nXFixed = nYFixed;
      break;
    case _Grid1:
      if (size == 0) size = 5;
      nPts = size * size + 4 * size;
      nXandYFixed = nXFixed = nYFixed = 4 * size;
      break;
    case _Grid2:
      if (size == 0) size = 5;
      nPts = size * size;
      nXandYFixed = 0;
      break;
  }

  switch (ex) {
    case _Random: {
      RandomPoint pt(bb);
      points = new Point[nPts];
      for (i = 0; i < size; i++) points[i] = pt;
      for (i = 0; i < nYFixed / 4; i++) points[i].x = 0;
      for (i = nYFixed / 4; i < nYFixed / 2; i++) points[i].y = 1;
      for (i = nYFixed / 2; i < 3 * nYFixed / 4; i++) points[i].x = 1;
      for (i = 3 * nYFixed / 4; i < nYFixed; i++) points[i].y = 0;
      return;
    }
    case _Grid1: {
      points = new Point[nPts];
      double xGridStep = (bb.xMax - bb.xMin) / (size + 1);
      double yGridStep = (bb.yMax - bb.yMin) / (size + 1);
      for (i = 0; i < size; i++)
        points[num++] = Point(bb.xMin, bb.yMin + yGridStep * (1 + i));
      for (i = 0; i < size; i++)
        points[num++] = Point(bb.xMax, bb.yMin + yGridStep * (1 + i));
      for (i = 0; i < size; i++)
        points[num++] = Point(bb.xMin + xGridStep * (1 + i), bb.yMin);
      for (i = 0; i < size; i++)
        points[num++] = Point(bb.xMin + xGridStep * (1 + i), bb.yMax);
      for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
          points[num++] = Point(bb.xMin + xGridStep * (1 + i),
                                bb.yMin + yGridStep * (1 + j));
      return;
    }
    case _Grid2: {
      points = new Point[nPts];
      double xGridStep = (bb.xMax - bb.xMin) / (size - 1);
      double yGridStep = (bb.yMax - bb.yMin) / (size - 1);
      for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
          points[num++] =
              Point(bb.xMin + xGridStep * i, bb.yMin + yGridStep * j);
      return;
    }
    default:
      abkfatal(0, "Unknown placement example");
  }
}

Placement::Placement(istream& in) {
  int lineNo = 1;
  string word1, word2, word3;
  in >> eathash(lineNo) >> word1 >> word2 >> word3;
  abkfatal(strcasecmp(word1.c_str(), "Total") == 0,
           "Placement files should start with \'Total points :\'\n");
  abkfatal(strcasecmp(word2.c_str(), "points") == 0,
           "Placements files should start with \'Total points :\'\n");
  abkfatal(strcasecmp(word3.c_str(), ":") == 0,
           "Placements files should start with \'Total points :\'\n");
  in >> my_isnumber(lineNo) >> nPts >> skiptoeol;

  abkfatal(nPts > 0, "Empty placement");
  points = new Point[nPts];
  for (unsigned i = 0; i < nPts; i++)
    in >> eathash(lineNo) >> points[i] >> skiptoeol;
}

Placement& Placement::operator=(
    const Placement& from) {  // this can be sped up with strcpy()
  if (from.points == points) return *this;
  delete[] points;
  nPts = from.getSize();
  points = new Point[nPts];
  for (unsigned i = 0; i < getSize(); i++) points[i] = from[i];
  return *this;
}

#ifndef AUTHMAR
#define AUTHMAR
static const char AuthorMarkov[50] =
    "062897, Igor Markov, VLSI CAD ABKGROUP UCLA";
#endif

double operator-(const Placement& arg1, const Placement& arg2)
/* returns the Linf distance between placements
   which you can compare then to epsilon     */
{
  abkfatal(arg1.getSize() == arg2.getSize(),
           " Comparing placements of different size");
  struct LeftAccessor {
    const Placement& pl;
    LeftAccessor(const Placement& p) : pl(p) {}
    const Point& operator()(unsigned i) const { return pl[i]; }
  };
  return placement_utils::linfDistance(arg1.getSize(), LeftAccessor(arg1),
                                       LeftAccessor(arg2));
}

ostream& operator<<(ostream& out, const Placement& arg) {
  TimeStamp tm;
  out << tm;
  out << "Total points : " << arg.nPts << endl;
  for (unsigned i = 0; i < arg.nPts; i++) out << arg.points[i] << " " << "\n";
  return out;
}

void Placement::flipXY() {
  for (unsigned i = getSize(); i != 0;) {
    i--;
    std::swap(points[i].x, points[i].y);
  }
}
