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

#ifndef _PLACEMENT_UTILS_H_
#define _PLACEMENT_UTILS_H_

#include <Geoms/plGeom.h>

#include <algorithm>

namespace placement_utils {

template <class PointGetter>
inline Point centerOfGravity(unsigned nPts, PointGetter get) {
  Point cog(0, 0);
  for (unsigned i = 0; i < nPts; ++i) {
    const Point& pt = get(i);
    cog.x += pt.x;
    cog.y += pt.y;
  }
  if (nPts) {
    cog.x /= nPts;
    cog.y /= nPts;
  }
  return cog;
}

template <class PointGetter>
inline Point weightedCenterOfGravity(unsigned nPts, PointGetter get,
                                     const double* weights) {
  Point cog(0, 0);
  for (unsigned i = 0; i < nPts; ++i) {
    const Point& pt = get(i);
    cog.x += weights[i] * pt.x;
    cog.y += weights[i] * pt.y;
  }
  if (nPts) {
    cog.x /= nPts;
    cog.y /= nPts;
  }
  return cog;
}

template <class PointGetter>
inline double polygonArea(unsigned nPts, PointGetter get) {
  if (nPts < 2) return 0.0;

  double area = 0.0;
  Point prev = get(0);

  for (unsigned i = 1; i < nPts; ++i) {
    const Point& curr = get(i);
    area += prev.x * curr.y - curr.x * prev.y;
    prev = curr;
  }

  const Point& first = get(0);
  area += prev.x * first.y - first.x * prev.y;
  return 0.5 * area;
}

template <class PointGetter>
inline bool isInsidePolygon(unsigned nPts, PointGetter get, const Point& p) {
  if (nPts < 3) return false;

  unsigned counter = 0;
  Point p1 = get(0);
  Point p2;
  double xinters;

  for (unsigned i = 0; i < nPts; ++i) {
    p2 = get(i % nPts);
    if ((p.y > std::min(p1.y, p2.y)) && (p.y <= std::max(p1.y, p2.y)) &&
        (p.x <= std::max(p1.x, p2.x)) && (p1.y != p2.y)) {
      xinters = (p.y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y) + p1.x;
      if (p1.x == p2.x || p.x <= xinters) ++counter;
    }
    p1 = p2;
  }

  return (counter % 2) != 0;
}

template <class LeftGetter, class RightGetter>
inline double linfDistance(unsigned nPts, LeftGetter left, RightGetter right) {
  double maxDst = 0.0;
  for (unsigned i = 0; i < nPts; ++i) {
    const double dst = mDist(left(i), right(i));
    if (maxDst < dst) maxDst = dst;
  }
  return maxDst;
}

template <class LeftGetter, class RightGetter>
inline double averageDistance(unsigned nPts, LeftGetter left,
                              RightGetter right) {
  if (!nPts) return 0.0;

  double sum = 0.0;
  for (unsigned i = 0; i < nPts; ++i) sum += mDist(left(i), right(i));
  return sum / nPts;
}

}  // namespace placement_utils

#endif
