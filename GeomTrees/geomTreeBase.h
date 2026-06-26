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

//! CONTACTS = "ABK, Stefanus, Igor"

#ifndef _GEOM_TREE_BASE_H
#define _GEOM_TREE_BASE_H

#include <utility>

#include "Geoms/bbox.h"
#include "Geoms/point.h"

//: Represent geometric tree edge by two points
class GeomTreeEdge : public std::pair<Point, Point> {
 public:
  GeomTreeEdge(const Point& f, const Point& s);
  GeomTreeEdge() : std::pair<Point, Point>() {}

  BBox getBBox() const;
  double getMDist() const;
};

//: Node of geometric tree
class GeomTreeNode {
 public:
  int parent, id;
  double weight;
  GeomTreeEdge edge;

  GeomTreeNode() : parent(-1), id(-1), weight(0), edge() {}
};

inline double GeomTreeEdge::getMDist() const { return mDist(first, second); }

#endif
