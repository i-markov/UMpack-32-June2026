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

#ifndef _MSTBASE_H_
#define _MSTBASE_H_

#include <GeomTrees/geomTreeBase.h>
#include <Geoms/point.h>

#include <vector>

//: Base class for minimum spanning tree
class MSTBase {
 protected:
  std::vector<GeomTreeNode> _mstNodes;
  double _mstCost;

  void print(std::ostream& out) const;

 public:
  MSTBase(const std::vector<Point>& points)
      : _mstNodes(points.size()), _mstCost(-1) {
    for (unsigned i = 0; i < points.size(); i++) {
      _mstNodes[i].id = i;
      _mstNodes[i].edge.first.x = _mstNodes[i].edge.second.x = points[i].x;
      _mstNodes[i].edge.first.y = _mstNodes[i].edge.second.y = points[i].y;
      //_mstNodes[i].edge.first = _mstNodes[i].edge.second = points[i];
    }
  }
  ~MSTBase() { _mstNodes.clear(); }

  double getCost() const { return _mstCost; }
  const std::vector<GeomTreeNode>& getTree() const { return _mstNodes; }
  std::vector<GeomTreeNode>& getTree() { return _mstNodes; }

  friend std::ostream& operator<<(std::ostream& out, const MSTBase& mst) {
    mst.print(out);
    return out;
  }
};

inline void MSTBase::print(std::ostream& out) const {
  out << "Root: 0 (weight = " << _mstNodes[0].weight << ")\n"
      << "  loc: " << _mstNodes[0].edge.first << std::endl;
  for (unsigned i = 1; i < _mstNodes.size(); i++) {
    out << "Node: " << i << " (weight = " << _mstNodes[i].weight
        << ", parent = " << _mstNodes[i].parent << ")\n"
        << "  loc: " << _mstNodes[i].edge.first << std::endl
        << "  dist = " << _mstNodes[i].edge.getMDist() << std::endl;
  }
}

#endif
