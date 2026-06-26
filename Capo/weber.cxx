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

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <algorithm>

#include "Ctainers/bitBoard.h"
#include "capoBin.h"
#include "capoPlacer.h"

using std::sort;
using std::vector;

// Created by Igor Markov on April 7 1999

Point CapoPlacer::getWeberLocation(const CapoBin& bin) {
  Point center = bin.getCenter(), closest(-DBL_MAX, -DBL_MAX);
  BBox wb = getWeberRegion(bin);

  if (center.x < wb.xMin)
    closest.x = wb.xMin;
  else if (center.x <= wb.xMax)
    closest.x = center.x;
  else
    closest.x = wb.xMax;

  if (center.y < wb.yMin)
    closest.y = wb.yMin;
  else if (center.y <= wb.yMax)
    closest.y = center.y;
  else
    closest.y = wb.yMax;

  return closest;
}

Point CapoPlacer::getWeberLocation(unsigned nodeIdx) {
  Point node = _placement[nodeIdx], closest(-DBL_MAX, -DBL_MAX);
  BBox wb = getWeberRegion(nodeIdx);

  if (node.x < wb.xMin)
    closest.x = wb.xMin;
  else if (node.x <= wb.xMax)
    closest.x = node.x;
  else
    closest.x = wb.xMax;

  if (node.y < wb.yMin)
    closest.y = wb.yMin;
  else if (node.y <= wb.yMax)
    closest.y = node.y;
  else
    closest.y = wb.yMax;

  return closest;
}

BBox CapoPlacer::getWeberRegion(const CapoBin& bin) {
  const HGraphWDimensions& hgraph = getNetlistHGraph();

  const vector<unsigned>& nodeIds = bin.getCellIds();
  bit_vector nodes(hgraph.getNumNodes(), false);  // in the bin
  BitBoard netBoard(hgraph.getNumEdges());        // insident to bin

  unsigned i;
  for (i = 0; i != nodeIds.size(); i++) {
    unsigned nodeIdx = nodeIds[i];
    nodes[nodeIdx] = true;
    const HGFNode& node = hgraph.getNodeByIdx(nodeIdx);

    for (itHGFEdgeLocal edgeIt = node.edgesBegin(); edgeIt != node.edgesEnd();
         edgeIt++)
      netBoard.setBit((*edgeIt)->getIndex());
  }

  const auto& nets = netBoard.getIndicesOfSetBits();

  // now find the Weber region
  vector<double> xEnds, yEnds;
  xEnds.reserve(nodeIds.size()), yEnds.reserve(nodeIds.size());
  for (auto netIt = nets.begin(); netIt != nets.end(); netIt++) {
    const HGFEdge& edge = hgraph.getEdgeByIdx(*netIt);
    BBox box;  // created empty
    bool hasFixed = false;
    for (itHGFNodeLocal nodeIt = edge.nodesBegin(); nodeIt != edge.nodesEnd();
         nodeIt++) {
      unsigned nodeIdx = (*nodeIt)->getIndex();
      if (!nodes[nodeIdx]) {
        hasFixed = true;
        box += _placement[nodeIdx];
      }
    }
    if (hasFixed) {
      xEnds.push_back(box.xMin);
      xEnds.push_back(box.xMax);
      yEnds.push_back(box.yMin);
      yEnds.push_back(box.yMax);
    }
  }
  if (xEnds.size()) {
    sort(xEnds.begin(), xEnds.end());
    sort(yEnds.begin(), yEnds.end());
    unsigned n = xEnds.size() / 2;
    double xLeft = xEnds[n - 1], xRight = xEnds[n];
    double yLeft = yEnds[n - 1], yRight = yEnds[n];
    return BBox(xLeft, yLeft, xRight, yRight);
  } else
    return bin.getBBox();
}

BBox CapoPlacer::getWeberRegion(unsigned nodeIdx) {
  const HGraphWDimensions& hgraph = getNetlistHGraph();
  bit_vector visited(hgraph.getNumEdges(), false);  // insident to bin

  const HGFNode& node = hgraph.getNodeByIdx(nodeIdx);

  vector<double> xEnds, yEnds;
  xEnds.reserve(2 * node.getDegree()), yEnds.reserve(2 * node.getDegree());
  for (itHGFEdgeLocal edgeIt = node.edgesBegin(); edgeIt != node.edgesEnd();
       edgeIt++) {
    BBox box;  // created empty
    const HGFEdge& edge = (**edgeIt);
    unsigned netIdx = edge.getIndex();
    if (!visited[netIdx]) {
      visited[netIdx] = true;
      bool hasFixed = false;
      for (itHGFNodeLocal nodeIt = edge.nodesBegin(); nodeIt != edge.nodesEnd();
           nodeIt++) {
        unsigned nodeNum = (*nodeIt)->getIndex();
        if (nodeNum != nodeIdx) {
          hasFixed = true;
          box += _placement[nodeNum];
        }
      }
      if (hasFixed) {
        xEnds.push_back(box.xMin);
        xEnds.push_back(box.xMax);
        yEnds.push_back(box.yMin);
        yEnds.push_back(box.yMax);
      }
    }
  }

  if (xEnds.empty()) return _placement[nodeIdx];

  sort(xEnds.begin(), xEnds.end());
  sort(yEnds.begin(), yEnds.end());
  unsigned n = xEnds.size() / 2;
  double xLeft = xEnds[n - 1], xRight = xEnds[n];
  double yLeft = yEnds[n - 1], yRight = yEnds[n];
  return BBox(xLeft, yLeft, xRight, yRight);
}

unsigned CapoPlacer::nodesToMove(const CapoBin& bin) {
  const vector<unsigned>& nodes = bin.getCellIds();
  unsigned count = 0;
  for (unsigned k = 0; k != nodes.size(); k++) {
    unsigned idx = nodes[k];
    BBox weberBBox = getWeberRegion(idx);
    //    if (!weberBBox.intersects(bin.getBBox())) count++;
    if (!bin.getBBox().contains(weberBBox)) count++;
  }
  return count;
}
