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

// Created 990435 by Andrew Caldwell

#ifndef _RBPL_PINLOC_H_
#define _RBPL_PINLOC_H_

#include "ABKCommon/abkassert.h"
#include <string>
#include <vector>
#include "HGraphWDims/hgWDims.h"
#include "Placement/placeWOri.h"

class PinLocCalculator {
  const HGraphWDimensions& _hgraph;
  const PlacementWOrient _placement;

 public:
  PinLocCalculator(const HGraphWDimensions& hgraph,
                   const std::vector<Orientation>& orient, const Placement& pl)
      : _hgraph(hgraph), _placement(pl, orient) {}

  PinLocCalculator(const HGraphWDimensions& hgraph, const PlacementWOrient& pl)
      : _hgraph(hgraph), _placement(pl) {}

  Point getPinOffset(unsigned cellId, unsigned netId) const;
  Point getPinOffset(unsigned cellId, unsigned netId, Orientation ori) const;
};

inline Point PinLocCalculator::getPinOffset(unsigned cellId,
                                            unsigned netId) const {
  return getPinOffset(cellId, netId, _placement.getOrient(cellId));
}

inline Point PinLocCalculator::getPinOffset(unsigned cellId, unsigned netId,
                                            Orientation ori) const {
  const HGFNode& node = _hgraph.getNodeByIdx(cellId);

  bool found = false;
  itHGFEdgeLocal e;
  for (e = node.edgesBegin(); e != node.edgesEnd(); e++)
    if ((*e)->getIndex() == netId) {
      found = true;
      break;
    }

  abkassert(found, "asked for pin offset for an edge not on the node");

  unsigned edgeOffset = e - node.edgesBegin();

  // pinOffset is the offset if the cell has orientation 'N'
  Point pinOffset = _hgraph.edgesOnNodePinOffset(edgeOffset, cellId);

  double h = _hgraph.getNodeHeight(cellId);
  double w = _hgraph.getNodeWidth(cellId);

  if (!ori.isFaceUp()) pinOffset.x = w - pinOffset.x;

  double tmp;
  switch (ori.getAngle()) {
    case 0:  //( x,  y)
      break;
    case 90:  //( y, -x)
      tmp = pinOffset.x;
      pinOffset.x = pinOffset.y;
      pinOffset.y = w - tmp;
      break;
    case 180:  //(-x, -y)
      pinOffset.x = w - pinOffset.x;
      pinOffset.y = h - pinOffset.y;
      break;
    case 270:  //(-y,  x)
      tmp = pinOffset.y;
      pinOffset.y = pinOffset.x;
      pinOffset.x = h - tmp;
      break;
    default:
      abkfatal(0, " Orientation corrupted ");
  }

  // abkassert(pinOffset.x >= 0 && pinOffset.x <= max(w,h),
  //		"pinOffset.x is out of range");
  // abkassert(pinOffset.y >= 0 && pinOffset.y <= max(w,h),
  //		"pinOffset.y is out of range");

  return pinOffset;
}

#endif
