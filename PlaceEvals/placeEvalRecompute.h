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

#ifndef __PLACE_EVAL_RECOMPUTE_H__
#define __PLACE_EVAL_RECOMPUTE_H__

#include <vector>

#include "ABKCommon/abkassert.h"
#include "HGraph/hgFixed.h"
#include "Placement/layoutBBoxes.h"
#include "Placement/placement.h"
#include "expectedBBox.h"
#include "placeEvalTables.h"

namespace PlaceEvalRecompute {

template <bool UseEdgeWeight>
inline void bboxCost(const HGraphFixed& hg, const Placement& placement,
                     std::vector<double>& netCost) {
  BBox netBBox;
  unsigned index;

  for (itHGFEdgeLocal e = hg.edgesBegin(); e != hg.edgesEnd(); e++) {
    netBBox.clear();
    index = (*e)->getIndex();
    for (itHGFNodeLocal n = (*e)->nodesBegin(); n != (*e)->nodesEnd(); n++)
      netBBox += placement[(*n)->getIndex()];
    abkfatal(index < netCost.size(), "index out of range");
    netCost[index] = netBBox.getHalfPerim();
    if (UseEdgeWeight) netCost[index] *= (*e)->getWeight();
  }
}

template <bool ZeroSameRegion, bool ShrinkToCenter, bool UseChengScale>
inline void hbboxCost(const HGraphFixed& hg, const LayoutBBoxes& regions,
                      const std::vector<unsigned>& assignment,
                      std::vector<double>& netCost) {
  BBox topLeft, bottomRight, saveBBox, currentBBox;
  unsigned index;

  for (itHGFEdgeLocal e = hg.edgesBegin(); e != hg.edgesEnd(); e++) {
    topLeft.clear();
    bottomRight.clear();
    saveBBox.clear();
    bool sameBBox = true;
    index = (*e)->getIndex();
    abkfatal(index < netCost.size(), "index out of range");

    for (itHGFNodeLocal n = (*e)->nodesBegin(); n != (*e)->nodesEnd(); n++) {
      unsigned regId = assignment[(*n)->getIndex()];
      abkassert(regId < regions.size(), "Assignment index too big");
      currentBBox = regions[regId];
      if (ShrinkToCenter) {
        HBBoxShrinkX2Center(currentBBox, 0.16666667);
        HBBoxShrinkY2Center(currentBBox, 0.16666667);
      }

      if (topLeft.isEmpty() && bottomRight.isEmpty())
        saveBBox = topLeft = bottomRight = currentBBox;
      else {
        topLeft = HBBoxUpdateTopLeft(topLeft, currentBBox);
        bottomRight = HBBoxUpdateBottomRight(bottomRight, currentBBox);
        if (ZeroSameRegion && sameBBox &&
            (saveBBox.xMin != currentBBox.xMin ||
             saveBBox.xMax != currentBBox.xMax ||
             saveBBox.yMin != currentBBox.yMin ||
             saveBBox.yMax != currentBBox.yMax))
          sameBBox = false;
      }
    }

    if (ZeroSameRegion && sameBBox) {
      netCost[index] = 0;
      continue;
    }

    BBox expBBox((topLeft.xMin + topLeft.xMax) / 2,
                 (bottomRight.yMin + bottomRight.yMax) / 2,
                 (bottomRight.xMin + bottomRight.xMax) / 2,
                 (topLeft.yMin + topLeft.yMax) / 2);
    netCost[index] = expBBox.getHalfPerim();
    if (UseChengScale) {
      unsigned degree = (*e)->getDegree();
      if (degree > MaxChengTableDegree) degree = MaxChengTableDegree;
      netCost[index] *= ChengTable[degree];
    }
  }
}

}  // namespace PlaceEvalRecompute

#endif
