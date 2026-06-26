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

// created by Stefanus Mantik on 06/06/98
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "hbboxWRSMTPlEval.h"

#include "expectedBBox.h"
#include "placeEvalTables.h"

void HBBoxWRSMTPlEval::_recomputeCost() {
  BBox topLeft, bottomRight;
  unsigned index;

  for (itHGFEdgeLocal e = _hg.edgesBegin(); e != _hg.edgesEnd(); e++) {
    topLeft.clear();
    bottomRight.clear();
    index = (*e)->getIndex();
    unsigned degree = (*e)->getDegree();
    for (itHGFNodeLocal n = (*e)->nodesBegin(); n != (*e)->nodesEnd(); n++) {
      unsigned regId = (*_assignment)[(*n)->getIndex()];
      abkassert(regId < _regions->size(), "Assignment index too big");
      if (topLeft.isEmpty() && bottomRight.isEmpty())
        topLeft = bottomRight = (*_regions)[regId];
      else {
        topLeft = HBBoxUpdateTopLeft(topLeft, (*_regions)[regId]);
        bottomRight = HBBoxUpdateBottomRight(bottomRight, (*_regions)[regId]);
      }
    }
    BBox expBBox((topLeft.xMin + topLeft.xMax) / 2,
                 (bottomRight.yMin + bottomRight.yMax) / 2,
                 (bottomRight.xMin + bottomRight.xMax) / 2,
                 (topLeft.yMin + topLeft.yMax) / 2);
    if (degree > MaxRSMTTableDegree) degree = MaxRSMTTableDegree;
    abkfatal(index < _netCost.size(), "index out of range");
    _netCost[index] = lookupTable(expBBox, degree);
  }
}
