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

// created by Stefanus Mantik on 10/13/98
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "hbboxNwChengPlEval.h"

#include <cmath>

#include "expectedBBox.h"
#include "placeEvalTables.h"

using std::vector;

void HBBoxNwChengPlEval::_computeCorrection() {
  _cutCorrection = std::vector<double>(_regions->size() + 1, 0);
  _uncutCorrection = std::vector<double>(_regions->size() + 1, 0);

  for (itHGFNodeGlobal n = _hg.nodesBegin(); n != _hg.nodesEnd(); n++) {
    _cutCorrection[_getRegionId(**n)] += 1;
    _uncutCorrection[_getRegionId(**n)] += 1;
  }

  _cutCorrection[_regions->size()] = 0;  // this is for cell that is not
                                         // in any block
  _uncutCorrection[_regions->size()] = 0;

  for (unsigned i = 0; i < _regions->size(); i++) {
    if (_cutCorrection[i] <= 1) _cutCorrection[i] = 1;
    if (_uncutCorrection[i] <= 1) _uncutCorrection[i] = 1;

    _cutCorrection[i] = 1.0 / pow(_cutCorrection[i], 0.10);
    _uncutCorrection[i] = 1.0 / pow(_uncutCorrection[i], 0.38);
    abkfatal(_cutCorrection[i] <= 1, "correction factor is too large");
  }
}

void HBBoxNwChengPlEval::_recomputeCost() {
  for (itHGFEdgeLocal e = _hg.edgesBegin(); e != _hg.edgesEnd(); e++)
    _netCost[(*e)->getIndex()] = _computeOneEdge(**e);
}

double HBBoxNwChengPlEval::_computeOneEdge(const HGFEdge& edge) {
  unsigned degree = edge.getDegree();
  if (degree < 2) return 0.0;

  Point centerOfMass = _getCenterOfMass(edge);
  BBox topLeft, bottomRight;

  BBox actualBBox;
  unsigned edgeId = edge.getIndex();

  _edgeIsCut[edgeId] = _isEdgeCut(edge);

  for (itHGFNodeLocal n = edge.nodesBegin(); n != edge.nodesEnd(); n++) {
    unsigned regId = _getRegionId(**n);
    BBox currentBBox = _getNodesBBox(**n);

    actualBBox += (*_placement)[(*n)->getIndex()];

    if (regId != _regions->size())  // a block, not just a point.
    {
      if (_edgeIsCut[edgeId])
        _scaleBBox(currentBBox, centerOfMass, _cutCorrection[regId] * _betaCut);
      else
        _scaleBBox(currentBBox, centerOfMass,
                   _uncutCorrection[regId] * _betaUnCut);
    }

    if (topLeft.isEmpty() && bottomRight.isEmpty())
      topLeft = bottomRight = currentBBox;
    else {
      topLeft = HBBoxUpdateTopLeft(topLeft, currentBBox);
      bottomRight = HBBoxUpdateBottomRight(bottomRight, currentBBox);
    }
  }

  _actualWl[edge.getIndex()] = actualBBox.getHalfPerim();

  BBox expBBox((topLeft.xMin + topLeft.xMax) / 2,
               (bottomRight.yMin + bottomRight.yMax) / 2,
               (bottomRight.xMin + bottomRight.xMax) / 2,
               (topLeft.yMin + topLeft.yMax) / 2);

  if (degree > MaxChengTableDegree) degree = MaxChengTableDegree;
  return expBBox.getHalfPerim() * ChengTable[degree];
}

void HBBoxNwChengPlEval::_scaleBBox(BBox& bbox, const Point& center,
                                    double scalingFactor) {
  // scale the x side
  double originalSize = bbox.xMax - bbox.xMin;
  double newSize = originalSize * scalingFactor;
  double delta = originalSize - newSize;

  if (delta < 0) {
    newSize = originalSize;
    delta = 0;
  }

  double left = center.x - bbox.xMin, right = bbox.xMax - center.x;

  if (left <= 0)  // the block is completely on the right of center line
    bbox.xMax = bbox.xMin + newSize;
  else if (right <= 0)  // the block is on the left of center line
    bbox.xMin = bbox.xMax - newSize;
  else {
    bbox.xMin = bbox.xMin + (left * delta / originalSize);
    bbox.xMax = bbox.xMin + newSize;
  }

  // scale the y side
  originalSize = bbox.yMax - bbox.yMin;
  newSize = originalSize * scalingFactor;
  delta = originalSize - newSize;

  if (delta < 0) {
    newSize = originalSize;
    delta = 0;
  }

  double bottom = center.y - bbox.yMin, top = bbox.yMax - center.y;
  if (bottom <= 0)  // the block is completely on the bottom of center line
    bbox.yMax = bbox.yMin + newSize;
  else if (top <= 0)  // the block is on the top of center line
    bbox.yMin = bbox.yMax - newSize;
  else {
    bbox.yMin = bbox.yMin + (bottom * delta / originalSize);
    bbox.yMax = bbox.yMin + newSize;
  }
}

Point HBBoxNwChengPlEval::_getCenterOfMass(const HGFEdge& edge) {
  Point centerOfMass(0, 0);

  for (itHGFNodeLocal n = edge.nodesBegin(); n != edge.nodesEnd(); n++)
    centerOfMass += _getNodesBBox(**n).getGeomCenter();

  centerOfMass.x /= edge.getDegree();
  centerOfMass.y /= edge.getDegree();

  return centerOfMass;
}

bool HBBoxNwChengPlEval::_isEdgeCut(const HGFEdge& edge) {
  itHGFNodeLocal n = edge.nodesBegin();
  unsigned firstRegion = _getRegionId(**n);
  if (firstRegion == _regions->size())  // a placed node..not in a box
    return true;
  n++;

  for (; n != edge.nodesEnd(); n++) {
    if (_getRegionId(**n) != firstRegion) return true;
  }

  return false;
}
