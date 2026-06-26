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

// created by Stefanus Mantik on 10/06/98

#ifndef __HBBOXNPLEVAL_H__
#define __HBBOXNPLEVAL_H__

#include <vector>

#include <ABKCommon/abktypes.h>

#include "placeEvalXFace.h"

class HBBoxNPlEval : public PlaceEvalXFace {
  double _betaCut;
  double _betaUnCut;

  std::vector<double> _cutCorrection;    // currently 1/pow(#cells, 0.15)
  std::vector<double> _uncutCorrection;  // currently 1/pow(#cells, 0.04)

  bit_vector _edgeIsCut;
  std::vector<double> _actualWl;

  virtual void _recomputeCost();
  void _computeCorrection();
  bool _isEdgeCut(const HGFEdge& edge);

  double _computeOneEdge(const HGFEdge& edge);

  Point _getCenterOfMass(
      const HGFEdge& edge);  // get the center of
                             // mass of the locations of the cells on the
                             // net, assuming each cell is in the center
                             // of it's bbox, and has wt 1.

  void _scaleBBox(BBox& bbox, const Point& center, double scalingFactor);
  // side effect - modifies bbox.
  // scale bbox according to the latest algorithm, given center and
  // correction scale.

  unsigned _getRegionId(const HGFNode& node) const {
    abkassert((*_assignment)[node.getIndex()] <= _regions->size(),
              "Assignment index too big");
    return (*_assignment)[node.getIndex()];
  }

  BBox _getNodesBBox(const HGFNode& node) const {
    unsigned regId = _getRegionId(node);
    const Point& cellLoc = (*_placement)[node.getIndex()];
    if (regId == _regions->size())
      return BBox(cellLoc);
    else
      return (*_regions)[_getRegionId(node)];
  }

 public:
  HBBoxNPlEval(const HGraphFixed& hg, const LayoutBBoxes& regions,
               std::vector<unsigned>& assignment, const Placement& placement)
      : PlaceEvalXFace(hg, regions, assignment),
        _betaCut(1.0),
        _betaUnCut(2.4),
        _edgeIsCut(hg.getNumEdges(), false),
        _actualWl(hg.getNumEdges(), 0) {
    _placement = &placement;
    _computeCorrection();
    _recomputeCost();  // computes the estimate and 'actual' WLs
  }

  HBBoxNPlEval(const HGraphFixed& hg, const LayoutBBoxes& regions,
               std::vector<unsigned>& assignment)
      : PlaceEvalXFace(hg, regions, assignment),
        _betaCut(1.0),
        _betaUnCut(2.4),
        _edgeIsCut(hg.getNumEdges(), false),
        _actualWl(hg.getNumEdges(), 0) {
    _computeCorrection();
    _recomputeCost();  // computes the estimate and 'actual' WLs
  }

  HBBoxNPlEval(const HGraphFixed& hg, const LayoutBBoxes& regions,
               const Placement& placement)
      : PlaceEvalXFace(hg, regions, placement),
        _betaCut(1.0),
        _betaUnCut(2.4),
        _edgeIsCut(hg.getNumEdges(), false),
        _actualWl(hg.getNumEdges(), 0) {
    _computeCorrection();
    _recomputeCost();
  }

  double getEstCutNetWl() {
    double estWl = 0;

    for (unsigned i = 0; i < _hg.getNumEdges(); i++) {
      if (!_edgeIsCut[i]) continue;
      estWl += _netCost[i];
    }
    return estWl;
  }

  double getCutNetWl() {
    double sumWl = 0;
    for (unsigned i = 0; i < _hg.getNumEdges(); i++) {
      if (!_edgeIsCut[i]) continue;
      sumWl += _actualWl[i];
    }

    return sumWl;
  }

  double getEstUnCutNetWl() {
    double estWl = 0;
    for (unsigned i = 0; i < _hg.getNumEdges(); i++) {
      if (_edgeIsCut[i]) continue;
      estWl += _netCost[i];
    }
    return estWl;
  }

  double getUnCutNetWl() {
    double sumWl = 0;
    for (unsigned i = 0; i < _hg.getNumEdges(); i++) {
      if (_edgeIsCut[i]) continue;
      sumWl += _actualWl[i];
    }
    return sumWl;
  }

  unsigned getNumCutNets() {
    unsigned num = 0;
    for (unsigned i = 0; i < _hg.getNumEdges(); i++)
      if (_edgeIsCut[i]) num++;

    return num;
  }

  unsigned getNumUnCutNets() {
    unsigned num = 0;
    for (unsigned i = 0; i < _hg.getNumEdges(); i++)
      if (!_edgeIsCut[i]) num++;

    return num;
  }
};

#endif
