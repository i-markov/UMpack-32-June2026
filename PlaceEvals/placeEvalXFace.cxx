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

// created by Stefanus Mantik on 09/25/98
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "placeEvalXFace.h"

using std::cout;
using std::endl;
using std::ostream;

void PlaceEvalXFace::_print(ostream& o) const {
  for (unsigned i = 0; i < _netCost.size(); i++) {
    o << "Net: " << i << " , cost: " << _netCost[i] << endl;
  }
}

void PlaceEvalXFace::average(ostream& out, unsigned degree,
                             const LayoutBBoxes& bboxList,
                             bool includeSmallerNets) {
  double aveSum = 0, aveFormula = 0, aveRatio = 0;
  unsigned counter = 0;
  for (itHGFEdgeGlobal e = _hg.edgesBegin(); e != _hg.edgesEnd(); e++) {
    if ((*e)->getDegree() == degree ||
        (includeSmallerNets && (*e)->getDegree() < degree)) {
      BBox netBBox;
      bool terminalNet = false;
      for (itHGFNodeLocal n = (*e)->nodesBegin(); n != (*e)->nodesEnd(); n++) {
        netBBox += (*_placement)[(*n)->getIndex()];
        if (_hg.isTerminal((*n)->getIndex())) terminalNet = true;
      }
      if (terminalNet) continue;
      bool contained = false;
      for (LayoutBBoxes::const_iterator b = bboxList.begin();
           b != bboxList.end(); b++)
        if ((*b).contains(netBBox)) contained = true;
      if (contained) {
        double halfPerim = netBBox.getHalfPerim();
        double cost = _netCost[(*e)->getIndex()];
        aveSum += halfPerim;
        aveFormula += cost;
        aveRatio += halfPerim / cost;
        counter++;
      }
    }
  }
  aveSum /= counter;
  aveFormula /= counter;
  aveRatio /= counter;
  out << "Ave Sum: " << aveSum << "  Ave Formula: " << aveFormula
      << "  Ave Ratio: " << aveRatio << "  Ratio: " << aveSum / aveFormula
      << endl;
}

void PlaceEvalXFace::plot(ostream& xyStream, unsigned degree, const BBox& bbox,
                          bool includeSmallerNets) {
  for (itHGFEdgeGlobal e = _hg.edgesBegin(); e != _hg.edgesEnd(); e++) {
    if ((*e)->getDegree() == degree ||
        (includeSmallerNets && (*e)->getDegree() < degree)) {
      BBox netBBox;
      bool terminalNet = false;
      itHGFNodeLocal n;
      for (n = (*e)->nodesBegin(); n != (*e)->nodesEnd(); n++) {
        netBBox += (*_placement)[(*n)->getIndex()];
        if (_hg.isTerminal((*n)->getIndex())) terminalNet = true;
      }
      if (terminalNet) continue;
      if (bbox.contains(netBBox))
        for (n = (*e)->nodesBegin(); n != (*e)->nodesEnd(); n++)
          xyStream << (*_placement)[(*n)->getIndex()] << endl;
    }
  }
}

void PlaceEvalXFace::measure(ostream& out, unsigned degree, const BBox& bbox,
                             bool includeSmallerNets) {
  for (itHGFEdgeGlobal e = _hg.edgesBegin(); e != _hg.edgesEnd(); e++) {
    if ((*e)->getDegree() == degree ||
        (includeSmallerNets && (*e)->getDegree() < degree)) {
      BBox netBBox;
      bool terminalNet = false;
      for (itHGFNodeLocal n = (*e)->nodesBegin(); n != (*e)->nodesEnd(); n++) {
        netBBox += (*_placement)[(*n)->getIndex()];
        if (_hg.isTerminal((*n)->getIndex())) terminalNet = true;
      }
      if (terminalNet) continue;
      if (bbox.contains(netBBox))
        out << "Sum: " << netBBox.getHalfPerim()
            << "  Formula: " << _netCost[(*e)->getIndex()] << "  Ratio: "
            << netBBox.getHalfPerim() / _netCost[(*e)->getIndex()] << endl;
    }
  }
}

void PlaceEvalXFace::measure(ostream& out, unsigned degree, const BBox& bbox,
                             const LayoutBBoxes& nextLevelBBoxes,
                             bool includeSmallerNets, bool useContained) {
  for (itHGFEdgeGlobal e = _hg.edgesBegin(); e != _hg.edgesEnd(); e++) {
    if ((*e)->getDegree() == degree ||
        (includeSmallerNets && (*e)->getDegree() < degree)) {
      BBox netBBox;
      bool terminalNet = false;
      for (itHGFNodeLocal n = (*e)->nodesBegin(); n != (*e)->nodesEnd(); n++) {
        netBBox += (*_placement)[(*n)->getIndex()];
        if (_hg.isTerminal((*n)->getIndex())) terminalNet = true;
      }
      if (terminalNet) continue;
      if (bbox.contains(netBBox)) {
        bool contained = false;
        for (LayoutBBoxes::const_iterator b = nextLevelBBoxes.begin();
             b != nextLevelBBoxes.end(); b++)
          if ((*b).contains(netBBox)) contained = true;
        if ((useContained && contained) || (!useContained && !contained))
          out << "Sum: " << netBBox.getHalfPerim()
              << "  Formula: " << _netCost[(*e)->getIndex()] << "  Ratio: "
              << netBBox.getHalfPerim() / _netCost[(*e)->getIndex()] << endl;
      }
    }
  }
}
