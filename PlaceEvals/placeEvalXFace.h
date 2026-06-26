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

// created by Stefanus Mantik on 06/05/98

#ifndef __PLACE_EVAL_XFACE_H__
#define __PLACE_EVAL_XFACE_H__

#include "HGraph/hgFixed.h"
#include "Placement/layoutBBoxes.h"
#include "Placement/placement.h"

class PlaceEvalXFace {
 protected:
  const Placement* _placement;  // for regular bbox

  const LayoutBBoxes* _regions;  // for heuristic bbox
  std::vector<unsigned>* _assignment;
  bool _localAssignment;
  // true if the assignment is declared locally for this class
  // so at the destructor, we need to remove it from memory

  const HGraphFixed& _hg;

  std::vector<double> _netCost;  // cost per net

  virtual void _recomputeCost() = 0;  //{ abkfatal(0,"wrong call"); }

  void _print(std::ostream& o) const;

 public:
  PlaceEvalXFace(const HGraphFixed& hg, const Placement& place)
      : _placement(&place),
        _regions(NULL),
        _assignment(NULL),
        _localAssignment(false),
        _hg(hg),
        _netCost(hg.getNumEdges(), 0) {
    abkfatal(_placement->getSize() == _hg.getNumNodes(),
             "Size mismatch between HGraph and Placement");
  }

  PlaceEvalXFace(const HGraphFixed& hg, const LayoutBBoxes& regions,
                 std::vector<unsigned>& assignment)
      : _placement(NULL),
        _regions(&regions),
        _assignment(&assignment),
        _localAssignment(false),
        _hg(hg),
        _netCost(hg.getNumEdges(), 0) {
    abkfatal(_assignment->size() == _hg.getNumNodes(),
             "Size mismatch between HGraph and Assignment");
  }

  PlaceEvalXFace(const HGraphFixed& hg, const LayoutBBoxes& regions,
                 const Placement& place)
      : _placement(&place),
        _regions(&regions),
        _localAssignment(true),
        _hg(hg),
        _netCost(hg.getNumEdges(), 0) {
    unsigned regSize = _regions->size();

    _assignment = new std::vector<unsigned>(_hg.getNumNodes(), regSize);

    //'discover' the assignment of nodes to blocks
    for (unsigned index = 0; index < _hg.getNumNodes(); index++) {
      for (unsigned i = 0; i < regSize; i++) {
        if ((*_regions)[i].contains((*_placement)[index])) {
          (*_assignment)[index] = i;
          break;
        }
      }
      // abkfatal((*_assignment)[index] < regSize || _hg.isTerminal(index),
      //      "non-terminal cell is not in any block");
      abkfatal((*_assignment)[index] <= regSize, "assignment out of range");
      // assignment == regSize if node is not contained in any block
    }
  }

  virtual ~PlaceEvalXFace() {
    if (_localAssignment) delete _assignment;
  }

  virtual double getCost() {
    double totalCost = 0;
    // for(unsigned i = 0; i < _netCost.size(); i++)
    for (itHGFEdgeGlobal e = _hg.edgesBegin(); e != _hg.edgesEnd(); e++)
      totalCost += _netCost[(*e)->getIndex()];
    return totalCost;
  }

  virtual double getCutCost() {
    double totalCost = 0;
    std::vector<unsigned> assignment(_hg.getNumNodes(), _regions->size());
    for (itHGFNodeGlobal n = _hg.nodesBegin(); n != _hg.nodesEnd(); n++) {
      unsigned index = (*n)->getIndex();
      for (unsigned i = 0; i < _regions->size(); i++)
        if ((*_regions)[i].contains((*_placement)[index])) {
          assignment[(*n)->getIndex()] = i;
          break;
        }
    }

    for (itHGFEdgeGlobal e = _hg.edgesBegin(); e != _hg.edgesEnd(); e++) {
      unsigned index = (*e)->getIndex();
      unsigned currentAssignment;
      bool cutNet = false;
      for (itHGFNodeLocal nn = (*e)->nodesBegin(); nn != (*e)->nodesEnd();
           nn++) {
        if (nn == (*e)->nodesBegin())
          currentAssignment = assignment[(*nn)->getIndex()];
        else if (currentAssignment != assignment[(*nn)->getIndex()] ||
                 currentAssignment == _regions->size())
          cutNet = true;
      }
      if (cutNet) totalCost += _netCost[index];
    }
    return totalCost;
  }

  virtual double getNetCost(unsigned k) { return _netCost[k]; }

  // Reinitializes the evaluator if the placement or region (with assignment)
  // changed "non-incrementally"
  virtual void resetTo(const Placement& newPlace) {
    abkfatal(_placement != NULL, "cannot reassign a NULL placement");
    _placement = &newPlace;
    _recomputeCost();
  }
  virtual void resetTo(std::vector<unsigned>& newAssignment) {
    abkfatal(_assignment != NULL, "cannot reassign a NULL assignment");
    _assignment = &newAssignment;
    _recomputeCost();
  }
  virtual void resetTo(const LayoutBBoxes& newRegions,
                       std::vector<unsigned>& newAssignment) {
    abkfatal(_regions != NULL && _assignment != NULL,
             "cannot reassign a NULL regions list or assignment");
    _regions = &newRegions;
    _assignment = &newAssignment;
    _recomputeCost();
  }

  virtual void average(std::ostream& out, unsigned degree,
                       const LayoutBBoxes& bboxList,
                       bool includeSmallerNets = false);

  virtual void plot(std::ostream& xyStream, unsigned degree, const BBox& bbox,
                    bool includeSmallerNets = false);

  virtual void measure(std::ostream& out, unsigned degree, const BBox& bbox,
                       bool includeSmallerNets = false);

  virtual void measure(std::ostream& out, unsigned degree, const BBox& bbox,
                       const LayoutBBoxes& nextLevelBBoxes,
                       bool includeSmallerNets = false,
                       bool useContained = true);

  friend std::ostream& operator<<(std::ostream& o, const PlaceEvalXFace& plXF) {
    plXF._print(o);
    return o;
  }
};

#endif
