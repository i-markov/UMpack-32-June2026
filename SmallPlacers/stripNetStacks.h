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

#ifndef _STRIP_NETSTACKS_SROW_H_
#define _STRIP_NETSTACKS_SROW_H_

#include <set>

#include "Ctainers/riskyStack.h"
#include "SmallPlacement/smallPlaceProb.h"
#include "SmallPlacement/xInterval.h"
#include "mbset.h"

class EdgeAndOffset {
 public:
  unsigned edge;
  BBox pinOffset;
};

class StripNetStacks {
 public:
  StripNetStacks(const SmallPlacementProblem& problem,
                 EdgeAndOffset**& firstNets);

  virtual ~StripNetStacks();

  void populate(const SmallPlacementProblem& problem);

  void copy(const StripNetStacks& t);

  double addNodeToLeft(unsigned idx);

  double getCost() const { return _cost; }

  double getLeftEdge() const { return _leftEdge; }

  void print();

 public:
  //_firstNets[k] points to the edges touching moveable cell k.
  EdgeAndOffset**& _firstNets;

  unsigned _numPins;
  unsigned _numNets;
  unsigned _numCells;

  XInterval* _netStacks;  // status (width) of an edge.
  std::vector<int> _nMoveable;
  mbitset _toLeft;   // trailing edges to the left.
  double _leftEdge;  // posn next cell added at (on the left).

  double _totalWidth;
  std::vector<double> _bloatWidths;  // bloat cells
                                      // for equal distribution of whitespace.

  double _cost;  // current cost so far.

  double _negWS;  // for support of negative width cells.
};

inline double StripNetStacks::addNodeToLeft(unsigned idx) {
  // add node to the left.

  double newCost = _cost;

  for (EdgeAndOffset* pt = _firstNets[idx]; pt != _firstNets[idx + 1]; pt++) {
    XInterval* netIvl = &_netStacks[pt->edge];

    if (netIvl->isEmpty()) {
      netIvl[0] += (_leftEdge + pt->pinOffset.xMin);
      netIvl[0] += (_leftEdge + pt->pinOffset.xMax);

      newCost += netIvl[0].getLength();  // maybe cost since equiv. pins.
    } else {
      double oldLen = netIvl[0].getLength();
      netIvl[0] +=
          (_leftEdge + pt->pinOffset.xMin);  // always safe to add both pts.
      netIvl[0] += (_leftEdge + pt->pinOffset.xMax);
      double newLen = netIvl[0].getLength();

      newCost += (newLen - oldLen);
    }

    _nMoveable[pt->edge]--;

    // track trailing edges.
    if (_nMoveable[pt->edge] != 0) {
      _toLeft.set(pt->edge);
    } else {
      _toLeft.unset(pt->edge);
    }
  }

  _negWS -= std::min(
      0.0, _bloatWidths[idx]);  // negWS doesn't increase for +ve width cells.

  _leftEdge += _bloatWidths[idx];  // next cell WILL go here.

  // edges on left frontier can update.
  unsigned k;
  mbitset tmp = _toLeft;
  for (k = 0; tmp.isNonZero(); tmp.next(), k++) {
    if (!tmp.lastBitIsOne()) continue;

    XInterval* netIvl = &_netStacks[k];

    // nb: trailing edge can't have empty interval.
    // nb: negative cells have an effect on how far a trailing edge can update.

    if (_leftEdge + _negWS >= netIvl[0].min) {
      double oldLen = netIvl[0].getLength();
      netIvl[0] += (_leftEdge + _negWS);
      double newLen = netIvl[0].getLength();

      newCost += (newLen - oldLen);
    }
  }

  _cost = newCost;

  return newCost;
}

#endif
