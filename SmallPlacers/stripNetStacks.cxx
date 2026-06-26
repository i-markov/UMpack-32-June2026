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

#include "stripNetStacks.h"

#include <iostream>

using std::cout;
using std::endl;
using std::min;
using std::setw;

StripNetStacks::StripNetStacks(const SmallPlacementProblem& problem,
                               EdgeAndOffset**& firstNets)
    : _numNets(problem.getNumNets()),
      _numCells(problem.getNumCells()),
      _firstNets(firstNets),
      _totalWidth(0),
      _bloatWidths(problem.getNumCells()) {
  _numPins = _numNets * (_numCells + 1);

  _netStacks = new XInterval[_numNets];  // pointers to current net costs

  _cost = 0.;

  unsigned c;  //, n ;
  for (c = 0; c != problem.getNumCells(); c++)
    _totalWidth += problem.getCellWidth(c);

  const SmallPlacementRow& row = problem.getRows()[0];
  _leftEdge = row.xMin;

  double space = row.xMax - row.xMin;
  double bloat = (space - _totalWidth) / _numCells;

  _negWS = 0.0;  // amount of negative WS in terms of bloated cell widths.
  for (c = 0; c < problem.getNumCells(); c++) {
    // if no WS, then bloat width is the cell width.
    _bloatWidths[c] = problem.getCellWidth(c) + bloat;

    _negWS += min(0.0, _bloatWidths[c]);  // as far as stacks concerned,
                                          // cell width is the bloat width.
  }

  _nMoveable.reserve(_numNets);

  _toLeft.reset();

  populate(problem);
}

void StripNetStacks::populate(const SmallPlacementProblem& problem) {
  double initCost = 0;
  //  unsigned nextEdgeStack = 0;

  XInterval rowIvl(problem.getRows()[0].xMin, problem.getRows()[0].xMax);

  const SmallPlacementNetlist& netlist = problem.getNetlist();

  unsigned n, c;
  for (n = _numNets; n--;) {
    _nMoveable[n] = 0;
    for (c = _numCells; c--;) {
      if (!netlist.getPinOffset(c, n).isEmpty()) {
        _nMoveable[n]++;
      }
    }
  }

  for (n = _numNets; n--;) {
    const BBox& netBBox = problem.getNetTerminalBBox(n);

    XInterval& netIvl = _netStacks[n];

    netIvl = XInterval(netBBox.xMin, netBBox.xMax);
    if (!netBBox.isEmpty()) {
      if (netIvl.max < rowIvl.min) {
        netIvl.max = rowIvl.min;
      } else if (netIvl.min > rowIvl.max) {
        netIvl.min = rowIvl.max;
      }

      initCost += netIvl.getLength();
    }
  }
  _cost = initCost;

  for (n = _numNets; n--;) {
    const XInterval& netIvl = _netStacks[n];
    if ((!netIvl.isEmpty()) && (_nMoveable[n] != 0)) {
      if (netIvl.min <= _leftEdge || netIvl.max <= _leftEdge) {
        _toLeft.set(n);
      }
    }
  }
}

StripNetStacks::~StripNetStacks() { delete[] _netStacks; }

void StripNetStacks::print() {
  cout << "Current WL for each net is:" << endl;

  unsigned n;
  for (n = 0; n < _numNets; n++)
    cout << setw(2) << n << " : " << setw(3) << _nMoveable[n] << " : "
         << _netStacks[n].getLength() << endl;

  cout << " bits: " << _toLeft << " ";
  cout << " cost: " << _cost << endl;
}

void StripNetStacks::copy(const StripNetStacks& t) {
  _leftEdge = t._leftEdge;

  _cost = t._cost;
  _toLeft = t._toLeft;

  _negWS = t._negWS;

  unsigned i;
  for (i = 0; i < _numNets; i++) {
    _nMoveable[i] = t._nMoveable[i];
    _netStacks[i] = t._netStacks[i];
  }
}
