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

// Created by Andrew Caldwell on 05/14/99
//  Based on netStacks.cxx by Igor Markov

#include "spNetStacks.h"

#include <iostream>

using std::cout;
using std::endl;
using std::max;
using std::setw;

SmallPlacerNetStacks::SmallPlacerNetStacks(const SmallPlacementProblem& problem,
                                           double minWL)
    : _numNets(problem.getNumNets()),
      _totalWidth(0),
      _bloatWidths(problem.getNumCells()) {
  unsigned numCells = problem.getNumCells();
  unsigned maxNumPins = _numNets * (numCells + 1);

  _netData = new XInterval[2 * maxNumPins];  // net costs
  _netStacks = new XInterval*[_numNets];     // pointers to current net costs

  _minWL.reserve(_numNets);
  for (unsigned n = 0; n < _numNets; n++) _minWL[n] = minWL;

  _adjNets = new EdgeAndOffset[maxNumPins];
  _firstNets = new EdgeAndOffset*[numCells + 1];

  _cost.reserve(numCells + 2);

  unsigned c;
  for (c = 0; c != problem.getNumCells(); c++)
    _totalWidth += problem.getCellWidth(c);

  const SmallPlacementRow& row = problem.getRows()[0];
  _leftEdge = row.xMin;
  _rightEdge = row.xMax;

  double space = _rightEdge - _leftEdge;
  double bloat = (space - _totalWidth) / numCells;

  for (c = 0; c < problem.getNumCells(); c++)
    _bloatWidths[c] = problem.getCellWidth(c) + bloat;

  _nMoveable.reserve(_numNets);
  _moveableWidth.reserve(_numNets);
  _leftPinCorrection.reserve(_numNets);
  _ritePinCorrection.reserve(_numNets);

  _toLeft = new mbitset[numCells + 1];
  _toRite = new mbitset[numCells + 1];

  _toLeftIdx = 0;
  _toRiteIdx = 0;

  populate(problem);
}

void SmallPlacerNetStacks::populate(const SmallPlacementProblem& problem) {
  double initCost = 0;
  unsigned nextEdgeStack = 0;

  XInterval rowIvl(problem.getRows()[0].xMin, problem.getRows()[0].xMax);

  unsigned n, c;
  for (n = 0; n < problem.getNumNets(); n++) {
    _nMoveable[n] = 0;

    _moveableWidth[n] = 0.;
    _leftPinCorrection[n] = 0.;
    _ritePinCorrection[n] = 0.;
  }

  const SmallPlacementNetlist& netlist = problem.getNetlist();

  unsigned idx = 0;
  for (c = 0; c < problem.getNumCells(); c++) {
    _firstNets[c] = _adjNets + idx;
    for (n = 0; n < _numNets; n++) {
      const BBox& pinOffset = netlist.getPinOffset(c, n);
      if (!pinOffset.isEmpty()) {
        _adjNets[idx].edge = n;
        _adjNets[idx].pinOffset = pinOffset;
        idx++;

        _nMoveable[n]++;
        _moveableWidth[n] += _bloatWidths[c];
        _leftPinCorrection[n] =
            max((_bloatWidths[c] - pinOffset.xMax), _leftPinCorrection[n]);
        _ritePinCorrection[n] = max(pinOffset.xMin, _ritePinCorrection[n]);
      }

      /************************
                  double pinOffset = netlist.getPinOffset(c,n).x;

                  if(pinOffset != DBL_MAX)	//they are connected...
                  {
                      _adjNets[idx].edge = n ;
                      _adjNets[idx].pinOffset = pinOffset ;
                      idx++ ;

                      _nMoveable[n]++ ;
                      _moveableWidth[n]     += _bloatWidths[c] ;
                      _leftPinCorrection[n]  =
       max((_bloatWidths[c]-pinOffset),_leftPinCorrection[n]) ;
                      _ritePinCorrection[n]  =
       max(pinOffset,_ritePinCorrection[n]) ;
                  }
       ************************/
    }
  }
  _firstNets[c] = _adjNets + idx;

  idx = 0;
  for (n = 0; n < problem.getNumNets(); n++) {
    // NB: need to compute minWL first while observing user spec.
    double minWL =
        _moveableWidth[n] - _leftPinCorrection[n] - _ritePinCorrection[n];
    _minWL[n] = max(_minWL[n], minWL);

    _netStacks[n] = _netData + nextEdgeStack;
    nextEdgeStack += (2 * problem.getNumCells());

    const BBox& netBBox = problem.getNetTerminalBBox(n);
    XInterval& netIvl = *_netStacks[n];
    netIvl = XInterval(netBBox.xMin, netBBox.xMax);

    if (!netIvl.isEmpty()) {
      double rmin = rowIvl.min;
      double& nmax = netIvl.max;
      if (nmax < rmin)
        nmax = rmin;  // same as += rowInv.min
      else if (netIvl.min > rowIvl.max) {
        netIvl.min = rowIvl.max;
      }  // same as += rowInv.max

      initCost += max(netIvl.getLength(), _minWL[n]);
    } else
      initCost += _minWL[n];
  }
  _cost.push_back(initCost);

  for (n = 0; n < _numNets; n++) {
    const XInterval& netIvl = *_netStacks[n];
    if ((!netIvl.isEmpty()) && (_nMoveable[n] != 0)) {
      if (netIvl.min <= _leftEdge || netIvl.max <= _leftEdge) {
        _toLeft[0].set(n);
      }
      if (netIvl.max >= _rightEdge || netIvl.min >= _rightEdge) {
        _toRite[0].set(n);
      }
    }
  }

  mbitset intersect = _toLeft[0] & _toRite[0];
  _toLeft[0] -= intersect;
  _toRite[0] -= intersect;
}

SmallPlacerNetStacks::~SmallPlacerNetStacks() {
  delete[] _netData;
  delete[] _netStacks;
  delete[] _adjNets;
  delete[] _firstNets;

  delete[] _toLeft;
  delete[] _toRite;
}

void SmallPlacerNetStacks::print() {
  cout << "Current WL for each net is:" << endl;

  unsigned n;
  for (n = 0; n < _numNets; n++)
    cout << setw(2) << n << ": " << _netStacks[n]->getLength() << endl;
}
