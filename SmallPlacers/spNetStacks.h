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
//  Derived from netStacks.h by Igor Markov

#ifndef _SPL_NETSTACKS_SROW_H_
#define _SPL_NETSTACKS_SROW_H_

#include <set>

#include "Ctainers/riskyStack.h"
#include "SmallPlacement/smallPlaceProb.h"
#include "SmallPlacement/xInterval.h"
#include "mbset.h"

typedef XInterval* ivlPt;

class SmallPlacerNetStacks {
 protected:
  struct EdgeAndOffset {
    unsigned edge;
    // double   pinOffset;
    BBox pinOffset;
  };

  std::vector<int> _nMoveable;

  std::vector<double> _moveableWidth;
  std::vector<double> _leftPinCorrection;
  std::vector<double> _ritePinCorrection;

  int _toLeftIdx;
  int _toRiteIdx;
  mbitset* _toLeft;
  mbitset* _toRite;

  unsigned _numNets;

  // each element of _netData stores the xInterval for one hyperedge
  // at some point in the recursion.
  XInterval* _netData;  // netData points to an array of XIntervals
                        // of size #edges*#nodes

  // _netStacks[k] is a pointer to the element of _netData
  // describing edge k at the current moment in the recursion.
  XInterval** _netStacks;

  //_adjNets gives, for each node, the list of adjacent nets and
  // their offsets
  EdgeAndOffset* _adjNets;

  // _firstNets[k] is a pointer to the first element of
  // _nets that refers to a net touching movable node k
  EdgeAndOffset** _firstNets;

  RiskyStack<double> _cost;

  double _leftEdge;  // when addNodeToLeft is called, the cell being placed
                     // will have its origin at leftEdge, and leftEdge
                     // will increase by the cell's width
  double _rightEdge;
  double _totalWidth;
  std::vector<double> _bloatWidths;  // bloat cells
                                      // s.t. whitespace is equally distributed

  std::vector<double> _minWL;  // min x-length any net can evaluate to

  void populate(const SmallPlacementProblem& problem);

 public:
  SmallPlacerNetStacks(const SmallPlacementProblem& problem, double minWL);

  virtual ~SmallPlacerNetStacks();

  // return value is new cost
  // the NetStacks maintain the left and right edges, at which
  // the cell will be placed
  double addNodeToLeft(unsigned idx);  // idx of the node to add
  double addNodeToRight(unsigned idx);

  double removeNodeFromLeft(unsigned idx);
  double removeNodeFromRight(unsigned idx);

  double getCost() const { return _cost.back(); }
  const RiskyStack<double>& getCostStack() const { return _cost; }

  double getLeftEdge() const { return _leftEdge; }
  double getRightEdge() const { return _rightEdge; }

  void print();
};

inline double SmallPlacerNetStacks::addNodeToLeft(unsigned idx) {
  double newCost = _cost.back();

  _toLeftIdx++;
  _toLeft[_toLeftIdx] = _toLeft[_toLeftIdx - 1];

  // NB: adding node to left can change the trailing edges on right.
  _toRiteIdx++;
  _toRite[_toRiteIdx] = _toRite[_toRiteIdx - 1];

  // cout << "adding node left: " << idx << endl ;

  double rightCellEnd = _leftEdge + _bloatWidths[idx];

  for (EdgeAndOffset* pt = _firstNets[idx]; pt != _firstNets[idx + 1]; pt++) {
    ivlPt& netIvl = _netStacks[pt->edge];
    netIvl[1] = netIvl[0];

    // if (netIvl->isEmpty()) netIvl[1]    += (_leftEdge + pt->pinOffset) ;
    if (netIvl->isEmpty())
      netIvl[1] += (_leftEdge + pt->pinOffset.xMax);
    else  // newCost does not change since the new segment has len 0
    {
      // netIvl[1]    += (_leftEdge + pt->pinOffset) ;
      netIvl[1] += (_leftEdge + pt->pinOffset.xMax);
      double oldLen = netIvl[0].getLength();
      double newLen = netIvl[1].getLength();
      if (oldLen >= _minWL[pt->edge])
        newCost += (newLen - oldLen);
      else if (newLen >= _minWL[pt->edge])
        newCost += (newLen - _minWL[pt->edge]);
      //         else newCost does not change ;-)
    }

    netIvl++;

    _nMoveable[pt->edge]--;
    _moveableWidth[pt->edge] -= _bloatWidths[idx];
    if (_nMoveable[pt->edge] != 0) {
      if (rightCellEnd > netIvl[0].max) _toLeft[_toLeftIdx].set(pt->edge);
    } else
      _toLeft[_toLeftIdx].unset(pt->edge);
  }

  mbitset intersection = _toLeft[_toLeftIdx] & _toRite[_toRiteIdx];
  _toLeft[_toLeftIdx] -= intersection;
  _toRite[_toRiteIdx] -= intersection;

  _leftEdge = rightCellEnd;  // += _bloatWidths[idx];
  // edges on left frontier can update.
  unsigned k;
  mbitset tmp = _toLeft[_toLeftIdx];
  for (k = 0; tmp.isNonZero(); tmp.next(), k++) {
    if (!tmp.lastBitIsOne()) continue;
    ivlPt& netIvl = _netStacks[k];
    netIvl[1] = netIvl[0];

    //      note: a trailing net cannot have empty interval
    //      note: add inevitable WL too.
    //      note: [MAR] above comment seems to be wrong, need to handle rare
    //      empty interval

    netIvl[1] +=
        (_leftEdge + std::max(0., (_moveableWidth[k] - _leftPinCorrection[k])));

    if (!netIvl->isEmpty()) {
      double oldLen = netIvl[0].getLength();
      double newLen = netIvl[1].getLength();

      if (oldLen >= _minWL[k])
        newCost += (newLen - oldLen);
      else if (newLen >= _minWL[k])
        newCost += (newLen - _minWL[k]);
      //          else newCost does not change ;-)
    }

    netIvl++;
  }

  _cost.push_back(newCost);

  return newCost;
}

inline double SmallPlacerNetStacks::removeNodeFromLeft(unsigned idx) {
  // cout << "deleting node from left: " << idx << endl ;

  unsigned k;
  mbitset tmp = _toLeft[_toLeftIdx];
  for (k = 0; tmp.isNonZero(); tmp.next(), k++)
    if (tmp.lastBitIsOne()) --_netStacks[k];

  _toLeftIdx--;
  _toRiteIdx--;

  for (EdgeAndOffset* pt = _firstNets[idx]; pt != _firstNets[idx + 1]; pt++) {
    --_netStacks[pt->edge];
    _nMoveable[pt->edge]++;
    _moveableWidth[pt->edge] += _bloatWidths[idx];
  }

  _leftEdge -= _bloatWidths[idx];

  _cost.pop_back();
  return _cost.back();
}

inline double SmallPlacerNetStacks::addNodeToRight(unsigned idx) {
  double newCost = _cost.back();

  _toRiteIdx++;
  _toRite[_toRiteIdx] = _toRite[_toRiteIdx - 1];

  // NB: add right can change trail left.
  _toLeftIdx++;
  _toLeft[_toLeftIdx] = _toLeft[_toLeftIdx - 1];

  // cout << "adding node rite: " << idx << endl ;

  _rightEdge -= _bloatWidths[idx];
  for (EdgeAndOffset* pt = _firstNets[idx]; pt != _firstNets[idx + 1]; pt++) {
    ivlPt& netIvl = _netStacks[pt->edge];
    netIvl[1] = netIvl[0];

    if (netIvl->isEmpty()) {
      // netIvl[1] += (_rightEdge + pt->pinOffset) ; // newCost is constant
      netIvl[1] += (_rightEdge + pt->pinOffset.xMin);  // newCost is constant
    } else {
      // netIvl[1] += (_rightEdge + pt->pinOffset) ;
      netIvl[1] += (_rightEdge + pt->pinOffset.xMin);
      double oldLen = netIvl[0].getLength();
      double newLen = netIvl[1].getLength();
      if (oldLen >= _minWL[pt->edge])
        newCost += (newLen - oldLen);
      else if (newLen >= _minWL[pt->edge])
        newCost += (newLen - _minWL[pt->edge]);
      //        else newCost is constant;
    }

    netIvl++;

    _nMoveable[pt->edge]--;
    _moveableWidth[pt->edge] -= _bloatWidths[idx];
    if (_nMoveable[pt->edge] != 0) {
      if (_rightEdge < netIvl[0].min) _toRite[_toRiteIdx].set(pt->edge);
    } else
      _toRite[_toRiteIdx].unset(pt->edge);
  }

  mbitset intersection = _toLeft[_toLeftIdx] & _toRite[_toRiteIdx];
  _toLeft[_toLeftIdx] -= intersection;
  _toRite[_toRiteIdx] -= intersection;

  // edges on right frontier can update.
  unsigned k;
  mbitset tmp = _toRite[_toRiteIdx];
  for (k = 0; tmp.isNonZero(); tmp.next(), k++) {
    if (!tmp.lastBitIsOne()) continue;
    ivlPt& netIvl = _netStacks[k];
    netIvl[1] = netIvl[0];

    //      note: a trailing net cannot have an empty interval
    //      note: add inevitable WL too.
    //      note: [MAR] above comment seems to be wrong, need to handle rare
    //      empty interval

    netIvl[1] += (_rightEdge -
                  std::max(0., (_moveableWidth[k] - _ritePinCorrection[k])));

    if (!netIvl->isEmpty()) {
      double oldLen = netIvl[0].getLength();
      double newLen = netIvl[1].getLength();

      if (oldLen >= _minWL[k])
        newCost += (newLen - oldLen);
      else if (newLen >= _minWL[k])
        newCost += (newLen - _minWL[k]);
      //          else newCost does not change ;-)
    }

    netIvl++;
  }

  _cost.push_back(newCost);
  return newCost;
}

inline double SmallPlacerNetStacks::removeNodeFromRight(unsigned idx) {
  // cout << "deleting node from rite: " << idx << endl ;

  unsigned k;
  mbitset tmp = _toRite[_toRiteIdx];
  for (k = 0; tmp.isNonZero(); tmp.next(), k++)
    if (tmp.lastBitIsOne()) --_netStacks[k];

  _toRiteIdx--;
  _toLeftIdx--;

  for (EdgeAndOffset* pt = _firstNets[idx]; pt != _firstNets[idx + 1]; pt++) {
    --_netStacks[pt->edge];
    _nMoveable[pt->edge]++;
    _moveableWidth[pt->edge] += _bloatWidths[idx];
  }

  _rightEdge += _bloatWidths[idx];
  _cost.pop_back();

  return _cost.back();
}

#endif
