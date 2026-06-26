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

// branchPl.h
// Base Small Placer
//  - Small Placement Problem
//    created by Andrew Caldwell  981206

#ifndef __BRANCH_PL_H__
#define __BRANCH_PL_H__

#include <deque>

#include "Combi/partitionIds.h"
#include "baseSmallPl.h"

class SmallPlacerNetStacks;

class BranchSmallPlacer : public BaseSmallPlacer {
  typedef std::deque<unsigned> bbpl_queue;
  typedef std::vector<unsigned> bbpl_stack;

  bbpl_queue _q;
  bbpl_stack _s;
  std::vector<unsigned> _c;
  unsigned _idx;
  PartitionIds _allocSet;
  double _allocSetWL[32800];  // 32768 = 2^15

  Placement _curSoln;
  double _curWL;
  Point _nxtLoc;

  std::vector<double> _bestSolnCostStack;
  SmallPlacerNetStacks* _netStacks;

  double _minNetCost;

  void addCell();  // move a cell from the queue to the stack (permutation).
                   // updates the placement, and all other structures.
  void addCellLeft();

  void removeCell();  // move a cell off the stack back to the queue.
                      // updates the placement, and all other structures.
  void removeCellLeft();

  void initSoln();  // pseudo-qp type initial soln.

  void print();

 public:
  BranchSmallPlacer(SmallPlacementProblem& problem, Parameters params,
                    MaxMem* maxMem);

  virtual ~BranchSmallPlacer();
};

class CompareNodesByInitPl {
  const Placement& pl;

 public:
  CompareNodesByInitPl(const Placement& place) : pl(place) {}

  bool operator()(unsigned n1, unsigned n2) const {
    return pl[n1].x < pl[n2].x;
  }
};

class CompareNodesByDegree {
 public:
  CompareNodesByDegree(const std::vector<int>& edgeDegrees)
      : m_edgeDegrees(edgeDegrees) {}
  bool operator()(unsigned n1, unsigned n2) const {
    return m_edgeDegrees[n1] < m_edgeDegrees[n2];
  }

 private:
  const std::vector<int>& m_edgeDegrees;
};

#endif
