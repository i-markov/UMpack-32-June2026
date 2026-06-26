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

#ifndef __DP_PL_H__
#define __DP_PL_H__

#include "baseSmallPl.h"

class EdgeAndOffset;   // forward decl.
class StripNetStacks;  // forward decl.

struct DPOneRowTableEntry {
  unsigned x[3];  // third ordinate for negWS cells.
  DPOneRowTableEntry* parent;

  StripNetStacks* _stacks;
  int _dir;
};

////////////////////////////////////////////////////////////////////////////////

class DPSmallPlacer : public BaseSmallPlacer {
  // list of adjacent edges and pin offsets for each cell.
  EdgeAndOffset* _adjNets;

  // _firstNets[k] is a pointer to the first net touching cell k.
  EdgeAndOffset** _firstNets;

  int _idx;

  double _curWL;
  Placement _curSoln;
  Placement _bestSolnOverPass;

  std::vector<unsigned> _q;

  std::vector<DPOneRowTableEntry*> _tablePool;

  void print();
  void buildConnectivity();
  void initialSoln();

 public:
  DPSmallPlacer(SmallPlacementProblem& problem, Parameters params,
                MaxMem* maxMem);
  double go(std::vector<unsigned>& _a, std::vector<unsigned>& _b);
  double go1(std::vector<unsigned>& _a, std::vector<unsigned>& _b);

  virtual ~DPSmallPlacer();

  DPOneRowTableEntry* getTableEntryFromPool();
  void releaseTableEntryToPool(DPOneRowTableEntry* t);
};

class DPCompareNodesByInitPl {
  const Placement& pl;

 public:
  DPCompareNodesByInitPl(const Placement& place) : pl(place) {}

  bool operator()(unsigned n1, unsigned n2) const {
    return pl[n1].x < pl[n2].x;
  }
};

////////////////////////////////////////////////////////////////////////////////

#endif
