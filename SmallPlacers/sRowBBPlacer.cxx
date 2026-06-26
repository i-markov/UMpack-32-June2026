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

// Base Small Placer
//  - Small Placement Problem
//    created by Andrew Caldwell  981206

#include "sRowBBPlacer.h"

#include <numeric>

#include "spNetStacks.h"

using std::cout;
using std::endl;
using std::setw;
using std::vector;

BranchSmallPlacer::BranchSmallPlacer(SmallPlacementProblem& problem,
                                     Parameters params, MaxMem* maxMem)
    : BaseSmallPlacer(problem, params, maxMem),
      _c(problem.getNumCells() + 2),
      _curSoln(problem.getNumCells()) {
  abkfatal(problem.getNumRows() == 1,
           "BB SmallPlacer can't deal with multiple rows");

  Placement& soln = _problem.getSolnBuffer()[0].placement;
  // vector<Orient>& ori     = _problem.getSolnBuffer()[0].orientations;
  const unsigned numCells = _problem.getNumCells();
  //  const unsigned numNets  = _problem.getNumNets() ;

  // setup net stacks.
  double minWL = 0;
  double rowHeight = problem.getRows()[0].getHeight();
  if (params.minWL == Parameters::HalfRowHeight)
    minWL = rowHeight / 2;
  else if (params.minWL == Parameters::RowHeight)
    minWL = rowHeight;
  else if (params.minWL == Parameters::TwiceRowHeight)
    minWL = rowHeight * 2;

  _netStacks = new SmallPlacerNetStacks(problem, minWL);

  // setup the stack and queue to produce the ordering
  _s.reserve(numCells);

  _c[numCells] = numCells;
  _c[numCells + 1] = numCells + 1;

  for (_idx = numCells; _idx != 0;) _q.push_back(--_idx);

  double bestWL = DBL_MAX;
  _curWL = DBL_MAX;

  unsigned n;
  for (n = numCells; n > 0;) soln[--n].y = problem.getRows()[0].yMin;

  bbpl_queue __q = _q;  // __q will be used for sorting

  // generate a bound from initial solution
  Placement __curSoln = _curSoln;
  if (_problem.getSolnBuffer()[0].populated && _params.useInitSoln) {
    std::sort(__q.begin(), __q.end(), CompareNodesByInitPl(soln));
    while (!__q.empty()) {
      _s.push_back(__q.front());
      __q.pop_front();
      __curSoln[_s.back()].x = _netStacks->getLeftEdge();
      _netStacks->addNodeToLeft(_s.back());
    }
    _curWL = _netStacks->getCost();
    while (!_s.empty()) {
      __q.push_front(_s.back());  // restores the orig. order
      _netStacks->removeNodeFromLeft(_s.back());
      _s.pop_back();
    }
  }

  initSoln();  // this populates _q with the order produced by QP

  if (1)  // ... that is if we want to run QP
  {       // generate bound from a QP solution
    while (!_q.empty()) {
      _s.push_back(_q.front());
      _q.pop_front();
      _curSoln[_s.back()].x = _netStacks->getLeftEdge();
      _netStacks->addNodeToLeft(_s.back());
    }
    bestWL = _netStacks->getCost();
    while (!_s.empty()) {
      _q.push_front(_s.back());  // restors the orig. order
      _netStacks->removeNodeFromLeft(_s.back());
      _s.pop_back();
    }
  }

  if (_params.verb.getForMajStats() > 5)
    cout << " Initial solution gave cost : " << _curWL
         << " QP gave cost : " << bestWL << endl;

  if (_curWL < bestWL) {
    //_q = __q; //nb: using "=" causes compiler errors in debug.
    for (n = 0; n < numCells; n++) _q[n] = __q[n];

    bestWL = _curWL;
    _curSoln = __curSoln;
    if (_params.verb.getForMajStats() > 5)
      cout << " Initial solution was better than QP" << endl;
  }

  // generate a bound from initial ordering
  unsigned c;
  for (c = 0; c != numCells; c++) {
    __curSoln[c].x = _netStacks->getLeftEdge();
    _netStacks->addNodeToLeft(c);
  }
  _curWL = _netStacks->getCost();
  for (c = numCells; c != 0;) _netStacks->removeNodeFromLeft(--c);

  if (_params.verb.getForMajStats() > 5)
    cout << " initial ordering gave cost : " << _curWL;

  if (_curWL < bestWL) {
    for (c = 0; c != numCells; c++) _q[c] = c;
    bestWL = _curWL;
    _curSoln = __curSoln;
    if (_params.verb.getForMajStats() > 5)
      cout << " Initial ordering was even better" << endl;
  }

  for (c = numCells; c != 0;) {
    c--;
    soln[c].x = _curSoln[c].x;
  }
  _problem.getSolnBuffer()[0].populated = true;

  if (_params.verb.getForMajStats() > 5)
    cout << " Initial Bound: " << bestWL << endl;

  // note: hey... _q is ordered assuming left to right, but if
  //       we are placing left AND right, then need to reorder _q.
  if (_params.assignMethod != Parameters::LeftOnly) {
    bbpl_queue t = _q;

    unsigned i = 0, j = numCells - 1, k = 0;
    for (_idx = numCells; _idx != 0; k++) {
      _idx--;

      if (_idx % 2)
        _q[k] = t[j--];  // would be add to right.
      else
        _q[k] = t[i++];  // would be add to left.
    }
  }

  _curWL = _netStacks->getCost();

  _idx = numCells - 1;

  unsigned numAdds = 0;

  //  cout << " curWL " << _curWL << " bestWL " << bestWL << endl;

  if (params.subsetBound)
    abkfatal(numCells <= 15, "Can't take more than 15 cells at the moment\n");
  bool poorPrefix = false;             // these 3 lines only needed
  unsigned k = 0, ns = 1 << numCells;  // for subset bounds
  for (k = 0; k != ns; k++) _allocSetWL[k] = DBL_MAX;

  maxMem->update("Branch and Bound Single Row SmallPlacer");

  while (_idx < numCells) {
    addCell();

    numAdds++;
    /*
            cout << " Stack: " <<  _s;
            cout << " Allocated set    " << _allocSet;
            cout << " Curr WL " << _curWL
                 << ", best for this set       "
                 << _allocSetWL[_allocSet.getUnsigned()]
                 << endl << endl;
    */

    if (params.subsetBound) {
      if (_allocSetWL[_allocSet.getUnsigned()] > _curWL) {
        _allocSetWL[_allocSet.getUnsigned()] = _curWL;
        poorPrefix = false;
      } else
        poorPrefix = true;
    }

    if (_c[_idx] == 0 || _curWL >= bestWL || poorPrefix) {
      _c[_idx] = 0;  // force a bound

      if (!poorPrefix &&
          _curWL <
              bestWL - 0.1e-8)  // got here if:
                                //  new best complete soln (curWL < best)
                                //  bounded partial soln (curWL > best)
                                //   so there is no need to additionally
                                //   check to ensure this is a complete soln
      {
        if (_params.verb.getForMajStats() > 5) {
          //                 cout<<" Previous best: " << bestWL
          //                     <<" improvement: " << bestWL-_curWL << endl;
          cout << " New Best: " << _curWL << " found after " << setw(5)
               << numAdds << " add operations" << endl;
        }

        bestWL = _curWL;

        unsigned cell;
        for (cell = numCells; cell != 0;) {
          cell--;
          soln[cell].x = _curSoln[cell].x;
        }
        _problem.getSolnBuffer()[0].populated = true;
      }

      while (_c[_idx] == 0)  //|| _netStacks->getCost()==_BestWL
      {
        if (_idx < numCells) removeCell();
        _c[++_idx]--;
      }
    }
    _idx--;
  }

  if (_params.verb.getForMajStats() > 5) {
    cout << " Total Add Operations: " << numAdds << endl;
    cout << "Final soln has cost: " << bestWL << endl;
  }
}

BranchSmallPlacer::~BranchSmallPlacer() {
  if (_netStacks) delete _netStacks;
}

void BranchSmallPlacer::initSoln() {
  unsigned numNets = _problem.getNumNets();
  unsigned numCells = _problem.getNumCells();
  unsigned maxNumPins = numNets * numCells;

  double xMin = _problem.getRows()[0].xMin;
  double xMax = _problem.getRows()[0].xMax;

  const SmallPlacementNetlist& netlist = _problem.getNetlist();

  vector<int> edgeDeg(numNets, 0);
  vector<int> cellDeg(numCells, 0);

  vector<double> fixedTerm(numNets, 0.);  // 0 won't influence cog.

  unsigned* adjCell = new unsigned[maxNumPins];
  unsigned** firstCell = new unsigned*[numNets + 1];

  unsigned n, c, idx;
  for (idx = 0, n = 0; n < numNets; n++) {
    const BBox& netBBox = _problem.getNetTerminalBBox(n);
    if (!netBBox.isEmpty()) {
      fixedTerm[n] = (netBBox.xMin + netBBox.xMax) / 2.;
      edgeDeg[n]++;
    }

    firstCell[n] = adjCell + idx;
    for (c = 0; c < numCells; c++) {
      // double pinOffset = netlist.getPinOffset(c,n).x ;
      // if (pinOffset != DBL_MAX)
      if (!netlist.getPinOffset(c, n).isEmpty()) {
        adjCell[idx++] = c;
        edgeDeg[n]++;
      }
    }
  }
  firstCell[n] = adjCell + idx;

  unsigned* adjNet = new unsigned[maxNumPins];
  unsigned** firstNet = new unsigned*[numCells + 1];

  for (idx = 0, c = 0; c < numCells; c++) {
    firstNet[c] = adjNet + idx;
    for (n = 0; n < numNets; n++) {
      // double pinOffset = netlist.getPinOffset(c,n).x ;
      // if (pinOffset != DBL_MAX)
      if (!netlist.getPinOffset(c, n).isEmpty()) {
        adjNet[idx++] = n;
        cellDeg[c]++;
      }
    }

    _curSoln[c].x = (xMin + xMax) / 2.;
  }
  firstNet[c] = adjNet + idx;

  // perform cog iterations.
  for (unsigned iter = 1; iter != 11; iter++) {
    for (unsigned curr = 0; curr < numCells; curr++)  // update cell i posn.
    {
      double cellPosn = 0.;
      for (unsigned* i = firstNet[curr]; i != firstNet[curr + 1]; i++) {
        n = *i;

        double netPosn = fixedTerm[n];
        for (unsigned* j = firstCell[n]; j != firstCell[n + 1]; j++) {
          c = *j;
          // netPosn += (_curSoln[c].x + netlist.getPinOffset(c,n).x) ;

          const BBox& pinOffset = netlist.getPinOffset(c, n);
          netPosn += (_curSoln[c].x + 0.5 * (pinOffset.xMin + pinOffset.xMax));
        }
        cellPosn += (netPosn / ((double)edgeDeg[n]));
      }
      if (cellDeg[curr] != 0)
        _curSoln[curr].x = cellPosn / ((double)cellDeg[curr]);
    }
  }

  // record the order.
  for (c = 0; c < numCells; c++) _q[c] = c;
  std::sort(_q.begin(), _q.end(), CompareNodesByInitPl(_curSoln));

  delete[] adjNet;
  delete[] firstNet;
  delete[] adjCell;
  delete[] firstCell;
}

void BranchSmallPlacer::addCell() {
  unsigned cellId = _q.front();
  _s.push_back(cellId);
  if (_params.subsetBound) _allocSet.addToPart(cellId);
  _q.pop_front();
  _c[_idx] = _idx;

  const unsigned cell = _s.back();

  if (_params.assignMethod == Parameters::LeftOnly) {
    _curSoln[cell].x = _netStacks->getLeftEdge();
    _curWL = _netStacks->addNodeToLeft(cell);
  } else {
    if (_idx % 2) {
      _curWL = _netStacks->addNodeToRight(cell);
      _curSoln[cell].x = _netStacks->getRightEdge();
    } else {
      _curSoln[cell].x = _netStacks->getLeftEdge();
      _curWL = _netStacks->addNodeToLeft(cell);
    }
  }
}

void BranchSmallPlacer::addCellLeft() {
  _s.push_back(_q.front());
  _q.pop_front();
  _c[_idx] = _idx;

  const unsigned cell = _s.back();

  _curSoln[cell].x = _netStacks->getLeftEdge();
  _curWL = _netStacks->addNodeToLeft(cell);
}

void BranchSmallPlacer::removeCell() {
  unsigned cell = _s.back();
  if (_params.subsetBound) _allocSet.removeFromPart(cell);
  _q.push_back(cell);
  _s.pop_back();

  if (_params.assignMethod == Parameters::LeftOnly) {
    _curWL = _netStacks->removeNodeFromLeft(cell);
  } else {
    if (_idx % 2) {
      _curWL = _netStacks->removeNodeFromRight(cell);
    } else {
      _curWL = _netStacks->removeNodeFromLeft(cell);
    }
  }
}

void BranchSmallPlacer::removeCellLeft() {
  unsigned cell = _s.back();

  _q.push_back(cell);
  _s.pop_back();

  _curWL = _netStacks->removeNodeFromLeft(cell);
}

void BranchSmallPlacer::print() {
  cout << "Q: ";
  unsigned i;
  for (i = 0; i < _q.size(); i++) cout << _q[i] << " ";
  cout << endl;

  cout << "S: ";
  for (i = 0; i < _s.size(); i++) cout << _s[i] << " ";
  cout << endl;

  cout << "C: ";
  for (i = _c.size(); i != 0;) cout << _c[--i] << " ";
  cout << endl;

  cout << "Idx:   " << _idx << endl;
  cout << "CurWL: " << _curWL << endl << endl;

  //_netStacks->print();

  cout << "__________________________" << endl << endl;
}
