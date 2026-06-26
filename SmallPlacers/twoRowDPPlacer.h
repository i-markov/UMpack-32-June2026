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

#ifndef __INCLUDED_TWO_ROW_DP_PLACER_H_
#define __INCLUDED_TWO_ROW_DP_PLACER_H_

#include "SmallPlacement/smallPlaceProb.h"
#include "SmallPlacement/xInterval.h"
#include "baseSmallPl.h"
#include "mbset.h"

class twoRowDPStacks;  // forward decl.

// everything DP for two row needs is in here.

////////////////////////////////////////////////////////////////////////////////

struct TwoRowEdgeAndOffset {
 public:
  unsigned edge;

  BBox pinOffset;
};

////////////////////////////////////////////////////////////////////////////////

struct tableEntry {
  unsigned x[8];
  twoRowDPStacks* st;
  tableEntry* parent;
};

////////////////////////////////////////////////////////////////////////////////

class twoRowDPSmallPlacer : public BaseSmallPlacer {
 public:
  twoRowDPSmallPlacer(SmallPlacementProblem& problem, Parameters params,
                      MaxMem* maxMem);
  virtual ~twoRowDPSmallPlacer();

 private:
  double go(std::vector<std::vector<int> >& _list);

  void print();
  void buildConnectivity();
  void initialSoln();

  void dump_open_queue(std::vector<tableEntry*>& op);

  tableEntry* getTableEntryFromPool();
  void releaseTableEntryToPool(tableEntry* p);

  class DPCompareNodesByInitPl {
    const Placement& pl;

   public:
    DPCompareNodesByInitPl(const Placement& place) : pl(place) {}

    bool operator()(unsigned n1, unsigned n2) const {
      return pl[n1].x < pl[n2].x;
    }
  };

  TwoRowEdgeAndOffset* _adjNets;     // adjacent nets and offsets for each cell.
  TwoRowEdgeAndOffset** _firstNets;  // ptr to first net touching cell k.

  int _idx;

  double _curWL;
  Placement _curSoln;

  std::vector<int> _q;

  bool _usingInitSoln;
  double _bestCost;
  double _minRowViol[2];  // row viol for the best soln.

  std::vector<tableEntry*> _tablePool;
};

////////////////////////////////////////////////////////////////////////////////

// Net stacks for the 2-row DP problem.
class twoRowDPStacks {
 public:
  twoRowDPStacks(const SmallPlacementProblem& problem,
                 TwoRowEdgeAndOffset**& firstNets)
      : _problem(problem),
        _firstNets(firstNets),
        _numNets(problem.getNumNets()),
        _numCells(problem.getNumCells()) {
    _numPins = _numNets * (_numCells + 1);

    _netStacksX = new XInterval[_numNets];
    _netStacksY = new XInterval[_numNets];

    _leftEdge[0] = problem.getRows()[0].xMin;  // working left edge.
    _leftEdge[1] = problem.getRows()[1].xMin;

    _riteEdge[0] = problem.getRows()[0]
                       .xMin;  // true left edge; i.e., left edge of cells >= 0.
    _riteEdge[1] = problem.getRows()[1].xMin;

    _bottomEdge[0] = problem.getRows()[0].yMin;
    _bottomEdge[1] = problem.getRows()[1].yMin;

    _cost = 0;
    _costX = 0;
    _costY = 0;
    _toLeft.reset();
    _nMoveable.resize(_numNets);

    unsigned c, n;
    for (n = 0; n != _numNets; n++) _nMoveable[n] = 0;

    _negWS = 0;  // nb: negWS is the amount of negative cell width (and it's
                 // negative).
    for (c = 0; c != _numCells; c++) {
      _negWS += std::min(0.0, problem.getCellWidth(c));
    }
    const SmallPlacementNetlist& netlist = problem.getNetlist();
    for (c = 0; c != _numCells; c++) {
      for (n = 0; n != _numNets; n++) {
        // double pinOffset = netlist.getPinOffset(c,n).x;
        // if (pinOffset != DBL_MAX) _nMoveable[n]++;

        if (!netlist.getPinOffset(c, n).isEmpty()) _nMoveable[n]++;
      }
    }

    double initCostX = 0., initCostY = 0.;
    double initCost = 0., leftEdge = std::min(_leftEdge[0], _leftEdge[1]);
    for (n = 0; n != _numNets; n++) {
      const BBox& netBBox = problem.getNetTerminalBBox(n);
      _netStacksX[n] = XInterval(netBBox.xMin, netBBox.xMax);
      _netStacksY[n] = XInterval(netBBox.yMin, netBBox.yMax);

      if (!_netStacksX[n].isEmpty()) {
        initCost += _netStacksX[n].getLength();
        initCostX += _netStacksX[n].getLength();
      }
      if (!_netStacksY[n].isEmpty()) {
        initCost += _netStacksY[n].getLength();
        initCostY += _netStacksY[n].getLength();
      }

      if (_nMoveable[n] != 0 && !_netStacksX[n].isEmpty()) {
        const XInterval& netIvl = _netStacksX[n];
        if (netIvl.min <= leftEdge || netIvl.max <= leftEdge) {
          _toLeft.set(n);
        }
      }
    }
    _cost = initCost;
    _costX = initCostX;
    _costY = initCostY;
  }

  ~twoRowDPStacks() {
    delete[] _netStacksX;
    delete[] _netStacksY;
  }

  void copy(const twoRowDPStacks& t) {
    // copy whatever is important.
    _leftEdge[0] = t._leftEdge[0];
    _leftEdge[1] = t._leftEdge[1];

    _riteEdge[0] = t._riteEdge[0];
    _riteEdge[1] = t._riteEdge[1];

    _cost = t._cost;
    _costX = t._costX;
    _costY = t._costY;
    _toLeft = t._toLeft;

    _negWS = t._negWS;

    unsigned i;
    for (i = 0; i != _numNets; i++) {
      _nMoveable[i] = t._nMoveable[i];
      _netStacksX[i] = t._netStacksX[i];
      _netStacksY[i] = t._netStacksY[i];
    }
  }

  void addNode(int cell, int row) {
    double newCost = _cost;
    double newCostX = _costX;
    double newCostY = _costY;

    double cellWidth = _problem.getCellWidth(cell);
    double leftEdge = _leftEdge[row];
    double bottomEdge = _problem.getRows()[row].yMin;

    for (TwoRowEdgeAndOffset* pt = _firstNets[cell]; pt != _firstNets[cell + 1];
         pt++) {
      XInterval* netIvlX = &_netStacksX[pt->edge];
      XInterval* netIvlY = &_netStacksY[pt->edge];

      if (netIvlX->isEmpty()) {
        // add both pts since maybe has width.
        *netIvlX += (leftEdge + pt->pinOffset.xMin);
        *netIvlX += (leftEdge + pt->pinOffset.xMax);

        newCost += netIvlX->getLength();
        newCostX += netIvlX->getLength();
      } else {
        double oldLen = netIvlX->getLength();
        *netIvlX += (leftEdge + pt->pinOffset.xMin);
        *netIvlX += (leftEdge + pt->pinOffset.xMax);
        double newLen = netIvlX->getLength();

        newCost += (newLen - oldLen);
        newCostX += (newLen - oldLen);
      }

      if (netIvlY->isEmpty()) {
        *netIvlY += (bottomEdge + pt->pinOffset.yMin);
        *netIvlY += (bottomEdge + pt->pinOffset.yMax);

        newCost += netIvlY->getLength();
        newCostY += netIvlY->getLength();
      } else {
        double oldLen = netIvlY->getLength();
        *netIvlY += (bottomEdge + pt->pinOffset.yMin);
        *netIvlY += (bottomEdge + pt->pinOffset.yMax);
        double newLen = netIvlY->getLength();

        newCost += (newLen - oldLen);
        newCostY += (newLen - oldLen);
      }

      _nMoveable[pt->edge]--;

      if (_nMoveable[pt->edge] != 0)
        _toLeft.set(pt->edge);
      else
        _toLeft.unset(pt->edge);
    }

    _negWS -= std::min(0.0, cellWidth);  // if cellWidth < 0, less negWS so add.

    _leftEdge[row] += cellWidth;  // can decrease with negative cell width.
    _riteEdge[row] = std::max(_leftEdge[row], _riteEdge[row]);

    // can't update past minimum of the row left edges.
    double minLeftEdge = std::min(_leftEdge[0], _leftEdge[1]);

    // nb: inevitable WL is only updated to the right... the minimum
    //     left edge is the minimum of the row left edges FURTHER REDUCED
    //     due to remaining negative cell width.

    minLeftEdge += _negWS;

    // left edges can update too.
    unsigned k;
    mbitset tmp = _toLeft;
    double farthestLeft;
    for (k = 0; tmp.isNonZero(); tmp.next(), k++) {
      if (!tmp.lastBitIsOne()) continue;

      farthestLeft = minLeftEdge + _problem.getNetlist().getMinXPinOffset(k);

      XInterval* netIvlX = &_netStacksX[k];

      if (farthestLeft >= netIvlX->min) {
        // nb: only do if farthestLeft is far enough right.
        double oldLen = netIvlX->getLength();
        *netIvlX += farthestLeft;
        double newLen = netIvlX->getLength();

        newCost += (newLen - oldLen);
        newCostX += (newLen - oldLen);
      }
    }
    _cost = newCost;
    _costX = newCostX;
    _costY = newCostY;
  }

  double getCost() { return _cost; }
  double getLeftEdge(int row) { return _leftEdge[row]; }
  double getRiteEdge(int row) { return _riteEdge[row]; }
  double getBottomEdge(int row) { return _bottomEdge[row]; }

  void print() {
    std::cout << "stacks:" << std::endl;
    unsigned n;
    for (n = 0; n != _numNets; n++) {
      std::cout << "net " << n << ": ";
      if (_netStacksX[n].isEmpty()) {
        std::cout << "[empty]";
      } else {
        std::cout << " xmin: " << _netStacksX[n].min
                  << " xmax: " << _netStacksX[n].max;
        //                    <<" wlx: " <<_netStacksX[n].getLength();
      }
      if (_netStacksY[n].isEmpty()) {
        std::cout << " [empty]" << std::endl;
      } else {
        std::cout << "  ymin: " << _netStacksY[n].min
                  << "  ymax: " << _netStacksY[n].max << std::endl;
        //                    <<"  wlx: "
        //                    <<_netStacksX[n].getLength()<<std::endl;
      }
    }
  }

 public:
  const SmallPlacementProblem& _problem;
  TwoRowEdgeAndOffset**& _firstNets;

  unsigned _numPins;
  unsigned _numNets;
  unsigned _numCells;

  XInterval* _netStacksX;
  XInterval* _netStacksY;
  std::vector<int> _nMoveable;
  double _leftEdge[2];    // two rows.
  double _riteEdge[2];    // two rows.
  double _bottomEdge[2];  // two rows.
  mbitset _toLeft;
  double _cost;
  double _costX;
  double _costY;

  double _negWS;
};

struct two_dp_table_compare
// : public binary_function<tableEntry*,tableEntry*,bool>
{
  bool operator()(tableEntry* x, tableEntry* y) {
    return x->st->getCost() > y->st->getCost();
  }
};
#endif
