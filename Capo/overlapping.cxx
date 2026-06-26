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

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <algorithm>

#include "Ctainers/bitBoard.h"
#include "FMPart/fmPart.h"
#include "capoBin.h"
#include "capoPlacer.h"
#include "partProbForCapo.h"

using std::cout;
using std::endl;
using std::ostream;
using std::vector;

// Created by Igor Markov on March 24, 1999
// the algorithm is gnarly, so don't try to modify for no reason

#define DEBUGINFO 0
// could be 1

class CapoBinMatchRec {
  CapoBin* _bin0;
  CapoBin* _bin1;
  double _deltaWS;

 public:
  bool used;

  CapoBinMatchRec()
      : _bin0(NULL), _bin1(NULL), _deltaWS(-DBL_MAX), used(false){};
  CapoBinMatchRec(CapoBin* b0, CapoBin* b1)
      : _bin0(b0), _bin1(b1), used(false) {
    recomputeDeltaWS();
  }
  CapoBin* getBin0() { return _bin0; }
  CapoBin* getBin1() { return _bin1; }
  unsigned getIndex0() { return _bin0->getIndex(); }
  unsigned getIndex1() { return _bin1->getIndex(); }

  double getDeltaWS() const { return _deltaWS; }
  void recomputeDeltaWS();
  bool hasSiblings() const { return _bin0->_sibling == _bin1; }
};

bool matchRecUsed(const CapoBinMatchRec& mr) { return mr.used; }

void CapoBinMatchRec::recomputeDeltaWS() {
  _deltaWS =
      fabs(_bin0->getRelativeWhitespace() - _bin1->getRelativeWhitespace()) *
      100;  // to have units in %
}

bool operator<(const CapoBinMatchRec& r1, const CapoBinMatchRec& r2) {
  return r1.getDeltaWS() > r2.getDeltaWS();
}
// so that sort() results in a *descending* order

ostream& operator<<(ostream& os, CapoBinMatchRec& curB) {
  os << " Indices: " << curB.getBin1()->getIndex() << ","
     << curB.getBin0()->getIndex() << " imbalance " << curB.getDeltaWS() << "%";
  return os;
}

void CapoPlacer::doOverlapping(vector<CapoBin*>& layout) {
  vector<CapoBinMatchRec> edges;
  vector<CapoBinMatchRec> siblings;
  edges.reserve(layout.size());     // may be several x more
  siblings.reserve(layout.size());  // may be several x more

  for (unsigned b = 0; b != layout.size(); b++) {
    CapoBin* curB = layout[b];

    if (curB->getNumCells() <= 10) continue;

    unsigned k;

    //      cout << "----------" << *curB << "-----" << endl;

    for (k = 0; k != curB->_neighborsAbove.size(); k++) {
      CapoBin* otherB = curB->_neighborsAbove[k];
      if (otherB->getNumCells() <= 10) continue;

      CapoBinMatchRec tmp(curB, otherB);
      if (tmp.hasSiblings())
        siblings.push_back(tmp);
      else
        edges.push_back(tmp);
    };

    for (k = 0; k != curB->_leftNeighbors.size(); k++) {
      CapoBin* otherB = curB->_leftNeighbors[k];
      if (otherB->getNumCells() <= 10) continue;

      CapoBinMatchRec tmp(curB, otherB);
      if (tmp.hasSiblings())
        siblings.push_back(tmp);
      else
        edges.push_back(tmp);
    };
  };

  std::sort(edges.begin(), edges.end());

  if (edges.empty()) return;

  edges.reserve(edges.size() + siblings.size());

  if (DEBUGINFO || _params.verb.getForMajStats() > 2) {
    cout << " --- Overlapping candidates : " << endl;
    if (edges.empty())
      cout << "      none" << endl;
    else
      for (unsigned k = 0; k != edges.size(); k++) cout << edges[k] << endl;

    cout << " --- Siblings : " << endl;
    if (siblings.empty())
      cout << "      none" << endl;
    else
      for (unsigned k = 0; k != siblings.size(); k++)
        cout << siblings[k] << endl;
    cout << endl;
  }

  vector<CapoBinMatchRec> matching;
  matching.reserve(layout.size());
  vector<CapoBinMatchRec>::iterator firstAvailable = edges.begin(),
                                    counter = edges.begin();
  BitBoard takenBins(layout.size());
  abkassert(layout.size(), " Empty layout ?");
  unsigned firstIndex = layout[0]->getIndex();

  /* Pass one out of four */

  while (firstAvailable < edges.end() && matching.size() < layout.size()) {
    unsigned binNum0 = firstAvailable->getIndex0() - firstIndex,
             binNum1 = firstAvailable->getIndex1() - firstIndex;
    if (takenBins.isBitSet(binNum0) || takenBins.isBitSet(binNum1))
      firstAvailable++;  // can't use firstAvailable for matching
    else                 // can use for matching
    {
      takenBins.setBit(binNum0), takenBins.setBit(binNum1);
      matching.push_back(*firstAvailable);
      firstAvailable++->used = true;
    }
  }

  edges.erase(std::remove_if(edges.begin(), edges.end(), matchRecUsed),
              edges.end());

  // now have matching, and potentially useful neighboring pairs in the range
  // [firstAvailable, edges.end()]

  if (DEBUGINFO || _params.verb.getForMajStats() > 2) {
    cout << " --- First matching : " << endl;
    if (matching.empty())
      cout << "      empty" << endl;
    else
      for (counter = matching.begin(); counter != matching.end(); counter++)
        cout << *counter << endl;
    cout << endl;
  }

  // apply matching here
  unsigned k;
  for (k = 0; k != matching.size(); k++)
    repartitionBins(matching[k].getBin0(), matching[k].getBin1());

  /* Passes two to four */

  // add siblings
  edges.insert(edges.end(), siblings.begin(), siblings.end());
  firstAvailable = edges.begin();

  unsigned passNum;
  for (passNum = 2; passNum != 5; passNum++) {
    matching.clear();
    takenBins.clear();

    for (counter = firstAvailable; counter != edges.end(); counter++)
      counter->recomputeDeltaWS();
    std::sort(firstAvailable, edges.end());

    if (DEBUGINFO || _params.verb.getForMajStats() > 2) {
      cout << " --- Candidates for matching #" << passNum << " : " << endl;
      if (firstAvailable >= edges.end())
        cout << "     none" << endl;
      else
        for (counter = firstAvailable; counter != edges.end(); counter++)
          cout << *counter << endl;
      cout << endl;
    }

    while (firstAvailable < edges.end() && matching.size() < layout.size()) {
      unsigned binNum0 = firstAvailable->getIndex0() - firstIndex,
               binNum1 = firstAvailable->getIndex1() - firstIndex;
      if (takenBins.isBitSet(binNum0) || takenBins.isBitSet(binNum1))
        firstAvailable++;  // can't use firstAvailable for matching
      else                 // can use for matching
      {
        takenBins.setBit(binNum0), takenBins.setBit(binNum1);
        matching.push_back(*firstAvailable);
        firstAvailable++->used = true;
      }
    }

    edges.erase(std::remove_if(edges.begin(), edges.end(), matchRecUsed),
                edges.end());
    firstAvailable = edges.begin();

    // now have matching, and potentially useful neighboring pairs in the range
    // [firstAvailable, edges.end()]

    if (DEBUGINFO || _params.verb.getForMajStats() > 2) {
      cout << " --- Matching # " << passNum << " : " << endl;
      if (matching.empty())
        cout << "      empty" << endl;
      else
        for (counter = matching.begin(); counter != matching.end(); counter++)
          cout << *counter << endl;
      cout << endl;
    }

    // apply matching here
    for (k = 0; k != matching.size(); k++)
      repartitionBins(matching[k].getBin0(), matching[k].getBin1());

    if (matching.empty()) break;
  }
}

void CapoPlacer::repartitionBins(CapoBin* bin0, CapoBin* bin1) {
  if (bin0->getNumCells() < 10 || bin1->getNumCells() < 10) return;

  if ((bin0->getBBox().yMin >= bin1->getBBox().yMax) ||
      (bin0->getBBox().xMax <= bin1->getBBox().xMin)) {
  }  // order is correct..do nothing
  else if ((bin1->getBBox().yMin >= bin0->getBBox().yMax) ||
           (bin1->getBBox().xMax <= bin0->getBBox().xMin)) {  // bins reversed
    std::swap(bin0, bin1);
  } else  // no clear seperating line....
  {
    Point p0Center = bin0->getCenter();
    Point p1Center = bin1->getCenter();

    if (p1Center.y > bin0->getBBox().yMax &&
        p0Center.y < bin1->getBBox().yMin) {
      std::swap(bin0, bin1);
    }
  }

  vector<const CapoBin*> bins(2);
  bins[0] = bin0;
  bins[1] = bin1;
  PartitioningProblemForCapo problem(bins, *this, _cellToBinMap,
                                     partProbEdgesVisited, thetoWeights,
                                     _params.verb);

  PartitionerParams repartParams(_params.mlParams);
  repartParams.verb = ("0_0_0");
  FMPartitioner repart(problem, _params.maxMem, repartParams);

  const SubHGraph& hgraph = problem.getSubHGraph();
  const Partitioning& soln = problem.getBestSoln();

  // produce lists of cells in bins and update _cellToBinMap
  vector<unsigned> cellsInBin0;
  cellsInBin0.reserve(soln.size() / 2);
  vector<unsigned> cellsInBin1;
  cellsInBin1.reserve(soln.size() / 2);
  Point p0Center = bin0->getCenter();
  Point p1Center = bin1->getCenter();

  unsigned nId;
  for (nId = hgraph.getNumTerminals(); nId != hgraph.getNumNodes(); nId++) {
    unsigned origId = hgraph.newNode2OrigIdx(nId);
    if (soln[nId].lowestNumPart() == 1) {
      cellsInBin1.push_back(origId);
      _cellToBinMap[origId] = bin1;
      _placement[origId] = p1Center;
    } else {
      cellsInBin0.push_back(origId);
      _cellToBinMap[origId] = bin0;
      _placement[origId] = p0Center;
    }
  }

  // update the bins and _cellToBinMap

  bin0->resetCellIds(cellsInBin0);
  bin1->resetCellIds(cellsInBin1);
}

void CapoPlacer::repartitionBins(CapoBin* bin0, CapoBin* bin1,
                                 double relCongestionBin0,
                                 double relCongestionBin1) {
  if (bin0->getNumCells() < 10 || bin1->getNumCells() < 10) return;

  if ((bin0->getBBox().yMin >= bin1->getBBox().yMax) ||
      (bin0->getBBox().xMax <= bin1->getBBox().xMin)) {
  }  // order is correct..do nothing
  else if ((bin1->getBBox().yMin >= bin0->getBBox().yMax) ||
           (bin1->getBBox().xMax <= bin0->getBBox().xMin)) {  // bins reversed
    std::swap(bin0, bin1);
  } else  // no clear seperating line....
  {
    Point p0Center = bin0->getCenter();
    Point p1Center = bin1->getCenter();

    if (p1Center.y > bin0->getBBox().yMax &&
        p0Center.y < bin1->getBBox().yMin) {
      std::swap(bin0, bin1);
    }
  }

  vector<const CapoBin*> bins(2);
  bins[0] = bin0;
  bins[1] = bin1;

  // TO DO
  unsigned numCells = bin0->getNumCells() + bin1->getNumCells();
  vector<double> partitionCapacities(2);
  partitionCapacities[0] = relCongestionBin1;
  partitionCapacities[1] = relCongestionBin0;

  PartitioningProblemForCapo problem(bins, *this, _cellToBinMap,
                                     partProbEdgesVisited, thetoWeights,
                                     partitionCapacities, _params.verb);

  const SubHGraph& hgraph = problem.getSubHGraph();

  if (numCells < 200) {
    PartitionerParams repartParams(_params.mlParams);
    repartParams.verb = ("0_0_0");
    FMPartitioner repart(problem, _params.maxMem, repartParams);
  } else {
    MLPart::Parameters mlParams;
    mlParams.verb = ("0_0_0");
    if (hgraph.getNumNodes() < 500) mlParams.clParams.sizeOfTop = 100;

    problem.reserveBuffers(2);  // sets are pairs of runs+a VC
    PartitioningBuffer refinedSolns(hgraph.getNumNodes(), 2);
    vector<double> partAreas(2);  // in the BSF soln
    unsigned bestSet;
    double origCut;
    origCut = DBL_MAX;

    vector<const char*> nodeNames;

    unsigned setNum;
    for (setNum = 0; setNum < 2; setNum++) {
      MLPart* partitioner1;

      partitioner1 = new MLPart(problem, mlParams, bbBitBoards, _params.maxMem);

      refinedSolns[setNum] = problem.getBestSoln();
      if (partitioner1->getBestSolnCost() < origCut) {
        origCut = partitioner1->getBestSolnCost();
        bestSet = setNum;
        partAreas[0] = partitioner1->getPartitionArea(0);
        partAreas[1] = partitioner1->getPartitionArea(1);
      }

      PartitionIds clearId;
      clearId.setToAll(2);
      Partitioning& buff0 = problem.getSolnBuffers()[0];
      Partitioning& buff1 = problem.getSolnBuffers()[1];
      std::fill(buff0.begin() + hgraph.getNumTerminals(), buff0.end(), clearId);
      std::fill(buff1.begin() + hgraph.getNumTerminals(), buff1.end(), clearId);
      delete partitioner1;
    }
    problem.reserveBuffers(2);
    for (setNum = 0; setNum < 2; setNum++)
      problem.getSolnBuffers()[setNum] = refinedSolns[setNum];
    problem.setBestSolnNum(bestSet);
  }

  const Partitioning& soln = problem.getBestSoln();

  // produce lists of cells in bins and update _cellToBinMap
  vector<unsigned> cellsInBin0;
  cellsInBin0.reserve(soln.size() / 2);
  vector<unsigned> cellsInBin1;
  cellsInBin1.reserve(soln.size() / 2);
  Point p0Center = bin0->getCenter();
  Point p1Center = bin1->getCenter();

  unsigned nId;
  for (nId = hgraph.getNumTerminals(); nId != hgraph.getNumNodes(); nId++) {
    unsigned origId = hgraph.newNode2OrigIdx(nId);
    if (soln[nId].lowestNumPart() == 1) {
      cellsInBin1.push_back(origId);
      _cellToBinMap[origId] = bin1;
      _placement[origId] = p1Center;
    } else {
      cellsInBin0.push_back(origId);
      _cellToBinMap[origId] = bin0;
      _placement[origId] = p0Center;
    }
  }

  // update the bins and _cellToBinMap

  bin0->resetCellIds(cellsInBin0);
  bin1->resetCellIds(cellsInBin1);
}
