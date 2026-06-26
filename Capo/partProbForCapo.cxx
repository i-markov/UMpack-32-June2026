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
#include <climits>
#include <cmath>
#include <cstdio>
#include <iomanip>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkfunc.h"
#include "AnalytPl/analytPl.h"
#ifdef USEFLUTE
#include "Flute/flute.h"
#endif
#include "GeomTrees/FastSteiner/global.h"
#include "GeomTrees/FastSteiner/mst2.h"
#include "GeomTrees/FastSteiner/stnr1.h"
#include "Partitioning/termiProp.h"
#include "capoBin.h"
#include "capoPlacer.h"
#include "partProbForCapo.h"
#include "subHGForCapo.h"

using std::cout;
using std::endl;
using std::log;
using std::max;
using std::min;
using std::vector;

void PartitioningProblemForCapo::setPartDims() {
  abkassert(_bins.size() == 2, "error: called setPartDims without 2 bins");

  _nDim = _mDim = 0;

  if (_bins[0]->getBBox().yMin >= _bins[1]->getBBox().yMax)  // hcut
  {
    _nDim = 2;
    _mDim = 1;
  }

  if (_bins[1]->getBBox().yMin >=
      _bins[0]->getBBox().yMax) {  // hcut, bins reversed
    cout << "Bins are reversed: HCut" << endl;
    cout << "P0: " << *_bins[0] << endl;
    cout << "P1: " << *_bins[1] << endl;
    abkfatal(0, "Bins are reversed: p0 should be above p1");
  }

  if (_bins[0]->getBBox().xMax <= _bins[1]->getBBox().xMin)  // vcut
  {
    _nDim = 1;
    _mDim = 2;
  }

  if (_bins[1]->getBBox().xMax <=
      _bins[0]->getBBox().xMin) {  // vcut, bins reversed
    cout << "Bins are reversed: VCut" << endl;
    cout << "P0: " << *_bins[0] << endl;
    cout << "P1: " << *_bins[1] << endl;
    abkfatal(0, "Bins are reversed: p0 should be left of p1");
  }

  if (_nDim == 0)  // no clear seperating line....
  {
    Point p0Center = _bins[0]->getCenter();
    Point p1Center = _bins[1]->getCenter();

    if (p1Center.y > _bins[0]->getBBox().yMax &&
        p0Center.y < _bins[1]->getBBox().yMin)
      abkfatal(0, "Bins are reversed: p0 should be above p1");

    if (p0Center.y > _bins[1]->getBBox().yMax &&
        p1Center.y < _bins[0]->getBBox().yMin) {
      _nDim = 2;
      _mDim = 1;
    } else {
      cout << "Bins with Horizontal Overlap" << endl;
      cout << "P0: " << *_bins[0] << endl;
      cout << "P1: " << *_bins[1] << endl;
      abkfatal(0, "bins w/ horizontal overlap");
    }
  }

  const CapoBin& bin0 = *_bins[0];
  const CapoBin& bin1 = *_bins[1];

  double C0 = bin0.getTotalCellArea();
  double C1 = bin1.getTotalCellArea();
  _totalWeight = new std::vector<double>(1, C0 + C1);

  _partitions = new std::vector<BBox>(2, BBox());
  (*_partitions)[0] = bin0.getBBox();
  (*_partitions)[1] = bin1.getBBox();

  _capacities =
      new std::vector<std::vector<double> >(2, std::vector<double>(1, 0));
  _maxCapacities =
      new std::vector<std::vector<double> >(2, std::vector<double>(1, 0));
  _minCapacities =
      new std::vector<std::vector<double> >(2, std::vector<double>(1, 0));

  double S0 = bin0.getCapacity();
  double S1 = bin1.getCapacity();
  double w0 = max(bin0.getRelativeWhitespace(), 0.0);
  double w1 = max(bin1.getRelativeWhitespace(), 0.0);

  double P0Max = max(S0 - (S0 * w0), S0 - (S0 * w1));
  double P1Max = max(S1 - (S1 * w0), S1 - (S1 * w1));
  double P0Min = min(S0 - (S0 * w0), S0 - (S0 * w1));
  double P1Min = min(S1 - (S1 * w0), S1 - (S1 * w1));

  double Target0 = 0.5 * (P0Max + P0Min);
  double Target1 = 0.5 * (P1Max + P1Min);

  (*_maxCapacities)[0][0] = P0Max;
  (*_maxCapacities)[1][0] = P1Max;
  (*_minCapacities)[0][0] = P0Min;
  (*_minCapacities)[1][0] = P1Min;
  (*_capacities)[0][0] = Target0;
  (*_capacities)[1][0] = Target1;

  if (_verb.getForActions() > 7) {
    cout << "Repartitioning " << endl;
    cout << "  P0 Capacity(S0)      " << S0 << endl;
    cout << "  P1 Capacity(S1)      " << S1 << endl;
    cout << "  Total Cell Area(C)   " << C0 + C1 << endl;
    cout << "  Original %WS         " << w0 * 100 << "/" << w1 * 100 << endl;
    cout << "  Targets              " << Target0 << "/" << Target1 << endl;
    cout << "  Max                  " << P0Max << "/" << P1Max << endl;
    cout << "  Min %WhiteSpace      " << ((S0 - P0Max) / S0) * 100 << "/"
         << ((S1 - P1Max) / S1) * 100 << endl;
    cout << "  Min                  " << P0Min << "/" << P1Min << endl;
    cout << "  Max %WhiteSpace      " << ((S0 - P0Min) / S0) * 100 << "/"
         << ((S1 - P1Min) / S1) * 100 << endl;
    cout << "  Tolerance %          " << 100.0 * (P0Max - Target0) / Target0
         << "/" << 100.0 * (P1Max - Target1) / Target1 << endl;
  }

  _setPartitionCenters();
}

void PartitioningProblemForCapo::setPartDims(const vector<double>& partDims) {
  abkassert(_bins.size() == 2, "error: called setPartDims without 2 bins");

  _nDim = _mDim = 0;

  if (_bins[0]->getBBox().yMin >= _bins[1]->getBBox().yMax)  // hcut
  {
    _nDim = 2;
    _mDim = 1;
  }

  if (_bins[1]->getBBox().yMin >=
      _bins[0]->getBBox().yMax) {  // hcut, bins reversed
    cout << "Bins are reversed: HCut" << endl;
    cout << "P0: " << *_bins[0] << endl;
    cout << "P1: " << *_bins[1] << endl;
    abkfatal(0, "Bins are reversed: p0 should be above p1");
  }

  if (_bins[0]->getBBox().xMax <= _bins[1]->getBBox().xMin)  // vcut
  {
    _nDim = 1;
    _mDim = 2;
  }

  if (_bins[1]->getBBox().xMax <=
      _bins[0]->getBBox().xMin) {  // vcut, bins reversed
    cout << "Bins are reversed: VCut" << endl;
    cout << "P0: " << *_bins[0] << endl;
    cout << "P1: " << *_bins[1] << endl;
    abkfatal(0, "Bins are reversed: p0 should be left of p1");
  }

  if (_nDim == 0)  // no clear seperating line....
  {
    Point p0Center = _bins[0]->getCenter();
    Point p1Center = _bins[1]->getCenter();

    if (p1Center.y > _bins[0]->getBBox().yMax &&
        p0Center.y < _bins[1]->getBBox().yMin)
      abkfatal(0, "Bins are reversed: p0 should be above p1");

    if (p0Center.y > _bins[1]->getBBox().yMax &&
        p1Center.y < _bins[0]->getBBox().yMin) {
      _nDim = 2;
      _mDim = 1;
    } else {
      cout << "Bins with Horizontal Overlap" << endl;
      cout << "P0: " << *_bins[0] << endl;
      cout << "P1: " << *_bins[1] << endl;
      abkfatal(0, "bins w/ horizontal overlap");
    }
  }

  const CapoBin& bin0 = *_bins[0];
  const CapoBin& bin1 = *_bins[1];

  double C0 = bin0.getTotalCellArea();
  double C1 = bin1.getTotalCellArea();
  _totalWeight = new std::vector<double>(1, C0 + C1);

  _partitions = new std::vector<BBox>(2, BBox());
  (*_partitions)[0] = bin0.getBBox();
  (*_partitions)[1] = bin1.getBBox();

  _capacities =
      new std::vector<std::vector<double> >(2, std::vector<double>(1, 0));
  _maxCapacities =
      new std::vector<std::vector<double> >(2, std::vector<double>(1, 0));
  _minCapacities =
      new std::vector<std::vector<double> >(2, std::vector<double>(1, 0));

  double S0 = bin0.getCapacity();
  double S1 = bin1.getCapacity();
  double w0 = max(bin0.getRelativeWhitespace(), 0.0);
  double w1 = max(bin1.getRelativeWhitespace(), 0.0);
  double totArea = C0 + C1;
  double totWS = S0 + S1 - C0 - C1;
  double ws0 = totWS * partDims[1];
  double ws1 = totWS * partDims[0];
  double Target0, Target1;
  if (1) {
    Target0 = 0.5 * totArea * (partDims[0] + 0.5);
    Target1 = 0.5 * totArea * (partDims[1] + 0.5);
  } else if (totWS < 0) {
    Target0 = 0.5 * totArea;
    Target1 = 0.5 * totArea;
  } else {
    Target0 = S0 - ws0;
    Target1 = S1 - ws1;
  }

  if (Target0 >= S0 || Target1 >= S1) {
    if (totWS < 0) {
      Target0 = 0.5 * totArea;
      Target1 = 0.5 * totArea;
    } else {
      Target0 = S0 - ws0;
      Target1 = S1 - ws1;
    }
  }

  if (Target0 > S0) {
    Target0 = S0;
    Target1 = totArea - Target0;
  }
  if (Target1 > S1) {
    Target1 = S1;
    Target0 = totArea - Target1;
  }

  if (Target0 == 0) {
    Target0 = 0.2 * totArea;
    Target1 = 0.8 * totArea;
  } else if (Target1 == 0) {
    Target0 = 0.8 * totArea;
    Target1 = 0.2 * totArea;
  }

  double tolPct = 0;
  double P0Max = min(Target0 * (1.0 + tolPct / 100.0), S0);
  double P0Min = Target0 * (1.0 - tolPct / 100.0);
  double P1Max = min(Target1 * (1.0 + tolPct / 100.0), S1);
  double P1Min = Target1 * (1.0 - tolPct / 100.0);

  (*_maxCapacities)[0][0] = P0Max;
  (*_maxCapacities)[1][0] = P1Max;
  (*_minCapacities)[0][0] = P0Min;
  (*_minCapacities)[1][0] = P1Min;
  (*_capacities)[0][0] = Target0;
  (*_capacities)[1][0] = Target1;

  if (1 /*_verb.getForActions() > 7*/) {
    cout << "Repartitioning " << endl;
    cout << "  P0 Capacity(S0)      " << S0 << endl;
    cout << "  P1 Capacity(S1)      " << S1 << endl;
    cout << "  Total Cell Area(C)   " << C0 + C1 << endl;
    cout << "  PartDims             " << partDims[0] << "/" << partDims[1]
         << endl;
    cout << "  Original %WS         " << w0 * 100 << "/" << w1 * 100 << endl;
    cout << "  Targets              " << Target0 << "/" << Target1 << endl;
    cout << "  Max                  " << P0Max << "/" << P1Max << endl;
    cout << "  Min %WhiteSpace      " << ((S0 - P0Max) / S0) * 100 << "/"
         << ((S1 - P1Max) / S1) * 100 << endl;
    cout << "  Min                  " << P0Min << "/" << P1Min << endl;
    cout << "  Max %WhiteSpace      " << ((S0 - P0Min) / S0) * 100 << "/"
         << ((S1 - P1Min) / S1) * 100 << endl;
    cout << "  Tolerance %          " << 100.0 * (P0Max - Target0) / Target0
         << "/" << 100.0 * (P1Max - Target1) / Target1 << endl;
  }

  _setPartitionCenters();
}

void PartitioningProblemForCapo::setBuffer() {
  unsigned numStarts = 1;

  if (_bins.size() == 1) {
    _solnBuffers = new PartitioningBuffer(_hgraph->getNumNodes(), 1);
    _solnBuffers->setBeginUsedSoln(0);
    _solnBuffers->setEndUsedSoln(1);

    PartitionIds clearId;
    clearId.setToAll(2);
    std::fill((*_solnBuffers)[0].begin(), (*_solnBuffers)[0].end(), clearId);
    // set the terminal locations
    if (_ownsHGraphs) {
      // this is only true when we build our own new hgraph
      (*_solnBuffers)[0][0].setToPart(0);
      (*_solnBuffers)[0][1].setToPart(1);
    }
  } else {
    numStarts = 1;

    _solnBuffers = new PartitioningBuffer(_hgraph->getNumNodes(), numStarts);
    _solnBuffers->setBeginUsedSoln(0);
    _solnBuffers->setEndUsedSoln(numStarts);

    for (unsigned s = 0; s < numStarts; s++) {
      // set the terminal locations
      if (_ownsHGraphs) {
        // this is only true when we build our own new hgraph
        (*_solnBuffers)[s][0].setToPart(0);
        (*_solnBuffers)[s][1].setToPart(1);
      }

      // put the movable nodes in their current bins
      itHGFNodeGlobal n = _subHG->nodesBegin();
      n += _subHG->getNumTerminals();
      for (; n != _subHG->nodesEnd(); n++) {
        unsigned origId = _subHG->newNode2OrigIdx((*n)->getIndex());

        if (_cellToBinMap[origId] == _bins[0])
          (*_solnBuffers)[s][(*n)->getIndex()].setToPart(0);
        else if (_cellToBinMap[origId] == _bins[1])
          (*_solnBuffers)[s][(*n)->getIndex()].setToPart(1);
        else
          abkfatal(_subHG->isTerminal(s), "non-terminal is not in P0 or P1");
      }
    }

    _bestSolnNum = 0;
  }

  _costs = std::vector<double>(numStarts);
  _violations = std::vector<double>(numStarts);
  _imbalances = std::vector<double>(numStarts);
}

void PartitioningProblemForCapo::addEdgePropagateTerminalsReweight(
    SubHGraphForCapo* newHG, const itHGFEdgeLocal& eItr,
    const vector<unsigned>& nonTermsOnNet, bool term0, bool term1,
    double reweightingFactor) {
  double newNetWeight = (**eItr).getWeight() * reweightingFactor;

  if (_capoPl.getParams().splitterParams.netDegreeWeighting) {
    double degree =
        pow(static_cast<double>(nonTermsOnNet.size()),
            _capoPl.getParams().splitterParams.netDegreeWeightExponent);

    newNetWeight *= degree;
  }

  if (equalDouble(newNetWeight, 0.)) return;  // a zero net weight is useless

  unsigned edgeDegree = nonTermsOnNet.size();
  if (term0) edgeDegree++;
  if (term1) edgeDegree++;

  if (edgeDegree <= 1) return;

  HGFEdge& newEdge = newHG->fastAddNewEdge(**eItr, edgeDegree);

  unsigned insertedPins = 0;
  if (term0) {
    newHG->fastAddSrcSnk(&newHG->getNodeByIdx(0), &newEdge);
    insertedPins++;
  }
  if (term1) {
    newHG->fastAddSrcSnk(&newHG->getNodeByIdx(1), &newEdge);
    insertedPins++;
  }

  vector<unsigned>::const_iterator c;
  for (c = nonTermsOnNet.begin(); c != nonTermsOnNet.end(); c++) {
    HGFNode& newNd = newHG->getNodeByIdx(newHG->origNode2NewIdx(*c));
    newHG->fastAddSrcSnk(&newNd, &newEdge);
    insertedPins++;
  }
  abkassert(insertedPins == edgeDegree,
            "reweighted Capo HGraph edge degree mismatch");

  newEdge.setWeight(newNetWeight);
}

void PartitioningProblemForCapo::buildUnweightedHGraph(bool honourGrpConstr) {
  _edgesVisited.clear();

  double p0Boundry, p1Boundry, p1Left, p0Right, p0Bottom, p1Top;
  const BBox& p0Box = (*_partitions)[0];
  const BBox& p1Box = (*_partitions)[1];

  if (_capoPl.getParams().splitterParams.newFuzzy &&
      (!_capoPl.getParams().splitterParams.newFuzzyOnlyRepart || _isRepart)) {
    double totalArea = (*_capacities)[0][0] + (*_capacities)[1][0];
    double p0MinFraction = (*_minCapacities)[0][0] / totalArea;
    double p1MinFraction = (*_minCapacities)[1][0] / totalArea;
    if (_nDim == 2)  // horizontal bisection
    {
      double ySpace = max(p0Box.yMin - p1Box.yMax, 0.0);
      p0Bottom = p0Box.yMin - 0.5 * ySpace;
      p1Top = p1Box.yMax + 0.5 * ySpace;
      double ySpan = p0Box.yMax - p1Box.yMin;
      p0Boundry = p0Box.yMax - p0MinFraction * ySpan;
      p1Boundry = p1Box.yMin + p1MinFraction * ySpan;
    } else {
      double xSpace = max(p1Box.xMin - p0Box.xMax, 0.0);
      p0Right = p0Box.xMax + 0.5 * xSpace;
      p1Left = p1Box.xMin - 0.5 * xSpace;
      double xSpan = p1Box.xMax - p0Box.xMin;
      p0Boundry = p0Box.xMin + p0MinFraction * xSpan;
      p1Boundry = p1Box.xMax - p1MinFraction * xSpan;
    }
  } else {
    const double fuzzyFactor = 0.3;
    if (_nDim == 2)  // horizontal bisection
    {
      // vert. space between the bins
      double ySpace = max(p0Box.yMin - p1Box.yMax, 0.0);

      p0Bottom = p0Box.yMin - 0.5 * ySpace;
      p1Top = p1Box.yMax + 0.5 * ySpace;

      // lowest y-loc that propagates into p0
      p0Boundry = p0Bottom + (p0Box.yMax - p0Bottom) * fuzzyFactor;
      // highest y-loc that propagates into p1
      p1Boundry = p1Top - (p1Top - p1Box.yMin) * fuzzyFactor;
    } else {
      double xSpace = max(p1Box.xMin - p0Box.xMax, 0.0);
      p0Right = p0Box.xMax + 0.5 * xSpace;
      p1Left = p1Box.xMin - 0.5 * xSpace;
      p0Boundry = p0Right - (p0Right - p0Box.xMin) * fuzzyFactor;
      p1Boundry = p1Left + (p1Box.xMax - p1Left) * fuzzyFactor;
    }
  }

  unsigned totalCells = 0;
  for (unsigned b = 0; b < _bins.size(); ++b)
    totalCells += _bins[b]->getNumCells();

  vector<unsigned> unionOfCells = _bins[0]->getCellIds();
  unionOfCells.reserve(totalCells);
  for (unsigned b = 1; b < _bins.size(); ++b)
    unionOfCells.insert(unionOfCells.end(), _bins[b]->cellIdsBegin(),
                        _bins[b]->cellIdsEnd());

  bool splitting_whole_netlist =
      (unionOfCells.size() == (_capoPl.getNetlistHGraph().getNumNodes() -
                               _capoPl.getNetlistHGraph().getNumTerminals()));

  double terminal_fraction =
      static_cast<double>(_capoPl.getNetlistHGraph().getNumTerminals()) /
      static_cast<double>(_capoPl.getNetlistHGraph().getNumNodes());

  if (splitting_whole_netlist &&
      (terminal_fraction <
       _capoPl.getParams().splitterParams.termMergeThreshold)) {
    if (_params.verb.getForActions())
      cout << "Splitting the entire netlist, saving memory by not building a "
              "new subgraph"
           << endl;

    _ownsHGraphs = false;

    _hgraph = static_cast<HGraphFixed*>(
        const_cast<HGraphWDimensions*>(&_capoPl.getNetlistHGraph()));
    _subHG = static_cast<SubHGraph*>(
        const_cast<HGraphWDimensions*>(&_capoPl.getNetlistHGraph()));

    _hgraph->temporarilyZeroOutTermWeights();

    _terminalToBlock =
        new std::vector<unsigned>(_hgraph->getNumTerminals(), UINT_MAX);
    _padBlocks = new std::vector<BBox>(_hgraph->getNumTerminals());
    PartitionIds movable;
    movable.setToAll(2);
    _fixedConstr = new Partitioning(_hgraph->getNumNodes(), movable);

    double cellWidth, cellHeight;
    for (unsigned i = 0; i < _hgraph->getNumTerminals(); ++i) {
      (*_terminalToBlock)[i] = i;

      unsigned angle = _capoPl.getPlacement().getOrient(i).getAngle();
      bool notRotated = angle == 0 || angle == 180;
      cellWidth = notRotated ? _capoPl.getNetlistHGraph().getNodeWidth(i)
                             : _capoPl.getNetlistHGraph().getNodeHeight(i);
      cellHeight = notRotated ? _capoPl.getNetlistHGraph().getNodeHeight(i)
                              : _capoPl.getNetlistHGraph().getNodeWidth(i);

      (*_padBlocks)[i] += _capoPl.getPlacement()[i];  // bottomleft
      (*_padBlocks)[i] += (_capoPl.getPlacement()[i] +
                           Point(cellWidth, cellHeight));  // topright

      // figure if this terminal goes in part 0 or 1
      const HGFNode& n = _hgraph->getNodeByIdx(i);
      const BBox& termLoc =
          (n.getDegree() > 0)
              ? _capoPl.getPinLocation(i, (*(n.edgesBegin()))->getIndex())
              : _capoPl.getPlacement()[i];
      if (_nDim == 2) {
        if (termLoc.yMin <= p1Boundry && termLoc.yMax >= p0Boundry) {
        } else if (termLoc.yMax >= p0Boundry)
          (*_fixedConstr)[i].setToPart(0);
        else if (termLoc.yMin <= p1Boundry)
          (*_fixedConstr)[i].setToPart(1);
      } else {
        if (termLoc.xMin <= p0Boundry && termLoc.xMax >= p1Boundry) {
        } else if (termLoc.xMin <= p0Boundry)
          (*_fixedConstr)[i].setToPart(0);
        else if (termLoc.xMax >= p1Boundry)
          (*_fixedConstr)[i].setToPart(1);
      }
    }

  } else {
    bool destroyNetlistHGraph =
        _size == Large &&
        static_cast<double>(unionOfCells.size()) >=
            _capoPl.getParams().keepHGraphFraction *
                static_cast<double>(_capoPl.getNetlistHGraph().getNumNodes());

    if (destroyNetlistHGraph) {
      // they are unnecessary and will come back later when the hgraph is
      // rebuilt
      _capoPl.clearNamesFromHGraph();
    }

    _ownsHGraphs = true;

    SubHGraphForCapo* newHG = new SubHGraphForCapo(2, unionOfCells.size());

    _hgraph = static_cast<HGraphFixed*>(newHG);
    _subHG = static_cast<SubHGraph*>(newHG);

    vector<unsigned>::const_iterator cItr;
    itHGFEdgeLocal eItr;

    newHG->_param.makeAllSrcSnk = true;
    unsigned newId = 2;
    for (cItr = unionOfCells.begin(); cItr != unionOfCells.end();
         ++cItr, ++newId)
      newHG->addNode(_capoPl.getNetlistHGraph(),
                     _capoPl.getNetlistHGraph().getNodeByIdx(*cItr), newId);

    vector<unsigned> nonTermsOnNet;  // ids of all nonTerminals on this edge
    unsigned compoundTerm;           // the id of the (1) terminal on this edge
    bool essential;                  // is the net being examined essential?

    unsigned numBins = _bins.size();

    for (cItr = unionOfCells.begin(); cItr != unionOfCells.end(); ++cItr) {
      const HGFNode& node = _capoPl.getNetlistHGraph().getNodeByIdx(*cItr);

      for (eItr = node.edgesBegin(); eItr != node.edgesEnd(); ++eItr) {
        if (_edgesVisited.isBitSet((*eItr)->getIndex())) continue;
        _edgesVisited.setBit((*eItr)->getIndex());
        if ((*eItr)->getDegree() > 100) continue;

        nonTermsOnNet.clear();
        compoundTerm = UINT_MAX;
        essential = true;

        PartitionIds termIds;
        itHGFNodeLocal adjNd;

        unsigned eId = (*eItr)->getIndex();
        for (adjNd = (*eItr)->nodesBegin(); adjNd != (*eItr)->nodesEnd();
             adjNd++) {
          unsigned adjIdx = (*adjNd)->getIndex();

          if (_cellToBinMap[adjIdx] == _bins[0])
            nonTermsOnNet.push_back(adjIdx);
          else if (numBins == 2 && _cellToBinMap[adjIdx] == _bins[1]) {
            nonTermsOnNet.push_back(adjIdx);
          } else {  // propagate this terminal.
            const BBox& termLoc = _capoPl.getPinLocation(adjIdx, eId);
            termIds.loadBitsFrom(0);

            if (_capoPl.getNetlistHGraph().isTerminal(adjIdx)) {
              if (_nDim > 1) {
                if (termLoc.yMax >= p0Bottom && termLoc.yMin <= p1Top) {
                  essential = false;
                  break;
                }
                if (termLoc.yMax >= p0Bottom) termIds.addToPart(0);
                if (termLoc.yMin <= p1Top) termIds.addToPart(1);
              } else {
                if (termLoc.xMax >= p1Left && termLoc.xMin <= p0Right) {
                  essential = false;
                  break;
                }
                if (termLoc.xMin <= p0Right) termIds.addToPart(0);
                if (termLoc.xMax >= p1Left) termIds.addToPart(1);
              }
            } else {
              if (_nDim > 1) {
                if (termLoc.yMax > p1Boundry && termLoc.yMin < p0Boundry) {
                  essential = false;
                  break;
                }
                if (termLoc.yMax > p1Boundry) termIds.addToPart(0);
                if (termLoc.yMin < p0Boundry) termIds.addToPart(1);
              } else {
                if (termLoc.xMax > p0Boundry && termLoc.xMin < p1Boundry) {
                  essential = false;
                  break;
                }
                if (termLoc.xMin < p1Boundry) termIds.addToPart(0);
                if (termLoc.xMax > p0Boundry) termIds.addToPart(1);
              }
            }
            if (termIds[0] && termIds[1]) continue;

            if ((termIds[0] && compoundTerm == 1) ||
                (termIds[1] && compoundTerm == 0)) {
              essential = false;
              break;
            }

            if (termIds[0])
              compoundTerm = 0;
            else if (termIds[1])
              compoundTerm = 1;
            else {
              cout << "Node " << adjIdx << endl;
              cout << "Located at " << _capoPl.getPinLocation(adjIdx, eId)
                   << endl;
              cout << "Propagated to " << termIds << endl;
              cout << "p0Boundry: " << p0Boundry << endl;
              cout << "p1Boundry: " << p1Boundry << endl;
              cout << "p0Bottom " << p0Bottom << endl;
              cout << "p1Top " << p1Top << endl;
              cout << "p0Right " << p0Right << endl;
              cout << "p1Left " << p1Left << endl;
              if (_capoPl.getNetlistHGraph().isTerminal(adjIdx))
                cout << "Node " << adjIdx << " is an actual terminal " << endl;
              if (_nDim > 1)
                cout << " Both in terms if y-loc" << endl;
              else
                cout << " Both in terms if x-loc" << endl;

              abkfatal(0, "terminal propagated nowhere");
            }
          }
        }

        if (!essential) continue;

        double factor = 1.0;
        bool t0 = compoundTerm == 0, t1 = compoundTerm == 1;
        addEdgePropagateTerminalsReweight(newHG, eItr, nonTermsOnNet, t0, t1,
                                          factor);
      }
    }

    _capoPl.getParams().maxMem->update(
        "PartProbForCapo construction before finalize");

    // destroy the original hypergraph if necessary
    if (destroyNetlistHGraph) {
      _capoPl.destroyHGraph();
    }

    newHG->finalize();

    _capoPl.getParams().maxMem->update(
        "PartProbForCapo construction after finalize");

    _terminalToBlock = new std::vector<unsigned>(2, UINT_MAX);

    (*_terminalToBlock)[0] = 0;
    (*_terminalToBlock)[1] = 1;

    _padBlocks = new std::vector<BBox>(2);

    (*_padBlocks)[0] += Point((*_partitions)[0].xMin, (*_partitions)[0].yMax);
    (*_padBlocks)[1] += Point((*_partitions)[1].xMax, (*_partitions)[1].yMin);

    PartitionIds movable;
    movable.setToAll(2);
    _fixedConstr = new Partitioning(_hgraph->getNumNodes(), movable);
    (*_fixedConstr)[0].setToPart(0);
    (*_fixedConstr)[1].setToPart(1);
  }

  if (honourGrpConstr && _capoPl.groupRegionConstr.size() > 0) {
    vector<unsigned> fixToP0;
    vector<unsigned> fixToP1;

    const BBox& p0Box = (*_partitions)[0];
    const BBox& p1Box = (*_partitions)[1];
    BBox intersectionP0;
    BBox intersectionP1;
    bool warnedOnce = false;
    bool warnedOnceMore = false;
    const RBPlacement& RBP = _capoPl.getRBP();

    vector<double> P0wGrpsIntAreas(_capoPl.groupRegionConstr.size(), -1.);
    vector<double> P1wGrpsIntAreas(_capoPl.groupRegionConstr.size(), -1.);

    double areaBlkBox = RBP.getContainedSiteAreaInBBox(p0Box) +
                        RBP.getContainedSiteAreaInBBox(p1Box);

    double area0, area1;
    for (unsigned i = 0; i < unionOfCells.size(); ++i) {
      unsigned cellId = unionOfCells[i];
      unsigned groupId = _capoPl.cellToGrpMapping[cellId];
      if (groupId != UINT_MAX) {
        const BBox& rgnBBox = _capoPl.groupRegionConstr[groupId].first;
        if (_nDim > 1)  // horizontal split
        {
          if (greaterOrEqualDouble(rgnBBox.yMin, p0Box.yMin)) {
            fixToP0.push_back(i + 2);
          } else if (lessOrEqualDouble(rgnBBox.yMax, p1Box.yMax)) {
            fixToP1.push_back(i + 2);
          } else  // bbox intersects
          {
            intersectionP0 = p0Box.intersect(rgnBBox);
            intersectionP1 = p1Box.intersect(rgnBBox);
            if (equalDouble(P0wGrpsIntAreas[groupId], -1.)) {
              area0 = RBP.getContainedSiteAreaInBBox(intersectionP0);
              P0wGrpsIntAreas[groupId] = area0;
            } else {
              area0 = P0wGrpsIntAreas[groupId];
            }

            if (equalDouble(P1wGrpsIntAreas[groupId], -1.)) {
              area1 = RBP.getContainedSiteAreaInBBox(intersectionP1);
              P1wGrpsIntAreas[groupId] = area1;
            } else {
              area1 = P1wGrpsIntAreas[groupId];
            }

            double area = area0 + area1;

            if ((area / areaBlkBox) > 0.99) {
              if (!warnedOnceMore && _verb.getForActions() > 2) {
                abkwarn(0, "Region w grp constraints completely in pl. bin ");
                cout << "Rgn BBox " << rgnBBox << " Pl. bin "
                     << _bins[0]->getBBox() << endl;
                warnedOnceMore = true;
              }
            } else if ((area0 / area) > 0.7)
              fixToP0.push_back(i + 2);
            else if ((area1 / area) > 0.7)
              fixToP1.push_back(i + 2);
            else if (!warnedOnce && _verb.getForActions() > 2) {
              abkwarn(0, "Cutting a group with regionConstraints \n");
              cout << "Rgn BBox " << rgnBBox << endl;
              warnedOnce = true;
            }

            /*
            if(intersectionP0.getArea() > intersectionP1.getArea())
              fixToP0.push_back(i+2);
            else
              fixToP1.push_back(i+2);
            */
          }
        } else  // vertical split
        {
          if (lessOrEqualDouble(rgnBBox.xMax, p0Box.xMax))
            fixToP0.push_back(i + 2);
          else if (greaterOrEqualDouble(rgnBBox.xMin, p1Box.xMin))
            fixToP1.push_back(i + 2);
          else  // bbox intersects
          {
            intersectionP0 = p0Box.intersect(rgnBBox);
            intersectionP1 = p1Box.intersect(rgnBBox);
            if (equalDouble(P0wGrpsIntAreas[groupId], -1.)) {
              area0 = RBP.getContainedSiteAreaInBBox(intersectionP0);
              P0wGrpsIntAreas[groupId] = area0;
            } else {
              area0 = P0wGrpsIntAreas[groupId];
            }
            if (equalDouble(P1wGrpsIntAreas[groupId], -1.)) {
              area1 = RBP.getContainedSiteAreaInBBox(intersectionP1);
              P1wGrpsIntAreas[groupId] = area1;
            } else {
              area1 = P1wGrpsIntAreas[groupId];
            }

            double area = area0 + area1;
            if ((area / areaBlkBox) > 0.99) {
              if (!warnedOnceMore && _verb.getForActions() > 2) {
                abkwarn(
                    0, "Region w grp constraints completely in placement bin ");
                cout << "Rgn BBox " << rgnBBox << " Pl. bin "
                     << _bins[0]->getBBox() << endl;
                warnedOnceMore = true;
              }
            } else if ((area0 / area) > 0.7)
              fixToP0.push_back(i + 2);
            else if ((area1 / area) > 0.7)
              fixToP1.push_back(i + 2);
            else if (!warnedOnce && _verb.getForActions() > 2) {
              abkwarn(0, "Cutting a group with regionConstraints ");
              cout << "Rgn BBox " << rgnBBox << endl;
              warnedOnce = true;
            }
          }
        }
      }
    }

    if (_verb.getForActions() > 1) {
      if (_nDim > 1)
        cout << "Horizontal Split " << endl;
      else
        cout << "Vertical Split " << endl;

      cout << "Num Nodes fixed to Partition 0: " << fixToP0.size() << endl;
      cout << "Num Nodes fixed to Partition 1: " << fixToP1.size() << endl;
    }

    for (unsigned i = 0; i < fixToP0.size(); ++i)
      (*_fixedConstr)[fixToP0[i]].setToPart(0);
    for (unsigned i = 0; i < fixToP1.size(); ++i)
      (*_fixedConstr)[fixToP1[i]].setToPart(1);
  }
}

void PartitioningProblemForCapo::setTolerance(double tolPct) {
  double absoluteMax = _bins[0]->getTotalCellArea();
  double minCapacity = (*_capacities)[0][0];
  for (unsigned p = 1; p < (*_capacities).size(); ++p) {
    minCapacity = std::min(minCapacity, (*_capacities)[p][0]);
  }
  // double halfdiff = 0.5*absoluteMax*tolPct;
  double halfdiff = minCapacity * tolPct;
  for (unsigned p = 0; p < (*_capacities).size(); ++p) {
    double target = (*_capacities)[p][0];
    (*_maxCapacities)[p][0] = min(absoluteMax, target + halfdiff);
    (*_minCapacities)[p][0] = max(0., target - halfdiff);
  }
}

void PartitioningProblemForCapo::increaseTolerance(double tolPct) {
  double curTol = getTolerance();
  if (lessOrEqualDouble(curTol, 0.)) {
    setTolerance(tolPct);
  } else {
    multiplyTolerance(tolPct / curTol);
  }
}

void PartitioningProblemForCapo::decreaseTolerance(double tolPct) {
  double curTol = getTolerance();
  if (lessOrEqualDouble(curTol, 0.)) {
    setTolerance(tolPct);
  } else {
    multiplyTolerance(tolPct / curTol);
  }
}

void PartitioningProblemForCapo::multiplyTolerance(double factor) {
  double absoluteMax = _bins[0]->getTotalCellArea();
  for (unsigned p = 0; p < (*_capacities).size(); ++p) {
    double multMaxDiff =
        factor * ((*_maxCapacities)[p][0] - (*_capacities)[p][0]);
    (*_maxCapacities)[p][0] =
        min(absoluteMax, (*_capacities)[p][0] + multMaxDiff);
    double multMinDiff =
        factor * ((*_capacities)[p][0] - (*_minCapacities)[p][0]);
    (*_minCapacities)[p][0] = max(0., (*_capacities)[p][0] - multMinDiff);
  }
}

double PartitioningProblemForCapo::getTolerance(unsigned wtNum) const {
  // double absoluteMax = _bins[0]->getTotalCellArea();
  double minCapacity = (*_capacities)[0][0];
  for (unsigned p = 1; p < (*_capacities).size(); ++p) {
    minCapacity = std::min(minCapacity, (*_capacities)[p][0]);
  }

  if (equalDouble(minCapacity, 0.)) return 0.;

  // double P0Tol = ( (*_maxCapacities)[0][wtNum] - (*_minCapacities)[0][wtNum]
  // ) / absoluteMax;
  double P0Tol = 0.5 *
                 ((*_maxCapacities)[0][wtNum] - (*_minCapacities)[0][wtNum]) /
                 minCapacity;

  // double P1Tol = ( (*_maxCapacities)[1][wtNum] - (*_minCapacities)[1][wtNum]
  // ) / absoluteMax;
  double P1Tol = 0.5 *
                 ((*_maxCapacities)[1][wtNum] - (*_minCapacities)[1][wtNum]) /
                 minCapacity;

  return min(P0Tol, P1Tol);
}

void PartitioningProblemForCapo::setSafeWSToleranceH(unsigned splitRow,
                                                     double safeWS) {
  const CapoBin& bin = *_bins[0];

  double p0SiteArea = 0;                          // S0
  double p1SiteArea = 0;                          // S1
  double totalCellArea = bin.getTotalCellArea();  // C
  double totalCap = bin.getCapacity();            // S

  if (splitRow != UINT_MAX && bin.getNumRows() > 1) {
    unsigned c;
    for (c = 0; c < splitRow; ++c)
      p1SiteArea += bin.getContainedAreaInCoreRow(c);
    for (c = splitRow; c < bin.getNumRows(); ++c)
      p0SiteArea += bin.getContainedAreaInCoreRow(c);
  } else {
    p0SiteArea = ceil(0.5 * totalCap);
    p1SiteArea = floor(0.5 * totalCap);
  }

  double MaxC0, MaxC1, MinC0, MinC1;

  MaxC0 = min(totalCellArea, (1 - safeWS) * p0SiteArea);
  MaxC1 = min(totalCellArea, (1 - safeWS) * p1SiteArea);
  MinC0 = totalCellArea - MaxC1;
  MinC1 = totalCellArea - MaxC0;

  (*_maxCapacities)[0][0] = MaxC0;
  (*_maxCapacities)[1][0] = MaxC1;
  (*_minCapacities)[0][0] = MinC0;
  (*_minCapacities)[1][0] = MinC1;
  (*_capacities)[0][0] = (p0SiteArea / totalCap) * totalCellArea;
  (*_capacities)[1][0] = totalCellArea - (*_capacities)[0][0];
}

void PartitioningProblemForCapo::setMinLocalWSToleranceH(unsigned splitRow,
                                                         double minLocalWS) {
  const CapoBin& bin = *_bins[0];

  double p0SiteArea = 0;                          // S0
  double p1SiteArea = 0;                          // S1
  double totalCellArea = bin.getTotalCellArea();  // C
  double totalCap = bin.getCapacity();            // S

  if (splitRow != UINT_MAX && bin.getNumRows() > 1) {
    unsigned c;
    for (c = 0; c < splitRow; ++c)
      p1SiteArea += bin.getContainedAreaInCoreRow(c);
    for (c = splitRow; c < bin.getNumRows(); ++c)
      p0SiteArea += bin.getContainedAreaInCoreRow(c);
  } else {
    p0SiteArea = ceil(0.5 * totalCap);
    p1SiteArea = floor(0.5 * totalCap);
  }

  double MaxC0, MaxC1, MinC0, MinC1;

  double relWS = bin.getRelativeWhitespace();
  double cellUtil = 1 - relWS;
  double n = ceil(log(static_cast<double>(bin.getNumRows())) / log(2.0));
  double oneOverN = 1 / n;
  double alpha = pow(minLocalWS / relWS, oneOverN);

  MaxC0 = min(totalCellArea,
              (1 - alpha) * p0SiteArea + alpha * cellUtil * p0SiteArea);
  MaxC1 = min(totalCellArea,
              (1 - alpha) * p1SiteArea + alpha * cellUtil * p1SiteArea);
  MinC0 = totalCellArea - MaxC1;
  MinC1 = totalCellArea - MaxC0;

  (*_maxCapacities)[0][0] = MaxC0;
  (*_maxCapacities)[1][0] = MaxC1;
  (*_minCapacities)[0][0] = MinC0;
  (*_minCapacities)[1][0] = MinC1;
  (*_capacities)[0][0] = (p0SiteArea / totalCap) * totalCellArea;
  (*_capacities)[1][0] = totalCellArea - (*_capacities)[0][0];
}

void PartitioningProblemForCapo::setPaperMethodToleranceH(unsigned splitRow) {
  const CapoBin& bin = *_bins[0];

  double p0SiteArea = 0;                          // S0
  double p1SiteArea = 0;                          // S1
  double totalCellArea = bin.getTotalCellArea();  // C
  double totalCap = bin.getCapacity();            // S

  if (splitRow != UINT_MAX && bin.getNumRows() > 1) {
    unsigned c;
    for (c = 0; c < splitRow; ++c)
      p1SiteArea += bin.getContainedAreaInCoreRow(c);
    for (c = splitRow; c < bin.getNumRows(); ++c)
      p0SiteArea += bin.getContainedAreaInCoreRow(c);
  } else {
    p0SiteArea = ceil(0.5 * totalCap);
    p1SiteArea = floor(0.5 * totalCap);
  }

  double MaxC0, MaxC1, MinC0, MinC1;

  if (totalCellArea < totalCap)  // not an overfilled bin
  {
    double relWS = bin.getRelativeWhitespace();
    double cellUtil = 1 - relWS;
    double n = ceil(log(static_cast<double>(bin.getNumRows())) / log(2.0));
    double oneOverNP1 = 1 / (n + 1);
    double alpha = (pow(cellUtil, oneOverNP1) - cellUtil) /
                   (relWS * pow(cellUtil, oneOverNP1));

    MaxC0 = min(totalCellArea,
                (1 - alpha) * p0SiteArea + alpha * cellUtil * p0SiteArea);
    MaxC1 = min(totalCellArea,
                (1 - alpha) * p1SiteArea + alpha * cellUtil * p1SiteArea);
    MinC0 = totalCellArea - MaxC1;
    MinC1 = totalCellArea - MaxC0;

    (*_capacities)[0][0] = (p0SiteArea / totalCap) * totalCellArea;
    (*_capacities)[1][0] = totalCellArea - (*_capacities)[0][0];
  } else  // set tight tolerance
  {
    (*_capacities)[0][0] = (p0SiteArea / totalCap) * totalCellArea;
    (*_capacities)[1][0] = totalCellArea - (*_capacities)[0][0];
    double halfMedianCellSize = 0.5 * bin.getMedianCellArea();
    MaxC0 = min(totalCellArea, (*_capacities)[0][0] + halfMedianCellSize);
    MinC0 = max(0., (*_capacities)[0][0] - halfMedianCellSize);
    MaxC1 = min(totalCellArea, (*_capacities)[1][0] + halfMedianCellSize);
    MinC1 = max(0., (*_capacities)[1][0] - halfMedianCellSize);
  }

  (*_maxCapacities)[0][0] = MaxC0;
  (*_maxCapacities)[1][0] = MaxC1;
  (*_minCapacities)[0][0] = MinC0;
  (*_minCapacities)[1][0] = MinC1;
}

void PartitioningProblemForCapo::setSafeWSToleranceV(double xSplit,
                                                     double safeWS) {
  const CapoBin& bin = *_bins[0];

  double p0SiteArea = 0;                          // S0
  double p1SiteArea = 0;                          // S1
  double totalCellArea = bin.getTotalCellArea();  // C
  double totalCap = bin.getCapacity();            // S
  vector<double> siteAreas(2, 0.0);

  if (bin.getNumRows() > 1) {
    bin.computePartAreas(xSplit, siteAreas);
  } else {
    siteAreas[0] = ceil(0.5 * totalCap);
    siteAreas[1] = floor(0.5 * totalCap);
  }
  p0SiteArea = siteAreas[0];
  p1SiteArea = siteAreas[1];

  double MaxC0, MaxC1, MinC0, MinC1;

  MaxC0 = min(totalCellArea, (1 - safeWS) * p0SiteArea);
  MaxC1 = min(totalCellArea, (1 - safeWS) * p1SiteArea);
  MinC0 = totalCellArea - MaxC1;
  MinC1 = totalCellArea - MaxC0;

  (*_maxCapacities)[0][0] = MaxC0;
  (*_maxCapacities)[1][0] = MaxC1;
  (*_minCapacities)[0][0] = MinC0;
  (*_minCapacities)[1][0] = MinC1;
  (*_capacities)[0][0] = (p0SiteArea / totalCap) * totalCellArea;
  (*_capacities)[1][0] = totalCellArea - (*_capacities)[0][0];
}

void PartitioningProblemForCapo::setMinLocalWSToleranceV(double xSplit,
                                                         double minLocalWS) {
  const CapoBin& bin = *_bins[0];

  double p0SiteArea = 0;                          // S0
  double p1SiteArea = 0;                          // S1
  double totalCellArea = bin.getTotalCellArea();  // C
  double totalCap = bin.getCapacity();            // S
  vector<double> siteAreas(2, 0.0);

  if (bin.getNumRows() > 1) {
    bin.computePartAreas(xSplit, siteAreas);
  } else {
    siteAreas[0] = ceil(0.5 * totalCap);
    siteAreas[1] = floor(0.5 * totalCap);
  }
  p0SiteArea = siteAreas[0];
  p1SiteArea = siteAreas[1];

  double MaxC0, MaxC1, MinC0, MinC1;

  double relWS = bin.getRelativeWhitespace();
  double cellUtil = 1 - relWS;
  double n = ceil(log(static_cast<double>(bin.getNumRows())) / log(2.0));
  double oneOverN = 1 / n;
  double alpha = pow(minLocalWS / relWS, oneOverN);

  MaxC0 = min(totalCellArea,
              (1 - alpha) * p0SiteArea + alpha * cellUtil * p0SiteArea);
  MaxC1 = min(totalCellArea,
              (1 - alpha) * p1SiteArea + alpha * cellUtil * p1SiteArea);
  MinC0 = totalCellArea - MaxC1;
  MinC1 = totalCellArea - MaxC0;

  (*_maxCapacities)[0][0] = MaxC0;
  (*_maxCapacities)[1][0] = MaxC1;
  (*_minCapacities)[0][0] = MinC0;
  (*_minCapacities)[1][0] = MinC1;
  (*_capacities)[0][0] = (p0SiteArea / totalCap) * totalCellArea;
  (*_capacities)[1][0] = totalCellArea - (*_capacities)[0][0];
}

void PartitioningProblemForCapo::setPaperMethodToleranceV(double xSplit) {
  const CapoBin& bin = *_bins[0];

  double p0SiteArea = 0;                          // S0
  double p1SiteArea = 0;                          // S1
  double totalCellArea = bin.getTotalCellArea();  // C
  double totalCap = bin.getCapacity();            // S
  vector<double> siteAreas(2, 0.0);

  if (bin.getNumRows() > 1) {
    bin.computePartAreas(xSplit, siteAreas);
  } else {
    siteAreas[0] = ceil(0.5 * totalCap);
    siteAreas[1] = floor(0.5 * totalCap);
  }
  p0SiteArea = siteAreas[0];
  p1SiteArea = siteAreas[1];

  double MaxC0, MaxC1, MinC0, MinC1;

  if (totalCellArea < totalCap)  // not an overfilled bin
  {
    double relWS = bin.getRelativeWhitespace();
    double cellUtil = 1 - relWS;
    double n = ceil(log(static_cast<double>(bin.getNumRows())) / log(2.0));
    double oneOverNP1 = 1 / (n + 1);
    double alpha = (pow(cellUtil, oneOverNP1) - cellUtil) /
                   (relWS * pow(cellUtil, oneOverNP1));

    MaxC0 = min(totalCellArea,
                (1 - alpha) * p0SiteArea + alpha * cellUtil * p0SiteArea);
    MaxC1 = min(totalCellArea,
                (1 - alpha) * p1SiteArea + alpha * cellUtil * p1SiteArea);
    MinC0 = totalCellArea - MaxC1;
    MinC1 = totalCellArea - MaxC0;

    (*_capacities)[0][0] = (p0SiteArea / totalCap) * totalCellArea;
    (*_capacities)[1][0] = totalCellArea - (*_capacities)[0][0];
  } else  // set tight tolerance
  {
    (*_capacities)[0][0] = (p0SiteArea / totalCap) * totalCellArea;
    (*_capacities)[1][0] = totalCellArea - (*_capacities)[0][0];
    double halfMedianCellSize = 0.5 * bin.getMedianCellArea();
    MaxC0 = min(totalCellArea, (*_capacities)[0][0] + halfMedianCellSize);
    MinC0 = max(0., (*_capacities)[0][0] - halfMedianCellSize);
    MaxC1 = min(totalCellArea, (*_capacities)[1][0] + halfMedianCellSize);
    MinC1 = max(0., (*_capacities)[1][0] - halfMedianCellSize);
  }

  (*_maxCapacities)[0][0] = MaxC0;
  (*_maxCapacities)[1][0] = MaxC1;
  (*_minCapacities)[0][0] = MinC0;
  (*_minCapacities)[1][0] = MinC1;
}

void PartitioningProblemForCapo::buildShredHGraph(
    PartitioningProblemForCapo& origPartProb) {
  _edgesVisited.clear();

  double totalCap = (*_capacities)[0][0] + (*_capacities)[1][0];
  double threshold = 0.1 * getTolerance(0) * totalCap;  // 10% of tolerance
  double shredWidth = 0.5 * threshold;
  double shredHeight = 0.5 * threshold;

  // any cell above the threshold size will be shredded

  const HGraphFixed& origHGraph = origPartProb.getHGraph();

  unsigned totalCells = 0;
  for (unsigned b = 0; b < _bins.size(); b++)
    totalCells += _bins[b]->getNumCells();

  vector<unsigned> origUnionOfCells = _bins[0]->getCellIds();
  origUnionOfCells.reserve(totalCells);
  for (unsigned b = 1; b < _bins.size(); b++)
    origUnionOfCells.insert(origUnionOfCells.end(), _bins[b]->cellIdsBegin(),
                            _bins[b]->cellIdsEnd());

  vector<bool> isNodeShredded(origUnionOfCells.size(), false);

  const HGraphWDimensions& nlhgraph = _capoPl.getNetlistHGraph();

  vector<int> unionOfCells;
  unionOfCells.reserve(origUnionOfCells.size());
  vector<unsigned>::const_iterator cItr;
  unsigned tempId = 0;
  for (cItr = origUnionOfCells.begin(); cItr != origUnionOfCells.end();
       cItr++, tempId++) {
    const HGFNode& thisNode = nlhgraph.getNodeByIdx(*cItr);
    double cellArea = nlhgraph.getWeight(thisNode.getIndex());
    double cellWidth = sqrt(cellArea);
    double cellHeight = sqrt(cellArea);
    if (cellArea > threshold) {
      isNodeShredded[tempId] = true;
      int numNewNodesHoriz = int(ceil(cellWidth / shredWidth));
      int numNewNodesVert = int(ceil(cellHeight / shredHeight));
      int numNewNodes = numNewNodesHoriz * numNewNodesVert;
      int newId = 0 - (*cItr);
      for (int i = 0; i < numNewNodes; ++i)
        unionOfCells.push_back(newId);  //-ve idx represents shredded node
    } else
      unionOfCells.push_back(*cItr);
  }

  SubHGraphForCapo* newHG = new SubHGraphForCapo(2, unionOfCells.size());

  _hgraph = static_cast<HGraphFixed*>(newHG);
  _subHG = static_cast<SubHGraph*>(newHG);

  newHG->_param.makeAllSrcSnk = true;

  OrigToNewMap orig2New;
  unsigned newId = 2;
  unsigned id = 0;
  for (cItr = origUnionOfCells.begin(); cItr != origUnionOfCells.end();
       cItr++) {
    if (!isNodeShredded[id]) {
      newHG->addNode(nlhgraph, nlhgraph.getNodeByIdx(*cItr), newId);
      newId++;
    } else {
      const HGFNode& origNode = nlhgraph.getNodeByIdx(*cItr);
      double cellArea = nlhgraph.getWeight(origNode.getIndex());
      double cellWidth = sqrt(cellArea);
      double cellHeight = sqrt(cellArea);
      int numNewNodesHoriz = int(ceil(cellWidth / shredWidth));
      int numNewNodesVert = int(ceil(cellHeight / shredHeight));
      // int numNewNodes = numNewNodesHoriz*numNewNodesVert;
      for (int i = 0; i < numNewNodesVert; ++i)
        for (int j = 0; j < numNewNodesHoriz; ++j) {
          double weight = shredWidth * shredHeight;
          newHG->addNode(weight, *cItr, newId);
          newId++;
        }
    }
    id++;
  }

  vector<unsigned> nonTermsOnNet;  // ids of all nonTerminals on this edge
  unsigned compoundTerm;           // the id of the (1) terminal on this edge
  // bool	     essential;     //is the net being examined essential?

  const SubHGraph& origSHGraph = origPartProb.getSubHGraph();
  // itHGFEdgeGlobal eItr;
  itHGFEdgeLocal eItr;
  for (eItr = origHGraph.edgesBegin(); eItr != origHGraph.edgesEnd(); eItr++) {
    nonTermsOnNet.clear();
    compoundTerm = UINT_MAX;

    // PartitionIds termIds;
    itHGFNodeLocal adjNd;

    // unsigned eId = (*eItr)->getIndex();
    for (adjNd = (*eItr)->nodesBegin(); adjNd != (*eItr)->nodesEnd(); adjNd++) {
      unsigned adjIdx = (*adjNd)->getIndex();
      if (adjIdx == 0 || adjIdx == 1) {
        compoundTerm = adjIdx;
      } else {
        unsigned origId = origSHGraph.newNode2OrigIdx(adjIdx);
        const HGFNode& origNode = nlhgraph.getNodeByIdx(origId);
        double cellArea = nlhgraph.getWeight(origNode.getIndex());
        if (cellArea < threshold) {
          nonTermsOnNet.push_back(adjIdx);
        } else {
          // double cellWidth  = sqrt(cellArea);
          // double cellHeight = sqrt(cellArea);
          // int numNewNodesHoriz = int(ceil(cellWidth/shredWidth));
          // int numNewNodesVert = int(ceil(cellHeight/shredHeight));
          // int numNewNodes = numNewNodesHoriz*numNewNodesVert;
          // shred here
          nonTermsOnNet.push_back(adjIdx);
        }
      }
    }

    unsigned edgeDegree = nonTermsOnNet.size();
    if (compoundTerm != UINT_MAX) edgeDegree++;

    HGFEdge& newEdge = newHG->fastAddNewEdge(**eItr, edgeDegree);

    unsigned insertedPins = 0;
    if (compoundTerm != UINT_MAX) {
      newHG->fastAddSrcSnk(&newHG->getNodeByIdx(compoundTerm), &newEdge);
      insertedPins++;
    }

    vector<unsigned>::iterator c;
    for (c = nonTermsOnNet.begin(); c != nonTermsOnNet.end(); c++) {
      HGFNode& newNd = newHG->getNodeByIdx(newHG->origNode2NewIdx(*c));
      newHG->fastAddSrcSnk(&newNd, &newEdge);
      insertedPins++;
    }
    abkassert(insertedPins == edgeDegree,
              "shredded Capo HGraph edge degree mismatch");
  }

  _capoPl.getParams().maxMem->update(
      "PartProbForCapo construction before finalize");

  // destroy the original hypergraph if necessary
  bool destroyNetlistHGraph =
      _size == Large &&
      static_cast<double>(unionOfCells.size()) >=
          _capoPl.getParams().keepHGraphFraction *
              static_cast<double>(_capoPl.getNetlistHGraph().getNumNodes());

  if (destroyNetlistHGraph) {
    _capoPl.destroyHGraph();
  }

  newHG->finalize();

  _capoPl.getParams().maxMem->update(
      "PartProbForCapo construction after finalize");

  _terminalToBlock = new std::vector<unsigned>(2, UINT_MAX);
  _padBlocks = new std::vector<BBox>(2);

  PartitionIds movable;
  movable.setToAll(2);
  _fixedConstr = new Partitioning(_hgraph->getNumNodes(), movable);
  (*_fixedConstr)[0].setToPart(0);
  (*_fixedConstr)[1].setToPart(1);

  (*_terminalToBlock)[0] = 0;
  (*_terminalToBlock)[1] = 1;

  (*_padBlocks)[0] += Point((*_partitions)[0].xMin, (*_partitions)[0].yMax);
  (*_padBlocks)[1] += Point((*_partitions)[1].xMax, (*_partitions)[1].yMin);
}

void PartitioningProblemForCapo::buildHGraphForPartitioning(
    bool honourGrpConstr) {
  if (_capoPl.getParams().wtCut) {
    // do sampling
    if (_capoPl.getParams().samples > 0 && _size != Small) {
      for (unsigned i = 0; i < _capoPl.getParams().samples; ++i) {
        // randomly distribute cells in all bins
        // and update the point sets

        // ugly hack
        const_cast<CapoPlacer*>(&_capoPl)->randomlyMoveCellsInTheirBins();

        // compute unified theto values (three values)
        computeUnifiedTheto();
      }

      // divide out the non zero weights
      const vector<unsigned>& visitedEdges = _edgesVisited.getIndicesOfSetBits();
      for (unsigned i = 0; i < visitedEdges.size(); ++i) {
        _weights[visitedEdges[i]][0] /=
            static_cast<double>(_capoPl.getParams().samples);
        _weights[visitedEdges[i]][1] /=
            static_cast<double>(_capoPl.getParams().samples);
        _weights[visitedEdges[i]][2] /=
            static_cast<double>(_capoPl.getParams().samples);
      }
    } else {
      // compute unified theto values (three values)
      computeUnifiedTheto();
    }

    // depending on the three values, HGraph construction starts
    constructHG(honourGrpConstr);

    // clear weights
    const vector<unsigned>& visitedEdges = _edgesVisited.getIndicesOfSetBits();
    for (unsigned i = 0; i < visitedEdges.size(); ++i) {
      _weights[visitedEdges[i]][0] = 0.;
      _weights[visitedEdges[i]][1] = 0.;
      _weights[visitedEdges[i]][2] = 0.;
    }
  } else {
    buildUnweightedHGraph(honourGrpConstr);
  }
}

double CapoPlacer::calculateSteinerFromPointSet(const vector<Point>& pointSet,
#ifdef USEFLUTE
                                                bool useFastSt, bool useFLUTE,
                                                bool useMST) const
#else
                                                bool useFastSt,
                                                bool useMST) const
#endif
{
  const unsigned size = pointSet.size();

  if (size <= 1) {
    return 0.;
  } else if (size == 2) {
    return (fabs(pointSet[0].x - pointSet[1].x) +
            fabs(pointSet[0].y - pointSet[1].y));
  } else if (useMST) {
    // call MST to get pin pairs
    vector<fastSteiner::Point> points(size);
    vector<long> parent(size);

    for (unsigned i = 0; i < size; ++i) {
      points[i].x = pointSet[i].x;
      points[i].y = pointSet[i].y;
    }

    fastSteiner::mst2(size, &points[0], &parent[0]);

    double length = 0.;

    for (unsigned long i = 0; i < size; ++i) {
      if (static_cast<unsigned long>(parent[i]) == i) continue;

      BBox pinPair;
      pinPair += pointSet[i];
      pinPair += pointSet[parent[i]];
      length += pinPair.getHalfPerim();
    }

    fastSteiner::mst2_package_done();

    return length;
  } else if (size == 3) {
    double x_span = 0;
    double y_span = 0;

    if (pointSet[0].x > pointSet[1].x) {
      if (pointSet[0].x > pointSet[2].x) {
        // pointSet[0].x is the max
        if (pointSet[1].x > pointSet[2].x) {
          // pointSet[2].x is the min
          x_span = pointSet[0].x - pointSet[2].x;
        } else {
          // pointSet[1].x is the min
          x_span = pointSet[0].x - pointSet[1].x;
        }
      } else {
        // pointSet[2] is the max
        // pointSet[1] is the min
        x_span = pointSet[2].x - pointSet[1].x;
      }
    } else {
      if (pointSet[1].x > pointSet[2].x) {
        // pointSet[1].x is the max
        if (pointSet[0].x > pointSet[2].x) {
          // pointSet[2].x is the min
          x_span = pointSet[1].x - pointSet[2].x;
        } else {
          // pointSet[0].x is the min
          x_span = pointSet[1].x - pointSet[0].x;
        }
      } else {
        // pointSet[2] is the max
        // pointSet[0] is the min
        x_span = pointSet[2].x - pointSet[0].x;
      }
    }

    if (pointSet[0].y > pointSet[1].y) {
      if (pointSet[0].y > pointSet[2].y) {
        // pointSet[0].y is the max
        if (pointSet[1].y > pointSet[2].y) {
          // pointSet[2].y is the min
          y_span = pointSet[0].y - pointSet[2].y;
        } else {
          // pointSet[1].y is the min
          y_span = pointSet[0].y - pointSet[1].y;
        }
      } else {
        // pointSet[2] is the max
        // pointSet[1] is the min
        y_span = pointSet[2].y - pointSet[1].y;
      }
    } else {
      if (pointSet[1].y > pointSet[2].y) {
        // pointSet[1].y is the max
        if (pointSet[0].y > pointSet[2].y) {
          // pointSet[2].y is the min
          y_span = pointSet[1].y - pointSet[2].y;
        } else {
          // pointSet[0].y is the min
          y_span = pointSet[1].y - pointSet[0].y;
        }
      } else {
        // pointSet[2] is the max
        // pointSet[0] is the min
        y_span = pointSet[2].y - pointSet[0].y;
      }
    }

    return (x_span + y_span);
  } else {
    double span = std::numeric_limits<double>::max();

#ifdef USEFLUTE
    if (useFastSt && (!useFLUTE || size > D))
#else
    if (useFastSt)
#endif
    {
      long n_terminals = size;
      long max_points = 4 * n_terminals;
      long flags = 0;
      long max_rounds = DEFAULT_MAX_ROUNDS;
      long max_phases_per_round = DEFAULT_PHASES_PER_ROUND;
      long max_stnr_per_round = DEFAULT_STNR_PER_ROUND;
      long max_stnr_per_phase = DEFAULT_STNR_PER_PHASE;
      double cut_off = DEFAULT_CUT_OFF;
      long n_rounds, n_phases;

      abkfatal(n_terminals <= 22, "didn't allocate enough mem for this");

      for (unsigned i = 0; i < pointSet.size(); ++i) {
        _fastStPt[i].x = pointSet[i].x;
        _fastStPt[i].y = pointSet[i].y;
      }

      long n_points = n_terminals;

      double stnr_len = fastSteiner::stnr1(
          &n_points, max_points, _fastStPt, _fastStParent, &n_rounds,
          max_rounds, &n_phases, max_phases_per_round, max_stnr_per_round,
          max_stnr_per_phase, cut_off, flags);

      span = std::min(span, stnr_len);

      n_points = n_terminals;

      std::reverse(_fastStPt, _fastStPt + n_points);

      stnr_len = fastSteiner::stnr1(
          &n_points, max_points, _fastStPt, _fastStParent, &n_rounds,
          max_rounds, &n_phases, max_phases_per_round, max_stnr_per_round,
          max_stnr_per_phase, cut_off, flags);

      span = std::min(span, stnr_len);
    }

#ifdef USEFLUTE
    if (useFLUTE) {
      double* x = new double[size];
      double* y = new double[size];
      for (unsigned i = 0; i < size; ++i) {
        x[i] = pointSet[i].x;
        y[i] = pointSet[i].y;
      }

      unsigned isLegal = 0;
      Tree fltree = flautist(size, x, y, ACCURACY, fixedObstacles,
                             relevantObstacles, isLegal);

      bool alwaysFlute = true;

      if (isLegal || alwaysFlute) {
        span = std::min(span, fltree.length);
        free(fltree.branch);
        delete[] x;
        delete[] y;
      } else {
        free(fltree.branch);
        delete[] x;
        delete[] y;
        vector<Point> pointsOnNet = pointSet;
        mst_simple steinerbuilder;
        vector<int> relevantobs(relevantObstacles.size());
        for (unsigned i = 0; i < relevantObstacles.size(); ++i) {
          relevantobs[i] = (int)relevantObstacles[i];
        }
        steinerbuilder.getStePoints(fixedObstacles, relevantobs, pointsOnNet);
        vector<int> Tree;
        double length = 0.;
        vector<double> allLength;
        steinerbuilder.calMST(fixedObstacles, pointsOnNet, Tree, length,
                              allLength);
        span = std::min(span, length);
      }
    }
#endif

    return span;
  }
}

unsigned combineFixedAndMovablePointSets(const BBox& binBBox,
                                         const vector<PtSetPoint>& fixedPS,
                                         const vector<PtSetPoint>& movablePS,
                                         vector<Point>& combinedPS) {
  vector<PtSetPoint>::const_iterator fixed = fixedPS.begin();
  vector<PtSetPoint>::const_iterator movable = movablePS.begin();
  unsigned movablesInBin = 0;

  combinedPS.clear();
  combinedPS.reserve(fixedPS.size() + movablePS.size());

  // join the two sorted vectors with a linear pass
  // but don't choose all of the movable points
  while (fixed != fixedPS.end() && movable != movablePS.end()) {
    if (fixed->pt < movable->pt) {
      combinedPS.push_back(fixed->pt);
      ++fixed;
    } else if (movable->pt < fixed->pt) {
      if (lessOrEqualDouble(binBBox.xMin, movable->pt.x) &&
          lessThanDouble(movable->pt.x, binBBox.xMax) &&
          lessOrEqualDouble(binBBox.yMin, movable->pt.y) &&
          lessThanDouble(movable->pt.y, binBBox.yMax)) {
        movablesInBin += movable->count;
      } else {
        combinedPS.push_back(movable->pt);
      }
      ++movable;
    } else {
      if (lessOrEqualDouble(binBBox.xMin, movable->pt.x) &&
          lessThanDouble(movable->pt.x, binBBox.xMax) &&
          lessOrEqualDouble(binBBox.yMin, movable->pt.y) &&
          lessThanDouble(movable->pt.y, binBBox.yMax)) {
        movablesInBin += movable->count;
      }
      combinedPS.push_back(movable->pt);
      ++fixed;
      ++movable;
    }
  }

  for (; fixed != fixedPS.end(); ++fixed) {
    combinedPS.push_back(fixed->pt);
  }

  for (; movable != movablePS.end(); ++movable) {
    if (lessOrEqualDouble(binBBox.xMin, movable->pt.x) &&
        lessThanDouble(movable->pt.x, binBBox.xMax) &&
        lessOrEqualDouble(binBBox.yMin, movable->pt.y) &&
        lessThanDouble(movable->pt.y, binBBox.yMax)) {
      movablesInBin += movable->count;
    } else {
      combinedPS.push_back(movable->pt);
    }
  }

  return movablesInBin;
}

void calculateSpansFromHPWL(const vector<Point>& pointSet,
                            const Point& p0CenterPoint,
                            const Point& p1CenterPoint,
                            const unsigned& movablesInBin, double& oneSpan,
                            double& twoSpan, double& threeSpan) {
  BBox net;
  if (pointSet.size() > 0) {
    // pointsets are sorted by x-coord
    net.xMin = pointSet.front().x;
    net.xMax = pointSet.back().x;
    net.yMin = pointSet.front().y;
    net.yMax = pointSet.front().y;
  }

  bool containsP0 = net.contains(p0CenterPoint);
  bool containsP1 = net.contains(p1CenterPoint);
  for (unsigned i = 1; i < pointSet.size(); ++i) {
    if (pointSet[i].y < net.yMin)
      net.yMin = pointSet[i].y;
    else if (pointSet[i].y > net.yMax)
      net.yMax = pointSet[i].y;
    containsP0 = containsP0 || net.contains(p0CenterPoint);
    containsP1 = containsP1 || net.contains(p1CenterPoint);
    if (containsP0 && containsP1) break;
  }

  if (containsP0) {
    oneSpan = net.getHalfPerim();
  } else {
    BBox temp(net);
    temp += p0CenterPoint;
    oneSpan = temp.getHalfPerim();
  }

  if (containsP1) {
    twoSpan = net.getHalfPerim();
  } else {
    net += p1CenterPoint;
    twoSpan = net.getHalfPerim();
  }

  if (movablesInBin > 1) {
    if (containsP0) {
      threeSpan = net.getHalfPerim();
    } else {
      net += p0CenterPoint;
      threeSpan = net.getHalfPerim();
    }
  } else {
    threeSpan = std::max(oneSpan, twoSpan);
  }
}

void PartitioningProblemForCapo::computeUnifiedTheto() {
  _edgesVisited.clear();

  const BBox& p0Box = (*_partitions)[0];
  const BBox& p1Box = (*_partitions)[1];
  Point p0CenterPoint, p1CenterPoint;
  double unitSpan;

  if (_nDim == 2)  // horizontal bisection
  {
    double p0Center = 0.5 * (p0Box.yMax - p0Box.yMin) + p0Box.yMin;
    double p1Center = 0.5 * (p1Box.yMax - p1Box.yMin) + p1Box.yMin;
    unitSpan = fabs(p0Center - p1Center);
    // BBox one includes the p0Center point
    p0CenterPoint = Point(p0Box.xMin + 0.5 * p0Box.getWidth(), p0Center);
    // BBox two includes the p1Center point
    p1CenterPoint = Point(p1Box.xMin + 0.5 * p1Box.getWidth(), p1Center);
  } else {
    double p0Center = 0.5 * (p0Box.xMax - p0Box.xMin) + p0Box.xMin;
    double p1Center = 0.5 * (p1Box.xMax - p1Box.xMin) + p1Box.xMin;
    unitSpan = fabs(p0Center - p1Center);
    // BBox one includes the p0Center point
    p0CenterPoint = Point(p0Center, p0Box.yMin + 0.5 * p0Box.getHeight());
    // BBox two includes the p1Center point
    p1CenterPoint = Point(p1Center, p1Box.yMin + 0.5 * p1Box.getHeight());
  }

  const vector<unsigned>& unionOfCells = _bins[0]->getCellIds();

  abkfatal(_bins.size() == 1, "There should only be 1 bin to be partitioned");

  vector<Point> points;

  const BBox& binBBox = _bins[0]->getBBox();

  vector<unsigned>::const_iterator cItr;
  itHGFEdgeLocal eItr;

  const HGraphWDimensions& nlhgraph = _capoPl.getNetlistHGraph();

  for (cItr = unionOfCells.begin(); cItr != unionOfCells.end(); ++cItr) {
    const HGFNode& node = nlhgraph.getNodeByIdx(*cItr);

    for (eItr = node.edgesBegin(); eItr != node.edgesEnd(); ++eItr) {
      unsigned eId = (*eItr)->getIndex();

      if (_edgesVisited.isBitSet(eId)) continue;

      _edgesVisited.setBit(eId);

      unsigned movablesInBin = combineFixedAndMovablePointSets(
          binBBox, _capoPl.getPointSet_fixed(eId),
          _capoPl.getPointSet_movable(eId), points);

      if (movablesInBin == 0) {
        cout << endl << "net Id " << eId << endl;
        cout << "bin id " << _bins[0]->getIndex() << endl;
        cout << "bin bbox " << binBBox << endl;
        cout << "fixed point set" << endl;
        for (unsigned i = 0; i < _capoPl.getPointSet_fixed(eId).size(); ++i) {
          cout << _capoPl.getPointSet_fixed(eId)[i].pt << " multiplicity "
               << _capoPl.getPointSet_fixed(eId)[i].count << endl;
        }
        cout << "movable point set" << endl;
        for (unsigned i = 0; i < _capoPl.getPointSet_movable(eId).size(); ++i) {
          cout << _capoPl.getPointSet_movable(eId)[i].pt << " multiplicity "
               << _capoPl.getPointSet_movable(eId)[i].count << endl;
        }
        cout << endl;
        cout << "net degree " << (*eItr)->getDegree() << endl;
        cout << "cells on the net" << endl;
        for (itHGFNodeLocal n = (*eItr)->nodesBegin(); n != (*eItr)->nodesEnd();
             ++n) {
          unsigned nID = (*n)->getIndex();

          cout << "Id " << nID << " \tplaced? "
               << (_capoPl.isPlaced(nID) ? "yes " : "no ") << "\tfixed? "
               << (_capoPl.getRBP().isFixed(nID) ? "yes " : "no ")
               << "\tpin pos "
               << _capoPl.getPinLocation(nID, eId).getGeomCenter()
               << "\tcell pos " << _capoPl.getPlacement()[nID];
          if (_capoPl.getCellToBinMap()[nID] == NULL)
            cout << "\tbin id NONE" << endl;
          else
            cout << "\tbin id " << _capoPl.getCellToBinMap()[nID]->getIndex()
                 << endl;
        }
      }

      abkfatal(movablesInBin > 0,
               "Why are we considering a net with no movables?");

      double oneSpan, twoSpan, threeSpan;

      bool useSteiner =
          (!_capoPl.getParams().useHPWL && (points.size() < 20)) ||
          _capoPl.getParams().useMST;

      if (useSteiner) {
        std::pair<vector<Point>::iterator, vector<Point>::iterator> p0Iters =
            std::equal_range(points.begin(), points.end(), p0CenterPoint);

        if (p0Iters.first == p0Iters.second) {
          unsigned position = p0Iters.first - points.begin();
          points.insert(p0Iters.first, p0CenterPoint);
          oneSpan = _capoPl.calculateSteinerFromPointSet(
              points,
#ifdef USEFLUTE
              _capoPl.getParams().useFastSt, _capoPl.getParams().useFLUTE,
              _capoPl.getParams().useMST);
#else
              _capoPl.getParams().useFastSt, _capoPl.getParams().useMST);
#endif
          points.erase(points.begin() + position);
        } else {
          oneSpan = _capoPl.calculateSteinerFromPointSet(
              points,
#ifdef USEFLUTE
              _capoPl.getParams().useFastSt, _capoPl.getParams().useFLUTE,
              _capoPl.getParams().useMST);
#else
              _capoPl.getParams().useFastSt, _capoPl.getParams().useMST);
#endif
        }

        std::pair<vector<Point>::iterator, vector<Point>::iterator> p1Iters =
            std::equal_range(points.begin(), points.end(), p1CenterPoint);

        if (p1Iters.first == p1Iters.second) {
          points.insert(p1Iters.first, p1CenterPoint);
          twoSpan = _capoPl.calculateSteinerFromPointSet(
              points,
#ifdef USEFLUTE
              _capoPl.getParams().useFastSt, _capoPl.getParams().useFLUTE,
              _capoPl.getParams().useMST);
#else
              _capoPl.getParams().useFastSt, _capoPl.getParams().useMST);
#endif
        } else {
          twoSpan = _capoPl.calculateSteinerFromPointSet(
              points,
#ifdef USEFLUTE
              _capoPl.getParams().useFastSt, _capoPl.getParams().useFLUTE,
              _capoPl.getParams().useMST);
#else
              _capoPl.getParams().useFastSt, _capoPl.getParams().useMST);
#endif
        }

        if (movablesInBin > 1) {
          // points already has p1CenterPoint
          p0Iters =
              std::equal_range(points.begin(), points.end(), p0CenterPoint);

          if (p0Iters.first == p0Iters.second) {
            points.insert(p0Iters.first, p0CenterPoint);
            threeSpan = _capoPl.calculateSteinerFromPointSet(
                points,
#ifdef USEFLUTE
                _capoPl.getParams().useFastSt, _capoPl.getParams().useFLUTE,
                _capoPl.getParams().useMST);
#else
                _capoPl.getParams().useFastSt, _capoPl.getParams().useMST);
#endif
          } else {
            threeSpan = _capoPl.calculateSteinerFromPointSet(
                points,
#ifdef USEFLUTE
                _capoPl.getParams().useFastSt, _capoPl.getParams().useFLUTE,
                _capoPl.getParams().useMST);
#else
                _capoPl.getParams().useFastSt, _capoPl.getParams().useMST);
#endif
          }
        } else {
          threeSpan = std::max(oneSpan, twoSpan);
        }

        // fidelity checking
        threeSpan = min(threeSpan, oneSpan + unitSpan);
        threeSpan = min(threeSpan, twoSpan + unitSpan);
        oneSpan = min(oneSpan, threeSpan);
        twoSpan = min(twoSpan, threeSpan);

        // averaging of three numbers
        double two_min = twoSpan * 0.99;
        double two_max = twoSpan * 1.01;
        double three_min = threeSpan * 0.99;
        double three_max = threeSpan * 1.01;

        if ((oneSpan > two_min) && (oneSpan < two_max)) {
          if ((oneSpan > three_min) && (oneSpan < three_max)) {
            double avg = (oneSpan + twoSpan + threeSpan) / 3;
            oneSpan = avg;
            twoSpan = avg;
            threeSpan = avg;
          } else {
            double avg = (oneSpan + twoSpan) * 0.5;
            oneSpan = avg;
            twoSpan = avg;
          }
        } else if ((oneSpan > three_min) && (oneSpan < three_max)) {
          double avg = (oneSpan + threeSpan) * 0.5;
          oneSpan = avg;
          threeSpan = avg;
        } else if ((twoSpan > three_min) && (twoSpan < three_max)) {
          double avg = (twoSpan + threeSpan) * 0.5;
          twoSpan = avg;
          threeSpan = avg;
        }
        // end averaging
      } else {
        // using HPWL to figure out the three spans
        calculateSpansFromHPWL(points, p0CenterPoint, p1CenterPoint,
                               movablesInBin, oneSpan, twoSpan, threeSpan);
      }

      _weights[eId][0] += oneSpan;
      _weights[eId][1] += twoSpan;
      _weights[eId][2] += threeSpan;
      // end unified theto
    }
  }
}

void PartitioningProblemForCapo::constructHG(bool honourGrpConstr) {
  _edgesVisited.clear();

  double p0Boundry, p1Boundry, p1Left, p0Right, p0Bottom, p1Top;
  double p0Center = 0.0, p1Center = 0.0;
  const BBox& p0Box = (*_partitions)[0];
  const BBox& p1Box = (*_partitions)[1];

  if (_capoPl.getParams().splitterParams.newFuzzy &&
      (!_capoPl.getParams().splitterParams.newFuzzyOnlyRepart || _isRepart)) {
    double totalArea = (*_capacities)[0][0] + (*_capacities)[1][0];
    double p0MinFraction = (*_minCapacities)[0][0] / totalArea;
    double p1MinFraction = (*_minCapacities)[1][0] / totalArea;
    if (_nDim == 2)  // horizontal bisection
    {
      double ySpace = max(p0Box.yMin - p1Box.yMax, 0.0);
      p0Bottom = p0Box.yMin - 0.5 * ySpace;
      p1Top = p1Box.yMax + 0.5 * ySpace;
      p0Center = 0.5 * (p0Box.yMax - p0Box.yMin) + p0Box.yMin;
      p1Center = 0.5 * (p1Box.yMax - p1Box.yMin) + p1Box.yMin;
      double ySpan = p0Box.yMax - p1Box.yMin;
      p0Boundry = p0Box.yMax - p0MinFraction * ySpan;
      p1Boundry = p1Box.yMin + p1MinFraction * ySpan;
    } else {
      double xSpace = max(p1Box.xMin - p0Box.xMax, 0.0);
      p0Right = p0Box.xMax + 0.5 * xSpace;
      p1Left = p1Box.xMin - 0.5 * xSpace;
      p0Center = 0.5 * (p0Box.xMax - p0Box.xMin) + p0Box.xMin;
      p1Center = 0.5 * (p1Box.xMax - p1Box.xMin) + p1Box.xMin;
      double xSpan = p1Box.xMax - p0Box.xMin;
      p0Boundry = p0Box.xMin + p0MinFraction * xSpan;
      p1Boundry = p1Box.xMax - p1MinFraction * xSpan;
    }
  } else {
    const double fuzzyFactor = 0.3;
    if (_nDim == 2)  // horizontal bisection
    {
      // vert. space between the bins
      double ySpace = max(p0Box.yMin - p1Box.yMax, 0.0);
      p0Bottom = p0Box.yMin - 0.5 * ySpace;
      p1Top = p1Box.yMax + 0.5 * ySpace;
      p0Center = 0.5 * (p0Box.yMax - p0Box.yMin) + p0Box.yMin;
      p1Center = 0.5 * (p1Box.yMax - p1Box.yMin) + p1Box.yMin;

      // lowest y-loc that propagates into p0
      p0Boundry = p0Bottom + (p0Box.yMax - p0Bottom) * fuzzyFactor;
      // highest y-loc that propagates into p1
      p1Boundry = p1Top - (p1Top - p1Box.yMin) * fuzzyFactor;
    } else {
      double xSpace = max(p1Box.xMin - p0Box.xMax, 0.0);
      p0Right = p0Box.xMax + 0.5 * xSpace;
      p1Left = p1Box.xMin - 0.5 * xSpace;
      p0Center = 0.5 * (p0Box.xMax - p0Box.xMin) + p0Box.xMin;
      p1Center = 0.5 * (p1Box.xMax - p1Box.xMin) + p1Box.xMin;
      p0Boundry = p0Right - (p0Right - p0Box.xMin) * fuzzyFactor;
      p1Boundry = p1Left + (p1Box.xMax - p1Left) * fuzzyFactor;
    }
  }

  const vector<unsigned>& unionOfCells = _bins[0]->getCellIds();

  abkfatal(_bins.size() == 1, "There should only be 1 bin to be partitioned");

  bool splitting_whole_netlist =
      (unionOfCells.size() == (_capoPl.getNetlistHGraph().getNumNodes() -
                               _capoPl.getNetlistHGraph().getNumTerminals()));

  double terminal_fraction =
      static_cast<double>(_capoPl.getNetlistHGraph().getNumTerminals()) /
      static_cast<double>(_capoPl.getNetlistHGraph().getNumNodes());

  if (splitting_whole_netlist &&
      (terminal_fraction <
       _capoPl.getParams().splitterParams.termMergeThreshold)) {
    if (_params.verb.getForActions())
      cout << "Splitting the entire netlist, saving memory by not building a "
              "new subgraph"
           << endl;

    _ownsHGraphs = false;

    _hgraph = static_cast<HGraphFixed*>(
        const_cast<HGraphWDimensions*>(&_capoPl.getNetlistHGraph()));
    _subHG = static_cast<SubHGraph*>(
        const_cast<HGraphWDimensions*>(&_capoPl.getNetlistHGraph()));

    _hgraph->temporarilyZeroOutTermWeights();

    _terminalToBlock =
        new std::vector<unsigned>(_hgraph->getNumTerminals(), UINT_MAX);
    _padBlocks = new std::vector<BBox>(_hgraph->getNumTerminals());
    PartitionIds movable;
    movable.setToAll(2);
    _fixedConstr = new Partitioning(_hgraph->getNumNodes(), movable);

    double cellWidth, cellHeight;
    for (unsigned i = 0; i < _hgraph->getNumTerminals(); ++i) {
      (*_terminalToBlock)[i] = i;

      unsigned angle = _capoPl.getPlacement().getOrient(i).getAngle();
      bool notRotated = angle == 0 || angle == 180;
      cellWidth = notRotated ? _capoPl.getNetlistHGraph().getNodeWidth(i)
                             : _capoPl.getNetlistHGraph().getNodeHeight(i);
      cellHeight = notRotated ? _capoPl.getNetlistHGraph().getNodeHeight(i)
                              : _capoPl.getNetlistHGraph().getNodeWidth(i);

      (*_padBlocks)[i] += _capoPl.getPlacement()[i];  // bottomleft
      (*_padBlocks)[i] += (_capoPl.getPlacement()[i] +
                           Point(cellWidth, cellHeight));  // topright

      // figure if this terminal goes in part 0 or 1
      const HGFNode& n = _hgraph->getNodeByIdx(i);
      const BBox& termLoc =
          (n.getDegree() > 0)
              ? _capoPl.getPinLocation(i, (*(n.edgesBegin()))->getIndex())
              : _capoPl.getPlacement()[i];
      if (_nDim == 2) {
        if (termLoc.yMin <= p1Boundry && termLoc.yMax >= p0Boundry) {
        } else if (termLoc.yMax >= p0Boundry)
          (*_fixedConstr)[i].setToPart(0);
        else if (termLoc.yMin <= p1Boundry)
          (*_fixedConstr)[i].setToPart(1);
      } else {
        if (termLoc.xMin <= p0Boundry && termLoc.xMax >= p1Boundry) {
        } else if (termLoc.xMin <= p0Boundry)
          (*_fixedConstr)[i].setToPart(0);
        else if (termLoc.xMax >= p1Boundry)
          (*_fixedConstr)[i].setToPart(1);
      }
    }
  } else {
    _ownsHGraphs = true;
    SubHGraphForCapo* newHG = new SubHGraphForCapo(2, unionOfCells.size());
    _hgraph = static_cast<HGraphFixed*>(newHG);
    _subHG = static_cast<SubHGraph*>(newHG);

    vector<unsigned>::const_iterator cItr;
    itHGFEdgeLocal eItr;

    const HGraphWDimensions& nlhgraph = _capoPl.getNetlistHGraph();

    newHG->_param.makeAllSrcSnk = true;
    unsigned newId = 2;
    for (cItr = unionOfCells.begin(); cItr != unionOfCells.end();
         ++cItr, ++newId) {
      newHG->addNode(nlhgraph, nlhgraph.getNodeByIdx(*cItr), newId);
    }

    const CapoBin& bin = *_bins[0];
    double relWS = bin.getRelativeWhitespace();
    double wsAdjustUtilization = max(0.05, 1 - relWS);
    double mult_factor = 15;
    double unitSpan = fabs(p0Center - p1Center);

    vector<unsigned> nonTermsOnNet;  // ids of all nonTerminals on this edge

    for (cItr = unionOfCells.begin(); cItr != unionOfCells.end(); ++cItr) {
      const HGFNode& node = nlhgraph.getNodeByIdx(*cItr);
      for (eItr = node.edgesBegin(); eItr != node.edgesEnd(); ++eItr) {
        unsigned eId = (*eItr)->getIndex();

        if (_edgesVisited.isBitSet(eId)) continue;

        _edgesVisited.setBit(eId);

        double oneSpan = _weights[eId][0];
        double twoSpan = _weights[eId][1];
        double threeSpan = _weights[eId][2];

        nonTermsOnNet.clear();
        for (itHGFNodeLocal adjNd = (*eItr)->nodesBegin();
             adjNd != (*eItr)->nodesEnd(); ++adjNd) {
          unsigned adjIdx = (*adjNd)->getIndex();
          if (_cellToBinMap[adjIdx] == _bins[0]) {
            nonTermsOnNet.push_back(adjIdx);
          }
        }

        // check for inessential net
        if (equalDouble(oneSpan, twoSpan) && equalDouble(twoSpan, threeSpan))
          continue;

        if (equalDouble(oneSpan, 0.) && equalDouble(twoSpan, 0.)) {
          // has no terminals
          double net_weight = 1.;
          bool t0 = false, t1 = false;
          addEdgePropagateTerminalsReweight(newHG, eItr, nonTermsOnNet, t0, t1,
                                            net_weight * mult_factor);
          continue;
        }

        if (lessThanDouble(oneSpan, twoSpan)) {
          double net_weight = (twoSpan - oneSpan) / unitSpan;
          if (equalDouble(twoSpan, threeSpan) &&
              _capoPl.getParams().WSadjustment) {
            // net is terminal bound from the p0 side
            double factor = 1.;
            factor /= pow(wsAdjustUtilization,
                          _capoPl.getParams().termWeightModifier);
            factor = std::min(factor, 3.);
            net_weight *= factor;
          }
          bool t0 = true, t1 = false;
          addEdgePropagateTerminalsReweight(newHG, eItr, nonTermsOnNet, t0, t1,
                                            net_weight * mult_factor);
        } else if (greaterThanDouble(oneSpan, twoSpan)) {
          double net_weight = (oneSpan - twoSpan) / unitSpan;
          if (equalDouble(oneSpan, threeSpan) &&
              _capoPl.getParams().WSadjustment) {
            // net is terminal bound from the p1 side
            double factor = 1.;
            factor /= pow(wsAdjustUtilization,
                          _capoPl.getParams().termWeightModifier);
            factor = std::min(factor, 3.);
            net_weight *= factor;
          }
          bool t0 = false, t1 = true;
          addEdgePropagateTerminalsReweight(newHG, eItr, nonTermsOnNet, t0, t1,
                                            net_weight * mult_factor);
        }

        if (greaterThanDouble(threeSpan, std::max(oneSpan, twoSpan))) {
          double net_weight =
              (threeSpan - std::max(oneSpan, twoSpan)) / unitSpan;
          bool t0 = false, t1 = false;
          addEdgePropagateTerminalsReweight(newHG, eItr, nonTermsOnNet, t0, t1,
                                            net_weight * mult_factor);
        }
      }
    }

    _capoPl.getParams().maxMem->update(
        "PartProbForCapo construction before finalize");

    // destroy the original hypergraph if necessary
    bool destroyNetlistHGraph =
        _size == Large &&
        static_cast<double>(unionOfCells.size()) >=
            _capoPl.getParams().keepHGraphFraction *
                static_cast<double>(_capoPl.getNetlistHGraph().getNumNodes());

    if (destroyNetlistHGraph) {
      _capoPl.destroyHGraph();
    }

    newHG->finalize();

    _capoPl.getParams().maxMem->update(
        "PartProbForCapo construction after finalize");

    _terminalToBlock = new std::vector<unsigned>(2, UINT_MAX);

    (*_terminalToBlock)[0] = 0;
    (*_terminalToBlock)[1] = 1;

    _padBlocks = new std::vector<BBox>(2);

    (*_padBlocks)[0] += Point((*_partitions)[0].xMin, (*_partitions)[0].yMax);
    (*_padBlocks)[1] += Point((*_partitions)[1].xMax, (*_partitions)[1].yMin);

    PartitionIds movable;
    movable.setToAll(2);
    _fixedConstr = new Partitioning(_hgraph->getNumNodes(), movable);
    (*_fixedConstr)[0].setToPart(0);
    (*_fixedConstr)[1].setToPart(1);
  }

  if (honourGrpConstr && _capoPl.groupRegionConstr.size() > 0) {
    vector<unsigned> fixToP0;
    vector<unsigned> fixToP1;

    const BBox& p0Box = (*_partitions)[0];
    const BBox& p1Box = (*_partitions)[1];
    BBox intersectionP0;
    BBox intersectionP1;
    bool warnedOnce = false;
    bool warnedOnceMore = false;
    const RBPlacement& RBP = _capoPl.getRBP();

    vector<double> P0wGrpsIntAreas(_capoPl.groupRegionConstr.size(), -1.);
    vector<double> P1wGrpsIntAreas(_capoPl.groupRegionConstr.size(), -1.);

    double areaBlkBox = RBP.getContainedSiteAreaInBBox(p0Box) +
                        RBP.getContainedSiteAreaInBBox(p1Box);

    double area0, area1;
    for (unsigned i = 0; i < unionOfCells.size(); ++i) {
      unsigned cellId = unionOfCells[i];
      unsigned groupId = _capoPl.cellToGrpMapping[cellId];

      if (groupId != UINT_MAX) {
        const BBox& rgnBBox = _capoPl.groupRegionConstr[groupId].first;
        if (_nDim > 1)  // horizontal split
        {
          if (greaterOrEqualDouble(rgnBBox.yMin, p0Box.yMin))
            fixToP0.push_back(i + 2);
          else if (lessOrEqualDouble(rgnBBox.yMax, p1Box.yMax))
            fixToP1.push_back(i + 2);
          else  // bbox intersects
          {
            intersectionP0 = p0Box.intersect(rgnBBox);
            intersectionP1 = p1Box.intersect(rgnBBox);
            if (equalDouble(P0wGrpsIntAreas[groupId], -1.)) {
              area0 = RBP.getContainedSiteAreaInBBox(intersectionP0);
              P0wGrpsIntAreas[groupId] = area0;
            } else {
              area0 = P0wGrpsIntAreas[groupId];
            }
            if (equalDouble(P1wGrpsIntAreas[groupId], -1.)) {
              area1 = RBP.getContainedSiteAreaInBBox(intersectionP1);
              P1wGrpsIntAreas[groupId] = area1;
            } else {
              area1 = P1wGrpsIntAreas[groupId];
            }

            double area = area0 + area1;

            if ((area / areaBlkBox) > 0.99) {
              if (!warnedOnceMore && _verb.getForActions() > 2) {
                abkwarn(0, "Region w grp constraints completely in pl. bin ");
                cout << "Rgn BBox " << rgnBBox << " Pl. bin "
                     << _bins[0]->getBBox() << endl;
                warnedOnceMore = true;
              }
            } else if ((area0 / area) > 0.7)
              fixToP0.push_back(i + 2);
            else if ((area1 / area) > 0.7)
              fixToP1.push_back(i + 2);
            else if (!warnedOnce && _verb.getForActions() > 2) {
              abkwarn(0, "Cutting a group with regionConstraints \n");
              cout << "Rgn BBox " << rgnBBox << endl;
              warnedOnce = true;
            }
          }
        } else  // vertical split
        {
          if (lessOrEqualDouble(rgnBBox.xMax, p0Box.xMax))
            fixToP0.push_back(i + 2);
          else if (greaterOrEqualDouble(rgnBBox.xMin, p1Box.xMin))
            fixToP1.push_back(i + 2);
          else  // bbox intersects
          {
            intersectionP0 = p0Box.intersect(rgnBBox);
            intersectionP1 = p1Box.intersect(rgnBBox);
            if (equalDouble(P0wGrpsIntAreas[groupId], -1.)) {
              area0 = RBP.getContainedSiteAreaInBBox(intersectionP0);
              P0wGrpsIntAreas[groupId] = area0;
            } else {
              area0 = P0wGrpsIntAreas[groupId];
            }
            if (equalDouble(P1wGrpsIntAreas[groupId], -1.)) {
              area1 = RBP.getContainedSiteAreaInBBox(intersectionP1);
              P1wGrpsIntAreas[groupId] = area1;
            } else {
              area1 = P1wGrpsIntAreas[groupId];
            }

            double area = area0 + area1;
            if ((area / areaBlkBox) > 0.99) {
              if (!warnedOnceMore && _verb.getForActions() > 2) {
                abkwarn(
                    0, "Region w grp constraints completely in placement bin ");
                cout << "Rgn BBox " << rgnBBox << " Pl. bin "
                     << _bins[0]->getBBox() << endl;
                warnedOnceMore = true;
              }
            } else if ((area0 / area) > 0.7)
              fixToP0.push_back(i + 2);
            else if ((area1 / area) > 0.7)
              fixToP1.push_back(i + 2);
            else if (!warnedOnce && _verb.getForActions() > 2) {
              abkwarn(0, "Cutting a group with regionConstraints \n");
              cout << "Rgn BBox " << rgnBBox << endl;
              warnedOnce = true;
            }
          }
        }
      }
    }

    if (_verb.getForActions() > 1) {
      if (_nDim > 1)
        cout << "Horizontal Split " << endl;
      else
        cout << "Vertical Split " << endl;

      cout << "Num Nodes fixed to Partition 0: " << fixToP0.size() << endl;
      cout << "Num Nodes fixed to Partition 1: " << fixToP1.size() << endl;
    }

    for (unsigned i = 0; i < fixToP0.size(); ++i)
      (*_fixedConstr)[fixToP0[i]].setToPart(0);
    for (unsigned i = 0; i < fixToP1.size(); ++i)
      (*_fixedConstr)[fixToP1[i]].setToPart(1);
  }
}

void PartitioningProblemForCapo::calculateCapacities_LargeH() {
  const CapoBin& bin = *_bins[0];

  double capacity = bin.getCapacity();
  bool changeDefaultFlow = false;

  if (_capoPl.getParams().splitterParams.useRgnConstr) {
    vector<double> partCaps(2);
    changeDefaultFlow = hPartCapsFrmRgnConstr(_splitRow, partCaps);
    if (changeDefaultFlow) {
      double p0Target = partCaps[0] * bin.getTotalCellArea();
      double p1Target = bin.getTotalCellArea() - p0Target;
      setCapacity(0, p0Target);
      setCapacity(1, p1Target);
    } else {
      double capRatio = _p0Capacity / capacity;
      double p0CellArea = capRatio * bin.getTotalCellArea();
      setCapacity(0, p0CellArea);
      setCapacity(1, bin.getTotalCellArea() - p0CellArea);
    }
  } else {
    double capRatio = _p0Capacity / capacity;
    double p0CellArea = capRatio * bin.getTotalCellArea();
    setCapacity(0, p0CellArea);
    setCapacity(1, bin.getTotalCellArea() - p0CellArea);
  }

  _splitPtFixed =
      _splitPtFixed ||
      (_capoPl.getParams().splitterParams.useRgnConstr && changeDefaultFlow);
}

void PartitioningProblemForCapo::calculateCapacities_LargeH_repart() {
  const CapoBin& bin = *_bins[0];

  if (_capoPl.getParams().splitterParams.useRgnConstr) {
    vector<double> partCaps(2);
    bool changeDefaultFlow = hPartCapsFrmRgnConstr(_splitRow, partCaps);
    if (changeDefaultFlow) {
      double p0Target = partCaps[0] * bin.getTotalCellArea();
      double p1Target = bin.getTotalCellArea() - p0Target;
      setCapacity(0, p0Target);
      setCapacity(1, p1Target);
    } else {
      double capRatio = _p0Capacity / bin.getCapacity();
      double p0CellArea = capRatio * bin.getTotalCellArea();
      setCapacity(0, p0CellArea);
      setCapacity(1, bin.getTotalCellArea() - p0CellArea);
    }
  } else {
    double p0Target =
        (_p0Capacity / bin.getCapacity()) * bin.getTotalCellArea();
    setCapacity(0, p0Target);
    setCapacity(1, bin.getTotalCellArea() - p0Target);
  }
}

void PartitioningProblemForCapo::calculateTolerance_LargeH() {
  const CapoSplitterParams& params = _capoPl.getParams().splitterParams;
  const CapoBin& bin = *_bins[0];
  const double& binWS = bin.getRelativeWhitespace();
  _repartSmallWS = params.repartSmallWS && binWS < params.smallWS;
  _changedTol = false;

  if (_splitPtFixed) {
    if (greaterOrEqualDouble(params.safeWS, 0.) && binWS >= params.safeWS) {
      setSafeWSToleranceH(_splitRow, params.safeWS);
    } else if (greaterOrEqualDouble(params.minLocalWS, 0.) &&
               binWS >= params.minLocalWS) {
      setMinLocalWSToleranceH(_splitRow, params.minLocalWS);
    } else {
      setTolerance(0.5 * params.defaultRepartTolerance);
    }
  } else if (params.doRepartitioning || _repartSmallWS ||
             !params.useWSTolMethod) {
    setTolerance(params.defaultTolerance);
  } else {
    if (greaterOrEqualDouble(params.safeWS, 0.) && binWS >= params.safeWS) {
      setSafeWSToleranceH(_splitRow, params.safeWS);
    } else if (greaterOrEqualDouble(params.minLocalWS, 0.) &&
               binWS >= params.minLocalWS) {
      setMinLocalWSToleranceH(_splitRow, params.minLocalWS);
    } else {
      setPaperMethodToleranceH(_splitRow);
      if (getTolerance() > params.defaultTolerance) {
        if (params.verb.getForActions() > 4) {
          cout << "Initial Tolerance " << getTolerance() * 100. << endl;
        }
        decreaseTolerance(params.defaultTolerance);
        if (params.verb.getForActions() > 4) {
          cout << "New Tolerance " << getTolerance() * 100. << endl;
        }
        //_changedTol = true;
      }
    }
    if (getTolerance() < params.minTolerance) {
      if (params.verb.getForActions() > 4) {
        cout << "Initial Tolerance " << getTolerance() * 100. << endl;
      }
      increaseTolerance(params.minTolerance);
      if (params.verb.getForActions() > 4) {
        cout << "New Tolerance " << getTolerance() * 100. << endl;
      }
      _changedTol = true;
    }
  }
}

void PartitioningProblemForCapo::calculateTolerance_LargeH_repart() {
  const CapoSplitterParams& params = _capoPl.getParams().splitterParams;
  const CapoBin& bin = *_bins[0];
  double binWS = bin.getRelativeWhitespace();
  double smallestTol = bin.getMedianCellArea() / bin.getTotalCellArea();
  _repartSmallWS = params.repartSmallWS && binWS < params.smallWS;

  if (_repartSmallWS) {
    double medianCellSize = 0.5 * bin.getMedianCellArea();
    const std::vector<std::vector<double> >& caps = getCapacities();
    double minTol =
        max(medianCellSize / caps[0][0], medianCellSize / caps[1][0]);
    double tol = max(0., binWS);
    tol = max(minTol, tol);
    setTolerance(tol);
  } else if (_splitPtFixed) {
    if (greaterOrEqualDouble(params.safeWS, 0.) && binWS >= params.safeWS) {
      setSafeWSToleranceH(_splitRow, params.safeWS);
    } else if (greaterOrEqualDouble(params.minLocalWS, 0.) &&
               binWS >= params.minLocalWS) {
      setMinLocalWSToleranceH(_splitRow, params.minLocalWS);
    } else {
      setTolerance(0.25 * params.defaultRepartTolerance);
    }
  } else if (params.useWSTolMethod) {
    if (greaterOrEqualDouble(params.safeWS, 0.) && binWS >= params.safeWS) {
      setSafeWSToleranceH(_splitRow, params.safeWS);
    } else if (greaterOrEqualDouble(params.minLocalWS, 0.) &&
               binWS >= params.minLocalWS) {
      setMinLocalWSToleranceH(_splitRow, params.minLocalWS);
    } else {
      setPaperMethodToleranceH(_splitRow);
      if (getTolerance() > params.defaultRepartTolerance) {
        decreaseTolerance(params.defaultRepartTolerance);
      }
    }
  } else {
    setTolerance(params.defaultRepartTolerance);
  }

  if (getTolerance() < smallestTol) {
    increaseTolerance(smallestTol);
  }
}

// largeV
void PartitioningProblemForCapo::calculateCapacities_LargeV() {
  const CapoBin& bin = *_bins[0];
  double capacity = bin.getCapacity();
  bool changeDefaultFlow = false;

  if (_capoPl.getParams().splitterParams.useRgnConstr) {
    vector<double> partCaps(2);
    changeDefaultFlow = vPartCapsFrmRgnConstr(_xSplit, partCaps);
    if (changeDefaultFlow) {
      double p0Target = partCaps[0] * bin.getTotalCellArea();
      double p1Target = bin.getTotalCellArea() - p0Target;
      setCapacity(0, p0Target);
      setCapacity(1, p1Target);
    } else {
      double capRatio = _p0Capacity / capacity;
      double p0CellArea = capRatio * bin.getTotalCellArea();
      setCapacity(0, p0CellArea);
      setCapacity(1, bin.getTotalCellArea() - p0CellArea);
    }
  } else {
    double capRatio = _p0Capacity / capacity;
    double p0CellArea = capRatio * bin.getTotalCellArea();
    setCapacity(0, p0CellArea);
    setCapacity(1, bin.getTotalCellArea() - p0CellArea);
  }

  _splitPtFixed =
      _splitPtFixed ||
      (_capoPl.getParams().splitterParams.useRgnConstr && changeDefaultFlow);
}

void PartitioningProblemForCapo::calculateCapacities_LargeV_repart() {
  const CapoBin& bin = *_bins[0];

  if (_capoPl.getParams().splitterParams.useRgnConstr) {
    vector<double> partCaps(2);
    bool changeDefaultFlow = vPartCapsFrmRgnConstr(_xSplit, partCaps);
    if (changeDefaultFlow) {
      double p0Target = partCaps[0] * bin.getTotalCellArea();
      double p1Target = bin.getTotalCellArea() - p0Target;
      setCapacity(0, p0Target);
      setCapacity(1, p1Target);
    } else {
      double capRatio = _p0Capacity / bin.getCapacity();
      double p0CellArea = capRatio * bin.getTotalCellArea();
      setCapacity(0, p0CellArea);
      setCapacity(1, bin.getTotalCellArea() - p0CellArea);
    }
  } else {
    double p0Target =
        (_p0Capacity / bin.getCapacity()) * bin.getTotalCellArea();
    setCapacity(0, p0Target);
    setCapacity(1, bin.getTotalCellArea() - p0Target);
  }
}

void PartitioningProblemForCapo::calculateTolerance_LargeV() {
  const CapoSplitterParams& params = _capoPl.getParams().splitterParams;
  const CapoBin& bin = *_bins[0];
  double binWS = bin.getRelativeWhitespace();
  _repartSmallWS = params.repartSmallWS && binWS < params.smallWS;
  bool _uniformWS = lessThanDouble(params.safeWS, 0.) &&
                    lessThanDouble(params.minLocalWS, 0.);

  if (_splitPtFixed) {
    if (greaterOrEqualDouble(params.safeWS, 0.) && binWS >= params.safeWS) {
      setSafeWSToleranceV(_xSplit, params.safeWS);
    } else if (greaterOrEqualDouble(params.minLocalWS, 0.) &&
               binWS >= params.minLocalWS) {
      setMinLocalWSToleranceV(_xSplit, params.minLocalWS);
    } else {
      setTolerance(0.5 * params.defaultTolerance);
    }
  } else if (params.useWSTolMethod && !_uniformWS) {
    if (greaterOrEqualDouble(params.safeWS, 0.) && binWS >= params.safeWS) {
      setSafeWSToleranceV(_xSplit, params.safeWS);
    } else if (greaterOrEqualDouble(params.minLocalWS, 0.) &&
               binWS >= params.minLocalWS) {
      setMinLocalWSToleranceV(_xSplit, params.minLocalWS);
    } else {
      setPaperMethodToleranceV(_xSplit);
    }
    if (getTolerance() < params.minTolerance) {
      increaseTolerance(params.minTolerance);
    }
  } else {
    setTolerance(params.defaultTolerance);
  }
}

void PartitioningProblemForCapo::calculateTolerance_LargeV_repart() {
  const CapoSplitterParams& params = _capoPl.getParams().splitterParams;
  const CapoBin& bin = *_bins[0];
  double binWS = bin.getRelativeWhitespace();
  double smallestTol = bin.getMedianCellArea() / bin.getTotalCellArea();
  _repartSmallWS = params.repartSmallWS && binWS < params.smallWS;

  if (_repartSmallWS) {
    double medianCellSize = 0.5 * bin.getMedianCellArea();
    const std::vector<std::vector<double> >& caps = getCapacities();
    double minTol =
        max(medianCellSize / caps[0][0], medianCellSize / caps[1][0]);
    double tol = max(0., binWS);
    tol = max(minTol, tol);
    setTolerance(tol);
  } else if (_splitPtFixed) {
    if (greaterOrEqualDouble(params.safeWS, 0.) && binWS >= params.safeWS) {
      setSafeWSToleranceV(_xSplit, params.safeWS);
    } else if (greaterOrEqualDouble(params.minLocalWS, 0.) &&
               binWS >= params.minLocalWS) {
      setMinLocalWSToleranceV(_xSplit, params.minLocalWS);
    } else {
      setTolerance(0.25 * params.defaultRepartTolerance);
    }
  } else if (params.useWSTolMethod) {
    if (greaterOrEqualDouble(params.safeWS, 0.) && binWS >= params.safeWS) {
      setSafeWSToleranceV(_xSplit, params.safeWS);
    } else if (greaterOrEqualDouble(params.minLocalWS, 0.) &&
               binWS >= params.minLocalWS) {
      setMinLocalWSToleranceV(_xSplit, params.minLocalWS);
    } else {
      setPaperMethodToleranceV(_xSplit);
      if (getTolerance() > params.defaultRepartTolerance) {
        decreaseTolerance(params.defaultRepartTolerance);
      }
    }
  } else {
    setTolerance(params.defaultRepartTolerance);
  }

  if (getTolerance() < smallestTol) {
    increaseTolerance(smallestTol);
  }
}

// split smallH
void PartitioningProblemForCapo::calculateCapacities_SmallH() {
  const CapoBin& bin = *_bins[0];
  double capRatio = _p0Capacity / bin.getCapacity();
  double p0CellArea = capRatio * bin.getTotalCellArea();
  setCapacity(0, p0CellArea);
  setCapacity(1, bin.getTotalCellArea() - p0CellArea);
}

void PartitioningProblemForCapo::calculateTolerance_SmallH(bool midSize) {
  const CapoSplitterParams& params = _capoPl.getParams().splitterParams;
  const CapoBin& bin = *_bins[0];
  double binWS = bin.getRelativeWhitespace();
  _changedTol = false;
  _repartSmallWS = midSize && params.repartSmallWS && binWS < params.smallWS;

  if (params.doRepartitioning || params.useWSTolMethod || _repartSmallWS) {
    if (greaterOrEqualDouble(params.safeWS, 0.) && binWS >= params.safeWS) {
      setSafeWSToleranceH(_splitRow, params.safeWS);
    } else if (greaterOrEqualDouble(params.minLocalWS, 0.) &&
               binWS >= params.minLocalWS) {
      setMinLocalWSToleranceH(_splitRow, params.minLocalWS);
    } else {
      setPaperMethodToleranceH(_splitRow);
      if (getTolerance() > params.defaultTolerance) {
        decreaseTolerance(params.defaultTolerance);
      }
    }
    if (midSize && getTolerance() < params.minTolerance) {
      increaseTolerance(params.minTolerance);
      _changedTol = true;
    }
  } else {
    setTolerance(params.defaultTolerance);
  }
}

void PartitioningProblemForCapo::calculateTolerance_SmallH_repart(
    bool midSize) {
  const CapoSplitterParams& params = _capoPl.getParams().splitterParams;
  const CapoBin& bin = *_bins[0];
  double binWS = bin.getRelativeWhitespace();
  _repartSmallWS = midSize && params.repartSmallWS && binWS < params.smallWS;

  if (midSize && _repartSmallWS) {
    double medianCellSize = bin.getMedianCellArea();
    const std::vector<std::vector<double> >& caps = getCapacities();
    double minTol =
        max(medianCellSize / caps[0][0], medianCellSize / caps[1][0]);
    double tol = max(0., binWS);
    tol = max(minTol, tol);
    setTolerance(tol);
  } else if (params.useWSTolMethod) {
    if (greaterOrEqualDouble(params.safeWS, 0.) && binWS >= params.safeWS) {
      setSafeWSToleranceH(_splitRow, params.safeWS);
    } else if (greaterOrEqualDouble(params.minLocalWS, 0.) &&
               binWS >= params.minLocalWS) {
      setMinLocalWSToleranceH(_splitRow, params.minLocalWS);
    } else {
      setPaperMethodToleranceH(_splitRow);
      if (getTolerance() > params.defaultRepartTolerance) {
        decreaseTolerance(params.defaultRepartTolerance);
      }
    }
  } else {
    setTolerance(params.defaultRepartTolerance);
  }
}

// split smallV
void PartitioningProblemForCapo::calculateCapacities_SmallV() {
  const CapoBin& bin = *_bins[0];
  double p0Target = (_p0Capacity / bin.getCapacity()) * bin.getTotalCellArea();
  double p1Target = bin.getTotalCellArea() - p0Target;
  setCapacity(0, p0Target);
  setCapacity(1, p1Target);
}

void PartitioningProblemForCapo::calculateTolerance_SmallV(bool midSize) {
  const CapoSplitterParams& params = _capoPl.getParams().splitterParams;
  const CapoBin& bin = *_bins[0];
  double binWS = bin.getRelativeWhitespace();
  _repartSmallWS = midSize && params.repartSmallWS && binWS < params.smallWS;
  bool _uniformWS = lessThanDouble(params.safeWS, 0.) &&
                    lessThanDouble(params.minLocalWS, 0.);

  if (params.useWSTolMethod && !_uniformWS) {
    if (greaterOrEqualDouble(params.safeWS, 0.) && binWS >= params.safeWS) {
      setSafeWSToleranceV(_xSplit, params.safeWS);
    } else if (greaterOrEqualDouble(params.minLocalWS, 0.) &&
               binWS >= params.minLocalWS) {
      setMinLocalWSToleranceV(_xSplit, params.minLocalWS);
    } else {
      setPaperMethodToleranceV(_xSplit);
    }
    if (getTolerance() < params.minTolerance) {
      increaseTolerance(params.minTolerance);
    }
  } else {
    setTolerance(params.defaultTolerance);
  }
}

void PartitioningProblemForCapo::calculateTolerance_SmallV_repart(
    bool midSize) {
  const CapoSplitterParams& params = _capoPl.getParams().splitterParams;
  const CapoBin& bin = *_bins[0];
  double binWS = bin.getRelativeWhitespace();
  _repartSmallWS = midSize && params.repartSmallWS && binWS < params.smallWS;
  bool _uniformWS = lessThanDouble(params.safeWS, 0.) &&
                    lessThanDouble(params.minLocalWS, 0.);

  if (midSize && _repartSmallWS) {
    double medianCellSize = 0.5 * bin.getMedianCellArea();
    // double medianCellSize =  bin.getLargestCellArea();
    const std::vector<std::vector<double> >& caps = getCapacities();
    double minTol =
        max(medianCellSize / caps[0][0], medianCellSize / caps[1][0]);
    double tol = max(0., binWS);
    tol = max(minTol, tol);
    setTolerance(tol);
  } else if (params.useWSTolMethod && !_uniformWS) {
    if (greaterOrEqualDouble(params.safeWS, 0.) && binWS >= params.safeWS) {
      setSafeWSToleranceV(_xSplit, params.safeWS);
    } else if (greaterOrEqualDouble(params.minLocalWS, 0.) &&
               binWS >= params.minLocalWS) {
      setMinLocalWSToleranceV(_xSplit, params.minLocalWS);
    } else {
      setPaperMethodToleranceV(_xSplit);
    }
  } else {
    setTolerance(params.defaultRepartTolerance);
  }
}

bool PartitioningProblemForCapo::hPartCapsFrmRgnConstr(
    unsigned splitRow, vector<double>& partCaps) {
  bool changeDefaultFlow = true;
  const BBox& bbox = _bins[0]->getBBox();
  BBox bbox0 = bbox;
  BBox bbox1 = bbox;
  const vector<const RBPCoreRow*>& rows = _bins[0]->getRows();
  double ySplit = rows[splitRow]->getYCoord();
  bbox0.yMin = ySplit;
  bbox1.yMax = ySplit;

  double p1InitSiteArea = 0;

  for (unsigned r = 0; r < splitRow; ++r) {
    p1InitSiteArea += _bins[0]->getContainedAreaInCoreRow(r);
  }

  double p0InitSiteArea = _bins[0]->getCapacity() - p1InitSiteArea;

  double reqdCellArea0 = 0;
  double reqdCellArea1 = 0;
  double containedArea0 = 0;
  double containedArea1 = 0;

  for (unsigned i = 0; i < _capoPl.regionUtilization.size(); ++i) {
    const BBox& rgnBBox = _capoPl.regionUtilization[i].first;
    double rgnUtil = _capoPl.regionUtilization[i].second;
    BBox intersect = bbox.intersect(rgnBBox);
    BBox intersect0 = bbox0.intersect(rgnBBox);
    BBox intersect1 = bbox1.intersect(rgnBBox);
    double temp = _bins[0]->getContainedAreaInBBox(intersect);
    double temp0 = _bins[0]->getContainedAreaInBBox(intersect0);
    double temp1 = _bins[0]->getContainedAreaInBBox(intersect1);
    containedArea0 += temp0;
    containedArea1 += temp1;
    reqdCellArea0 += rgnUtil * temp0;
    reqdCellArea1 += rgnUtil * temp1;
    /*One region covers almost the entire bin. Disregard this region*/
    if (temp / _bins[0]->getCapacity() > 0.97) changeDefaultFlow = false;
  }
  if (containedArea0 == 0 && containedArea1 == 0)
    changeDefaultFlow = false;
  else if (containedArea0 / p0InitSiteArea < 0.005 &&
           containedArea1 / p1InitSiteArea < 0.005)
    changeDefaultFlow = false;

  double remainingCellArea =
      _bins[0]->getTotalCellArea() - reqdCellArea0 - reqdCellArea1;

  double p0SiteArea = p0InitSiteArea - containedArea0;
  double p1SiteArea = p1InitSiteArea - containedArea1;

  if ((p0SiteArea + p1SiteArea) <= 0) {
    // partCaps[0] = reqdCellArea0/(reqdCellArea0+reqdCellArea1);
    // partCaps[1] = reqdCellArea1/(reqdCellArea0+reqdCellArea1);
    partCaps[0] = p0InitSiteArea / (p0InitSiteArea + p1InitSiteArea);
    partCaps[1] = p1InitSiteArea / (p0InitSiteArea + p1InitSiteArea);
    return (false);
  }

  // TO DO minimize rel WS here
  reqdCellArea0 += remainingCellArea * p0SiteArea / (p0SiteArea + p1SiteArea);
  reqdCellArea1 += remainingCellArea * p1SiteArea / (p0SiteArea + p1SiteArea);

  if (reqdCellArea0 < 0) reqdCellArea0 = 0;
  if (reqdCellArea1 < 0) reqdCellArea1 = 0;

  double totalCellArea = reqdCellArea0 + reqdCellArea1;
  if (totalCellArea == 0) {
    partCaps[0] = 0.5;
    partCaps[1] = 0.5;
  } else {
    partCaps[0] = reqdCellArea0 / totalCellArea;
    partCaps[1] = reqdCellArea1 / totalCellArea;
  }
  return (changeDefaultFlow);
}

bool PartitioningProblemForCapo::vPartCapsFrmRgnConstr(
    double xSplit, vector<double>& partCaps) {
  bool changeDefaultFlow = true;
  const BBox& bbox = _bins[0]->getBBox();
  BBox bbox0 = bbox;
  BBox bbox1 = bbox;
  bbox0.xMax = xSplit;
  bbox1.xMin = xSplit;

  double p0InitSiteArea = _bins[0]->getContainedAreaInBBox(bbox0);
  double p1InitSiteArea = _bins[0]->getCapacity() - p0InitSiteArea;

  double reqdCellArea0 = 0;
  double reqdCellArea1 = 0;
  double containedArea0 = 0;
  double containedArea1 = 0;

  for (unsigned i = 0; i < _capoPl.regionUtilization.size(); ++i) {
    const BBox& rgnBBox = _capoPl.regionUtilization[i].first;
    double rgnUtil = _capoPl.regionUtilization[i].second;
    BBox intersect = bbox.intersect(rgnBBox);
    BBox intersect0 = bbox0.intersect(rgnBBox);
    BBox intersect1 = bbox1.intersect(rgnBBox);
    double temp = _bins[0]->getContainedAreaInBBox(intersect);
    double temp0 = _bins[0]->getContainedAreaInBBox(intersect0);
    double temp1 = _bins[0]->getContainedAreaInBBox(intersect1);
    containedArea0 += temp0;
    containedArea1 += temp1;
    reqdCellArea0 += rgnUtil * temp0;
    reqdCellArea1 += rgnUtil * temp1;
    /*One region covers almost the entire bin. Disregard this region*/
    if (temp / _bins[0]->getCapacity() > 0.97) changeDefaultFlow = false;
  }
  if (containedArea0 == 0 && containedArea1 == 0)
    changeDefaultFlow = false;
  else if (containedArea0 / p0InitSiteArea < 0.005 &&
           containedArea1 / p1InitSiteArea < 0.005)
    changeDefaultFlow = false;

  double remainingCellArea =
      _bins[0]->getTotalCellArea() - reqdCellArea0 - reqdCellArea1;

  // cout<<"p0SiteArea "<<p0SiteArea<<" p1SiteArea "<<p1SiteArea<<endl;
  // cout<<"containedArea0 "<<containedArea0<<" containedArea1
  // "<<containedArea1<<endl;
  double p0SiteArea = p0InitSiteArea - containedArea0;
  double p1SiteArea = p1InitSiteArea - containedArea1;

  if ((p0SiteArea + p1SiteArea) <= 0) {
    // partCaps[0] = reqdCellArea0/(reqdCellArea0+reqdCellArea1);
    // partCaps[1] = reqdCellArea1/(reqdCellArea0+reqdCellArea1);
    partCaps[0] = p0InitSiteArea / (p0InitSiteArea + p1InitSiteArea);
    partCaps[1] = p1InitSiteArea / (p0InitSiteArea + p1InitSiteArea);
    return (false);
  }

  // TO DO minimize rel WS here
  reqdCellArea0 += remainingCellArea * p0SiteArea / (p0SiteArea + p1SiteArea);
  reqdCellArea1 += remainingCellArea * p1SiteArea / (p0SiteArea + p1SiteArea);

  if (reqdCellArea0 < 0) reqdCellArea0 = 0;
  if (reqdCellArea1 < 0) reqdCellArea1 = 0;

  double totalCellArea = reqdCellArea0 + reqdCellArea1;
  if (totalCellArea == 0) {
    partCaps[0] = 0.5;
    partCaps[1] = 0.5;
  } else {
    partCaps[0] = reqdCellArea0 / totalCellArea;
    partCaps[1] = reqdCellArea1 / totalCellArea;
  }
  return (changeDefaultFlow);
}

class SortByLocs {
 public:
  SortByLocs(vector<unsigned>& idxToNodeID, const PlacementWOrient& placement,
             bool sortByX)
      : _idxToNodeID(idxToNodeID), _placement(placement), _sortByX(sortByX) {}

  bool operator()(unsigned left, unsigned right) {
    if (_sortByX)
      return _placement[_idxToNodeID[left]].x <
             _placement[_idxToNodeID[right]].x;
    else
      return _placement[_idxToNodeID[left]].y <
             _placement[_idxToNodeID[right]].y;
  }

  vector<unsigned>& _idxToNodeID;
  const PlacementWOrient& _placement;
  bool _sortByX;
};

void PartitioningProblemForCapo::fixCellsUsingIntermediatePlacement()
/*
 * jflu & aaronnn:
 * run SOR to get hints of which partitions cells should be in,
 * then fix cells in the partition
 */
{
  double pctCellsToFix = _capoPl.getParams().prePartFixLimit;
  double fixMargin = _capoPl.getParams().prePartFixMargin;

  bool isHCut = (_nDim > 1) ? true : false;
  bool isVCut = (isHCut) ? false : true;

  unsigned b, numBins = _bins.size();

  vector<vector<unsigned> > nodeIds(numBins);
  vector<BBox> binBBox(numBins);

  for (b = 0; b < numBins; b++) {
    const CapoBin& curBin = *(_bins[b]);
    const BBox& binBox = curBin.getBBox();
    binBBox[b] = binBox;
    vector<unsigned>::const_iterator nItr;
    for (nItr = curBin.cellIdsBegin(); nItr != curBin.cellIdsEnd(); ++nItr)
      nodeIds[b].push_back(*nItr);
  }

  for (b = 0; b < numBins; b++) {
    // use whatever location has already been set for the cells
    const PlacementWOrient& placement = _capoPl.getPlacement();

    vector<unsigned>& nodesInBin = nodeIds[b];
    unsigned numNodesInBin = nodeIds[b].size();
    vector<unsigned> sortedNodeIdx(numNodesInBin);
    // fill up indices for sorting
    for (unsigned ni = 0; ni < numNodesInBin; ni++) sortedNodeIdx[ni] = ni;

    bool sortByX = (isVCut) ? true : false;
    sort(sortedNodeIdx.begin(), sortedNodeIdx.end(),
         SortByLocs(nodesInBin, placement, sortByX));

    // check for validity
    if (isVCut && placement[nodesInBin[sortedNodeIdx[numNodesInBin - 1]]].x -
                          placement[nodesInBin[sortedNodeIdx[0]]].x <
                      0.1 * binBBox[b].getWidth()) {
      if (_verb.getForActions() > 7) {
        cout << "j";
        cout.flush();
      }
      continue;
    }
    if (isHCut && placement[nodesInBin[sortedNodeIdx[numNodesInBin - 1]]].y -
                          placement[nodesInBin[sortedNodeIdx[0]]].y <
                      0.1 * binBBox[b].getHeight()) {
      if (_verb.getForActions() > 7) {
        cout << "j";
        cout.flush();
      }
      continue;
    }

    PartitionIds freeToMove;
    freeToMove.setToAll(2);

    unsigned numCellsToFix = unsigned(ceil(numNodesInBin * pctCellsToFix));

    // fix from left to right

    double cutOff;
    if (isVCut)
      cutOff = binBBox[b].xMin + (fixMargin * binBBox[b].getWidth());
    else
      cutOff = binBBox[b].yMin + (fixMargin * binBBox[b].getHeight());

    for (unsigned i = 0, numCellsFixed = 0;
         i < numNodesInBin && numCellsFixed < numCellsToFix; ++i) {
      if (isVCut && placement[nodesInBin[sortedNodeIdx[i]]].x > cutOff) break;
      if (isHCut && placement[nodesInBin[sortedNodeIdx[i]]].y > cutOff) break;

      if (_verb.getForActions() > 7) {
        cout << "c";
        cout.flush();
      }
      // fix this cell i
      unsigned nodeID = nodesInBin[sortedNodeIdx[i]];
      unsigned virtualNodeID = _subHG->origNode2NewIdx(nodeID);

      abkfatal(virtualNodeID < (*_fixedConstr).size(), "nodeID out of range");

      PartitionIds constrainedTo = (*_fixedConstr)[virtualNodeID];
      if (constrainedTo != freeToMove) continue;

      (*_fixedConstr)[virtualNodeID].setToPart(0);
      if (_verb.getForActions() > 7) {
        cout << "F";
        cout.flush();
      }

      ++numCellsFixed;
    }

    // fix from right to left

    if (isVCut)
      cutOff = binBBox[b].xMax - (fixMargin * binBBox[b].getWidth());
    else
      cutOff = binBBox[b].yMax - (fixMargin * binBBox[b].getHeight());

    for (unsigned i = numNodesInBin - 1, numCellsFixed = 0;
         i < numNodesInBin && numCellsFixed < numCellsToFix; --i) {
      if (isVCut && placement[nodesInBin[sortedNodeIdx[i]]].x < cutOff) break;
      if (isHCut && placement[nodesInBin[sortedNodeIdx[i]]].y < cutOff) break;

      if (_verb.getForActions() > 7) {
        cout << "C";
        cout.flush();
      }
      // fix this cell i
      unsigned nodeID = nodesInBin[sortedNodeIdx[i]];
      unsigned virtualNodeID = _subHG->origNode2NewIdx(nodeID);

      abkfatal(virtualNodeID < (*_fixedConstr).size(), "nodeID out of range");

      PartitionIds constrainedTo = (*_fixedConstr)[virtualNodeID];
      if (constrainedTo != freeToMove) continue;

      (*_fixedConstr)[virtualNodeID].setToPart(1);
      if (_verb.getForActions() > 7) {
        cout << "F";
        cout.flush();
      }

      ++numCellsFixed;
    }
  }
}

void PartitioningProblemForCapo::cellBloating() {
  const CapoBin& curBin = *(_bins[0]);

  double availableWS = 0.15 * curBin.getWhitespace();

  if (lessOrEqualDouble(availableWS, 0.)) return;

  double total = 0;
  for (vector<unsigned>::const_iterator nItr = curBin.cellIdsBegin();
       nItr != curBin.cellIdsEnd(); ++nItr) {
    if (greaterOrEqualDouble(_capoPl._congestion[(*nItr)],
                             _capoPl._congestionCutoff)) {
      total += _capoPl._congestion[(*nItr)];
    }
  }

  for (vector<unsigned>::const_iterator nItr = curBin.cellIdsBegin();
       nItr != curBin.cellIdsEnd(); ++nItr) {
    unsigned originalCellID = (*nItr);

    if (lessThanDouble(_capoPl._congestion[originalCellID],
                       _capoPl._congestionCutoff)) {
      continue;
    }

    const HGFNode& node = _subHG->getNewNodeByOrigIdx(originalCellID);

    double oldWeight = _subHG->getWeight(node.getIndex());

    // calculate amount to bloat by
    double newWeight =
        (_capoPl._congestion[originalCellID] / total) * availableWS + oldWeight;

    // cap the amount by 1.5x of the original area
    newWeight = min(1.5 * oldWeight, newWeight);

    // bloat the cell, update the weights
    _subHG->setWeight(node.getIndex(), newWeight);
  }
}

void PartitioningProblemForCapo::destroyHGraph() {
  if (_hgraph != NULL) {
    if (_ownsHGraphs)
      delete _hgraph;
    else
      _hgraph->reinstateTermWeights();
    _hgraph = _subHG = NULL;
  }
}

void PartitioningProblemForCapo::swapHGraphs(
    PartitioningProblemForCapo& newProb) {
  std::swap(_hgraph, newProb._hgraph);
  std::swap(_subHG, newProb._subHG);
  std::swap(_ownsHGraphs, newProb._ownsHGraphs);
}
