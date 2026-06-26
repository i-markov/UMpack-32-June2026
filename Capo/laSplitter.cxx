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

// created by Andrew Caldwell on 10/17/99 caldwell@cs.ucla.edu

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "laSplitter.h"

#include "PlaceEvals/bboxPlEval.h"
#include "splitLarge.h"

using std::cout;
using std::endl;
using std::vector;

bool splitHorizontally(const CapoBin& bin) {
  unsigned numRows = bin.getNumRows();
  unsigned numCells = bin.getNumCells();

  return (numRows > 1 &&
          (numCells < 15 * numRows || bin.getHeight() > bin.getWidth()));
}

LookAheadSplitter::LookAheadSplitter(CapoBin& binToSplit, CapoPlacer& capo,
                                     unsigned branchingFactor,
                                     unsigned lookAheadLvls, Verbosity verb)
    : BaseBinSplitter(binToSplit, capo, verb),
      _branchingFactor(branchingFactor),
      _lookAheadLvls(lookAheadLvls),
      _estWL(_branchingFactor, 0.0),
      _pinCounts(_branchingFactor),
      _netIsInBin(capo.partProbEdgesVisited) {
  if (_params.verb.getForActions() > 1)
    cout << " Using LookAhead Splitter:  BF == " << _branchingFactor
         << "  LA == " << _lookAheadLvls << endl;

  _params.numMLSets = branchingFactor;

  setNetsBitBoard();

  // this is so that the temp bins which are created don't
  // screw up adjacent bin's neighbor pointers.
  // 1) disconect binToSplit from its neighbors
  // 2) store binToSplit's vectors of neighbors
  // 3) clear bTS's vectors of neighbors (so it thinks it doesn't have any)
  // 4) split the bin as many times/levels as desired
  // 5) from the saved neighbor vectors, hook-up the new bins!

  _bin.unLinkNeighbors();
  _bTSLeftNeighbors = _bin._leftNeighbors;
  _bTSRightNeighbors = _bin._rightNeighbors;
  _bTSNeighborsAbove = _bin._neighborsAbove;
  _bTSNeighborsBelow = _bin._neighborsBelow;
  _bin._leftNeighbors.clear();
  _bin._rightNeighbors.clear();
  _bin._neighborsAbove.clear();
  _bin._neighborsBelow.clear();

  // unsigned numMLSets = 1;

  Verbosity sVerb = _params.verb;
  sVerb.setForActions(sVerb.getForActions() / 4);
  sVerb.setForMajStats(sVerb.getForMajStats() / 4);
  sVerb.setForSysRes(sVerb.getForSysRes() / 4);

  unsigned b;
  for (b = 0; b < _branchingFactor; b++) {
    if (_params.verb.getForActions() > 2) cout << "  Branch  " << b << endl;
    vector<CapoBin*> miniLayout;
    vector<CapoBin*> newMiniLayout;

    miniLayout.push_back(&_bin);
    const_cast<CapoPlacer&>(_capo).updateInfoAboutBin(&_bin);

    for (unsigned lvl = 0; lvl <= _lookAheadLvls; lvl++) {
      if (_params.verb.getForActions() > 3)
        cout << "    LookAhead level " << lvl << endl;
      for (unsigned bin = 0; bin < miniLayout.size(); bin++) {
        BaseBinSplitter* subSplitter;
        CapoBin& subBin = *miniLayout[bin];
        if (splitHorizontally(subBin))
          subSplitter = new SplitLargeBinHorizontally(subBin, _capo, _params);
        else
          subSplitter = new SplitLargeBinVertically(subBin, _capo, _params);

        if (miniLayout[bin] == &_bin)  // top level bin..
        {
          // this is one of the splits we're deciding between
          _topLevelBins.push_back(subSplitter->getNewBins());
        } else if (miniLayout[bin] != _topLevelBins.back()[0] &&
                   miniLayout[bin] != _topLevelBins.back()[1])
          delete miniLayout[bin];

        // yes..this is a Bad Thing...but I don't have a better soln
        // immediately w/o changing a lot of const->non-cost, and
        // making this public
        const_cast<CapoPlacer&>(_capo).updateInfoAboutBin(
            subSplitter->getNewBins()[0]);
        const_cast<CapoPlacer&>(_capo).updateInfoAboutBin(
            subSplitter->getNewBins()[1]);

        newMiniLayout.push_back(subSplitter->getNewBins()[0]);
        newMiniLayout.push_back(subSplitter->getNewBins()[1]);

        delete subSplitter;
      }
      miniLayout = newMiniLayout;
      newMiniLayout.clear();
    }

    for (unsigned bin = 0; bin < miniLayout.size(); bin++) {
      if (miniLayout[bin] != _topLevelBins.back()[0] &&
          miniLayout[bin] != _topLevelBins.back()[1]) {
        delete miniLayout[bin];
        miniLayout[bin] = NULL;
      }
    }

    _estWL[b] = evalBinsWL(_capo.getPlacement());

    if (_params.verb.getForMajStats() > 2) {
      cout << "  Measured HPWL:  " << _estWL[b] << endl;
      cout << endl;
    }
  }

  double bestWL = _estWL[0];
  unsigned bestSoln = 0;

  for (b = 1; b < _branchingFactor; b++)
    if (_estWL[b] < bestWL) {
      bestSoln = b;
      bestWL = _estWL[b];
    }

  if (_params.verb.getForMajStats() > 2)
    cout << "Best solution was " << bestSoln << endl;

  _newBins = _topLevelBins[bestSoln];
  _topLevelBins[bestSoln][0] = NULL;
  _topLevelBins[bestSoln][1] = NULL;

  // now, we have to hookup the neighbors!!
  attachNeighbors();
  _newBins[0]->_index = _bin.getIndex() + 1;
  _newBins[1]->_index = _bin.getIndex() + 2;
  _capo.nextBinNum = _bin.getIndex() + 3;

  _netIsInBin.clear();
}

void LookAheadSplitter::setNetsBitBoard() {
  const HGraphFixed& hgraph = _capo.getNetlistHGraph();

  if (_netIsInBin.getSize() < hgraph.getNumEdges())
    _netIsInBin = BitBoard(hgraph.getNumEdges());

  vector<unsigned>::const_iterator cIt;
  for (cIt = _bin.cellIdsBegin(); cIt != _bin.cellIdsEnd(); cIt++) {
    const HGFNode& node = hgraph.getNodeByIdx(*cIt);
    for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); e++)
      _netIsInBin.setBit((*e)->getIndex());
  }
}

double LookAheadSplitter::evalBinsWL(const Placement& pl) {
  double totalWL = 0;
  BBox netsBBox;

  const HGraphFixed& netlistHG = _capo.getNetlistHGraph();

  const auto& setBits = _netIsInBin.getIndicesOfSetBits();
  for (auto eId = setBits.begin(); eId != setBits.end(); eId++) {
    netsBBox.clear();
    const HGFEdge& edge = netlistHG.getEdgeByIdx(*eId);
    for (itHGFNodeLocal n = edge.nodesBegin(); n != edge.nodesEnd(); n++)
      netsBBox += pl[(*n)->getIndex()];

    totalWL += netsBBox.getHalfPerim();
  }

  return totalWL;
}

LookAheadSplitter::~LookAheadSplitter() {
  for (unsigned b = 0; b < _topLevelBins.size(); b++) {
    if (_topLevelBins[b][0]) delete _topLevelBins[b][0];
    if (_topLevelBins[b][1]) delete _topLevelBins[b][1];
    _topLevelBins[b][0] = NULL;
    _topLevelBins[b][1] = NULL;
  }
}

void LookAheadSplitter::attachNeighbors() {
  if (_newBins[0]->getCenter().x < _newBins[1]->getCenter().x) {
    abkfatal(_newBins[0]->getCenter().y == _newBins[1]->getCenter().y,
             "bins are off center-y");

    // the bins will be attached to eachother, and we
    // don't want to link this 2x...so..
    _newBins[0]->_rightNeighbors.clear();
    _newBins[1]->_leftNeighbors.clear();

    double xSplit =
        (_newBins[0]->getBBox().xMax + _newBins[0]->getBBox().xMin) / 2.0;

    vector<CapoBin*> p0NeighborsAbove, p1NeighborsAbove;
    unsigned k;
    for (k = 0; k != _bTSNeighborsAbove.size(); k++) {
      CapoBin* tempB = _bTSNeighborsAbove[k];
      double xMax = tempB->getBBox().xMax;
      double xMin = tempB->getBBox().xMin;
      if (xMin < xSplit) p0NeighborsAbove.push_back(tempB);
      if (xMax > xSplit) p1NeighborsAbove.push_back(tempB);
    }

    vector<CapoBin*> p0NeighborsBelow, p1NeighborsBelow;
    for (k = 0; k != _bTSNeighborsBelow.size(); k++) {
      CapoBin* tempB = _bTSNeighborsBelow[k];
      double xMax = tempB->getBBox().xMax;
      double xMin = tempB->getBBox().xMin;
      if (xMin < xSplit) p0NeighborsBelow.push_back(tempB);
      if (xMax > xSplit) p1NeighborsBelow.push_back(tempB);
    }

    for (k = 0; k < p0NeighborsAbove.size(); k++)
      _newBins[0]->_addNeighborAbove(p0NeighborsAbove[k]);
    for (k = 0; k < p0NeighborsBelow.size(); k++)
      _newBins[0]->_addNeighborBelow(p0NeighborsBelow[k]);
    for (k = 0; k < _bTSLeftNeighbors.size(); k++)
      _newBins[0]->_addLeftNeighbor(_bTSLeftNeighbors[k]);

    for (k = 0; k < p1NeighborsAbove.size(); k++)
      _newBins[1]->_addNeighborAbove(p1NeighborsAbove[k]);
    for (k = 0; k < p1NeighborsBelow.size(); k++)
      _newBins[1]->_addNeighborBelow(p1NeighborsBelow[k]);
    for (k = 0; k < _bTSRightNeighbors.size(); k++)
      _newBins[1]->_addRightNeighbor(_bTSRightNeighbors[k]);

    _newBins[0]->linkNeighbors();
    _newBins[1]->linkNeighbors();

    _newBins[0]->_addRightNeighbor(_newBins[1]);
    _newBins[1]->_addLeftNeighbor(_newBins[0]);
  } else  // horizontal cut
  {
    abkfatal(_newBins[0]->getCenter().x == _newBins[1]->getCenter().x,
             "bins are off center-x");

    // the bins will be attached to eachother, and we
    // don't want to link this 2x...so..
    _newBins[0]->_neighborsBelow.clear();
    _newBins[1]->_neighborsAbove.clear();

    double ySplit =
        (_newBins[0]->getBBox().yMin + _newBins[1]->getBBox().yMax) / 2.0;

    vector<CapoBin*> p0LeftNeighbors, p1LeftNeighbors;

    unsigned k;
    for (k = 0; k != _bTSLeftNeighbors.size(); k++) {
      CapoBin* tempB = _bTSLeftNeighbors[k];
      double yMax = tempB->getBBox().yMax;
      double yMin = tempB->getBBox().yMin;
      if (yMax > ySplit) p0LeftNeighbors.push_back(tempB);
      if (yMin < ySplit) p1LeftNeighbors.push_back(tempB);
    }

    vector<CapoBin*> p0RightNeighbors, p1RightNeighbors;
    for (k = 0; k != _bTSRightNeighbors.size(); k++) {
      CapoBin* tempB = _bTSRightNeighbors[k];
      double yMax = tempB->getBBox().yMax;
      double yMin = tempB->getBBox().yMin;
      if (yMax > ySplit) p0RightNeighbors.push_back(tempB);
      if (yMin < ySplit) p1RightNeighbors.push_back(tempB);
    }

    for (k = 0; k < p0LeftNeighbors.size(); k++)
      _newBins[0]->_addLeftNeighbor(p0LeftNeighbors[k]);
    for (k = 0; k < p0RightNeighbors.size(); k++)
      _newBins[0]->_addRightNeighbor(p0RightNeighbors[k]);
    for (k = 0; k < _bTSNeighborsAbove.size(); k++)
      _newBins[0]->_addNeighborAbove(_bTSNeighborsAbove[k]);

    for (k = 0; k < p1LeftNeighbors.size(); k++)
      _newBins[1]->_addLeftNeighbor(p1LeftNeighbors[k]);
    for (k = 0; k < p1RightNeighbors.size(); k++)
      _newBins[1]->_addRightNeighbor(p1RightNeighbors[k]);
    for (k = 0; k < _bTSNeighborsBelow.size(); k++)
      _newBins[1]->_addNeighborBelow(_bTSNeighborsBelow[k]);

    _newBins[0]->linkNeighbors();
    _newBins[1]->linkNeighbors();

    _newBins[0]->_addNeighborBelow(_newBins[1]);
    _newBins[1]->_addNeighborAbove(_newBins[0]);
  }
}
