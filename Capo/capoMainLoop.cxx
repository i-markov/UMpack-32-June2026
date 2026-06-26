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

// Created: May 10, 1999 by Andrew Caldwell,
//	taken from capoPlacer.cxx by Igor Markov and Andrew Caldwell

// CHANGES

// 001218  ilm   empty bins are now removed from layout

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <cstdio>
#include <sstream>
#include <string>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkfunc.h"
#include "ABKCommon/abktimer.h"
#include "ABKCommon/infolines.h"
#include "PlaceEvals/placeEvals.h"
#include "Placement/xyDraw.h"
#include "RBPlace/RBPlacement.h"
#include "Stats/trivStats.h"
#include "capoBin.h"
#include "capoPlacer.h"
// #include "../CapoWrapper/doSTA.h"

using std::cerr;
using std::cout;
using std::endl;
using std::max;
using std::ostringstream;
using std::pair;
using std::vector;

bool CapoPlacer::doOneLayer(AllowedPartDir partDir, bool& fpdSomething) {
  Timer capoLoopTimer;

  vector<CapoBin*>& activeLayout = _layout[_totalLayerNum % 2];
  vector<CapoBin*>& newLayout = _layout[(_totalLayerNum + 1) % 2];
  _curLayout = &(_layout[_totalLayerNum % 2]);

  newLayout.clear();
  _splits.clear();

  // refine each bin in the active layout, and put the
  // newly created bins into newLayout

  _totalOverfill = 0;

  if (_params.mixedSize) getBinsWMacrosInFront(activeLayout);

  /*
    //NOTE:  put back to do alternating RB
    if(partDir == VOnly)
    doVOnlyLayer(activeLayout, newLayout, _layerNum);
    else if(partDir == HOnly)
    doHOnlyLayer(activeLayout, newLayout, _layerNum);
    else
  */

  fpdSomething = doHAndVLayer(activeLayout, newLayout);

  _curLayout = &(_layout[(_totalLayerNum + 1) % 2]);  // switch to the new layer

  if (_params.verb.getForMajStats() > 0) cerr << endl;

  capoLoopTimer.stop();

  if (_params.verb.getForMajStats() > 0) {
    const HGraphWDimensions& hgraph = getNetlistHGraph();

    Timer layerStats;

    cout << endl << "# Current Time: " << TimeStamp(true);
    cout << "Stats for Layer " << _mainLayerNum << "." << _feedbackLayerNum
         << endl;

    double bboxWL = estimateWL();
    double bboxWLWt =
        _rbplace.evalTraditionalWeightedHPWL(_placement, hgraph, true);

    cout << "Layer " << _mainLayerNum << "." << _feedbackLayerNum
         << " Stats:  WL   (HP): " << bboxWL << endl;
    cout << "Layer " << _mainLayerNum << "." << _feedbackLayerNum
         << " Stats:  WtWL (HP): " << bboxWLWt << endl;
    if (_params.verb.getForMajStats() > 1) {
      double bboxChengWL = _rbplace.evalWeightedWL(_placement, hgraph, true);
      cout << "Layer " << _mainLayerNum << "." << _feedbackLayerNum
           << " Stats: WWL(Cheng): " << bboxChengWL << endl;
      //         BBoxWRSMTPlEval       rsmtWLEval (hgraph, _placement);
      //         cout<<"Layer "<<_mainLayerNum<<"."<<_feedbackLayerNum<<" Stats:
      //         WWL (RSMT): "<<rsmtWLEval.getCost()<<endl;
    }

    cout << "Layer " << _mainLayerNum << "." << _feedbackLayerNum
         << " Stats:  Time:  " << capoLoopTimer.getUserTime() << endl;

    cout << "Layer " << _mainLayerNum << "." << _feedbackLayerNum
         << " Stats:  " << MemUsage() << endl;
    cout << "Layer " << _mainLayerNum << "." << _feedbackLayerNum
         << " Stats:  " << *_params.maxMem << endl;

    _totalLayoutArea = _totalOverfill = _maxOverfill = 0;
    unsigned b;
    for (b = 0; b < _curLayout->size(); b++) {
      double cap = (*_curLayout)[b]->getCapacity();
      double ovr = (*_curLayout)[b]->getOverfill();
      _totalLayoutArea += cap;
      _totalOverfill += ovr;

      _maxOverfill = max(_maxOverfill, ovr / cap);
    }
    for (b = 0; b < _placedBins.size(); b++) {
      double cap = _placedBins[b]->getCapacity();
      double ovr = _placedBins[b]->getOverfill();
      _totalLayoutArea += cap;
      _totalOverfill += ovr;

      _maxOverfill = max(_maxOverfill, ovr / cap);
    }

    cout << "Layer " << _mainLayerNum << "." << _feedbackLayerNum
         << " Ave Overfill: " << 100 * (_totalOverfill / _totalLayoutArea)
         << "%" << endl;
    cout << "Layer " << _mainLayerNum << "." << _feedbackLayerNum
         << " Max Overfill: " << _maxOverfill * 100 << "%" << endl
         << endl;

    layerStats.stop();
    CapoPlacer::capoTimer::WLCalcTime += layerStats.getUserTime();
  }

  if (_params.doOverlapping &&
      (_mainLayerNum >= _params.startOverlappingLayer) &&
      (_mainLayerNum <= _params.endOverlappingLayer)) {
    capoLoopTimer.start();
    doOverlapping(*_curLayout);
    capoLoopTimer.stop();

    if (true || _params.verb.getForMajStats() > 1) {
      double hpWL = estimateWL();
      cout << "Layer " << _mainLayerNum << "." << _feedbackLayerNum
           << " Stats:  Post-Overlapping HP WL: " << hpWL << "\n";
      cout << "Layer " << _mainLayerNum << "." << _feedbackLayerNum
           << " Stats:  Overlapping Time:  " << capoLoopTimer.getUserTime()
           << endl;
      cout << endl;

      _totalOverfill = 0;
      for (unsigned i = 0; i < _curLayout->size(); i++)
        _totalOverfill += (*_curLayout)[i]->getOverfill();
      cout << "Layer " << _mainLayerNum << "." << _feedbackLayerNum
           << " Stats:  Ave Relative Overfill: "
           << _totalOverfill / _totalLayoutArea << endl;
    }
  }

  if (_params.plotBins) plotBins(*_curLayout);

  if (_params.saveLayerBBoxes) {
    LayoutBBoxes regions;
    for (unsigned box = 0; box < _curLayout->size(); box++)
      regions.push_back((*_curLayout)[box]->getBBox());

    ostringstream layerStr;
    layerStr << "LevelBBoxes_" << _mainLayerNum << "." << _feedbackLayerNum;
    const std::string layerFile = layerStr.str();
    regions.save(layerFile.c_str());
  }

  _totalLayerNum++;
  _feedbackLayerNum++;

  return _curLayout->size() != 0;
}

void CapoPlacer::doVOnlyLayer(vector<CapoBin*>& curLayout,
                              vector<CapoBin*>& newLayout) {
  if (_params.verb.getForMajStats() > 0) {
    cerr << "Layer{V}: " << _mainLayerNum << "." << _feedbackLayerNum
         << "  Bins: " << curLayout.size() << endl;
    cout << "Layer{V}: " << _mainLayerNum << "." << _feedbackLayerNum
         << "  Bins: " << curLayout.size() << endl;
  }

  unsigned totalNumBins = _curLayout->size() + _placedBins.size();

  abkfatal(totalNumBins > 0, "There are no bins at the start of doVOnlyLayer");

  unsigned minBinSize = getNetlistHGraph().getNumNodes() / totalNumBins;
  minBinSize = minBinSize < 50 ? 0 : (minBinSize * 2) / 3;
  // min becomes 0 when we get to small bin partitioning
  // so that we can't get stuck, unable to partition a very small
  // multi-row bin

  bool didVSplit = true;

  newLayout.clear();

  while (didVSplit) {
    didVSplit = false;
    if (_params.verb.getForMajStats() > 1) cout << "Begin V-Split pass" << endl;

    unsigned b;
    for (b = 0; b < curLayout.size(); b++) curLayout[b]->_canSplitBin = true;

    cerr << "   [";

    for (b = 0; b < curLayout.size(); b++) {
      CapoBin& bin = *curLayout[b];

      if (_params.verb.getForMajStats() > 3) {
        cout << "Bin: " << b << " [" << bin.getNumCells() << " cells]" << endl;
        if (_params.verb.getForMajStats() > 4) cout << bin << endl;
      }

      if (bin.getNumCells() == 0) {
        unsigned num = bin.getIndex();
        // bin.unLinkNeighbors();
        // delete curLayout[b];
        _placedBins.push_back(curLayout[b]);
        curLayout[b] = 0;
        if (_params.verb.getForActions() > 0)
          cout << " - moved empty bin [" << num << "] to placed bins" << endl;
      } else {
        bool binWasSplit = refineBin(bin, newLayout, VOnly, minBinSize,
                                     _params.tdCongestionCtl);
        didVSplit |= binWasSplit;
        if (binWasSplit) delete curLayout[b];
      }

      curLayout[b] = 0;
    }

    cerr << "]\n";
    if (didVSplit)  // prepare to do it again...
    {
      curLayout = newLayout;
      newLayout.clear();
    }
  }
  curLayout.clear();  // everything is in newLayout
}

void CapoPlacer::doHOnlyLayer(vector<CapoBin*>& curLayout,
                              vector<CapoBin*>& newLayout) {
  if (_params.verb.getForMajStats() > 0) {
    cerr << "Layer{H}: " << _mainLayerNum << "." << _feedbackLayerNum
         << "  Bins: " << curLayout.size() << endl;
    cout << "Layer{H}: " << _mainLayerNum << "." << _feedbackLayerNum
         << "  Bins: " << curLayout.size() << endl;
  }

  unsigned totalNumBins = _curLayout->size() + _placedBins.size();
  unsigned minBinSize = getNetlistHGraph().getNumNodes() / totalNumBins;
  minBinSize = minBinSize < 50 ? 0 : minBinSize / 2;
  // min becomes 0 when we get to small bin partitioning
  // so that we can't get stuck, unable to partition a very small
  // multi-row bin

  cerr << "   [";

  unsigned b;
  for (b = 0; b < curLayout.size(); b++) {
    curLayout[b]->_canSplitBin = true;
    curLayout[b]->_isNewBin = false;
  }

  if (_params.verb.getForMajStats() > 1) cout << "H-Split pass" << endl;

  vector<CapoBin*> tempNewLayout;
  tempNewLayout.reserve(curLayout.size() * 2);

  // now, do horizontal splits
  for (b = 0; b < curLayout.size(); b++) {
    CapoBin& bin = *curLayout[b];
    _totalOverfill += bin.getOverfill();

    if (!bin.canSplitBin()) {
      if (_params.verb.getForActions() > 5)
        cout << "Already partitioned: bin " << b << "(" << bin.getIndex() << ")"
             << endl;
      // will all be deleted after this loop
      continue;
    }

    if (_params.verb.getForMajStats() > 2) {
      cout << "Bin: " << b << " [" << bin.getNumCells() << " cells]" << endl;
      if (_params.verb.getForMajStats() > 4) cout << bin << endl;
    }

    if (bin.getNumCells() == 0) {
      unsigned num = bin.getIndex();
      // bin.unLinkNeighbors();
      // delete curLayout[b];
      _placedBins.push_back(curLayout[b]);
      curLayout[b] = 0;
      if (_params.verb.getForActions() > 0)
        cout << " - moved empty bin [" << num << "] to placed bins" << endl;
    } else
      refineBin(bin, tempNewLayout, HOnly, minBinSize, _params.tdCongestionCtl);
  }

  cerr << "]";

  // bins may have been passed over initially, then
  // split as a member of a chain..in which case they woul
  // already have been transfered to 'tempNewLayout'..
  for (b = 0; b < tempNewLayout.size(); b++) {
    if (tempNewLayout[b]->isNewBin() || tempNewLayout[b]->canSplitBin())
      newLayout.push_back(tempNewLayout[b]);
  }

  for (b = 0; b < curLayout.size(); b++) {
    if (!curLayout[b]->canSplitBin())  // was split....delete it
      delete curLayout[b];

    curLayout[b] = 0;
  }
  curLayout.clear();
}

void CapoPlacer::boost(vector<BBox>& bboxes) {
  itHGFEdgeGlobal edge;

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  bboxes.clear();
  for (edge = hgraph.edgesBegin(); edge != hgraph.edgesEnd(); edge++) {
    int eindex = (*edge)->getIndex();
    itHGFNodeLocal node;
    BBox box;
    box.clear();
    for (node = (*edge)->nodesBegin(); node != (*edge)->nodesEnd(); node++) {
      unsigned nodeId = (*node)->getIndex();
      box += _placement[nodeId];
    }
    bboxes[eindex] = box;
  }

  /*  for(edge = hgraph.edgesBegin(); edge != hgraph.edgesEnd(); edge++) {
      int eindex=(*edge)->getIndex();
      printf("Edge %d from %.2f %.2f in xrange and %.2f %.2f in yrange\n",
     eindex, bboxes[eindex].xMin, bboxes[eindex].xMax, bboxes[eindex].yMin,
     bboxes[eindex].yMax);
      }*/
}

int CapoPlacer::probeBin(CapoBin& bin) {
  unsigned numLeaves = bin.getNumCells();

  if (numLeaves == 1)
    return 0;
  else if (bin.getNumRows() == 1 && numLeaves < _params.smallPlaceThreshold)
    return 0;
  else {
    unsigned numRows = bin.getNumRows();
    Verbosity sVerb = _params.verb;
    sVerb.setForActions(sVerb.getForActions() / 2);
    sVerb.setForMajStats(sVerb.getForMajStats() / 2);
    sVerb.setForSysRes(sVerb.getForSysRes() / 2);

    vector<CapoBin*> newBins;
    updateInfoAboutBin(&bin);

    unsigned lAhead = _params.lookAhead;
    if (numLeaves < 500) lAhead = 0;

    double stretchFactor =
        (bin.getAvgRowSpacing() > 1.5 * bin.getAvgCellHeight() ? 2 : 1);
    bool splitHoriz = stretchFactor * bin.getHeight() > bin.getWidth();

    if (!splitHoriz)
      splitHoriz = (2 * bin.getHeight() > bin.getWidth() &&
                    numLeaves <= _params.smallPlaceThreshold * numRows);

    if (numRows > 1 && splitHoriz)
      return 1;
    else
      return 2;
  }
  return 0;
}

void CapoPlacer::boostNets(CapoBin& bin, vector<BBox>& bboxes) {
  const vector<unsigned>& cellids = bin.getCellIds();
  int dir = probeBin(bin);
  int BF = _params.boostFactor;
  int boost = 0;

  if (dir == 0) return;

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  for (unsigned i = 0; i < cellids.size(); i++) {
    HGFNode node = hgraph.getNodeByIdx(cellids[i]);
    unsigned nid = node.getIndex();
    itHGFEdgeGlobal edge;
    for (edge = node.edgesBegin(); edge != node.edgesEnd(); edge++) {
      if ((*edge)->getWeight() == BF) continue;
      int eindex = (*edge)->getIndex();

      if (dir == 2) {
        // vertical
        if (bboxes[eindex].xMin < _placement[nid].x &&
            _placement[nid].x < bboxes[eindex].xMax)
          continue;
        else if (bboxes[eindex].xMin == _placement[nid].x &&
                 bboxes[eindex].xMax != _placement[nid].x) {
          (*edge)->setWeight(BF);
          boost++;
        } else if (bboxes[eindex].xMax == _placement[nid].x &&
                   bboxes[eindex].xMin != _placement[nid].x) {
          (*edge)->setWeight(BF);
          boost++;
        }
      } else if (dir == 1) {
        // horizontal
        if (bboxes[eindex].yMin < _placement[nid].y &&
            _placement[nid].y < bboxes[eindex].yMax)
          continue;
        else if (bboxes[eindex].yMin == _placement[nid].y &&
                 bboxes[eindex].yMax != _placement[nid].y) {
          (*edge)->setWeight(BF);
          boost++;
        } else if (bboxes[eindex].yMax == _placement[nid].y &&
                   bboxes[eindex].yMin != _placement[nid].y) {
          (*edge)->setWeight(BF);
          boost++;
        }
      }
    }
  }
  printf("BOOSTED %d nets\n", boost);
}

void CapoPlacer::deboostNets(CapoBin& bin) {
  const vector<unsigned>& cellids = bin.getCellIds();

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  for (unsigned i = 0; i < cellids.size(); i++) {
    HGFNode node = hgraph.getNodeByIdx(cellids[i]);
    itHGFEdgeGlobal edge;
    for (edge = node.edgesBegin(); edge != node.edgesEnd(); edge++)
      (*edge)->setWeight(1);
  }
}

double CapoPlacer::squareCost() {
  double square = 0;
  itHGFEdgeGlobal edge;

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  for (edge = hgraph.edgesBegin(); edge != hgraph.edgesEnd(); edge++) {
    itHGFNodeLocal node;
    BBox box;
    box.clear();
    for (node = (*edge)->nodesBegin(); node != (*edge)->nodesEnd(); node++) {
      unsigned nodeId = (*node)->getIndex();
      box += _placement[nodeId];
    }
    square += pow(box.xMax - box.xMin, 2) + pow(box.yMax - box.yMin, 2);
  }
  return square;
}

bool CapoPlacer::doHAndVLayer(vector<CapoBin*>& curLayout,
                              vector<CapoBin*>& newLayout) {
  vector<BBox> bboxes;
  if (_params.boost) {
    cout << "Calculating the bounding boxes of nets" << endl;
    bboxes = vector<BBox>(getNetlistHGraph().getNumEdges());
  }

  if (_params.verb.getForMajStats() > 0) {
    cerr << "Layer:    " << _mainLayerNum << "." << _feedbackLayerNum
         << "  Bins: " << curLayout.size() << endl;
    cout << "Layer:    " << _mainLayerNum << "." << _feedbackLayerNum
         << "  Bins: " << curLayout.size() << endl;
  }

  cerr << "   [";

  for (unsigned b = 0; b < curLayout.size(); b++)
    if (curLayout[b] != NULL) curLayout[b]->_canSplitBin = true;

  bool atleast1FPed = false;
  if (_params.mixedSize && _rbplace.getNumMacros() > 0) {
    Timer fpTimer;

    // <aaronnn> reset the bitvector that marks the nodes/bins needing FP
    //_nodesNeedingFP.clear(); // <aaronnn> don't clear - remember all nodes
    //that needed FP

#ifdef USEPATOMA
    if (_params.usePatoma) {
      // use PATOMA instead of the standard FP procedure
      atleast1FPed = doPatomaLayer(curLayout, newLayout);

      // put the PATOMAed bins to the new layout SO THAT
      // they are not refined (in code below)
      getPatomaedBinsIntoNewLayout(curLayout, newLayout);
      processMacrosAfterlayer(newLayout, PATOMA_LAYER);
    } else
#endif
    {
      // do the standard floorplanning
      atleast1FPed = doFPLayer(curLayout, newLayout);

      // put the FPed bins to the new layout SO THAT
      // they are not refined (in code below)
      getFPedBinsIntoNewLayout(curLayout, newLayout);
      processMacrosAfterlayer(newLayout, FP_LAYER);
    }

    // can skip if no bins were floorplanned
    if (atleast1FPed) {
      for (unsigned i = 0; i < newLayout.size(); ++i) {
        if (newLayout[i] != NULL) {
          // reset the core rows and row offsets
          // and recalculate the BBox (may have shrunk)
          newLayout[i]->resetRows();
          // the BBox may have changed so reposition
          // cells that are left in the fpd bins
          updateInfoAboutBin(newLayout[i]);
        }
      }
      for (unsigned i = 0; i < curLayout.size(); ++i) {
        if (curLayout[i] != NULL) {
          // reset the core rows and row offsets
          // and recalculate the BBox (may have shrunk)
          curLayout[i]->resetRows();
          // the BBox may have changed so reposition
          // cells that are left in the bins to be parted
          updateInfoAboutBin(curLayout[i]);
        }
      }
      for (unsigned i = 0; i < _placedBins.size(); ++i) {
        // reset the core rows and row offsets
        // and recalculate the BBox (may have shrunk)
        _placedBins[i]->resetRows();
      }
    }

    fpTimer.stop();
    CapoPlacer::capoTimer::FloorplanTime += fpTimer.getUserTime();
  }

  // REMINDER: after this point, all floorplanned bins are not in
  // "curLayout" but in "newLayout" instead

  vector<CapoBin*> splitLayout;

  unsigned minBinSize = 0;  // will partition everything
  for (unsigned b = 0; b < curLayout.size(); b++) {
    if (curLayout[b] != NULL) {
      curLayout[b]->_canSplitBin = true;
      CapoBin& bin = *curLayout[b];
      _totalOverfill += bin.getOverfill();

      abkwarn(bin.getNumRows() > 0 || bin.getNumCells() == 0,
              "non-empty bin with no core rows.");  //-----
      if (_params.verb.getForMajStats() > 1) {
        cout << "Bin: " << b << " [" << bin.getNumCells() << " cells]" << endl;
        if (_params.verb.getForMajStats() > 1) cout << bin << endl;
      }

      if (bin.getNumCells() == 0) {
        // treat empty bins as placed
        unsigned num = bin.getIndex();
        _placedBins.push_back(curLayout[b]);
        curLayout[b] = NULL;
        if (_params.verb.getForActions() > 0)
          cout << " - moved empty bin [" << num << "] to placed bins" << endl;
      }
      /*
        else if(bin.getNumSites() == 0)
        {
        _placedBins.push_back(curLayout[b]);
        curLayout[b] = 0;
        }
      */
      else if (bin.getNumRows() == 0) {
        // treat empty bins as placed
        unsigned num = bin.getIndex();
        _placedBins.push_back(curLayout[b]);
        curLayout[b] = NULL;
        if (_params.verb.getForActions() > 0)
          cout << " - moved non-empty bin [" << num
               << "] with no core rows to placed bins" << endl;
      } else {
        // the general case when there are some cells in the bin
        if (_params.boost)
          if (_mainLayerNum <= _params.boostLayer) boostNets(bin, bboxes);

        bool deleteBin = refineBin(bin, splitLayout, HAndV, minBinSize,
                                   _params.tdCongestionCtl);

        if (deleteBin) {
          bin._sibling->_sibling = bin._sibling;
        } else {
          // zero out the reference to this bin
          // so it will not be deleted afterwards
          curLayout[b] = NULL;
        }
      }
    }
  }  // done considering each non-NULL bin in "curLayout"

  if (_ispd04CongMap != NULL) {
    Timer congMapTimer;

    if (_ispd04CongMapPlacement != NULL) {
      // fill in this placement with information from the old and new placements
      // depending on the state of ECO for each bin
      for (unsigned i = 0; i < _ispd04CongMapPlacement->size(); ++i) {
        const CapoBin* bin = _cellToBinMap[i];

        if (bin == NULL || !bin->isEcoAllowed()) {
          (*_ispd04CongMapPlacement)[i] = _placement[i];
          _ispd04CongMapPlacement->setOrient(i, _placement.getOrient(i));
        } else {
          (*_ispd04CongMapPlacement)[i] = getOraclePlacement()[i];
          _ispd04CongMapPlacement->setOrient(i,
                                             getOraclePlacement().getOrient(i));
        }
      }
    }

    const HGraphWDimensions& hgraph = getNetlistHGraph();

    _ispd04CongMap->calculateDemand(hgraph);

    pair<double, double> totalDemand =
        _ispd04CongMap->getDemand(_rbplace.getBBox(false));
    double avgDemand = (totalDemand.first + totalDemand.second) /
                       _rbplace.getBBox(false).getArea();
    double fracAvgDemand = 1.5 * avgDemand;

    // shift cutlines of all the bins that were split this layer
    for (unsigned i = 0; i < _splits.size(); ++i) {
      if (_splits[i].p1 == NULL) {
        newLayout.push_back(_splits[i].p0);
        continue;
      }

      BaseBinSplitter splitter(*_splits[i].parent, *this,
                               _params.splitterParams);

      vector<vector<unsigned> > cellsInBins(2);
      cellsInBins[0] = _splits[i].p0->getCellIds();
      cellsInBins[1] = _splits[i].p1->getCellIds();

      double p0MaxCap = _splits[i].p0->getMaxRecommendedCellArea();
      double p1MaxCap = _splits[i].p1->getMaxRecommendedCellArea();

      vector<double> partAreas(2, 0.);
      for (unsigned c = 0; c < cellsInBins[0].size(); ++c)
        partAreas[0] += hgraph.getWeight(cellsInBins[0][c]);

      for (unsigned c = 0; c < cellsInBins[1].size(); ++c)
        partAreas[1] += hgraph.getWeight(cellsInBins[1][c]);

      if (_splits[i].splitHoriz) {
        if (splitter.shouldCongestionShiftH(
                _splits[i].splitRow, *getISPD04CongMap(), fracAvgDemand)) {
          pair<double, double> WSWeight = splitter.getCongestionRatiosH(
              _splits[i].splitRow, *getISPD04CongMap());
          if (equalDouble(WSWeight.first, WSWeight.second) ||
              equalDouble(WSWeight.first, 0.) ||
              equalDouble(WSWeight.second, 0.)) {
            WSWeight = std::make_pair(1., 1.);
          }
          vector<double> actualCaps(2, 0.);
          unsigned newSplitRow = splitter.findBestSplitRow(
              partAreas, WSWeight.first, WSWeight.second, actualCaps);

          if (_params.splitterParams.eco && _params.splitterParams.ecoShift &&
              _splits[i].p0->isEcoAllowed() && _splits[i].p1->isEcoAllowed() &&
              _splits[i].splitRow != newSplitRow) {
            splitter.shiftOraclePlacementH(
                _splits[i].parent->getRows()[_splits[i].splitRow]->getYCoord(),
                _splits[i].parent->getRows()[newSplitRow]->getYCoord());
          }

          splitter.createHSplitBins(cellsInBins, newSplitRow, p0MaxCap,
                                    p1MaxCap);

          vector<CapoBin*>& newBins = splitter.getNewBins();

          for (unsigned b = 0; b < newBins.size(); ++b) {
            newLayout.push_back(newBins[b]);
            updateInfoAboutBin(newBins[b]);
          }

          delete _splits[i].p0;
          _splits[i].p0 = NULL;
          delete _splits[i].p1;
          _splits[i].p1 = NULL;
        } else {  // do not shift for congestion, use the old bins
          newLayout.push_back(_splits[i].p0);
          newLayout.push_back(_splits[i].p1);
        }
      } else {  // Vert split
        if (splitter.shouldCongestionShiftV(
                _splits[i].xSplit, *getISPD04CongMap(), fracAvgDemand)) {
          pair<double, double> WSWeight = splitter.getCongestionRatiosV(
              _splits[i].xSplit, *getISPD04CongMap());
          if (equalDouble(WSWeight.first, WSWeight.second) ||
              equalDouble(WSWeight.first, 0.) ||
              equalDouble(WSWeight.second, 0.)) {
            WSWeight = std::make_pair(1., 1.);
          }
          vector<double> actualCaps(2, 0.);
          double newXSplit = _splits[i].parent->findXSplit(
              partAreas[0], partAreas[1], 0, WSWeight.first, WSWeight.second,
              actualCaps);

          if (_params.splitterParams.eco && _params.splitterParams.ecoShift &&
              _splits[i].p0->isEcoAllowed() && _splits[i].p1->isEcoAllowed() &&
              _splits[i].xSplit != newXSplit) {
            splitter.shiftOraclePlacementV(_splits[i].xSplit, newXSplit);
          }

          splitter.createVSplitBins(cellsInBins, newXSplit, p0MaxCap, p1MaxCap);

          vector<CapoBin*>& newBins = splitter.getNewBins();

          for (unsigned b = 0; b < newBins.size(); ++b) {
            newLayout.push_back(newBins[b]);
            updateInfoAboutBin(newBins[b]);
          }

          delete _splits[i].p0;
          _splits[i].p0 = NULL;
          delete _splits[i].p1;
          _splits[i].p1 = NULL;
        } else {  // do not shift for congestion, use the old bins
          newLayout.push_back(_splits[i].p0);
          newLayout.push_back(_splits[i].p1);
        }
      }
    }
  } else {
    // otherwise put the split bins onto newLayout
    newLayout.insert(newLayout.end(), splitLayout.begin(), splitLayout.end());
  }

  for (unsigned i = 0; i < curLayout.size(); ++i) {
    if (curLayout[i] != NULL) delete curLayout[i];
  }
  curLayout.clear();

  cerr << "]";

  return atleast1FPed;
}

bool CapoPlacer::doFPLayer(vector<CapoBin*>& curLayout,
                           vector<CapoBin*>& newLayout) {
  if (_params.verb.getForMajStats() > 5) {
    cout << "Layer:    " << _mainLayerNum << "." << _feedbackLayerNum
         << "  Bins: " << curLayout.size() << endl;
    cout << "Doing Floorplanning Pass " << endl;
  }

  bool atleast1FPed = false;

  bool atBot = atBottomLayers();

  for (unsigned b = 0; b < curLayout.size(); b++) {
    if (curLayout[b] != 0) {
      const HGraphWDimensions& hgraph = getNetlistHGraph();

      CapoBin& bin = *curLayout[b];
      unsigned numLeaves = bin.getNumCells();

      if (numLeaves == 0) continue;

      vector<unsigned> macrosInBin;
      for (unsigned c = 0; c < bin.getNumCells(); c++) {
        unsigned idx = bin.getCellIds()[c];
        if (hgraph.isNodeMacro(idx)) {
          macrosInBin.push_back(idx);
        }
      }
      if (_params.lookAheadFP &&
		     _mainLayerNum > 0 &&
             macrosInBin.size() > 0 &&
             !atBot &&
             !bin._needFP /*&&
             (_mainLayerNum < 6 || bin.getNumObstacles())*/)
         {
        // <aaronnn> run a 'faster' FP on the current bin to see if a legal
        // solution is possible
        if (!lookAheadFPSuccess(b, macrosInBin, curLayout, newLayout)) {
          // TODO: rip stuff out of lookAheadFPSuccess() and call a function
          // here like mergeBinsForBacktrack()
          continue;
        }
      }

      bin._floorplanned = false;

      bool desiredDir = getDefaultSplitDir(bin);
      bool needFP = bin._needFP;
      double splitPt;

      if (!needFP)
        needFP = bin.fpCondMet(desiredDir, splitPt);
      else {
        //(void)bin.fpCondMet(desiredDir, splitPt);
        if (_params.verb.getForActions() > 0)
          cout << endl << "FPing hierarchical bin " << bin.getName() << endl;
      }

      if (!needFP) continue;

      // the bin needs floorplanning before cut can continue

      if (_params.verb.getForMajStats() > 5) cout << bin << endl;

      // <aaronnn> find the 2 largest macros in the bin
      int maxAreaId1 = -1;
      int maxAreaId2 = -1;
      double maxArea1 = -1;
      double maxArea2 = -1;
      for (vector<unsigned>::const_iterator nItr = bin.cellIdsBegin();
           nItr != bin.cellIdsEnd(); ++nItr) {
        unsigned idx = *nItr;
        if (hgraph.isNodeMacro(idx)) {
          double cellWidth = hgraph.getNodeWidth(idx);
          double cellHeight = hgraph.getNodeHeight(idx);
          double cellArea = cellWidth * cellHeight;
          if (cellArea > maxArea1) {
            maxArea2 = maxArea1;
            maxAreaId2 = maxAreaId1;
            maxArea1 = cellArea;
            maxAreaId1 = int(idx);
          }
        }
      }
      // decide if we're floorplanning because of a giant macro
      bool floorplanningGiantMacro = false;
      unsigned numNodesNeedingFP = 0;
      for (unsigned i = 0; i < _nodesNeedingFP.size(); i++) {
        if (_nodesNeedingFP[i] == true) numNodesNeedingFP++;
      }
      if (numNodesNeedingFP == 1) {
        if (maxAreaId1 != -1 && maxAreaId2 != -1 && maxArea1 > 8 * maxArea2) {
          // floorplanning more than one macro,
          // largest macro is much larger than next macro
          floorplanningGiantMacro = true;
        } else if (maxAreaId2 == -1) {
          // only one macro in bin,
          // I consider macros always 'giant' wrt std cells
          floorplanningGiantMacro = true;
        }
      }

      bool mergedSibling = false;
      bool success = false;

      if (_params.useAnalytFP && floorplanningGiantMacro) {
        success = floorplanBinWithSOR(bin);
        //           success = floorplanBinWithFlows(bin);
      } else {
        success = floorplanBin(bin);
      }

      if (!success) {
        // jflu
        const BBox& binBBox = bin.getBBox();
        // coerce all the cells in the bin
        for (vector<unsigned>::const_iterator nItr = bin.cellIdsBegin();
             nItr != bin.cellIdsEnd(); ++nItr) {
          unsigned nodeId = (*nItr);

          unsigned angle = _placement.getOrient(nodeId).getAngle();
          bool notRotated = angle == 0 || angle == 180;
          double nodeWidth = notRotated ? hgraph.getNodeWidth(nodeId)
                                        : hgraph.getNodeHeight(nodeId);
          double nodeHeight = notRotated ? hgraph.getNodeHeight(nodeId)
                                         : hgraph.getNodeWidth(nodeId);

          Point oldPlace(_placement[nodeId]);

          if (greaterThanDouble(_placement[nodeId].x + nodeWidth,
                                binBBox.xMax)) {
            _placement[nodeId].x = binBBox.xMax - nodeWidth;
          }
          if (lessThanDouble(_placement[nodeId].x, binBBox.xMin)) {
            _placement[nodeId].x = binBBox.xMin;
          }

          if (greaterThanDouble(_placement[nodeId].y + nodeHeight,
                                binBBox.yMax)) {
            _placement[nodeId].y = binBBox.yMax - nodeHeight;
          }
          if (lessThanDouble(_placement[nodeId].y, binBBox.yMin)) {
            _placement[nodeId].y = binBBox.yMin;
          }

          if (oldPlace != _placement[nodeId]) {
            if (_params.wtCut || _params.doFastPlaceMoves) {
              Timer psTimer;
              updatePointSetsForMovedCell(nodeId, oldPlace, _placement[nodeId]);
              psTimer.stop();
              CapoPlacer::capoTimer::PointSetTime += psTimer.getUserTime();
            }
          }
        }
      }

      unsigned binMacroNum = 0;
      for (vector<unsigned>::const_iterator cellIdPtr = bin.cellIdsBegin();
           cellIdPtr != bin.cellIdsEnd(); cellIdPtr++) {
        if (hgraph.isNodeMacro(*cellIdPtr)) binMacroNum++;
      }

      FloorplanStats::numInstances++;
      FloorplanStats::instancesMaxMacroNum =
          max(FloorplanStats::instancesMaxMacroNum, binMacroNum);
      FloorplanStats::totalFPcells += bin.getNumCells();
      FloorplanStats::cellsInInstance.push_back(bin.getNumCells());
      FloorplanStats::totalFPmacros += binMacroNum;
      FloorplanStats::macrosInInstance.push_back(binMacroNum);
      if (_params.FPEvadeObstacles) {
        FloorplanStats::totalFPobstacles += bin.getNumObstacles();
        FloorplanStats::obstaclesInInstance.push_back(bin.getNumObstacles());
      }

      if (success)
        FloorplanStats::numSuccessfulInstances++;
      else {
        FloorplanStats::numFailedInstances++;
        FloorplanStats::failedInstancesMaxMacroNum =
            max(FloorplanStats::failedInstancesMaxMacroNum, binMacroNum);
      }

      if (!success) {
        // now go up a hierarchy and solve again
        CapoBin* mergedBin = mergeBinWithSibling(bin);
        if (mergedBin != NULL) {
          unsigned numMacros = 0;
          for (unsigned c = 0; c < mergedBin->getNumCells(); c++) {
            unsigned idx = mergedBin->getCellIds()[c];
            if (hgraph.isNodeMacro(idx)) numMacros++;
          }

          // if bin doesn't have too many macros or we aren't doing selective
          // floorplanning, mark all for FP to not create problems
          if (numMacros <= 50 || (!_params.selectiveFloorplanning &&
                                  numMacros <= _params.scampiThresh)) {
            markNodesNeedFP(mergedBin->getCellIds());
          } else if (bin._sibling != NULL) {
            // check FP conditions on sibling bin as well (mark nodes for FP)
            double splitPt;
            bool desiredDir = getDefaultSplitDir(*(bin._sibling));
            (void)bin._sibling->fpCondMet(desiredDir, splitPt);
          }

          mergedBin->_needFP = true;
          updateCell2Bin(mergedBin);
          newLayout.push_back(mergedBin);
          if (_params.verb.getForActions() > 0)
            cout << endl
                 << "Created Hierarchical Bin " << mergedBin->getName()
                 << " for FPing" << endl;

          cerr << "i";
          bin._needFP = true;
          bin._floorplanned = false;
          mergedSibling = true;
        } else {
          // sibling not found, then do not
          // refloorplan. use current illegal solution

          // floorplan with only area objective
          // success = floorplanBin(bin, false);
          cerr << "o";
          bin._needFP = false;
          bin._floorplanned = true;
          // don't push this yet. push this after layer is done
          // newLayout.push_back(&bin);
        }
      } else {
        bin._needFP = false;
        bin._floorplanned = true;
        cerr << "o";
      }

      if (bin._floorplanned) atleast1FPed = true;

      if (mergedSibling) {
        bool foundSibling = false;
        CapoBin* bin2 = bin._sibling;
        if (bin._needFP) {
          for (unsigned b2 = 0; b2 < curLayout.size(); b2++) {
            if (curLayout[b2] != 0)
              if (curLayout[b2] == bin2) {
                // bin2->unLinkNeighbors();
                delete bin2;
                curLayout[b2] = 0;
                foundSibling = true;
                break;
              }
          }
          // bin.unLinkNeighbors();
        }
        if (!foundSibling) bin2->_sibling = bin2;
        delete curLayout[b];
        curLayout[b] = 0;
      }
    }
  }
  return (atleast1FPed);
}

#ifdef USEPATOMA
bool CapoPlacer::doPatomaLayer(vector<CapoBin*>& curLayout,
                               vector<CapoBin*>& newLayout) {
  if (_params.verb.getForMajStats() > 5) {
    cout << "Layer:    " << _mainLayerNum << "." << _feedbackLayerNum
         << "  Bins: " << curLayout.size() << endl;
    cout << "Doing PATOMA Pass " << endl;
  }

  bool atleast1FPed = false;

  // for each non-NULL bin
  for (unsigned b = 0; b < curLayout.size(); b++)
    if (curLayout[b] != NULL) {
      CapoBin& bin = *curLayout[b];
      unsigned numLeaves = bin.getNumCells();
      unsigned numRows = bin.getNumRows();

      if (numLeaves > 0) {
        // determine whether FP is needed
        double stretchFactor =
            (bin.getAvgRowSpacing() > 1.5 * bin.getAvgCellHeight() ? 2 : 1);
        bool splitHoriz = stretchFactor * bin.getHeight() > bin.getWidth();

        double vertFactor = _params.vertFactor;  // 7.0;
        bool foundMacros = false;
        if (!splitHoriz && _params.mixedSize) {
          const vector<unsigned>& cellIds = bin.getCellIds();
          unsigned i;
          for (i = 0; i < cellIds.size(); i++) {
            unsigned idx = cellIds[i];
            if (hgraph.isNodeMacro(idx)) {
              vertFactor = 0;
              foundMacros = true;
              break;
            }
          }
        }

        if (!splitHoriz)
          splitHoriz = (vertFactor * bin.getHeight() > bin.getWidth() &&
                        numLeaves <= _params.smallPlaceThreshold * numRows &&
                        bin.getRelativeWhitespace() <= 0.3 && numRows <= 5);

        bool desiredDir = splitHoriz;
        double splitPt;
        bool needFP = (bin._needFP || bin.fpCondMet(desiredDir, splitPt));

        //             //if(!foundMacros && !needFP)
        //             //needFP = false;
        //             if(!needFP)
        //                needFP = bin.fpCondMet(desiredDir, splitPt);
        //             else
        //                if(_params.verb.getForActions() > 0)
        //                   cout<<endl<<"FPing hierarchical bin
        //                   "<<bin.getName()<<endl;

        // determine whether Patoma is needed here
        bool needPatoma = bin.isPatomaEndCase();
        //             needPatoma = needPatoma || bin.patomaCondMet();
        needPatoma = needFP;  // || bin.isPatomaEndCase();

        // some bins can be unfloorplanned since
        // there are no macros
        if (needPatoma) {
          if (_params.verb.getForMajStats() > 5) cout << bin << endl;

          // -----floorplan the bin-----
          bool success = false;
          if (bin.getNumCells() == 1 && false) {
            // special cases with at most 1 cell
            unsigned cellId = *(bin.cellIdsBegin());
            if (hgraph.isNodeMacro(cellId)) {
              _placement[cellId] = bin.getBBox().getBottomLeft();
              double nodeWidth = hgraph.getNodeWidth(cellId);
              double nodeHeight = hgraph.getNodeHeight(cellId);
              if (nodeWidth <= bin.getBBox().getWidth() &&
                  nodeHeight <= bin.getBBox().getHeight())
                success = true;
            }
          } else if (bin.isPatomaPacked()) {
            // if previously proven to be successful,
            // then skip it (NOTE: if previously failed,
            // then the bin wont' exist in curLayout in the
            // first place!)
            success = true;
          } else {
            // general cases with at least 2 cells
            success = patomaBin(bin);
          }
          // -----finish floorplanning-----

          // -----statistics only-----
          int binMacroNum = 0;
          for (vector<unsigned>::const_iterator cellIdPtr = bin.cellIdsBegin();
               cellIdPtr != bin.cellIdsEnd(); cellIdPtr++) {
            if (hgraph.isNodeMacro(*cellIdPtr)) binMacroNum++;
          }

          FloorplanStats::numInstances++;
          FloorplanStats::instancesMaxMacroNum =
              max(FloorplanStats::instancesMaxMacroNum, binMacroNum);

          if (success)
            FloorplanStats::numSuccessfulInstances++;
          else {
            FloorplanStats::numFailedInstances++;
            FloorplanStats::failedInstancesMaxMacroNum =
                max(FloorplanStats::failedInstancesMaxMacroNum, binMacroNum);
          }
          // ----done with statistics----

          bool mergedSibling = false;
          if (!success) {
            // the floorplanning wasn't successful
            // now go up a hierarchy and solve again
            vector<CapoBin*> mix;
            CapoBin* bin1 = &bin;
            CapoBin* bin2 = bin._sibling;
            if (bin1 != bin2) {
              mix.push_back(bin1);
              mix.push_back(bin2);

              int dirc = 2;
              if (bin1->getBBox().xMin == bin2->getBBox().xMin &&
                  bin1->getBBox().xMax == bin2->getBBox().xMax)
                dirc = 1;
              if (bin1->getBBox().yMin == bin2->getBBox().yMin &&
                  bin1->getBBox().yMax == bin2->getBBox().yMax)
                dirc = 0;
              if (dirc == 1 || dirc == 0) {
                // the merged bins are indeed siblings
                // then merge the bins and mark it as
                // Patoma end case
                CapoBin* bin = new CapoBin(mix, dirc);
                bin->_needFP = true;
                bin->setPatomaEndCase(true);

                updateCell2Bin(bin);
                newLayout.push_back(bin);
                mergedSibling = true;

                if (_params.verb.getForActions() > 0)
                  cout << "Created Hierarchical Bin " << bin->getName()
                       << " as PATOMA end case." << endl;
              }
            }

            if (!mergedSibling) {
              // if sibling not found, then do not
              // refloorplan. use current illegal solution
              cerr << "o";
              bin._needFP = false;
              bin._floorplanned = true;
            } else {
              // if sibling is found, then go up and reflooplan
              // WARNING: the merged bin "*bin" (NOT THIS ONE)
              //  is always added to "newLayout" above
              cerr << "i";
              bin._needFP = true;
              bin._floorplanned = false;
            }
          }  // done handling the non-successful case
          else {
            // successful or previously proven to be successful
            bin._needFP = false;
            bin._floorplanned = true;
            cerr << ((bin.isPatomaPacked()) ? "+" : "o");

          }  // done handling the successful case

          // deal with PATOMA flags here
          bin.setPatomaPacked(true);
          if (success) {
            // PatomaEndCase remains the same, i.e.
            //   endCase => endCase and commit macros
            //   non-endCase =>non-endCase and not commit macros
            bin.setPatomaMacrosCommitted(bin.isPatomaEndCase());

          } else {
            if (bin.isPatomaEndCase()) {
              // commit if it is an end case
              bin.setPatomaMacrosCommitted(true);
            } else {
              // forced to commit since it can't be merged
              bin.setPatomaMacrosCommitted(!mergedSibling);
            }
            bin.setPatomaEndCase(true);
          }
          // done dealing with PATOMA flags

          if (bin._floorplanned) atleast1FPed = true;

          // further handling when siblings are merged
          if (mergedSibling) {
            // if a bin is merged with its sibling,
            // the merged bin is already added to "newLayout" up there
            // its slibling and itself are deleted
            bool foundSibling = false;
            CapoBin* bin2 = bin._sibling;
            if (bin._needFP) {
              for (unsigned b2 = 0; b2 < curLayout.size(); b2++)
                if (curLayout[b2] != NULL) {
                  // foreach non-NULL bin in curLayout, check if it is
                  // the sibling of "bin", if so, delete it
                  if (curLayout[b2] == bin2) {
                    delete bin2;
                    curLayout[b2] = NULL;
                    foundSibling = true;
                    break;
                  }
                }
            }
            if (!foundSibling) {
              // if the sibling of of "bin" is the not in "curLayout",
              // then its sibling's slibling is itself
              bin2->_sibling = bin2;
            }

            delete curLayout[b];
            curLayout[b] = NULL;
          }  // done handling the merging of siblings
        }
      }
    }  // done with all non-NULL bins
  return (atleast1FPed);
}
#endif

CapoBin* CapoPlacer::mergeBinWithSibling(CapoBin& bin) {
  vector<CapoBin*> mix;
  CapoBin* bin1 = &bin;
  CapoBin* bin2 = bin._sibling;

  if (bin1 != bin2) {
    mix.push_back(bin1);
    mix.push_back(bin2);
    int dirc = 2;
    if (equalDouble(bin1->getBBox().xMin, bin2->getBBox().xMin) &&
        equalDouble(bin1->getBBox().xMax, bin2->getBBox().xMax))
      dirc = 1;
    if (equalDouble(bin1->getBBox().yMin, bin2->getBBox().yMin) &&
        equalDouble(bin1->getBBox().yMax, bin2->getBBox().yMax))
      dirc = 0;
    if (dirc == 1 || dirc == 0) {
      CapoBin* mergedBin = new CapoBin(mix, dirc);
      return mergedBin;
    }
  }

  if (_params.verb.getForActions() > 1) {
    cout << endl
         << "Could not find a sibling bin! Going to keep bad floorplan."
         << endl;
  }
  return NULL;
}

class FPCheckOrderByDecreasingArea {
  const HGraphWDimensions& _hgraph;

 public:
  FPCheckOrderByDecreasingArea(const HGraphWDimensions& hgraph)
      : _hgraph(hgraph) {}
  bool operator()(unsigned i, unsigned j) {
    return (_hgraph.getNodeWidth(i) * _hgraph.getNodeHeight(i) >
            _hgraph.getNodeWidth(j) * _hgraph.getNodeHeight(j));
  }
};

bool CapoPlacer::lookAheadFPSuccess(unsigned b, vector<unsigned> macrosInBin,
                                    vector<CapoBin*>& curLayout,
                                    vector<CapoBin*>& newLayout) {
  const HGraphWDimensions& hgraph = getNetlistHGraph();

  CapoBin& bin = *curLayout[b];

  // <aaronnn> run a 'faster' FP on the current bin to see if a legal solution
  // is possible
  bool lookAheadFPSuccess = true;

  // mark nodes FPed to perform a fuzzy check for legality
  const vector<unsigned>& macrosToCheck = macrosInBin;
  // sort(macrosToCheck.begin(), macrosToCheck.end(),
  // FPCheckOrderByDecreasingArea(hgraph));

  vector<unsigned> tmpFPmarks;
  vector<unsigned> oldNodeIDs;
  vector<Point> oldNodeLocations;
  // check the largest macros
  unsigned nodesSelected = 0;
  for (unsigned c = 0; bin.getNumCells() < 10000 && c < macrosToCheck.size();
       c++) {
    unsigned idx = macrosToCheck[c];

    double nodeArea = hgraph.getNodeHeight(idx) * hgraph.getNodeWidth(idx);

    if (nodeArea < 0.07 * bin.getTotalCellArea()) continue;

    nodesSelected++;

    oldNodeIDs.push_back(idx);
    oldNodeLocations.push_back(_placement[idx]);

    if (nodeNeedsFP(idx)) continue;

    markNodeNeedsFP(idx);
    tmpFPmarks.push_back(idx);
  }

  if (nodesSelected) {
    Verbosity origVerb = _params.verb;
    _params.verb = Verbosity("0 0 0");
    lookAheadFPSuccess =
        floorplanBin(bin, /*minWL=*/false, /*lookAheadFP=*/true);
    _params.verb = origVerb;
  }

  // restore original node locations
  for (unsigned i = 0; i < oldNodeIDs.size(); i++) {
    _placement[oldNodeIDs[i]] = oldNodeLocations[i];
  }

  if (!lookAheadFPSuccess) {
    // partitioner gave us bad bins, we shall go up and fix this problem
    // ourselves
    CapoBin* mergedBin = mergeBinWithSibling(bin);
    if (mergedBin != NULL) {
      // successfully merged bin with its sibling

      unsigned numMacros = 0;
      for (unsigned c = 0; c < mergedBin->getNumCells(); c++) {
        unsigned idx = mergedBin->getCellIds()[c];
        if (hgraph.isNodeMacro(idx)) numMacros++;
      }

      if (numMacros <= 50 || (!_params.selectiveFloorplanning &&
                              numMacros <= _params.scampiThresh)) {
        // if bin doesn't have too many macros, mark all for FP to not create
        // problems
        markNodesNeedFP(mergedBin->getCellIds());
      } else if (bin._sibling != NULL) {
        // check FP conditions on sibling bin as well (mark nodes for FP)
        double splitPt;
        bool desiredDir = getDefaultSplitDir(*(bin._sibling));
        (void)bin._sibling->fpCondMet(desiredDir, splitPt);
      }

      mergedBin->_needFP = true;
      updateCell2Bin(mergedBin);
      newLayout.push_back(mergedBin);
      bin._needFP = true;
      bin._floorplanned = false;

      if (_params.verb.getForActions() > 0)
        cout << endl
             << "Unfloorplannable bin - created Hierarchical Bin "
             << mergedBin->getName() << " for FPing" << endl;

      bool foundSibling = false;
      CapoBin* bin2 = bin._sibling;
      if (bin._needFP) {
        for (unsigned b2 = 0; b2 < curLayout.size(); b2++) {
          if (curLayout[b2] != 0)
            if (curLayout[b2] == bin2) {
              // bin2->unLinkNeighbors();
              delete bin2;
              curLayout[b2] = 0;
              foundSibling = true;
              break;
            }
        }
        // bin.unLinkNeighbors();
      }
      if (!foundSibling) bin2->_sibling = bin2;
      delete curLayout[b];
      curLayout[b] = 0;

      return false;
    }
  }

  // unmark nodes temporarily marked for lookAheadFP,
  for (unsigned i = 0; i < tmpFPmarks.size(); i++)
    unmarkNodeNeedsFP(tmpFPmarks[i]);

  return true;
}
