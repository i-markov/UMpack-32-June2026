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

// Created: Sep 11, 1997 by Igor Markov & Andy Caldwell

// CHANGES
// 970923 aec moving to new Capo structure: major changes include removal pf
//       padBin, addition of a persistent slice, move away from
//       the 'changeNetlist' 'changeLayout' design into a unified
//       'refineBin' structure.
// 971113 aec added class CapoPlacer::Parameters
// 971118 aec added SA multi-start and multi-level params to
//           CapoPlacer::Paremeters
// 980215 aec moved parameter code to capoParams.cxx
// 980305 ilm removed fixedConstraint handling for Partitioner::Parameters
// 980401 ilm split into refineBin.cxx, splitBin.cxx and the rest
// 990202 aec capoPlacer now operates on an RBPlacement and an HGraphWDims,
//       rather than DB
// 990319 ilm the main Capo loop now copies all placed bins to the new
//            level and stops when there are no bins to refine

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkio.h"
#include "ABKCommon/abkfunc.h"
#include "ABKCommon/abkrand.h"
#include "ABKCommon/abktimer.h"
#include "ABKCommon/infolines.h"
#include "AnalytPl/analytPl.h"
#include "ClusteredHGraph/clustHGraph.h"
#include "Constraints/allConstr.h"
#ifdef USEFLUTE
#include "Flute/flute.h"
#endif
#include "GeomTrees/FastSteiner/global.h"
#include "GeomTrees/FastSteiner/random.h"
#include "GeomTrees/FastSteiner/stnr1.h"
#include "ParquetDBFromRBP/ParquetDBFromRBP.h"
#ifdef USEPATOMA
#include "Patoma/patoma.h"
#endif
#include "PlaceEvals/pinPlEval.h"
#include "PlaceEvals/placeEvals.h"
#include "Placement/xyDraw.h"
#include "RBPlace/RBPlacement.h"
#include "Stats/trivStats.h"
#ifdef USETRAFFIC
#include "Traffic/traffic.h"
#endif

#include "capoBin.h"
#include "capoPlacer.h"
#include "thunderPart.h"

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::fstream;
using std::ifstream;
using std::max;
using std::min;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::pair;
using std::setw;
using std::string;
using std::vector;

namespace {

double medianUnsignedValues(vector<unsigned>& values, unsigned count) {
  if (count == 0) return 0.;

  if ((count % 2) == 0) {
    unsigned middlehigh = count / 2;
    unsigned middlelow = middlehigh - 1;
    nth_element(values.begin(), values.begin() + middlelow, values.end());
    nth_element(values.begin() + middlehigh, values.begin() + middlehigh,
                values.end());
    return 0.5 * (static_cast<double>(values[middlelow]) +
                  static_cast<double>(values[middlehigh]));
  }

  unsigned middle = count / 2;
  nth_element(values.begin(), values.begin() + middle, values.end());
  return static_cast<double>(values[middle]);
}

}  // namespace

double CapoPlacer::capoTimer::FMTime = 0.;
double CapoPlacer::capoTimer::MLTime = 0.;
#ifdef USEHMETIS
double CapoPlacer::capoTimer::HMetisTime = 0.;
#endif
double CapoPlacer::capoTimer::BBTime = 0.;
double CapoPlacer::capoTimer::ECOTime = 0.;
double CapoPlacer::capoTimer::SmTime = 0.;
double CapoPlacer::capoTimer::SmPlaceProb = 0.;
double CapoPlacer::capoTimer::SetupTime = 0.;
double CapoPlacer::capoTimer::FloorplanTime = 0.;
double CapoPlacer::capoTimer::FloorClusterTime = 0.;
double CapoPlacer::capoTimer::FloorAnnealTime = 0.;
double CapoPlacer::capoTimer::MacroOverlapTime = 0.;
double CapoPlacer::capoTimer::WLCalcTime = 0.;
double CapoPlacer::capoTimer::FeedbackProcessing = 0.;
double CapoPlacer::capoTimer::ThunderTime = 0.;
double CapoPlacer::capoTimer::SORTime = 0.;
double CapoPlacer::capoTimer::PointSetTime = 0.;
double CapoPlacer::capoTimer::FastPlaceMoveTime = 0.;
double CapoPlacer::capoTimer::BloatTime = 0.;
double CapoPlacer::capoTimer::CongMapTime = 0.;
double CapoPlacer::capoTimer::HGraphRebuildTime = 0.;
#ifdef USEFLOW
double CapoPlacer::capoTimer::FlowTime = 0.;
#endif
unsigned CapoPlacer::FloorplanStats::numInstances = 0;
unsigned CapoPlacer::FloorplanStats::numSuccessfulInstances = 0;
unsigned CapoPlacer::FloorplanStats::numFailedInstances = 0;
unsigned CapoPlacer::FloorplanStats::instancesMaxMacroNum = 0;
unsigned CapoPlacer::FloorplanStats::failedInstancesMaxMacroNum = 0;
unsigned CapoPlacer::FloorplanStats::totalFPcells = 0;
unsigned CapoPlacer::FloorplanStats::totalFPmacros = 0;
unsigned CapoPlacer::FloorplanStats::totalFPobstacles = 0;
vector<unsigned> CapoPlacer::FloorplanStats::cellsInInstance;
vector<unsigned> CapoPlacer::FloorplanStats::macrosInInstance;
vector<unsigned> CapoPlacer::FloorplanStats::obstaclesInInstance;
unsigned CapoPlacer::EcoStats::ecoPartsUsed = 0;
unsigned CapoPlacer::EcoStats::ecoPartsRejected = 0;
unsigned CapoPlacer::EcoStats::ecoFPsUsed = 0;

// forward declarations of helper functions
static double getRealNodeWidth(unsigned index,
                               const PlacementWOrient& placement,
                               const HGraphWDimensions& hgraph);
static double getRealNodeHeight(unsigned index,
                                const PlacementWOrient& placement,
                                const HGraphWDimensions& hgraph);

void CapoPlacer::printNetStat() {
  int first = 1;
  double len = 0.0, mlen = 0.0, ncost = 0.0;
  itHGFEdgeGlobal edge;
  double alen[10];
  BBox coreBBox;
  vector<const RBPCoreRow*> rows;
  vector<CROffset> endOffsets;

  itRBPCoreRow r;
  for (r = _rbplace.coreRowsBegin(); r != _rbplace.coreRowsEnd(); ++r) {
    rows.push_back(&(*r));
    const RBPSubRow& lastSR = *(r->subRowsBegin() + (r->getNumSubRows() - 1));
    CROffset tmpOffset(r->getNumSubRows() - 1, lastSR.getNumSites() - 1);
    endOffsets.push_back(tmpOffset);
    coreBBox += Point((*r)[0].getXStart(), r->getYCoord());
    coreBBox += Point((*r)[endOffsets.back().first].getXEnd(),
                      r->getYCoord() + r->getHeight());
  }

  // bzero(alen, sizeof(double)*10);
  memset(alen, 0, sizeof(double) * 10);

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  for (edge = hgraph.edgesBegin(); edge != hgraph.edgesEnd(); edge++) {
    itHGFNodeLocal node;
    BBox box;
    box.clear();
    for (node = (*edge)->nodesBegin(); node != (*edge)->nodesEnd(); node++) {
      unsigned nodeId = (*node)->getIndex();
      box += _placement[nodeId];
    }
    len = fabs(box.xMax - box.xMin) + fabs(box.yMax - box.yMin);
    if (len > mlen || first) {
      mlen = len;
      first = 0;
    }
  }

  printf("Max edge len is %.2f\n", mlen);

  ncost = coreBBox.xMax - coreBBox.xMin;
  ncost += coreBBox.yMax - coreBBox.yMin;

  for (edge = hgraph.edgesBegin(); edge != hgraph.edgesEnd(); edge++) {
    itHGFNodeLocal node;
    BBox box;
    box.clear();
    for (node = (*edge)->nodesBegin(); node != (*edge)->nodesEnd(); node++) {
      unsigned nodeId = (*node)->getIndex();
      box += _placement[nodeId];
    }

    len = fabs(box.xMax - box.xMin) + fabs(box.yMax - box.yMin);
    //    printf("%.2f %.2f\n", len, ncost);
    int bin = (int)(len / (ncost / 9));
    //    printf("%d bin\n\n", bin);
    alen[bin]++;
  }

  for (int i = 0; i < 10; i++) {
    printf("%d %.2f\n", i, alen[i]);
  }
}

// extern int AMB;

void CapoPlacer::updateCell2Bin(const CapoBin* bin) {
  vector<unsigned>::const_iterator cIt;
  for (cIt = bin->cellIdsBegin(); cIt != bin->cellIdsEnd(); cIt++) {
    _cellToBinMap[*cIt] = bin;
  }
}

void CapoPlacer::undoLayer(vector<CapoBin*>& newLayout,
                           vector<CapoBin*>& extraLayout) {
  int i = 0;
  int csize = _curLayout->size();

  //  cerr << "layout of size " << csize << endl;
  if (_params.feedback) {
    newLayout.clear();
    while (i < csize) {
      vector<CapoBin*> mix;
      mix.clear();
      if (i != (csize - 1)) {
        mix.push_back((*_curLayout)[i]);
        mix.push_back((*_curLayout)[i + 1]);
        CapoBin* bin1 = (*_curLayout)[i];
        CapoBin* bin2 = (*_curLayout)[i + 1];
        int dirc = 2;
        if (equalDouble(bin1->getBBox().xMin, bin2->getBBox().xMin) &&
            equalDouble(bin1->getBBox().xMax, bin2->getBBox().xMax)) {
          dirc = 1;
        }
        if (equalDouble(bin1->getBBox().yMin, bin2->getBBox().yMin) &&
            equalDouble(bin1->getBBox().yMax, bin2->getBBox().yMax)) {
          dirc = 0;
        }

        /*
          if((dirc==1 && !(bin1->getBBox().yMin == bin2->getBBox().yMax ||
          bin1->getBBox().yMax == bin2->getBBox().yMin)) ||
          (dirc==0 && !(bin1->getBBox().xMin == bin2->getBBox().xMax ||
          bin1->getBBox().xMax == bin2->getBBox().xMin)))
          {
          cout<<bin1->getBBox()<<" : "<<bin2->getBBox()<<endl;
          //abkfatal(0, "bug in sheriefs code\n");
          }
        */

        if ((dirc == 1 || dirc == 0)) {
          if (bin2 == bin1->_sibling) {
            // CapoBin *bin=new CapoBin(CapoBin(mix, dirc));
            CapoBin* bin = new CapoBin(mix, dirc);
            newLayout.push_back(bin);
            updateCell2Bin(bin);
            bin1->_sibling->_sibling = bin1->_sibling;
            bin2->_sibling->_sibling = bin2->_sibling;
            bin1->_sibling = bin1;
            bin2->_sibling = bin2;
            // bin1->unLinkNeighbors();
            // bin2->unLinkNeighbors();
            delete (bin1);
            delete (bin2);
            i += 2;
          } else {
            extraLayout.push_back(bin1);
            i++;
          }
        } else {
          extraLayout.push_back(bin1);
          // newLayout.push_back(bin1);
          i++;
        }
      } else {
        CapoBin* bin1 = (*_curLayout)[i];
        extraLayout.push_back(bin1);
        i++;
      }
    }

    if (_params.verb.getForActions() > 0) {
      cout << "HPWL after undoing layer for feedback: " << estimateWL() << endl;
    }
  }
}

void CapoPlacer::capoPreMainLoop(unsigned& maxBins, unsigned& saveAtBins,
                                 bool& haveUnplacedBins, bool& doVCuts,
                                 bool& saveBins) {
  if (_params.verb.getForActions() > 0)
    cerr << endl << "Launching the UCLA/Michigan Capo placer" << endl;

  setupAndCheck();

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  if (_params.wtCut || _params.doFastPlaceMoves) {
    Timer psTimer;
    // reserve the size of the arrays to num edges
    _pointSet_fixed.resize(hgraph.getNumEdges());
    _pointSet_movable.resize(hgraph.getNumEdges());
    buildPointSetsFromScratch();
    psTimer.stop();
    CapoPlacer::capoTimer::PointSetTime += psTimer.getUserTime();
  }

  // printNetStat();
  if (_params.mixedSize) {
    float maxCoreHeight = _rbplace.getMaxHeightCoreRow();
    float maxCellWidth = std::numeric_limits<float>::max();
    const_cast<HGraphWDimensions&>(hgraph).markNodesAsMacro(maxCellWidth,
                                                            maxCoreHeight);

    if (_params.allMacroSoft) {
      cout << "-allMacroSoft option found," << endl
           << " it forces aspect ratios of all macros to be in " << "["
           << _params.minAR << "," << _params.maxAR << "]" << endl;
      const_cast<HGraphWDimensions&>(hgraph).markMacrosAsSoft(_params.minAR,
                                                              _params.maxAR);
    } else if (!_params.softFileName.empty()) {
      cout << "reading soft block information from " << _params.softFileName
           << endl;
      const_cast<HGraphWDimensions&>(hgraph).markSelectedAsSoftFromFile(
          _params.softFileName);
    }
  }

  maxBins = _params.stopAtBins ? _params.stopAtBins : UINT_MAX;

  saveAtBins = _params.saveAtBins ? _params.saveAtBins : UINT_MAX;

  haveUnplacedBins = true;
  doVCuts = true;
  saveBins = true;
}

bool CapoPlacer::capoMainLoopTest(bool haveUnplacedBins, unsigned maxBins) {
  return (haveUnplacedBins && _curLayout->size() < maxBins);
}

void CapoPlacer::capoMainLoop(unsigned& saveAtBins, bool& haveUnplacedBins,
                              bool& saveBins, bool& doVCuts,
                              vector<CapoBin*>& newLayout,
                              vector<CapoBin*>& extraLayout,
                              vector<CapoBin*>& bestLayout) {
  // <aaronnn> spin random numbers
  if (_mainLayerNum < _params.layerSeeds.size() && _mainLayerNum > 0)
    ClusteredHGraph::spinRandomNumber(_params.layerSeeds[_mainLayerNum]);

  AllowedPartDir partDir;

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  // if bloating, update the congestion information
  if (_params.useBloating) {
    updateCongestionEstimation();
  }

  // build new congestion map if necessary
  if (_params.splitterParams.congestionShifting /* && _mainLayerNum > 1 */ &&
      (!_params.splitterParams.eco || _params.splitterParams.ecoShift)) {
    Timer congMapTimer;

    unsigned gridCellsPerSide =
        2 * static_cast<unsigned>(ceil(sqrt(static_cast<double>(
                _curLayout->size() + _placedBins.size())))) +
        1;
    unsigned maxVGridCells = _rbplace.getNumCoreRows();
    if (_params.splitterParams.eco && _params.splitterParams.ecoShift) {
      _ispd04CongMapPlacement = new PlacementWOrient(_placement.size());
      _ispd04CongMap = new ISPD04CongMap(
          _rbplace, *_ispd04CongMapPlacement, gridCellsPerSide,
          std::min(gridCellsPerSide, maxVGridCells));
    } else {
      _ispd04CongMap =
          new ISPD04CongMap(_rbplace, _placement, gridCellsPerSide,
                            std::min(gridCellsPerSide, maxVGridCells));
    }
    _ispd04CongMap->setObstacles(&fixedObstacles);
    congMapTimer.stop();
    CapoPlacer::capoTimer::CongMapTime += congMapTimer.getUserTime();
  }

#ifdef USEFLOW
  // if using flows, solve the flow within the allowed bins
  if (_params.splitterParams.flowECO) {
    Timer flowTimer;

    // equalize the orients and the positions of fixed objects
    _oraclePlacement = _placement;

    vector<unsigned> movableCells;
    vector<BBox> cellBBoxes;

    for (unsigned i = hgraph.getNumTerminals(); i < hgraph.getNumNodes(); ++i) {
      const CapoBin* bin = _cellToBinMap[i];
      if (bin != NULL) {
        movableCells.push_back(i);
        cellBBoxes.push_back(bin->getBBox());
      }
    }

    if (_params.verb.getForActions()) {
      cout << "Generating flow-based solution" << endl;
    }
    _rbplace.doUnconstrainedXFlow(_placement, hgraph, movableCells, cellBBoxes,
                                  _oraclePlacement, _params.maxMem);
    _rbplace.doUnconstrainedYFlow(_placement, hgraph, movableCells, cellBBoxes,
                                  _oraclePlacement, _params.maxMem);
    if (_params.verb.getForActions()) {
      cout << endl;
    }

    flowTimer.stop();
    CapoPlacer::capoTimer::FlowTime += flowTimer.getUserTime();
  }
#endif

  if (_curLayout->size() >= saveAtBins && saveBins) {
    // take a snapshot of bins
    saveBinsCopy();
    saveBins = false;
  }

  if (_mainLayerNum < 3)
    partDir = HAndV;
  else if (_mainLayerNum == 3 || doVCuts) {
    partDir = VOnly;
    doVCuts = false;  // next time will be horizontal
  } else {
    partDir = HOnly;
    doVCuts = true;
  }

  if (_params.capoNE || _params.useQuad) {
    Timer sorTimer;

    solveQuadraticMinSoln();
    //_rbplace.updatePlacementWOri(_placement);

    sorTimer.stop();

    CapoPlacer::capoTimer::SORTime += sorTimer.getUserTime();

    // Uncomment to plot every SOR layer
    // static unsigned numSORs=0;
    // stringstream fname;
    // fname << "out" << static_cast<int>(numSORs++);
    //_rbplace.updatePlacementWOri(_placement);
    // Plotters::RBPlacePlotter::Parameters params;
    // params.plotNodes = true;
    //_rbplace.saveAsPlot(params,fname.str().c_str());
  }

  // <aaronnn> save intermediate placement before partitioning
  if (_params.prePartSavePl && _mainLayerNum > 0) {
    char plFileName[] = "out-tmp.pl";
    _rbplace.updatePlacementWOri(_placement);
    cout << "Saving " << plFileName << " ..." << endl;
    _rbplace.savePlacement(plFileName);
  }

  bool dummy;
  _feedbackLayerNum = 0;
  haveUnplacedBins = doOneLayer(partDir, dummy);

  /*
    if(_curLayout->size() >= maxBins)
    {
    _rbplace.updatePlacementWOri(_placement);
    _rbplace.saveAsPlot("plotNodes",-DBL_MAX,DBL_MAX,-DBL_MAX,DBL_MAX,"out");
    exit(0);
    }
  */
  // if replaceSmallBins == AtEveryLayer, the placed bins
  // will still be in curLayout.  Otherwise, they will have
  // been moved into _placedBins
  //      cerr<<"Attempting to write to OA..."<<endl;
  //      writePlacement(singletonCurrDes::getCurrDesign(),hgraph,_placement);

  Timer feedback;

  bool atBot = atBottomLayers();
  if (_params.verb.getForActions() > 0 && _params.feedback &&
      _curLayout->size() > 0 && atBot) {
    cout << "Skipping feedback at the bottom layers" << endl;
  }

  bool allECO = true;
  for (unsigned i = 0; i < _curLayout->size(); ++i) {
    if (!(*_curLayout)[i]->isEcoAllowed()) {
      allECO = false;
      break;
    }
  }

  if (_params.verb.getForActions() > 0 && _params.feedback &&
      _curLayout->size() > 0 && _params.splitterParams.eco && allECO) {
    cout << "Skipping feedback during ECO" << endl;
  }

  if (_params.verb.getForActions() > 0 && _params.feedback &&
      _curLayout->size() == 0) {
    cout << "Skipping feedback due to lack of bins" << endl;
  }

  feedback.stop();
  CapoPlacer::capoTimer::FeedbackProcessing += feedback.getUserTime();

  if (_params.feedback && !atBot && !allECO && _curLayout->size() > 0) {
    feedback.start();
    unsigned origNumMLSets = _params.splitterParams.numMLSets;
    double origBFactorAdjCut = _params.splitterParams.bFactorAdjustCutoff;
    MLPartParams::VcyclType origVcycling = _params.mlParams.Vcycling;

    if (_params.acc) {
      _params.splitterParams.numMLSets = 1;
      _params.mlParams.Vcycling = MLPartParams::NoVcycles;
      _params.splitterParams.bFactorAdjustCutoff = 0.;
    } else {
      _params.splitterParams.numMLSets = 2;
    }

    double bestBBox = estimateWL();
    if (_params.verb.getForActions() > 0)
      cout << "Feedback Iter 0:" << " HPWL: " << bestBBox << endl;
    bestLayout.clear();
    for (unsigned k = 0; k < _curLayout->size(); k++) {
      // CapoBin *nbin=new CapoBin((*_curLayout)[k]);
      CapoBin* nbin = (*_curLayout)[k];
      bestLayout.push_back(nbin);
    }

    int iter = _params.termIter;
    extraLayout.clear();

    feedback.stop();
    CapoPlacer::capoTimer::FeedbackProcessing += feedback.getUserTime();

    while (iter > 0) {
      feedback.start();
      if (_params.verb.getForActions() > 0)
        cout << "Undoing layer for feedback" << endl;
      undoLayer(newLayout, extraLayout);
      _curLayout = &newLayout;
      _layout[0] = newLayout;
      _layout[1] = newLayout;

      if (_params.acc && iter == 1) {
        // reinstate parameters for final iteration of accelerated feedback
        _params.splitterParams.numMLSets = origNumMLSets;
        _params.splitterParams.bFactorAdjustCutoff = origBFactorAdjCut;
        _params.mlParams.Vcycling = origVcycling;
      }

      bool fpdSomething = false;

      feedback.stop();
      CapoPlacer::capoTimer::FeedbackProcessing += feedback.getUserTime();

      if (_curLayout->size() > 0) doOneLayer(partDir, fpdSomething);
      feedback.start();

      if (fpdSomething) {
        for (unsigned k = 0; k < extraLayout.size(); ++k) {
          extraLayout[k]->resetRows();
          // reset the core rows and row offsets
          // and recalculate the BBox (may have shrunk)
          extraLayout[k]->resetRows();
          // the BBox may have changed so reposition
          // cells that are left in the fpd bins
          updateInfoAboutBin(extraLayout[k]);
          // cells may have moved so keep the
          // point sets up to date
        }
      }

      double bboxWL4 = estimateWL();
      if (_params.verb.getForActions() > 0)
        cout << "Feedback Iter: " << _params.termIter - iter + 1
             << " HPWL: " << bboxWL4 << endl;
      if (_params.contMode == 2) {
        if (bboxWL4 > bestBBox) break;
        bestBBox = bboxWL4;
        bestLayout.clear();
        for (unsigned k = 0; k < _curLayout->size(); k++) {
          // CapoBin *nbin=new CapoBin((*_curLayout)[k]);
          CapoBin* nbin = (*_curLayout)[k];
          bestLayout.push_back(nbin);
        }
      }

      if (_params.contMode == 1 ||
          (_params.contMode == 3 && bboxWL4 < bestBBox)) {
        bestBBox = bboxWL4;
        bestLayout.clear();
        for (unsigned k = 0; k < _curLayout->size(); k++) {
          // CapoBin *nbin=new CapoBin((*_curLayout)[k]);
          CapoBin* nbin = (*_curLayout)[k];
          bestLayout.push_back(nbin);
        }
      }
      iter--;
      feedback.stop();
      CapoPlacer::capoTimer::FeedbackProcessing += feedback.getUserTime();
    }
    feedback.start();
    // TODO by sadya. bestLayout is shady. the way it should be implemented
    // causes memory leaks. need to figure this out properly. For now
    // the other modes of feeback except unconstrained do not work properly.

    if (_params.verb.getForActions() > 1)
      if (_params.contMode == 1)
        cout << "Unconstrained Controller Feedback" << endl;
      else if (_params.contMode == 2)
        cout << "Non-Increasing Controller Feedback" << endl;
      else if (_params.contMode == 3)
        cout << "Best Controller Feedback" << endl;

    _curLayout->clear();
    for (unsigned k = 0; k < bestLayout.size(); k++)
      _curLayout->push_back(bestLayout[k]);
    for (unsigned k = 0; k < extraLayout.size(); k++)
      _curLayout->push_back(extraLayout[k]);

    double bboxWL4 = estimateWL();
    if (_params.verb.getForMajStats() > 0)
      cout << " Best HPWL from feedback: " << bboxWL4 << endl;

    _layout[0] = *_curLayout;
    _layout[1] = *_curLayout;
    haveUnplacedBins = (_curLayout->size() != 0);

    // reinstate old parameters
    _params.splitterParams.numMLSets = origNumMLSets;
    _params.splitterParams.bFactorAdjustCutoff = origBFactorAdjCut;
    _params.mlParams.Vcycling = origVcycling;
    feedback.stop();
    CapoPlacer::capoTimer::FeedbackProcessing += feedback.getUserTime();
  }

  // <aaronnn> launch ThunderPart (global placement refinement)
  if (_params.thunder && !atBot && _curLayout->size() > 2) {
    Timer thunderTimer;

    cout << endl;
    ThunderPartParams thunderParams;
    thunderParams.verb = _params.verb;
    ThunderPart thunderPart(*_curLayout, hgraph, _placement, _rbplace,
                            thunderParams);

    double initialHPWL = thunderPart.getHPWL();

    if (_params.verb.getForActions() > 0)
      cout << "HPWL before: " << initialHPWL << endl;

    thunderPart.partitionize(
        unsigned(_params.thunderMoveMultiplier * _rbplace.getNumCells()));

    double finalHPWL = thunderPart.getHPWL();

    if (_params.verb.getForActions() > 0) {
      cout << "HPWL after: " << finalHPWL << endl;

      double deltaHPWL = finalHPWL - initialHPWL;
      if (deltaHPWL < 0) {
        deltaHPWL = -deltaHPWL;
        double pctDeltaHPWL = deltaHPWL / initialHPWL * 100;
        deltaHPWL = floor((deltaHPWL * 1000) + 0.5) / 1000;
        pctDeltaHPWL = floor((pctDeltaHPWL * 1000) + 0.5) / 1000;

        cout << "Improvement in HPWL: " << deltaHPWL << " (" << pctDeltaHPWL
             << " %)" << endl
             << endl;
      } else if (deltaHPWL < 1e-6 && deltaHPWL > 1e-6) {
        cout << "No change in HPWL" << endl << endl;
      } else {
        // should never happen
        double pctDeltaHPWL = deltaHPWL / initialHPWL * 100;
        deltaHPWL = floor((deltaHPWL * 1000) + 0.5) / 1000;
        pctDeltaHPWL = floor((pctDeltaHPWL * 1000) + 0.5) / 1000;

        cout << "Degradation in HPWL: " << deltaHPWL << " (" << pctDeltaHPWL
             << " %)" << endl
             << endl;
      }
    }

    // save new bin assignments
    const vector<vector<unsigned> >& resultBins = thunderPart.getResultBins();
    for (unsigned b = 0; b < _curLayout->size(); b++) {
      if ((*_curLayout)[b] == NULL) continue;

      (*_curLayout)[b]->resetCellIds(resultBins[b]);
      updateInfoAboutBin((*_curLayout)[b]);
    }
    _layout[0] = *_curLayout;
    _layout[1] = *_curLayout;
    thunderTimer.stop();
    CapoPlacer::capoTimer::ThunderTime += thunderTimer.getUserTime();

    double totalLayoutArea = 0;
    double totalOverfill = 0;
    double totalWhitespace = 0;
    double maxOverfill = 0;
    for (unsigned b = 0; b < _curLayout->size(); b++) {
      double ws = (*_curLayout)[b]->getWhitespace();
      double cap = (*_curLayout)[b]->getCapacity();
      double ovr = (*_curLayout)[b]->getOverfill();
      totalWhitespace += ws;
      totalLayoutArea += cap;
      totalOverfill += ovr;
      maxOverfill = max(maxOverfill, ovr / cap);
    }
    for (unsigned b = 0; b < _placedBins.size(); b++) {
      double ws = _placedBins[b]->getWhitespace();
      double cap = _placedBins[b]->getCapacity();
      double ovr = _placedBins[b]->getOverfill();
      totalWhitespace += ws;
      totalLayoutArea += cap;
      totalOverfill += ovr;
      maxOverfill = max(maxOverfill, ovr / cap);
    }

    cout << "Layer " << _mainLayerNum << " Ave Overfill after ThunderPart: "
         << 100 * (totalOverfill / totalLayoutArea) << "%" << endl;
    cout << "Layer " << _mainLayerNum
         << " Max Overfill after ThunderPart: " << maxOverfill * 100 << "%"
         << endl;
    cout << "Layer " << _mainLayerNum << " Ave Whitespace after ThunderPart: "
         << 100 * (totalWhitespace / totalLayoutArea) << "%" << endl
         << endl;
  }

  if (_params.netStatsFileName != "") {
    fstream netStatLog;
    netStatLog.open(_params.netStatsFileName.c_str(),
                    fstream::out | fstream::app);
    netStatLog << "Statistics for layer " << _mainLayerNum << " :" << endl;
    vector<HGFEdge*> netsCutThisLayer, netsCutBefore;
    collectNetStats(netsCutThisLayer, netsCutBefore);
    netStatLog << "\tNumber of nets cut this layer: " << netsCutThisLayer.size()
               << endl;
    netStatLog << "\tNumber of nets cut before: " << netsCutBefore.size()
               << endl;
    netStatLog << "\tNumber of total nets cut: "
               << netsCutThisLayer.size() + netsCutBefore.size() << endl;

    bool usePinLocs = true, useWts = true;
    netStatLog << "\tNets cut this layer: " << endl;
    for (vector<HGFEdge*>::iterator it = netsCutThisLayer.begin();
         it != netsCutThisLayer.end(); ++it) {
      unsigned thisNetIndex = (*it)->getIndex();
      string thisNetName = hgraph.getNetNameByIndex(thisNetIndex).c_str();
      double cutNet_length = pinPlEval::oneNetHPWL(
          _placement, hgraph, thisNetIndex, usePinLocs, useWts);
      netStatLog << "\t\tName: \"" << thisNetName
                 << "\"\tIndex: " << thisNetIndex << "\tHPWL: " << cutNet_length
                 << endl;
    }

    netStatLog << "\tNets cut before: " << endl;
    for (vector<HGFEdge*>::iterator it = netsCutBefore.begin();
         it != netsCutBefore.end(); ++it) {
      unsigned thisNetIndex = (*it)->getIndex();
      string thisNetName = hgraph.getNetNameByIndex(thisNetIndex).c_str();
      double cutNet_length = pinPlEval::oneNetHPWL(
          _placement, hgraph, thisNetIndex, usePinLocs, useWts);
      netStatLog << "\t\tName: \"" << thisNetName
                 << "\"\tIndex: " << thisNetIndex << "\tHPWL: " << cutNet_length
                 << endl;
    }

    netStatLog << endl << endl;
    netStatLog.close();
  }

  if (_params.verb.getForMajStats() > 1) {
    cout << "\nTotal cut nets: " << getTotalNetCut() << endl;
    vector<unsigned> externalCutNets, containedNets;
    getNetCutInfo(externalCutNets, containedNets);

    if (externalCutNets.size() == _curLayout->size()) {
      cout << "Per Bin Stats:  " << "<numCells><cellArea><cutEdges>" << endl;

      for (unsigned ijk = 0; ijk < externalCutNets.size(); ijk++) {
        cout << "Bin " << ijk << " has " << (*_curLayout)[ijk]->getNumCells()
             << " cells, " << (*_curLayout)[ijk]->getTotalCellArea()
             << " cell area, " << externalCutNets[ijk]
             << " cut nets, "
             //              <<containedNets[ijk] << " contained nets, "
             << (*_curLayout)[ijk]->getNumAdjacentCells() << " adj. cells and "
             << (*_curLayout)[ijk]->getNumAdjacentCellsWithDuplicates()
             << " adj. cells with duplicates" << endl;
      }
    }
  }

  if (_params.plotCutLines && _mainLayerNum < _params.cutLineLayerCap) {
    BBox layout = getRBP().getBBox(true);
    double xthresh = 0.05 * (layout.xMax - layout.xMin);
    double ythresh = 0.05 * (layout.yMax - layout.yMin);
    std::set<unsigned> seenbins;

    recordCutLayer(_mainLayerNum);

    for (unsigned k = 0; k < _curLayout->size(); ++k) {
      CapoBin* curbin = (*_curLayout)[k];
      if (curbin == NULL) continue;
      if (curbin->_needFP || curbin->_floorplanned) continue;
      // only want to plot each line once
      if (seenbins.find(curbin->_index) != seenbins.end()) continue;
      CapoBin* sib = (*_curLayout)[k]->_sibling;
      if (sib == NULL || sib == curbin) continue;
      if (sib->_needFP || sib->_floorplanned) continue;
      seenbins.insert(curbin->_index);
      seenbins.insert(sib->_index);

      if (sib->_bbox.xMin == curbin->_bbox.xMax) {
        const double dist = curbin->_bbox.yMax - curbin->_bbox.yMin;
        if (dist > ythresh || _params.cutLineLayerCap != UINT_MAX) {
          Point bot(curbin->_bbox.xMax, curbin->_bbox.yMin);
          Point top(curbin->_bbox.xMax, curbin->_bbox.yMax);
          recordCutLine(bot, top, sib->_index, curbin->_index);
        }
      } else if (sib->_bbox.xMax == curbin->_bbox.xMin) {
        const double dist = curbin->_bbox.yMax - curbin->_bbox.yMin;
        if (dist > ythresh || _params.cutLineLayerCap != UINT_MAX) {
          Point bot(curbin->_bbox.xMin, curbin->_bbox.yMin);
          Point top(curbin->_bbox.xMin, curbin->_bbox.yMax);
          recordCutLine(bot, top, sib->_index, curbin->_index);
        }
      } else if (sib->_bbox.yMin == curbin->_bbox.yMax) {
        const double dist = curbin->_bbox.xMax - curbin->_bbox.xMin;
        if (dist > xthresh || _params.cutLineLayerCap != UINT_MAX) {
          Point left(curbin->_bbox.xMin, curbin->_bbox.yMax);
          Point right(curbin->_bbox.xMax, curbin->_bbox.yMax);
          recordCutLine(left, right, sib->_index, curbin->_index);
        }
      } else  // sib->_bbox.yMax == curbin->_bbox.yMin
      {
        const double dist = curbin->_bbox.xMax - curbin->_bbox.xMin;
        if (dist > xthresh || _params.cutLineLayerCap != UINT_MAX) {
          Point left(curbin->_bbox.xMin, curbin->_bbox.yMin);
          Point right(curbin->_bbox.xMax, curbin->_bbox.yMin);
          recordCutLine(left, right, sib->_index, curbin->_index);
        }
      }
    }
  }

  ++_mainLayerNum;

  if (_params.doFastPlaceMoves && !haveUnplacedBins &&
      fastPlaceMovesMade == 0) {
    haveUnplacedBins = doFastPlaceMoves();
  }

  // destroy the congestion map if necessary
  if (_ispd04CongMap != NULL) {
    Timer congMapTimer;

    delete _ispd04CongMap;
    _ispd04CongMap = NULL;

    delete _ispd04CongMapPlacement;
    _ispd04CongMapPlacement = NULL;

    congMapTimer.stop();
    CapoPlacer::capoTimer::CongMapTime += congMapTimer.getUserTime();
  }
}

void CapoPlacer::collectNetStats(vector<HGFEdge*>& netsCutThisLayer,
                                 vector<HGFEdge*>& netsCutBefore) {
  netsCutThisLayer.clear();
  netsCutBefore.clear();
  vector<int> nodeBin;
  getBinMembership(nodeBin);

  bool usePinLocs = true;
  bool useWts = true;
  bool netCut = false;
  unsigned cutNet_count = 0;
  double cutNet_length = 0;
  double notCutNet_length = 0;
  unsigned notCutNet_count = 0;
  const unsigned UNASSIGNED = std::numeric_limits<unsigned>::max();

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  // for every net in the design
  for (itHGFEdgeLocal netIterator = hgraph.edgesBegin();
       netIterator != hgraph.edgesEnd(); ++netIterator) {
    unsigned thisNetIndex = (*netIterator)->getIndex();
    netCut = false;
    unsigned binIDu = UNASSIGNED;
    itHGFNodeLocal node;
    // for every node on the net
    for (node = (*netIterator)->nodesBegin();
         node != (*netIterator)->nodesEnd(); node++) {
      unsigned nodeId = (*node)->getIndex();

      if (binIDu == UNASSIGNED) {
        binIDu = static_cast<unsigned>(nodeBin[nodeId]);
        continue;
      }
      unsigned nodeBinOfNodeIDu = static_cast<unsigned>(nodeBin[nodeId]);
      if (binIDu != nodeBinOfNodeIDu) {
        if (_layerCut.find(thisNetIndex) == _layerCut.end()) {
          netsCutThisLayer.push_back(*netIterator);
          _layerCut[thisNetIndex].push_back(binIDu);
          _layerCut[thisNetIndex].push_back(nodeBin[nodeId]);
          ++cutNet_count;
          cutNet_length += pinPlEval::oneNetHPWL(
              _placement, hgraph, thisNetIndex, usePinLocs, useWts);
          netCut = true;
          break;
        } else {
          // if binID is not in the _layerCut[thisNetIndex] vector
          if (find(_layerCut[thisNetIndex].begin(),
                   _layerCut[thisNetIndex].end(),
                   binIDu) == _layerCut[thisNetIndex].end()) {
            _layerCut[thisNetIndex].push_back(binIDu);
            netCut = true;
          }
          // if nodeBin[nodeId] is not in the _layerCut[thisNetIndex] vector
          if (find(_layerCut[thisNetIndex].begin(),
                   _layerCut[thisNetIndex].end(),
                   nodeBinOfNodeIDu) == _layerCut[thisNetIndex].end()) {
            _layerCut[thisNetIndex].push_back(nodeBinOfNodeIDu);
            netCut = true;
          }
          if (netCut) {
            ++cutNet_count;
            netsCutBefore.push_back(*netIterator);
            cutNet_length += pinPlEval::oneNetHPWL(
                _placement, hgraph, thisNetIndex, usePinLocs, useWts);
            break;
          }
        }
      }
    }  // end for every node on the net

    if (!netCut) {
      notCutNet_length += pinPlEval::oneNetHPWL(
          _placement, hgraph, thisNetIndex, usePinLocs, useWts);
      ++notCutNet_count;
    }  // end if thisNet is not cut
  }  // end for every net in the design
}

void CapoPlacer::capoPostMainLoop(Timer& capoTimer, RBPlacement& rbplace,
                                  bool& saveBins) {
  if (_params.saveBins) {
    cout << "Saving bin membership in bin_mem.dat" << endl;
    printBinMembership("bin_mem.dat");
    cout << "Saving hierarchical cell names cell_names.dat" << endl;
    printHierCellNames("cell_names.dat");
  }

  if (_params.replaceSmallBins == AtTheEnd) {
    Timer replaceTimer;
    cout << " Replacing " << _placedBins.size() << " Small Bins" << endl;
    replaceSmallBins();
    replaceTimer.stop();
    cout << " Replacing Small Bins took: " << replaceTimer.getUserTime()
         << endl;
  }

  if (_params.boost) printf("SQUARE COST %.2f\n", squareCost());
  capoTimer.stop();

  if (_params.verb.getForMajStats() > 0) {
    cout << "  CapoPlacer took : " << capoTimer << endl;
    double realtot = capoTimer.getUserTime();
    cout << "  Breakdown by component - \n";
    cout << "    MLPart took:        " << setw(11)
         << floor(capoTimer::MLTime * 100 + 0.5) / 100 << "sec ("
         << floor(capoTimer::MLTime / capoTimer.getUserTime() * 10000 + 0.5) /
                100
         << "%)" << endl;
    cout << "    FMPart took:        " << setw(11)
         << floor(capoTimer::FMTime * 100 + 0.5) / 100 << "sec ("
         << floor(capoTimer::FMTime / capoTimer.getUserTime() * 10000 + 0.5) /
                100
         << "%)" << endl;
#ifdef USEHMETIS
    if (capoTimer::HMetisTime > 0.) {
      cout << "    HMetis took:        " << setw(11)
           << floor(capoTimer::HMetisTime * 100 + 0.5) / 100 << "sec ("
           << floor(capoTimer::HMetisTime / capoTimer.getUserTime() * 10000 +
                    0.5) /
                  100
           << "%)" << endl;
    }
#endif
    cout << "    BBPart took:        " << setw(11)
         << floor(capoTimer::BBTime * 100 + 0.5) / 100 << "sec ("
         << floor(capoTimer::BBTime / capoTimer.getUserTime() * 10000 + 0.5) /
                100
         << "%)" << endl;
    if (_params.splitterParams.eco) {
      cout << "    ECO partitioning took: " << setw(8)
           << floor(capoTimer::ECOTime * 100 + 0.5) / 100 << "sec ("
           << floor(capoTimer::ECOTime / capoTimer.getUserTime() * 10000 +
                    0.5) /
                  100
           << "%)" << endl;
    }
    if (_params.useQuad) {
      cout << "    SOR took:   " << setw(19)
           << floor(capoTimer::SORTime * 100 + 0.5) / 100 << "sec ("
           << floor(capoTimer::SORTime / capoTimer.getUserTime() * 10000 +
                    0.5) /
                  100
           << "%)" << endl;
    }
    if (_params.thunder) {
      cout << "    ThunderPart took:   " << setw(11)
           << floor(capoTimer::ThunderTime * 100 + 0.5) / 100 << "sec ("
           << floor(capoTimer::ThunderTime / capoTimer.getUserTime() * 10000 +
                    0.5) /
                  100
           << "%)" << endl;
    }
    cout << "    SmPlace took:       " << setw(11)
         << floor(capoTimer::SmTime * 100 + 0.5) / 100 << "sec ("
         << floor(capoTimer::SmTime / capoTimer.getUserTime() * 10000 + 0.5) /
                100
         << "%)" << endl;
    cout << "    ProblemSetup took:  " << setw(11)
         << floor(capoTimer::SetupTime * 100 + 0.5) / 100 << "sec ("
         << floor(capoTimer::SetupTime / capoTimer.getUserTime() * 10000 +
                  0.5) /
                100
         << "%)" << endl;
    cout << "    SmPlProbSetup took: " << setw(11)
         << floor(capoTimer::SmPlaceProb * 100 + 0.5) / 100 << "sec ("
         << floor(capoTimer::SmPlaceProb / capoTimer.getUserTime() * 10000 +
                  0.5) /
                100
         << "%)" << endl;
    cout << "    Level Stats took:   " << setw(11)
         << floor(capoTimer::WLCalcTime * 100 + 0.5) / 100 << "sec ("
         << floor(capoTimer::WLCalcTime / capoTimer.getUserTime() * 10000 +
                  0.5) /
                100
         << "%)" << endl;
    if (_params.wtCut) {
      cout << "    PointSet Maintenance:" << setw(10)
           << floor(capoTimer::PointSetTime * 100 + 0.5) / 100 << "sec ("
           << floor(capoTimer::PointSetTime / capoTimer.getUserTime() * 10000 +
                    0.5) /
                  100
           << "%)" << endl;
    }
    if (_params.doFastPlaceMoves) {
      cout << "    Fast Place Moves:   " << setw(11)
           << floor(capoTimer::FastPlaceMoveTime * 100 + 0.5) / 100 << "sec ("
           << floor(capoTimer::FastPlaceMoveTime / capoTimer.getUserTime() *
                        10000 +
                    0.5) /
                  100
           << "%)" << endl;
    }
    if (_params.useBloating) {
      cout << "    Bloating took:      " << setw(11)
           << floor(capoTimer::BloatTime * 100 + 0.5) / 100 << "sec ("
           << floor(capoTimer::BloatTime / capoTimer.getUserTime() * 10000 +
                    0.5) /
                  100
           << "%)" << endl;
    }
    if (_params.splitterParams.congestionShifting) {
      cout << "    Congestion Maps took:" << setw(10)
           << floor(capoTimer::CongMapTime * 100 + 0.5) / 100 << "sec ("
           << floor(capoTimer::CongMapTime / capoTimer.getUserTime() * 10000 +
                    0.5) /
                  100
           << "%)" << endl;
    }
#ifdef USEFLOW
    if (_params.splitterParams.flowECO) {
      cout << "    Flows took:          " << setw(10)
           << floor(capoTimer::FlowTime * 100 + 0.5) / 100 << "sec ("
           << floor(capoTimer::FlowTime / capoTimer.getUserTime() * 10000 +
                    0.5) /
                  100
           << "%)" << endl;
    }
#endif
    if (_params.feedback) {
      cout << "    Feedback Processing took: " << setw(5)
           << floor(capoTimer::FeedbackProcessing * 100 + 0.5) / 100 << "sec ("
           << floor(capoTimer::FeedbackProcessing / capoTimer.getUserTime() *
                        10000 +
                    0.5) /
                  100
           << "%)" << endl;
    }
    if (_params.keepHGraphFraction < 1.) {
      cout << "    Rebuilding Hypergraphs took: " << setw(5)
           << floor(capoTimer::HGraphRebuildTime * 100 + 0.5) / 100 << "sec ("
           << floor(capoTimer::HGraphRebuildTime / capoTimer.getUserTime() *
                        10000 +
                    0.5) /
                  100
           << "%)" << endl;
    }
    if (_params.mixedSize && _rbplace.getNumMacros() > 0) {
      cout << "    Floorplanning took: " << setw(11)
           << floor(capoTimer::FloorplanTime * 100 + 0.5) / 100 << "sec ("
           << floor(capoTimer::FloorplanTime / capoTimer.getUserTime() * 10000 +
                    0.5) /
                  100
           << "%)" << endl;
      cout << "      Floorplanning Clustering took: " << setw(7)
           << floor(capoTimer::FloorClusterTime * 100 + 0.5) / 100 << "sec ("
           << floor(capoTimer::FloorClusterTime / capoTimer.getUserTime() *
                        10000 +
                    0.5) /
                  100
           << "%)" << endl;
      cout << "      Floorplanning Annealing took:  " << setw(7)
           << floor(capoTimer::FloorAnnealTime * 100 + 0.5) / 100 << "sec ("
           << floor(capoTimer::FloorAnnealTime / capoTimer.getUserTime() *
                        10000 +
                    0.5) /
                  100
           << "%)" << endl;
      if (_params.removeMacroOverlap) {
        cout << "      Macro Overlap Removal took:    " << setw(7)
             << floor(capoTimer::MacroOverlapTime * 100 + 0.5) / 100 << "sec ("
             << floor(capoTimer::MacroOverlapTime / capoTimer.getUserTime() *
                          10000 +
                      0.5) /
                    100
             << "%)" << endl;
      }
      cout << "      No. Floorplanning Instances: "
           << FloorplanStats::numInstances
           << " (successful: " << FloorplanStats::numSuccessfulInstances
           << ", failed: " << FloorplanStats::numFailedInstances << ")" << endl;
      cout << "      The largest floorplanning instance has "
           << FloorplanStats::instancesMaxMacroNum << " macros." << endl;
      cout << "      The largest failed instance has "
           << FloorplanStats::failedInstancesMaxMacroNum << " macros." << endl;
      cout << "      Average cells per instance: "
           << static_cast<double>(FloorplanStats::totalFPcells) /
                  static_cast<double>(FloorplanStats::numInstances)
           << endl;
      double cellMedian =
          medianUnsignedValues(FloorplanStats::cellsInInstance,
                               FloorplanStats::numInstances);
      cout << "      Median cells per instance: " << cellMedian << endl;
      cout << "      Average macros per instance: "
           << static_cast<double>(FloorplanStats::totalFPmacros) /
                  static_cast<double>(FloorplanStats::numInstances)
           << endl;
      double macroMedian =
          medianUnsignedValues(FloorplanStats::macrosInInstance,
                               FloorplanStats::numInstances);
      cout << "      Median macros per instance: " << macroMedian << endl;
      if (_params.FPEvadeObstacles) {
        cout << "      Average obstacles per instance: "
             << static_cast<double>(FloorplanStats::totalFPobstacles) /
                    static_cast<double>(FloorplanStats::numInstances)
             << endl;
        double obstacleMedian =
            medianUnsignedValues(FloorplanStats::obstaclesInInstance,
                                 FloorplanStats::numInstances);
        cout << "      Median obstacles per instance: " << obstacleMedian
             << endl;
      }
    }
    if (_params.splitterParams.eco) {
      cout << "    Eco Partitioning solutions used: " << EcoStats::ecoPartsUsed
           << endl;
      cout << "    Eco Partitioning solutions rejected: "
           << EcoStats::ecoPartsRejected << endl;
      cout << "    Eco Floorplanning solutions used: " << EcoStats::ecoFPsUsed
           << endl;
    }

    double tot = capoTimer::MLTime + capoTimer::FMTime +
#ifdef USEHMETIS
                 capoTimer::HMetisTime +
#endif
                 capoTimer::BBTime + capoTimer::ECOTime +
#ifdef USEFLOW
                 capoTimer::FlowTime +
#endif
                 capoTimer::SmTime + capoTimer::SetupTime +
                 capoTimer::SmPlaceProb + capoTimer::WLCalcTime +
                 capoTimer::FeedbackProcessing + capoTimer::FloorplanTime +
                 capoTimer::SORTime + capoTimer::ThunderTime +
                 capoTimer::PointSetTime + capoTimer::FastPlaceMoveTime +
                 capoTimer::BloatTime + capoTimer::CongMapTime +
                 capoTimer::HGraphRebuildTime;

    cout << "    Total runtime of measured components: "
         << floor(tot * 100 + 0.5) / 100 << "sec ("
         << floor(tot / realtot * 10000 + 0.5) / 100 << "%)" << endl;
    _params.maxMem->setPrintExtra(true);
    cout << endl << "  " << *_params.maxMem << endl;
    _params.maxMem->setPrintExtra(false);
  }

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  if (_params.netStatsFileName != "") {
    fstream netStatLog;
    netStatLog.open(_params.netStatsFileName.c_str(),
                    fstream::out | fstream::app);
    netStatLog << "Final wirelength of all nets: " << endl;
    for (itHGFEdgeLocal netIterator = hgraph.edgesBegin();
         netIterator != hgraph.edgesEnd(); ++netIterator) {
      bool usePinLocs = true, useWts = true;
      unsigned thisNetIndex = (*netIterator)->getIndex();
      string thisNetName = hgraph.getNetNameByIndex(thisNetIndex).c_str();
      double cutNet_length = pinPlEval::oneNetHPWL(
          _placement, hgraph, thisNetIndex, usePinLocs, useWts);
      netStatLog << "\t\tName: \"" << thisNetName
                 << "\"\tIndex: " << thisNetIndex << "\tHPWL: " << cutNet_length
                 << endl;
    }

    netStatLog << endl << endl;
    netStatLog.close();
  }

  if (_params.verb.getForMajStats() > 0) {
    cout << "  Final Average Relative Overfill: "
         << 100 * (_totalOverfill / _totalLayoutArea) << "%" << endl;
    cout << "  Max 1bin Overfill: " << 100 * _maxOverfill << "%" << endl;
  }

  rbplace.resetPlacementWOri(_placement);

  if (_params.mixedSize) {
    rbplace.alignCellsToRows();
  }

  if (_params.saveBinsFloorplan) {
    if (!saveBins) {
      cout << "Saving Bins as Floorplanning instance with base FP" << endl;
      saveBinsAsFloorplan("FP", _saveBins);
    } else {
      cout << " Could not save bins because capo run ended earlier" << endl;
    }
  }
}

CapoPlacer::CapoPlacer(RBPlacement& rbplace, const CapoParameters& params,
                       CongestionMaps* congestionMap)
    //                     const vector<char*>* altCellNames)
    : _params(params),
      _rbplace(rbplace),
      _coreBBox(/*BBoxFromRBPlace*/),
      _placement(static_cast<PlacementWOrient>(rbplace)),
      _oraclePlacement(
          (_params.splitterParams.WSOracle || _params.splitterParams.eco)
              ? static_cast<PlacementWOrient>(rbplace)
              : PlacementWOrient(0)),
      _congestionMap(congestionMap),
      _ispd04CongMap(NULL),
      _ispd04CongMapPlacement(NULL),
      _placed(rbplace.getFixed()),
      _totalLayerNum(0),
      _mainLayerNum(0),
      _feedbackLayerNum(0),
      nextBinNum(0),
      _cellToBinMap(rbplace.getHGraph().getNumNodes(),
                    static_cast<const CapoBin*>(NULL)),
      _nodeSeen(rbplace.getHGraph().getNumNodes()),
      _edgeSeen(rbplace.getHGraph().getNumEdges()),
      _nodesNeedingFP(rbplace.getHGraph().getNumNodes(), false),
      partProbEdgesVisited(rbplace.getHGraph().getNumEdges()),
      thetoWeights(params.wtCut ? rbplace.getHGraph().getNumEdges() : 0,
                   vector<double>(3, 0.)),
      _numOrigNets(100, 0),
      _numEssentialNets(100, 0),
      _origNetPins(100, 0),
      _essentialNetPins(100, 0),
      _numProblemsOfSize(100, 0),
      _netHasBeenCut(rbplace.getHGraph().getNumEdges(), false),
      _numNotSolved(0),
      _numSmallPartProbs(0),
      //      _slicingTree(rbplace.getBBox(true)),
      constraints(_rbplace.constraints),
      regionUtilization(_rbplace.regionUtilization),
      groupRegionConstr(_rbplace.groupRegionConstr),
      cellToGrpMapping(_rbplace.cellToGrpMapping),
      _congestion(rbplace.getHGraph().getNumNodes(), 0.),
      _congestionCutoff(1.),
      fastPlaceMovesMade(0),
      _fastStPt(NULL),
      _fastStParent(NULL) {
  compilerCheck();

  Timer capoTimer;

  CapoPlacer::capoTimer::clearTimer();
  CapoPlacer::FloorplanStats::clearFPStats();
  CapoPlacer::EcoStats::clearEcoStats();

  _params.maxMem->update("Beginning of Capo");

  _params.maxMem->setPrintExtra(true);
  cout << endl << *_params.maxMem << endl;
  _params.maxMem->setPrintExtra(false);

  if (!_params.automatic) return;

  rbplace.clearStdCellsFromCoreRows();

#ifdef USEFLUTE
  if (_params.useFLUTE || _params.fastPlaceFLUTE) {
    Timer lutTime;
    readLUT();
    lutTime.stop();
    CapoPlacer::capoTimer::SetupTime += lutTime.getUserTime();
  }
#endif

  if (_params.useFastSt || _params.fastPlaceFastSt) {
    Timer fastStSetup;
    fastSteiner::stnr1_package_init(4 * 22, 0);
    _fastStPt = new fastSteiner::Point[4 * 22];
    _fastStParent = new long[4 * 22];
    fastStSetup.stop();
    CapoPlacer::capoTimer::SetupTime += fastStSetup.getUserTime();
  }

  unsigned numObstacles = _rbplace.getNumFixedCore() +
                          _rbplace.getNumFixedMacros() +
                          _rbplace.getHGraph().getNumTerminals();
  unsigned gridSize = max(
      1u, static_cast<unsigned>(ceil(sqrt(static_cast<double>(numObstacles)))));

  obstacleGrid = vector<vector<vector<unsigned> > >(
      gridSize, vector<vector<unsigned> >(gridSize));
  BBox coreBBox = _rbplace.getBBox(false);
  double gridXSize = coreBBox.getWidth() / static_cast<double>(gridSize);
  double gridYSize = coreBBox.getHeight() / static_cast<double>(gridSize);

  for (unsigned i = 0; i < _rbplace.getHGraph().getNumNodes(); ++i) {
    if (_rbplace.isFixed(i)) {
      unsigned angle = _placement.getOrient(i).getAngle();
      bool notRotated = angle == 0 || angle == 180;
      double height = notRotated ? _rbplace.getHGraph().getNodeHeight(i)
                                 : _rbplace.getHGraph().getNodeWidth(i);
      double width = notRotated ? _rbplace.getHGraph().getNodeWidth(i)
                                : _rbplace.getHGraph().getNodeHeight(i);

      BBox fixedObs(_placement[i].x, _placement[i].y, _placement[i].x + width,
                    _placement[i].y + height);

      fixedObs = coreBBox.intersect(fixedObs);

      if (equalDouble(fixedObs.getArea(), 0.)) continue;

      unsigned startGridCol = static_cast<unsigned>(
          floor((fixedObs.xMin - coreBBox.xMin) / gridXSize));
      unsigned startGridRow = static_cast<unsigned>(
          floor((fixedObs.yMin - coreBBox.yMin) / gridYSize));
      unsigned endGridCol = static_cast<unsigned>(
          ceil((fixedObs.xMax - coreBBox.xMin) / gridXSize));
      unsigned endGridRow = static_cast<unsigned>(
          ceil((fixedObs.yMax - coreBBox.yMin) / gridYSize));

      for (unsigned rIdx = startGridRow; rIdx < endGridRow; ++rIdx)
        for (unsigned cIdx = startGridCol; cIdx < endGridCol; ++cIdx)
          obstacleGrid[rIdx][cIdx].push_back(fixedObstacles.size());

      relevantObstacles.push_back(fixedObstacles.size());
      fixedObstacles.push_back(fixedObs);
    }
  }

  unsigned maxBins, saveAtBins;
  bool haveUnplacedBins, doVCuts, saveBins;

  //    int totalAMB1=0, totalAMB2=0;
  vector<CapoBin*> newLayout;
  vector<CapoBin*> extraLayout;
  vector<CapoBin*> bestLayout;

  capoPreMainLoop(maxBins, saveAtBins, haveUnplacedBins, doVCuts, saveBins);

  while (capoMainLoopTest(haveUnplacedBins, maxBins)) {
    capoMainLoop(saveAtBins, haveUnplacedBins, saveBins, doVCuts, newLayout,
                 extraLayout, bestLayout);
  }
  capoPostMainLoop(capoTimer, rbplace, saveBins);
}

// jflu: pointSet to speedup net traversal

void sortAndUniquefyPointSet(vector<PtSetPoint>& ptset) {
  unsigned len = ptset.size();
  unsigned i = 0, j = 1;

  sort(ptset.begin(), ptset.end());

  while (j < ptset.size()) {
    if (ptset[i].pt == ptset[j].pt) {
      ptset[i].count += ptset[j].count;
      --len;
    } else {
      ++i;
      ptset[i] = ptset[j];
    }
    ++j;
  }

  ptset.resize(len);
}

void CapoPlacer::buildPointSetsFromScratch() {
  const HGraphWDimensions& hgraph = getNetlistHGraph();

  for (unsigned eId = 0; eId < hgraph.getNumEdges(); ++eId) {
    _pointSet_fixed[eId].clear();
    _pointSet_movable[eId].clear();

    for (itHGFNodeLocal nItr = hgraph.getEdgeByIdx(eId).nodesBegin();
         nItr != hgraph.getEdgeByIdx(eId).nodesEnd(); ++nItr) {
      if (_placed[(*nItr)->getIndex()]) {
        _pointSet_fixed[eId].push_back(
            getPinLocation((*nItr)->getIndex(), eId).getGeomCenter());
      } else {
        _pointSet_movable[eId].push_back(
            getPinLocation((*nItr)->getIndex(), eId).getGeomCenter());
      }
    }

    // sort and uniquefy
    sortAndUniquefyPointSet(_pointSet_fixed[eId]);
    sortAndUniquefyPointSet(_pointSet_movable[eId]);
  }
}

void removePtFromPtSet(vector<PtSetPoint>& ptSet, const PtSetPoint& pt) {
  std::pair<vector<PtSetPoint>::iterator, vector<PtSetPoint>::iterator> iters =
      std::equal_range(ptSet.begin(), ptSet.end(), pt);

  abkfatal(iters.first != iters.second, "pt to remove not found in point set");
  abkfatal((iters.first + 1) == iters.second,
           "duplicates found in point set when removing");
  abkfatal(iters.first->count >= 1, "incorrect multiplicity when removing");

  if (iters.first->count > 1) {
    --iters.first->count;
  } else {
    ptSet.erase(iters.first);
  }
}

void addPtToPtSet(vector<PtSetPoint>& ptSet, const PtSetPoint& pt) {
  std::pair<vector<PtSetPoint>::iterator, vector<PtSetPoint>::iterator> iters =
      std::equal_range(ptSet.begin(), ptSet.end(), pt);

  if (iters.first == iters.second) {
    ptSet.insert(iters.first, pt);
  } else {
    abkfatal((iters.first + 1) == iters.second,
             "duplicates found in point set when adding");
    abkfatal(iters.first->count >= 1, "incorrect multiplicity when adding");
    ++iters.first->count;
  }
}

void CapoPlacer::updatePointSetsForMovedCell(unsigned cellId,
                                             const Point& oldPos,
                                             const Point& newPos) {
  const HGFNode& node = getNetlistHGraph().getNodeByIdx(cellId);

  if (_placed[cellId]) {
    for (itHGFEdgeLocal eItr = node.edgesBegin(); eItr != node.edgesEnd();
         ++eItr) {
      unsigned eId = (*eItr)->getIndex();

      Point offset = getPinOffset(cellId, eId).getGeomCenter();

      PtSetPoint oldPinPtSetPoint(oldPos + offset, 1);
      PtSetPoint newPinPtSetPoint(newPos + offset, 1);

      removePtFromPtSet(_pointSet_fixed[eId], oldPinPtSetPoint);
      addPtToPtSet(_pointSet_fixed[eId], newPinPtSetPoint);
    }
  } else {
    for (itHGFEdgeLocal eItr = node.edgesBegin(); eItr != node.edgesEnd();
         ++eItr) {
      unsigned eId = (*eItr)->getIndex();

      Point offset = getPinOffset(cellId, eId).getGeomCenter();

      PtSetPoint oldPinPtSetPoint(oldPos + offset, 1);
      PtSetPoint newPinPtSetPoint(newPos + offset, 1);

      removePtFromPtSet(_pointSet_movable[eId], oldPinPtSetPoint);
      addPtToPtSet(_pointSet_movable[eId], newPinPtSetPoint);
    }
  }
}

void CapoPlacer::updatePointSetsForNewlyPlacedCell(unsigned cellId,
                                                   const Point& oldPos,
                                                   const Point& newPos) {
  const HGFNode& node = getNetlistHGraph().getNodeByIdx(cellId);

  abkfatal(_placed[cellId], "wrong function called");

  for (itHGFEdgeLocal eItr = node.edgesBegin(); eItr != node.edgesEnd();
       ++eItr) {
    unsigned eId = (*eItr)->getIndex();

    // ugly hack to keep offsets consistent,
    // as offset can change based on the cell being placed or not
    _placed[cellId] = false;
    Point oldOffset = getPinOffset(cellId, eId).getGeomCenter();
    _placed[cellId] = true;
    Point newOffset = getPinOffset(cellId, eId).getGeomCenter();

    PtSetPoint oldPinPtSetPoint(oldPos + oldOffset, 1);
    PtSetPoint newPinPtSetPoint(newPos + newOffset, 1);

    removePtFromPtSet(_pointSet_movable[eId], oldPinPtSetPoint);
    addPtToPtSet(_pointSet_fixed[eId], newPinPtSetPoint);
  }
}

const vector<PtSetPoint>& CapoPlacer::getPointSet_fixed(unsigned netID) const {
  return _pointSet_fixed[netID];
}

const vector<PtSetPoint>& CapoPlacer::getPointSet_movable(
    unsigned netID) const {
  return _pointSet_movable[netID];
}

double CapoPlacer::estimateWL() {
  return _rbplace.evalHPWL(_placement, getNetlistHGraph(), true);
}

void CapoPlacer::replaceSmallBins() {
  vector<CapoBin*>::iterator bin;
  double prevWL, newWL;

  // to disable point set maintenance during the replacement
  bool oldWtCut = _params.wtCut;
  bool oldFastPlaceMoves = _params.doFastPlaceMoves;
  bool oldCongWS = _params.splitterParams.congestionShifting;
  _params.wtCut = false;
  _params.doFastPlaceMoves = false;
  _params.splitterParams.congestionShifting = false;

  if (_params.verb.getForMajStats() > 0) {
    prevWL = estimateWL();
    cout << "\nEst WL before re-placing small bins " << prevWL << endl;
  }

  vector<CapoBin*> placedBinsCpy = _placedBins;
  _placedBins.clear();

  for (bin = placedBinsCpy.begin(); bin != placedBinsCpy.end(); bin++)
    refineBin(**bin, *_curLayout, HAndV, 0);
  // curLayout is a dummy var here..nothing will be
  // added to it, as these are all small bins

  if (_params.verb.getForMajStats() > 0) {
    newWL = estimateWL();
    cout << "\nEst WL after re-placing small bins " << newWL << endl;
  }

  _params.wtCut = oldWtCut;
  _params.doFastPlaceMoves = oldFastPlaceMoves;
  _params.splitterParams.congestionShifting = oldCongWS;
}

void CapoPlacer::setupAndCheck() {
  const HGraphWDimensions& hgraph = getNetlistHGraph();
  vector<unsigned> coreCells;
  coreCells.reserve(hgraph.getNumNodes());
  vector<unsigned> obstacleCells;
  obstacleCells.reserve(hgraph.getNumNodes());
  const BBox coreBBox = _rbplace.getBBox(/*includeTerms=*/false);

  Orientation nOrient;  // the 'default' N orientation
  itHGFNodeGlobal cell;
  for (cell = hgraph.nodesBegin(); cell != hgraph.nodesEnd(); cell++) {
    unsigned cellId = (*cell)->getIndex();

    if (_rbplace.isCoreCell(cellId)) {
      if (_rbplace.isFixed(cellId)) continue;

      coreCells.push_back(cellId);

      _placement.setOrient(cellId, nOrient);
      abkassert(hgraph.getWeight(cellId) != 0, "core cell has 0 area ");
    } else  // must be fixed
    {
      if (!hgraph.isTerminal(cellId) && _params.verb.getForMajStats() > 3) {
        cout << "Non-Terminal is not a movable core cell" << endl;
        cout << "Cell " << cellId << " located at " << _placement[cellId]
             << endl;
      }
      abkassert(_rbplace.isFixed(cellId), "unfixed non-core cell");
    }

    if (_rbplace.isFixed(cellId)) {
      // <aaronnn> whatever is fixed inside corebox is a possible obstacle
      BBox cellBBox;
      cellBBox.add(_placement[cellId].x, _placement[cellId].y);
      cellBBox.add(
          _placement[cellId].x + getRealNodeWidth(cellId, _placement, hgraph),
          _placement[cellId].y + getRealNodeHeight(cellId, _placement, hgraph));

      BBox intersection = cellBBox.intersect(coreBBox);

      if (greaterThanDouble(intersection.getArea(), 0.)) {
        obstacleCells.push_back(cellId);
      }
    }
  }

  vector<const RBPCoreRow*> rows;
  vector<CROffset> endOffsets;
  itRBPCoreRow r;
  for (r = _rbplace.coreRowsBegin(); r != _rbplace.coreRowsEnd(); r++) {
    rows.push_back(&(*r));
    const RBPSubRow& lastSR = *(r->subRowsBegin() + (r->getNumSubRows() - 1));
    CROffset tmpOffset(r->getNumSubRows() - 1, lastSR.getNumSites() - 1);
    endOffsets.push_back(tmpOffset);
    _coreBBox += Point((*r)[0].getXStart(), r->getYCoord());
    _coreBBox += Point((*r)[endOffsets.back().first].getXEnd(),
                       r->getYCoord() + r->getHeight());
  }

  vector<CROffset> startOffsets(rows.size(), CROffset(0, 0));

  if (_params.verb.getForMajStats() > 1)
    cout << endl << " Core BBox: " << _coreBBox << endl;

  vector<CapoBin*> noNeighbors;
  _layout[0].push_back(new CapoBin(
      coreCells, rows.begin(), rows.end(), startOffsets, endOffsets,
      noNeighbors, noNeighbors, noNeighbors, noNeighbors, *this, "BIN"));

  if (_params.splitterParams.WSOracle)
    _layout[0].back()->setOracleCellIds(coreCells);

  if (_params.splitterParams.eco) _layout[0].back()->setEcoAllowed(true);

  // <aaronnn> give the new bin some obstacles
  _layout[0].back()->setObstacleCellIds(obstacleCells);

  CapoBin* newBin = _layout[0].back();
  Point center = newBin->getCenter();

  for (unsigned c = 0; c < coreCells.size(); c++) {
    _cellToBinMap[coreCells[c]] = newBin;
    _placement[coreCells[c]] = center;
  }

  _curLayout = &(_layout[0]);

  if (_params.verb.getForMajStats() > 4) printNetlistStats();

  double layoutArea = _rbplace.getSitesArea();
  double totalNodesArea = hgraph.getNodesArea();
  double util = 100 * totalNodesArea / layoutArea;
  cout << "Utilization of given design is " << util << " %" << endl
       << "  Total placeable cell area: " << totalNodesArea << endl
       << "  Total placeable site area: " << layoutArea << endl
       << endl;

  if (_params.usePredeterminedCutlines)
    readCutlineFile("preDeterminedCutlines.txt");
}

void CapoPlacer::updateInfoAboutBin(const CapoBin* bin) {
  Timer psTimer;
  vector<unsigned>::const_iterator cIt;
  Point binCenter = bin->getCenter();
  for (cIt = bin->cellIdsBegin(); cIt != bin->cellIdsEnd(); ++cIt) {
    _cellToBinMap[*cIt] = bin;
    if (_placement[*cIt] != binCenter) {
      if (_params.wtCut || _params.doFastPlaceMoves) {
        updatePointSetsForMovedCell(*cIt, _placement[*cIt], binCenter);
      }
      _placement[*cIt] = binCenter;
    }
  }
  psTimer.stop();
  if (_params.wtCut || _params.doFastPlaceMoves) {
    CapoPlacer::capoTimer::PointSetTime += psTimer.getUserTime();
  }
}

void CapoPlacer::updatePlInfoAboutBin(const CapoBin* bin) {
  vector<unsigned>::const_iterator cIt;
  Point binCenter = bin->getCenter();
  for (cIt = bin->cellIdsBegin(); cIt != bin->cellIdsEnd(); cIt++) {
    _placement[*cIt] = binCenter;
    _rbplace.updatePlacement(*cIt, binCenter);
  }
}

void CapoPlacer::updateMapInfoAboutBin(const CapoBin* bin) {
  vector<unsigned>::const_iterator cIt;
  for (cIt = bin->cellIdsBegin(); cIt != bin->cellIdsEnd(); cIt++) {
    _cellToBinMap[*cIt] = bin;
  }
}

void CapoPlacer::plotBins(vector<CapoBin*>& layout) {
  double sum = 0.0;
  unsigned b;
  vector<const CapoBin*> tmpBins;
  tmpBins.reserve(layout.size() + _placedBins.size());
  for (b = 0; b != layout.size(); b++) {
    sum += layout[b]->getNumCells();
    tmpBins.push_back(layout[b]);
  }

  for (b = 0; b != _placedBins.size(); b++) {
    sum += _placedBins[b]->getNumCells();
    tmpBins.push_back(_placedBins[b]);
  }

  if ((layout.size() + _placedBins.size()) > 200 &&
      sum / (layout.size() + _placedBins.size()) < 10)
    return;

  if (_mainLayerNum > 15) return;
  ostringstream fileNameStream;
  fileNameStream << "capoBins-layer" << _mainLayerNum << "."
                 << _feedbackLayerNum << ".gpl";
  const string fileName = fileNameStream.str();
  ofstream out(fileName.c_str());
  out << "# " << layout.size() << " Capo bins from layer " << _mainLayerNum
      << "." << _feedbackLayerNum << endl
      << "# Need gnuplot 3.7 or later version \n"
      << "#set noxtics \n#set noytics \nset noborder \nset nokey\n"
      << "set title 'Capo Layer " << _mainLayerNum << "." << _feedbackLayerNum
      << "'\n";
  for (b = 0; b != tmpBins.size(); b++) {
    const CapoBin& curBin = *tmpBins[b];
    Point ctr = curBin.getBBox().getGeomCenter();
    out << "set label '"
        << int(curBin.getRelativeWhitespace() * 100)
        // << "%/" << int(nodesToMove(curBin)*100.0/curBin.getCellIds().size())
        << "%' at " << ctr.x << ", " << ctr.y << "\n";
    Point weberLoc = getWeberLocation(curBin);
    out << "set arrow from " << ctr.x << ", " << ctr.y << " to " << weberLoc.x
        << "," << weberLoc.y << endl;
  }
  out << "# set terminal postscript\n# set output '" << "capoBins-layer"
      << _mainLayerNum << "." << _feedbackLayerNum << ".eps'" << endl;
  out << " plot '-' w l\n";
  for (b = 0; b != tmpBins.size(); b++) {
    xyDrawRectangle(out, tmpBins[b]->getBBox());
    out << endl;
  }
  out << "EOF\n";
  out << "pause -1 'Press any key' " << endl;

  cout << " Saved " << fileName << endl;
}

BBox CapoPlacer::getPinOffset(unsigned cellId, unsigned netId) const {
  if (_params.useActualPinLocs && _rbplace.isFixed(cellId)) {
    const HGraphWDimensions& hgraph = getNetlistHGraph();
    const HGFNode& node = hgraph.getNodeByIdx(cellId);
    itHGFEdgeLocal e;
    for (e = node.edgesBegin(); e != node.edgesEnd(); ++e) {
      if ((*e)->getIndex() == netId) break;
    }

    abkassert(e != node.edgesEnd(),
              "asked for pin offset for an edge not on the node");

    unsigned edgeOffset = e - node.edgesBegin();

    return hgraph.edgesOnNodePinOffsetBBox(edgeOffset, cellId,
                                           _placement.getOrient(cellId));
  } else {
    return BBox(0., 0., 0., 0.);
  }
}

BBox CapoPlacer::getPinLocation(unsigned cellId, unsigned netId) const {
  BBox pinOffset = getPinOffset(cellId, netId);

  pinOffset.TranslateBy(_placement[cellId]);

  return pinOffset;
}

unsigned CapoPlacer::getTotalNetCut() const {
  const HGraphWDimensions& hgraph = getNetlistHGraph();
  vector<int> nodeBin;
  getBinMembership(nodeBin);
  unsigned numTotalCut = 0;
  itHGFEdgeGlobal edge;
  for (edge = hgraph.edgesBegin(); edge != hgraph.edgesEnd(); edge++) {
    int binID = -2;
    itHGFNodeLocal node;
    for (node = (*edge)->nodesBegin(); node != (*edge)->nodesEnd(); node++) {
      unsigned nodeId = (*node)->getIndex();

      // uncomment the following line if we don't want to count
      // nets that are connected to pads
      // if(nodeBin[nodeId] == -1) continue;

      if (binID == -2) {
        binID = nodeBin[nodeId];
        continue;
      }
      if (binID != nodeBin[nodeId]) {
        numTotalCut++;
        break;
      }
    }
  }
  return numTotalCut;
}

void CapoPlacer::getNetCutInfo(vector<unsigned>& externalCutNets,
                               vector<unsigned>& containedNets) const {
  vector<int> nodeBin;
  getBinMembership(nodeBin);
  unsigned numBins = _placedBins.size() + _curLayout->size();
  externalCutNets.clear();
  externalCutNets.reserve(numBins);
  containedNets.clear();
  containedNets.reserve(numBins);
  for (unsigned i = 0; i < numBins; i++) {
    externalCutNets.push_back(0);
    containedNets.push_back(0);
  }

  itHGFEdgeGlobal edge;
  BitBoard binUsed(numBins);
  const HGraphWDimensions& hgraph = getNetlistHGraph();
  for (edge = hgraph.edgesBegin(); edge != hgraph.edgesEnd(); edge++) {
    int binID = -2;
    itHGFNodeLocal node;
    binUsed.clear();
    for (node = (*edge)->nodesBegin(); node != (*edge)->nodesEnd(); node++) {
      unsigned nodeId = (*node)->getIndex();

      // uncomment the following line if we don't want to count
      // nets that are connected to pads
      // if(nodeBin[nodeId] == -1) continue;

      if (binID == -2) {
        binID = nodeBin[nodeId];
        //              if (binID>=0) binUsed.setBit(binID);
        continue;
      }
      if (binID != nodeBin[nodeId]) {
        if (binID >= 0 && !binUsed.isBitSet(binID)) {
          binUsed.setBit(binID);
          externalCutNets[binID]++;
        }
        if (nodeBin[nodeId] >= 0 && !binUsed.isBitSet(nodeBin[nodeId])) {
          binUsed.setBit(nodeBin[nodeId]);
          externalCutNets[nodeBin[nodeId]]++;
        }
      }
    }
    if (binUsed.getNumBitsSet() == 1 && binID >= 0) {
      const vector<unsigned>& usedBins = binUsed.getIndicesOfSetBits();
      containedNets[usedBins[0]]++;
    }
  }
}

void CapoPlacer::getBinMembership(vector<int>& nodeBin) const {
  const HGraphWDimensions& hgraph = getNetlistHGraph();
  nodeBin.assign(hgraph.getNumNodes(), -1);
  // Initialize the bin number to -1 so that cells that are not
  // in any bin will have -1 as its bin id

  // cells that are in placed bins
  int b, vectorSize = static_cast<int>(_placedBins.size());
  for (b = 0; b < vectorSize; b++) {
    const CapoBin& curBin = *(_placedBins[b]);
    vector<unsigned>::const_iterator nItr;
    for (nItr = curBin.cellIdsBegin(); nItr != curBin.cellIdsEnd(); nItr++)
      nodeBin[*nItr] = b;
  }

  // cells that are in unplaced bins
  // The bin id continues from the placed bins
  int baseID = static_cast<int>(_placedBins.size());
  vectorSize = static_cast<int>(_curLayout->size());
  for (b = 0; b < vectorSize; b++) {
    const CapoBin& curBin = *(*_curLayout)[b];
    vector<unsigned>::const_iterator nItr;
    for (nItr = curBin.cellIdsBegin(); nItr != curBin.cellIdsEnd(); nItr++)
      nodeBin[*nItr] = b + baseID;
  }
}

void CapoPlacer::getBinMembership(vector<int>& nodeBin,
                                  const vector<CapoBin*>& bins) const {
  const HGraphWDimensions& hgraph = getNetlistHGraph();
  nodeBin.assign(hgraph.getNumNodes(), -1);
  // Initialize the bin number to -1 so that cells that are not
  // in any bin will have -1 as its bin id

  // cells that are in bins
  int b, vectorSize = static_cast<int>(bins.size());
  for (b = 0; b < vectorSize; b++) {
    const CapoBin& curBin = *(bins[b]);
    vector<unsigned>::const_iterator nItr;
    for (nItr = curBin.cellIdsBegin(); nItr != curBin.cellIdsEnd(); nItr++)
      nodeBin[*nItr] = b;
  }
}

void CapoPlacer::getHierCellNames(vector<string>& cellNames) const {
  const HGraphWDimensions& hgraph = getNetlistHGraph();
  cellNames.clear();
  cellNames.insert(cellNames.end(), hgraph.getNumNodes(), string());

  // cells that are in placed bins
  unsigned vectorSize = _placedBins.size();
  ::bit_vector seenCell(hgraph.getNumNodes(), false);

  for (unsigned b = 0; b < vectorSize; b++) {
    const CapoBin& curBin = *(_placedBins[b]);
    abkfatal(curBin.getName() != string(),
             " NULL bin name in CapoPlacer::getHierCellName()"
             " --- make sure that small bins were not replaced earlier\n");

    // cout<<"Shortened Bin name: "<<shortBinName<<endl;

    vector<unsigned>::const_iterator nItr;
    for (nItr = curBin.cellIdsBegin(); nItr != curBin.cellIdsEnd(); nItr++) {
      abkfatal(!seenCell[*nItr], "already found this cell");
      seenCell[*nItr] = true;
      ostringstream nameSuffix;
      nameSuffix << "_" << *nItr;
      cellNames[*nItr] = curBin.getName() + nameSuffix.str();
    }
  }

  // cells that are in unplaced bins
  // The bin id continues from the placed bins
  for (unsigned k = 0; k != cellNames.size(); k++) {
    if (cellNames[k] == string()) cellNames[k] = hgraph.getNodeNameByIndex(k).c_str();
  }
}

void CapoPlacer::printBinMembership(ostream& os) const {
  os << "# every line below is of the form  \"<CellName> <BinNumber>\"\n"
     << "# cells that do not belong to any bin will have "
     << "BinNumber = -1\n";
  if (_placedBins.size() > 0)
    os << "# bins with BinNumber from 0 to " << _placedBins.size() - 1
       << " are placed bins\n";

  vector<int> nodeBin;
  getBinMembership(nodeBin);

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  itHGFNodeGlobal node;
  for (node = hgraph.nodesBegin(); node != hgraph.nodesEnd(); ++node)
    os << setw(50) << hgraph.getNodeNameByIndex((*node)->getIndex()) << "  "
       << setw(4) << nodeBin[(*node)->getIndex()] << endl;
}

void CapoPlacer::printBinMembership(const char* fileName) const {
  ofstream outFile(fileName);
  printBinMembership(outFile);
  outFile.close();
}

void CapoPlacer::printHierCellNames(ostream& os) const {
  os << "# every line below is of the form  \"<CellName> <NewCellName>\"\n"
     << "# cells that do not belong to any bin keep original names \n";
  //     << "# cells that do not belong to any bin will have "
  //     << " NoBin|<id> as name\n";

  vector<string> cellNames;
  getHierCellNames(cellNames);

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  itHGFNodeGlobal node;
  for (node = hgraph.nodesBegin(); node != hgraph.nodesEnd(); node++)
    os << setw(45) << hgraph.getNodeNameByIndex((*node)->getIndex()) << "  "
       << setw(25) << cellNames[(*node)->getIndex()] << endl;
}

void CapoPlacer::printHierCellNames(const char* fileName) const {
  ofstream outFile(fileName);
  printHierCellNames(outFile);
}

void CapoPlacer::saveBinsCopy() {
  int b, vectorSize = static_cast<int>(_placedBins.size());

  for (b = 0; b < vectorSize; b++) {
    CapoBin* newBin = new CapoBin(_placedBins[b]);
    _saveBins.push_back(newBin);
  }

  // cells that are in unplaced bins
  vectorSize = static_cast<int>(_curLayout->size());
  for (b = 0; b < vectorSize; b++) {
    CapoBin* newBin = new CapoBin((*_curLayout)[b]);
    _saveBins.push_back(newBin);
  }
}

void CapoPlacer::saveBinsAsFloorplan(const char* fileName,
                                     vector<CapoBin*> bins) {
  const string binsFileName = string(fileName) + ".bins";
  const string netsFileName = string(fileName) + ".nets";
  const string plFileName = string(fileName) + ".pl";

  ofstream binsFile(binsFileName.c_str());
  ofstream netsFile(netsFileName.c_str());
  ofstream plFile(plFileName.c_str());

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  // write the .bins file and .pl file
  binsFile << "UCSC bins 1.0" << endl;
  binsFile << TimeStamp() << User() << Platform() << endl;
  binsFile << "NumSoftRectangularBins : " << bins.size() << endl;
  binsFile << "NumHardRectilinearBins : 0" << endl;
  binsFile << "NumTerminals : " << hgraph.getNumTerminals() << endl << endl;

  plFile << "UCSC pl 1.0" << endl;
  plFile << TimeStamp() << User() << Platform() << endl;

  int b, vectorSize = static_cast<int>(bins.size());
  vector<string> binNames;

  for (b = 0; b < vectorSize; b++) {
    const CapoBin& curBin = *(bins[b]);
    const BBox& binBox = curBin.getBBox();
    binNames.push_back(curBin.getName());
    // double binArea = (binBox.xMax-binBox.xMin)*(binBox.yMax-binBox.yMin);
    double binArea = 0;
    vector<unsigned>::const_iterator nItr;
    for (nItr = curBin.cellIdsBegin(); nItr != curBin.cellIdsEnd(); nItr++)
      binArea += hgraph.getNodeWidth(*nItr) * hgraph.getNodeHeight(*nItr);

    double binAR = (binBox.xMax - binBox.xMin) / (binBox.yMax - binBox.yMin);
    double binWidth = sqrt(binArea * binAR);
    double binHeight = binWidth / binAR;

    binsFile << curBin.getName() << " softrectangular " << binArea
             << " 0.5 2.0 " << endl;

    /*
      binsFile<<binName<<" hardrectilinear 4 (0, 0) ("
    <<"0, "<<binBox.yMax-binBox.yMin<<") ("
    <<binBox.xMax-binBox.xMin<<", "<<binBox.yMax-binBox.yMin<<") ("
    <<binBox.xMax-binBox.xMin<<", 0)"<<endl;
      plFile<<setw(8) <<binName<<" "<<setw(10)<<binBox.xMin
    <<setw(10)<<binBox.yMin<<endl;
    */

    plFile << setw(8) << curBin.getName() << " " << setw(10) << binBox.xMin
           << setw(10) << binBox.yMin << " DIMS = (" << binWidth << ", "
           << binHeight << ")" << endl;
  }

  itHGFNodeGlobal nIt;
  for (nIt = hgraph.terminalsBegin(); nIt != hgraph.terminalsEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned termIdx = node.getIndex();
    if (!hgraph.getNodeNameByIndex(termIdx).empty()) {
      binsFile << hgraph.getNodeNameByIndex(termIdx) << " ";
      plFile << setw(8) << hgraph.getNodeNameByIndex(termIdx) << " ";
    } else {
      binsFile << "p" << termIdx + 1 << " ";
      plFile << "p" << termIdx + 1 << " ";
    }
    binsFile << "terminal " << endl;

    const Point termLoc = _rbplace[termIdx];
    plFile << setw(10) << termLoc.x << setw(10) << termLoc.y << endl;
  }

  // write the .nets file
  vector<int> netDegrees;
  vector<int> nodeBin;
  getBinMembership(nodeBin, bins);  // mapping from cells to bins

  unsigned numEdges = 0;
  unsigned numPins = 0;
  itHGFEdgeGlobal e;
  for (e = hgraph.edgesBegin(); e != hgraph.edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    bool needNet = 0;
    if (edge.getDegree() > 1) {
      itHGFNodeLocal firstNode = (*e)->nodesBegin();
      if (hgraph.isTerminal((*firstNode)->getIndex()))
        needNet = 1;
      else {
        unsigned firstBinId = nodeBin[(*firstNode)->getIndex()];
        itHGFNodeLocal v;
        for (v = (*e)->nodesBegin() + 1; v != (*e)->nodesEnd(); v++) {
          if (hgraph.isTerminal((*v)->getIndex())) {
            needNet = 1;
            break;
          }
          unsigned binId = nodeBin[(*v)->getIndex()];
          if (binId != firstBinId)  // atleast 1 different bin
          {
            needNet = 1;
            break;
          }
        }
      }
    }
    if (needNet) {
      ++numEdges;
      unsigned netDegree = 0;
      ::bit_vector seenBins(bins.size(), 0);
      itHGFNodeLocal v;
      for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); v++) {
        const HGFNode& node = (**v);
        unsigned idx = node.getIndex();
        if (hgraph.isTerminal(idx)) {
          ++netDegree;
          continue;
        }

        unsigned binId = nodeBin[idx];
        if (!seenBins[binId]) {
          seenBins[binId] = 1;
          ++netDegree;
        }
      }
      netDegrees.push_back(netDegree);
      numPins += netDegree;
    }
  }

  netsFile << TimeStamp() << User() << Platform() << endl;
  netsFile << "NumNets : " << numEdges << endl;
  netsFile << "NumPins : " << numPins << endl << endl;

  unsigned netIdx = 0;
  for (e = hgraph.edgesBegin(); e != hgraph.edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    bool needNet = 0;
    if (edge.getDegree() > 1) {
      itHGFNodeLocal firstNode = (*e)->nodesBegin();
      if (hgraph.isTerminal((*firstNode)->getIndex()))
        needNet = 1;
      else {
        unsigned firstBinId = nodeBin[(*firstNode)->getIndex()];
        itHGFNodeLocal v;
        for (v = (*e)->nodesBegin() + 1; v != (*e)->nodesEnd(); v++) {
          if (hgraph.isTerminal((*v)->getIndex())) {
            needNet = 1;
            break;
          }
          unsigned binId = nodeBin[(*v)->getIndex()];
          if (binId != firstBinId)  // atleast 1 different bin
          {
            needNet = 1;
            break;
          }
        }
      }
    }

    if (needNet) {
      netsFile << "NetDegree : " << netDegrees[netIdx] << " " << endl;
      ++netIdx;
      ::bit_vector seenBins(bins.size(), 0);
      itHGFNodeLocal v;
      unsigned nodeOffset = 0;
      unsigned edgeId = (**e).getIndex();
      for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); v++, nodeOffset++) {
        const HGFNode& node = (**v);
        unsigned idx = node.getIndex();
        if (hgraph.isTerminal(idx)) {
          netsFile << setw(10) << hgraph.getNodeNameByIndex(idx) << " B "
                   << endl;
          continue;
        }

        Point pOffset = hgraph.nodesOnEdgePinOffset(nodeOffset, edgeId);

        Point absOffset = _rbplace[idx];
        absOffset.x = absOffset.x + pOffset.x;
        absOffset.y = absOffset.y + pOffset.y;

        unsigned binId = nodeBin[idx];
        CapoBin* bin = bins[binId];
        Point binCenter = bin->getCenter();
        double binWidth = bin->getWidth();
        double binHeight = bin->getHeight();
        double xOff = 100 * (absOffset.x - binCenter.x) / binWidth;
        double yOff = 100 * (absOffset.y - binCenter.y) / binHeight;

        if (!seenBins[binId]) {
          seenBins[binId] = 1;
          netsFile << binNames[binId] << " B  : %" << xOff << " %" << yOff
                   << endl;
        }
      }
    }
  }

  binsFile.close();
  netsFile.close();
  plFile.close();
}

void CapoPlacer::solveQuadraticMinSoln() {
  if (_params.verb.getForActions())
    cout << endl << "SOR Solver Invoked ..." << endl;

  vector<CapoBin*> bins;
  int b, numBins = static_cast<int>(_curLayout->size());
  for (b = 0; b < numBins; b++) {
    CapoBin* newBin = new CapoBin((*_curLayout)[b]);
    bins.push_back(newBin);
  }

  vector<vector<unsigned> > nodeIds(numBins);
  vector<BBox> binBBox(numBins);

  for (b = 0; b < numBins; b++) {
    const CapoBin& curBin = *(bins[b]);
    const BBox& binBox = curBin.getBBox();
    binBBox[b] = binBox;
    vector<unsigned>::const_iterator nItr;
    for (nItr = curBin.cellIdsBegin(); nItr != curBin.cellIdsEnd(); nItr++)
      nodeIds[b].push_back(*nItr);
  }

  AnalyticalSolver solver(_rbplace, _placement, nodeIds, binBBox);
  if (_params.quadSkipNetsLargerThan != UINT_MAX)
    solver.setSkipNetsLargerThan(_params.quadSkipNetsLargerThan);

  BBox layoutBBox = _rbplace.getBBox(true);
  double maxDimension = std::max(layoutBBox.getHeight(), layoutBBox.getWidth());
  double epsilon = maxDimension / 100;
  if (!_params.noCOG) {
    if (_params.useLinearQP) {
      solver.solveLinearMin(epsilon);
    } else {
      solver.solveQuadraticMin(epsilon);
    }
  } else {
    if (_params.useLinearQP) {
      solver.solveLinearSOR(epsilon, 0);
    } else {
      solver.solveSOR(epsilon, 0);
    }
  }

  vector<vector<Point> >& nodeLocs = solver.getNodeLocs();
  for (b = 0; b < numBins; ++b) {
    for (unsigned i = 0; i < nodeIds[b].size(); ++i) {
      _placement[nodeIds[b][i]] = nodeLocs[b][i];
    }

    delete bins[b];
  }
}

unsigned CapoPlacer::getNumRegionedCellsInBin(CapoBin* bin) {
  unsigned numConstrCells = 0;
  if (constraints.getSize() > 0) {
    const vector<unsigned>& nodeIds = bin->getCellIds();
    for (unsigned i = 0; i != nodeIds.size(); i++) {
      unsigned cellId = nodeIds[i];
      unsigned groupId = cellToGrpMapping[cellId];
      if (groupId != UINT_MAX) numConstrCells++;
    }
  }
  return (numConstrCells);
}

#ifdef USEFLOW
bool CapoPlacer::floorplanBinWithFlows(CapoBin& bin) {
  vector<unsigned> movableCells;
  vector<BBox> cellBBoxes;

  for (unsigned i = 0; i < bin.getCellIds().size(); ++i) {
    movableCells.push_back(bin.getCellIds()[i]);
    cellBBoxes.push_back(bin.getBBox());
  }

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  if (_params.verb.getForActions()) {
    cout << "Floorplanning with flows" << endl;
  }
  _rbplace.doUnconstrainedXFlow(_placement, hgraph, movableCells, cellBBoxes,
                                _placement, _params.maxMem);
  _rbplace.doUnconstrainedYFlow(_placement, hgraph, movableCells, cellBBoxes,
                                _placement, _params.maxMem);
  if (_params.verb.getForActions()) {
    cout << endl;
  }

  return true;
}
#endif

bool CapoPlacer::floorplanBinWithSOR(CapoBin& bin) {
  vector<vector<unsigned> > nodeIds(1);
  vector<BBox> binBBox(1);

  for (unsigned b = 0; b < 1; b++) {
    const BBox& binBox = bin.getBBox();
    binBBox[b] = binBox;
    vector<unsigned>::const_iterator nItr;
    for (nItr = bin.cellIdsBegin(); nItr != bin.cellIdsEnd(); nItr++)
      nodeIds[b].push_back(*nItr);
  }

  AnalyticalSolver solver(_rbplace, _placement, nodeIds, binBBox);
  if (_params.quadSkipNetsLargerThan != UINT_MAX)
    solver.setSkipNetsLargerThan(_params.quadSkipNetsLargerThan);

  BBox layoutBBox = _rbplace.getBBox(true);
  double maxDimension = std::max(layoutBBox.getHeight(), layoutBBox.getWidth());
  double epsilon = maxDimension / 100;
  if (!_params.noCOG) {
    if (_params.useLinearQP) {
      solver.solveLinearMin(epsilon);
    } else {
      solver.solveQuadraticMin(epsilon);
    }
  } else {
    if (_params.useLinearQP) {
      solver.solveLinearSOR(epsilon, 0);
    } else {
      solver.solveSOR(epsilon, 0);
    }
  }

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  vector<vector<Point> >& nodeLocs = solver.getNodeLocs();
  for (unsigned i = 0; i < nodeIds[0].size(); ++i) {
    unsigned cellId = nodeIds[0][i];

    if (hgraph.isNodeMacro(cellId) && nodeNeedsFP(cellId)) {
      if (_params.wtCut || _params.doFastPlaceMoves) {
        Timer psTimer;
        updatePointSetsForMovedCell(cellId, _placement[cellId], nodeLocs[0][i]);
        psTimer.stop();
        CapoPlacer::capoTimer::PointSetTime += psTimer.getUserTime();
      }
      _placement[cellId] = nodeLocs[0][i];
    }
  }

  return true;
}

bool CapoPlacer::floorplanBin(CapoBin& bin, bool minWL, bool lookAheadMode) {
  const HGraphWDimensions& hgraph = getNetlistHGraph();
  if (_params.splitterParams.eco && bin._ecoFP) {
    ++CapoPlacer::EcoStats::ecoFPsUsed;
    for (vector<unsigned>::const_iterator cellIt = bin.cellIdsBegin();
         cellIt != bin.cellIdsEnd(); ++cellIt) {
      if (hgraph.isNodeMacro(*cellIt)) {
        if (_params.wtCut || _params.doFastPlaceMoves) {
          Timer psTimer;
          updatePointSetsForMovedCell(*cellIt, _placement[*cellIt],
                                      _oraclePlacement[*cellIt]);
          psTimer.stop();
          CapoPlacer::capoTimer::PointSetTime += psTimer.getUserTime();
        }
        _placement[*cellIt] = _oraclePlacement[*cellIt];
        _placement.setOrient(*cellIt, _oraclePlacement.getOrient(*cellIt));
      }
    }
    return true;
  }

  std::vector<unsigned> nodeIds;
  double nodesArea = 0, macroArea = 0;
  for (vector<unsigned>::const_iterator cellIt = bin.cellIdsBegin();
       cellIt != bin.cellIdsEnd(); ++cellIt) {
    if (hgraph.getNodeByIdx(*cellIt).getDegree() > 0 ||
        hgraph.isNodeMacro(*cellIt))  // don't include any disconnected cells
    {
      nodeIds.push_back(*cellIt);
      double nodeArea =
          hgraph.getNodeWidth(*cellIt) * hgraph.getNodeHeight(*cellIt);
      nodesArea += nodeArea;
      if (hgraph.isNodeMacro(*cellIt)) macroArea += nodeArea;
    }
  }
  // double macroPercent = 100.0*macroArea/nodesArea;

  const BBox& binBBox = bin.getBBox();
  double siteSpacing = _rbplace.coreRowsBegin()->getSpacing();
  double rowSpacing = (_rbplace.coreRowsBegin() + 1)->getYCoord() -
                      _rbplace.coreRowsBegin()->getYCoord();
#ifdef USETRAFFIC
  BBox trafficBox = bin.getBBox();
  // this will likely need to be rewritten to use an arbitrary placement,
  // and not an rbplacement
  if (_params.useTraffic) {
    trafficMain(_rbplace, trafficBox, nodeIds, _params.trafficOdds,
                _params.trafficThresh);
  }
#endif

  std::unique_ptr<ParquetDBFromRBP> fpDB(new ParquetDBFromRBP(
      _placement, hgraph, nodeIds, binBBox, _params.noParquetRotation,
      _params.snapMacrosX || _params.snapMacrosY,
      static_cast<float>(siteSpacing), static_cast<float>(rowSpacing)));

  double reqdAR = binBBox.getWidth() / binBBox.getHeight();
  double reqdArea = binBBox.getArea();
  double maxWS = 100 * (reqdArea - nodesArea) / nodesArea;

  double coreCellArea = nodesArea - macroArea;
  double newMaxWS = 30;

  double maxWSHier = (100 * reqdArea - 100 * macroArea - newMaxWS * macroArea) /
                         (coreCellArea * (100 + newMaxWS)) -
                     1.0;
  // double newMaxWS = 100*(reqdArea-(macroArea +
  // (1+maxWSHier)*coreCellArea))/(macroArea + (1+maxWSHier)*coreCellArea);

  maxWSHier *= 100;
  if (maxWSHier < -100) maxWSHier = -80;

  if (maxWS < 0) maxWS = 0;

  std::unique_ptr<parquetfp::Command_Line> fpParams(
      new parquetfp::Command_Line());

  fpParams->verb = _params.verb;

  fpParams->getSeed = false;
  // fpParams->seed = 0;
  // fpParams->seed = 1082073072;
  // fpParams->getSeed = true;
  // fpParams->setSeed();
  fpParams->seed = SeedGen;
  fpParams->reqdAR = reqdAR;
  fpParams->maxWS = maxWS;
  bool scaleTerms = false;
  fpParams->scaleTerms = scaleTerms;
  fpParams->softBlocks = true;
  fpParams->solveMulti = true;
  fpParams->solveTop = true;

  // <aaronnn> use whatever we have as an initial solution I have this
  // because I want to reuse lookaheadFP's solution.  but it shouldn't be
  // worse than random even if lookahead is not enabled
  // fpParams->takePl = true;

  // <aaronnn> put FP in lookahead mode?
  fpParams->lookAheadFP = lookAheadMode;

  if (!lookAheadMode) {
    fpParams->shrinkToSize = static_cast<float>(_params.parquetShrink);
  }

  // <aaronnn> cluster macros
  if (_params.clusterMacros)
    fpParams->dontClusterMacros = false;
  else
    fpParams->dontClusterMacros = true;

  if (lookAheadMode)
    fpParams->FPrep = "BTree";
  else
    fpParams->FPrep = _params.FPrep;

  if (maxWSHier < -10)
    fpParams->maxWSHier = maxWSHier;
  else
    fpParams->maxWSHier = -10;

  fpParams->wireWeight = 0.4f;
  fpParams->areaWeight = 0.35f;
  fpParams->minWL = minWL;
  fpParams->initCompact = true;
  fpParams->compact = false;
  // fpParams->initQP = true;
  // fpParams->iterations = 10;
  fpParams->maxIterHier = 5;
  if (fpDB->getNumNodes() > 2000 || fpDB->getNumNodes() < 60)
    fpParams->maxTopLevelNodes = 100;

  if (_params.verb.getForMajStats() > 0)
    cout << endl << "REQD AR " << reqdAR << " maxWS " << maxWS << endl;
  fpDB->markTallNodesAsMacros((*bin.rowsBegin())->getHeight());

  // <aaronnn> tell FP the nodes we want FP-ed
  fpDB->markNodesNeedingFP(_nodesNeedingFP);

  // <aaronnn> tell FP what obstacles it has to deal with
  if (_params.FPEvadeObstacles && bin.getNumObstacles() > 0) {
    const std::vector<unsigned> obstacleIds(bin.getObstacleCellIds().begin(),
                                            bin.getObstacleCellIds().end());
    fpDB->addObstacles(_placement, hgraph, obstacleIds,
                       _params.FPEvadeObstaclesThresh);
  }

  // fpDB->save("temp");
  // fpParams->printAnnealerParams();

  fpParams->packleft = true;
  fpParams->packbot = true;

  if (_params.annealToCutLines) {
    const string::size_type len = bin.getName().length();
    if (len > 3) {
      if (bin.getName()[len - 2] == 'H') {
        if (bin.getName()[len - 1] == '1') fpParams->packbot = false;
        if (len > 6 && bin.getName()[len - 5] == 'V') {
          if (bin.getName()[len - 4] == '0') fpParams->packleft = false;
        }
      } else if (bin.getName()[len - 2] == 'V') {
        if (bin.getName()[len - 1] == '0') fpParams->packleft = false;
        if (len > 6 && bin.getName()[len - 5] == 'H') {
          if (bin.getName()[len - 4] == '1') fpParams->packbot = false;
        }
      }
    }

    if (!fpParams->packleft && _params.verb.getForActions() > 0)
      cout << "Forcing Parquet to pack to the right" << endl;
    if (!fpParams->packbot && _params.verb.getForActions() > 0)
      cout << "Forcing Parquet to pack to the top" << endl;
  }

  if (fpDB->getNumNodes() > 100 ||
      1)  // <aaronnn> lol || 1.
          // WARNING: it better always enter here or obstacles won't be handled
  {
    parquetfp::SolveMulti solveMulti(fpDB.get(), fpParams.get(),
                                     _params.maxMem);
    solveMulti.go();
    CapoPlacer::capoTimer::FloorClusterTime += solveMulti.clusterTime;
    CapoPlacer::capoTimer::FloorAnnealTime += solveMulti.annealTime;
  } else {
    parquetfp::Annealer annealer(fpParams.get(), fpDB.get(), _params.maxMem);
    annealer.go();
    CapoPlacer::capoTimer::FloorAnnealTime += annealer.annealTime;
  }

  // fpDB->expandDesign(binBBox.getWidth(), binBBox.getHeight());
  parquetfp::Point offset;

  double xMaxWMacro = static_cast<double>(fpDB->getXMaxWMacroOnly());
  double yMaxWMacro = static_cast<double>(fpDB->getYMaxWMacroOnly());

  // the following code pushes the whole design down&left to satisfy
  // top-right corner constraints
  offset.x = offset.y = 0.f;
  if (xMaxWMacro > (binBBox.getWidth()))
    offset.x = static_cast<float>(binBBox.getWidth() - xMaxWMacro);
  if (yMaxWMacro > (binBBox.getHeight()))
    offset.y = static_cast<float>(binBBox.getHeight() - yMaxWMacro);
  if (!_params.FPEvadeObstacles || bin.getNumObstacles() == 0) {
    // shift design only if FP obstacle evasion is disabled, or there are no
    // obstacles
    fpDB->shiftDesign(offset);
  }

  // following code shifts the entire floorplan to bin's lower left corner
  offset.x = static_cast<float>(binBBox.xMin);
  offset.y = static_cast<float>(binBBox.yMin);
  fpDB->shiftDesign(offset);

  if (!_params.FPEvadeObstacles || bin.getNumObstacles() == 0) {
    // shiftOptimize only if FP obstacle evasion is disabled, or there are no
    // obstacles
    parquetfp::Point outlineBottomLeft(offset);
    parquetfp::Point outlineTopRight;
    outlineTopRight.x = static_cast<float>(binBBox.xMax);
    outlineTopRight.y = static_cast<float>(binBBox.yMax);
#ifdef USEFLUTE
    bool useSteiner = false;
    fpDB->shiftOptimizeDesign(outlineBottomLeft, outlineTopRight, scaleTerms,
                              useSteiner, _params.verb);
#else
    fpDB->shiftOptimizeDesign(outlineBottomLeft, outlineTopRight, scaleTerms,
                              _params.verb);
#endif
  }
  // fpDB->cornerOptimizeDesign();

  // fpDB->plot("out.gpl", 0, 0, 0, 0, 0, 0, 0, 0);
  // exit(0);
  bool success = true;
  const std::vector<Point>& nodeLocs = fpDB->getNodeLocs();
  const std::vector<Orient>& nodeOrients = fpDB->getNodeOrients();

  for (unsigned i = 0; i < nodeIds.size(); ++i) {
    // updates the locations and dimensions of macros only
    // if(hgraph.isNodeMacro(nodeIds[i]))

    // <aaronnn> modified: updates the locations and dimensions of macros
    // that we wanted floorplanned only
    if (hgraph.isNodeMacro(nodeIds[i]) && nodeNeedsFP(nodeIds[i])) {
      unsigned cellId = nodeIds[i];

      if (!lookAheadMode && (_params.wtCut || _params.doFastPlaceMoves)) {
        Timer psTimer;
        updatePointSetsForMovedCell(cellId, _placement[cellId], nodeLocs[i]);
        psTimer.stop();
        CapoPlacer::capoTimer::PointSetTime += psTimer.getUserTime();
      }

      _placement[cellId] = nodeLocs[i];

      Orient savedOrient = _placement.getOrient(cellId);

      _placement.setOrient(cellId, nodeOrients[i]);

      BBox macroOutline;
      double nodeWidth, nodeHeight;

      unsigned angle = nodeOrients[i].getAngle();
      if (angle == 0 || angle == 180) {
        if (hgraph.isNodeSoft(cellId)) {
          if (!lookAheadMode) {
            // update the width and height if the node is soft
            const_cast<HGraphWDimensions&>(hgraph).setNodeWidth(
                fpDB->getNodes()->getNodeWidth(i), cellId);

            const_cast<HGraphWDimensions&>(hgraph).setNodeHeight(
                fpDB->getNodes()->getNodeHeight(i), cellId);
          }
        }
        nodeWidth = hgraph.getNodeWidth(cellId);
        nodeHeight = hgraph.getNodeHeight(cellId);
      } else {
        if (hgraph.isNodeSoft(cellId)) {
          if (!lookAheadMode) {
            // update the width and height, but need to flip it
            const_cast<HGraphWDimensions&>(hgraph).setNodeWidth(
                fpDB->getNodes()->getNodeHeight(i), cellId);

            const_cast<HGraphWDimensions&>(hgraph).setNodeHeight(
                fpDB->getNodes()->getNodeWidth(i), cellId);
          }
        }
        nodeWidth = hgraph.getNodeHeight(cellId);
        nodeHeight = hgraph.getNodeWidth(cellId);
      }

      if (hgraph.isNodeSoft(cellId)) {
        const_cast<HGraphWDimensions&>(hgraph).makeNodePinsCenter(cellId);

        //             printf("soft-block[%d]: width: %.2f height: %.2f AR: %.2f
        //             minAR: %.2f maxAR: %.2f\n",
        //                    cellId,
        //                    hgraph.getNodeWidth(cellId),
        //                    hgraph.getNodeHeight(cellId),
        //                    hgraph.getNodeWidth(cellId) /
        //                    hgraph.getNodeHeight(cellId),
        //                    hgraph.getNodeMinAR(cellId),
        //                    hgraph.getNodeMaxAR(cellId));
      }

      macroOutline.add(nodeLocs[i].x, nodeLocs[i].y);
      macroOutline.add(nodeLocs[i].x + nodeWidth, nodeLocs[i].y + nodeHeight);
      if (!bin._needFP) {
        if (greaterThanDouble(macroOutline.xMax, binBBox.xMax) ||
            greaterThanDouble(macroOutline.yMax, binBBox.yMax) ||
            lessThanDouble(macroOutline.xMin, binBBox.xMin) ||
            lessThanDouble(macroOutline.yMin, binBBox.yMin)) {
          if (lessThanFloat(fpParams->shrinkToSize, 0.f)) success = false;
          if (_params.verb.getForActions() > 5) {
            cout << "Macro out of specified region" << endl;
            cout << " Macro dims " << macroOutline << endl;
            cout << " Region dims " << binBBox << endl;
          }
        }
      } else  // for second time we have more relaxed constraints for
              // determining success
      {
        if (greaterThanDouble(macroOutline.xMax, binBBox.xMax) ||
            greaterThanDouble(macroOutline.yMax, binBBox.yMax) ||
            lessThanDouble(macroOutline.xMin, binBBox.xMin) ||
            lessThanDouble(macroOutline.yMin, binBBox.yMin)) {
          if (lessThanFloat(fpParams->shrinkToSize, 0.f)) success = false;
          if (_params.verb.getForActions() > 5) {
            cout << "Macro out of specified region" << endl;
            cout << " Macro dims " << macroOutline << endl;
            cout << " Region dims " << binBBox << endl;
          }
        }
      }

      // <aaronnn> check for overlaps too, when determining FP success
      if (lessThanFloat(fpParams->shrinkToSize, 0.f) &&
          _params.FPEvadeObstacles && bin.getNumObstacles() > 0) {
        const vector<unsigned>& obstacleCellIds = bin.getObstacleCellIds();
        for (unsigned i = 0; i < obstacleCellIds.size(); ++i) {
          unsigned obstacleCellId = obstacleCellIds[i];
          double obstacleWidth, obstacleHeight;
          unsigned obstacleAngle =
              _placement.getOrient(obstacleCellId).getAngle();
          bool obstacleNotRotated = obstacleAngle == 0 || obstacleAngle == 180;
          if (obstacleNotRotated) {
            obstacleWidth = hgraph.getNodeWidth(obstacleCellId);
            obstacleHeight = hgraph.getNodeHeight(obstacleCellId);
          } else {
            obstacleWidth = hgraph.getNodeHeight(obstacleCellId);
            obstacleHeight = hgraph.getNodeWidth(obstacleCellId);
          }
          BBox obstacleOutline;
          obstacleOutline.add(_placement[obstacleCellId].x,
                              _placement[obstacleCellId].y);
          obstacleOutline.add(_placement[obstacleCellId].x + obstacleWidth,
                              _placement[obstacleCellId].y + obstacleHeight);

          if (macroOutline.intersects(obstacleOutline)) {
            success = false;
            break;
          }
        }
      }

      if (lookAheadMode) {
        _placement.setOrient(cellId, savedOrient);
      }
    } else if (!lookAheadMode) {
      Point oldPlace(_placement[nodeIds[i]]);

      _placement[nodeIds[i]] = nodeLocs[i];

      double nodeWidth, nodeHeight;
      unsigned angle = _placement.getOrient(nodeIds[i]).getAngle();
      bool notRotated = angle == 0 || angle == 180;
      if (notRotated) {
        nodeWidth = hgraph.getNodeWidth(nodeIds[i]);
        nodeHeight = hgraph.getNodeHeight(nodeIds[i]);
      } else {
        nodeWidth = hgraph.getNodeHeight(nodeIds[i]);
        nodeHeight = hgraph.getNodeWidth(nodeIds[i]);
      }

      if (greaterThanDouble(_placement[nodeIds[i]].x + nodeWidth,
                            binBBox.xMax)) {
        _placement[nodeIds[i]].x = binBBox.xMax - nodeWidth;
      }
      if (lessThanDouble(_placement[nodeIds[i]].x, binBBox.xMin)) {
        _placement[nodeIds[i]].x = binBBox.xMin;
      }

      if (greaterThanDouble(_placement[nodeIds[i]].y + nodeHeight,
                            binBBox.yMax)) {
        _placement[nodeIds[i]].y = binBBox.yMax - nodeHeight;
      }
      if (lessThanDouble(_placement[nodeIds[i]].y, binBBox.yMin)) {
        _placement[nodeIds[i]].y = binBBox.yMin;
      }

      if (_params.wtCut || _params.doFastPlaceMoves) {
        Timer psTimer;
        updatePointSetsForMovedCell(nodeIds[i], oldPlace,
                                    _placement[nodeIds[i]]);
        psTimer.stop();
        CapoPlacer::capoTimer::PointSetTime += psTimer.getUserTime();
      }
    }
  }

  bool floorplanningSuccess = success;

  if (lookAheadMode) {
    return (success);
  }

  // if already gone up the hierarchy once, floorplan with only area and accept
  //    bool finalSuccess=true;
  if (!success && bin._needFP) success = true;

  // the following code disables going up the hierarchy if the layer number
  // is close to the top. Unsuccessfull floorplan is also accepted.
  if (_mainLayerNum < 1) {
    if (!success && _params.verb.getForActions() > 0)
      cout << "Not going up the hierarchy and accepting illegal floorplanning "
              "solution in early layers"
           << endl;
    success = true;
  }
  // uncomment following line, to disable going up the hierarchy all together
  // success = true;

  if (_params.plotFPOutlines) {
    ofstream dat;

    if (FloorplanStats::numInstances == 0)
      dat.open(_params.FPOutlinesFileName.c_str());
    else
      dat.open(_params.FPOutlinesFileName.c_str(), std::ios::app);

    dat << "\n# Floorplanning instance number "
        << FloorplanStats::numInstances + 1 << endl;
    dat << "# Bin name : " << bin.getName() << endl;
    dat << "# Bin id : " << bin.getIndex() << endl;
    dat << "# Packed to the " << (fpParams->packleft ? "left" : "right")
        << endl;
    dat << "# Packed to the " << (fpParams->packbot ? "bottom" : "top") << endl;
    dat << "# Floorplanning " << (floorplanningSuccess ? "Success!" : "Failure")
        << endl;
    dat << "# Final placement " << (success ? "accepted" : "rejected") << endl;
    dat << endl;
    BBox binBBox = bin.getBBox();
    double offset = (binBBox.xMax - binBBox.xMin < binBBox.yMax - binBBox.yMin)
                        ? binBBox.xMax - binBBox.xMin
                        : binBBox.yMax - binBBox.yMin;
    offset *= 0.15;
    BBox layout = _rbplace.getBBox();
    double cap = (layout.xMax - layout.xMin < layout.yMax - layout.yMin)
                     ? layout.xMax - layout.xMin
                     : layout.yMax - layout.yMin;
    cap *= 0.02;
    if (offset > cap) offset = cap;
    dat << binBBox.xMin << "  " << binBBox.yMin + offset << endl
        << binBBox.xMin << "  " << binBBox.yMax - offset << endl
        << binBBox.xMin + offset << "  " << binBBox.yMax << endl
        << binBBox.xMax - offset << "  " << binBBox.yMax << endl
        << binBBox.xMax << "  " << binBBox.yMax - offset << endl
        << binBBox.xMax << "  " << binBBox.yMin + offset << endl
        << binBBox.xMax - offset << "  " << binBBox.yMin << endl
        << binBBox.xMin + offset << "  " << binBBox.yMin << endl
        << binBBox.xMin << "  " << binBBox.yMin + offset << endl;
    dat.close();
  }

  return (success);
}

#ifdef USEPATOMA
bool CapoPlacer::patomaBin(CapoBin& bin, bool minWL) {
  const HGraphWDimensions& hgraph = getNetlistHGraph();
  std::vector<unsigned> nodeIds;
  double nodesArea = 0, macroArea = 0;
  vector<unsigned>::const_iterator cellIt;
  for (cellIt = bin.cellIdsBegin(); cellIt != bin.cellIdsEnd(); ++cellIt) {
    nodeIds.push_back(*cellIt);
    double nodeArea =
        hgraph.getNodeWidth(*cellIt) * hgraph.getNodeHeight(*cellIt);
    nodesArea += nodeArea;
    if (hgraph.isNodeMacro(*cellIt)) macroArea += nodeArea;
  }
  // double macroPercent = 100.0*macroArea/nodesArea;

  const BBox binBBox = bin.getBBox();
  const std::unique_ptr<ParquetDBFromRBP> fpDB(new ParquetDBFromRBP(
      _placement, hgraph, nodeIds, const_cast<BBox&>(binBBox),
      _params.noParquetRotation));

  const double reqdAR = binBBox.getWidth() / binBBox.getHeight();
  const double reqdArea = binBBox.getArea();
  const double maxWS = max(100 * (reqdArea - nodesArea) / nodesArea, 0.0);
  //    if(maxWS < 0)
  //       maxWS = 0;

  const double coreCellArea = nodesArea - macroArea;
  const double newMaxWS = 30;

  double maxWSHier =
      ((100 * reqdArea - 100 * macroArea - newMaxWS * macroArea) /
           (coreCellArea * (100 + newMaxWS)) -
       1.0);
  // double newMaxWS = 100*(reqdArea-(macroArea +
  // (1+maxWSHier)*coreCellArea))/(macroArea + (1+maxWSHier)*coreCellArea);
  maxWSHier *= 100;
  if (maxWSHier < -100) maxWSHier = -80;

  const std::unique_ptr<parquetfp::Command_Line> fpParams(
      new parquetfp::Command_Line());

  fpParams->verb = _params.verb;

  fpParams->getSeed = false;
  fpParams->seed = SeedGen;
  // fpParams->seed = 0;
  // fpParams->seed = 1082073072;
  // fpParams->getSeed = true;
  // fpParams->setSeed();
  fpParams->reqdAR = reqdAR;
  bool scaleTerms = false;
  fpParams->scaleTerms = scaleTerms;
  fpParams->maxWS = maxWS;
  fpParams->softBlocks = true;
  fpParams->solveMulti = true;
  fpParams->solveTop = true;
  // fpParams->dontClusterMacros = true;

  fpParams->FPrep = _params.FPrep;

  if (maxWSHier < -10)
    fpParams->maxWSHier = maxWSHier;
  else
    fpParams->maxWSHier = -10;

  //    fpParams->wireWeight = 0.4;
  //    fpParams->areaWeight = 0.35;
  //    if(minWL)
  //       fpParams->minWL = true;
  //    else
  //       fpParams->minWL = false;
  fpParams->initCompact = true;

  // area-only floorplanning
  fpParams->wireWeight = 0.0;
  fpParams->minWL = false;

  fpParams->compact = true;
  // fpParams->initQP = true;
  // fpParams->iterations = 10;
  fpParams->maxIterHier = 5;
  if (fpDB->getNumNodes() > 10000 || fpDB->getNumNodes() < 60)
    fpParams->maxTopLevelNodes = 100;

  if (_params.verb.getForMajStats() > 0)
    cout << endl << "REQD AR " << reqdAR << " maxWS " << maxWS << endl;
  fpDB->markTallNodesAsMacros((*bin.rowsBegin())->getHeight());

  // ----count # soft nodes and # macros-----
  unsigned numSoftNodes = 0;
  unsigned numMacros = 0;
  for (unsigned i = 0; i < fpDB->getNumNodes(); i++) {
    if (fpDB->isNodeSoft(i)) numSoftNodes++;

    if (fpDB->isMacro(i)) numMacros++;
  }

  // must set it to an end case, otherwise, infinite loop
  const bool terminalCase = (numMacros == 1);
  const bool allMacroSoft = numSoftNodes == numMacros;
  cout << "# Nodes: " << fpDB->getNumNodes()
       << " # Soft Nodes: " << numSoftNodes << " # Macros: " << numMacros;
  if (bin.isPatomaEndCase()) cout << " (PATOMA end case)";
  cout << endl;
  // -----------------------------------------

  // fpDB->save("temp");
  // fpParams->printAnnealerParams();

  if (allMacroSoft) {
    parquetfp::SolveMulti solveMulti(fpDB.get(), fpParams.get());
    std::unique_ptr<parquetfp::DB> clusteredDB(solveMulti.clusterOnly());

    vector<unsigned> allNodes(clusteredDB->getNumNodes());
    for (unsigned i = 0; i < clusteredDB->getNumNodes(); i++) allNodes[i] = i;

    Timer T;
    Patoma patoma(bin.getBBox(), allNodes, *clusteredDB);
    patoma.go(Patoma::LOOK_AHEAD);
    T.stop();

    // indirectly update *fpDB, may consider putting this fch
    // to fpDB:: instead
    solveMulti.updatePlaceUnCluster(clusteredDB.get());

    CapoPlacer::capoTimer::FloorClusterTime += solveMulti.clusterTime;
    CapoPlacer::capoTimer::FloorAnnealTime += T.getUserTime();
  } else if (fpDB->getNumNodes() > 100 || true) {
    parquetfp::SolveMulti solveMulti(fpDB.get(), fpParams.get());
    solveMulti.go();
    CapoPlacer::capoTimer::FloorClusterTime += solveMulti.clusterTime;
    CapoPlacer::capoTimer::FloorAnnealTime += solveMulti.annealTime;
  } else {
    parquetfp::Annealer annealer(fpParams.get(), fpDB.get());
    annealer.go();
    CapoPlacer::capoTimer::FloorAnnealTime += annealer.annealTime;
  }

  // fpDB->expandDesign(binBBox.getWidth(), binBBox.getHeight());
  parquetfp::Point offset;

  const double xMaxWMacro = fpDB->getXMaxWMacroOnly();
  const double yMaxWMacro = fpDB->getYMaxWMacroOnly();

  // the following code pushes the whole design down&left to satisfy
  // top-right corner constraints
  offset.x = offset.y = 0;
  if (xMaxWMacro > (binBBox.getWidth()))
    offset.x = binBBox.getWidth() - xMaxWMacro;
  if (yMaxWMacro > (binBBox.getHeight()))
    offset.y = binBBox.getHeight() - yMaxWMacro;
  fpDB->shiftDesign(offset);

  // following code shifts the entire floorplan to bin's lower left corner
  offset.x = binBBox.xMin;
  offset.y = binBBox.yMin;
  fpDB->shiftDesign(offset);

  parquetfp::Point outlineBottomLeft(offset);
  parquetfp::Point outlineTopRight;
  outlineTopRight.x = binBBox.xMax;
  outlineTopRight.y = binBBox.yMax;
  fpDB->shiftOptimizeDesign(outlineBottomLeft, outlineTopRight, scaleTerms,
                            _params.verb);

  // fpDB->plot("out.gpl", 0, 0, 0, 0, 0, 0, 0, 0);
  // exit(0);
  bool success = true;
  const std::vector<Point>& nodeLocs = fpDB->getNodeLocs();
  const std::vector<Orient>& nodeOrients = fpDB->getNodeOrients();

  for (unsigned i = 0; i < nodeIds.size(); ++i) {
    // updates the locations and dimensions of macros only
    // from *fpDB to _placement
    if (hgraph.isNodeMacro(nodeIds[i])) {
      unsigned cellId = nodeIds[i];

      _placement[cellId] = nodeLocs[i];
      _placement.setOrient(cellId, nodeOrients[i]);

      BBox obstacle;
      double nodeWidth, nodeHeight;

      unsigned angle = nodeOrients[i].getAngle();
      if (angle == 0 || angle == 180) {
        if (hgraph.isNodeSoft(cellId)) {
          // update the width and height if the node is soft
          const_cast<HGraphWDimensions&>(hgraph).setNodeWidth(
              fpDB->getNodes()->getNodeWidth(i), cellId);

          const_cast<HGraphWDimensions&>(hgraph).setNodeHeight(
              fpDB->getNodes()->getNodeHeight(i), cellId);
        }
        nodeWidth = hgraph.getNodeWidth(cellId);
        nodeHeight = hgraph.getNodeHeight(cellId);
      } else {
        if (hgraph.isNodeSoft(cellId)) {
          // update the width and height, but need to flip it
          const_cast<HGraphWDimensions&>(hgraph).setNodeWidth(
              fpDB->getNodes()->getNodeHeight(i), cellId);

          const_cast<HGraphWDimensions&>(hgraph).setNodeHeight(
              fpDB->getNodes()->getNodeWidth(i), cellId);
        }
        nodeWidth = hgraph.getNodeHeight(cellId);
        nodeHeight = hgraph.getNodeWidth(cellId);
      }

      if (hgraph.isNodeSoft(cellId)) {
        const_cast<HGraphWDimensions&>(hgraph).makeNodePinsCenter(cellId);

        //             printf("soft-block[%d]: width: %.2f height: %.2f AR: %.2f
        //             minAR: %.2f maxAR: %.2f\n",
        //                    cellId,
        //                    hgraph.getNodeWidth(cellId),
        //                    hgraph.getNodeHeight(cellId),
        //                    hgraph.getNodeWidth(cellId) /
        //                    hgraph.getNodeHeight(cellId),
        //                    hgraph.getNodeMinAR(cellId),
        //                    hgraph.getNodeMaxAR(cellId));
      }

      obstacle.add(nodeLocs[i].x, nodeLocs[i].y);
      obstacle.add(nodeLocs[i].x + nodeWidth, nodeLocs[i].y + nodeHeight);
      if (!bin._needFP) {
        if (obstacle.xMax > 1.0 * binBBox.xMax ||
            obstacle.yMax > 1.0 * binBBox.yMax ||
            obstacle.xMin < 1.0 * binBBox.xMin ||
            obstacle.yMin < 1.0 * binBBox.yMin) {
          success = false;
          if (_params.verb.getForActions() > 5) {
            cout << "BIN out of specified region " << endl;
            cout << " BIN dims " << obstacle << endl;
            cout << " RGN dims " << binBBox << endl;
          }
        }
      } else  // for second time we have more relaxed constraints for
              // determining success
      {
        if (obstacle.xMax > 1.0 * binBBox.xMax ||
            obstacle.yMax > 1.0 * binBBox.yMax ||
            obstacle.xMin < 1.0 * binBBox.xMin ||
            obstacle.yMin < 1.0 * binBBox.yMin) {
          success = false;
          if (_params.verb.getForActions() > 5) {
            cout << "BIN out of specified region " << endl;
            cout << " BIN dims " << obstacle << endl;
            cout << " RGN dims " << binBBox << endl;
          }
        }
      }
    }  // done with the case for macros
    else {
      // if the current node is not a macro (i.e. cell),
      // still use its actual location though it may
      // be outside of "binBBox"

      //       _placement[nodeIds[i]] = binBBox.getGeomCenter();
      _placement[nodeIds[i]] = nodeLocs[i];
    }
  }  // done updating every nodes (macros || cells)

  // if allready gone up the hierarchy once, floorplan with only area and accept
  //    bool finalSuccess=true;
  //     if(!success && bin._needFP)

  // if already gone up one, i.e. it is an end-case for Patoma
  // then accept the potentially illegal solution
  if (bin.isPatomaEndCase() && !terminalCase) {
    cout << "accept the bad solution." << endl;
    success = true;
  }

  // DONE if it's a terminal case, do not split anymore
  if (success && terminalCase) bin.setPatomaEndCase(true);

  // the following code disables going up the hierarchy if the layer number
  // is close to the top. Unsuccessfull floorplan is also accepted.
  if (_mainLayerNum < 2) {
    if (!success && _params.verb.getForActions() > 0)
      cout << "Not going up the hierarchy and accepting illegal floorplanning "
              "solution in early layers"
           << endl;
    success = true;
  }
  // uncomment following line, to disable going up the hierarchy all together
  // success = true;

  cout << "PATOMA: " << ((success) ? string("success") : string("failed"))
       << endl;
  return (success);
}
#endif

void CapoPlacer::processMacrosAfterlayer(vector<CapoBin*>& layout,
                                         CapoPlacer::LayerType layerType) {
  const HGraphWDimensions& hgraph = getNetlistHGraph();

  bool atLeast1FPed = false;
  for (unsigned i = 0; i < layout.size(); ++i) {
    if (layout[i] != 0) {
      CapoBin* bin = layout[i];
      bool shouldCommit = false;
      if (layerType == FP_LAYER) shouldCommit = bin->_floorplanned;
#ifdef USEPATOMA
      else if (layerType == PATOMA_LAYER)
        shouldCommit = bin->_floorplanned || bin->isPatomaMacrosCommitted();
#endif
      else
        abkwarn(false, "invalid layer type specified");

      if (shouldCommit) {
        BBox binBBox = bin->getBBox();
        bin->_floorplanned = false;
        vector<unsigned>::const_iterator cellIt;
        for (cellIt = bin->cellIdsBegin(); cellIt != bin->cellIdsEnd();
             ++cellIt) {
          unsigned cellId = *cellIt;

          // <aaronnn> don't block out every site with macros,
          // just the sites under macros that we wanted floorplanned.
          if (hgraph.isNodeMacro(cellId) && nodeNeedsFP(cellId)) {
            atLeast1FPed = true;
            _placed[cellId] = true;
            _cellToBinMap[cellId] = NULL;

            if (_params.wtCut || _params.doFastPlaceMoves) {
              Timer psTimer;
              updatePointSetsForNewlyPlacedCell(cellId, _placement[cellId],
                                                _placement[cellId]);
              psTimer.stop();
              CapoPlacer::capoTimer::PointSetTime += psTimer.getUserTime();
            }

            _rbplace.splitCoreRowsByCell(_placement, cellId);
          }
        }
        // <aaronnn> we don't want to remove all macros,
        // just the macros that we wanted floorplanned
        bin->removeFloorplannedMacros();

        if (_params.splitterParams.WSOracle) {
          bin->clearOracleCellIds();
        }
        if (_params.splitterParams.eco && !bin->_ecoFP) {
          bin->setEcoAllowed(false);
        }
#ifdef USEPATOMA
        bin->setPatomaEndCase(false);
        bin->setPatomaMacrosCommitted(false);
#endif
      }
    }
  }

  // remove any overlaps between macros
  if (_params.removeMacroOverlap && atLeast1FPed) {
    Timer macroOverlapSetup;
    vector<unsigned> cellIds;
    vector<unsigned> potentiallyMovedCells;
    vector<Point> potentiallyMovedLocs;
    cellIds.reserve(hgraph.getNumNodes());
    potentiallyMovedCells.reserve(hgraph.getNumNodes());
    potentiallyMovedLocs.reserve(hgraph.getNumNodes());
    for (unsigned i = 0; i < hgraph.getNumNodes(); ++i) {
      if ((_placed[i] && hgraph.isNodeMacro(i)) || _rbplace.isFixed(i) ||
          hgraph.isTerminal(i)) {
        cellIds.push_back(i);
        if (!_rbplace.isFixed(i) && !hgraph.isTerminal(i)) {
          potentiallyMovedCells.push_back(i);
          potentiallyMovedLocs.push_back(_placement[i]);
        }
      }
    }
    macroOverlapSetup.stop();
    CapoPlacer::capoTimer::MacroOverlapTime += macroOverlapSetup.getUserTime();
    Timer macroOverlapRem;
    bool snapMacrosX = _params.snapMacrosX;
    bool snapMacrosY = _params.snapMacrosY;
    bool putCellsInRows = false;
    bool print = true;
    unsigned badMoveCap = 1; /* maybe a commandline option? */
    _rbplace.removeMacroOverlaps(_params.maxMem, _placement, cellIds,
                                 snapMacrosX, snapMacrosY, putCellsInRows,
                                 print, badMoveCap);
    macroOverlapRem.stop();
    CapoPlacer::capoTimer::MacroOverlapTime += macroOverlapRem.getUserTime();

    // do flows if everything is legal
#ifdef USEFLOW
    if (_params.flowAfterParquet) {
      double totalMacroOverlap = _rbplace.calcOverlapSpecific(
          _placement, _rbplace.getBBox(true), cellIds, OnlyMacro, false);

      if (equalDouble(totalMacroOverlap, 0.)) {
        double hpwlBefore = estimateWL();
        BBox core = _rbplace.getBBox(true);
        bool flowX = true;
        double Xdist = max(_rbplace.getRowSpacing(),
                           _params.flowLimitRatio * core.getWidth());
        bool flowY = true;
        double Ydist = max(_rbplace.getRowHeight(),
                           _params.flowLimitRatio * core.getHeight());

        cout << "HPWL before flows " << hpwlBefore << endl;
        _rbplace.doFlowHelper(NULL, flowX, Xdist, flowY, Ydist, _placement,
                              &cellIds, _params.maxMem);
        double hpwlAfter = estimateWL();
        cout << "HPWL after flows " << hpwlAfter << endl;
      } else {
        cout << "Skipping flow due to overlap." << endl;
      }
    }
#endif

    // update the point sets
    if (_params.wtCut || _params.doFastPlaceMoves) {
      Timer psTimer;
      for (unsigned i = 0; i < potentiallyMovedCells.size(); ++i) {
        if (potentiallyMovedLocs[i] != _placement[potentiallyMovedCells[i]]) {
          updatePointSetsForMovedCell(potentiallyMovedCells[i],
                                      potentiallyMovedLocs[i],
                                      _placement[potentiallyMovedCells[i]]);
        }
      }
      psTimer.stop();
      CapoPlacer::capoTimer::PointSetTime += psTimer.getUserTime();
    }
  }
}

void CapoPlacer::getBinsWMacrosInFront(vector<CapoBin*>& layout) {
  vector<CapoBin*> newLayout;
  newLayout.reserve(layout.size());

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  for (unsigned i = 0; i < layout.size(); ++i) {
    if (layout[i] != 0) {
      CapoBin* bin = layout[i];
      bool foundMacro = false;
      vector<unsigned>::const_iterator cellIt;
      for (cellIt = bin->cellIdsBegin(); cellIt != bin->cellIdsEnd();
           ++cellIt) {
        unsigned cellId = *cellIt;
        if (hgraph.isNodeMacro(cellId)) {
          foundMacro = true;
          break;
        }
      }
      if (foundMacro)
        newLayout.insert(newLayout.begin(), layout[i]);
      else
        newLayout.push_back(layout[i]);
    } else
      newLayout.push_back(layout[i]);
  }
  for (unsigned i = 0; i < layout.size(); ++i) layout[i] = newLayout[i];
}

void CapoPlacer::getFPedBinsIntoNewLayout(vector<CapoBin*>& layout,
                                          vector<CapoBin*>& newLayout) {
  // foreach non-NULL bin in "layout", add it to "newLayout"
  // only move floorplanned or hierarchical bins to newLayout
  for (unsigned i = 0; i < layout.size(); ++i)
    if (layout[i] != NULL) {
      if (layout[i]->_floorplanned || layout[i]->_needFP) {
        newLayout.push_back(layout[i]);
        layout[i] = NULL;
      }
    }
}

#ifdef USEPATOMA
void CapoPlacer::getPatomaedBinsIntoNewLayout(vector<CapoBin*>& layout,
                                              vector<CapoBin*>& newLayout) {
  for (unsigned i = 0; i < layout.size(); ++i)
    if (layout[i] != NULL) {
      // only move Patoma end-cases to the new layout
      // and further partition the successful cases
      if (layout[i]->isPatomaEndCase()) {
        newLayout.push_back(layout[i]);
        layout[i] = NULL;
      }
    }
}
#endif

void CapoPlacer::swapCells(CapoBin* b1, unsigned b1Cell, CapoBin* b2,
                           unsigned b2Cell) {
  vector<unsigned> cellIds;

  cellIds = b1->getCellIds();
  if (b1Cell != UINT_MAX) {
    vector<unsigned>::iterator pos =
        std::find(cellIds.begin(), cellIds.end(), b1Cell);
    abkfatal(pos != cellIds.end(), "cell has to be in the source to remove it");
    cellIds.erase(pos);
  }
  if (b2Cell != UINT_MAX) {
    cellIds.push_back(b2Cell);
    _placed[b2Cell] = false;
  }
  b1->resetCellIds(cellIds);
  updateInfoAboutBin(b1);

  cellIds = b2->getCellIds();
  if (b2Cell != UINT_MAX) {
    vector<unsigned>::iterator pos =
        std::find(cellIds.begin(), cellIds.end(), b2Cell);
    abkfatal(pos != cellIds.end(), "cell has to be in the source to remove it");
    cellIds.erase(pos);
  }
  if (b1Cell != UINT_MAX) {
    cellIds.push_back(b1Cell);
    _placed[b1Cell] = true;
  }
  b2->resetCellIds(cellIds);
  updateInfoAboutBin(b2);
}

double calcDistance(const CapoBin& b1, const CapoBin& b2) {
  Point b1center = b1.getCenter(), b2center = b2.getCenter();
  return fabs(b1center.x - b2center.x) + fabs(b1center.y - b2center.y);
}

CapoBin* CapoPlacer::findNearestBin(vector<CapoBin*> haveSpace,
                                    const CapoBin& start, double spaceneeded,
                                    unsigned* swapCell, double* bestdistance) {
  double spaceavailable, distance, bestCellSize;
  unsigned bestCell;
  CapoBin* closest = NULL;
  *swapCell = UINT_MAX;
  *bestdistance = DBL_MAX;

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  for (vector<CapoBin*>::iterator b = haveSpace.begin(); b != haveSpace.end();
       ++b) {
    if (*b == 0) continue;
    spaceavailable = (**b).getWhitespace();
    if (spaceavailable < 1.e-6) {
      *b = 0;
      continue;
    }
    distance = calcDistance(start, **b);
    if (*bestdistance < distance) continue;
    bestCell = UINT_MAX;
    bestCellSize = DBL_MAX;
    if (spaceavailable < spaceneeded) {
      double target = spaceneeded - spaceavailable;
      if (target >= spaceneeded) continue;
      for (unsigned i = 0; i < (**b).getNumCells(); ++i) {
        unsigned tmp = (**b).getCellIds()[i];
        if (!_rbplace.isCoreCell(tmp) || _rbplace.isFixed(tmp)) continue;
        double size = hgraph.getWeight(tmp);
        if (size >= target && size < spaceneeded && size < bestCellSize) {
          bestCell = tmp;
          bestCellSize = size;
        }
      }
      if (bestCell == UINT_MAX) continue;
    }
    closest = *b;
    *swapCell = bestCell;
    *bestdistance = distance;
  }

  return closest;
}

// Begin Upper and Lower bound functions
double CapoPlacer::getWorstCaseHPWL(void) {
  double totalWL = 0.0;
  unsigned numEdges = _netUpperBounds.size();
  for (unsigned i = 0; i < numEdges; ++i) {
    totalWL += _netUpperBounds[i].getHalfPerim();
  }
  return totalWL;
}

double CapoPlacer::getBestCaseHPWL(void) {
  double totalWL = 0.0;
  unsigned numEdges = _netLowerBounds.size();
  for (unsigned i = 0; i < numEdges; ++i) {
    if (_netLowerBounds[i].xMin != -DBL_MAX &&
        _netLowerBounds[i].yMin != -DBL_MAX &&
        _netLowerBounds[i].xMax != DBL_MAX &&
        _netLowerBounds[i].yMax != DBL_MAX)
      totalWL += _netLowerBounds[i].getHalfPerim();
  }
  return totalWL;
}

void CapoPlacer::updateNetUpperAndLowerBounds() {
  unsigned numEdges = getNetlistHGraph().getNumEdges();
  for (unsigned i = 0; i < numEdges; ++i) {
    _netUpperBounds[i] = evalNetUpperBoundBB(i);
    _netLowerBounds[i] = evalNetLowerBoundBB(i);
  }
}

BBox CapoPlacer::evalNetLowerBoundBB(unsigned netToEval) const {
  // the BBox rval will be used to store intervals which are independent
  // for both y and x lower bounds.  We may have an x lower bound and not
  // y and vice versa.  Prior to establishing a lower bound interval
  // we will keep optimistic top bottom and left right estimates.
  BBox rval(-DBL_MAX, -DBL_MAX, DBL_MAX, DBL_MAX);
  const HGraphWDimensions& hgraph = getNetlistHGraph();
  const HGFEdge& e = hgraph.getEdgeByIdx(netToEval);
  itHGFNodeLocal nIt = e.nodesBegin();

  double optL = -DBL_MAX;  // greatest left of all bboxes.
  double optR = DBL_MAX;   // least right of all bboxes.
  double optT = DBL_MAX;   // least top of all bboxes.
  double optB = -DBL_MAX;  // greatest bottom of all bboxes.

  // look for the first cell with a bin or placement
  for (; nIt != e.nodesEnd(); ++nIt) {
    HGFNode* n = *nIt;
    unsigned nIdx = n->getIndex();
    abkfatal(nIdx < _cellToBinMap.size(), "Node Index out of range");
    const CapoBin* bin = _cellToBinMap[nIdx];
    if (bin) {
      BBox bin_bbox = bin->getBBox();
      if (bin_bbox.xMin > optL) optL = bin_bbox.xMin;
      if (bin_bbox.xMax < optR) optR = bin_bbox.xMax;
      if (bin_bbox.yMin > optB) optB = bin_bbox.yMin;
      if (bin_bbox.yMax < optT) optT = bin_bbox.yMax;
    } else if (_placed[nIdx]) {
      Point plNpt = _placement[nIdx];
      if (plNpt.x < optR) optR = plNpt.x;
      if (plNpt.x > optL) optL = plNpt.x;
      if (plNpt.y < optT) optT = plNpt.y;
      if (plNpt.y > optB) optB = plNpt.y;
    } else {
      abkfatal(0, "Where is this cell?");
    }
  }

  if (optL > optR) {
    rval.xMin = optR;
    rval.xMax = optL;
  }
  if (optB > optT) {
    rval.yMin = optT;
    rval.yMax = optB;
  }
  return rval;
}

BBox CapoPlacer::evalNetUpperBoundBB(unsigned netToEval) const {
  const HGraphWDimensions& hgraph = getNetlistHGraph();
  BBox rval(_rbplace.getBBox());
  const HGFEdge& e = hgraph.getEdgeByIdx(netToEval);
  itHGFNodeLocal nIt = e.nodesBegin();
  // initialize rval to the first the bin of the first cell
  // that has a bin or, if the cell is placed, take its location
  while (nIt != e.nodesEnd()) {
    unsigned nIdx1 = (*nIt)->getIndex();
    abkassert(nIdx1 < _cellToBinMap.size(), "Node Index out of range");
    const CapoBin* bin1 = _cellToBinMap[nIdx1];
    if (_placed[nIdx1]) {
      rval.setTopRight(_placement[nIdx1]);
      rval.setBottomLeft(_placement[nIdx1]);
      break;
    } else if (bin1) {
      rval = bin1->getBBox();
      // cout<<"In Upper Bound, found cell "<<nIdx1<<" in a bin"<<endl;
      // cout<<"   Bin has BB " << rval << endl;
      break;
    } else {
      abkassert(0, "Where is this cell?");
    }
  }
  for (; nIt != e.nodesEnd(); ++nIt) {
    unsigned nIdx = (*nIt)->getIndex();
    abkassert(nIdx < _cellToBinMap.size(), "Node Index out of range");
    const CapoBin* thisbin = _cellToBinMap[nIdx];
    if (_placed[nIdx]) {
      rval += _placement[nIdx];

    } else if (thisbin) {
      BBox cellBBox = thisbin->getBBox();
      rval.expandToInclude(cellBBox);
    } else {
      abkassert(0, "Where is this cell?");
    }
  }
  return rval;
}
// End upper and lower bound functions

void CapoPlacer::recordCutLine(const Point& a, const Point& b, unsigned bin1,
                               unsigned bin2) {
  if (bin2 < bin1) {
    unsigned tmp = bin1;
    bin1 = bin2;
    bin2 = tmp;
  }

  ofstream dat(_params.CutLinesFileName.c_str(), std::ios::app);

  dat << "# cut line between bins " << bin1 << " and " << bin2 << endl;
  dat << a << endl << b << endl << endl;

  dat.close();
}

void CapoPlacer::recordCutLayer(unsigned layer) {
  ofstream dat;
  if (layer == 0) {
    dat.open(_params.CutLinesFileName.c_str());
  } else {
    dat.open(_params.CutLinesFileName.c_str(), std::ios::app);
  }

  dat << "# Writing cut lines for layer number " << layer << endl << endl;

  dat.close();
}

bool CapoPlacer::nodeNeedsFP(unsigned nodeID) const {
  abkassert(nodeID >= 0 && nodeID < _nodesNeedingFP.size(),
            "nodeNeedsFP OOB");  // XXX remove later
  return _nodesNeedingFP[nodeID];
}

void CapoPlacer::markNodeNeedsFP(unsigned nodeID) {
  abkassert(nodeID >= 0 && nodeID < _nodesNeedingFP.size(),
            "markNodeNeedsFP OOB");  // XXX remove later
  _nodesNeedingFP[nodeID] = true;
}

void CapoPlacer::unmarkNodeNeedsFP(unsigned nodeID) {
  abkassert(nodeID >= 0 && nodeID < _nodesNeedingFP.size(),
            "markNodeNeedsFP OOB");  // XXX remove later
  _nodesNeedingFP[nodeID] = false;
}

void CapoPlacer::markBinNeedsFP(unsigned binIdx) {
  if ((*_curLayout)[binIdx] == NULL) return;

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  for (unsigned i = 0; i < (*_curLayout)[binIdx]->getNumCells(); ++i) {
    unsigned cellId = (*_curLayout)[binIdx]->getCellIds()[i];

    if (!hgraph.isNodeMacro(cellId)) {
      // do nothing for non-macros
      continue;
    }

    _nodesNeedingFP[cellId] = true;
  }
}

void CapoPlacer::markNodesNeedFP(const vector<unsigned>& cellIds) {
  const HGraphWDimensions& hgraph = getNetlistHGraph();

  for (unsigned i = 0; i < cellIds.size(); ++i) {
    if (!hgraph.isNodeMacro(cellIds[i])) {
      // do nothing for non-macros
      continue;
    }

    _nodesNeedingFP[cellIds[i]] = true;
  }
}

bool CapoPlacer::getDefaultSplitDir(CapoBin& bin) {
  bool splitHoriz = false;

  double stretchFactor =
      (bin.getAvgRowSpacing() > 1.5 * bin.getAvgCellHeight() ? 2 : 1);

  splitHoriz = stretchFactor * bin.getHeight() > bin.getWidth();

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  double vertFactor = _params.vertFactor;  // 7.0;
  if (!splitHoriz && _params.mixedSize) {
    const vector<unsigned>& cellIds = bin.getCellIds();
    for (unsigned i = 0; i < cellIds.size(); i++) {
      unsigned idx = cellIds[i];
      if (hgraph.isNodeMacro(idx)) {
        vertFactor = 0;
        break;
      }
    }
  }

  unsigned numLeaves = bin.getNumCells();
  unsigned numRows = bin.getNumRows();

  if (!splitHoriz)
    splitHoriz = (vertFactor * bin.getHeight() > bin.getWidth() &&
                  numLeaves <= _params.smallPlaceThreshold * numRows &&
                  bin.getRelativeWhitespace() <= 0.3 && numRows <= 5);

  return splitHoriz;
}

bool CapoPlacer::atBottomLayers() {
  unsigned binswithfewrows = 0;
  for (unsigned k = 0; k < _curLayout->size(); k++) {
    CapoBin* bin = (*_curLayout)[k];

    if (bin == NULL) continue;

    if (bin->getNumRows() <= 3) binswithfewrows++;
  }

  bool atBot = ((double)binswithfewrows / (double)_curLayout->size()) >= 0.8;
  return atBot;
}

CapoPlacer::~CapoPlacer() {
  vector<CapoBin*> allBins;
  allBins.reserve(_placedBins.size() + _layout[0].size() + _layout[1].size());
  allBins.insert(allBins.end(), _placedBins.begin(), _placedBins.end());
  _placedBins.clear();
  allBins.insert(allBins.end(), _layout[0].begin(), _layout[0].end());
  _layout[0].clear();
  allBins.insert(allBins.end(), _layout[1].begin(), _layout[1].end());
  _layout[1].clear();
  std::sort(allBins.begin(), allBins.end());
  vector<CapoBin*>::iterator new_end =
      std::unique(allBins.begin(), allBins.end());
  allBins.erase(new_end, allBins.end());
  for (unsigned i = 0; i < allBins.size(); ++i) {
    delete allBins[i];
  }
  allBins.clear();
  if (_params.useFastSt || _params.fastPlaceFastSt) {
    fastSteiner::stnr1_package_done();
    if (_fastStPt != NULL) {
      delete[] _fastStPt;
      _fastStPt = NULL;
    }
    if (_fastStParent != NULL) {
      delete[] _fastStParent;
      _fastStParent = NULL;
    }
  }
}

bool CapoPlacer::doFastPlaceMoves() {
  Timer fastPlaceMovement;

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  double beginSteinerLen = 0.;
  if (!_params.fastPlaceHPWL) {
    for (unsigned i = 0; i < hgraph.getNumEdges(); ++i) {
      const bool usePins = true;
      double netlen = std::numeric_limits<double>::max();

      if (_params.fastPlaceFastSt)
        netlen = std::min(netlen, pinPlEval::oneNetHPWL_FastSteiner(
                                      _placement, hgraph, i, usePins));

#ifdef USEFLUTE
      if (_params.fastPlaceFLUTE)
        netlen = std::min(netlen, pinPlEval::oneNetHPWL_Flute(
                                      _placement, hgraph, i, usePins));
#endif

      beginSteinerLen += netlen;
    }
  }

  double beginHPWL = estimateWL();
  cout << "Beginning HPWL: " << beginHPWL << endl;
  if (!_params.fastPlaceHPWL)
    cout << "Beginning Estimated Steiner Length: " << beginSteinerLen << endl;

  cout << endl << "Starting Fast Place Moves" << endl;
  cout << "    [" << flush;

  BitBoard needsUpdate(_rbplace.getNumCells());
  vector<FastPlaceMove> moveList(_rbplace.getNumCells());
  std::map<unsigned, bool> affectedBins;
  vector<double> curNetLens(hgraph.getNumEdges(), -1.);

  for (unsigned i = 0; i < _curLayout->size(); ++i) {
    CapoBin* bin = (*_curLayout)[i];
    if (bin == NULL) continue;

    for (unsigned j = 0; j < bin->getCellIds().size(); ++j) {
      unsigned cellId = bin->getCellIds()[j];

      moveList[cellId].cellId = cellId;
      moveList[cellId].source = bin;
      chooseFastPlaceMove(moveList[cellId], curNetLens);
    }
  }

  for (unsigned i = 0; i < _placedBins.size(); ++i) {
    CapoBin* bin = _placedBins[i];
    if (bin == NULL) continue;

    for (unsigned j = 0; j < bin->getCellIds().size(); ++j) {
      unsigned cellId = bin->getCellIds()[j];

      moveList[cellId].cellId = cellId;
      moveList[cellId].source = bin;
      chooseFastPlaceMove(moveList[cellId], curNetLens);
    }
  }

  vector<unsigned> moveOrder(_rbplace.getNumCells());
  for (unsigned i = 0; i < moveOrder.size(); ++i) {
    moveOrder[i] = i;
  }

  sort(moveOrder.begin(), moveOrder.end(), DecreasingByCost(moveList));

  double sumOfCosts = 0.;

  while (lessThanDouble(moveList[moveOrder.back()].cost, 0.)) {
    makeFastPlaceMove(moveList[moveOrder.back()]);

    sumOfCosts += moveList[moveOrder.back()].cost;

    ++fastPlaceMovesMade;
    cout << "+" << flush;

    needsUpdate.clear();

    addCellsNeedingUpdate(moveList[moveOrder.back()].source, needsUpdate);
    addCellsNeedingUpdate(moveList[moveOrder.back()].dest, needsUpdate);
    addCellsNeedingUpdate(moveList[moveOrder.back()].cellId, needsUpdate);
    invalidateNetLens(moveList[moveOrder.back()].cellId, curNetLens);
    affectedBins[moveList[moveOrder.back()].source->getIndex()] = true;
    affectedBins[moveList[moveOrder.back()].dest->getIndex()] = true;

    const vector<unsigned>& updatedCells = needsUpdate.getIndicesOfSetBits();
    for (unsigned i = 0; i < updatedCells.size(); ++i) {
      chooseFastPlaceMove(moveList[updatedCells[i]], curNetLens);
    }

    sort(moveOrder.begin(), moveOrder.end(), DecreasingByCost(moveList));
  }

  fastPlaceMovesMade = 1;

  cout << "]" << endl;

  // update the bin lists
  vector<CapoBin*> tempBins;
  tempBins.reserve(_placedBins.size());
  for (unsigned i = 0; i < _placedBins.size(); ++i) {
    if (_placedBins[i] == NULL) continue;

    if (affectedBins.find(_placedBins[i]->getIndex()) == affectedBins.end()) {
      tempBins.push_back(_placedBins[i]);
    } else {
      _curLayout->push_back(_placedBins[i]);
      for (unsigned j = 0; j < _placedBins[i]->getCellIds().size(); ++j) {
        _placed[_placedBins[i]->getCellIds()[j]] = false;
      }
    }
  }
  _placedBins = tempBins;
  _layout[0] = *_curLayout;
  _layout[1] = *_curLayout;

  // don't update the congestion map,
  // it will be rebuilt at the next layer
  bool congWSsave = _params.splitterParams.congestionShifting;
  _params.splitterParams.congestionShifting = false;
  for (unsigned i = 0; i < _curLayout->size(); ++i) {
    updateInfoAboutBin((*_curLayout)[i]);
  }
  _params.splitterParams.congestionShifting = congWSsave;

  bool haveUnplacedBins = (_curLayout->size() > 0);

  buildPointSetsFromScratch();

  double endSteinerLen = 0.;
  if (!_params.fastPlaceHPWL) {
    for (unsigned i = 0; i < hgraph.getNumEdges(); ++i) {
      const bool usePins = true;
      double netlen = std::numeric_limits<double>::max();

      if (_params.fastPlaceFastSt)
        netlen = std::min(netlen, pinPlEval::oneNetHPWL_FastSteiner(
                                      _placement, hgraph, i, usePins));

#ifdef USEFLUTE
      if (_params.fastPlaceFLUTE)
        netlen = std::min(netlen, pinPlEval::oneNetHPWL_Flute(
                                      _placement, hgraph, i, usePins));
#endif

      endSteinerLen += netlen;
    }
  }

  double endHPWL = estimateWL();

  fastPlaceMovement.stop();

  CapoPlacer::capoTimer::FastPlaceMoveTime += fastPlaceMovement.getUserTime();

  cout << "Fast Place Movement took " << fastPlaceMovement << endl << endl;
  cout << "Ending HPWL: " << endHPWL << endl;
  if (!_params.fastPlaceHPWL)
    cout << "Ending Estimated Steiner Length: " << endSteinerLen << endl
         << endl;

  return haveUnplacedBins;
}

double CapoPlacer::calcFastPlaceCost(unsigned cellId, const Point& newCellPos,
                                     vector<double>& curNetLens) {
  double cost = 0.;

  const HGraphWDimensions& hgraph = getNetlistHGraph();
  const HGFNode& node = hgraph.getNodeByIdx(cellId);

  for (itHGFEdgeLocal edge = node.edgesBegin(); edge != node.edgesEnd();
       ++edge) {
    unsigned eId = (*edge)->getIndex();
    Point currPos = getPinLocation(cellId, eId).getGeomCenter();
    Point newPos = currPos - _placement[cellId] + newCellPos;

    vector<Point> points;

    for (unsigned i = 0; i < _pointSet_fixed[eId].size(); ++i) {
      if (_placed[cellId] && _pointSet_fixed[eId][i].pt == currPos) {
        if (_pointSet_fixed[eId][i].count > 1) {
          points.push_back(_pointSet_fixed[eId][i].pt);
        }
      } else {
        points.push_back(_pointSet_fixed[eId][i].pt);
      }
    }

    for (unsigned i = 0; i < _pointSet_movable[eId].size(); ++i) {
      if (!_placed[cellId] && _pointSet_movable[eId][i].pt == currPos) {
        if (_pointSet_movable[eId][i].count > 1) {
          points.push_back(_pointSet_movable[eId][i].pt);
        }
      } else {
        points.push_back(_pointSet_movable[eId][i].pt);
      }
    }

    sort(points.begin(), points.end());
    vector<Point>::iterator new_end = unique(points.begin(), points.end());
    points.erase(new_end, points.end());

    bool useSteiner = !_params.fastPlaceHPWL && points.size() < 20;

    double beforeLen, afterLen;
    if (useSteiner) {
      beforeLen = curNetLens[eId];
      if (beforeLen < 0) {
        if (binary_search(points.begin(), points.end(), currPos)) {
          beforeLen = calculateSteinerFromPointSet(
              points,
#ifdef USEFLUTE
              _params.fastPlaceFastSt, _params.fastPlaceFLUTE, false);
#else
              _params.fastPlaceFastSt, false);
#endif
        } else {
          vector<Point> temp = points;
          temp.push_back(currPos);
          sort(temp.begin(), temp.end());
          beforeLen = calculateSteinerFromPointSet(
              temp,
#ifdef USEFLUTE
              _params.fastPlaceFastSt, _params.fastPlaceFLUTE, false);
#else
              _params.fastPlaceFastSt, false);
#endif
        }
        curNetLens[eId] = beforeLen;
      }

      if (binary_search(points.begin(), points.end(), newPos)) {
        afterLen = calculateSteinerFromPointSet(points,
#ifdef USEFLUTE
                                                _params.fastPlaceFastSt,
                                                _params.fastPlaceFLUTE, false);
#else
                                                _params.fastPlaceFastSt, false);
#endif
      } else {
        points.push_back(newPos);
        sort(points.begin(), points.end());
        afterLen = calculateSteinerFromPointSet(points,
#ifdef USEFLUTE
                                                _params.fastPlaceFastSt,
                                                _params.fastPlaceFLUTE, false);
#else
                                                _params.fastPlaceFastSt, false);
#endif
      }
    } else {
      bool containsBefore = false;
      bool containsAfter = false;
      BBox net;
      for (unsigned i = 0; i < points.size(); ++i) {
        net += points[i];
        containsBefore = containsBefore || net.contains(currPos);
        containsAfter = containsAfter || net.contains(newPos);
        if (containsBefore && containsAfter) break;
      }

      if (containsBefore) {
        beforeLen = net.getHalfPerim();
      } else {
        BBox temp = net;
        temp += currPos;
        beforeLen = temp.getHalfPerim();
      }

      if (containsAfter) {
        afterLen = net.getHalfPerim();
      } else {
        net += newPos;
        afterLen = net.getHalfPerim();
      }
    }

    cost += (afterLen - beforeLen);
  }

  return cost;
}

void CapoPlacer::chooseFastPlaceMove(FastPlaceMove& move,
                                     vector<double>& curNetLens) {
  move.cost = 0.;
  move.pos = _placement[move.cellId];
  move.dest = move.source;

  const HGraphWDimensions& hgraph = getNetlistHGraph();
  double cellArea = hgraph.getWeight(move.cellId);

  Point currentOffsetRatio =
      _placement[move.cellId] - move.source->getBBox().getBottomLeft();
  currentOffsetRatio.x /= move.source->getBBox().getWidth();
  currentOffsetRatio.y /= move.source->getBBox().getHeight();

  // loop over bins neighbors above
  for (unsigned i = 0; i < move.source->_neighborsAbove.size(); ++i) {
    CapoBin* neighbor = move.source->_neighborsAbove[i];
    if (neighbor == NULL) continue;

    if (lessThanDouble(neighbor->getWhitespace(), cellArea)) continue;

    // calculate a new position in this new bin
    Point newOffset = currentOffsetRatio;
    newOffset.x *= neighbor->getBBox().getWidth();
    newOffset.y *= neighbor->getBBox().getHeight();

    Point newPos = neighbor->getBBox().getBottomLeft() + newOffset;

    double Cost = calcFastPlaceCost(move.cellId, newPos, curNetLens);

    if (lessThanDouble(Cost, move.cost)) {
      move.cost = Cost;
      move.dest = neighbor;
      move.pos = newPos;
    }
  }

  // loop over bins neighbors below
  for (unsigned i = 0; i < move.source->_neighborsBelow.size(); ++i) {
    CapoBin* neighbor = move.source->_neighborsBelow[i];
    if (neighbor == NULL) continue;

    if (lessThanDouble(neighbor->getWhitespace(), cellArea)) continue;

    // calculate a new position in this new bin
    Point newOffset = currentOffsetRatio;
    newOffset.x *= neighbor->getBBox().getWidth();
    newOffset.y *= neighbor->getBBox().getHeight();

    Point newPos = neighbor->getBBox().getBottomLeft() + newOffset;

    double Cost = calcFastPlaceCost(move.cellId, newPos, curNetLens);

    if (lessThanDouble(Cost, move.cost)) {
      move.cost = Cost;
      move.dest = neighbor;
      move.pos = newPos;
    }
  }

  // loop over bins left neighbors
  for (unsigned i = 0; i < move.source->_leftNeighbors.size(); ++i) {
    CapoBin* neighbor = move.source->_leftNeighbors[i];
    if (neighbor == NULL) continue;

    if (lessThanDouble(neighbor->getWhitespace(), cellArea)) continue;

    // calculate a new position in this new bin
    Point newOffset = currentOffsetRatio;
    newOffset.x *= neighbor->getBBox().getWidth();
    newOffset.y *= neighbor->getBBox().getHeight();

    Point newPos = neighbor->getBBox().getBottomLeft() + newOffset;

    double Cost = calcFastPlaceCost(move.cellId, newPos, curNetLens);

    if (lessThanDouble(Cost, move.cost)) {
      move.cost = Cost;
      move.dest = neighbor;
      move.pos = newPos;
    }
  }

  // loop over bins right neighbors
  for (unsigned i = 0; i < move.source->_rightNeighbors.size(); ++i) {
    CapoBin* neighbor = move.source->_rightNeighbors[i];
    if (neighbor == NULL) continue;

    if (lessThanDouble(neighbor->getWhitespace(), cellArea)) continue;

    // calculate a new position in this new bin
    Point newOffset = currentOffsetRatio;
    newOffset.x *= neighbor->getBBox().getWidth();
    newOffset.y *= neighbor->getBBox().getHeight();

    Point newPos = neighbor->getBBox().getBottomLeft() + newOffset;

    double Cost = calcFastPlaceCost(move.cellId, newPos, curNetLens);

    if (lessThanDouble(Cost, move.cost)) {
      move.cost = Cost;
      move.dest = neighbor;
      move.pos = newPos;
    }
  }
}

void CapoPlacer::makeFastPlaceMove(FastPlaceMove& move) {
  vector<unsigned> cellsTemp = move.source->getCellIds();

  cellsTemp.erase(find(cellsTemp.begin(), cellsTemp.end(), move.cellId));

  move.source->resetCellIds(cellsTemp);

  CapoBin* temp = move.source;
  move.source = move.dest;
  move.dest = temp;

  cellsTemp = move.source->getCellIds();
  cellsTemp.push_back(move.cellId);

  move.source->resetCellIds(cellsTemp);

  updatePointSetsForMovedCell(move.cellId, _placement[move.cellId], move.pos);

  _placement[move.cellId] = move.pos;
}

void CapoPlacer::invalidateNetLens(unsigned cellId,
                                   vector<double>& curNetLens) {
  // invalidate (set to -1) the length all the nets attached
  // to this cell
  const HGraphWDimensions& hgraph = getNetlistHGraph();
  const HGFNode& node = hgraph.getNodeByIdx(cellId);

  for (itHGFEdgeLocal edge = node.edgesBegin(); edge != node.edgesEnd();
       ++edge) {
    curNetLens[(*edge)->getIndex()] = -1.;
  }
}

void CapoPlacer::addCellsNeedingUpdate(CapoBin* bin, BitBoard& needsUpdate) {
  // loop over bins neighbors above
  for (unsigned i = 0; i < bin->_neighborsAbove.size(); ++i) {
    CapoBin* neighbor = bin->_neighborsAbove[i];
    if (neighbor == NULL) continue;

    for (unsigned j = 0; j < neighbor->getCellIds().size(); ++j) {
      needsUpdate.setBit(neighbor->getCellIds()[j]);
    }
  }

  // loop over bins neighbors below
  for (unsigned i = 0; i < bin->_neighborsBelow.size(); ++i) {
    CapoBin* neighbor = bin->_neighborsBelow[i];
    if (neighbor == NULL) continue;

    for (unsigned j = 0; j < neighbor->getCellIds().size(); ++j) {
      needsUpdate.setBit(neighbor->getCellIds()[j]);
    }
  }

  // loop over bins left neighbors
  for (unsigned i = 0; i < bin->_leftNeighbors.size(); ++i) {
    CapoBin* neighbor = bin->_leftNeighbors[i];
    if (neighbor == NULL) continue;

    for (unsigned j = 0; j < neighbor->getCellIds().size(); ++j) {
      needsUpdate.setBit(neighbor->getCellIds()[j]);
    }
  }

  // loop over bins right neighbors
  for (unsigned i = 0; i < bin->_rightNeighbors.size(); ++i) {
    CapoBin* neighbor = bin->_rightNeighbors[i];
    if (neighbor == NULL) continue;

    for (unsigned j = 0; j < neighbor->getCellIds().size(); ++j) {
      needsUpdate.setBit(neighbor->getCellIds()[j]);
    }
  }
}

void CapoPlacer::addCellsNeedingUpdate(unsigned cellId, BitBoard& needsUpdate) {
  const HGraphWDimensions& hgraph = getNetlistHGraph();
  const HGFNode& node = hgraph.getNodeByIdx(cellId);

  for (itHGFEdgeLocal edge = node.edgesBegin(); edge != node.edgesEnd();
       ++edge) {
    for (itHGFNodeLocal adjNode = (*edge)->nodesBegin();
         adjNode != (*edge)->nodesEnd(); ++adjNode) {
      needsUpdate.setBit((*adjNode)->getIndex());
    }
  }
}

void CapoPlacer::updateCongestionEstimation() {
  Timer bloatTimer;

  double totalPinDensity = 0;
  const double size_factor = 50.;
  unsigned nonzeronodes = 0;

  const HGraphWDimensions& hgraph = getNetlistHGraph();
  BitBoard nodesVisited(hgraph.getNumNodes());

  for (itHGFNodeLocal nItr = hgraph.nodesBegin(); nItr != hgraph.nodesEnd();
       ++nItr) {
    nodesVisited.clear();
    const unsigned nodeId = (*nItr)->getIndex();

    if (hgraph.isTerminal(nodeId) || _cellToBinMap[nodeId] == NULL) {
      _congestion[nodeId] = 0;
      continue;
    }

    ++nonzeronodes;

    unsigned currentBinIndex = _cellToBinMap[nodeId]->getIndex();
    double currentCellArea = hgraph.getWeight(nodeId);
    double total_area = currentCellArea;
    unsigned total_pin = (*nItr)->getDegree();
    nodesVisited.setBit(nodeId);

    // get the one hop neighbors and their areas
    vector<unsigned> oneHopNeighbors;
    for (itHGFEdgeLocal eItr = (*nItr)->edgesBegin();
         eItr != (*nItr)->edgesEnd(); ++eItr) {
      for (itHGFNodeLocal hopItr = (*eItr)->nodesBegin();
           hopItr != (*eItr)->nodesEnd(); ++hopItr) {
        unsigned hopIndex = (*hopItr)->getIndex();

        if (nodesVisited.isBitSet(hopIndex)) continue;

        nodesVisited.setBit(hopIndex);

        unsigned hopBinIndex = (_cellToBinMap[hopIndex] == NULL)
                                   ? UINT_MAX
                                   : _cellToBinMap[hopIndex]->getIndex();

        if (hopBinIndex != currentBinIndex) continue;

        oneHopNeighbors.push_back(hopIndex);

        double hopArea = hgraph.getWeight(hopIndex);

        if (greaterOrEqualDouble(size_factor * currentCellArea, hopArea)) {
          total_pin += (*hopItr)->getDegree();
          total_area += hopArea;
        }
      }
    }

    // now we have a list of the one hop neighbors
    // we have also counted their area and pins
    // lets find the two hop neighbors and count their area and pins
    for (unsigned i = 0; i < oneHopNeighbors.size(); ++i) {
      const HGFNode& oneHop = hgraph.getNodeByIdx(oneHopNeighbors[i]);

      for (itHGFEdgeLocal eItr = oneHop.edgesBegin(); eItr != oneHop.edgesEnd();
           ++eItr) {
        for (itHGFNodeLocal hopItr = (*eItr)->nodesBegin();
             hopItr != (*eItr)->nodesEnd(); ++hopItr) {
          unsigned hopIndex = (*hopItr)->getIndex();

          if (nodesVisited.isBitSet(hopIndex)) continue;

          nodesVisited.setBit(hopIndex);

          unsigned hopBinIndex = (_cellToBinMap[hopIndex] == NULL)
                                     ? UINT_MAX
                                     : _cellToBinMap[hopIndex]->getIndex();

          if (hopBinIndex != currentBinIndex) continue;

          double hopArea = hgraph.getWeight(hopIndex);

          if (greaterOrEqualDouble(size_factor * currentCellArea, hopArea)) {
            total_pin += (*hopItr)->getDegree();
            total_area += hopArea;
          }
        }
      }
    }

    _congestion[nodeId] = static_cast<double>(total_pin) / total_area;
    totalPinDensity += _congestion[nodeId];
  }  // end node traversal

  double avgPinDensity = totalPinDensity / static_cast<double>(nonzeronodes);

  // normalize to the average pin density
  for (unsigned i = 0; i < hgraph.getNumNodes(); ++i) {
    _congestion[i] /= avgPinDensity;
  }

  vector<double> cong_temp = _congestion;
  unsigned position = cong_temp.size() - nonzeronodes +
                      static_cast<unsigned>(static_cast<double>(nonzeronodes) *
                                            _params.bloatCongPercentile) -
                      1;
  nth_element(cong_temp.begin(), cong_temp.begin() + position, cong_temp.end());
  _congestionCutoff = cong_temp[position];

  bloatTimer.stop();
  CapoPlacer::capoTimer::BloatTime += bloatTimer.getUserTime();
}

// direction: 1 is horizontal, 0 is vertical
bool CapoPlacer::isCutlinePredetermined(const string name, double& loc,
                                        bool& direction) {
  string::size_type constrStrLen = 0, nameLen = 0;
  bool foundAConstraint = false;
  nameLen = name.size();
  for (unsigned i = 0; i < _requiredSlicingTree.size(); ++i) {
    constrStrLen = _requiredSlicingTree[i].first.size();
    if (nameLen + 1 == constrStrLen &&
        _requiredSlicingTree[i].first.compare(0, nameLen, name) == 0) {
      foundAConstraint = true;
      loc = _requiredSlicingTree[i].second;
      if (_requiredSlicingTree[i].first[nameLen] == 'V')
        direction = 0;
      else
        direction = 1;
      break;
    }
  }
  return (foundAConstraint);
}

void CapoPlacer::readCutlineFile(const string cutlineFileName) {
  string binName;
  double location;
  ifstream cutlineFile(cutlineFileName.c_str());
  // abkwarn2(cutlineFile, "Could not open cutline file ",
  // cutlineFileName.c_str());

  if (cutlineFile) {
    while (!cutlineFile.eof()) {
      eatblank(cutlineFile);
      if (cutlineFile.eof())
        break;
      else {
        if (cutlineFile.peek() == '\n' || cutlineFile.peek() == '\r' ||
            cutlineFile.peek() == EOF) {
          cutlineFile.get();
          continue;
        }
      }
      cutlineFile >> binName >> location;
      _requiredSlicingTree.push_back(pair<string, double>(binName, location));
    }
    cutlineFile.close();
  }
}

double getRealNodeWidth(unsigned index, const PlacementWOrient& placement,
                        const HGraphWDimensions& hgraph) {
  unsigned angle = placement.getOrient(index).getAngle();
  bool notRotated = angle == 0 || angle == 180;
  double nodeWidth =
      notRotated ? hgraph.getNodeWidth(index) : hgraph.getNodeHeight(index);

  return nodeWidth;
}

double getRealNodeHeight(unsigned index, const PlacementWOrient& placement,
                         const HGraphWDimensions& hgraph) {
  unsigned angle = placement.getOrient(index).getAngle();
  bool notRotated = angle == 0 || angle == 180;
  double nodeHeight =
      notRotated ? hgraph.getNodeHeight(index) : hgraph.getNodeWidth(index);

  return nodeHeight;
}

void CapoPlacer::randomlyMoveCellsInTheirBins(void) {
  RandomRawDouble rng;

  BBox coreBBox = _rbplace.getBBox(false);
  double gridXSize =
      coreBBox.getWidth() / static_cast<double>(obstacleGrid.size());
  double gridYSize =
      coreBBox.getHeight() / static_cast<double>(obstacleGrid.size());

  for (unsigned i = 0; i < _placement.size(); ++i) {
    const CapoBin* bin = _cellToBinMap[i];

    if (bin == NULL) continue;

    for (unsigned j = 0; j < 3; ++j) {
      Point newPos;
      newPos.x = bin->getBBox().xMin +
                 rng * max(0., bin->getBBox().getWidth() -
                                   _rbplace.getHGraph().getNodeWidth(i));
      newPos.y = bin->getBBox().yMin +
                 rng * max(0., bin->getBBox().getHeight() -
                                   _rbplace.getHGraph().getNodeHeight(i));

      unsigned gridColIdx = min(
          obstacleGrid.size() - 1,
          static_cast<size_t>(floor((newPos.x - coreBBox.xMin) / gridXSize)));
      unsigned gridRowIdx = min(
          obstacleGrid.size() - 1,
          static_cast<size_t>(floor((newPos.y - coreBBox.yMin) / gridYSize)));

      bool legal = true;
      for (unsigned k = 0; k < obstacleGrid[gridColIdx][gridRowIdx].size();
           ++k) {
        unsigned obsIdx = obstacleGrid[gridColIdx][gridRowIdx][k];

        if (fixedObstacles[obsIdx].contains(newPos)) {
          legal = false;
          break;
        }
      }

      if (legal || j == 2) {
        updatePointSetsForMovedCell(i, _placement[i], newPos);
        _placement[i] = newPos;
        break;
      }
    }
  }
}
