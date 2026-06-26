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

#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <memory>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkfunc.h"
#include "ABKCommon/abktimer.h"
#include "FMPart/fmPart.h"
#include "FMPart/fmPartPlus.h"
#include "MLPart/mlPart.h"
#include "cutPosition.h"
#include "partProbForCapo.h"
#include "splitLarge.h"
#ifdef USEHMETIS
#include "HMetis/hmetisPart.h"
#endif

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

SplitLargeBinVertically::SplitLargeBinVertically(
    CapoBin& binToSplit, CapoPlacer& capo, CapoSplitterParams params,
    PlacementWOrient* placement, SplitPtFixedType splitPtFixedType,
    double xSplit)
    : BaseBinSplitter(binToSplit, capo, params) {
  if (_params.verb.getForActions() > 2)
    cout << "Splitting a large bin Vertically" << endl;

  Timer tm22;

  double p0Capacity = DBL_MAX;
  double origXSplit = xSplit;

  bool splitPtFixed = splitPtFixedType == SPLITPT_FIXED;
  if (splitPtFixedType != SPLITPT_UNFIXED) {
    BBox binBBox = _bin.getBBox();
    double siteSpacing = _bin.getRow(0).getSpacing();
    double givenXSplit =
        binBBox.xMin +
        rint((xSplit - binBBox.xMin) / siteSpacing) * siteSpacing;

    if (lessOrEqualDouble(givenXSplit, binBBox.xMin) ||
        greaterOrEqualDouble(givenXSplit, binBBox.xMax)) {
      abkwarn(0, "xSplit out of bin BBox while splitting vertically\n");
      cout << "xSplit " << xSplit << " Aligned xSplit : " << givenXSplit
           << ", Bin BBox " << binBBox;
      givenXSplit = rint((binBBox.xMin + binBBox.xMax) * 0.5);
      splitPtFixed = false;
    }

    BBox p0BBox = binBBox;
    p0BBox.xMax = givenXSplit;
    BBox p1BBox = binBBox;
    p1BBox.xMin = givenXSplit;
    p0Capacity = _bin.getContainedAreaInBBox(p0BBox);
    xSplit = givenXSplit;
  } else {
    double capacity = _bin.getCapacity();
    double targetCap0, targetCap1;

    vector<double> partCaps(2, 0.);

    CutLinePositionOracle* oracle = 0;
    bool oracleHasSpoken = false;
    double oracle_xSplit = DBL_MAX;
    if (_params.useCutOracle) {
      oracle = new CutLinePositionOracle(binToSplit);
      oracle_xSplit = oracle->getVerticalCut(oracleHasSpoken);
      double LeftWidth = (oracle_xSplit - binToSplit.getBBox().xMin);
      targetCap0 = capacity * LeftWidth / _bin.getBBox().getWidth();
      targetCap1 = capacity - targetCap0;
      partCaps[0] = targetCap0;
      partCaps[1] = targetCap1;
      p0Capacity = partCaps[0];
      // oracleHasSpoken = false; //HARD WIRE IGNORE ORACLE!!!
    }

    if (!_params.useCutOracle ||
        !oracleHasSpoken)  // Find the cutline the old way
    {
      if (_params.useCutOracle && _params.verb.getForActions() > 4)
        cerr << "Not using oracle's vertical cut suggestion" << endl;
      targetCap0 = 0.5 * capacity;
      targetCap1 = targetCap0;

      if (_params.fillerCellsPartCaps) {
        vector<double> targetCellAreas(2, 0.5 * _bin.getTotalCellArea());
        bool changeCaps = partAreasForFillerCells(1, targetCellAreas);
        if (changeCaps) {
          targetCap0 = capacity * targetCellAreas[0] /
                       (targetCellAreas[0] + targetCellAreas[1]);
          targetCap1 = capacity * targetCellAreas[1] /
                       (targetCellAreas[0] + targetCellAreas[1]);
        }
      }

      if (_params.verb.getForActions() > 4)
        cout << "Targets for calculating xSplit point are " << targetCap0
             << " : " << targetCap1 << endl;

      double p0WSWeight = 1.;
      double p1WSWeight = 1.;
      xSplit = _bin.findXSplit(targetCap0, targetCap1, 0, p0WSWeight,
                               p1WSWeight, partCaps);

      if (_params.verb.getForActions() > 4)
        cout << "Potential xSplit Pt: " << xSplit << endl;

      p0Capacity = partCaps[0];
    } else {
      xSplit = oracle_xSplit;
      if (_params.useCutOracle && _params.verb.getForActions() > 4)
        cerr << "Using oracle's cut suggestion. " << endl;
    }

    if (_params.useCutOracle && _params.verb.getForActions() > 3) {
      cerr << "Cutline Oracle Summary: " << endl;
      cerr << "\tSplitting large bin vertically" << endl;
      cerr << "\tCapo chose cutline position: " << xSplit << endl;
      cerr << "\tWhich has cut length: " << oracle->cutLineLengthV(xSplit)
           << endl;
      cerr << "\tOracle chose cutline position: " << oracle_xSplit << endl;
      cerr << "\tWhich has cut length: "
           << oracle->cutLineLengthV(oracle_xSplit) << endl;
    }
  }

  BaseBinSplitter& splitter = *this;
  bool isHCut = false;
  PartProbSizeType sizeLarge = Large;
  bool isRepart = false;
  double splitPt = xSplit;
  Verbosity verb = Verbosity("1 1 1");
  bool honorGrpConstr = _capo.getParams().useGrpConstr;
  unsigned splitRow = UINT_MAX;
  PartitioningProblemForCapo problem(
      _bin, _capo, _capo.partProbEdgesVisited, _capo.thetoWeights, splitter,
      isHCut, sizeLarge, isRepart, splitPt, splitPtFixed, splitRow, p0Capacity,
      verb, honorGrpConstr);
  tm22.stop();
  CapoPlacer::capoTimer::SetupTime += tm22.getUserTime();

  if (_params.eco && _bin.isEcoAllowed() &&
      (_params.ecoOverfillPct > _bin.getRelativeWhitespace())) {
    _bin.setEcoAllowed(false);
  }

  if (_params.eco && _bin.isEcoAllowed()) {
    // if it wasn't actually destroyed, this is a no-op
    _capo.rebuildHGraph();

    _capo.getParams().maxMem->update(
        "Split Large V after rebuilding the HGraph for eco");

    Timer ECOtimer;

    CapoBin::vCutLineInfo ecoSearchResults = _bin.scanCutLinesV(
        problem.getSubHGraph(), problem.getMinCapacities()[0][0],
        problem.getMaxCapacities()[0][0], problem.getMinCapacities()[1][0],
        problem.getMaxCapacities()[1][0]);

    double binWS = _bin.getRelativeWhitespace();
    double ecoXSplit = ecoSearchResults.xSplit;
    double shiftedXSplit = ecoSearchResults.xSplit;
    vector<double> partAreas(2, 0.), siteAreas(2, 0.);
    partAreas[0] = ecoSearchResults.p0Area;
    partAreas[1] = ecoSearchResults.p1Area;

    if (_params.ecoShift) {
      if (greaterOrEqualDouble(_params.safeWS, 0.) && binWS >= _params.safeWS) {
        // make sure each side of the split has as least this much WS
        shiftedXSplit =
            _bin.shiftXSplitAndSnap(ecoXSplit, _params.safeWS, partAreas[0],
                                    partAreas[1], 5.0, siteAreas);
      } else if (greaterOrEqualDouble(_params.minLocalWS, 0.) &&
                 binWS >= _params.minLocalWS) {
        // make sure each side of the split has as least this much WS
        shiftedXSplit =
            _bin.shiftXSplitAndSnap(ecoXSplit, _params.minLocalWS, partAreas[0],
                                    partAreas[1], 5.0, siteAreas);
      } else  // uniform WS
      {
        double p0WSWeight = 1.;
        double p1WSWeight = 1.;
        shiftedXSplit = _bin.findXSplit(partAreas[0], partAreas[1], 0,
                                        p0WSWeight, p1WSWeight, siteAreas);
      }
    } else {
      // populate siteAreas
      BBox p0BBox = _bin.getBBox();
      p0BBox.xMax = ecoXSplit;
      siteAreas[0] = _bin.getContainedAreaInBBox(p0BBox);
      siteAreas[1] = _bin.getCapacity() - siteAreas[0];
    }

    double ecoSplitCost = ecoSearchResults.cutCost;

    double p0WS = 1. - (partAreas[0] / siteAreas[0]);
    double p1WS = 1. - (partAreas[1] / siteAreas[1]);

    if (_params.ecoOverfillPct > p0WS || _params.ecoOverfillPct > p1WS) {
      // it looks like the eco split will cause WS issues,
      // so go from scratch
#ifdef USEFLOW
      if (!_params.flowECO) {
        // eco with flows is different, we can still use a later
        // solution if this one isn't all that good
        _bin.setEcoAllowed(false);
      }
#else
      _bin.setEcoAllowed(false);
#endif

      ECOtimer.stop();
      CapoPlacer::capoTimer::ECOTime += ECOtimer.getUserTime();
      ++CapoPlacer::EcoStats::ecoPartsRejected;

      callPartitioner(problem, placement);
    } else {
      Partitioning soln(problem.getHGraph().getNumNodes());

      for (unsigned n = 0; n < problem.getHGraph().getNumTerminals(); ++n) {
        unsigned fixed = problem.getFixedConstr()[n].getUnsigned();
        soln[n].loadBitsFrom(fixed & 3u);
      }

      for (unsigned n = problem.getHGraph().getNumTerminals();
           n < problem.getHGraph().getNumNodes(); ++n) {
        unsigned origID = problem.getSubHGraph().newNode2OrigIdx(n);
        unsigned angle = capo.getOraclePlacement().getOrient(origID).getAngle();
        bool notRotated = angle == 0 || angle == 180;
        double halfCellWidth =
            notRotated ? 0.5 * capo.getNetlistHGraph().getNodeWidth(origID)
                       : 0.5 * capo.getNetlistHGraph().getNodeHeight(origID);
        if (lessOrEqualDouble(
                capo.getOraclePlacement()[origID].x + halfCellWidth,
                ecoXSplit)) {
          soln[n].setToPart(0);
        } else {
          soln[n].setToPart(1);
        }
      }

      problem.reserveBuffers(1);
      problem.getSolnBuffers()[0] = soln;
      problem.setBestSolnNum(0);

      problem.setXSplit(shiftedXSplit);

      // lets see if we can improve upon the part we just generated
      FMPartitioner::Parameters fmParams(_capo.getParams().mlParams);
      fmParams.useEarlyStop = true;
      fmParams.maxNumPasses = 1;
      fmParams.useClip = false;
      fmParams.verb.setForActions(_params.verb.getForActions() / 2);
      fmParams.verb.setForMajStats(_params.verb.getForMajStats() / 2);
      fmParams.verb.setForSysRes(_params.verb.getForSysRes() / 2);
      fmParams.relaxedTolerancePasses = 0;

      std::unique_ptr<MultiStartPartitioner> partitioner;
      if (_params.useFMPartPlus) {
        partitioner.reset(new FMPartitionerPlus(
            problem, _capo.getParams().maxMem, fmParams, true));
      } else {
        partitioner.reset(new FMPartitioner(problem, _capo.getParams().maxMem,
                                            fmParams, true));
      }

      double newCut = partitioner->getBestSolnCost();

      if (newCut < (1. - _params.ecoImprovePct) * ecoSplitCost) {
        // it looks like the old split wasn't all that great
        // this means that the new one probably isn't all that
        // good either, so go from scratch
#ifdef USEFLOW
        if (!_params.flowECO) {
          // eco with flows is different, we can still use a later
          // solution if this one isn't all that good
          _bin.setEcoAllowed(false);
        }
#else
        _bin.setEcoAllowed(false);
#endif
        problem.setXSplit(splitPt);

        ECOtimer.stop();
        CapoPlacer::capoTimer::ECOTime += ECOtimer.getUserTime();
        ++CapoPlacer::EcoStats::ecoPartsRejected;

        callPartitioner(problem, placement);
      } else {
        // the old split probably wasn't so bad
        // keep the old and move on with it
        problem.getSolnBuffers()[0] = soln;

        if (ecoXSplit != shiftedXSplit) {
          // if we are keeping the ECO split and we had to shift
          // the split row, we need to scale the locations of
          // movables in the oracle placement
          shiftOraclePlacementV(ecoXSplit, shiftedXSplit);
        }

        ECOtimer.stop();
        CapoPlacer::capoTimer::ECOTime += ECOtimer.getUserTime();
        ++CapoPlacer::EcoStats::ecoPartsUsed;
      }
    }
  } else {
    callPartitioner(problem, placement);
    // <aaronnn> adjust cut of partitioning solution to resize bins by
    // utilization
    if (_params.adjustWSByDensity && !problem.getSplitPtFixed()) {
      problem.setXSplit(adjustSplitByUtilizationDensity(
          problem, problem.getXSplit(), isHCut));
    }
  }

  if (problem.getSplitPtFixed() && _params.verb.getForActions() > 4) {
    cout << "Requested xSplit " << origXSplit << " Aligned xSplit " << xSplit
         << " Final xSplit " << problem.getXSplit() << endl;
  }

  bool isLarge = true;
  createVSplitBins(isLarge, problem, problem.getXSplit());
}

void SplitLargeBinVertically::callPartitioner(
    PartitioningProblemForCapo& problem, PlacementWOrient* placement) {
  bool foundGoodSplit = false;
  unsigned numTries = 0;

  while (!foundGoodSplit && numTries++ < 2) {
    foundGoodSplit = callMLPart(problem, placement);

    if (foundGoodSplit) {
      if (_params.verb.getForActions() > 9) {
        cout << " Success.  Found a good split" << endl;
      }
    } else if (numTries < 2) {
      if (_params.verb.getForActions() > 9)
        cout << "  Failed. clearing part buffers" << endl;

      PartitioningBuffer& buff = problem.getSolnBuffers();
      unsigned numTerms = problem.getHGraph().getNumTerminals();

      unsigned sol;
      for (sol = buff.beginUsedSoln(); sol != buff.endUsedSoln(); sol++) {
        Partitioning& part = buff[sol];
        for (unsigned i = numTerms; i < part.size(); i++) part[i].setToAll(2);
      }

      // consider changing the problem here to attempt to
      // get a better solution
    } else if (_params.verb.getForActions() > 9)
      cout << "  Failed too many times. Taking the previous soln." << endl;
  }
}

bool SplitLargeBinVertically::callMLPart(PartitioningProblemForCapo& problem,
                                         PlacementWOrient* placement) {
  double binWS = _bin.getRelativeWhitespace();

  if (_params.verb.getForMajStats() > 4) {
    cout << " First MLPart" << endl;
    cout << " Targets: " << problem.getCapacities()[0][0] << "/"
         << problem.getCapacities()[1][0] << endl;
    cout << " Max:     " << problem.getMaxCapacities()[0][0] << "/"
         << problem.getMaxCapacities()[1][0] << endl;
    cout << " Min:     " << problem.getMinCapacities()[0][0] << "/"
         << problem.getMinCapacities()[1][0] << endl;
    cout << " Tolerance: " << 100 * problem.getTolerance(0) << endl;
  }

  MLPart::Parameters mlParams(_capo.getParams().mlParams);

  mlParams.verb.setForActions(_params.verb.getForActions() / 2);
  mlParams.verb.setForMajStats(_params.verb.getForMajStats() / 2);
  mlParams.verb.setForSysRes(_params.verb.getForSysRes() / 2);

  if (problem.getHGraph().getNumNodes() < 500)
    mlParams.clParams.sizeOfTop = 100;

  problem.reserveBuffers(2);  // sets are pairs of runs+a VC
  PartitioningBuffer refinedSolns(problem.getHGraph().getNumNodes(),
                                  _params.numMLSets + 1);
  vector<double> partAreas(2);  // in the BSF soln
  unsigned bestSet = 0;
  double origCut, secondBest, finalCut;
  origCut = secondBest = DBL_MAX;

  vector<string> nodeNames;
  if (_params.useNameHierarchy) {
    if (_params.verb.getForMajStats() > 2)
      cout << "Using name hierarchy to construct clustering" << endl;

    const HGraphFixed& netlistHG = _capo.getNetlistHGraph();
    vector<unsigned>::const_iterator cItr;

    for (cItr = _bin.cellIdsBegin(); cItr != _bin.cellIdsEnd(); cItr++)
      nodeNames.push_back(netlistHG.getNodeNameByIndex(*cItr).c_str());
  }
  FillableHierarchyParameters fillHParams;

  // problem.saveAsNodesNets((string("problem")+string(_bin.getName())).c_str());

  mlParams.Vcycling = MLPartParams::NoVcycles;
  for (unsigned setNum = 0; setNum <= _params.numMLSets; ++setNum) {
    if (setNum == _params.numMLSets) {
      if (setNum == 1) break;
      if (_params.verb.getForActions() > 1)
        cout << "Finished regularly scheduled ML sets: best cut is " << origCut
             << ", second best is " << secondBest << endl;
      if (lessOrEqualDouble(origCut, 0.)) {
        if (_params.verb.getForActions() > 1)
          cout << "Best cut optimal, skipping extra ML set" << endl;
        break;
      }
      double diffpercent = 100 * (secondBest - origCut) / origCut;
      if (_params.verb.getForActions() > 1)
        cout << "Percent difference between best 2 solns is " << diffpercent
             << "%" << endl;
      if (diffpercent > _params.extraBranchCutoff) {
        if (_params.verb.getForActions() > 1)
          cout << "Running one more ML set" << endl;
      } else {
        if (_params.verb.getForActions() > 1)
          cout << "Difference not greater than " << _params.extraBranchCutoff
               << "%, skipping extra ML set" << endl;
        break;
      }
    }

    double thisCut = -DBL_MAX;
    vector<double> thesePartAreas(2);
#ifdef USEHMETIS
    if (_params.useHMetis) {
      Timer hmetisTimer;

      HMetisParams hMetisParams;

      hMetisParams.verb.setForActions(_params.verb.getForActions() / 2);
      hMetisParams.verb.setForMajStats(_params.verb.getForMajStats() / 2);
      hMetisParams.verb.setForSysRes(_params.verb.getForSysRes() / 2);
      hMetisParams.startsPerSolution = 3;

      std::unique_ptr<HMetisPartitioner> metisPart(new HMetisPartitioner(
          problem, _capo.getParams().maxMem, hMetisParams));

      thisCut = metisPart->getBestSolnCost();
      thesePartAreas[0] = metisPart->getPartSizes()[0];
      thesePartAreas[1] = metisPart->getPartSizes()[1];

      hmetisTimer.stop();
      CapoPlacer::capoTimer::HMetisTime += hmetisTimer.getUserTime();
    } else
#endif
    {
      Timer mlPartTime;

      std::unique_ptr<MLPart> partitioner1;

      if (_params.useNameHierarchy) {
        std::vector<const char*> tmpNodeNames;
        string_vec_to_char_star_view(nodeNames, tmpNodeNames);
        FillableHierarchy fillH(
            tmpNodeNames, const_cast<char*>(_params.delimiters.c_str()),
            problem.getHGraph(), fillHParams);
        partitioner1.reset(new MLPart(problem, mlParams, _capo.bbBitBoards,
                                      _capo.getParams().maxMem, &fillH));
      } else if (_params.useQuadCluster &&
                 problem.getHGraph().getNumNodes() > 200) {
        partitioner1.reset(new MLPart(problem, mlParams, placement,
                                      _capo.bbBitBoards,
                                      _capo.getParams().maxMem));
      } else {
        partitioner1.reset(new MLPart(problem, mlParams, _capo.bbBitBoards,
                                      _capo.getParams().maxMem));
      }

      thisCut = partitioner1->getBestSolnCost();
      thesePartAreas[0] = partitioner1->getPartitionArea(0);
      thesePartAreas[1] = partitioner1->getPartitionArea(1);

      mlPartTime.stop();
      CapoPlacer::capoTimer::MLTime += mlPartTime.getUserTime();
    }

    refinedSolns[setNum] = problem.getBestSoln();
    if (thisCut < origCut) {
      secondBest = origCut;
      origCut = thisCut;
      bestSet = setNum;
      partAreas[0] = thesePartAreas[0];
      partAreas[1] = thesePartAreas[1];
    } else if (thisCut < secondBest) {
      secondBest = thisCut;
    }

    PartitionIds clearId;
    clearId.setToAll(2);
    Partitioning& buff0 = problem.getSolnBuffers()[0];
    Partitioning& buff1 = problem.getSolnBuffers()[1];
    std::fill(buff0.begin() + problem.getHGraph().getNumTerminals(),
              buff0.end(), clearId);
    std::fill(buff1.begin() + problem.getHGraph().getNumTerminals(),
              buff1.end(), clearId);
  }

  problem.reserveBuffers(1);
  problem.getSolnBuffers()[0] = refinedSolns[bestSet];
  problem.setBestSolnNum(0);

  // do some VCycling if useful (i.e. problem large enough)
  unsigned sizeOfTop = mlParams.clParams.sizeOfTop;
  if (_capo.getParams().mlParams.Vcycling != MLPartParams::NoVcycles &&
      sizeOfTop != 0 &&
      problem.getHGraph().getNumNodes() -
              problem.getHGraph().getNumTerminals() >
          sizeOfTop &&
      origCut > 100.) {
    if (_params.verb.getForActions() > 4) {
      cout << " Vcycling best solution" << endl;
      cout << " Targets: " << problem.getCapacities()[0][0] << "/"
           << problem.getCapacities()[1][0] << endl;
      cout << " Max:     " << problem.getMaxCapacities()[0][0] << "/"
           << problem.getMaxCapacities()[1][0] << endl;
      cout << " Min:     " << problem.getMinCapacities()[0][0] << "/"
           << problem.getMinCapacities()[1][0] << endl;
      cout << " Tolerance: " << 100 * problem.getTolerance(0) << endl;
    }
    MLPart::Parameters vCycleParams = mlParams;
    vCycleParams.Vcycling = MLPartParams::Initial;
    Timer vcyclingtime;
    MLPart vcycler(problem, vCycleParams, _capo.bbBitBoards,
                   _capo.getParams().maxMem);
    origCut = vcycler.getBestSolnCost();
    partAreas[0] = vcycler.getPartitionArea(0);
    partAreas[1] = vcycler.getPartitionArea(1);
    vcyclingtime.stop();
    CapoPlacer::capoTimer::MLTime += vcyclingtime.getUserTime();
  }

  vector<double> siteAreas(2, 0.);
  double oldXSplit = problem.getXSplit();

  if (!problem.getSplitPtFixed()) {
    if (_params.WSOracle) {
      problem.setXSplit(_bin.shiftXSplit(oldXSplit, 0, partAreas[0],
                                         partAreas[1], siteAreas));
    } else if (greaterOrEqualDouble(_params.safeWS, 0.) &&
               binWS >= _params.safeWS) {
      // make sure each side of the split has as least this much WS
      // problem.setXSplit(_bin.shiftXSplit(oldXSplit, _params.safeWS,
      // partAreas[0], partAreas[1], siteAreas));
      problem.setXSplit(_bin.shiftXSplitAndSnap(oldXSplit, _params.safeWS,
                                                partAreas[0], partAreas[1], 5.0,
                                                siteAreas));
    } else if (greaterOrEqualDouble(_params.minLocalWS, 0.) &&
               binWS >= _params.minLocalWS) {
      // make sure each side of the split has as least this much WS
      // problem.setXSplit(_bin.shiftXSplit(oldXSplit, _params.minLocalWS,
      // partAreas[0], partAreas[1], siteAreas));
      problem.setXSplit(_bin.shiftXSplitAndSnap(oldXSplit, _params.minLocalWS,
                                                partAreas[0], partAreas[1], 5.0,
                                                siteAreas));
    } else  // uniform WS
    {
      double p0WSWeight = 1.;
      double p1WSWeight = 1.;
      problem.setXSplit(_bin.findXSplit(partAreas[0], partAreas[1], 0,
                                        p0WSWeight, p1WSWeight, siteAreas));
    }
  }

  if (!_params.doRepartitioning &&
      ((_params.repartShiftCL && problem.getXSplit() == oldXSplit) ||
       !_params.repartShiftCL) &&
      !problem.getRepartSmallWS() && !problem.getSplitPtFixed())  // we're done
  {
    finalCut = origCut;
    if (_params.verb.getForMajStats() > 4)
      cout << " Best Solution had cut " << finalCut << endl;
    return true;
  }

  // double idealXSplit = problem.getXSplit();
  if (!problem.getSplitPtFixed()) {
    if (!(greaterOrEqualDouble(_params.safeWS, 0.) &&
          binWS >= _params.safeWS) &&
        !(greaterOrEqualDouble(_params.minLocalWS, 0.) &&
          binWS >= _params.minLocalWS)) {
      double p0WSWeight = 1.;
      double p1WSWeight = 1.;
      problem.setXSplit(_bin.findXSplit(partAreas[0], partAreas[1], 5,
                                        p0WSWeight, p1WSWeight, siteAreas));
    }
  } else {
    // the following is here just to populate siteAreas
    _bin.computePartAreas(problem.getXSplit(), siteAreas);
  }

  if (problem.getXSplit() == oldXSplit && !problem.getRepartSmallWS() &&
      !problem.getSplitPtFixed())  // no change...
  {
    finalCut = origCut;
    if (_params.verb.getForMajStats() > 4)
      cout << " Best Solution had cut " << finalCut << endl;
    return true;
  }

  // cutline was shifted, repartition w/ the new
  // targets and new tolerance.

  if (_params.verb.getForMajStats() > 1)
    cout << "Doing Repartitioning phase" << endl;

  problem.destroyHGraph();

  // if it wasn't actually destroyed, this is a no-op
  _capo.rebuildHGraph();

  _capo.getParams().maxMem->update(
      "Split Large V after rebuilding the HGraph for repartitioning");

  Timer repartProb;
  bool honorGrpConstraints = _capo.getParams().useGrpConstr;
  unsigned splitRow = UINT_MAX;
  PartitioningProblemForCapo problemNew(
      _bin, _capo, _capo.partProbEdgesVisited, _capo.thetoWeights, *this, false,
      Large, true, problem.getXSplit(), problem.getSplitPtFixed(), splitRow,
      siteAreas[0], Verbosity("1 1 1"), honorGrpConstraints);

  PartitioningBuffer& bufOld = problem.getSolnBuffers();
  unsigned bestSolnNum = problem.getBestSolnNum();
  PartitioningBuffer* pNewBuf = problemNew.swapOutSolnBuffers(&bufOld);
  Partitioning& pOldFC = const_cast<Partitioning&>(problem.getFixedConstr());
  Partitioning* pNewFC = problemNew.swapOutFixedConst(&pOldFC);
  problemNew.setBestSolnNum(bestSolnNum);
  problemNew.getSolnBuffers().setBeginUsedSoln(bestSolnNum);
  problemNew.getSolnBuffers().setEndUsedSoln(bestSolnNum + 1);

  repartProb.stop();
  CapoPlacer::capoTimer::SetupTime += repartProb.getUserTime();

  FMPartitioner::Parameters fmParams(_capo.getParams().mlParams);
  fmParams.useEarlyStop = false;
  fmParams.useClip = false;
  fmParams.verb.setForActions(_params.verb.getForActions() / 2);
  fmParams.verb.setForMajStats(_params.verb.getForMajStats() / 2);
  fmParams.verb.setForSysRes(_params.verb.getForSysRes() / 2);
  fmParams.relaxedTolerancePasses = 0;

  PartitioningBuffer bestSolns(problemNew.getHGraph().getNumNodes(), 2);
  double bestCuts[2];
  double bestSplits[2];
  bool satisfiesWS[2];

  bestSolns[0] = problem.getBestSoln();
  bestSolns[1] = problem.getBestSoln();
  bestCuts[0] = origCut;
  bestCuts[1] = origCut;
  bestSplits[0] = problem.getXSplit();
  bestSplits[1] = problem.getXSplit();

  double p0WS = (siteAreas[0] - partAreas[0]) / siteAreas[0];
  double p1WS = (siteAreas[1] - partAreas[1]) / siteAreas[1];
  if (greaterOrEqualDouble(_params.safeWS, 0.) && binWS >= _params.safeWS) {
    satisfiesWS[0] = p0WS >= _params.safeWS && p1WS >= _params.safeWS;
  } else if (greaterOrEqualDouble(_params.minLocalWS, 0.) &&
             binWS >= _params.minLocalWS) {
    satisfiesWS[0] = p0WS >= _params.minLocalWS && p1WS >= _params.minLocalWS;
  } else {
    satisfiesWS[0] = p0WS >= 0. && p1WS >= 0.;
  }

  // override the problemNew tolerance calculation
  double oldTol = problem.getTolerance();
  double newTol = problemNew.getTolerance();
  if (oldTol != newTol) {
    if (lessOrEqualDouble(newTol, 0.)) {
      problemNew.setTolerance(oldTol);
    } else {
      problemNew.multiplyTolerance(oldTol / newTol);
    }
  }

  double smallestTol = _bin.getMedianCellArea() / _bin.getTotalCellArea(),
         lastP0Max, lastP0Min, lastP1Max, lastP1Min;

  if (_params.verb.getForMajStats() > 4)
    cout << "Minimum repartitioning tolerance: " << 100. * smallestTol << "%"
         << endl;

  bool doTighten = false;
  bool hitBottom = false;
  if (problemNew.getTolerance() <= smallestTol) {
    lastP0Max = problemNew.getMaxCapacities()[0][0];
    lastP0Min = problemNew.getMinCapacities()[0][0];
    lastP1Max = problemNew.getMaxCapacities()[1][0];
    lastP1Min = problemNew.getMinCapacities()[1][0];
    // tighten at least once
    problemNew.multiplyTolerance(0.5);
    doTighten = true;
  } else {
    do {
      lastP0Max = problemNew.getMaxCapacities()[0][0];
      lastP0Min = problemNew.getMinCapacities()[0][0];
      lastP1Max = problemNew.getMaxCapacities()[1][0];
      lastP1Min = problemNew.getMinCapacities()[1][0];
      problemNew.multiplyTolerance(0.5);
      if (problemNew.getTolerance() <= smallestTol) {
        problemNew.increaseTolerance(smallestTol);
        hitBottom = true;
      }
      if ((problemNew.getMinCapacities()[0][0] > partAreas[0]) ||
          (partAreas[0] > problemNew.getMaxCapacities()[0][0]) ||
          (problemNew.getMinCapacities()[1][0] > partAreas[1]) ||
          (partAreas[1] > problemNew.getMaxCapacities()[1][0])) {
        // current solution does not satisfy these constraints
        // so we should try this
        doTighten = true;
      }
    } while (!doTighten && !hitBottom);
  }
  doTighten = doTighten || !satisfiesWS[0];

  if (doTighten) {
    hitBottom = false;
    bool firstPass = true;
    while (true) {
      if (_params.verb.getForMajStats() > 4) {
        cout << " Tightening pass" << endl;
        cout << " Targets: " << problemNew.getCapacities()[0][0] << "/"
             << problemNew.getCapacities()[1][0] << endl;
        cout << " Max:     " << problemNew.getMaxCapacities()[0][0] << "/"
             << problemNew.getMaxCapacities()[1][0] << endl;
        cout << " Min:     " << problemNew.getMinCapacities()[0][0] << "/"
             << problemNew.getMinCapacities()[1][0] << endl;
        cout << " Tolerance: " << 100 * problemNew.getTolerance(0) << endl;
      }

      Timer tm9;
      std::unique_ptr<MultiStartPartitioner> partitioner2;
      if (_params.useFMPartPlus) {
        partitioner2.reset(new FMPartitionerPlus(
            problemNew, _capo.getParams().maxMem, fmParams, true));
      } else {
        partitioner2.reset(new FMPartitioner(
            problemNew, _capo.getParams().maxMem, fmParams, true));
      }
      tm9.stop();
      CapoPlacer::capoTimer::FMTime += tm9.getUserTime();

      double newCut = partitioner2->getBestSolnCost();
      partAreas[0] = 0.;
      partAreas[1] = 0.;
      for (unsigned c = 0; c < problemNew.getHGraph().getNumNodes(); ++c) {
        unsigned part = problemNew.getBestSoln()[c].isInPart(0) ? 0 : 1;
        partAreas[part] += problemNew.getHGraph().getWeight(c);
      }

      if (newCut <= 1.02 * origCut || (firstPass && !satisfiesWS[0])) {
        // check to see if we have gotten to the min tolerance
        if (problemNew.getTolerance() <= (smallestTol + 1e-6) || hitBottom) {
          if (_params.verb.getForActions() > 4) {
            cout << " Vcycling good tightened solution" << endl;
            cout << " Targets: " << problemNew.getCapacities()[0][0] << "/"
                 << problemNew.getCapacities()[1][0] << endl;
            cout << " Max:     " << problemNew.getMaxCapacities()[0][0] << "/"
                 << problemNew.getMaxCapacities()[1][0] << endl;
            cout << " Min:     " << problemNew.getMinCapacities()[0][0] << "/"
                 << problemNew.getMinCapacities()[1][0] << endl;
            cout << " Tolerance: " << 100 * problemNew.getTolerance(0) << endl;
          }
          MLPart::Parameters vCycleParams = mlParams;
          vCycleParams.Vcycling = MLPartParams::Initial;
          Timer vcyclingtime;
          MLPart vcycler(problemNew, vCycleParams, _capo.bbBitBoards,
                         _capo.getParams().maxMem);
          bestCuts[1] = vcycler.getBestSolnCost();
          bestSolns[1] = problemNew.getBestSoln();
          vcyclingtime.stop();
          CapoPlacer::capoTimer::MLTime += vcyclingtime.getUserTime();
          break;
        }
        // if not, reduce the tolerance further and try again
        else {
          doTighten = false;
          do {
            lastP0Max = problemNew.getMaxCapacities()[0][0];
            lastP0Min = problemNew.getMinCapacities()[0][0];
            lastP1Max = problemNew.getMaxCapacities()[1][0];
            lastP1Min = problemNew.getMinCapacities()[1][0];
            problemNew.multiplyTolerance(0.5);
            if (problemNew.getTolerance() <= smallestTol) {
              problemNew.increaseTolerance(smallestTol);
              hitBottom = true;
            }
            if ((problemNew.getMinCapacities()[0][0] > partAreas[0]) ||
                (partAreas[0] > problemNew.getMaxCapacities()[0][0]) ||
                (problemNew.getMinCapacities()[1][0] > partAreas[1]) ||
                (partAreas[1] > problemNew.getMaxCapacities()[1][0])) {
              // current solution does not satisfy these constraints
              // so we should try this
              doTighten = true;
            }
          } while (!doTighten && !hitBottom);
          if (!doTighten) {
            if (_params.verb.getForActions() > 4) {
              cout << " Vcycling good tightened solution" << endl;
              cout << " Targets: " << problemNew.getCapacities()[0][0] << "/"
                   << problemNew.getCapacities()[1][0] << endl;
              cout << " Max:     " << problemNew.getMaxCapacities()[0][0] << "/"
                   << problemNew.getMaxCapacities()[1][0] << endl;
              cout << " Min:     " << problemNew.getMinCapacities()[0][0] << "/"
                   << problemNew.getMinCapacities()[1][0] << endl;
              cout << " Tolerance: " << 100 * problemNew.getTolerance(0)
                   << endl;
            }
            MLPart::Parameters vCycleParams = mlParams;
            vCycleParams.Vcycling = MLPartParams::Initial;
            Timer vcyclingtime;
            MLPart vcycler(problemNew, vCycleParams, _capo.bbBitBoards,
                           _capo.getParams().maxMem);
            bestCuts[1] = vcycler.getBestSolnCost();
            bestSolns[1] = problemNew.getBestSoln();
            vcyclingtime.stop();
            CapoPlacer::capoTimer::MLTime += vcyclingtime.getUserTime();
            break;
          }
        }
      } else {
        if (_params.verb.getForActions() > 4) {
          cout << " Vcycling tightened solution with 2% increased cut" << endl;
          cout << " Targets: " << problemNew.getCapacities()[0][0] << "/"
               << problemNew.getCapacities()[1][0] << endl;
          cout << " Max:     " << problemNew.getMaxCapacities()[0][0] << "/"
               << problemNew.getMaxCapacities()[1][0] << endl;
          cout << " Min:     " << problemNew.getMinCapacities()[0][0] << "/"
               << problemNew.getMinCapacities()[1][0] << endl;
          cout << " Tolerance: " << 100 * problemNew.getTolerance(0) << endl;
        }
        MLPart::Parameters vCycleParams = mlParams;
        vCycleParams.Vcycling = MLPartParams::Initial;
        Timer vcyclingtime;
        MLPart vcycler(problemNew, vCycleParams, _capo.bbBitBoards,
                       _capo.getParams().maxMem);
        double improvedCut = vcycler.getBestSolnCost();
        vcyclingtime.stop();
        CapoPlacer::capoTimer::MLTime += vcyclingtime.getUserTime();

        if (improvedCut <= 1.04 * origCut) {
          bestCuts[1] = improvedCut;
          bestSolns[1] = problemNew.getBestSoln();
          break;
        } else {
          problemNew.getSolnBuffers()[bestSolnNum] = bestSolns[1];
          problemNew.setBestSolnNum(bestSolnNum);
          problemNew.setMaxCapacity(0, lastP0Max);
          problemNew.setMinCapacity(0, lastP0Min);
          problemNew.setMaxCapacity(1, lastP1Max);
          problemNew.setMinCapacity(1, lastP1Min);
          if (_params.verb.getForActions() > 4) {
            cout << " Vcycling previous tightened solution" << endl;
            cout << " Targets: " << problemNew.getCapacities()[0][0] << "/"
                 << problemNew.getCapacities()[1][0] << endl;
            cout << " Max:     " << problemNew.getMaxCapacities()[0][0] << "/"
                 << problemNew.getMaxCapacities()[1][0] << endl;
            cout << " Min:     " << problemNew.getMinCapacities()[0][0] << "/"
                 << problemNew.getMinCapacities()[1][0] << endl;
            cout << " Tolerance: " << 100 * problemNew.getTolerance(0) << endl;
          }
          MLPart::Parameters vCycleParams = mlParams;
          vCycleParams.Vcycling = MLPartParams::Initial;
          Timer vcyclingtime;
          MLPart vcycler(problemNew, vCycleParams, _capo.bbBitBoards,
                         _capo.getParams().maxMem);
          bestCuts[1] = vcycler.getBestSolnCost();
          bestSolns[1] = problemNew.getBestSoln();
          vcyclingtime.stop();
          CapoPlacer::capoTimer::MLTime += vcyclingtime.getUserTime();
          break;
        }
      }
      bestCuts[1] = newCut;
      bestSolns[1] = problemNew.getBestSoln();
      firstPass = false;
    }
  }

  partAreas[0] = 0;
  partAreas[1] = 0;
  for (unsigned c = 0; c < problemNew.getHGraph().getNumNodes(); ++c) {
    unsigned part = bestSolns[1][c].isInPart(0) ? 0 : 1;
    partAreas[part] += problemNew.getHGraph().getWeight(c);
  }

  if (!problemNew.getSplitPtFixed()) {
    if (_params.WSOracle) {
      bestSplits[1] = _bin.shiftXSplit(bestSplits[1], 0, partAreas[0],
                                       partAreas[1], siteAreas);
    } else if (greaterOrEqualDouble(_params.safeWS, 0.) &&
               binWS >= _params.safeWS) {
      // make sure each side of the split has as least this much WS
      // bestSplits[1] = _bin.shiftXSplit(bestSplits[1], _params.safeWS,
      // partAreas[0], partAreas[1], siteAreas);
      bestSplits[1] =
          _bin.shiftXSplitAndSnap(bestSplits[1], _params.safeWS, partAreas[0],
                                  partAreas[1], 5.0, siteAreas);
    } else if (greaterOrEqualDouble(_params.minLocalWS, 0.) &&
               binWS >= _params.minLocalWS) {
      // make sure each side of the split has as least this much WS
      // bestSplits[1] = _bin.shiftXSplit(bestSplits[1], _params.minLocalWS,
      // partAreas[0], partAreas[1], siteAreas);
      bestSplits[1] =
          _bin.shiftXSplitAndSnap(bestSplits[1], _params.minLocalWS,
                                  partAreas[0], partAreas[1], 5.0, siteAreas);
    } else {
      double p0WSWeight = 1.;
      double p1WSWeight = 1.;
      bestSplits[1] = _bin.findXSplit(partAreas[0], partAreas[1], 0, p0WSWeight,
                                      p1WSWeight, siteAreas);
    }
  }

  // check to see if this satisfies the min local constraints
  p0WS = (siteAreas[0] - partAreas[0]) / siteAreas[0];
  p1WS = (siteAreas[1] - partAreas[1]) / siteAreas[1];
  if (greaterOrEqualDouble(_params.safeWS, 0.) && binWS >= _params.safeWS) {
    satisfiesWS[1] = p0WS >= _params.safeWS && p1WS >= _params.safeWS;
  } else if (greaterOrEqualDouble(_params.minLocalWS, 0.) &&
             binWS >= _params.minLocalWS) {
    satisfiesWS[1] = p0WS >= _params.minLocalWS && p1WS >= _params.minLocalWS;
  } else {
    satisfiesWS[1] = p0WS >= 0. && p1WS >= 0.;
  }

  /*
  // possible decision making scheme
  if(satisfiesWS[0])
  {
  if(satisfiesWS[1])
  {
  if(bestCuts[0] < bestCuts[1])
  {
  abkwarn(0,"Dropping tightened solution because both ok wrt local ws and orig
  had better cut"); problemNew.getSolnBuffers()[0] = bestSolns[0];
  problemNew.setBestSolnNum(0);
  problemNew.setXSplit(bestSplits[0]);
  finalCut = bestCuts[0];
  }
  else
  {
  problemNew.setXSplit(bestSplits[1]);
  finalCut = bestCuts[1];
  }
  }
  else
  {
  abkwarn(0,"Original solution satisfied WS, but tightened did not. Probably a
  bug."); problemNew.getSolnBuffers()[0] = bestSolns[0];
  problemNew.setBestSolnNum(0);
  problemNew.setXSplit(bestSplits[0]);
  finalCut = bestCuts[0];
  }
  }
  else
  {
  if(satisfiesWS[1])
  {
  problemNew.setXSplit(bestSplits[1]);
  finalCut = bestCuts[1];
  }
  else
  {
  if(bestCuts[0] < bestCuts[1])
  {
  abkwarn(0,"Dropping tightened solution because both bad wrt local ws and orig
  had better cut"); problemNew.getSolnBuffers()[0] = bestSolns[0];
  problemNew.setBestSolnNum(0);
  problemNew.setXSplit(bestSplits[0]);
  finalCut = bestCuts[0];
  }
  else
  {
  problemNew.setXSplit(bestSplits[1]);
  finalCut = bestCuts[1];
  }
  }
  }
   */

  problemNew.setXSplit(bestSplits[1]);
  finalCut = bestCuts[1];

  unsigned bestSolnNumNew = problemNew.getBestSolnNum();
  problemNew.swapOutSolnBuffers(pNewBuf);
  problemNew.swapOutFixedConst(pNewFC);
  problem.setBestSolnNum(bestSolnNumNew);
  problem.setXSplit(problemNew.getXSplit());
  problem.swapHGraphs(problemNew);

  if (_params.verb.getForMajStats() > 4)
    cout << " Best Solution had cut " << finalCut << endl;

  return true;
}
