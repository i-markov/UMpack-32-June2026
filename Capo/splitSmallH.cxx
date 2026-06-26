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

#include "FMPart/fmPart.h"
#include "FMPart/fmPartPlus.h"
#include "MLPart/mlPart.h"
#include "ShellPart/bbFMPart.h"
#include "partProbForCapo.h"
#include "splitSmall.h"

#include <memory>

using std::cout;
using std::endl;
using std::string;
using std::vector;

SplitSmallBinHorizontally::SplitSmallBinHorizontally(CapoBin& binToSplit,
                                                     CapoPlacer& capo,
                                                     CapoSplitterParams params)
    : BaseBinSplitter(binToSplit, capo, params) {
  abkfatal(_bin.getNumRows() > 1, "can't split a 1 row bin horizontally");

  if (_params.verb.getForActions() > 2)
    cout << "Splitting a small bin horizontally" << endl;

  // notes:

  // shifting the cutline will be difficult on such small
  // problems, so we won't allow as much tolerance.
  // also, with the small number of rows small problems have,
  // shifting by even 1 row would mean a large change in
  // tolerance

  // P0 is the top partition

  Timer tm24;

  const vector<const RBPCoreRow*>& rows = _bin.getRows();
  double totalCellArea = _bin.getTotalCellArea();
  vector<double> targetCellAreas(2, 0.5 * totalCellArea);
  vector<double> actualCaps(2, 0);

  double p0WSWeight = 1.;
  double p1WSWeight = 1.;
  unsigned splitRow =
      findBestSplitRow(targetCellAreas, p0WSWeight, p1WSWeight, actualCaps);

  if (_params.verb.getForActions() > 4) {
    cout << "P0/P1 rows: " << rows.size() - splitRow << "/" << splitRow << endl;
  }

  bool honorGrpConstraints = _capo.getParams().useGrpConstr;

  PartProbSizeType Size =
      (_bin.getNumCells() > _capo.getParams().smallPartThreshold) ? Medium
                                                                  : Small;
  bool splitPtFixed = false;
  PartitioningProblemForCapo problem(
      _bin, _capo, _capo.partProbEdgesVisited, _capo.thetoWeights, *this, true,
      Size, false, rows[splitRow]->getYCoord(), splitPtFixed, splitRow,
      actualCaps[0], Verbosity("1 1 1"), honorGrpConstraints);
  tm24.stop();
  CapoPlacer::capoTimer::SetupTime += tm24.getUserTime();

  if (_params.eco && _bin.isEcoAllowed() &&
      (_params.ecoOverfillPct > _bin.getRelativeWhitespace())) {
    _bin.setEcoAllowed(false);
  }

  bool shouldPartByGroups = shouldPartitionByGroups(problem);
  if (shouldPartByGroups) _bin.setEcoAllowed(false);

  if (_params.eco && _bin.isEcoAllowed()) {
    Timer ECOtimer;

    CapoBin::hCutLineInfo ecoSearchResults = _bin.scanCutLinesH(
        problem.getSubHGraph(), problem.getMinCapacities()[0][0],
        problem.getMaxCapacities()[0][0], problem.getMinCapacities()[1][0],
        problem.getMaxCapacities()[1][0]);

    double binWS = _bin.getRelativeWhitespace();
    unsigned ecoSplitRow = ecoSearchResults.splitRow;
    unsigned shiftedSplitRow = ecoSplitRow;
    vector<double> partAreas(2, 0.), actualCaps(2, 0.);
    partAreas[0] = ecoSearchResults.p0Area;
    partAreas[1] = ecoSearchResults.p1Area;

    if (_params.ecoShift) {
      if (greaterOrEqualDouble(_params.safeWS, 0.) && binWS >= _params.safeWS) {
        // make sure each side of the split has as least this much WS
        shiftedSplitRow = shiftSplitRowAndSnap(ecoSplitRow, partAreas,
                                               actualCaps, _params.safeWS, 5.0);
      } else if (greaterOrEqualDouble(_params.minLocalWS, 0.) &&
                 binWS >= _params.minLocalWS) {
        // make sure each side of the split has as least this much WS
        shiftedSplitRow = shiftSplitRowAndSnap(
            ecoSplitRow, partAreas, actualCaps, _params.minLocalWS, 5.0);
      } else  // uniform WS
      {
        double p0WSWeight = 1.;
        double p1WSWeight = 1.;
        shiftedSplitRow =
            findBestSplitRow(partAreas, p0WSWeight, p1WSWeight, actualCaps);
      }
    } else {
      // populate actualCaps
      BBox p0BBox = _bin.getBBox();
      p0BBox.yMin = _bin.getRows()[ecoSplitRow]->getYCoord();
      actualCaps[0] = _bin.getContainedAreaInBBox(p0BBox);
      actualCaps[1] = _bin.getCapacity() - actualCaps[0];
    }

    double ecoSplitRowCoord = _bin.getRows()[ecoSplitRow]->getYCoord();
    double ecoSplitCost = ecoSearchResults.cutCost;

    double p0WS = 1. - (partAreas[0] / actualCaps[0]);
    double p1WS = 1. - (partAreas[1] / actualCaps[1]);

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

      callPartitioner(problem);
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
        double halfCellHeight =
            notRotated ? 0.5 * capo.getNetlistHGraph().getNodeHeight(origID)
                       : 0.5 * capo.getNetlistHGraph().getNodeWidth(origID);
        if (greaterOrEqualDouble(
                capo.getOraclePlacement()[origID].y + halfCellHeight,
                ecoSplitRowCoord)) {
          soln[n].setToPart(0);
        } else {
          soln[n].setToPart(1);
        }
      }

      problem.reserveBuffers(1);
      problem.getSolnBuffers()[0] = soln;
      problem.setBestSolnNum(0);

      problem.setSplitRow(shiftedSplitRow);

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
        problem.setSplitRow(splitRow);

        ECOtimer.stop();
        CapoPlacer::capoTimer::ECOTime += ECOtimer.getUserTime();
        ++CapoPlacer::EcoStats::ecoPartsRejected;

        callPartitioner(problem);
      } else {
        // the old split probably wasn't so bad
        // keep the old and move on with it
        problem.getSolnBuffers()[0] = soln;

        if (ecoSplitRow != shiftedSplitRow) {
          // if we are keeping the ECO split and we had to shift
          // the split row, we need to scale the locations of
          // movables in the oracle placement
          shiftOraclePlacementH(ecoSplitRowCoord,
                                _bin.getRows()[shiftedSplitRow]->getYCoord());
        }

        ECOtimer.stop();
        CapoPlacer::capoTimer::ECOTime += ECOtimer.getUserTime();
        ++CapoPlacer::EcoStats::ecoPartsUsed;
      }
    }
  } else {
    callPartitioner(problem);
    if (_params.adjustWSByDensity) {
      bool isHCut = true;
      double splitPt = adjustSplitByUtilizationDensity(
          problem, _bin.getRows()[splitRow]->getYCoord(), isHCut);

      for (unsigned newSplitRow = 0; newSplitRow < rows.size() - 1;
           newSplitRow++) {
        double rowTop = _bin.getRows()[newSplitRow]->getYCoord() +
                        _bin.getRows()[newSplitRow]->getHeight();
        if (rowTop > splitPt) {
          splitRow = newSplitRow;
          problem.setSplitRow(splitRow);
          break;
        }
      }
    }
  }

  if (_params.verb.getForActions() > 8)
    cout << "Done partitioning...creating bins" << endl;

  bool isLarge = false;
  createHSplitBins(isLarge, problem, problem.getSplitRow());
}

void SplitSmallBinHorizontally::callPartitioner(
    PartitioningProblemForCapo& problem) {
  bool foundGoodSplit = false;
  unsigned numTries = 0;

  const Partitioning& fixedConstr = problem.getFixedConstr();
  const HGraphFixed& hgraph = problem.getHGraph();
  bool atleast1Movable = false;
  PartitionIds part01;
  part01.setToAll(2);
  for (unsigned k = 0; k < hgraph.getNumNodes(); ++k) {
    if (fixedConstr[k].isEmpty() || fixedConstr[k] == part01) {
      atleast1Movable = true;
      break;
    }
  }

  while (!foundGoodSplit && numTries++ < 2) {
    if (_bin.getNumCells() <= _capo.getParams().smallPartThreshold &&
        atleast1Movable) {
      if (_params.verb.getForActions() > 8) cout << "Using BBPart" << endl;

      foundGoodSplit = callBBFM(problem);
    } else {
      if (_params.verb.getForActions() > 8) cout << "Using FMPart" << endl;

      foundGoodSplit = callFMPart(problem);
    }

    if (foundGoodSplit) {
      if (_params.verb.getForActions() > 9)
        cout << "  Success.  Found a good partitioning" << endl;
    } else if (numTries < 2) {
      if (_params.verb.getForActions() > 9)
        cout << "  Failed. clearing part buffers" << endl;

      PartitioningBuffer& buff = problem.getSolnBuffers();
      unsigned numTerms = problem.getHGraph().getNumTerminals();

      for (unsigned sol = buff.beginUsedSoln(); sol != buff.endUsedSoln();
           ++sol) {
        Partitioning& part = buff[sol];
        for (unsigned i = numTerms; i < part.size(); ++i) part[i].setToAll(2);
      }

      // consider changing the problem here to attempt to
      // get a better solution
    } else {
      if (_params.verb.getForActions() > 9)
        cout << "  Failed too many times. Taking the best soln" << endl;
    }
  }
}

bool SplitSmallBinHorizontally::callBBFM(PartitioningProblemForCapo& problem) {
  double binWS = _bin.getRelativeWhitespace();

  if (_params.verb.getForActions() > 4) {
    cout << " BBFMPart" << endl;
    cout << " Targets: " << problem.getCapacities()[0][0] << "/"
         << problem.getCapacities()[1][0] << endl;
    cout << " Max    : " << problem.getMaxCapacities()[0][0] << "/"
         << problem.getMaxCapacities()[1][0] << endl;
    cout << " Min    : " << problem.getMinCapacities()[0][0] << "/"
         << problem.getMinCapacities()[1][0] << endl;
    cout << " Tolerance : " << 100 * problem.getTolerance(0) << endl;
    cout << endl;
  }

  double origCut = partitionByGroups(problem);
  bool partitionedByGroups = (origCut != -1.);

  if (!partitionedByGroups) {
    Timer smallPartTimer;

    problem.reserveBuffers(2);

    BBFMPart::Parameters bbfmparams(_capo.getParams().bbfmParams);
    bbfmparams.fmStarts = 5;
    bbfmparams.bbParams.pushLimit = 1000000;  //~2 seconds
    bbfmparams.verb = _params.verb;
    bbfmparams.verb.setForActions(bbfmparams.verb.getForActions() / 4);
    bbfmparams.verb.setForMajStats(bbfmparams.verb.getForMajStats() / 4);
    bbfmparams.verb.setForSysRes(bbfmparams.verb.getForSysRes() / 4);

    problem.getSolnBuffers().setBeginUsedSoln(0);
    problem.getSolnBuffers().setEndUsedSoln(1);

    // problem.saveAsNodesNets((string("problem")+string(_bin.getName())).c_str());

    BBFMPart partitioner1(problem, _capo.bbBitBoards, _capo.getParams().maxMem,
                          bbfmparams);
    smallPartTimer.stop();
    CapoPlacer::capoTimer::BBTime += smallPartTimer.getUserTime();

    unsigned origSplitRow = problem.getSplitRow();

    // otherwise, try rounding the split the other way..
    double origP0Capacity = problem.getP0Capacity();
    if ((_bin.getNumRows() % 2 != 0)) {
      if (problem.getSplitRow() <=
          floor(0.5 * static_cast<double>(_bin.getNumRows()))) {
        if (problem.getSplitRow() < (_bin.getNumRows() - 1)) {
          problem.setP0Capacity(
              problem.getP0Capacity() -
              _bin.getContainedAreaInCoreRow(problem.getSplitRow()));
          problem.setSplitRow(problem.getSplitRow() + 1);
        }
      } else {
        if (problem.getSplitRow() > 1) {
          problem.setP0Capacity(
              problem.getP0Capacity() +
              _bin.getContainedAreaInCoreRow(problem.getSplitRow() - 1));
          problem.setSplitRow(problem.getSplitRow() - 1);
        }
      }
    }

    abkfatal(
        problem.getSplitRow() != 0 && problem.getSplitRow() < _bin.getNumRows(),
        "new bin has no rows");

    if (origSplitRow != problem.getSplitRow()) {
      Timer tm16;

      double origP0Target = problem.getCapacities()[0][0];
      double origP1Target = problem.getCapacities()[1][0];
      double origP0MaxTarget = problem.getMaxCapacities()[0][0];
      double origP1MaxTarget = problem.getMaxCapacities()[1][0];
      double origP0MinTarget = problem.getMinCapacities()[0][0];
      double origP1MinTarget = problem.getMinCapacities()[1][0];
      problem.calculateCapacities_SmallH();
      problem.calculateTolerance_SmallH(
          false);  // false because this is a small instance
      problem.getSolnBuffers().setBeginUsedSoln(1);
      problem.getSolnBuffers().setEndUsedSoln(2);

      BBFMPart partitioner2(problem, _capo.bbBitBoards,
                            _capo.getParams().maxMem, bbfmparams);

      problem.getSolnBuffers().setBeginUsedSoln(0);

      // note to self...could include area-balances
      // and pin-balances in this decision as well
      if (problem.getCosts()[0] <= problem.getCosts()[1]) {
        problem.setBestSolnNum(0);
        problem.setP0Capacity(origP0Capacity);
        problem.calculateCapacities_SmallH();
        problem.setSplitRow(origSplitRow);
        problem.setCapacity(0, origP0Target);
        problem.setCapacity(1, origP1Target);
        problem.setMaxCapacity(0, origP0MaxTarget);
        problem.setMaxCapacity(1, origP1MaxTarget);
        problem.setMinCapacity(0, origP0MinTarget);
        problem.setMinCapacity(1, origP1MinTarget);
      } else {
        problem.setBestSolnNum(1);
      }

      tm16.stop();
      CapoPlacer::capoTimer::BBTime += tm16.getUserTime();
    }
  }

  vector<double> partAreas(2, 0.), actualCaps(2, 0.);
  const Partitioning& soln = problem.getBestSoln();
  const HGraphFixed& hgraph = problem.getHGraph();

  for (unsigned c = hgraph.getNumTerminals(); c < hgraph.getNumNodes(); ++c) {
    unsigned part = soln[c].isInPart(0) ? 0 : 1;
    partAreas[part] += hgraph.getWeight(c);
  }

  unsigned oldSplitRow = problem.getSplitRow();
  if (_params.WSOracle) {
    problem.setSplitRow(
        shiftSplitRow(problem.getSplitRow(), partAreas, actualCaps, 0.));
  } else if (greaterOrEqualDouble(_params.safeWS, 0.) &&
             binWS >= _params.safeWS) {
    // make sure each side of the split has as least this much WS
    problem.setSplitRow(shiftSplitRow(problem.getSplitRow(), partAreas,
                                      actualCaps, _params.safeWS));
  } else if (greaterOrEqualDouble(_params.minLocalWS, 0.) &&
             binWS >= _params.minLocalWS) {
    // make sure each side of the split has as least this much WS
    problem.setSplitRow(shiftSplitRow(problem.getSplitRow(), partAreas,
                                      actualCaps, _params.minLocalWS));
  } else  // uniform WS
  {
    double p0WSWeight = 1.;
    double p1WSWeight = 1.;
    problem.setSplitRow(
        findBestSplitRow(partAreas, p0WSWeight, p1WSWeight, actualCaps));
  }
  bool splitRowChanged = problem.getSplitRow() != oldSplitRow;

  if (!partitionedByGroups && (_params.repartShiftCL && splitRowChanged) ||
      problem.getRepartSmallWS()) {
    if (_params.verb.getForMajStats() > 1)
      cout << "Doing Repartitioning phase" << endl;

    Timer repartProb;
    const vector<const RBPCoreRow*>& rows = _bin.getRows();
    bool honorGrpConstraints = _capo.getParams().useGrpConstr;
    bool splitPtFixed = false;
    unsigned splitRow = problem.getSplitRow();
    PartitioningProblemForCapo problemNew(
        _bin, _capo, _capo.partProbEdgesVisited, _capo.thetoWeights, *this,
        true, Small, true, rows[splitRow]->getYCoord(), splitPtFixed, splitRow,
        actualCaps[0], Verbosity("1 1 1"), honorGrpConstraints);

    PartitioningBuffer& bufOld = problem.getSolnBuffers();
    unsigned bestSolnNumOld = problem.getBestSolnNum();
    PartitioningBuffer* pNewBuf = problemNew.swapOutSolnBuffers(&bufOld);
    Partitioning& pOldFC = const_cast<Partitioning&>(problem.getFixedConstr());
    Partitioning* pNewFC = problemNew.swapOutFixedConst(&pOldFC);
    problemNew.setBestSolnNum(bestSolnNumOld);

    unsigned bestSolnNum = problemNew.getBestSolnNum();
    problemNew.getSolnBuffers().setBeginUsedSoln(bestSolnNum);
    problemNew.getSolnBuffers().setEndUsedSoln(bestSolnNum + 1);

    repartProb.stop();
    CapoPlacer::capoTimer::SetupTime += repartProb.getUserTime();

    if (_params.verb.getForActions() > 4) {
      cout << " Second (smaller tol) FMPart" << endl;
      cout << " Targets: " << problemNew.getCapacities()[0][0] << "/"
           << problemNew.getCapacities()[1][0] << endl;
      cout << " Max    : " << problemNew.getMaxCapacities()[0][0] << "/"
           << problemNew.getMaxCapacities()[1][0] << endl;
      cout << " Min    : " << problemNew.getMinCapacities()[0][0] << "/"
           << problemNew.getMinCapacities()[1][0] << endl;
      cout << " Tolerance : " << 100 * problemNew.getTolerance(0) << endl;
      cout << endl;
    }

    Timer tm4;

    FMPartitioner::Parameters fmParams(_capo.getParams().mlParams);
    fmParams.useEarlyStop = false;
    fmParams.verb = _params.verb;
    fmParams.verb.setForActions(fmParams.verb.getForActions() / 4);
    fmParams.verb.setForMajStats(fmParams.verb.getForMajStats() / 4);
    fmParams.verb.setForSysRes(fmParams.verb.getForSysRes() / 4);
    fmParams.useClip = false;
    fmParams.relaxedTolerancePasses = 0;

    if (_params.useFMPartPlus) {
      FMPartitionerPlus partitioner2(problemNew, _capo.getParams().maxMem,
                                     fmParams, true);
    } else {
      FMPartitioner partitioner2(problemNew, _capo.getParams().maxMem, fmParams,
                                 true);
    }

    unsigned bestSolnNumNew = problemNew.getBestSolnNum();
    bool retValue = true;
    if (problemNew.getViolations()[bestSolnNumNew]) {
      if (_params.verb.getForActions() > 9)
        cout << " ->FM Failed.  No legal solutions" << endl;
      // retValue = false;
    }
    if (_params.verb.getForMajStats() > 4)
      cout << " Best Solution had cut " << problemNew.getCosts()[bestSolnNumNew]
           << endl;

    problemNew.swapOutSolnBuffers(pNewBuf);
    problemNew.swapOutFixedConst(pNewFC);
    problem.setBestSolnNum(bestSolnNumNew);
    problem.setSplitRow(problemNew.getSplitRow());

    partAreas[0] = 0;
    partAreas[1] = 0;
    const Partitioning& soln = problem.getBestSoln();
    const HGraphFixed& hgraph = problem.getHGraph();

    for (unsigned c = hgraph.getNumTerminals(); c < hgraph.getNumNodes(); ++c) {
      unsigned part = soln[c].isInPart(0) ? 0 : 1;
      partAreas[part] += hgraph.getWeight(c);
    }

    if (_params.WSOracle) {
      problem.setSplitRow(
          shiftSplitRow(oldSplitRow, partAreas, actualCaps, 0.));
    } else if (greaterOrEqualDouble(_params.safeWS, 0.) &&
               binWS >= _params.safeWS) {
      // make sure each side of the split has as least this much WS
      problem.setSplitRow(
          shiftSplitRow(oldSplitRow, partAreas, actualCaps, _params.safeWS));
    } else if (greaterOrEqualDouble(_params.minLocalWS, 0.) &&
               binWS >= _params.minLocalWS) {
      // make sure each side of the split has as least this much WS
      problem.setSplitRow(shiftSplitRow(oldSplitRow, partAreas, actualCaps,
                                        _params.minLocalWS));
    } else  // uniform WS
    {
      double p0WSWeight = 1.;
      double p1WSWeight = 1.;
      problem.setSplitRow(
          findBestSplitRow(partAreas, p0WSWeight, p1WSWeight, actualCaps));
    }

    tm4.stop();
    CapoPlacer::capoTimer::FMTime += tm4.getUserTime();
    return retValue;
  }

  unsigned bestSolnNum = problem.getBestSolnNum();
  if (_params.verb.getForMajStats() > 4) {
    cout << " Best Solution had cut " << problem.getCosts()[bestSolnNum]
         << endl;
  }

  if (problem.getViolations()[bestSolnNum]) {
    if (_params.verb.getForActions() > 9)
      cout << " ->FM Failed.  No legal solutions" << endl;
    return false;
  }

  return true;
}

bool SplitSmallBinHorizontally::callFMPart(
    PartitioningProblemForCapo& problem) {
  double binWS = _bin.getRelativeWhitespace();

  if (_params.verb.getForActions() > 4) {
    cout << " First FMPart" << endl;
    cout << " Targets: " << problem.getCapacities()[0][0] << "/"
         << problem.getCapacities()[1][0] << endl;
    cout << " Max    : " << problem.getMaxCapacities()[0][0] << "/"
         << problem.getMaxCapacities()[1][0] << endl;
    cout << " Min    : " << problem.getMinCapacities()[0][0] << "/"
         << problem.getMinCapacities()[1][0] << endl;
    cout << " Tolerance : " << 100 * problem.getTolerance(0) << endl;
    cout << endl;
  }

  FMPartitioner::Parameters fmParams(_capo.getParams().mlParams);

  double origCut = partitionByGroups(problem);
  bool partitionedByGroups = (origCut != -1.);

  if (!partitionedByGroups) {
    fmParams.useEarlyStop = false;
    fmParams.useClip = (_bin.getNumCells() > 75) ? 1 : 0;
    fmParams.verb = _params.verb;
    fmParams.verb.setForActions(fmParams.verb.getForActions() / 4);
    fmParams.verb.setForMajStats(fmParams.verb.getForMajStats() / 4);
    fmParams.verb.setForSysRes(fmParams.verb.getForSysRes() / 4);
    fmParams.relaxedTolerancePasses =
        _capo.getParams().mlParams.relaxedTolerancePasses;

    unsigned numFMStarts = 4;
    problem.reserveBuffers(numFMStarts);
    problem.getSolnBuffers().setBeginUsedSoln(0);
    problem.getSolnBuffers().setEndUsedSoln(numFMStarts);

    Timer tm3;
    std::unique_ptr<MultiStartPartitioner> partitioner1;
    if (_params.useFMPartPlus) {
      partitioner1.reset(
          new FMPartitionerPlus(problem, _capo.getParams().maxMem, fmParams));
    } else {
      partitioner1.reset(
          new FMPartitioner(problem, _capo.getParams().maxMem, fmParams));
    }
    tm3.stop();
    CapoPlacer::capoTimer::FMTime += tm3.getUserTime();
    origCut = partitioner1->getBestSolnCost();
  }

  vector<double> partAreas(2);
  const Partitioning& soln = problem.getBestSoln();
  const HGraphFixed& hgraph = problem.getHGraph();
  double finalCut = DBL_MAX;

  for (unsigned c = hgraph.getNumTerminals(); c < hgraph.getNumNodes(); c++) {
    unsigned part = soln[c].isInPart(0) ? 0 : 1;
    partAreas[part] += hgraph.getWeight(c);
  }

  vector<double> actualCaps(2, 0);
  unsigned oldSplitRow = problem.getSplitRow();
  if (_params.WSOracle) {
    problem.setSplitRow(
        shiftSplitRow(problem.getSplitRow(), partAreas, actualCaps, 0.));
  } else if (greaterOrEqualDouble(_params.safeWS, 0.) &&
             binWS >= _params.safeWS) {
    // make sure each side of the split has as least this much WS
    problem.setSplitRow(shiftSplitRow(problem.getSplitRow(), partAreas,
                                      actualCaps, _params.safeWS));
  } else if (greaterOrEqualDouble(_params.minLocalWS, 0.) &&
             binWS >= _params.minLocalWS) {
    // make sure each side of the split has as least this much WS
    problem.setSplitRow(shiftSplitRow(problem.getSplitRow(), partAreas,
                                      actualCaps, _params.minLocalWS));
  } else  // uniform WS
  {
    double p0WSWeight = 1.;
    double p1WSWeight = 1.;
    problem.setSplitRow(
        findBestSplitRow(partAreas, p0WSWeight, p1WSWeight, actualCaps));
  }

  if (!partitionedByGroups &&
      ((_params.repartShiftCL && problem.getSplitRow() != oldSplitRow) ||
       problem.getRepartSmallWS() || problem.getChangedTol())) {
    if (_params.verb.getForMajStats() > 1)
      cout << "Doing Repartitioning phase" << endl;

    Timer repartProb;
    const vector<const RBPCoreRow*>& rows = _bin.getRows();
    bool splitPtFixed = false;
    unsigned splitRow = problem.getSplitRow();
    Verbosity verb = Verbosity("1 1 1");
    bool honorGrpConstraints = false;
    PartitioningProblemForCapo problemNew(
        _bin, _capo, _capo.partProbEdgesVisited, _capo.thetoWeights, *this,
        true, Medium, true, rows[splitRow]->getYCoord(), splitPtFixed, splitRow,
        actualCaps[0], verb, honorGrpConstraints);

    PartitioningBuffer& bufOld = problem.getSolnBuffers();
    unsigned bestSolnNum = problem.getBestSolnNum();
    PartitioningBuffer* pNewBuf = problemNew.swapOutSolnBuffers(&bufOld);
    problemNew.setBestSolnNum(bestSolnNum);
    problemNew.getSolnBuffers().setBeginUsedSoln(bestSolnNum);
    problemNew.getSolnBuffers().setEndUsedSoln(bestSolnNum + 1);
    repartProb.stop();
    CapoPlacer::capoTimer::SetupTime += repartProb.getUserTime();

    fmParams.useEarlyStop = false;
    fmParams.useClip = false;
    fmParams.verb.setForActions(_params.verb.getForActions() / 4);
    fmParams.verb.setForMajStats(_params.verb.getForMajStats() / 4);
    fmParams.verb.setForSysRes(_params.verb.getForSysRes() / 4);
    fmParams.relaxedTolerancePasses = 0;

    PartitioningBuffer bestSolns(hgraph.getNumNodes(), 2);
    double bestCuts[2];
    unsigned bestSplits[2];
    bool satisfiesWS[2];

    bestSolns[0] = problem.getBestSoln();
    bestSolns[1] = problem.getBestSoln();
    bestCuts[0] = origCut;
    bestCuts[1] = origCut;
    bestSplits[0] = problemNew.getSplitRow();
    bestSplits[1] = problemNew.getSplitRow();

    double p0WS = (actualCaps[0] - partAreas[0]) / actualCaps[0];
    double p1WS = (actualCaps[1] - partAreas[1]) / actualCaps[1];
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

    if (_params.verb.getForActions() > 4)
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
        if (_params.verb.getForActions() > 4) {
          cout << " Tightening pass" << endl;
          cout << " Targets: " << problemNew.getCapacities()[0][0] << "/"
               << problemNew.getCapacities()[1][0] << endl;
          cout << " Max:     " << problemNew.getMaxCapacities()[0][0] << "/"
               << problemNew.getMaxCapacities()[1][0] << endl;
          cout << " Min:     " << problemNew.getMinCapacities()[0][0] << "/"
               << problemNew.getMinCapacities()[1][0] << endl;
          cout << " Tolerance : " << 100 * problemNew.getTolerance(0) << endl;
        }

        Timer tm4;
        std::unique_ptr<MultiStartPartitioner> partitioner2;
        if (_params.useFMPartPlus) {
          partitioner2.reset(new FMPartitionerPlus(
              problemNew, _capo.getParams().maxMem, fmParams, true));
        } else {
          partitioner2.reset(new FMPartitioner(
              problemNew, _capo.getParams().maxMem, fmParams, true));
        }
        tm4.stop();
        CapoPlacer::capoTimer::FMTime += tm4.getUserTime();

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
            bestCuts[1] = newCut;
            bestSolns[1] = problemNew.getBestSoln();
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
              bestCuts[1] = newCut;
              bestSolns[1] = problemNew.getBestSoln();
              break;
            }
          }
        } else {
          if (newCut <= 1.04 * origCut) {
            bestCuts[1] = newCut;
            bestSolns[1] = problemNew.getBestSoln();
            break;
          } else {
            problemNew.getSolnBuffers()[bestSolnNum] = bestSolns[1];
            problemNew.setBestSolnNum(bestSolnNum);
            problemNew.setMaxCapacity(0, lastP0Max);
            problemNew.setMinCapacity(0, lastP0Min);
            problemNew.setMaxCapacity(1, lastP1Max);
            problemNew.setMinCapacity(1, lastP1Min);
            break;
          }
        }
        bestCuts[1] = newCut;
        bestSolns[1] = problemNew.getBestSoln();
        firstPass = false;
      }
    }

    partAreas[0] = 0.;
    partAreas[1] = 0.;
    for (unsigned c = 0; c < problemNew.getHGraph().getNumNodes(); ++c) {
      unsigned part = bestSolns[1][c].isInPart(0) ? 0 : 1;
      partAreas[part] += problemNew.getHGraph().getWeight(c);
    }

    if (_params.WSOracle) {
      bestSplits[1] = shiftSplitRow(bestSplits[1], partAreas, actualCaps, 0.);
    } else if (greaterOrEqualDouble(_params.safeWS, 0.) &&
               binWS >= _params.safeWS) {
      // make sure each side of the split has as least this much WS
      bestSplits[1] =
          shiftSplitRow(bestSplits[1], partAreas, actualCaps, _params.safeWS);
    } else if (greaterOrEqualDouble(_params.minLocalWS, 0.) &&
               binWS >= _params.minLocalWS) {
      // make sure each side of the split has as least this much WS
      bestSplits[1] = shiftSplitRow(bestSplits[1], partAreas, actualCaps,
                                    _params.minLocalWS);
    } else  // uniform WS
    {
      double p0WSWeight = 1.;
      double p1WSWeight = 1.;
      bestSplits[1] =
          findBestSplitRow(partAreas, p0WSWeight, p1WSWeight, actualCaps);
    }

    // check to see if this satisfies the min local constraints
    p0WS = (actualCaps[0] - partAreas[0]) / actualCaps[0];
    p1WS = (actualCaps[1] - partAreas[1]) / actualCaps[1];
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
    problemNew.setSplitRow(bestSplits[0]);
    finalCut = bestCuts[0];
    }
    else
    {
    problemNew.setSplitRow(bestSplits[1]);
    finalCut = bestCuts[1];
    }
    }
    else
    {
    abkwarn(0,"Original solution satisfied WS, but tightened did not. Probably a
    bug."); problemNew.getSolnBuffers()[0] = bestSolns[0];
    problemNew.setBestSolnNum(0);
    problemNew.setSplitRow(bestSplits[0]);
    finalCut = bestCuts[0];
    }
    }
    else
    {
    if(satisfiesWS[1])
    {
    problemNew.setSplitRow(bestSplits[1]);
    finalCut = bestCuts[1];
    }
    else
    {
    if(bestCuts[0] < bestCuts[1])
    {
    abkwarn(0,"Dropping tightened solution because both bad wrt local ws and
    orig had better cut"); problemNew.getSolnBuffers()[0] = bestSolns[0];
    problemNew.setBestSolnNum(0);
    _splitRow = bestSplits[0];
    finalCut = bestCuts[0];
    }
    else
    {
    problemNew.setSplitRow(bestSplits[1]);
    finalCut = bestCuts[1];
    }
    }
    }
     */

    problemNew.setSplitRow(bestSplits[1]);
    finalCut = bestCuts[1];

    unsigned bestSolnNumNew = problemNew.getBestSolnNum();
    problemNew.swapOutSolnBuffers(pNewBuf);
    problem.setBestSolnNum(bestSolnNumNew);
    problem.setSplitRow(problemNew.getSplitRow());

    if (problem.getViolations()[bestSolnNumNew]) {
      if (_params.verb.getForActions() > 9)
        cout << " ->FM Failed.  No legal solutions" << endl;
    }

    if (_params.verb.getForMajStats() > 4)
      cout << " Best Solution had cut " << finalCut << endl;

    return true;
  }

  unsigned bestSolnNum = problem.getBestSolnNum();
  if (_params.verb.getForMajStats() > 4)
    cout << " Best Solution had cut " << problem.getCosts()[bestSolnNum]
         << endl;

  if (problem.getViolations()[bestSolnNum]) {
    if (_params.verb.getForActions() > 9)
      cout << " ->FM Failed.  No legal solutions" << endl;
    return false;
  }

  return true;
}
