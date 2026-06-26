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

SplitSmallBinVertically::SplitSmallBinVertically(CapoBin& binToSplit,
                                                 CapoPlacer& capo,
                                                 CapoSplitterParams params)
    : BaseBinSplitter(binToSplit, capo, params) {
  if (_params.verb.getForActions() > 2)
    cout << "Splitting a small bin vertically" << endl;

  // 1)construct the partitioning problem
  // 2)while(!foundAGoodSplit)
  //     decide which partitioner to use, and partition the bin
  //     if needed, change tolerances and repartition
  //	  see if this is a good split (legal, etc)
  // 3)construct new bins

  Timer tm25;

  double capacity = _bin.getCapacity();
  double targetCap0, targetCap1;
  targetCap0 = 0.5 * capacity;
  targetCap1 = targetCap0;

  vector<double> partCaps;

  double p0WSWeight = 1.;
  double p1WSWeight = 1.;
  double xSplit = _bin.findXSplit(targetCap0, targetCap1, 0, p0WSWeight,
                                  p1WSWeight, partCaps);

  bool honorGrpConstraints = _capo.getParams().useGrpConstr;

  PartProbSizeType Size =
      (_bin.getNumCells() > _capo.getParams().smallPartThreshold) ? Medium
                                                                  : Small;
  bool splitPtFixed = false;
  unsigned splitRow = UINT_MAX;
  PartitioningProblemForCapo problem(
      _bin, _capo, _capo.partProbEdgesVisited, _capo.thetoWeights, *this, false,
      Size, false, xSplit, splitPtFixed, splitRow, partCaps[0],
      Verbosity("1 1 1"), honorGrpConstraints);
  tm25.stop();
  CapoPlacer::capoTimer::SetupTime += tm25.getUserTime();

  if (_params.eco && _bin.isEcoAllowed() &&
      (_params.ecoOverfillPct > _bin.getRelativeWhitespace())) {
    _bin.setEcoAllowed(false);
  }

  bool shouldPartByGroups = shouldPartitionByGroups(problem);
  if (shouldPartByGroups) _bin.setEcoAllowed(false);

  if (_params.eco && _bin.isEcoAllowed()) {
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
        problem.setXSplit(xSplit);

        ECOtimer.stop();
        CapoPlacer::capoTimer::ECOTime += ECOtimer.getUserTime();
        ++CapoPlacer::EcoStats::ecoPartsRejected;

        callPartitioner(problem);
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
    callPartitioner(problem);
    // <aaronnn> adjust cut of partitioning solution to resize bins by
    // utilization
    bool isHCut = false;
    if (_params.adjustWSByDensity) {
      problem.setXSplit(adjustSplitByUtilizationDensity(
          problem, problem.getXSplit(), isHCut));
    }
  }

  if (_params.verb.getForActions() > 8)
    cout << "Done partitioning...creating bins" << endl;

  bool isLarge = false;
  createVSplitBins(isLarge, problem, problem.getXSplit());
}

void SplitSmallBinVertically::callPartitioner(
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

bool SplitSmallBinVertically::callBBFM(PartitioningProblemForCapo& problem) {
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
    BBFMPart::Parameters bbfmparams(_capo.getParams().bbfmParams);
    bbfmparams.fmStarts = 5;
    bbfmparams.bbParams.pushLimit = 1000000;  //~2 seconds
    bbfmparams.verb = _params.verb;
    bbfmparams.verb.setForActions(bbfmparams.verb.getForActions() / 4);
    bbfmparams.verb.setForMajStats(bbfmparams.verb.getForMajStats() / 4);
    bbfmparams.verb.setForSysRes(bbfmparams.verb.getForSysRes() / 4);

    problem.reserveBuffers(1);
    problem.getSolnBuffers().setBeginUsedSoln(0);
    problem.getSolnBuffers().setEndUsedSoln(1);

    Timer tm17;
    BBFMPart partitioner1(problem, _capo.bbBitBoards, _capo.getParams().maxMem,
                          bbfmparams);
    tm17.stop();
    CapoPlacer::capoTimer::BBTime += tm17.getUserTime();
  }

  double oldXSplit = problem.getXSplit();
  vector<double> partAreas(2, 0.);
  vector<double> siteAreas(2, 0.);
  const Partitioning& soln = problem.getBestSoln();
  for (unsigned c = 0; c < problem.getHGraph().getNumNodes(); c++) {
    unsigned part = soln[c].isInPart(0) ? 0 : 1;
    partAreas[part] += problem.getHGraph().getWeight(c);
  }

  if (_params.WSOracle) {
    problem.setXSplit(
        _bin.shiftXSplit(oldXSplit, 0., partAreas[0], partAreas[1], siteAreas));
  } else if (greaterOrEqualDouble(_params.safeWS, 0.) &&
             binWS >= _params.safeWS) {
    // make sure each side of the split has as least this much WS
    problem.setXSplit(_bin.shiftXSplit(oldXSplit, _params.safeWS, partAreas[0],
                                       partAreas[1], siteAreas));
  } else if (greaterOrEqualDouble(_params.minLocalWS, 0.) &&
             binWS >= _params.minLocalWS) {
    // make sure each side of the split has as least this much WS
    problem.setXSplit(_bin.shiftXSplit(oldXSplit, _params.minLocalWS,
                                       partAreas[0], partAreas[1], siteAreas));
  } else  // uniform WS
  {
    double p0WSWeight = 1.;
    double p1WSWeight = 1.;
    problem.setXSplit(_bin.findXSplit(partAreas[0], partAreas[1], 10,
                                      p0WSWeight, p1WSWeight, siteAreas));
  }

  // really ought to be more pro-active here.......
  //     if(binWS > 0 && ((siteAreas[0] < partAreas[0] || siteAreas[1] <
  //     partAreas[1])))
  //     {
  //	problem.setXSplit(_bin.findXSplit(partAreas[0], partAreas[1], 0,
  //siteAreas));
  //     }

  // cutline was shifted
  if (!partitionedByGroups && _params.repartShiftCL &&
      problem.getXSplit() != oldXSplit) {
    if (_params.verb.getForMajStats() > 1)
      cout << "Doing Repartitioning phase" << endl;

    Timer repartProb;

    bool honorGrpConstraints = _capo.getParams().useGrpConstr;
    bool splitPtFixed = false;
    unsigned splitRow = UINT_MAX;
    PartitioningProblemForCapo problemNew(
        _bin, _capo, _capo.partProbEdgesVisited, _capo.thetoWeights, *this,
        false, Small, true, problem.getXSplit(), splitPtFixed, splitRow,
        siteAreas[0], Verbosity("1 1 1"), honorGrpConstraints);

    PartitioningBuffer& bufOld = problem.getSolnBuffers();
    unsigned bestSolnNumOld = problem.getBestSolnNum();

    PartitioningBuffer* pNewBuf = problemNew.swapOutSolnBuffers(&bufOld);
    Partitioning& pOldFC = const_cast<Partitioning&>(problem.getFixedConstr());
    Partitioning* pNewFC = problemNew.swapOutFixedConst(&pOldFC);
    problemNew.setBestSolnNum(bestSolnNumOld);

    repartProb.stop();
    CapoPlacer::capoTimer::SetupTime += repartProb.getUserTime();

    if (_params.verb.getForMajStats() > 4) {
      cout << " Repartitioning FMPart" << endl;
      cout << "Targets:  " << problemNew.getCapacities()[0][0] << "/"
           << problemNew.getCapacities()[1][0] << endl;
      cout << "Max:      " << problemNew.getMaxCapacities()[0][0] << "/"
           << problemNew.getMaxCapacities()[1][0] << endl;
      cout << "Min:      " << problemNew.getMinCapacities()[0][0] << "/"
           << problemNew.getMinCapacities()[1][0] << endl;
    }

    unsigned bestSolnNum = problemNew.getBestSolnNum();
    problemNew.getSolnBuffers().setBeginUsedSoln(bestSolnNum);
    problemNew.getSolnBuffers().setEndUsedSoln(bestSolnNum + 1);

    FMPartitioner::Parameters fmParams(_capo.getParams().mlParams);
    fmParams.useEarlyStop = false;
    fmParams.verb = _params.verb;
    fmParams.verb.setForActions(fmParams.verb.getForActions() / 4);
    fmParams.verb.setForMajStats(fmParams.verb.getForMajStats() / 4);
    fmParams.verb.setForSysRes(fmParams.verb.getForSysRes() / 4);
    fmParams.useClip = false;
    fmParams.relaxedTolerancePasses = 0;

    Timer tm2;
    if (_params.useFMPartPlus) {
      FMPartitionerPlus partitioner2(problemNew, _capo.getParams().maxMem,
                                     fmParams, true);
    } else {
      FMPartitioner partitioner2(problemNew, _capo.getParams().maxMem, fmParams,
                                 true);
    }
    tm2.stop();
    CapoPlacer::capoTimer::FMTime += tm2.getUserTime();

    unsigned bestSolnNumNew = problemNew.getBestSolnNum();
    if (_params.verb.getForMajStats() > 4)
      cout << " Best Solution had cut " << problemNew.getCosts()[bestSolnNumNew]
           << endl;

    problemNew.swapOutSolnBuffers(pNewBuf);
    problemNew.swapOutFixedConst(pNewFC);
    problem.setBestSolnNum(bestSolnNumNew);
    problem.setXSplit(problemNew.getXSplit());

    const Partitioning& solnNew = problem.getBestSoln();
    partAreas[0] = 0;
    partAreas[1] = 0;
    for (unsigned c = 0; c < problem.getHGraph().getNumNodes(); c++) {
      unsigned part = solnNew[c].isInPart(0) ? 0 : 1;
      partAreas[part] += problem.getHGraph().getWeight(c);
    }

    if (_params.WSOracle) {
      problem.setXSplit(_bin.shiftXSplit(problem.getXSplit(), 0., partAreas[0],
                                         partAreas[1], siteAreas));
    } else if (greaterOrEqualDouble(_params.safeWS, 0.) &&
               binWS >= _params.safeWS) {
      // make sure each side of the split has as least this much WS
      problem.setXSplit(_bin.shiftXSplit(problem.getXSplit(), _params.safeWS,
                                         partAreas[0], partAreas[1],
                                         siteAreas));
    } else if (greaterOrEqualDouble(_params.minLocalWS, 0.) &&
               binWS >= _params.minLocalWS) {
      // make sure each side of the split has as least this much WS
      problem.setXSplit(_bin.shiftXSplit(problem.getXSplit(),
                                         _params.minLocalWS, partAreas[0],
                                         partAreas[1], siteAreas));
    } else  // uniform WS
    {
      double p0WSWeight = 1.;
      double p1WSWeight = 1.;
      problem.setXSplit(_bin.findXSplit(partAreas[0], partAreas[1], 5,
                                        p0WSWeight, p1WSWeight, siteAreas));
    }

    return true;
  }

  unsigned bestSolnNum = problem.getBestSolnNum();
  if (_params.verb.getForMajStats() > 4)
    cout << " Best Solution had cut " << problem.getCosts()[bestSolnNum]
         << endl;

  return true;
}

bool SplitSmallBinVertically::callFMPart(PartitioningProblemForCapo& problem) {
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

  double binWS = _bin.getRelativeWhitespace();

  FMPartitioner::Parameters fmParams(_capo.getParams().mlParams);
  double origCut = partitionByGroups(problem);
  bool partitionedByGroups = (origCut != -1.);

  if (!partitionedByGroups) {
    fmParams.useEarlyStop = false;
    fmParams.useClip = (_bin.getNumCells() > 50) ? 1 : 0;
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

    Timer tm1;
    std::unique_ptr<MultiStartPartitioner> partitioner1;
    if (_params.useFMPartPlus) {
      partitioner1.reset(
          new FMPartitionerPlus(problem, _capo.getParams().maxMem, fmParams));
    } else {
      partitioner1.reset(
          new FMPartitioner(problem, _capo.getParams().maxMem, fmParams));
    }
    tm1.stop();
    CapoPlacer::capoTimer::FMTime += tm1.getUserTime();
    origCut = partitioner1->getBestSolnCost();
  }

  vector<double> partAreas(2, 0.);
  vector<double> siteAreas(2, 0.);
  const Partitioning& soln = problem.getBestSoln();
  const HGraphFixed& hgraph = problem.getHGraph();
  double finalCut = DBL_MAX;

  for (unsigned c = hgraph.getNumTerminals(); c < hgraph.getNumNodes(); c++) {
    unsigned part = soln[c].isInPart(0) ? 0 : 1;
    partAreas[part] += hgraph.getWeight(c);
  }

  // for small bins this isn't really the best
  // way to do this.  Really, we ought to compute the
  // range of splits that give us some min ws% in
  // each new bin, and then pick the best xSplit
  // within that range.
  double oldXSplit = problem.getXSplit();
  double p0WSWeight = 1.;
  double p1WSWeight = 1.;
  double idealXSplit = _bin.findXSplit(partAreas[0], partAreas[1], 0,
                                       p0WSWeight, p1WSWeight, siteAreas);

  if (_params.WSOracle) {
    problem.setXSplit(
        _bin.shiftXSplit(oldXSplit, 0., partAreas[0], partAreas[1], siteAreas));
  } else if (greaterOrEqualDouble(_params.safeWS, 0.) &&
             binWS >= _params.safeWS) {
    // make sure each side of the split has as least this much WS
    problem.setXSplit(_bin.shiftXSplit(oldXSplit, _params.safeWS, partAreas[0],
                                       partAreas[1], siteAreas));
  } else if (greaterOrEqualDouble(_params.minLocalWS, 0.) &&
             binWS >= _params.minLocalWS) {
    // make sure each side of the split has as least this much WS
    problem.setXSplit(_bin.shiftXSplit(oldXSplit, _params.minLocalWS,
                                       partAreas[0], partAreas[1], siteAreas));
  } else  // uniformWS
  {
    double p0WSWeight = 1.;
    double p1WSWeight = 1.;
    problem.setXSplit(_bin.findXSplit(partAreas[0], partAreas[1], 5, p0WSWeight,
                                      p1WSWeight, siteAreas));
  }

  if (_params.verb.getForMajStats() > 9) {
    cout << " CellAreas: " << partAreas[0] << "/" << partAreas[1] << endl;
    cout << " Total:     " << _bin.getTotalCellArea() << endl;
    cout << " SiteAreas: " << siteAreas[0] << "/" << siteAreas[1] << endl;
  }

  // cutline was shifted
  if (!partitionedByGroups &&
      (_params.doRepartitioning ||
       (_params.repartShiftCL && problem.getXSplit() != oldXSplit) ||
       problem.getRepartSmallWS())) {
    if (_params.verb.getForMajStats() > 1)
      cout << "Doing Repartitioning phase" << endl;

    Timer repartProb;

    bool honorGrpConstraints = _capo.getParams().useGrpConstr;
    bool splitPtFixed = false;
    unsigned splitRow = UINT_MAX;
    PartitioningProblemForCapo problemNew(
        _bin, _capo, _capo.partProbEdgesVisited, _capo.thetoWeights, *this,
        false, Medium, true, problem.getXSplit(), splitPtFixed, splitRow,
        siteAreas[0], Verbosity("1 1 1"), honorGrpConstraints);

    PartitioningBuffer& bufOld = problem.getSolnBuffers();
    unsigned bestSolnNum = problem.getBestSolnNum();
    PartitioningBuffer* pNewBuf = problemNew.swapOutSolnBuffers(&bufOld);
    problemNew.setBestSolnNum(bestSolnNum);
    problemNew.getSolnBuffers().setBeginUsedSoln(bestSolnNum);
    problemNew.getSolnBuffers().setEndUsedSoln(bestSolnNum + 1);

    repartProb.stop();
    CapoPlacer::capoTimer::SetupTime += repartProb.getUserTime();

    fmParams.relaxedTolerancePasses = 0;
    fmParams.useClip = false;

    PartitioningBuffer bestSolns(hgraph.getNumNodes(), 2);
    double bestCuts[2];
    double bestSplits[2];
    bool satisfiesWS[2];

    bestSolns[0] = problem.getBestSoln();
    bestSolns[1] = problem.getBestSoln();
    bestCuts[0] = origCut;
    bestCuts[1] = origCut;
    bestSplits[0] = problemNew.getXSplit();
    bestSplits[1] = problemNew.getXSplit();

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

        Timer tm2;
        std::unique_ptr<MultiStartPartitioner> partitioner2;
        if (_params.useFMPartPlus) {
          partitioner2.reset(new FMPartitionerPlus(
              problemNew, _capo.getParams().maxMem, fmParams, true));
        } else {
          partitioner2.reset(new FMPartitioner(
              problemNew, _capo.getParams().maxMem, fmParams, true));
        }
        tm2.stop();
        CapoPlacer::capoTimer::FMTime += tm2.getUserTime();

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
      bestSplits[1] = _bin.shiftXSplit(bestSplits[1], 0., partAreas[0],
                                       partAreas[1], siteAreas);
    } else if (greaterOrEqualDouble(_params.safeWS, 0.) &&
               binWS >= _params.safeWS) {
      // make sure each side of the split has as least this much WS
      bestSplits[1] = _bin.shiftXSplit(bestSplits[1], _params.safeWS,
                                       partAreas[0], partAreas[1], siteAreas);
    } else if (greaterOrEqualDouble(_params.minLocalWS, 0.) &&
               binWS >= _params.minLocalWS) {
      // make sure each side of the split has as least this much WS
      bestSplits[1] = _bin.shiftXSplit(bestSplits[1], _params.minLocalWS,
                                       partAreas[0], partAreas[1], siteAreas);
    } else  // uniform WS
    {
      double p0WSWeight = 1.;
      double p1WSWeight = 1.;
      bestSplits[1] = _bin.findXSplit(partAreas[0], partAreas[1], 0, p0WSWeight,
                                      p1WSWeight, siteAreas);
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
    abkwarn(0,"Dropping tightened solution because both bad wrt local ws and
    orig had better cut"); problemNew.getSolnBuffers()[0] = bestSolns[0];
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
    problem.setBestSolnNum(bestSolnNumNew);
    problem.setXSplit(problemNew.getXSplit());

    if (problem.getViolations()[bestSolnNumNew]) {
      if (_params.verb.getForActions() > 9)
        cout << " ->FM Failed.  No legal solutions" << endl;
    }

    if (_params.verb.getForMajStats() > 4)
      cout << " Best Solution had cut " << finalCut << endl;

    return true;
  } else {
    problem.setXSplit(idealXSplit);
  }

  unsigned bestSolnNum = problem.getBestSolnNum();
  if (_params.verb.getForMajStats() > 4)
    cout << " Best Solution had cut " << problem.getCosts()[bestSolnNum]
         << endl;

  if (problem.getViolations()[bestSolnNum]) {
    if (_params.verb.getForActions() > 9)
      cout << " ->FM Failed.  No legal solutions" << endl;
    problem.setXSplit(idealXSplit);
    return false;
  }

  return true;
}
