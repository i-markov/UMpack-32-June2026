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

// Created: April 1, 1998 by I.Markov by splitting A.Caldwell's capoPlacer.cxx

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

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
//        ilm split capoPlacer.cxx into this file, splitBin.cxx and the rest
//        ilm introduced "CapoPlacer::splitBin()"
// 990319 ilm refineBin() now returns bool "this bin has children that
//            need to be partitioned/placed"
// 001203 phm Fixed refine bin to use aspect ratio, in the same manner
//            as Feng Shui.  This should be a cmd line arg at some point
//            in the future, or calculated within the code.//001206 ilm The
//            above kicks in only when average row spacing is at least 1.5x the
//            row height

#include <climits>
#include <sstream>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkfunc.h"
#include "ABKCommon/abktimer.h"
#include "CongestionMaps/CongestionMaps.h"
#include "Placement/placement.h"
#include "RBPlace/RBPlacement.h"
#include "SmallPlacers/sRowSmPlacer.h"
#include "capoBin.h"
#include "capoPlacer.h"
#include "cutPosition.h"
#include "laSplitter.h"
#include "partProbForCapo.h"
#include "splitLarge.h"
#include "splitRow.h"
#include "splitSmall.h"

using std::cerr;
using std::cout;
using std::endl;
using std::max;
using std::ostringstream;
using std::pair;
using std::vector;

bool CapoPlacer::refineBin(CapoBin& bin, vector<CapoBin*>& newLayout,
                           AllowedPartDir partDir, unsigned minNumCellsToPart,
                           bool useCongestionInfo) {
  const unsigned numLeaves = bin.getNumCells();
  double bFactorAdjustCutoff = _params.splitterParams.bFactorAdjustCutoff;

  const HGraphWDimensions& hgraph = getNetlistHGraph();

  // if(_params.trackRentNums && numLeaves >= 10)
  //   storeRentNumbers(bin);
  bool isHCut = 0;

  if (numLeaves > 0 && bin.getRows().size() == 0) {
    bin._canSplitBin = false;
    _placedBins.push_back(&bin);

    if (_params.verb.getForActions() > 0)
      cout << " - moved non-empty bin [" << bin.getIndex()
           << "] with no core rows to placed bins" << endl;

    return false;
  } else if (numLeaves == 0) {
    bin._canSplitBin = false;
    _placedBins.push_back(&bin);

    if (_params.verb.getForActions() > 0)
      cout << " - moved empty bin [" << bin.getIndex() << "] to placed bins"
           << endl;

    return false;
  } else if (numLeaves == 1) {
    Timer weberTime;

    unsigned nodeId = *(bin.cellIdsBegin());
    const BBox& binBBox = bin.getBBox();
    Point oldPlace(_placement[nodeId]);
    _placement[nodeId] = Point(binBBox.xMin, (*bin.rowsBegin())->getYCoord());

    // find weber loc and place closest to this loc as possible
    Point weberLoc = getWeberLocation(nodeId);
    double nodeWidth = getRBP().getCoreRow(0).getSpacing();
    double nodeHeight = hgraph.getNodeHeight(nodeId);
    if (greaterThanDouble(weberLoc.x, (binBBox.xMax - nodeWidth)))
      weberLoc.x = binBBox.xMax - nodeWidth;
    if (lessThanDouble(weberLoc.x, binBBox.xMin)) weberLoc.x = binBBox.xMin;
    if (greaterThanDouble(weberLoc.y, (binBBox.yMax - nodeHeight)))
      weberLoc.y = binBBox.yMax - nodeHeight;
    if (lessThanDouble(weberLoc.y, binBBox.yMin)) weberLoc.y = binBBox.yMin;

    itCBCoreRow cr;
    itRBPSubRow sr;
    unsigned crIdx = 0;
    for (cr = bin.rowsBegin(); cr != (bin.rowsEnd() - 1); cr++) {
      double rowLowYCoord = (*cr)->getYCoord();
      double rowHighYCoord = (*(cr + 1))->getYCoord();
      if (weberLoc.y >= rowLowYCoord && weberLoc.y < rowHighYCoord) break;
      ++crIdx;
    }
    weberLoc.y = (*cr)->getYCoord();

    pair<unsigned, unsigned> stOffset = bin.getStartOffset(crIdx);
    pair<unsigned, unsigned> eOffset = bin.getEndOffset(crIdx);

    unsigned startOffset = stOffset.first;
    unsigned endOffset = eOffset.first;

    bool foundValidLoc = false;
    for (sr = (*cr)->subRowsBegin() + startOffset;
         sr != (*cr)->subRowsBegin() + endOffset + 1; ++sr) {
      if ((weberLoc.x >= sr->getXStart()) && (weberLoc.x <= sr->getXEnd())) {
        if (sr->getXEnd() < (weberLoc.x + nodeWidth))
          weberLoc.x = sr->getXEnd() - nodeWidth;
        unsigned numSitesFrmStart = static_cast<unsigned>(
            (weberLoc.x - sr->getXStart()) / sr->getSpacing());

        weberLoc.x = sr->getXStart() + numSitesFrmStart * sr->getSpacing();
        foundValidLoc = true;
      }
    }
    if (!foundValidLoc) {
      if (weberLoc.x <= (*cr)->getXStart())
        weberLoc.x = (*cr)->getXStart();
      else if (weberLoc.x >= (*cr)->getXEnd())
        weberLoc.x = (*cr)->getXEnd() - nodeWidth;
      else {
        for (sr = (*cr)->subRowsBegin() + startOffset;
             sr != (*cr)->subRowsBegin() + endOffset; ++sr) {
          if (sr->getXEnd() < weberLoc.x &&
              (sr + 1)->getXStart() > weberLoc.x) {
            if ((sr->getXEnd() - nodeWidth - weberLoc.x) <
                ((sr + 1)->getXStart() - weberLoc.x))
              weberLoc.x = sr->getXEnd() - nodeWidth;
            else
              weberLoc.x = (sr + 1)->getXStart();
            foundValidLoc = true;
          }
        }
      }
    }

    // final check

    if (greaterThanDouble(weberLoc.x, (binBBox.xMax - nodeWidth)))
      weberLoc.x = binBBox.xMax - nodeWidth;
    if (lessThanDouble(weberLoc.x, binBBox.xMin)) weberLoc.x = binBBox.xMin;
    if (greaterThanDouble(weberLoc.y, (binBBox.yMax - nodeHeight)))
      weberLoc.y = binBBox.yMax - nodeHeight;
    if (lessThanDouble(weberLoc.y, binBBox.yMin)) weberLoc.y = binBBox.yMin;

    weberTime.stop();
    CapoPlacer::capoTimer::SmTime += weberTime.getUserTime();

    _placed[nodeId] = true;
    if (_params.wtCut) {
      Timer psTimer;
      updatePointSetsForNewlyPlacedCell(nodeId, oldPlace, weberLoc);
      psTimer.stop();
      CapoPlacer::capoTimer::PointSetTime += psTimer.getUserTime();
    }

    _placement[nodeId] = weberLoc;
    _cellToBinMap[nodeId] = NULL;

    cerr << "1";
    bin._canSplitBin = false;
    _placedBins.push_back(&bin);

    return false;
  }
  /*sadya:take care of special case where xmin == xmax or 0 capacity bins*/
  else if ((bin.getBBox().xMin >= bin.getBBox().xMax &&
            bin.getNumRows() == 1) ||
           bin.getCapacity() == 0 || bin.getNumSites() == 1 ||
           bin.getNumSites() == 0) {
    Timer specialCaseTime;

    if (_params.verb.getForActions() > 1) {
      if (bin.getCapacity() == 0) {
        abkwarn(0, "skipping partitioning for 0 Capacity bin");
      } else if (bin.getNumSites() == 1) {
        abkwarn(0, "skipping partitioning for bin w 1 site");
      } else {
        abkwarn(0, "skipping partitioning for 0-xspan bin w/ 1 subrow");
      }
    }

    specialCaseTime.stop();
    CapoPlacer::capoTimer::SmTime += specialCaseTime.getUserTime();

    const vector<unsigned>& cellIds = bin.getCellIds();
    for (unsigned i = 0; i < cellIds.size(); i++) {
      _placed[cellIds[i]] = true;
      _cellToBinMap[cellIds[i]] = NULL;
      Point newPos(bin.getBBox().xMin, bin.getBBox().yMin);

      if (_params.wtCut) {
        Timer psTimer;
        updatePointSetsForNewlyPlacedCell(cellIds[i], _placement[cellIds[i]],
                                          newPos);
        psTimer.stop();
        CapoPlacer::capoTimer::PointSetTime += psTimer.getUserTime();
      }

      _placement[cellIds[i]] = newPos;
    }
    cerr << "-";

    bin._canSplitBin = false;
    _placedBins.push_back(&bin);

    return false;
  } else if (bin.getNumRows() == 1 && numLeaves < _params.smallPlaceThreshold &&
             bin._startOffsets[0].first == bin._endOffsets[0].first) {
    Timer smallPlaceProbTime;

    if (_params.verb.getForActions() > 1) {
      abkwarn(bin._startOffsets[0].first == bin._endOffsets[0].first,
              "smallPlacer called on a bin w/ subrows");
    }

    const RBPCoreRow& onlyRow = **bin.rowsBegin();
    double spacing = onlyRow.getSpacing();
    double yCoord = onlyRow.getYCoord();

    SmallPlacementRow splRow(bin.getBBox(), spacing);

    // added by sadya to take care of white space by fake cells
    vector<unsigned> movables;
    vector<double> whiteSpaceWidths;
    vector<unsigned> cellIds;
    // [MAR] Note: original code was buggy, didn't handle sparse cellrows in
    // which site width != spacing.
    int totalCellSites =
        int(floor(bin.getTotalCellArea() / onlyRow.getSiteArea()));
    int totalWSSites = bin.getNumSites() - totalCellSites;
    double totalWSWidth = (double)totalWSSites * spacing;
    totalWSWidth =
        max(totalWSWidth,
            0.0);  // don't allow negative totalWSWidth in overfull bins

    unsigned numSitesWS = 0;  // default is 0

    if (totalWSWidth > 0) numSitesWS = unsigned(floor(totalWSWidth / spacing));

    unsigned sitesPerWS = 0;

    if (numSitesWS >= 3) {
      movables.push_back(UINT_MAX);
      sitesPerWS = unsigned(floor(numSitesWS / 3 + 0.5));
      numSitesWS -= sitesPerWS;
      whiteSpaceWidths.push_back(sitesPerWS * spacing);

      movables.push_back(UINT_MAX);
      sitesPerWS = unsigned(floor(numSitesWS / 2 + 0.5));
      numSitesWS -= sitesPerWS;
      whiteSpaceWidths.push_back(sitesPerWS * spacing);

      movables.push_back(UINT_MAX);
      sitesPerWS = numSitesWS;
      whiteSpaceWidths.push_back(sitesPerWS * spacing);
    } else if (numSitesWS == 2) {
      movables.push_back(UINT_MAX);
      sitesPerWS = unsigned(floor(numSitesWS / 2 + 0.5));
      numSitesWS -= sitesPerWS;
      whiteSpaceWidths.push_back(sitesPerWS * spacing);

      movables.push_back(UINT_MAX);
      sitesPerWS = numSitesWS;
      whiteSpaceWidths.push_back(sitesPerWS * spacing);
    } else if (numSitesWS == 1) {
      movables.push_back(UINT_MAX);
      sitesPerWS = numSitesWS;
      whiteSpaceWidths.push_back(sitesPerWS * spacing);
    }

    cellIds = bin.getCellIds();
    for (unsigned i = 0; i < bin.getNumCells(); ++i) {
      movables.push_back(cellIds[i]);
    }

    vector<SmallPlacementRow> smallPlRows;
    smallPlRows.push_back(splRow);
    // assume all cells are not fixed
    bit_vector fixed(movables.size(), 0);
    SmallPlacementProblem problem(hgraph, _placement, _placed, fixed, movables,
                                  smallPlRows, whiteSpaceWidths, spBitBoards);

    if (_params.saveSmallPlProb) {
      static unsigned numInstance = 0;
      if (numInstance == 0)
        cout << "Begin saving small placement problems..." << endl;
      ostringstream baseFileName;
      baseFileName << "smallPla-r" << problem.getNumRows() << "-c"
                   << problem.getNumCells() << "-n" << problem.getNumNets()
                   << "-" << numInstance++;
      /*                           static_cast<int>(problem.getNumRows()),
                       static_cast<int>(problem.getNumCells()),
                       static_cast<int>(problem.getNumNets()),
                       static_cast<int>(numInstance++));
      */
      problem.save(baseFileName.str().c_str());
    }

    smallPlaceProbTime.stop();
    CapoPlacer::capoTimer::SmPlaceProb += smallPlaceProbTime.getUserTime();

    Timer tm18;

    SmallPlParams splParams = _params.smplParams;
    splParams.verb.setForMajStats(splParams.verb.getForMajStats() / 10);
    splParams.verb.setForActions(splParams.verb.getForActions() / 10);
    splParams.verb.setForSysRes(splParams.verb.getForSysRes() / 10);

    SmallPlacer placer(problem, splParams, _params.maxMem);

    const Placement& soln = problem.getBestSoln().placement;
    abkassert(soln.getSize() == movables.size(), "undersized placement");

    double leftMostSite = onlyRow[bin._startOffsets[0].first].getXStart() +
                          bin._startOffsets[0].second * onlyRow.getSpacing();

    tm18.stop();
    CapoPlacer::capoTimer::SmTime += tm18.getUserTime();

    for (unsigned n = 0; n < movables.size(); n++) {
      if (movables[n] != UINT_MAX) {
        _placed[movables[n]] = true;
        _cellToBinMap[movables[n]] = NULL;
        double locOffset = soln[n].x - leftMostSite;
        unsigned numSites;
        if (locOffset <= 0)
          numSites = 0;
        else
          numSites = static_cast<unsigned>(floor(locOffset / spacing));

        Point newPos(numSites * spacing + leftMostSite, yCoord);

        if (_params.wtCut) {
          Timer psTimer;
          updatePointSetsForNewlyPlacedCell(movables[n],
                                            _placement[movables[n]], newPos);
          psTimer.stop();
          CapoPlacer::capoTimer::PointSetTime += psTimer.getUserTime();
        }

        _placement[movables[n]] = newPos;
      }
    }
    cerr << "-";

    bin.setIndexToNextIndex();

    if (_params.replaceSmallBins == AtEveryLayer)
      newLayout.push_back(&bin);
    else
      _placedBins.push_back(&bin);

    return false;
    // without deleting the bin
  } else {
    // Begin Non-trivial bin refinement
    unsigned numRows = bin.getNumRows();
    Verbosity sVerb = _params.verb;
    sVerb.setForActions(sVerb.getForActions() / 2);
    sVerb.setForMajStats(sVerb.getForMajStats() / 2);
    sVerb.setForSysRes(sVerb.getForSysRes() / 2);

    if (numLeaves < minNumCellsToPart) {
      if (_params.verb.getForActions() > 4)
        cout << "Bin below current min to split...skipping" << endl;
      newLayout.push_back(&bin);
      cerr << "x";
      return false;
    }

    BaseBinSplitter* splitter;
    vector<CapoBin*> newBins;
    updateMapInfoAboutBin(&bin);

    unsigned bFactor = _params.splitterParams.numMLSets;
    unsigned lAhead = _params.lookAhead;
    if (numLeaves < 500) lAhead = 0;

    if (_mainLayerNum == 0)
      bFactor *= 3;
    else if (_mainLayerNum == 1)
      bFactor *= 2;

    bool splitHoriz = false;
    if (bin.getWidth() < bin.getHeight() * 5.0 &&
        bin.getWidth() > bin.getHeight() * 1 / 5.0 &&
        _params.splitterParams.useOracleDirection) {
      if (_params.verb.getForActions() > 4)
        cerr << "Oracle choosing cut direction" << endl;
      CutLinePositionOracle oracle(bin);
      bool hCertain = false;
      bool vCertain = false;
      unsigned hCut = oracle.getHorizontalCut(hCertain);
      double vCut = oracle.getVerticalCut(vCertain);
      if (_params.verb.getForActions() > 4)
        cerr << " hCut Idx : " << hCut << endl;
      if (_params.verb.getForActions() > 4)
        cerr << " hCut Coord : " << bin.getRow(hCut).getYCoord() << endl;
      if (_params.verb.getForActions() > 4) cerr << " vCut : " << vCut << endl;
      double hLen = oracle.cutLineLengthH(hCut);
      double vLen = oracle.cutLineLengthV(vCut);
      if (_params.verb.getForActions() > 4) cerr << " hLen : " << hLen << endl;
      if (_params.verb.getForActions() > 4) cerr << " vLen : " << vLen << endl;

      splitHoriz = hLen < vLen;
      double vertFactor = _params.vertFactor;
      const RBPCoreRow& firstRow = **bin.rowsBegin();
      if (!splitHoriz) {
        splitHoriz = (bin.getBBox().getWidth() < 3 * firstRow.getSpacing()) ||
                     (vertFactor * bin.getHeight() > bin.getWidth() &&
                      numLeaves <= _params.smallPlaceThreshold * numRows &&
                      bin.getRelativeWhitespace() <= 0.3
                      //                           && numRows <= 5
                     );
      }
      if (_params.verb.getForActions() > 4) {
        if (splitHoriz)
          cerr << "Oracle chose Horiz" << endl;
        else
          cerr << "Oracle chose Vert" << endl;
      }
    } else {
      double stretchFactor =
          bin.getAvgRowSpacing() > 1.5 * bin.getAvgCellHeight() ? 2 : 1;

      if (_params.ARStretchFactor != -1 && _params.ARStretchFactor >= 0)
        stretchFactor = _params.ARStretchFactor;

      splitHoriz = stretchFactor * bin.getHeight() > bin.getWidth();

      double vertFactor = _params.vertFactor;
      if (!splitHoriz && _params.mixedSize) {
        const vector<unsigned>& cellIds = bin.getCellIds();
        for (unsigned i = 0; i < cellIds.size(); ++i) {
          unsigned idx = cellIds[i];
          if (hgraph.isNodeMacro(idx)) {
            vertFactor = 0.;
            break;
          }
        }
      }

      const RBPCoreRow& firstRow = **bin.rowsBegin();
      if (!splitHoriz) {
        splitHoriz = (bin.getBBox().getWidth() < 3 * firstRow.getSpacing()) ||
                     (vertFactor * bin.getHeight() > bin.getWidth() &&
                      numLeaves <= _params.smallPlaceThreshold * numRows &&
                      bin.getRelativeWhitespace() <= 0.3
                      //                            && numRows <= 5
                     );
      }
    }

    // in presence of groups in regions. need to decide here how/where to
    // partition
    bool changeDefaultFlow = false;
    double splitPt, splitPtBak;
    SplitPtFixedType splitPtFixedType = SPLITPT_UNFIXED;
    bool splitDir = splitHoriz;

    if (_params.alignCutline2Obs) {
      if (bin.shouldAlignCL2Obs(splitPt, splitHoriz)) {
        changeDefaultFlow = true;
        if (lessThanDouble(_params.splitterParams.safeWS, 0.) &&
            lessThanDouble(_params.splitterParams.minLocalWS,
                           0.)) {  // uniformWS
          // splitPtFixedType = SPLITPT_FIXED;
          splitPtFixedType = SPLITPT_HINT;
        } else {
          splitPtFixedType = SPLITPT_HINT;
        }
      }
    }
    if (_params.useGrpConstr) {
      unsigned numRegionedCells = getNumRegionedCellsInBin(&bin);

      if (numRegionedCells > 0) {
        bool desiredDir = splitHoriz;
        splitDir = splitHoriz;

        if (bin.findSplitPtDir(desiredDir, groupRegionConstr, splitDir,
                               splitPtBak)) {
          changeDefaultFlow = true;
          splitHoriz = splitDir;
          splitPt = splitPtBak;
          splitPtFixedType = SPLITPT_FIXED;
        }

        if (_params.verb.getForActions() > 2) {
          cout << "changeDefaultFlow = " << changeDefaultFlow
               << " orig splitDir " << splitHoriz << " new splitDir "
               << splitDir << " splitPoint " << splitPt << " binBBox "
               << bin.getBBox() << endl;
        }
      }
    }

    if (_params.usePredeterminedCutlines) {
      if (isCutlinePredetermined(bin.getName(), splitPtBak, splitDir)) {
        changeDefaultFlow = true;
        splitHoriz = splitDir;
        splitPt = splitPtBak;
        splitPtFixedType = SPLITPT_FIXED;
      }
    }

    if (changeDefaultFlow && _params.verb.getForActions() > 1) {
      cout << "Trying to split the bin ";
      if (splitHoriz)
        cout << "horizontally";
      else
        cout << "vertically";
      cout << " at location " << splitPt << endl;
    }

    isHCut = splitHoriz;

    numRows = bin.getNumRows();
    if (numRows == 0) {
      vector<unsigned>::const_iterator cellIt;
      for (cellIt = bin.cellIdsBegin(); cellIt != bin.cellIdsEnd(); ++cellIt) {
        _placed[*cellIt] = true;
        Point newPos(bin.getBBox().xMin, bin.getBBox().yMin);
        if (_params.wtCut) {
          Timer psTimer;
          updatePointSetsForNewlyPlacedCell(*cellIt, _placement[*cellIt],
                                            newPos);
          psTimer.stop();
          CapoPlacer::capoTimer::PointSetTime += psTimer.getUserTime();
        }
        _placement[*cellIt] = newPos;
      }
      bin._canSplitBin = false;
      _placedBins.push_back(&bin);

      return false;
    }

    double origSafeWS = _params.splitterParams.safeWS;
    double origMinLocalWS = _params.splitterParams.minLocalWS;
    if (_params.mixedSize) {
      double largestMacroArea = bin.largestMacroArea();
      double totalMacroArea = bin.getTotalMacroArea();
      double totalCellArea = bin.getTotalCellArea();
      if (greaterOrEqualDouble(totalMacroArea, 0.95 * totalCellArea) ||
          greaterOrEqualDouble(largestMacroArea, 0.05 * totalCellArea)) {
        _params.splitterParams.safeWS = -1;
        _params.splitterParams.minLocalWS = -1;
      } else if (greaterThanDouble(totalMacroArea, 0.)) {
        if (_params.splitterParams.safeWS > 0.) {
          _params.splitterParams.safeWS =
              max(_params.splitterParams.safeWS, 0.10);
        }
        if (_params.splitterParams.minLocalWS > 0.) {
          _params.splitterParams.minLocalWS =
              max(_params.splitterParams.minLocalWS, 0.05);
        }
      }
    }

    splitData dataForThisSplit;
    dataForThisSplit.parent = &bin;

    if (numRows > 1 && splitHoriz) {  // Horizontal partitioning

      dataForThisSplit.splitHoriz = true;

      if (partDir == VOnly) {
        if (_params.verb.getForActions() > 4)
          cout << "Bin wants to split the other way...skipping" << endl;
        newLayout.push_back(&bin);
        cerr << "x";
        if (_params.mixedSize) {
          _params.splitterParams.safeWS = origSafeWS;
          _params.splitterParams.minLocalWS = origMinLocalWS;
        }
        return false;
      }

      if (numLeaves <= _params.smallSplitThreshold) {
        CapoSplitterParams splitterParams(_params.splitterParams);
        splitterParams.verb = sVerb;
        splitter = new SplitSmallBinHorizontally(bin, *this, splitterParams);
      } else {
        if (lAhead > 0) {
          splitter = new LookAheadSplitter(bin, *this, bFactor, lAhead, sVerb);
        } else if (partDir == HAndV || !_params.allowRowSplits) {
          CapoSplitterParams splitterParams(_params.splitterParams);
          if (bFactorAdjustCutoff >= 1. &&
              (double)numLeaves >= bFactorAdjustCutoff) {
            splitterParams.numMLSets =
                splitterParams.numMLSets *
                (unsigned)(2. + log((double)numLeaves / bFactorAdjustCutoff));
          }
          if (!changeDefaultFlow) {
            splitter = new SplitLargeBinHorizontally(bin, *this, splitterParams,
                                                     &_placement);
          } else {
            splitter = new SplitLargeBinHorizontally(bin, *this, splitterParams,
                                                     &_placement,
                                                     splitPtFixedType, splitPt);
          }
        } else {
          // collect a row of bins
          vector<CapoBin*> rowOfBins;
          rowOfBins.push_back(&bin);
          double rowTop = rowOfBins[0]->getBBox().yMax;
          double rowBot = rowOfBins[0]->getBBox().yMin;

          CapoBin* adjBin = &bin;

          while (adjBin->_leftNeighbors.size() == 1) {
            CapoBin* lNbor = adjBin->_leftNeighbors[0];

            if (lNbor->getNumCells() < 100 ||
                lNbor->getNumCells() < minNumCellsToPart ||
                !lNbor->canSplitBin())
              break;

            if ((lNbor->getBBox().yMin != rowBot) ||
                (lNbor->getBBox().yMax != rowTop))
              break;

            adjBin = lNbor;
            rowOfBins.push_back(adjBin);
          }

          adjBin = &bin;
          while (adjBin->_rightNeighbors.size() == 1) {
            CapoBin* rNbor = adjBin->_rightNeighbors[0];

            if (rNbor->getNumCells() < 100 ||
                rNbor->getNumCells() < minNumCellsToPart ||
                !rNbor->canSplitBin())
              break;

            if ((rNbor->getBBox().yMin != rowBot) ||
                (rNbor->getBBox().yMax != rowTop))
              break;

            adjBin = rNbor;
            rowOfBins.push_back(adjBin);
          }

          if (rowOfBins.size() > 1) {
            CapoSplitterParams splitterParams(_params.splitterParams);
            if (bFactorAdjustCutoff >= 1. &&
                (double)numLeaves > bFactorAdjustCutoff) {
              splitterParams.numMLSets =
                  splitterParams.numMLSets *
                  (unsigned)(2. + log((double)numLeaves / bFactorAdjustCutoff));
            }
            splitter = new SplitRowOfBins(rowOfBins, *this, splitterParams);
            if (static_cast<SplitRowOfBins*>(splitter)->wasGoodSplit()) {
              for (unsigned sb = 0; sb < rowOfBins.size(); sb++)
                rowOfBins[sb]->_canSplitBin = false;
            } else {
              cout << "Failed Row Split" << endl;
              delete splitter;
              splitter =
                  new SplitLargeBinHorizontally(bin, *this, splitterParams);
            }
          } else {
            CapoSplitterParams splitterParams(_params.splitterParams);
            if (bFactorAdjustCutoff >= 1. &&
                (double)numLeaves > bFactorAdjustCutoff) {
              splitterParams.numMLSets =
                  splitterParams.numMLSets *
                  (unsigned)(2. + log((double)numLeaves / bFactorAdjustCutoff));
            }
            splitter =
                new SplitLargeBinHorizontally(bin, *this, splitterParams);
          }
        }
      }
    } else {  // Vertical partitioning

      dataForThisSplit.splitHoriz = false;

      if (partDir == HOnly) {
        if (_params.verb.getForActions() > 4)
          cout << "Bin wants to split the other way...skipping" << endl;
        newLayout.push_back(&bin);
        cerr << "x";
        if (_params.mixedSize) {
          _params.splitterParams.safeWS = origSafeWS;
          _params.splitterParams.minLocalWS = origMinLocalWS;
        }
        return false;
      }

      if (numLeaves <= _params.smallSplitThreshold) {
        CapoSplitterParams splitterParams(_params.splitterParams);
        splitterParams.verb = sVerb;
        splitter = new SplitSmallBinVertically(bin, *this, splitterParams);
      } else {
        if (lAhead > 0) {
          splitter = new LookAheadSplitter(bin, *this, bFactor, lAhead, sVerb);
        } else {
          CapoSplitterParams splitterParams(_params.splitterParams);
          if (bFactorAdjustCutoff >= 1. &&
              (double)numLeaves > bFactorAdjustCutoff) {
            splitterParams.numMLSets =
                splitterParams.numMLSets *
                (unsigned)(2. + log((double)numLeaves / bFactorAdjustCutoff));
          }
          if (!changeDefaultFlow) {
            splitter = new SplitLargeBinVertically(bin, *this, splitterParams,
                                                   &_placement);
          } else {
            splitter = new SplitLargeBinVertically(bin, *this, splitterParams,
                                                   &_placement,
                                                   splitPtFixedType, splitPt);
          }
        }
      }
    }

    newBins = splitter->getNewBins();

    dataForThisSplit.splitRow = splitter->getSplitRow();
    dataForThisSplit.xSplit = splitter->getXSplit();

    delete splitter;

    /*
    //The following block of code was added by DAP to early detect children that
    are impossible
    //thus implying that we must floorplan the current bin.  If that
    floorplanning fails, we
    //can still merge the two bins and go up one more time.  This should prevent
    most failures.

    bool childImpossible = false;
    for(unsigned i = 0; i < newBins.size(); ++i)
    {
      if( newBins[i]->fpCondMet_Strict() )
        childImpossible = true;
    }
    if(childImpossible)
    {
      vector<CapoBin*> backtrackLayout;
      vector<CapoBin*> failedBacktracks;
      newBins.clear();
      backtrackLayout.push_back(&bin);

      if(_params.verb.getForActions() > 0)
        cout<<"Found impossible to place child bins inside refineBin,
    backtracking... "<<endl;

      doFPLayer(backtrackLayout, failedBacktracks);

      if(failedBacktracks.size() > 0)
      {
        if(_params.verb.getForActions() > 0)
          cout<<"Failed a backtracking, doing last chance backtracking...
    "<<endl; newLayout.insert( newLayout.end(), failedBacktracks.begin(),
    failedBacktracks.end());
      }
    }
    //End Early detection block
    */

    if (newBins.size() < 2) {
      if (numLeaves <= _params.smallPartThreshold)
        cerr << "e";
      else
        cerr << "E";
    } else {
      if (numLeaves <= _params.smallPartThreshold)
        cerr << "*";
      else
        cerr << ".";
    }

    if (_params.capoNE || _params.useQuad) {
      for (unsigned b = 0; b < newBins.size(); b++)
        updatePlInfoAboutBin(newBins[b]);
    }

    // by sadya to incorporate top-down congestion control
    if (useCongestionInfo && _params.tdCongestionCtl) {
      double congestionPart0, congestionPart1;
      double relCongestionP0, relCongestionP1;
      const vector<unsigned>& binCellIds = bin.getCellIds();

      _congestionMap->remCongByCells(binCellIds, true);
      for (unsigned b = 0; b < newBins.size(); b++) {
        newLayout.push_back(newBins[b]);
        updateInfoAboutBin(newBins[b]);
        updatePlInfoAboutBin(newBins[b]);
      }
      _congestionMap->addCongByCells(binCellIds, true);

      bin._canSplitBin = false;

      unsigned binNumCells = bin.getNumCells();
      unsigned threshold = 500 * _congestionMap->getAvgNumCellsPerGrid();

      if (bin.getNumCells() > 30 && binNumCells < threshold) {
        cout << "Threshhold : " << threshold << " numCells " << binNumCells
             << endl;
        congestionPart0 = _congestionMap->getCongestion(newBins[0]->getBBox());
        congestionPart1 = _congestionMap->getCongestion(newBins[1]->getBBox());

        relCongestionP0 = congestionPart0 / (congestionPart0 + congestionPart1);
        relCongestionP1 = congestionPart1 / (congestionPart0 + congestionPart1);

        cout << "Before " << congestionPart0 << "  " << congestionPart1 << endl;
        cout << relCongestionP0 << "  " << relCongestionP1 << endl;

        if (1) {
          repartitionBins(newBins[0], newBins[1], relCongestionP0,
                          relCongestionP1);
        } else  // TO DO change areas of partitions according to congestion
        {
          /*
          double cap0 = newBins[0]->getCapacity();
          double cap1 = newBins[1]->getCapacity();
          double nodeArea0 = newBins[0]->getTotalCellArea();
          double nodeArea1 = newBins[1]->getTotalCellArea();
          double totWS = cap0+cap1 - nodeArea0-nodeArea1;
          if(totWS > 0)
          {
          double ws0 = relCongestionP0*totWS;
          double ws1 = relCongestionP1*totWS;
          double reqdCap0 = nodeArea0 + ws0;
          double reqdCap1 = nodeArea1 + ws1;
          if(isHCut)  //bin is horizontally cut
          {
          unsigned splitRow    = 1;
          double   p1Cap =  bin.getContainedAreaInCoreRow(0);
          unsigned r;
          for(r = 1; r < bin.getNumRows()-1; r++)
          {
          if(p1Cap >= reqdCap1)
          break;
          p1Cap += bin.getContainedAreaInCoreRow(r);
          }
          if(r == bin.getNumRows()-1)
          splitRow = r-1;
          else
          splitRow = r;
          //reSplitBinH(bin,*newBins[0],*newBins[1],splitRow);
          }
          else
          {
          vector<double> partCaps;
          double xSplit =
          bin.findXSplit(reqdCap0,reqdCap1,0,partCaps);
          //reSplitBinV(bin,*newBins[0],*newBins[1],xSplit);
          }
          }
          */
        }
        _congestionMap->remCongByCells(binCellIds, true);
        for (unsigned b = 0; b < newBins.size(); b++) {
          updateInfoAboutBin(newBins[b]);
          updatePlInfoAboutBin(newBins[b]);
        }
        _congestionMap->addCongByCells(binCellIds, true);
        congestionPart0 = _congestionMap->getCongestion(newBins[0]->getBBox());
        congestionPart1 = _congestionMap->getCongestion(newBins[1]->getBBox());

        relCongestionP0 = congestionPart0 / (congestionPart0 + congestionPart1);
        relCongestionP1 = congestionPart1 / (congestionPart0 + congestionPart1);

        cout << "After " << congestionPart0 << "  " << congestionPart1 << endl;
        cout << relCongestionP0 << "  " << relCongestionP1 << endl << endl;
      }
      if (_params.mixedSize) {
        _params.splitterParams.safeWS = origSafeWS;
        _params.splitterParams.minLocalWS = origMinLocalWS;
      }

      return true;
    } else {
      for (unsigned b = 0; b < newBins.size(); b++) {
        newLayout.push_back(newBins[b]);
        /*
        if(_params.capoNE || _params.useQuad)
          updateMapInfoAboutBin(newBins[b]);
        else
        */
        updateInfoAboutBin(newBins[b]);
      }

      if (newBins.size() == 2) {
        dataForThisSplit.p0 = newBins[0];
        dataForThisSplit.p1 = newBins[1];
        _splits.push_back(dataForThisSplit);
      } else if (newBins.size() == 1) {
        dataForThisSplit.p0 = newBins[0];
        dataForThisSplit.p1 = NULL;
        _splits.push_back(dataForThisSplit);
      }

      bin._canSplitBin = false;

      if (_params.mixedSize) {
        _params.splitterParams.safeWS = origSafeWS;
        _params.splitterParams.minLocalWS = origMinLocalWS;
      }

      return true;
    }
  }
  abkfatal(0, "Bin slipped through the cracks");
  return false;
}
