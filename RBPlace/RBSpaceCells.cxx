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

#include <cmath>
#include <map>

#include "ABKCommon/abkassert.h"
#include "CongestionMaps/CongestionMaps.h"
#include "RBPlacement.h"

using std::map;
using std::max;
using std::min;
using std::vector;

void RBPlacement::spaceCellsEqually(RBPSubRow& sr, unsigned firstCellOffset,
                                    unsigned numCells) {
  abkassert(sr.getNumCells() - firstCellOffset >= numCells,
            "spaceCells called with numCells > #cells to right of startCell");
  abkassert(numCells > 1, "spaceCells called for a single cell");

  std::vector<unsigned>::iterator cellsBegin =
      sr._cellIds.begin() + firstCellOffset;
  std::vector<unsigned>::iterator cellsEnd = cellsBegin + numCells;

  double rangeBegin, rangeEnd;
  rangeBegin = _placement[*cellsBegin].x;  // left-edge of first cell
  unsigned numSpaces;

  if (cellsEnd == sr._cellIds.end()) {
    // right justify the last cell w/ the end of the subrow

    cellsEnd--;  // we will now only place n-1 cells...
    unsigned lastCell = *(cellsEnd);
    double lastCellW = _hgWDims->getNodeWidth(lastCell);
    _placement[lastCell].x =
        sr.getXEnd() - sr.getCoreRow().siteRoundUp(lastCellW);
    numCells--;
  }

  rangeEnd = _placement[*cellsEnd].x;  // left-edge of next cell
  numSpaces = numCells;                // a space following each cell
                         // including the right one

  double totalCellWidth = 0;
  std::vector<unsigned>::iterator srCell;
  for (srCell = cellsBegin; srCell != cellsEnd; srCell++) {
    double cellWidth = _hgWDims->getNodeWidth(*srCell);
    totalCellWidth += sr.getCoreRow().siteRoundUp(cellWidth);
  }

  // deltaSpace is the average space to put between each cell
  double deltaSpace = ((rangeEnd - rangeBegin) - totalCellWidth) / numSpaces;

  double leftEdge = rangeBegin;
  double cumExtraSpace = 0;  // extra space is accumulated as cells are
                             // placed. when the extra space is > k*sitewidth,
                             // k sites are skipped (for k any non-zero
  // integer)

  for (srCell = cellsBegin; srCell != cellsEnd; srCell++) {
    unsigned cellId = *srCell;
    _placement[cellId].x = leftEdge;

    double cellWidth = _hgWDims->getNodeWidth(cellId);
    leftEdge += sr.getCoreRow().siteRoundUp(cellWidth);

    cumExtraSpace += deltaSpace;

    int numSitesOverlap =
        static_cast<int>(floor(cumExtraSpace / sr.getSpacing()));
    leftEdge += sr.getSpacing() * numSitesOverlap;
    cumExtraSpace -= sr.getSpacing() * numSitesOverlap;
  }
}

void RBPlacement::spaceCellsWithPinDimensionsAlg1(RBPSubRow& sr,
                                                  unsigned firstCellOffset,
                                                  unsigned numCells) {
  abkassert(sr.getNumCells() - firstCellOffset >= numCells,
            "spaceCells called with numCells > #cells to right of startCell");
  abkassert(numCells > 1, "spaceCells called for a single cell");

  std::vector<unsigned>::iterator cellsBegin =
      sr._cellIds.begin() + firstCellOffset;
  std::vector<unsigned>::iterator cellsEnd = cellsBegin + numCells;

  double rangeBegin, rangeEnd;
  rangeBegin = _placement[*cellsBegin].x;  // left-edge of first cell

  if (cellsEnd == sr._cellIds.end()) {
    // right justify the last cell w/ the end of the subrow

    cellsEnd--;  // we will now only place n-1 cells...
    unsigned lastCell = *(cellsEnd);
    double lastCellW = _hgWDims->getNodeWidth(lastCell);
    _placement[lastCell].x =
        sr.getXEnd() - sr.getCoreRow().siteRoundUp(lastCellW);
    numCells--;
  }

  rangeEnd = _placement[*cellsEnd].x +
             sr.getCoreRow().siteRoundUp(_hgWDims->getNodeWidth(*cellsEnd));
  // right-edge of next cell

  std::vector<unsigned>::iterator srCell;
  vector<double> xLeft(numCells + 1), xRight(numCells + 1);
  vector<unsigned> leftPins(numCells + 1), rightPins(numCells + 1),
      extremalPins(numCells + 1);
  vector<double> rightOffset(numCells + 1);
  vector<unsigned> cellCombined(numCells + 1, 1);
  double usedSpace = 0;

  unsigned cellOrder = 0, nextCellOrder = 0, cellId, totalExtremalPins = 0;
  for (srCell = cellsBegin; srCell != cellsEnd + 1; srCell++, cellOrder++) {
    cellId = (*srCell);
    const HGFNode& cell = _hgWDims->getNodeByIdx(cellId);

    itHGFEdgeLocal e;
    unsigned edgeOffset = 0;  // grab the extremal pin coordinates
    for (e = cell.edgesBegin(); e != cell.edgesEnd(); e++, edgeOffset++) {
      Point pinOffset = _hgWDims->edgesOnNodePinOffset(
          edgeOffset, cellId, _placement.getOrient(cellId));
      if (edgeOffset == 0)
        xLeft[cellOrder] = xRight[cellOrder] = pinOffset.x;
      else {
        xLeft[cellOrder] = min(pinOffset.x, xLeft[cellOrder]);
        xRight[cellOrder] = max(pinOffset.x, xRight[cellOrder]);
      }
    }

    rightOffset[cellOrder] =
        sr.getCoreRow().siteRoundUp(_hgWDims->getNodeWidth(cellId)) -
        xRight[cellOrder];
    // leftOffset is the same as xLeft

    if (cellOrder == 0) xLeft[cellOrder] = 0;
    if (cellOrder == numCells) {
      xRight[cellOrder] += rightOffset[cellOrder];
      rightOffset[cellOrder] = 0;
    }

    edgeOffset = 0;  // calculate the number of extremal pins
    leftPins[cellOrder] = rightPins[cellOrder] = extremalPins[cellOrder] = 0;
    for (e = cell.edgesBegin(); e != cell.edgesEnd(); e++, edgeOffset++) {
      Point pinOffset = _hgWDims->edgesOnNodePinOffset(
          edgeOffset, cellId, _placement.getOrient(cellId));
      if (pinOffset.x == xLeft[cellOrder] && cellOrder > 0)
        leftPins[cellOrder]++;
      if (pinOffset.x == xRight[cellOrder] && cellOrder < numCells)
        rightPins[cellOrder]++;
      if (pinOffset.x == xLeft[cellOrder] && cellOrder > 0 ||
          pinOffset.x == xRight[cellOrder] && cellOrder < numCells)
        extremalPins[cellOrder]++;
    }
    usedSpace += xRight[cellOrder] - xLeft[cellOrder];
    totalExtremalPins += extremalPins[cellOrder];
  }

  bool needToRecompute = true;
  double deltaSpace = 0;

  while (needToRecompute) {
    needToRecompute = false;
    deltaSpace = ((rangeEnd - rangeBegin) - usedSpace) / totalExtremalPins;

    for (cellOrder = 0; cellOrder < numCells;) {
      if (cellCombined[cellOrder] == 0) {
        cellOrder++;
        continue;
      }
      nextCellOrder = cellOrder + 1;
      while (nextCellOrder < numCells && cellCombined[nextCellOrder] == 0)
        nextCellOrder++;
      if (nextCellOrder >= numCells) break;
      if (rightOffset[cellOrder] + xLeft[nextCellOrder] < deltaSpace) {
        cellOrder++;
        continue;
      }

      needToRecompute = true;
      usedSpace += rightOffset[cellOrder];
      if (cellCombined[cellOrder] == 1)
        totalExtremalPins += leftPins[cellOrder] - extremalPins[cellOrder];
      else
        totalExtremalPins -= rightPins[cellOrder];
      if (cellCombined[nextCellOrder] == 1)
        totalExtremalPins +=
            rightPins[nextCellOrder] - extremalPins[nextCellOrder];
      else
        totalExtremalPins -= leftPins[nextCellOrder];
      cellCombined[cellOrder] += cellCombined[nextCellOrder];
      cellCombined[nextCellOrder] = 0;
      xRight[cellOrder] += rightOffset[cellOrder] + xRight[nextCellOrder];
      rightOffset[cellOrder] = rightOffset[nextCellOrder];
      rightPins[cellOrder] = rightPins[nextCellOrder];
      extremalPins[cellOrder] = leftPins[cellOrder] + rightPins[cellOrder];
    }
  }

  double leftEdge = rangeBegin;
  double cumExtraSpace = 0;  // extra space is accumulated as cells are
                             // placed. when the extra space is > k*sitewidth,
                             // k sites are skipped (for k any non-zero
                             // integer)

  cellOrder = 0;
  double addedSpace = 0;

  for (srCell = cellsBegin; srCell != cellsEnd; srCell++, cellOrder++) {
    cellId = *srCell;
    if (cellCombined[cellOrder] > 0) leftEdge += addedSpace;
    _placement[cellId].x = leftEdge;

    double cellWidth = _hgWDims->getNodeWidth(cellId);
    leftEdge += sr.getCoreRow().siteRoundUp(cellWidth);

    if (cellCombined[cellOrder] == 0) continue;
    nextCellOrder = cellOrder + cellCombined[cellOrder];
    double interpinSpace =
        ((rightPins[cellOrder] * extremalPins[cellOrder]) /
             (leftPins[cellOrder] + rightPins[cellOrder]) +
         (leftPins[nextCellOrder] * extremalPins[nextCellOrder]) /
             (leftPins[nextCellOrder] + rightPins[nextCellOrder])) *
        deltaSpace;
    cumExtraSpace +=
        interpinSpace - rightOffset[cellOrder] - xLeft[nextCellOrder];

    int numSitesOverlap =
        static_cast<unsigned>(floor(cumExtraSpace / sr.getSpacing()));
    addedSpace = sr.getSpacing() * numSitesOverlap;
    cumExtraSpace -= addedSpace;
  }
}

void RBPlacement::spaceCellsWithPinDimensionsAlg2(RBPSubRow& sr,
                                                  unsigned firstCellOffset,
                                                  unsigned numCells) {
  abkassert(sr.getNumCells() - firstCellOffset >= numCells,
            "spaceCells called with numCells > #cells to right of startCell");
  abkassert(numCells > 1, "spaceCells called for a single cell");

  std::vector<unsigned>::iterator cellsBegin =
      sr._cellIds.begin() + firstCellOffset;
  std::vector<unsigned>::iterator cellsEnd = cellsBegin + numCells;

  double rangeBegin, rangeEnd;
  rangeBegin = _placement[*cellsBegin].x;  // left-edge of first cell

  if (cellsEnd == sr._cellIds.end()) {
    // right justify the last cell w/ the end of the subrow

    cellsEnd--;  // we will now only place n-1 cells...
    unsigned lastCell = *(cellsEnd);
    double lastCellW = _hgWDims->getNodeWidth(lastCell);
    _placement[lastCell].x =
        sr.getXEnd() - sr.getCoreRow().siteRoundUp(lastCellW);
    numCells--;
  }

  rangeEnd = _placement[*cellsEnd].x;  // left-edge of next cell

  double totalCellWidth = 0;
  vector<double> minFeasibleX(numCells), maxFeasibleX(numCells),
      xLeft(numCells), xRight(numCells);
  vector<unsigned> numPins(numCells, 0), numLeftPins(numCells);
  std::vector<unsigned>::iterator srCell;
  unsigned cellOrder = 0, cellId, totalPins = 0;
  for (srCell = cellsBegin; srCell != cellsEnd; srCell++, cellOrder++) {
    cellId = (*srCell);
    minFeasibleX[cellOrder] = totalCellWidth;
    double cellWidth = _hgWDims->getNodeWidth(cellId);
    totalCellWidth += sr.getCoreRow().siteRoundUp(cellWidth);
    const HGFNode& cell = _hgWDims->getNodeByIdx(cellId);

    itHGFEdgeLocal e;
    unsigned edgeOffset = 0;  // grab the extremal pin coordinates
    for (e = cell.edgesBegin(); e != cell.edgesEnd(); e++, edgeOffset++) {
      numPins[cellOrder]++;
      Point pinOffset = _hgWDims->edgesOnNodePinOffset(
          edgeOffset, cellId, _placement.getOrient(cellId));
      if (edgeOffset == 0)
        xLeft[cellOrder] = xRight[cellOrder] = pinOffset.x;
      else {
        xLeft[cellOrder] = min(pinOffset.x, xLeft[cellOrder]);
        xRight[cellOrder] = max(pinOffset.x, xRight[cellOrder]);
      }
    }
    numLeftPins[cellOrder] = totalPins;
    totalPins += numPins[cellOrder];
  }

  double feasibleMax = rangeEnd - totalCellWidth;
  cellOrder = 0;
  for (srCell = cellsBegin; srCell != cellsEnd; srCell++, cellOrder++) {
    maxFeasibleX[cellOrder] = feasibleMax;
    double cellWidth = _hgWDims->getNodeWidth(*srCell);
    feasibleMax += sr.getCoreRow().siteRoundUp(cellWidth);
  }
  // deltaSpace is the average space to put between each pin
  double deltaSpace = (rangeEnd - rangeBegin) / (totalPins - 1);

  double leftEdge = rangeBegin;
  cellOrder = 0;
  for (srCell = cellsBegin; srCell != cellsEnd; srCell++, cellOrder++) {
    cellId = *srCell;
    double leftPinCoord = numLeftPins[cellOrder] * deltaSpace,
           rightPinCoord =
               (numLeftPins[cellOrder] + numPins[cellOrder] - 1) * deltaSpace;
    double idealX =
        (leftPinCoord - xLeft[cellOrder] + rightPinCoord - xRight[cellOrder]) /
        2;

    int numSitesOverlap =
        static_cast<unsigned>(floor(idealX / sr.getSpacing()));
    idealX = sr.getSpacing() * numSitesOverlap;  // snap to site

    if (idealX < minFeasibleX[cellOrder]) idealX = minFeasibleX[cellOrder];
    if (idealX > maxFeasibleX[cellOrder]) idealX = maxFeasibleX[cellOrder];
    if (idealX > leftEdge) leftEdge = idealX;
    _placement[cellId].x = leftEdge;

    double cellWidth = _hgWDims->getNodeWidth(cellId);
    leftEdge += sr.getCoreRow().siteRoundUp(cellWidth);
  }
}

// by sadya. works on entire placement
void RBPlacement::spaceCellsEquallyInRows(void) {
  itRBPCoreRow itc;
  itRBPSubRow its;

  for (itc = _coreRows.begin(); itc != _coreRows.end(); ++itc) {
    for (its = itc->subRowsBegin(); its != itc->subRowsEnd(); ++its) {
      RBPSubRow& sr = const_cast<RBPSubRow&>(*its);
      unsigned numCells = its->getNumCells();
      spaceCellsEqually(sr, 0, numCells);
    }
  }
}

// by sadya to allocate whitespace in all rows according to congestion info
void RBPlacement::spaceCellsWCongInfoInRows() {
  CongestionMaps congMap(this);
  congMap.initializeFMatrix();
  congMap.createCongestionMapProb();
  double horizGridSize = congMap.getHorizGridSize() * 5;

  itRBPCoreRow itc;
  itRBPSubRow its;
  itRBPCellIds itn;
  vector<BBox> binBBoxes;
  vector<double> binCongestion;

  for (itc = _coreRows.begin(); itc != _coreRows.end(); ++itc) {
    binBBoxes.clear();
    binCongestion.clear();
    double sitePitch = itc->getSpacing();
    //    double rowHeight = itc->getHeight();
    for (its = itc->subRowsBegin(); its != itc->subRowsEnd(); ++its) {
      RBPSubRow& sr = const_cast<RBPSubRow&>(*its);
      double srLength = its->getLength();
      double widthPerBin = floor(horizGridSize / sitePitch) * sitePitch;
      double srXStart = its->getXStart();
      double srYMin = its->getYCoord();
      double srYMax = srYMin + its->getHeight();
      double leftEdge = srXStart;
      double totalCellWidth = 0;
      for (itn = its->cellIdsBegin(); itn != its->cellIdsEnd(); ++itn)
        totalCellWidth += _hgWDims->getNodeWidth(*itn);
      double availWS = srLength - totalCellWidth;

      while (srLength > 0) {
        double thisBinWidth = widthPerBin;
        if (thisBinWidth > srLength) thisBinWidth = srLength;

        BBox thisBinBBox;
        thisBinBBox.add(leftEdge, srYMin);
        leftEdge += thisBinWidth;
        thisBinBBox.add(leftEdge, srYMax);
        binBBoxes.push_back(thisBinBBox);
        srLength -= thisBinWidth;
      }
      unsigned numHorizBins = binBBoxes.size();
      // calculate congestion for each bin
      double totalCongestion = 0;
      for (unsigned i = 0; i < numHorizBins; ++i) {
        double thisBinCong = congMap.getCongestion(binBBoxes[i]);
        binCongestion.push_back(thisBinCong);
        totalCongestion += thisBinCong;
      }
      // now allocate WS to each horiz bin
      itn = its->cellIdsBegin();
      unsigned leftCellOffset = 0;
      for (unsigned i = 0; i < numHorizBins; ++i) {
        double thisBinWS = availWS * binCongestion[i] / totalCongestion;
        double availCellSpace = binBBoxes[i].getWidth() - thisBinWS;
        int numCells = 0;
        double thisBinNodesWidth = 0;
        unsigned leftCell = *itn;
        for (; itn != its->cellIdsEnd(); ++itn) {
          if ((thisBinNodesWidth + _hgWDims->getNodeWidth(*itn)) >
              availCellSpace)
            break;
          thisBinNodesWidth += _hgWDims->getNodeWidth(*itn);
          ++numCells;
        }
        // expand nextBin to take up unused space
        double extraspace = 0;
        if (i != (numHorizBins - 1)) {
          extraspace = availCellSpace - thisBinNodesWidth;
          extraspace = floor(extraspace / sitePitch) * sitePitch;
          binBBoxes[i + 1].xMin -= extraspace;
        }
        unsigned rightCell;
        if (numCells == 1) {
          _placement[leftCell].x = binBBoxes[i].xMin;
          leftCellOffset += numCells;
        } else if (numCells != 0) {
          rightCell = *(itn - 1);

          _placement[leftCell].x = binBBoxes[i].xMin;
          _placement[rightCell].x = binBBoxes[i].xMax - extraspace -
                                    _hgWDims->getNodeWidth(rightCell);

          spaceCellsEqually(sr, leftCellOffset, numCells - 1);
          leftCellOffset += numCells;
        }
      }
    }
  }
}

// by sadya to allocate whitespace in all rows according to congestion info
// this function works by blowing up cells and then legalizing
void RBPlacement::spaceCellsWCongInfoInRows1() {
  CongestionMaps congMap(this);
  congMap.initializeFMatrix();
  congMap.createCongestionMapProb();
  double horizGridSize = congMap.getHorizGridSize() * 10;

  itRBPCoreRow itc;
  itRBPSubRow its;
  itRBPCellIds itn;
  vector<BBox> binBBoxes;
  vector<double> binCongestion;

  for (itc = _coreRows.begin(); itc != _coreRows.end(); ++itc) {
    binBBoxes.clear();
    binCongestion.clear();
    double sitePitch = itc->getSpacing();
    //    double rowHeight = itc->getHeight();
    for (its = itc->subRowsBegin(); its != itc->subRowsEnd(); ++its) {
      map<unsigned, double> oldWidths;
      vector<unsigned> nodesChanged;
      //        RBPSubRow& sr = const_cast<RBPSubRow&>(*its);
      double srLength = its->getLength();
      double widthPerBin = floor(horizGridSize / sitePitch) * sitePitch;
      double srXStart = its->getXStart();
      double srYMin = its->getYCoord();
      double srYMax = srYMin + its->getHeight();
      double leftEdge = srXStart;
      double totalCellWidth = 0;
      for (itn = its->cellIdsBegin(); itn != its->cellIdsEnd(); ++itn)
        totalCellWidth += _hgWDims->getNodeWidth(*itn);
      double availWS = srLength - totalCellWidth;

      while (srLength > 0) {
        double thisBinWidth = widthPerBin;
        if (thisBinWidth > srLength) thisBinWidth = srLength;

        BBox thisBinBBox;
        thisBinBBox.add(leftEdge, srYMin);
        leftEdge += thisBinWidth;
        thisBinBBox.add(leftEdge, srYMax);
        binBBoxes.push_back(thisBinBBox);
        srLength -= thisBinWidth;
      }
      unsigned numHorizBins = binBBoxes.size();
      // calculate congestion for each bin
      double totalCongestion = 0;
      for (unsigned i = 0; i < numHorizBins; ++i) {
        double thisBinCong = congMap.getCongestion(binBBoxes[i]);
        binCongestion.push_back(thisBinCong);
        totalCongestion += thisBinCong;
      }
      // now allocate WS to each horiz bin
      for (unsigned i = 0; i < numHorizBins; ++i) {
        double thisBinWS = availWS * binCongestion[i] / totalCongestion;
        std::vector<unsigned> cellsInBin;
        getCellsInBBox(binBBoxes[i], cellsInBin);
        double oldCellWidths = 0;
        for (unsigned j = 0; j < cellsInBin.size(); ++j) {
          unsigned cellid = cellsInBin[j];
          Point pl = _placement[cellid];
          double width = _hgWDims->getNodeWidth(cellid);
          if (pl.x <= binBBoxes[i].xMin)
            width -= (binBBoxes[i].xMin - pl.x);
          else if ((pl.x + width) > binBBoxes[i].xMax)
            width = binBBoxes[i].xMax - pl.x;

          oldCellWidths += width;
        }
        double thisBinCurrWS = binBBoxes[i].getWidth() - oldCellWidths;
        double diffWS = thisBinWS - thisBinCurrWS;
        if (diffWS > 0) {
          double newCellWidths = oldCellWidths + diffWS;
          double expansionFactor = newCellWidths / oldCellWidths;

          for (unsigned j = 0; j < cellsInBin.size(); ++j) {
            unsigned cellid = cellsInBin[j];
            double oldWidth = _hgWDims->getNodeWidth(cellid);
            oldWidths[cellid] = oldWidth;
            double newWidth = expansionFactor * oldWidth;
            nodesChanged.push_back(cellid);
            _hgWDims->setNodeWidth(newWidth, cellid);
          }
        }
      }
      shuffleSR(its);  // legalize
      for (unsigned i = 0; i < nodesChanged.size(); ++i) {
        double oldWidth = oldWidths[nodesChanged[i]];
        _hgWDims->setNodeWidth(oldWidth, nodesChanged[i]);
      }
      shuffleSR(its);  // legalize
    }
  }
}

void RBPlacement::shiftCellsLeft() {
  itRBPCoreRow itc;
  itRBPSubRow its;
  itRBPCellIds itn;

  for (itc = _coreRows.begin(); itc != _coreRows.end(); ++itc) {
    for (its = itc->subRowsBegin(); its != itc->subRowsEnd(); ++its) {
      double srXStart = its->getXStart();
      double leftEdge = srXStart;
      for (itn = its->cellIdsBegin(); itn != its->cellIdsEnd(); ++itn) {
        unsigned cellId = *itn;
        _placement[cellId].x = leftEdge;
        leftEdge += _hgWDims->getNodeWidth(cellId);
      }
    }
  }
}
