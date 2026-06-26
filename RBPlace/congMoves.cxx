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

#include "CongestionMaps/ISPD04CongMap.h"
#include "GeomTrees/FastSteiner/random.h"
#include "RBPlacement.h"

using std::cout;
using std::endl;

void RBPlacement::doCongestionMoves(void) {
  BBox core = getBBox(false);
  unsigned numVertTiles = getNumCoreRows();
  unsigned numHorizTiles = getNumCells() / (10 * numVertTiles);

  double tileHeight = core.getHeight() / static_cast<double>(numVertTiles);
  double tileWidth = core.getWidth() / static_cast<double>(numHorizTiles);

  std::vector<std::vector<std::vector<unsigned> > > cellMap(
      numHorizTiles, std::vector<std::vector<unsigned> >(numVertTiles));
  std::vector<std::vector<double> > tileAreaSupply(
      numHorizTiles, std::vector<double>(numVertTiles, 0.));
  std::vector<std::vector<double> > tileAreaDemand(
      numHorizTiles, std::vector<double>(numVertTiles, 0.));

  double siteSpacing = _coreRows[0].getSpacing();

  for (unsigned i = getHGraph().getNumTerminals(); i < getNumCells(); ++i) {
    if (isFixed(i) || getHGraph().isNodeMacro(i)) continue;

    unsigned xIndex = static_cast<unsigned>(
        floor((getPlacement()[i].x - core.xMin) / tileWidth));
    if (xIndex == numHorizTiles) --xIndex;

    unsigned yIndex = static_cast<unsigned>(
        floor((getPlacement()[i].y - core.yMin) / tileHeight));
    if (yIndex == numVertTiles) --yIndex;

    cellMap[xIndex][yIndex].push_back(i);

    tileAreaDemand[xIndex][yIndex] +=
        ceil(getHGraph().getNodeWidth(i) / siteSpacing);
  }

  for (unsigned j = 0; j < numVertTiles; ++j) {
    unsigned k = 0;
    for (unsigned i = 0; i < numHorizTiles; ++i) {
      double tileXMin = core.xMin + tileWidth * static_cast<double>(i);
      double tileXMax = tileXMin + tileWidth;

      while (k < _coreRows[j].getNumSubRows() &&
             _coreRows[j][k].getXEnd() <= tileXMin) {
        ++k;
      }

      while (k < _coreRows[j].getNumSubRows() &&
             tileXMax >= _coreRows[j][k].getXStart()) {
        double segmentXMin = std::max(tileXMin, _coreRows[j][k].getXStart());
        double segmentXMax = std::min(tileXMax, _coreRows[j][k].getXEnd());

        tileAreaSupply[i][j] += rint((segmentXMax - segmentXMin) / siteSpacing);

        ++k;
      }
    }
  }

  unsigned mapVertTiles = std::max(400u, numVertTiles);
  unsigned mapHorizTiles = mapVertTiles;

  ISPD04CongMap theCongMap(*this, getPlacement(), mapHorizTiles, mapVertTiles);
  theCongMap.calculateDemand(getHGraph());

  cout << "begin: total overfill " << theCongMap.getTotalOverfill() << endl;
  cout << "begin: max overfill " << theCongMap.getMaxOverfill() << endl;
  cout << "begin: overfull tiles " << theCongMap.getNumOverfullTiles() << endl;

  bool swappedSomething = false;
  unsigned swapsPerformed = 0;
  unsigned passes = 0;

  while (passes < 5) {
    double maxXDemand = 0.;

    BBox tileBox;
    tileBox.yMin = core.yMin;
    tileBox.yMax = tileBox.yMin + tileHeight;
    for (unsigned j = 0; j < numVertTiles; ++j) {
      tileBox.xMin = core.xMin;
      tileBox.xMax = tileBox.xMin + tileWidth;
      for (unsigned i = 0; i < numHorizTiles; ++i) {
        double currXDemand = theCongMap.getMaxXDemand(tileBox);
        maxXDemand = std::max(maxXDemand, currXDemand);

        tileBox.xMin += tileWidth;
        tileBox.xMax += tileWidth;
      }
      tileBox.yMin += tileHeight;
      tileBox.yMax += tileHeight;
    }

    cout << "maxXDemand is " << maxXDemand << endl;

    ++passes;
    swappedSomething = false;

    tileBox.yMin = core.yMin;
    tileBox.yMax = tileBox.yMin + tileHeight;
    for (unsigned j = 0; j < numVertTiles; ++j) {
      tileBox.xMin = core.xMin;
      tileBox.xMax = tileBox.xMin + tileWidth;
      for (unsigned i = 0; i < numHorizTiles; ++i) {
        double currXDemand = theCongMap.getMaxXDemand(tileBox);

        if (greaterOrEqualDouble(currXDemand, 0.9 * maxXDemand)) {
          double bestCurrXDemand = currXDemand;
          std::vector<unsigned> swapCells(2, UINT_MAX);
          bool swapUp = false;

          // look to the tile above it there is one for swaps
          if (j < (numVertTiles - 1)) {
            BBox aboveTile = tileBox;
            aboveTile.yMin += tileHeight;
            aboveTile.yMax += tileHeight;

            for (unsigned m = 0; m < cellMap[i][j].size(); ++m) {
              unsigned currCell = cellMap[i][j][m];
              std::vector<unsigned> bothCells(2, currCell);
              for (unsigned n = 0; n < cellMap[i][j + 1].size(); ++n) {
                unsigned aboveCell = cellMap[i][j + 1][n];
                bothCells.back() = aboveCell;

                // check to see if these cells can be safely swapped
                double currOverfill =
                    std::max(0., tileAreaDemand[i][j] - tileAreaSupply[i][j]);
                double aboveOverfill = std::max(
                    0., tileAreaDemand[i][j + 1] - tileAreaSupply[i][j + 1]);
                double newOverfill = std::max(
                    0.,
                    tileAreaDemand[i][j] -
                        ceil(getHGraph().getNodeWidth(currCell) / siteSpacing) +
                        ceil(getHGraph().getNodeWidth(aboveCell) /
                             siteSpacing) -
                        tileAreaSupply[i][j]);
                double newAboveOverfill =
                    std::max(0., tileAreaDemand[i][j + 1] -
                                     ceil(getHGraph().getNodeWidth(aboveCell) /
                                          siteSpacing) +
                                     ceil(getHGraph().getNodeWidth(aboveCell) /
                                          siteSpacing) -
                                     tileAreaSupply[i][j + 1]);

                if (lessOrEqualDouble(newOverfill, currOverfill) &&
                    lessOrEqualDouble(newAboveOverfill, aboveOverfill)) {
                  theCongMap.removeDemandForAttachedNets(getHGraph(),
                                                         bothCells);
                  std::swap(_placement[currCell], _placement[aboveCell]);
                  theCongMap.addDemandForAttachedNets(getHGraph(), bothCells);

                  double newCurrXDemand = theCongMap.getMaxXDemand(tileBox);
                  double newAboveDemand = theCongMap.getMaxXDemand(aboveTile);
                  double newMaxXDemand = theCongMap.getMaxXDemand(core);

                  theCongMap.removeDemandForAttachedNets(getHGraph(),
                                                         bothCells);
                  std::swap(_placement[currCell], _placement[aboveCell]);
                  theCongMap.addDemandForAttachedNets(getHGraph(), bothCells);

                  if (lessThanDouble(newCurrXDemand, currXDemand) &&
                      lessThanDouble(newAboveDemand, currXDemand) &&
                      lessOrEqualDouble(newMaxXDemand, maxXDemand) &&
                      lessThanDouble(newCurrXDemand, bestCurrXDemand)) {
                    // set this move to the best
                    bestCurrXDemand = newCurrXDemand;
                    swapCells = bothCells;
                    swapUp = true;
                  }
                }
              }
            }
          }
          // look to the tile below for swaps
          if (j > 0) {
            BBox belowTile = tileBox;
            belowTile.yMin -= tileHeight;
            belowTile.yMax -= tileHeight;

            for (unsigned m = 0; m < cellMap[i][j].size(); ++m) {
              unsigned currCell = cellMap[i][j][m];
              std::vector<unsigned> bothCells(2, currCell);
              for (unsigned n = 0; n < cellMap[i][j - 1].size(); ++n) {
                unsigned belowCell = cellMap[i][j - 1][n];
                bothCells.back() = belowCell;

                // check to see if these cells can be safely swapped
                double currOverfill =
                    std::max(0., tileAreaDemand[i][j] - tileAreaSupply[i][j]);
                double belowOverfill = std::max(
                    0., tileAreaDemand[i][j - 1] - tileAreaSupply[i][j - 1]);
                double newOverfill = std::max(
                    0.,
                    tileAreaDemand[i][j] -
                        ceil(getHGraph().getNodeWidth(currCell) / siteSpacing) +
                        ceil(getHGraph().getNodeWidth(belowCell) /
                             siteSpacing) -
                        tileAreaSupply[i][j]);
                double newBelowOverfill = std::max(
                    0.,
                    tileAreaDemand[i][j - 1] -
                        ceil(getHGraph().getNodeWidth(belowCell) /
                             siteSpacing) +
                        ceil(getHGraph().getNodeWidth(currCell) / siteSpacing) -
                        tileAreaSupply[i][j - 1]);

                if (lessOrEqualDouble(newOverfill, currOverfill) &&
                    lessOrEqualDouble(newBelowOverfill, belowOverfill)) {
                  theCongMap.removeDemandForAttachedNets(getHGraph(),
                                                         bothCells);
                  std::swap(_placement[currCell], _placement[belowCell]);
                  theCongMap.addDemandForAttachedNets(getHGraph(), bothCells);

                  double newCurrXDemand = theCongMap.getMaxXDemand(tileBox);
                  double newBelowDemand = theCongMap.getMaxXDemand(belowTile);
                  double newMaxXDemand = theCongMap.getMaxXDemand(core);

                  theCongMap.removeDemandForAttachedNets(getHGraph(),
                                                         bothCells);
                  std::swap(_placement[currCell], _placement[belowCell]);
                  theCongMap.addDemandForAttachedNets(getHGraph(), bothCells);

                  if (lessThanDouble(newCurrXDemand, currXDemand) &&
                      lessThanDouble(newBelowDemand, currXDemand) &&
                      lessOrEqualDouble(newMaxXDemand, maxXDemand) &&
                      lessThanDouble(newCurrXDemand, bestCurrXDemand)) {
                    // set this move to the best
                    bestCurrXDemand = newCurrXDemand;
                    swapCells = bothCells;
                    swapUp = false;
                  }
                }
              }
            }
          }

          // examine the swap we found, if any, and apply it
          if (lessThanDouble(bestCurrXDemand, currXDemand)) {
            theCongMap.removeDemandForAttachedNets(getHGraph(), swapCells);
            Point temp = getPlacement()[swapCells.front()];
            setLocation(swapCells.front(), getPlacement()[swapCells.back()]);
            setLocation(swapCells.back(), temp);
            theCongMap.addDemandForAttachedNets(getHGraph(), swapCells);
            tileAreaDemand[i][j] -=
                ceil(getHGraph().getNodeWidth(swapCells.front()) / siteSpacing);
            tileAreaDemand[i][j] +=
                ceil(getHGraph().getNodeWidth(swapCells.back()) / siteSpacing);
            if (swapUp) {
              tileAreaDemand[i][j + 1] += ceil(
                  getHGraph().getNodeWidth(swapCells.front()) / siteSpacing);
              tileAreaDemand[i][j + 1] -= ceil(
                  getHGraph().getNodeWidth(swapCells.back()) / siteSpacing);
            } else {
              tileAreaDemand[i][j - 1] += ceil(
                  getHGraph().getNodeWidth(swapCells.front()) / siteSpacing);
              tileAreaDemand[i][j - 1] -= ceil(
                  getHGraph().getNodeWidth(swapCells.back()) / siteSpacing);
            }
            ++swapsPerformed;
            swappedSomething = true;
          }
        }
        tileBox.xMin += tileWidth;
        tileBox.xMax += tileWidth;
      }
      tileBox.yMin += tileHeight;
      tileBox.yMax += tileHeight;
    }

    if (swappedSomething) {
      // legalization may be necessary
      for (itRBPCoreRow itc = _coreRows.begin(); itc != _coreRows.end();
           ++itc) {
        for (itRBPSubRow its = itc->subRowsBegin(); its != itc->subRowsEnd();
             ++its) {
          shuffleSR(its);
        }
      }
    } else {
      break;
    }
  }

  cout << "swaps: " << swapsPerformed << endl;
  cout << "end: total overfill " << theCongMap.getTotalOverfill() << endl;
  cout << "end: max overfill " << theCongMap.getMaxOverfill() << endl;
  cout << "end: overfull tiles " << theCongMap.getNumOverfullTiles() << endl;
}
