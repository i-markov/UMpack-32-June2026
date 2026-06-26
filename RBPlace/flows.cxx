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

#include <cmath>
#include <fstream>

#include "ABKCommon/abkfunc.h"
#include "ABKCommon/abktimer.h"
#include "CongestionMaps/ISPD06DensityMap.h"
#include "Ctainers/bitBoard.h"
#ifdef USEFLOW
#include "CS2/flowSolver.h"
#endif
#include "Placement/placeWOri.h"
#include "RBPlacement.h"

using std::cout;
using std::endl;
using std::make_pair;
using std::max;
using std::min;
using std::ofstream;
using std::pair;
using std::size_t;
using std::vector;

void addToCoreCols(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                   vector<vector<unsigned> > &coreCols, const BBox &coreArea,
                   double colSpacing, const unsigned cellId) {
  unsigned angle = pl.getOrient(cellId).getAngle();
  bool notRotated = angle == 0 || angle == 180;
  double nodeWidth =
      notRotated ? hg.getNodeWidth(cellId) : hg.getNodeHeight(cellId);
  double nodeHeight =
      notRotated ? hg.getNodeHeight(cellId) : hg.getNodeWidth(cellId);
  unsigned colsWide = static_cast<unsigned>(ceil(nodeWidth / colSpacing));

  if (equalDouble(nodeWidth, 0.) || equalDouble(nodeHeight, 0.)) return;

  size_t xMinIdx = 0, xMaxIdx = 0;
  if (pl[cellId].x < coreArea.xMin) {
    xMinIdx = 0;
    unsigned colsBelow = static_cast<unsigned>(
        ceil((coreArea.xMin - pl[cellId].x) / colSpacing));
    xMaxIdx = (colsBelow >= colsWide) ? 0 : colsWide - colsBelow;
  } else {
    xMinIdx = static_cast<unsigned>(
        floor((pl[cellId].x - coreArea.xMin) / colSpacing));
    xMaxIdx = xMinIdx + colsWide;
    unsigned xMaxIdx2 = static_cast<unsigned>(
        floor((pl[cellId].x + nodeWidth - coreArea.xMin) / colSpacing));
    if (xMaxIdx <= xMaxIdx2) {
      ++xMaxIdx;
    }
  }
  xMinIdx = min(xMinIdx, coreCols.size());
  xMaxIdx = min(xMaxIdx, coreCols.size());

  for (size_t i = xMinIdx; i < xMaxIdx; ++i) {
    // remove cell from coreCols[i]
    pair<vector<unsigned>::iterator, vector<unsigned>::iterator> range =
        equal_range(coreCols[i].begin(), coreCols[i].end(), cellId,
                    CompareYWOri(pl));

    coreCols[i].insert(range.first, cellId);
  }
}

void removeFromCoreCols(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                        vector<vector<unsigned> > &coreCols,
                        const BBox &coreArea, double colSpacing,
                        const unsigned cellId) {
  unsigned angle = pl.getOrient(cellId).getAngle();
  bool notRotated = angle == 0 || angle == 180;
  double nodeWidth =
      notRotated ? hg.getNodeWidth(cellId) : hg.getNodeHeight(cellId);
  double nodeHeight =
      notRotated ? hg.getNodeHeight(cellId) : hg.getNodeWidth(cellId);
  unsigned colsWide = static_cast<unsigned>(ceil(nodeWidth / colSpacing));

  if (equalDouble(nodeWidth, 0.) || equalDouble(nodeHeight, 0.)) return;

  size_t xMinIdx = 0, xMaxIdx = 0;
  if (pl[cellId].x < coreArea.xMin) {
    xMinIdx = 0;
    unsigned colsBelow = static_cast<unsigned>(
        ceil((coreArea.xMin - pl[cellId].x) / colSpacing));
    xMaxIdx = (colsBelow >= colsWide) ? 0 : colsWide - colsBelow;
  } else {
    xMinIdx = static_cast<unsigned>(
        floor((pl[cellId].x - coreArea.xMin) / colSpacing));
    xMaxIdx = xMinIdx + colsWide;
    unsigned xMaxIdx2 = static_cast<unsigned>(
        floor((pl[cellId].x + nodeWidth - coreArea.xMin) / colSpacing));
    if (xMaxIdx <= xMaxIdx2) {
      ++xMaxIdx;
    }
  }
  xMinIdx = min(xMinIdx, coreCols.size());
  xMaxIdx = min(xMaxIdx, coreCols.size());

  for (size_t i = xMinIdx; i < xMaxIdx; ++i) {
    // remove cell from coreCols[i]
    pair<vector<unsigned>::iterator, vector<unsigned>::iterator> range =
        equal_range(coreCols[i].begin(), coreCols[i].end(), cellId,
                    CompareYWOri(pl));

    for (vector<unsigned>::iterator j = range.first; j != range.second; ++j) {
      if (*j == cellId) {
        coreCols[i].erase(j);
        break;
      }
    }
  }
}

#ifdef USEFLOW
void RBPlacement::doYFlowInBBox(ISPD06DensityMap *map,
                                vector<vector<unsigned> > &coreCols,
                                const BBox &theBox,
                                vector<unsigned> *allowedMovables, double Ydist,
                                MaxMem *maxMem) {
  if (_params.verb.getForActions() > 0) {
    cout << " Y direction" << endl;
  }

  Timer buildProblemTime;

  BBox core = getBBox(false);
  BBox realBBox = core.intersect(theBox);

  double rowHeight = _coreRows[0].getHeight();

  if (greaterThanDouble(Ydist, 0.) && lessThanDouble(Ydist, rowHeight)) {
    cout << "  maximum Y distance (" << Ydist << ") too small" << endl;
    return;
  }

  // align the bbox to row and site boundaries
  realBBox.xMin = core.xMin + ceil((realBBox.xMin - core.xMin) /
                                   _coreRows[0].getSpacing()) *
                                  _coreRows[0].getSpacing();
  realBBox.xMax = core.xMin + floor((realBBox.xMax - core.xMin) /
                                    _coreRows[0].getSpacing()) *
                                  _coreRows[0].getSpacing();
  realBBox.yMin =
      core.yMin + ceil((realBBox.yMin - core.yMin) / rowHeight) * rowHeight;
  realBBox.yMax =
      core.yMin + floor((realBBox.yMax - core.yMin) / rowHeight) * rowHeight;

  vector<bool> allowed(getHGraph().getNumNodes(), false);
  if (allowedMovables != NULL) {
    for (unsigned i = 0; i < allowedMovables->size(); ++i) {
      allowed[(*allowedMovables)[i]] = true;
    }
  }

  // pick out all the movable cells within the "realBBox"
  BitBoard movableCellsBoard(getHGraph().getNumNodes());

  for (unsigned i = 0; i < _coreRows.size(); ++i) {
    if (lessThanDouble(_coreRows[i].getYCoord(), realBBox.yMin)) continue;
    if (greaterThanDouble(_coreRows[i].getYCoord(), realBBox.yMax)) break;

    for (unsigned j = 0; j < _coreRows[i].getAllCells().size(); ++j) {
      unsigned cellId = _coreRows[i].getAllCells()[j];

      if (allowedMovables != NULL && !allowed[cellId]) continue;

      if (lessThanDouble(_placement[cellId].x, realBBox.xMin)) continue;
      if (greaterThanDouble(_placement[cellId].x, realBBox.xMax)) break;

      if (lessOrEqualDouble(realBBox.yMin, _placement[cellId].y) &&
          lessOrEqualDouble(_placement[cellId].y, realBBox.yMax) &&
          !_isFixed[cellId]) {
        movableCellsBoard.setBit(cellId);
      }
    }
  }

  const auto &movableCells = movableCellsBoard.getIndicesOfSetBits();

  if (movableCells.size() == 0) return;

  // now that we have all the movable cells, find all of the nets that touch
  // them
  BitBoard *edgesBoard = new BitBoard(getHGraph().getNumEdges());
  for (unsigned i = 0; i < movableCells.size(); ++i) {
    const HGFNode &node = getHGraph().getNodeByIdx(movableCells[i]);
    for (itHGFEdgeLocal edge = node.edgesBegin(); edge != node.edgesEnd();
         ++edge) {
      edgesBoard->setBit((*edge)->getIndex());
    }
  }

  const auto &affectedEdges = edgesBoard->getIndicesOfSetBits();

  // build a list of necessaryEdges,
  // those that will change HPWL based on the movement
  vector<unsigned> *necessaryEdges = new vector<unsigned>(0);
  for (unsigned i = 0; i < affectedEdges.size(); ++i) {
    double minY = DBL_MAX;
    double maxY = -DBL_MAX;
    const HGFEdge &edge = getHGraph().getEdgeByIdx(affectedEdges[i]);
    unsigned offset = 0;
    for (itHGFNodeLocal node = edge.nodesBegin(); node != edge.nodesEnd();
         ++node, ++offset) {
      unsigned nodeId = (*node)->getIndex();
      // only concerned with fixed cells
      if (movableCellsBoard.isBitSet(nodeId)) continue;

      Point pinOffset = getHGraph().nodesOnEdgePinOffset(
          offset, affectedEdges[i], _placement.getOrient(nodeId));
      Point pinPos = _placement[nodeId] + pinOffset;

      minY = min(minY, pinPos.y);
      maxY = max(maxY, pinPos.y);
    }

    // now we have the span of the fixed cells on this net
    // if this span contains the span of the bbox, this net will not change
    // given movement of cells in the bbox
    if (greaterThanDouble(minY, realBBox.yMin) ||
        lessThanDouble(maxY, realBBox.yMax)) {
      necessaryEdges->push_back(affectedEdges[i]);
    }
  }

  delete edgesBoard;

  // make a map so we can go from cellId to its position in movableCells
  vector<unsigned> *cellIdToPos =
      new vector<unsigned>(getHGraph().getNumNodes(), UINT_MAX - 1);
  for (unsigned i = 0; i < movableCells.size(); ++i) {
    (*cellIdToPos)[movableCells[i]] = i;
  }

  // now we have lists of the movable cells as well as
  // all the necessary edges, so we can make all the constraints

  unsigned long sourceId = movableCells.size() + 1;
  unsigned long sinkId = sourceId + 1;
  unsigned long bottomEdge = sinkId + 1;
  unsigned long topEdge = bottomEdge + 1;
  unsigned long firstEdge = topEdge + 1;
  long unlimited = necessaryEdges->size();

  // map (from,to) -> (cap,cost)
  edgeHashMap *directedEdges = new edgeHashMap();

  // add the edges constraining the left and right sides of the entire area
  (*directedEdges)[make_pair(bottomEdge, topEdge)] =
      make_pair(unlimited, static_cast<long>(realBBox.getHeight() / rowHeight));
  (*directedEdges)[make_pair(topEdge, bottomEdge)] = make_pair(
      unlimited, -static_cast<long>(realBBox.getHeight() / rowHeight));

  // add the edges attaching the left edges to the sink
  // and the right edges to the source
  for (unsigned long i = 0; i < necessaryEdges->size(); ++i) {
    (*directedEdges)[make_pair(firstEdge + 2L * i, sinkId)] = make_pair(1L, 0L);
    (*directedEdges)[make_pair(sourceId, firstEdge + 2L * i + 1L)] =
        make_pair(1L, 0L);
  }

  // add the two edges for every pin
  for (unsigned long i = 0; i < necessaryEdges->size(); ++i) {
    unsigned edgeId = (*necessaryEdges)[i];
    const HGFEdge &edge = getHGraph().getEdgeByIdx(edgeId);
    unsigned offset = 0;
    unsigned fixedPins = 0;
    double fixedPinYMin = DBL_MAX;
    double fixedPinYMax = -DBL_MAX;
    for (itHGFNodeLocal n = edge.nodesBegin(); n != edge.nodesEnd();
         ++n, ++offset) {
      unsigned nodeId = (*n)->getIndex();
      Point pinOffset = getHGraph().nodesOnEdgePinOffset(
          offset, edgeId, _placement.getOrient(nodeId));
      if (movableCellsBoard.isBitSet(nodeId)) {
        // edge constraining left side
        (*directedEdges)[make_pair((*cellIdToPos)[nodeId] + 1L,
                                   firstEdge + 2L * i)] =
            make_pair(unlimited,
                      static_cast<long>(floor(pinOffset.y / rowHeight)));
        // edge constraining right side
        (*directedEdges)[make_pair(firstEdge + 2L * i + 1L,
                                   (*cellIdToPos)[nodeId] + 1L)] =
            make_pair(unlimited,
                      -static_cast<long>(floor(pinOffset.y / rowHeight)));
      } else {
        // fixed pin, collect all and have at most 2 constraints
        ++fixedPins;
        double yPinLoc = _placement[nodeId].y + pinOffset.y;
        fixedPinYMin = min(fixedPinYMin, yPinLoc);
        fixedPinYMax = max(fixedPinYMax, yPinLoc);
      }
    }
    if (fixedPins > 0) {
      // edge constraining bottom side
      (*directedEdges)[make_pair(bottomEdge, firstEdge + 2L * i)] = make_pair(
          unlimited,
          static_cast<long>(floor((fixedPinYMin - realBBox.yMin) / rowHeight)));
      // edge constraining top side
      (*directedEdges)[make_pair(firstEdge + 2L * i + 1L, topEdge)] = make_pair(
          unlimited,
          static_cast<long>(ceil((realBBox.yMax - fixedPinYMax) / rowHeight)));
    }
  }

  // for each cell, constrain it between the cells below and above it
  for (unsigned i = 0; i < coreCols.size(); ++i) {
    double coreColXMin =
        core.xMin + static_cast<double>(i) * _coreRows[0].getSpacing();
    double coreColXMax = coreColXMin + _coreRows[0].getSpacing();
    RBPSubRow targetSubRow(coreColXMin, coreColXMax);

    for (unsigned j = 0; j < coreCols[i].size(); ++j) {
      unsigned cellId = coreCols[i][j];

      if (!movableCellsBoard.isBitSet(cellId)) continue;

      // macros don't fit into the column structures perfectly, so we need to
      // check if relationships that look like constraints are real first before
      // adding them to the problem

      unsigned angle = _placement.getOrient(cellId).getAngle();
      bool notRotated = angle == 0 || angle == 180;
      double cellWidth = notRotated ? getHGraph().getNodeWidth(cellId)
                                    : getHGraph().getNodeHeight(cellId);
      double cellHeight = notRotated ? getHGraph().getNodeHeight(cellId)
                                     : getHGraph().getNodeWidth(cellId);
      cellHeight = ceil(cellHeight / rowHeight) * rowHeight;
      abkfatal(cellHeight >= rowHeight, "bad cell height");
      double cellXMin = _placement[cellId].x;
      double cellXMax = cellXMin + cellWidth;

      if (!getHGraph().isNodeMacro(cellId)) {
        // must search up and down
        // stop when we find a cell, a missing site, or a fixed object
        double yUpper = DBL_MAX;
        for (unsigned k = j + 1; k < coreCols[i].size(); ++k) {
          // stop when we see a cell or a fixed object
          if (!getHGraph().isNodeMacro(coreCols[i][k]) ||
              !movableCellsBoard.isBitSet(coreCols[i][k])) {
            yUpper = _placement[coreCols[i][k]].y;
            break;
          }
        }
        double yLower = -DBL_MAX;
        for (unsigned k = j; k > 0;) {
          --k;
          // stop when we see a cell or a fixed object
          if (!getHGraph().isNodeMacro(coreCols[i][k]) ||
              !movableCellsBoard.isBitSet(coreCols[i][k])) {
            unsigned prevAngle =
                _placement.getOrient(coreCols[i][k]).getAngle();
            bool prevNotRotated = prevAngle == 0 || prevAngle == 180;
            double prevCellHeight =
                prevNotRotated ? getHGraph().getNodeHeight(coreCols[i][k])
                               : getHGraph().getNodeWidth(coreCols[i][k]);
            prevCellHeight = ceil(prevCellHeight / rowHeight) * rowHeight;
            yLower = _placement[coreCols[i][k]].y + prevCellHeight;
            break;
          }
        }
        // found the upper and lower bounds for the cell based on std cells and
        // fixed objects, now look for missing sites
        unsigned coreRowIdx = findCoreRowIdx(_placement[cellId]);
        if (coreRowIdx != UINT_MAX) {
          // going up
          unsigned idx = coreRowIdx + 1;
          while (idx < _coreRows.size()) {
            if (greaterOrEqualDouble(_coreRows[idx].getYCoord(), yUpper)) break;
            pair<itRBPSubRow, itRBPSubRow> subRowIts =
                equal_range(_coreRows[idx].subRowsBegin(),
                            _coreRows[idx].subRowsEnd(), targetSubRow);
            if (subRowIts.first == subRowIts.second) {
              // site missing, must add a constraint
              double yMaximum =
                  min(realBBox.yMax, _coreRows[idx].getYCoord() - rowHeight);
              pair<unsigned long, unsigned long> fromTo =
                  make_pair(topEdge, (*cellIdToPos)[cellId] + 1L);
              if (directedEdges->find(fromTo) == directedEdges->end()) {
                (*directedEdges)[fromTo] = make_pair(
                    unlimited, -static_cast<long>(ceil(
                                   (realBBox.yMax - yMaximum) / rowHeight)));
              } else {
                long oldCost = (*directedEdges)[fromTo].second;
                (*directedEdges)[fromTo].second = min(
                    oldCost, -static_cast<long>(
                                 ceil((realBBox.yMax - yMaximum) / rowHeight)));
              }
              break;
            }
            ++idx;
          }
          // going down
          idx = coreRowIdx;
          while (idx > 0) {
            --idx;
            if (lessThanDouble(_coreRows[idx].getYCoord(), yLower)) break;
            pair<itRBPSubRow, itRBPSubRow> subRowIts =
                equal_range(_coreRows[idx].subRowsBegin(),
                            _coreRows[idx].subRowsEnd(), targetSubRow);
            if (subRowIts.first == subRowIts.second) {
              // site missing, must add a constraint
              double yMinimum =
                  max(realBBox.yMin, _coreRows[idx].getYCoord() + rowHeight);
              pair<unsigned long, unsigned long> fromTo =
                  make_pair((*cellIdToPos)[cellId] + 1L, bottomEdge);
              if (directedEdges->find(fromTo) == directedEdges->end()) {
                (*directedEdges)[fromTo] = make_pair(
                    unlimited, -static_cast<long>(ceil(
                                   (yMinimum - realBBox.yMin) / rowHeight)));
              } else {
                long oldCost = (*directedEdges)[fromTo].second;
                (*directedEdges)[fromTo].second = min(
                    oldCost, -static_cast<long>(
                                 ceil((yMinimum - realBBox.yMin) / rowHeight)));
              }
              break;
            }
          }
        }
      }

      unsigned k = j;
      // when a cell is being bound by a macro in this core column,
      // it could possibly be bound by several
      bool boundedBySomething = false;
      double boundedXMin = DBL_MAX;
      double boundedXMax = -DBL_MAX;
      while (true) {
        if (k == 0) {
          if (!boundedBySomething) {
            // if we have gotten here, the only thing constraining
            // cellId to the bottom is the bottom  of the core
            // bound on the bottom by the bottom of the chip
            double yMinimum = max(realBBox.yMin, core.yMin);
            pair<unsigned long, unsigned long> fromTo =
                make_pair((*cellIdToPos)[cellId] + 1L, bottomEdge);
            if (directedEdges->find(fromTo) == directedEdges->end()) {
              (*directedEdges)[fromTo] = make_pair(
                  unlimited, -static_cast<long>(
                                 ceil((yMinimum - realBBox.yMin) / rowHeight)));
            } else {
              long oldCost = (*directedEdges)[fromTo].second;
              (*directedEdges)[fromTo].second = min(
                  oldCost, -static_cast<long>(
                               ceil((yMinimum - realBBox.yMin) / rowHeight)));
            }
          }
          break;
        } else {
          // see if the cell that is at coreCols[i][k-1] could block cellId
          unsigned prevCell = coreCols[i][k - 1];
          unsigned prevAngle = _placement.getOrient(prevCell).getAngle();
          bool prevNotRotated = prevAngle == 0 || prevAngle == 180;
          double prevWidth = prevNotRotated
                                 ? getHGraph().getNodeWidth(prevCell)
                                 : getHGraph().getNodeHeight(prevCell);
          double prevHeight = prevNotRotated
                                  ? getHGraph().getNodeHeight(prevCell)
                                  : getHGraph().getNodeWidth(prevCell);
          prevHeight = ceil(prevHeight / rowHeight) * rowHeight;
          double prevXMin = _placement[prevCell].x;
          double prevXMax = prevXMin + prevWidth;

          if (greaterOrEqualDouble(prevXMin, cellXMax) ||
              greaterOrEqualDouble(cellXMin, prevXMax) ||
              (allowedMovables != NULL && !allowed[prevCell])) {
            // these two cells could never overlap, try the cell before this one
            --k;
          } else {
            // this is a true constraint
            if (greaterOrEqualDouble(max(coreColXMin, prevXMin), boundedXMin) &&
                lessOrEqualDouble(min(coreColXMax, prevXMax), boundedXMax)) {
              // no need to bound against this object,
              // as it is obscured by the object(s) we already bounded against
              --k;
            } else {
              if (movableCellsBoard.isBitSet(prevCell)) {
                // bound on the bottom by the movable object
                pair<unsigned long, unsigned long> fromTo = make_pair(
                    (*cellIdToPos)[cellId] + 1L, (*cellIdToPos)[prevCell] + 1L);
                if (directedEdges->find(fromTo) == directedEdges->end()) {
                  (*directedEdges)[fromTo] = make_pair(
                      unlimited,
                      -static_cast<long>(ceil(prevHeight / rowHeight)));
                } else {
                  long oldCost = (*directedEdges)[fromTo].second;
                  (*directedEdges)[fromTo].second =
                      min(oldCost,
                          -static_cast<long>(ceil(prevHeight / rowHeight)));
                }
                boundedBySomething = true;
                boundedXMin = max(coreColXMin, min(boundedXMin, prevXMin));
                boundedXMax = min(coreColXMax, max(boundedXMax, prevXMax));
              } else if (lessThanDouble(_placement[prevCell].y,
                                        _placement[cellId].y)) {
                // bound on the bottom by the fixed object
                double yMinimum =
                    max(realBBox.yMin, _placement[prevCell].y + prevHeight);
                yMinimum = max(yMinimum, core.yMin);
                pair<unsigned long, unsigned long> fromTo =
                    make_pair((*cellIdToPos)[cellId] + 1L, bottomEdge);
                if (directedEdges->find(fromTo) == directedEdges->end()) {
                  (*directedEdges)[fromTo] = make_pair(
                      unlimited, -static_cast<long>(ceil(
                                     (yMinimum - realBBox.yMin) / rowHeight)));
                } else {
                  long oldCost = (*directedEdges)[fromTo].second;
                  (*directedEdges)[fromTo].second = min(
                      oldCost, -static_cast<long>(ceil(
                                   (yMinimum - realBBox.yMin) / rowHeight)));
                }
                boundedBySomething = true;
                boundedXMin = max(coreColXMin, min(boundedXMin, prevXMin));
                boundedXMax = min(coreColXMax, max(boundedXMax, prevXMax));
              }
              --k;
              if (lessOrEqualDouble(boundedXMin, coreColXMin) &&
                  greaterOrEqualDouble(boundedXMax, coreColXMax)) {
                break;
              }
            }
          }
        }
      }

      k = j;
      boundedBySomething = false;
      boundedXMin = DBL_MAX;
      boundedXMax = -DBL_MAX;
      while (true) {
        if ((k + 1) == coreCols[i].size()) {
          if (!boundedBySomething) {
            // if we have gotten here, the only thing constraining
            // cellId above is the top of the core
            // bound above by the top of the chip
            double yMaximum = min(realBBox.yMax, core.yMax - cellHeight);
            pair<unsigned long, unsigned long> fromTo =
                make_pair(topEdge, (*cellIdToPos)[cellId] + 1L);
            if (directedEdges->find(fromTo) == directedEdges->end()) {
              (*directedEdges)[fromTo] = make_pair(
                  unlimited, -static_cast<long>(
                                 ceil((realBBox.yMax - yMaximum) / rowHeight)));
            } else {
              long oldCost = (*directedEdges)[fromTo].second;
              (*directedEdges)[fromTo].second = min(
                  oldCost, -static_cast<long>(
                               ceil((realBBox.yMax - yMaximum) / rowHeight)));
            }
          }
          break;
        } else {
          // see if the cell that is at coreCols[i][k+1] could block cellId
          unsigned nextCell = coreCols[i][k + 1];
          unsigned nextAngle = _placement.getOrient(nextCell).getAngle();
          bool nextNotRotated = nextAngle == 0 || nextAngle == 180;
          double nextWidth = nextNotRotated
                                 ? getHGraph().getNodeWidth(nextCell)
                                 : getHGraph().getNodeHeight(nextCell);
          double nextXMin = _placement[nextCell].x;
          double nextXMax = nextXMin + nextWidth;

          if (greaterOrEqualDouble(nextXMin, cellXMax) ||
              greaterOrEqualDouble(cellXMin, nextXMax) ||
              (allowedMovables != NULL && !allowed[nextCell])) {
            // these two cells could never overlap, try the cell after this one
            ++k;
          } else {
            // this is a true constraint
            if (greaterOrEqualDouble(max(coreColXMin, nextXMin), boundedXMin) &&
                lessOrEqualDouble(min(coreColXMax, nextXMax), boundedXMax)) {
              // no need to bound against this object,
              // as it is obscured by the object(s) we already bounded against
              ++k;
            } else {
              if (movableCellsBoard.isBitSet(nextCell)) {
                // bound above by the movable
                pair<unsigned long, unsigned long> fromTo = make_pair(
                    (*cellIdToPos)[nextCell] + 1L, (*cellIdToPos)[cellId] + 1L);
                if (directedEdges->find(fromTo) == directedEdges->end()) {
                  (*directedEdges)[fromTo] = make_pair(
                      unlimited,
                      -static_cast<long>(ceil(cellHeight / rowHeight)));
                } else {
                  long oldCost = (*directedEdges)[fromTo].second;
                  (*directedEdges)[fromTo].second =
                      min(oldCost,
                          -static_cast<long>(ceil(cellHeight / rowHeight)));
                }
                boundedBySomething = true;
                boundedXMin = max(coreColXMin, min(boundedXMin, nextXMin));
                boundedXMax = min(coreColXMax, max(boundedXMax, nextXMax));
              } else if (greaterThanDouble(_placement[nextCell].y,
                                           _placement[cellId].y)) {
                // bound above by the fixed object
                double yMaximum =
                    min(realBBox.yMax, _placement[nextCell].y - cellHeight);
                yMaximum = min(yMaximum, core.yMax - cellHeight);
                pair<unsigned long, unsigned long> fromTo =
                    make_pair(topEdge, (*cellIdToPos)[cellId] + 1L);
                if (directedEdges->find(fromTo) == directedEdges->end()) {
                  (*directedEdges)[fromTo] = make_pair(
                      unlimited, -static_cast<long>(ceil(
                                     (realBBox.yMax - yMaximum) / rowHeight)));
                } else {
                  long oldCost = (*directedEdges)[fromTo].second;
                  (*directedEdges)[fromTo].second = min(
                      oldCost, -static_cast<long>(ceil(
                                   (realBBox.yMax - yMaximum) / rowHeight)));
                }
                boundedBySomething = true;
                boundedXMin = max(coreColXMin, min(boundedXMin, nextXMin));
                boundedXMax = min(coreColXMax, max(boundedXMax, nextXMax));
              }
              ++k;
              if (lessOrEqualDouble(boundedXMin, coreColXMin) &&
                  greaterOrEqualDouble(boundedXMax, coreColXMax)) {
                break;
              }
            }
          }
        }
      }
    }
  }

  /*
     if(map != NULL)
     {
       // add additional constraints for each cell if we have the grid available
       for(unsigned long i = 0; i < movableCells.size(); ++i)
       {
         BBox grid = map->getGridCell(_placement[movableCells[i]]);

         pair< unsigned long, unsigned long > fromTo = make_pair(i + 1L,
     bottomEdge); double yMinimum = grid.yMin; if(directedEdges->find(fromTo) ==
     directedEdges->end())
         {
           (*directedEdges)[fromTo] = make_pair(unlimited,
     -static_cast<long>(ceil((yMinimum - realBBox.yMin)/rowHeight)));
         }
         else
         {
           long oldCost = (*directedEdges)[fromTo].second;
           (*directedEdges)[fromTo].second = min(oldCost,
     -static_cast<long>(ceil((yMinimum - realBBox.yMin)/rowHeight)));
         }

         fromTo = make_pair(topEdge, i + 1L);
         double yMaximum = grid.yMax;
         if(directedEdges->find(fromTo) == directedEdges->end())
         {
           (*directedEdges)[fromTo] = make_pair(unlimited,
     -static_cast<long>(ceil((realBBox.yMax - yMaximum)/rowHeight)));
         }
         else
         {
           long oldCost = (*directedEdges)[fromTo].second;
           (*directedEdges)[fromTo].second = min(oldCost,
     -static_cast<long>(ceil((realBBox.yMax - yMaximum)/rowHeight)));
         }
       }
     }
  */

  if (greaterThanDouble(Ydist, 0.)) {
    for (unsigned long i = 0; i < movableCells.size(); ++i) {
      pair<unsigned long, unsigned long> fromTo = make_pair(i + 1L, bottomEdge);
      double yMinimum = _placement[movableCells[i]].y - Ydist;
      if (directedEdges->find(fromTo) == directedEdges->end()) {
        (*directedEdges)[fromTo] = make_pair(
            unlimited,
            -static_cast<long>(ceil((yMinimum - realBBox.yMin) / rowHeight)));
      } else {
        long oldCost = (*directedEdges)[fromTo].second;
        (*directedEdges)[fromTo].second = min(
            oldCost,
            -static_cast<long>(ceil((yMinimum - realBBox.yMin) / rowHeight)));
      }

      fromTo = make_pair(topEdge, i + 1L);
      double yMaximum = _placement[movableCells[i]].y + Ydist;
      if (directedEdges->find(fromTo) == directedEdges->end()) {
        (*directedEdges)[fromTo] = make_pair(
            unlimited,
            -static_cast<long>(ceil((realBBox.yMax - yMaximum) / rowHeight)));
      } else {
        long oldCost = (*directedEdges)[fromTo].second;
        (*directedEdges)[fromTo].second = min(
            oldCost,
            -static_cast<long>(ceil((realBBox.yMax - yMaximum) / rowHeight)));
      }
    }
  }

  delete cellIdToPos;

  long numFlowNodes = firstEdge + 2L * necessaryEdges->size() - 1L;
  long numFlowEdges = directedEdges->size();

  delete necessaryEdges;

  flowSolver solver;
  solver.initProblem(numFlowNodes, numFlowEdges, unlimited, sourceId, sinkId,
                     directedEdges, maxMem);

#ifndef YDEBUG
  delete directedEdges;
#endif

  buildProblemTime.stop();

  if (_params.verb.getForMajStats() > 0) {
    cout << "  Building the flow problem took "
         << buildProblemTime.getUserTime() << " seconds " << endl;
  }

  Timer flowTime;

  if (_params.verb.getForActions() > 0) {
    cout << "  Calling the flow solver with " << numFlowNodes << " nodes and "
         << numFlowEdges << " edges" << endl;
  }

  solver.solve();

  flowTime.stop();

  if (_params.verb.getForMajStats() > 0) {
    cout << "  The flow solver took " << flowTime.getUserTime() << " seconds"
         << endl;
  }

  price_t bottomPrice = solver.getPrice(bottomEdge);

#ifdef YDEBUG
  for (edgeHashMap::const_iterator i = directedEdges->begin();
       i != directedEdges->end(); ++i) {
    unsigned long from = i->first.first;
    unsigned long to = i->first.second;
    long cap = i->second.first;
    price_t cost = i->second.second;

    if (cap != unlimited) continue;

    price_t fromPrice = solver.getPrice(from);
    price_t toPrice = solver.getPrice(to);

    if (toPrice - fromPrice > cost) {
      cout << "constraint violated! " << endl
           << " from " << from << " to " << to << " cost " << cost << endl
           << " from price " << fromPrice << " to price " << toPrice << " diff "
           << toPrice - fromPrice << endl;
    }
  }

  delete directedEdges;
#endif

  for (unsigned i = 0; i < movableCells.size(); ++i) {
#ifdef YDEBUG
    price_t oldRowHeights = static_cast<price_t>(
        (_placement[movableCells[i]].y - realBBox.yMin) / rowHeight);
#endif
    double newYLoc =
        realBBox.yMin +
        static_cast<double>(solver.getPrice(i + 1) - bottomPrice) * rowHeight;
    removeFromCoreCols(_placement, getHGraph(), coreCols, core,
                       _coreRows[0].getSpacing(), movableCells[i]);
    if (getHGraph().isNodeMacro(movableCells[i])) {
      removeFromCoreRows(_placement, movableCells[i]);
      _placement[movableCells[i]].y = newYLoc;
      addToCoreRows(_placement, movableCells[i]);
    } else {
      setLocation(movableCells[i],
                  Point(_placement[movableCells[i]].x, newYLoc));
    }
    addToCoreCols(_placement, getHGraph(), coreCols, core,
                  _coreRows[0].getSpacing(), movableCells[i]);
#ifdef YDEBUG
    solver.setPrice(i + 1, bottomPrice + oldRowHeights);
#endif
  }

#ifdef YDEBUG
  unsigned numCellsNotplaced = _cellsNotInRows;
  double overlap = calcOverlapGeneric(All, true);
  if (numCellsNotPlaced > 0 || greaterThanDouble(overlap, 0.)) {
    cout << "There was a problem, reverting placement" << endl;
    for (unsigned i = 0; i < movableCells.size(); ++i) {
      double newYLoc =
          realBBox.yMin +
          static_cast<double>(solver.getPrice(i + 1) - bottomPrice) * rowHeight;
      removeFromCoreCols(_placement, getHGraph(), coreCols, core,
                         _coreRows[0].getSpacing(), movableCells[i]);
      if (getHGraph().isNodeMacro(movableCells[i])) {
        removeFromCoreRows(_placement, movableCells[i]);
        _placement[movableCells[i]].y = newYLoc;
        addToCoreRows(_placement, movableCells[i]);
      } else {
        setLocation(movableCells[i],
                    Point(_placement[movableCells[i]].x, newYLoc));
      }
      addToCoreCols(_placement, getHGraph(), coreCols, core,
                    _coreRows[0].getSpacing(), movableCells[i]);
    }
  }
#endif
}
#endif

#ifdef USEFLOW
void RBPlacement::doUnconstrainedYFlow(const PlacementWOrient &origPl,
                                       const HGraphWDimensions &hg,
                                       const vector<unsigned> &cellIds,
                                       const vector<BBox> &cellBBoxes,
                                       PlacementWOrient &changedPl,
                                       MaxMem *maxMem) const {
  if (_params.verb.getForActions() > 0) {
    cout << " Unconstrained Y direction" << endl;
  }

  Timer buildProblemTime;

  BBox core = getBBox(false);
  double rowHeight = _coreRows[0].getHeight();

  // find all of the nets that touch these cells
  vector<bool> movableCells(hg.getNumNodes(), false);
  BitBoard *edgesBoard = new BitBoard(hg.getNumEdges());
  for (unsigned i = 0; i < cellIds.size(); ++i) {
    movableCells[cellIds[i]] = true;
    const HGFNode &node = hg.getNodeByIdx(cellIds[i]);
    for (itHGFEdgeLocal edge = node.edgesBegin(); edge != node.edgesEnd();
         ++edge) {
      edgesBoard->setBit((*edge)->getIndex());
    }
  }

  const auto &affectedEdges = edgesBoard->getIndicesOfSetBits();

  // build a list of necessaryEdges,
  // those that will change HPWL based on the movement
  vector<unsigned> *necessaryEdges = new vector<unsigned>(0);
  for (unsigned i = 0; i < affectedEdges.size(); ++i) {
    double minY = DBL_MAX;
    double maxY = -DBL_MAX;
    const HGFEdge &edge = hg.getEdgeByIdx(affectedEdges[i]);
    unsigned offset = 0;
    for (itHGFNodeLocal node = edge.nodesBegin(); node != edge.nodesEnd();
         ++node, ++offset) {
      unsigned nodeId = (*node)->getIndex();
      // only concerned with fixed cells
      if (movableCells[nodeId]) continue;

      Point pinOffset = hg.nodesOnEdgePinOffset(offset, affectedEdges[i],
                                                origPl.getOrient(nodeId));
      Point pinPos = origPl[nodeId] + pinOffset;

      minY = min(minY, pinPos.y);
      maxY = max(maxY, pinPos.y);
    }

    // now we have the span of the fixed cells on this net
    // if this span contains the span of the bbox, this net will not change
    // given movement of cells in the bbox
    if (greaterThanDouble(minY, core.yMin) || lessThanDouble(maxY, core.yMax)) {
      necessaryEdges->push_back(affectedEdges[i]);
    }
  }

  delete edgesBoard;

  // make a map so we can go from cellId to its position in movableCells
  vector<unsigned> *cellIdToPos =
      new vector<unsigned>(hg.getNumNodes(), UINT_MAX - 1);
  for (unsigned i = 0; i < cellIds.size(); ++i) {
    (*cellIdToPos)[cellIds[i]] = i;
  }

  // now we have lists of the movable cells as well as
  // all the necessary edges, so we can make all the constraints

  unsigned long sourceId = cellIds.size() + 1;
  unsigned long sinkId = sourceId + 1;
  unsigned long bottomEdge = sinkId + 1;
  unsigned long topEdge = bottomEdge + 1;
  unsigned long firstEdge = topEdge + 1;
  long unlimited = necessaryEdges->size();

  // map (from,to) -> (cap,cost)
  edgeHashMap *directedEdges = new edgeHashMap();

  // add the edges constraining the left and right sides of the entire area
  (*directedEdges)[make_pair(bottomEdge, topEdge)] =
      make_pair(unlimited, static_cast<long>(core.getHeight() / rowHeight));
  (*directedEdges)[make_pair(topEdge, bottomEdge)] =
      make_pair(unlimited, -static_cast<long>(core.getHeight() / rowHeight));

  // add the edges attaching the left edges to the sink
  // and the right edges to the source
  for (unsigned long i = 0; i < necessaryEdges->size(); ++i) {
    (*directedEdges)[make_pair(firstEdge + 2L * i, sinkId)] = make_pair(1L, 0L);
    (*directedEdges)[make_pair(sourceId, firstEdge + 2L * i + 1L)] =
        make_pair(1L, 0L);
  }

  // add the two edges for every pin

  for (unsigned long i = 0; i < necessaryEdges->size(); ++i) {
    unsigned edgeId = (*necessaryEdges)[i];
    const HGFEdge &edge = hg.getEdgeByIdx(edgeId);
    unsigned offset = 0;
    unsigned fixedPins = 0;
    double fixedPinYMin = DBL_MAX;
    double fixedPinYMax = -DBL_MAX;
    double cellAreaOnEdge = 0.;
    for (itHGFNodeLocal n = edge.nodesBegin(); n != edge.nodesEnd();
         ++n, ++offset) {
      unsigned nodeId = (*n)->getIndex();
      cellAreaOnEdge += hg.getWeight(nodeId);
      Point pinOffset =
          hg.nodesOnEdgePinOffset(offset, edgeId, origPl.getOrient(nodeId));
      if (movableCells[nodeId]) {
        // edge constraining left side
        (*directedEdges)[make_pair((*cellIdToPos)[nodeId] + 1L,
                                   firstEdge + 2L * i)] =
            make_pair(unlimited,
                      static_cast<long>(floor(pinOffset.y / rowHeight)));
        // edge constraining right side
        (*directedEdges)[make_pair(firstEdge + 2L * i + 1L,
                                   (*cellIdToPos)[nodeId] + 1L)] =
            make_pair(unlimited,
                      -static_cast<long>(floor(pinOffset.y / rowHeight)));
      } else {
        // fixed pin, collect all and have at most 2 constraints
        ++fixedPins;
        double yPinLoc = origPl[nodeId].y + pinOffset.y;
        fixedPinYMin = min(fixedPinYMin, yPinLoc);
        fixedPinYMax = max(fixedPinYMax, yPinLoc);
      }
    }
    if (fixedPins > 0) {
      // edge constraining bottom side
      (*directedEdges)[make_pair(bottomEdge, firstEdge + 2L * i)] = make_pair(
          unlimited,
          static_cast<long>(floor((fixedPinYMin - core.yMin) / rowHeight)));
      // edge constraining top side
      (*directedEdges)[make_pair(firstEdge + 2L * i + 1L, topEdge)] = make_pair(
          unlimited,
          static_cast<long>(ceil((core.yMax - fixedPinYMax) / rowHeight)));
    }

    // add constraint for minimum net length
    double minLength = 0.5 * sqrt(cellAreaOnEdge);
    (*directedEdges)[make_pair(firstEdge + 2L * i + 1L, firstEdge + 2L * i)] =
        make_pair(unlimited, -static_cast<long>(floor(minLength / rowHeight)));
  }

  // constrain each cell to be within its given bbox
  for (unsigned long i = 0; i < cellIds.size(); ++i) {
    if (!cellBBoxes[i].isEmpty() &&
        greaterThanDouble(cellBBoxes[i].getArea(), 0.)) {
      pair<unsigned long, unsigned long> fromTo = make_pair(i + 1L, bottomEdge);
      if (directedEdges->find(fromTo) == directedEdges->end()) {
        (*directedEdges)[fromTo] = make_pair(
            unlimited, -static_cast<long>(
                           ceil((cellBBoxes[i].yMin - core.yMin) / rowHeight)));
      } else {
        long oldCost = (*directedEdges)[fromTo].second;
        (*directedEdges)[fromTo].second = min(
            oldCost, -static_cast<long>(
                         ceil((cellBBoxes[i].yMin - core.yMin) / rowHeight)));
      }

      fromTo = make_pair(topEdge, i + 1L);
      if (directedEdges->find(fromTo) == directedEdges->end()) {
        (*directedEdges)[fromTo] = make_pair(
            unlimited, -static_cast<long>(
                           ceil((core.yMax - cellBBoxes[i].yMax) / rowHeight)));
      } else {
        long oldCost = (*directedEdges)[fromTo].second;
        (*directedEdges)[fromTo].second = min(
            oldCost, -static_cast<long>(
                         ceil((core.yMax - cellBBoxes[i].yMax) / rowHeight)));
      }
    }
  }

  delete cellIdToPos;

  long numFlowNodes = firstEdge + 2L * necessaryEdges->size() - 1L;
  long numFlowEdges = directedEdges->size();

  delete necessaryEdges;

  flowSolver solver;

  solver.initProblem(numFlowNodes, numFlowEdges, unlimited, sourceId, sinkId,
                     directedEdges, maxMem);

  delete directedEdges;

  buildProblemTime.stop();

  if (_params.verb.getForMajStats() > 0) {
    cout << "  Building the flow problem took "
         << buildProblemTime.getUserTime() << " seconds " << endl;
  }

  Timer flowTime;

  if (_params.verb.getForActions() > 0) {
    cout << "  Calling the flow solver with " << numFlowNodes << " nodes and "
         << numFlowEdges << " edges" << endl;
  }

  solver.solve();

  flowTime.stop();

  if (_params.verb.getForMajStats() > 0) {
    cout << "  The flow solver took " << flowTime.getUserTime() << " seconds"
         << endl;
  }

  price_t bottomPrice = solver.getPrice(bottomEdge);

  for (unsigned i = 0; i < cellIds.size(); ++i) {
    double newYLoc =
        core.yMin +
        static_cast<double>(solver.getPrice(i + 1) - bottomPrice) * rowHeight;

    changedPl[cellIds[i]].y = newYLoc;
  }
}
#endif

#ifdef USEFLOW
void RBPlacement::doXFlowInBBox(ISPD06DensityMap *map,
                                vector<vector<unsigned> > &coreCols,
                                const BBox &theBox,
                                vector<unsigned> *allowedMovables, double Xdist,
                                MaxMem *maxMem) {
  if (_params.verb.getForActions() > 0) {
    cout << " X direction" << endl;
  }

  Timer buildProblemTime;

  BBox core = getBBox(false);
  BBox realBBox = core.intersect(theBox);

  double rowSpacing = _coreRows[0].getSpacing();

  if (greaterThanDouble(Xdist, 0.) && lessThanDouble(Xdist, rowSpacing)) {
    cout << "  maximum X distance (" << Xdist << ") too small" << endl;
    return;
  }

  // align the bbox to row and site boundaries
  realBBox.xMin =
      core.xMin + ceil((realBBox.xMin - core.xMin) / rowSpacing) * rowSpacing;
  realBBox.xMax =
      core.xMin + floor((realBBox.xMax - core.xMin) / rowSpacing) * rowSpacing;
  realBBox.yMin =
      core.yMin + ceil((realBBox.yMin - core.yMin) / _coreRows[0].getHeight()) *
                      _coreRows[0].getHeight();
  realBBox.yMax = core.yMin + floor((realBBox.yMax - core.yMin) /
                                    _coreRows[0].getHeight()) *
                                  _coreRows[0].getHeight();

  vector<bool> allowed(getHGraph().getNumNodes(), false);
  if (allowedMovables != NULL) {
    for (unsigned i = 0; i < allowedMovables->size(); ++i) {
      allowed[(*allowedMovables)[i]] = true;
    }
  }

  // pick out all the movable cells within the "realBBox"
  BitBoard movableCellsBoard(getHGraph().getNumNodes());
  for (unsigned i = 0; i < _coreRows.size(); ++i) {
    if (lessThanDouble(_coreRows[i].getYCoord(), realBBox.yMin)) continue;
    if (greaterThanDouble(_coreRows[i].getYCoord(), realBBox.yMax)) break;

    for (unsigned j = 0; j < _coreRows[i].getAllCells().size(); ++j) {
      unsigned cellId = _coreRows[i].getAllCells()[j];

      if (allowedMovables != NULL && !allowed[cellId]) continue;

      if (lessThanDouble(_placement[cellId].x, realBBox.xMin)) continue;
      if (greaterThanDouble(_placement[cellId].x, realBBox.xMax)) break;

      if (!_isFixed[cellId]) movableCellsBoard.setBit(cellId);
    }
  }

  const auto &movableCells = movableCellsBoard.getIndicesOfSetBits();

  if (movableCells.size() == 0) return;

  bool usePins = true;
  double origScaledHPWL = evalHPWL(usePins);
  if (map != NULL) {
    map->compute();
    double penalty = map->getScaledOverflow();
    origScaledHPWL *= (1. + 0.01 * penalty);
  }

  // now that we have all the movable cells, find all of the nets that touch
  // them
  BitBoard *edgesBoard = new BitBoard(getHGraph().getNumEdges());
  for (unsigned i = 0; i < movableCells.size(); ++i) {
    const HGFNode &node = getHGraph().getNodeByIdx(movableCells[i]);
    for (itHGFEdgeLocal edge = node.edgesBegin(); edge != node.edgesEnd();
         ++edge) {
      edgesBoard->setBit((*edge)->getIndex());
    }
  }

  const auto &affectedEdges = edgesBoard->getIndicesOfSetBits();

  // build a list of necessaryEdges,
  // those that will change HPWL based on the movement
  vector<unsigned> *necessaryEdges = new vector<unsigned>(0);
  for (unsigned i = 0; i < affectedEdges.size(); ++i) {
    double minX = DBL_MAX;
    double maxX = -DBL_MAX;
    const HGFEdge &edge = getHGraph().getEdgeByIdx(affectedEdges[i]);
    unsigned offset = 0;
    for (itHGFNodeLocal node = edge.nodesBegin(); node != edge.nodesEnd();
         ++node, ++offset) {
      unsigned nodeId = (*node)->getIndex();
      // only concerned with fixed cells
      if (movableCellsBoard.isBitSet(nodeId)) continue;

      Point pinOffset = getHGraph().nodesOnEdgePinOffset(
          offset, affectedEdges[i], _placement.getOrient(nodeId));
      Point pinPos = _placement[nodeId] + pinOffset;

      minX = min(minX, pinPos.x);
      maxX = max(maxX, pinPos.x);
    }

    // now we have the span of the fixed cells on this net
    // if this span contains the span of the bbox, this net will not change
    // given movement of cells in the bbox
    if (greaterThanDouble(minX, realBBox.xMin) ||
        lessThanDouble(maxX, realBBox.xMax)) {
      necessaryEdges->push_back(affectedEdges[i]);
    }
  }

  delete edgesBoard;

  // make a map so we can go from cellId to its position in movableCells
  vector<unsigned> *cellIdToPos =
      new vector<unsigned>(getHGraph().getNumNodes(), UINT_MAX - 1);
  for (unsigned i = 0; i < movableCells.size(); ++i) {
    (*cellIdToPos)[movableCells[i]] = i;
  }

  // now we have lists of the movable cells as well as
  // all the necessary edges, so we can make all the constraints

  unsigned long sourceId = movableCells.size() + 1;
  unsigned long sinkId = sourceId + 1;
  unsigned long leftEdge = sinkId + 1;
  unsigned long rightEdge = leftEdge + 1;
  unsigned long firstEdge = rightEdge + 1;
  long unlimited = necessaryEdges->size();

  // map (from,to) -> (cap,cost)
  edgeHashMap *directedEdges = new edgeHashMap();

  // add the edges constraining the left and right sides of the entire area
  (*directedEdges)[make_pair(leftEdge, rightEdge)] =
      make_pair(unlimited, static_cast<long>(realBBox.getWidth() / rowSpacing));
  (*directedEdges)[make_pair(rightEdge, leftEdge)] = make_pair(
      unlimited, -static_cast<long>(realBBox.getWidth() / rowSpacing));

  // add the edges attaching the left edges to the sink
  // and the right edges to the source
  for (unsigned long i = 0; i < necessaryEdges->size(); ++i) {
    (*directedEdges)[make_pair(firstEdge + 2L * i, sinkId)] = make_pair(1L, 0L);
    (*directedEdges)[make_pair(sourceId, firstEdge + 2L * i + 1L)] =
        make_pair(1L, 0L);
  }

  // add the two edges for every pin
  for (unsigned long i = 0; i < necessaryEdges->size(); ++i) {
    unsigned edgeId = (*necessaryEdges)[i];
    const HGFEdge &edge = getHGraph().getEdgeByIdx(edgeId);
    unsigned offset = 0;
    unsigned fixedPins = 0;
    double fixedPinXMin = DBL_MAX;
    double fixedPinXMax = -DBL_MAX;
    for (itHGFNodeLocal n = edge.nodesBegin(); n != edge.nodesEnd();
         ++n, ++offset) {
      unsigned nodeId = (*n)->getIndex();
      Point pinOffset = getHGraph().nodesOnEdgePinOffset(
          offset, edgeId, _placement.getOrient(nodeId));
      if (movableCellsBoard.isBitSet(nodeId)) {
        // edge constraining left side
        (*directedEdges)[make_pair((*cellIdToPos)[nodeId] + 1L,
                                   firstEdge + 2L * i)] =
            make_pair(unlimited,
                      static_cast<long>(floor(pinOffset.x / rowSpacing)));
        // edge constraining right side
        (*directedEdges)[make_pair(firstEdge + 2L * i + 1L,
                                   (*cellIdToPos)[nodeId] + 1L)] =
            make_pair(unlimited,
                      -static_cast<long>(floor(pinOffset.x / rowSpacing)));
      } else {
        // fixed pin, collect all and have at most 2 constraints
        ++fixedPins;
        double xPinLoc = _placement[nodeId].x + pinOffset.x;
        fixedPinXMin = min(fixedPinXMin, xPinLoc);
        fixedPinXMax = max(fixedPinXMax, xPinLoc);
      }
    }
    if (fixedPins > 0) {
      // edge constraining left side
      (*directedEdges)[make_pair(leftEdge, firstEdge + 2L * i)] = make_pair(
          unlimited, static_cast<long>(
                         floor((fixedPinXMin - realBBox.xMin) / rowSpacing)));
      // edge constraining right side
      (*directedEdges)[make_pair(firstEdge + 2L * i + 1L, rightEdge)] =
          make_pair(unlimited,
                    static_cast<long>(
                        ceil((realBBox.xMax - fixedPinXMax) / rowSpacing)));
    }
  }

  // for each cell, constrain it between the cells in front and back of it
#ifdef XDEBUG
  bool bad = false;
#endif
  for (unsigned i = 0; i < _coreRows.size(); ++i) {
    double coreRowYMin = _coreRows[i].getYCoord();
    double coreRowYMax = coreRowYMin + _coreRows[i].getHeight();

    for (unsigned j = 0; j < _coreRows[i].getAllCells().size(); ++j) {
      unsigned cellId = _coreRows[i].getAllCells()[j];

#ifdef XDEBUG
      if (j != 0) {
        unsigned prev = _coreRows[i].getAllCells()[j - 1];
        if (_placement[prev].x > _placement[cellId].x) {
          cout << "inconsistency in the all cells structure" << endl;
          cout << "prev " << prev << "  curr " << cellId << endl;
          cout << "prev: is macro? " << getHGraph().isNodeMacro(prev) << endl;
          cout << "prev: is fixed? " << _isFixed[prev] << endl;
          cout << "prev: placement " << _placement[prev] << endl;
          cout << "curr: is macro? " << getHGraph().isNodeMacro(cellId) << endl;
          cout << "curr: is fixed? " << _isFixed[cellId] << endl;
          cout << "curr: placement " << _placement[cellId] << endl;
          bad = true;
        }
      }
#endif

      if (!movableCellsBoard.isBitSet(cellId)) continue;

      unsigned angle = _placement.getOrient(cellId).getAngle();
      bool notRotated = angle == 0 || angle == 180;
      double cellWidth = notRotated ? getHGraph().getNodeWidth(cellId)
                                    : getHGraph().getNodeHeight(cellId);
      cellWidth = ceil(cellWidth / rowSpacing) * rowSpacing;
      double cellHeight = notRotated ? getHGraph().getNodeHeight(cellId)
                                     : getHGraph().getNodeWidth(cellId);
      double cellYMin = _placement[cellId].y;
      double cellYMax = cellYMin + cellHeight;

      if (!getHGraph().isNodeMacro(cellId)) {
        const RBCellCoord &cellCoord = _cellCoords[cellId];

        if (cellCoord.cellOffset == 0) {
          // bound on the left by the subrow boundary
          double xMinimum = max(
              realBBox.xMin,
              _coreRows[cellCoord.coreRowId][cellCoord.subRowId].getXStart());
          (*directedEdges)[make_pair((*cellIdToPos)[cellId] + 1L, leftEdge)] =
              make_pair(unlimited,
                        -static_cast<long>(
                            ceil((xMinimum - realBBox.xMin) / rowSpacing)));
        }

        if ((cellCoord.cellOffset + 1) ==
            _coreRows[cellCoord.coreRowId][cellCoord.subRowId].getNumCells()) {
          // the cell is the last in its subrow
          // bound on the right by the subrow boundary
          double xMaximum =
              min(realBBox.xMax,
                  _coreRows[cellCoord.coreRowId][cellCoord.subRowId].getXEnd() -
                      cellWidth);
          (*directedEdges)[make_pair(rightEdge, (*cellIdToPos)[cellId] + 1L)] =
              make_pair(unlimited,
                        -static_cast<long>(
                            ceil((realBBox.xMax - xMaximum) / rowSpacing)));
        }
      }

      // macros don't fit into the row structures perfectly, so we need to check
      // if relationships that look like constraints are real first
      // before adding them to the problem

      unsigned k = j;
      // when a macro is being bound by another macro in this core row,
      // it could possibly be bound by 2 of them (no more, I think)
      bool boundedBySomething = false;
      double boundedYMin = DBL_MAX;
      double boundedYMax = -DBL_MAX;
      while (true) {
        if (k == 0) {
          if (!boundedBySomething) {
            // if we have gotten here, the only thing constraining
            // cellId to the left is the left side of the core
            // bound on the left by the left edge of the chip
            double xMinimum = max(realBBox.xMin, core.xMin);
            pair<unsigned long, unsigned long> fromTo =
                make_pair((*cellIdToPos)[cellId] + 1L, leftEdge);
            if (directedEdges->find(fromTo) == directedEdges->end()) {
              (*directedEdges)[fromTo] = make_pair(
                  unlimited, -static_cast<long>(ceil(
                                 (xMinimum - realBBox.xMin) / rowSpacing)));
            } else {
              long oldCost = (*directedEdges)[fromTo].second;
              (*directedEdges)[fromTo].second = min(
                  oldCost, -static_cast<long>(
                               ceil((xMinimum - realBBox.xMin) / rowSpacing)));
            }
          }
          break;
        } else {
          // see if the cell that is at _coreRows[i].getAllCells()[k-1] could
          // block cellId
          unsigned prevCell = _coreRows[i].getAllCells()[k - 1];
          unsigned prevAngle = _placement.getOrient(prevCell).getAngle();
          bool prevNotRotated = prevAngle == 0 || prevAngle == 180;
          double prevWidth = prevNotRotated
                                 ? getHGraph().getNodeWidth(prevCell)
                                 : getHGraph().getNodeHeight(prevCell);
          prevWidth = ceil(prevWidth / rowSpacing) * rowSpacing;
          double prevHeight = prevNotRotated
                                  ? getHGraph().getNodeHeight(prevCell)
                                  : getHGraph().getNodeWidth(prevCell);
          double prevYMin = _placement[prevCell].y;
          double prevYMax = prevYMin + prevHeight;

          if (greaterOrEqualDouble(prevYMin, cellYMax) ||
              greaterOrEqualDouble(cellYMin, prevYMax) ||
              (allowedMovables != NULL && !allowed[prevCell])) {
            // these two cells could never overlap, try the cell before this one
            --k;
          } else {
            // this is a true constraint
            if (greaterOrEqualDouble(max(coreRowYMin, prevYMin), boundedYMin) &&
                lessOrEqualDouble(min(coreRowYMax, prevYMax), boundedYMax)) {
              // no need to bound against this object,
              // as it is obscured by the object(s) we already bounded against
              --k;
            } else {
              if (movableCellsBoard.isBitSet(prevCell)) {
                // bound on the left by the movable object
                pair<unsigned long, unsigned long> fromTo = make_pair(
                    (*cellIdToPos)[cellId] + 1L, (*cellIdToPos)[prevCell] + 1L);
                if (directedEdges->find(fromTo) == directedEdges->end()) {
                  (*directedEdges)[fromTo] = make_pair(
                      unlimited,
                      -static_cast<long>(ceil(prevWidth / rowSpacing)));
                } else {
                  long oldCost = (*directedEdges)[fromTo].second;
                  (*directedEdges)[fromTo].second =
                      min(oldCost,
                          -static_cast<long>(ceil(prevWidth / rowSpacing)));
                }
                boundedBySomething = true;
                boundedYMin = max(coreRowYMin, min(boundedYMin, prevYMin));
                boundedYMax = min(coreRowYMax, max(boundedYMax, prevYMax));
              } else if (lessThanDouble(_placement[prevCell].x,
                                        _placement[cellId].x)) {
                // bound on the left by the fixed object
                double xMinimum =
                    max(realBBox.xMin, _placement[prevCell].x + prevWidth);
                xMinimum = max(xMinimum, core.xMin);
                pair<unsigned long, unsigned long> fromTo =
                    make_pair((*cellIdToPos)[cellId] + 1L, leftEdge);
                if (directedEdges->find(fromTo) == directedEdges->end()) {
                  (*directedEdges)[fromTo] = make_pair(
                      unlimited, -static_cast<long>(ceil(
                                     (xMinimum - realBBox.xMin) / rowSpacing)));
                } else {
                  long oldCost = (*directedEdges)[fromTo].second;
                  (*directedEdges)[fromTo].second = min(
                      oldCost, -static_cast<long>(ceil(
                                   (xMinimum - realBBox.xMin) / rowSpacing)));
                }
                boundedBySomething = true;
                boundedYMin = max(coreRowYMin, min(boundedYMin, prevYMin));
                boundedYMax = min(coreRowYMax, max(boundedYMax, prevYMax));
              }
              --k;
              if (lessOrEqualDouble(boundedYMin, coreRowYMin) &&
                  greaterOrEqualDouble(boundedYMax, coreRowYMax)) {
                break;
              }
            }
          }
        }
      }

      k = j;
      boundedBySomething = false;
      boundedYMin = DBL_MAX;
      boundedYMax = -DBL_MAX;
      while (true) {
        if ((k + 1) == _coreRows[i].getAllCells().size()) {
          if (!boundedBySomething) {
            // if we have gotten here, the only thing constraining
            // cellId to the right is the right side of the core
            // bound on the right by the right edge of the chip
            double xMaximum = min(realBBox.xMax, core.xMax - cellWidth);
            pair<unsigned long, unsigned long> fromTo =
                make_pair(rightEdge, (*cellIdToPos)[cellId] + 1L);
            if (directedEdges->find(fromTo) == directedEdges->end()) {
              (*directedEdges)[fromTo] = make_pair(
                  unlimited, -static_cast<long>(ceil(
                                 (realBBox.xMax - xMaximum) / rowSpacing)));
            } else {
              long oldCost = (*directedEdges)[fromTo].second;
              (*directedEdges)[fromTo].second = min(
                  oldCost, -static_cast<long>(
                               ceil((realBBox.xMax - xMaximum) / rowSpacing)));
            }
          }
          break;
        } else {
          // see if the cell that is at _coreRows[i].getAllCells()[k+1] could
          // block cellId
          unsigned nextCell = _coreRows[i].getAllCells()[k + 1];
          unsigned nextAngle = _placement.getOrient(nextCell).getAngle();
          bool nextNotRotated = nextAngle == 0 || nextAngle == 180;
          double nextHeight = nextNotRotated
                                  ? getHGraph().getNodeHeight(nextCell)
                                  : getHGraph().getNodeWidth(nextCell);
          double nextYMin = _placement[nextCell].y;
          double nextYMax = nextYMin + nextHeight;

          if (greaterOrEqualDouble(nextYMin, cellYMax) ||
              greaterOrEqualDouble(cellYMin, nextYMax) ||
              (allowedMovables != NULL && !allowed[nextCell])) {
            // these two cells could never overlap, try the cell after this one
            ++k;
          } else {
            // this is a true constraint
            if (greaterOrEqualDouble(max(coreRowYMin, nextYMin), boundedYMin) &&
                lessOrEqualDouble(min(coreRowYMax, nextYMax), boundedYMax)) {
              // no need to bound against this object,
              // as it is obscured by the object(s) we already bounded against
              ++k;
            } else {
              if (movableCellsBoard.isBitSet(nextCell)) {
                // bound on the right by the movable
                pair<unsigned long, unsigned long> fromTo = make_pair(
                    (*cellIdToPos)[nextCell] + 1L, (*cellIdToPos)[cellId] + 1L);
                if (directedEdges->find(fromTo) == directedEdges->end()) {
                  (*directedEdges)[fromTo] = make_pair(
                      unlimited,
                      -static_cast<long>(ceil(cellWidth / rowSpacing)));
                } else {
                  long oldCost = (*directedEdges)[fromTo].second;
                  (*directedEdges)[fromTo].second =
                      min(oldCost,
                          -static_cast<long>(ceil(cellWidth / rowSpacing)));
                }
                boundedBySomething = true;
                boundedYMin = max(coreRowYMin, min(boundedYMin, nextYMin));
                boundedYMax = min(coreRowYMax, max(boundedYMax, nextYMax));
              } else if (greaterThanDouble(_placement[nextCell].x,
                                           _placement[cellId].x)) {
                // bound on the right by the fixed object
                double xMaximum =
                    min(realBBox.xMax, _placement[nextCell].x - cellWidth);
                xMaximum = min(xMaximum, core.xMax - cellWidth);
                pair<unsigned long, unsigned long> fromTo =
                    make_pair(rightEdge, (*cellIdToPos)[cellId] + 1L);
                if (directedEdges->find(fromTo) == directedEdges->end()) {
                  (*directedEdges)[fromTo] = make_pair(
                      unlimited, -static_cast<long>(ceil(
                                     (realBBox.xMax - xMaximum) / rowSpacing)));
                } else {
                  long oldCost = (*directedEdges)[fromTo].second;
                  (*directedEdges)[fromTo].second = min(
                      oldCost, -static_cast<long>(ceil(
                                   (realBBox.xMax - xMaximum) / rowSpacing)));
                }
                boundedBySomething = true;
                boundedYMin = max(coreRowYMin, min(boundedYMin, nextYMin));
                boundedYMax = min(coreRowYMax, max(boundedYMax, nextYMax));
              }
              ++k;
              if (lessOrEqualDouble(boundedYMin, coreRowYMin) &&
                  greaterOrEqualDouble(boundedYMax, coreRowYMax)) {
                break;
              }
            }
          }
        }
      }
    }
  }

#ifdef XDEBUG
  abkfatal(!bad, "badness occurred");
#endif

  /*
     if(map != NULL)
     {
       // add additional constraints for each cell if we have the grid available
       for(unsigned long i = 0; i < movableCells.size(); ++i)
       {
         BBox grid = map->getGridCell(_placement[movableCells[i]]);

         pair< unsigned long, unsigned long > fromTo = make_pair(i + 1L,
     leftEdge); double xMinimum = grid.xMin; if(directedEdges->find(fromTo) ==
     directedEdges->end())
         {
           (*directedEdges)[fromTo] = make_pair(unlimited,
     -static_cast<long>(ceil((xMinimum - realBBox.xMin)/rowSpacing)));
         }
         else
         {
           long oldCost = (*directedEdges)[fromTo].second;
           (*directedEdges)[fromTo].second = min(oldCost,
     -static_cast<long>(ceil((xMinimum - realBBox.xMin)/rowSpacing)));
         }

         fromTo = make_pair(rightEdge, i + 1L);
         double xMaximum = grid.xMax;
         if(directedEdges->find(fromTo) == directedEdges->end())
         {
           (*directedEdges)[fromTo] = make_pair(unlimited,
     -static_cast<long>(ceil((realBBox.xMax - xMaximum)/rowSpacing)));
         }
         else
         {
           long oldCost = (*directedEdges)[fromTo].second;
           (*directedEdges)[fromTo].second = min(oldCost,
     -static_cast<long>(ceil((realBBox.xMax - xMaximum)/rowSpacing)));
         }
       }
     }
  */

  if (greaterThanDouble(Xdist, 0.)) {
    for (unsigned long i = 0; i < movableCells.size(); ++i) {
      pair<unsigned long, unsigned long> fromTo = make_pair(i + 1L, leftEdge);
      double xMinimum = _placement[movableCells[i]].x - Xdist;
      if (directedEdges->find(fromTo) == directedEdges->end()) {
        (*directedEdges)[fromTo] = make_pair(
            unlimited,
            -static_cast<long>(ceil((xMinimum - realBBox.xMin) / rowSpacing)));
      } else {
        long oldCost = (*directedEdges)[fromTo].second;
        (*directedEdges)[fromTo].second = min(
            oldCost,
            -static_cast<long>(ceil((xMinimum - realBBox.xMin) / rowSpacing)));
      }

      fromTo = make_pair(rightEdge, i + 1L);
      double xMaximum = _placement[movableCells[i]].x + Xdist;
      if (directedEdges->find(fromTo) == directedEdges->end()) {
        (*directedEdges)[fromTo] = make_pair(
            unlimited,
            -static_cast<long>(ceil((realBBox.xMax - xMaximum) / rowSpacing)));
      } else {
        long oldCost = (*directedEdges)[fromTo].second;
        (*directedEdges)[fromTo].second = min(
            oldCost,
            -static_cast<long>(ceil((realBBox.xMax - xMaximum) / rowSpacing)));
      }
    }
  }

  delete cellIdToPos;

  // build the flow problem into their data structures

  long numFlowNodes = firstEdge + 2L * necessaryEdges->size() - 1L;
  long numFlowEdges = directedEdges->size();

  delete necessaryEdges;

  flowSolver solver;
  solver.initProblem(numFlowNodes, numFlowEdges, unlimited, sourceId, sinkId,
                     directedEdges, maxMem);

#ifndef XDEBUG
  delete directedEdges;
#endif

  buildProblemTime.stop();

  if (_params.verb.getForMajStats() > 0) {
    cout << "  Building the flow problem took "
         << buildProblemTime.getUserTime() << " seconds " << endl;
  }

  Timer flowTime;

  if (_params.verb.getForActions() > 0) {
    cout << "  Calling the flow solver with " << numFlowNodes << " nodes and "
         << numFlowEdges << " edges" << endl;
  }

  solver.solve();

  flowTime.stop();

  if (_params.verb.getForMajStats() > 0) {
    cout << "  The flow solver took " << flowTime.getUserTime() << " seconds"
         << endl;
  }

  price_t leftPrice = solver.getPrice(leftEdge);

#ifdef XDEBUG
  for (edgeHashMap::const_iterator i = directedEdges->begin();
       i != directedEdges->end(); ++i) {
    unsigned long from = i->first.first;
    unsigned long to = i->first.second;
    long cap = i->second.first;
    price_t cost = i->second.second;

    if (cap != unlimited) continue;

    price_t fromPrice = solver.getPrice(from);
    price_t toPrice = solver.getPrice(to);

    if (toPrice - fromPrice > cost) {
      cout << "constraint violated! " << endl
           << " from " << from << " to " << to << " cost " << cost << endl
           << " from price " << fromPrice << " to price " << toPrice << " diff "
           << toPrice - fromPrice << endl;
    }
  }

  delete directedEdges;
#endif

  for (unsigned i = 0; i < movableCells.size(); ++i) {
    price_t oldRowSpaces = static_cast<price_t>(
        (_placement[movableCells[i]].x - realBBox.xMin) / rowSpacing);
    double newXLoc =
        realBBox.xMin +
        static_cast<double>(solver.getPrice(i + 1) - leftPrice) * rowSpacing;
    removeFromCoreCols(_placement, getHGraph(), coreCols, core, rowSpacing,
                       movableCells[i]);
    if (getHGraph().isNodeMacro(movableCells[i])) {
      removeFromCoreRows(_placement, movableCells[i]);
      _placement[movableCells[i]].x = newXLoc;
      addToCoreRows(_placement, movableCells[i]);
    } else {
      setLocation(movableCells[i],
                  Point(newXLoc, _placement[movableCells[i]].y));
    }
    addToCoreCols(_placement, getHGraph(), coreCols, core, rowSpacing,
                  movableCells[i]);
    solver.setPrice(i + 1, leftPrice + oldRowSpaces);
  }

  double finalScaledHPWL = evalHPWL(usePins);
  if (map != NULL) {
    map->compute();
    double penalty = map->getScaledOverflow();
    finalScaledHPWL *= (1. + 0.01 * penalty);
  }

#ifdef XDEBUG
  unsigned numCellsNotPlaced = _cellsNotInRows;
  double overlap = calcOverlapGeneric(All, true);
  if (numCellsNotPlaced > 0 || greaterThanDouble(overlap, 0.))
#endif
    if (lessThanDouble(origScaledHPWL, finalScaledHPWL)) {
      cout << "HPWL increased, reverting placement" << endl;
      for (unsigned i = 0; i < movableCells.size(); ++i) {
        double newXLoc =
            realBBox.xMin +
            static_cast<double>(solver.getPrice(i + 1) - leftPrice) *
                rowSpacing;
        removeFromCoreCols(_placement, getHGraph(), coreCols, core, rowSpacing,
                           movableCells[i]);
        if (getHGraph().isNodeMacro(movableCells[i])) {
          removeFromCoreRows(_placement, movableCells[i]);
          _placement[movableCells[i]].x = newXLoc;
          addToCoreRows(_placement, movableCells[i]);
        } else {
          setLocation(movableCells[i],
                      Point(newXLoc, _placement[movableCells[i]].y));
        }
        addToCoreCols(_placement, getHGraph(), coreCols, core, rowSpacing,
                      movableCells[i]);
      }
    }
}
#endif

#ifdef USEFLOW
void RBPlacement::doUnconstrainedXFlow(const PlacementWOrient &origPl,
                                       const HGraphWDimensions &hg,
                                       const vector<unsigned> &cellIds,
                                       const vector<BBox> &cellBBoxes,
                                       PlacementWOrient &changedPl,
                                       MaxMem *maxMem) const {
  if (_params.verb.getForActions() > 0) {
    cout << " Unconstrained X direction" << endl;
  }

  Timer buildProblemTime;

  BBox core = getBBox(false);
  double rowSpacing = _coreRows[0].getSpacing();

  // find all of the nets that touch these cells
  vector<bool> movableCells(hg.getNumNodes(), false);
  BitBoard *edgesBoard = new BitBoard(hg.getNumEdges());
  for (unsigned i = 0; i < cellIds.size(); ++i) {
    movableCells[cellIds[i]] = true;
    const HGFNode &node = hg.getNodeByIdx(cellIds[i]);
    for (itHGFEdgeLocal edge = node.edgesBegin(); edge != node.edgesEnd();
         ++edge) {
      edgesBoard->setBit((*edge)->getIndex());
    }
  }

  const auto &affectedEdges = edgesBoard->getIndicesOfSetBits();

  // build a list of necessaryEdges,
  // those that will change HPWL based on the movement
  vector<unsigned> *necessaryEdges = new vector<unsigned>(0);
  for (unsigned i = 0; i < affectedEdges.size(); ++i) {
    double minX = DBL_MAX;
    double maxX = -DBL_MAX;
    const HGFEdge &edge = hg.getEdgeByIdx(affectedEdges[i]);
    unsigned offset = 0;
    for (itHGFNodeLocal node = edge.nodesBegin(); node != edge.nodesEnd();
         ++node, ++offset) {
      unsigned nodeId = (*node)->getIndex();
      // only concerned with fixed cells
      if (movableCells[nodeId]) continue;

      Point pinOffset = hg.nodesOnEdgePinOffset(offset, affectedEdges[i],
                                                origPl.getOrient(nodeId));
      Point pinPos = origPl[nodeId] + pinOffset;

      minX = min(minX, pinPos.x);
      maxX = max(maxX, pinPos.x);
    }

    // now we have the span of the fixed cells on this net
    // if this span contains the span of the bbox, this net will not change
    // given movement of cells in the bbox
    if (greaterThanDouble(minX, core.xMin) || lessThanDouble(maxX, core.xMax)) {
      necessaryEdges->push_back(affectedEdges[i]);
    }
  }

  delete edgesBoard;

  // make a map so we can go from cellId to its position in movableCells
  vector<unsigned> *cellIdToPos =
      new vector<unsigned>(hg.getNumNodes(), UINT_MAX - 1);
  for (unsigned i = 0; i < cellIds.size(); ++i) {
    (*cellIdToPos)[cellIds[i]] = i;
  }

  // now we have lists of the movable cells as well as
  // all the necessary edges, so we can make all the constraints

  unsigned long sourceId = cellIds.size() + 1;
  unsigned long sinkId = sourceId + 1;
  unsigned long leftEdge = sinkId + 1;
  unsigned long rightEdge = leftEdge + 1;
  unsigned long firstEdge = rightEdge + 1;
  long unlimited = necessaryEdges->size();

  // map (from,to) -> (cap,cost)
  edgeHashMap *directedEdges = new edgeHashMap();

  // add the edges constraining the left and right sides of the entire area
  (*directedEdges)[make_pair(leftEdge, rightEdge)] =
      make_pair(unlimited, static_cast<long>(core.getWidth() / rowSpacing));
  (*directedEdges)[make_pair(rightEdge, leftEdge)] =
      make_pair(unlimited, -static_cast<long>(core.getWidth() / rowSpacing));

  // add the edges attaching the left edges to the sink
  // and the right edges to the source
  for (unsigned long i = 0; i < necessaryEdges->size(); ++i) {
    (*directedEdges)[make_pair(firstEdge + 2L * i, sinkId)] = make_pair(1L, 0L);
    (*directedEdges)[make_pair(sourceId, firstEdge + 2L * i + 1L)] =
        make_pair(1L, 0L);
  }

  // add the two edges for every pin
  for (unsigned long i = 0; i < necessaryEdges->size(); ++i) {
    unsigned edgeId = (*necessaryEdges)[i];
    const HGFEdge &edge = hg.getEdgeByIdx(edgeId);
    unsigned offset = 0;
    unsigned fixedPins = 0;
    double fixedPinXMin = DBL_MAX;
    double fixedPinXMax = -DBL_MAX;
    double cellAreaOnEdge = 0.;
    for (itHGFNodeLocal n = edge.nodesBegin(); n != edge.nodesEnd();
         ++n, ++offset) {
      unsigned nodeId = (*n)->getIndex();
      cellAreaOnEdge += hg.getWeight(nodeId);
      Point pinOffset =
          hg.nodesOnEdgePinOffset(offset, edgeId, origPl.getOrient(nodeId));
      if (movableCells[nodeId]) {
        // edge constraining left side
        (*directedEdges)[make_pair((*cellIdToPos)[nodeId] + 1L,
                                   firstEdge + 2L * i)] =
            make_pair(unlimited,
                      static_cast<long>(floor(pinOffset.x / rowSpacing)));
        // edge constraining right side
        (*directedEdges)[make_pair(firstEdge + 2L * i + 1L,
                                   (*cellIdToPos)[nodeId] + 1L)] =
            make_pair(unlimited,
                      -static_cast<long>(floor(pinOffset.x / rowSpacing)));
      } else {
        // fixed pin, collect all and have at most 2 constraints
        ++fixedPins;
        double xPinLoc = origPl[nodeId].x + pinOffset.x;
        fixedPinXMin = min(fixedPinXMin, xPinLoc);
        fixedPinXMax = max(fixedPinXMax, xPinLoc);
      }
    }
    if (fixedPins > 0) {
      // edge constraining left side
      (*directedEdges)[make_pair(leftEdge, firstEdge + 2L * i)] = make_pair(
          unlimited,
          static_cast<long>(floor((fixedPinXMin - core.xMin) / rowSpacing)));
      // edge constraining right side
      (*directedEdges)[make_pair(firstEdge + 2L * i + 1L, rightEdge)] =
          make_pair(
              unlimited,
              static_cast<long>(ceil((core.xMax - fixedPinXMax) / rowSpacing)));
    }

    // add constraint for minimum net length
    double minLength = 0.5 * sqrt(cellAreaOnEdge);
    (*directedEdges)[make_pair(firstEdge + 2L * i + 1L, firstEdge + 2L * i)] =
        make_pair(unlimited, -static_cast<long>(floor(minLength / rowSpacing)));
  }

  // constrain each cell to be within its given bbox
  for (unsigned long i = 0; i < cellIds.size(); ++i) {
    if (!cellBBoxes[i].isEmpty() &&
        greaterThanDouble(cellBBoxes[i].getArea(), 0.)) {
      pair<unsigned long, unsigned long> fromTo = make_pair(i + 1L, leftEdge);
      if (directedEdges->find(fromTo) == directedEdges->end()) {
        (*directedEdges)[fromTo] = make_pair(
            unlimited, -static_cast<long>(ceil(
                           (cellBBoxes[i].xMin - core.xMin) / rowSpacing)));
      } else {
        long oldCost = (*directedEdges)[fromTo].second;
        (*directedEdges)[fromTo].second = min(
            oldCost, -static_cast<long>(
                         ceil((cellBBoxes[i].xMin - core.xMin) / rowSpacing)));
      }

      fromTo = make_pair(rightEdge, i + 1L);
      if (directedEdges->find(fromTo) == directedEdges->end()) {
        (*directedEdges)[fromTo] = make_pair(
            unlimited, -static_cast<long>(ceil(
                           (core.xMax - cellBBoxes[i].xMax) / rowSpacing)));
      } else {
        long oldCost = (*directedEdges)[fromTo].second;
        (*directedEdges)[fromTo].second = min(
            oldCost, -static_cast<long>(
                         ceil((core.xMax - cellBBoxes[i].xMax) / rowSpacing)));
      }
    }
  }

  delete cellIdToPos;

  long numFlowNodes = firstEdge + 2L * necessaryEdges->size() - 1L;
  long numFlowEdges = directedEdges->size();

  delete necessaryEdges;

  flowSolver solver;

  solver.initProblem(numFlowNodes, numFlowEdges, unlimited, sourceId, sinkId,
                     directedEdges, maxMem);

  delete directedEdges;

  buildProblemTime.stop();

  if (_params.verb.getForMajStats() > 0) {
    cout << "  Building the flow problem took "
         << buildProblemTime.getUserTime() << " seconds " << endl;
  }

  Timer flowTime;

  if (_params.verb.getForActions() > 0) {
    cout << "  Calling the flow solver with " << numFlowNodes << " nodes and "
         << numFlowEdges << " edges" << endl;
  }

  solver.solve();

  flowTime.stop();

  if (_params.verb.getForMajStats() > 0) {
    cout << "  The flow solver took " << flowTime.getUserTime() << " seconds"
         << endl;
  }

  price_t leftPrice = solver.getPrice(leftEdge);

  for (unsigned i = 0; i < cellIds.size(); ++i) {
    double newXLoc =
        core.xMin +
        static_cast<double>(solver.getPrice(i + 1) - leftPrice) * rowSpacing;

    changedPl[cellIds[i]].x = newXLoc;
  }
}
#endif

#ifdef USEFLOW
void RBPlacement::doFlow(ISPD06DensityMap *map, bool Xflow, bool Yflow,
                         MaxMem *maxMem) {
  double maxXDist = -1., maxYDist = -1.;
  vector<unsigned> *allowedMovables = NULL;
  doFlowHelper(map, Xflow, maxXDist, Yflow, maxYDist, _placement,
               allowedMovables, maxMem);
}

void RBPlacement::doFlowHelper(ISPD06DensityMap *map, bool Xflow, double Xdist,
                               bool Yflow, double Ydist, PlacementWOrient &pl,
                               vector<unsigned> *allowedMovables,
                               MaxMem *maxMem) {
  if (!Xflow && !Yflow) return;

  BBox coreArea = getBBox(false);

  unsigned numCells = getHGraph().getNumNodes();

  unsigned sitesWide = static_cast<unsigned>(
      ceil(coreArea.getWidth() / _coreRows[0].getSpacing()));

  vector<vector<unsigned> > coreCols(sitesWide, vector<unsigned>(0));

  vector<bool> allowed(getHGraph().getNumNodes(), false);
  if (allowedMovables != NULL) {
    for (unsigned i = 0; i < allowedMovables->size(); ++i) {
      allowed[(*allowedMovables)[i]] = true;
    }
  }

  reinstateCoreRows();
  resetPlacementWOri(pl);
  for (unsigned i = 0; i < getHGraph().getNumNodes(); ++i) {
    if (!_isFixed[i] && getHGraph().isNodeMacro(i)) {
      addToCoreRows(_placement, i);
    }

    unsigned angle = pl.getOrient(i).getAngle();
    bool notRotated = angle == 0 || angle == 180;
    double cellWidth =
        notRotated ? getHGraph().getNodeWidth(i) : getHGraph().getNodeHeight(i);
    double cellHeight =
        notRotated ? getHGraph().getNodeHeight(i) : getHGraph().getNodeWidth(i);

    if (equalDouble(cellWidth, 0.) || equalDouble(cellHeight, 0.)) continue;

    double minX = _placement[i].x;
    // change me to unsigned
    int minIndex = static_cast<int>(floor(minX - coreArea.xMin) /
                                    _coreRows[0].getSpacing());
    int spacesWide =
        static_cast<unsigned>(ceil(cellWidth / _coreRows[0].getSpacing()));
    int maxIndex = minIndex + spacesWide;
    int maxIndex2 = static_cast<int>(floor(minX + cellWidth - coreArea.xMin) /
                                     _coreRows[0].getSpacing());
    if (maxIndex <= maxIndex2) {
      ++maxIndex;
    }

    for (int j = minIndex; j < maxIndex; ++j) {
      if (j >= 0 && static_cast<unsigned>(j) < coreCols.size()) {
        coreCols[j].push_back(i);
      }
    }
  }

  for (unsigned i = 0; i < coreCols.size(); ++i) {
    sort(coreCols[i].begin(), coreCols[i].end(), CompareYWOri(_placement));
  }

  unsigned tiles =
      static_cast<unsigned>(ceil(sqrt(static_cast<double>(numCells) / 50000.)));

  if (_params.verb.getForActions() > 0) {
    cout << " Using a " << tiles << "x" << tiles << " grid" << endl;
  }

  double width = coreArea.getWidth() / static_cast<double>(tiles);
  double height = coreArea.getHeight() / static_cast<double>(tiles);

  BBox tile(coreArea.xMin, coreArea.yMin, coreArea.xMin + width,
            coreArea.yMin + height);

  for (unsigned j = 0; j < tiles; ++j) {
    tile.xMin = coreArea.xMin;
    tile.xMax = tile.xMin + width;
    for (unsigned i = 0; i < tiles; ++i) {
      if (_params.verb.getForActions() > 0) {
        cout << " (" << i << "," << j << ")" << endl;
      }
      if (Xflow)
        doXFlowInBBox(map, coreCols, tile, allowedMovables, Xdist, maxMem);
      if (Yflow)
        doYFlowInBBox(map, coreCols, tile, allowedMovables, Ydist, maxMem);
      tile.xMin += width;
      tile.xMax += width;
    }
    tile.yMin += height;
    tile.yMax += height;
  }

  reinstateCoreRows();
  // remove sites below those macros that are allowed
  if (allowedMovables == NULL) {
    removeSitesBelowMacros();
  } else {
    for (unsigned i = 0; i < allowedMovables->size(); ++i) {
      if (!getHGraph().isNodeMacro((*allowedMovables)[i]) ||
          _isFixed[(*allowedMovables)[i]])
        continue;

      splitCoreRowsByCell(_placement, (*allowedMovables)[i]);
    }
    resetPlacement();
  }

  if (allowedMovables != NULL) {
    for (unsigned i = 0; i < allowedMovables->size(); ++i) {
      if (_isFixed[(*allowedMovables)[i]]) continue;

      pl[(*allowedMovables)[i]] = _placement[(*allowedMovables)[i]];
    }
  } else {
    pl = _placement;
  }
}
#endif
