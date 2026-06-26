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
#include "CongestionMaps/ISPD06DensityMap.h"
#include "RBPlacement.h"

using std::cout;
using std::endl;
using std::lower_bound;
using std::make_pair;
using std::max;
using std::min;
using std::pair;
using std::sort;
using std::upper_bound;

/*
bool lessThanImprove(const pair<double,double> &a, const pair<double,double> &b)
{
  // if wirelength doesn't get worse with either, priority is cong then WL
  if(greaterOrEqualDouble(a.second,0.) && greaterOrEqualDouble(b.second,0.))
  {
    if(lessThanDouble(a.first, b.first)) return true;
    else if(equalDouble(a.first, b.first)) return lessThanDouble(a.second,
b.second); else return false;
  }
  // else if neither cong gets worse, priority is on WL then cong
  else if(greaterOrEqualDouble(a.first,0.) && greaterOrEqualDouble(b.first,0.))
  {
    if(lessThanDouble(a.second, b.second)) return true;
    else if(equalDouble(a.second, b.second)) return lessThanDouble(a.first,
b.first); else return false;
  }
  // else priority is on cong first then WL
  else
  {
    if(lessThanDouble(a.first, b.first)) return true;
    else if(equalDouble(a.first, b.first)) return lessThanDouble(a.second,
b.second); else return false;
  }
}
*/

bool lessThanImprove(const pair<double, double> &a,
                     const pair<double, double> &b) {
  //  cout << "comparing (" << a.first << "," << a.second << ") with (" <<
  //  b.first << "," << b.second << ")" << endl;

  //  return (a.first + a.second) < (b.first + b.second);
  return (a < b);
}

double uglyHackHPWLDiv = 1.;
double uglyHackCongDiv = 1.;
double uglyHackHPWLCut = -0.00001;
// double uglyHackHPWLCut = 0;

pair<double, double> RBPlacement::calculateGreedyMoveImprove(
    BaseGeneric2DResourceMap *theMap, unsigned cell, const Point &newPos) {
  pair<double, double> improve;

  bool skipLargeNets = false;
  improve.second = -CalcSwapCost(cell, newPos, skipLargeNets) / uglyHackHPWLDiv;

  if (theMap == NULL) {
    improve.first = 0.;
  } else {
    //    if(greaterOrEqualDouble(improve.second, 0.))
    //    {
    improve.first =
        -theMap->estimateOverfullnessChangeFromCellMove(cell, newPos) /
        uglyHackCongDiv;
    //    }
    //    else
    //    {
    //      improve.first = -DBL_MAX;
    //    }
  }

  if (lessThanDouble(improve.second, uglyHackHPWLCut)) {
    improve.first = -DBL_MAX;
    improve.second = -DBL_MAX;
  }

  return improve;
}

void RBPlacement::findBestGreedyMove(BaseGeneric2DResourceMap *theMap,
                                     double cutoffRadius, const BBox &layout,
                                     unsigned cell, Point &newPos,
                                     pair<double, double> &improve) {
  newPos = _placement[cell];
  //  improve = make_pair(0.,0.);

  // calculate the x and y coords to find the Weber region
  std::vector<double> xcoords, ycoords;

  const HGFNode &node = _hgWDims->getNodeByIdx(cell);
  const Orient &cellOrient = _placement.getOrient(cell);

  unsigned eOffset = 0;
  for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd();
       ++e, ++eOffset) {
    Point cellOffset =
        _hgWDims->edgesOnNodePinOffset(eOffset, cell, cellOrient);

    unsigned nOffset = 0;
    for (itHGFNodeLocal n = (*e)->nodesBegin(); n != (*e)->nodesEnd();
         ++n, ++nOffset) {
      unsigned neighborCell = (*n)->getIndex();
      if (neighborCell == cell) continue;

      Point neighborOffset = _hgWDims->nodesOnEdgePinOffset(
          nOffset, (*e)->getIndex(), _placement.getOrient(neighborCell));

      xcoords.push_back(_placement[neighborCell].x + neighborOffset.x -
                        cellOffset.x);
      ycoords.push_back(_placement[neighborCell].y + neighborOffset.y -
                        cellOffset.y);
    }
  }

  BBox allowableBox(
      _placement[cell].x - cutoffRadius, _placement[cell].y - cutoffRadius,
      _placement[cell].x + cutoffRadius, _placement[cell].y + cutoffRadius);

  if (xcoords.size() < 2) {
    xcoords.push_back(allowableBox.xMin);
    xcoords.push_back(allowableBox.xMax);
    ycoords.push_back(allowableBox.yMin);
    ycoords.push_back(allowableBox.yMax);
  }

  sort(xcoords.begin(), xcoords.end());
  sort(ycoords.begin(), ycoords.end());

  unsigned minIdx = 0, maxIdx = 0;

  maxIdx = xcoords.size() / 2;
  if (xcoords.size() % 2 == 0) {
    minIdx = maxIdx - 1;
  } else {
    minIdx = maxIdx;
  }

  // binary search for the box that includes the cells current pos, so we don't
  // look in too much area
  unsigned xdistfromcenter = 0;
  pair<std::vector<double>::const_iterator,
       std::vector<double>::const_iterator>
      xpair = equal_range(xcoords.begin(), xcoords.end(), _placement[cell].x);
  if (xpair.first == xcoords.end() || xpair.second == xcoords.begin()) {
    xdistfromcenter = minIdx;
  } else {
    // index of biggest element <= target
    unsigned lowest = xpair.first - xcoords.begin();
    // index of smallest element >= target
    unsigned highest = (xpair.first != xpair.second)
                           ? xpair.second - xcoords.begin() - 1
                           : xpair.second - xcoords.begin();
    xdistfromcenter = max((minIdx >= lowest) ? minIdx - lowest : 0,
                          (highest >= maxIdx) ? highest - maxIdx : 0);
  }

  unsigned ydistfromcenter = 0;
  pair<std::vector<double>::const_iterator,
       std::vector<double>::const_iterator>
      ypair = equal_range(ycoords.begin(), ycoords.end(), _placement[cell].y);
  if (ypair.first == ycoords.end() || ypair.second == ycoords.begin()) {
    ydistfromcenter = minIdx;
  } else {
    // index of biggest element <= target
    unsigned lowest = ypair.first - ycoords.begin();
    // index of smallest element >= target
    unsigned highest = (ypair.first != ypair.second)
                           ? ypair.second - ycoords.begin() - 1
                           : ypair.second - ycoords.begin();
    ydistfromcenter = max((minIdx >= lowest) ? minIdx - lowest : 0,
                          (highest >= maxIdx) ? highest - maxIdx : 0);
  }

  // the necessary increase in idx so as to cover the smallest BBox that
  // contains _placement[cell]
  unsigned necessaryDist = max(xdistfromcenter, ydistfromcenter);
  unsigned biggestIdx = maxIdx + necessaryDist;

  // want to target 5 boxes at most
  unsigned increment = max(1u, necessaryDist / 5);

  BBox curBox, prevBox;
  Point weberCenter;

  while (maxIdx <= biggestIdx) {
    curBox.xMin = xcoords[minIdx];
    curBox.xMax = xcoords[maxIdx];
    curBox.yMin = ycoords[minIdx];
    curBox.yMax = ycoords[maxIdx];

    if (curBox != prevBox) {
      if (prevBox.isEmpty()) {
        // this is the first time, so we have the weber region
        weberCenter = curBox.getGeomCenter();

        BBox searchBox = allowableBox.intersect(curBox);

        if (!searchBox.isEmpty())
          findBestGreedyMoveInBox(theMap, cell, newPos, improve, searchBox,
                                  weberCenter);
      } else if (equalDouble(prevBox.xMin, prevBox.xMax) ||
                 equalDouble(prevBox.yMin, prevBox.yMax)) {
        // if the last box had a zero dimension (or two)
        // its more efficient to search the whole box instead of 4 separate
        // since they will cover the exact same area
        BBox searchBox = allowableBox.intersect(curBox);

        if (!searchBox.isEmpty())
          findBestGreedyMoveInBox(theMap, cell, newPos, improve, searchBox,
                                  weberCenter);
      } else {
        // curBox is an expansion of prevBox, so we search in the
        // 4 boxes in curBox not in prevBox

        BBox above;
        above.xMin = curBox.xMin;
        above.xMax = curBox.xMax;
        above.yMin = prevBox.yMax;
        above.yMax = curBox.xMax;

        if (!equalDouble(above.yMin, above.yMax)) {
          BBox searchBox = allowableBox.intersect(above);

          if (!searchBox.isEmpty())
            findBestGreedyMoveInBox(theMap, cell, newPos, improve, searchBox,
                                    weberCenter);
        }

        BBox below;
        below.xMin = curBox.xMin;
        below.xMax = curBox.xMax;
        below.yMin = curBox.yMin;
        below.yMax = prevBox.yMin;

        if (!equalDouble(below.yMin, below.yMax)) {
          BBox searchBox = allowableBox.intersect(below);

          if (!searchBox.isEmpty())
            findBestGreedyMoveInBox(theMap, cell, newPos, improve, searchBox,
                                    weberCenter);
        }

        BBox left;
        left.xMin = curBox.xMin;
        left.xMax = prevBox.xMin;
        left.yMin = prevBox.yMin;
        left.yMax = prevBox.yMax;

        if (!equalDouble(left.xMin, left.xMax)) {
          BBox searchBox = allowableBox.intersect(left);

          if (!searchBox.isEmpty())
            findBestGreedyMoveInBox(theMap, cell, newPos, improve, searchBox,
                                    weberCenter);
        }

        BBox right;
        right.xMin = prevBox.xMax;
        right.xMax = curBox.xMax;
        right.yMin = prevBox.yMin;
        right.yMax = prevBox.yMax;

        if (!equalDouble(right.xMin, right.xMax)) {
          BBox searchBox = allowableBox.intersect(right);

          if (!searchBox.isEmpty())
            findBestGreedyMoveInBox(theMap, cell, newPos, improve, searchBox,
                                    weberCenter);
        }
      }

      // take the first improvement we get
      if (newPos != _placement[cell]) break;

      prevBox = curBox;
    }

    if (minIdx > 0 && increment > minIdx) {
      minIdx = 0;
      maxIdx = xcoords.size() - 1;
    } else {
      minIdx -= increment;
      maxIdx += increment;
    }
  }
}

void RBPlacement::findBestGreedyMoveInBox(BaseGeneric2DResourceMap *theMap,
                                          unsigned cell, Point &newPos,
                                          pair<double, double> &improve,
                                          const BBox &box,
                                          const Point &weberCenter) {
  double width = _hgWDims->getNodeWidth(cell);
  double weberDist =
      fabs(newPos.x - weberCenter.x) + fabs(newPos.y - weberCenter.y);

  RBPCoreRow cr;
  cr._y = box.yMin;
  std::vector<RBPCoreRow>::const_iterator lowerY =
      lower_bound(_coreRows.begin(), _coreRows.end(), cr, compareCoreRowY());
  unsigned smallestCoreRowIdx = lowerY - _coreRows.begin();

  cr._y = box.yMax;
  std::vector<RBPCoreRow>::const_iterator upperY =
      upper_bound(_coreRows.begin(), _coreRows.end(), cr, compareCoreRowY());
  unsigned largestCoreRowIdx = upperY - _coreRows.begin();

  Point newPoint;
  for (unsigned i = smallestCoreRowIdx; i < largestCoreRowIdx; ++i) {
    newPoint.y = _coreRows[i].getYCoord();

    double spacing = _coreRows[i].getSpacing();
    double invSpacing = 1.0 / spacing;

    const std::vector<pair<double, double> > &wsLocs =
        _coreRows[i].getWhitespaceLocs(getHGraph());

    std::vector<pair<double, double> >::const_iterator lowest = lower_bound(
        wsLocs.begin(), wsLocs.end(), make_pair(box.xMin, box.xMin));
    unsigned lowestIdx = lowest - wsLocs.begin();
    if (lowest != wsLocs.begin() &&
        greaterOrEqualDouble(box.xMin, wsLocs[lowestIdx - 1].first) &&
        greaterOrEqualDouble(wsLocs[lowestIdx - 1].second, box.xMin)) {
      --lowestIdx;
    }

    std::vector<pair<double, double> >::const_iterator highest = upper_bound(
        wsLocs.begin(), wsLocs.end(), make_pair(box.xMax, box.xMax));
    unsigned highestIdx = highest - wsLocs.begin();
    if (highest != wsLocs.end() &&
        greaterOrEqualDouble(box.xMax, wsLocs[highestIdx].first) &&
        greaterOrEqualDouble(wsLocs[highestIdx].second, box.xMax)) {
      ++highestIdx;
    }

    double bestWLThisRow = -DBL_MAX;

    for (unsigned w = lowestIdx; w < highestIdx; ++w) {
      newPoint.x = wsLocs[w].first;
      newPoint.x = max(newPoint.x, box.xMin);
      newPoint.x = min(newPoint.x, box.xMax);
      newPoint.x =
          _coreRows[i].getXStart() +
          ceil((newPoint.x - _coreRows[i].getXStart()) * invSpacing) * spacing;

      double wsWidth = wsLocs[w].second - newPoint.x;

      if (greaterOrEqualDouble(wsWidth, width)) {
        double otherX = wsLocs[w].second - width;
        otherX = max(otherX, box.xMin);
        otherX = min(otherX, box.xMax);
        otherX =
            _coreRows[i].getXStart() +
            floor((otherX - _coreRows[i].getXStart()) * invSpacing) * spacing;

        double otherWsWidth = otherX - wsLocs[w].first;
        if (greaterOrEqualDouble(otherWsWidth, width) &&
            lessThanDouble(fabs(otherX - weberCenter.x),
                           fabs(newPoint.x - weberCenter.x))) {
          newPoint.x = otherX;
        }

        pair<double, double> newImprove =
            calculateGreedyMoveImprove(theMap, cell, newPoint);

        if (lessThanImprove(improve, newImprove)) {
          improve = newImprove;
          newPos = newPoint;
          weberDist =
              fabs(newPos.x - weberCenter.x) + fabs(newPos.y - weberCenter.y);
        } else if (improve == newImprove &&
                   lessThanDouble(fabs(newPoint.x - weberCenter.x) +
                                      fabs(newPoint.y - weberCenter.y),
                                  weberDist)) {
          improve = newImprove;
          newPos = newPoint;
          weberDist =
              fabs(newPos.x - weberCenter.x) + fabs(newPos.y - weberCenter.y);
        }
        if (theMap == NULL) {
          if (lessOrEqualDouble(newImprove.second, bestWLThisRow)) {
            bestWLThisRow = newImprove.second;
          } else {
            break;
          }
        }
      }
    }  // for each ws chunk in the row
  }
}

void RBPlacement::greedyMovement(greedyMoveMode optimizeMode,
                                 BaseGeneric2DResourceMap *theMap) {
  if (_params.verb.getForActions())
    cout << " Running greedy movement ..." << endl;

  Timer tm;

  uglyHackHPWLCut = -0.00001;

  double startWL = evalHPWL(true);
  uglyHackHPWLDiv = startWL;
  if (_params.verb.getForMajStats()) {
    cout << "  HPWL before greedy movement is " << startWL << endl;
  }

  std::vector<Point> newPos(getNumCells());
  std::vector<pair<double, double> > improvement(
      getNumCells(), make_pair(-DBL_MAX, -DBL_MAX));
  unsigned bestCell = UINT_MAX;
  pair<double, double> bestImprove = make_pair(-DBL_MAX, -DBL_MAX);

  std::vector<unsigned> allcells;
  allcells.reserve(getNumCells());

  std::vector<bool> moved(getNumCells(), false);

  BBox layout = getBBox(false);

  const unsigned numNodes = getNumCells() - _hgWDims->getNumTerminals();
  double distCutoff = 0.5 * max(layout.getWidth(), layout.getHeight());
  if (numNodes > _params.greedyMoveCells) {
    distCutoff = 0.5 * sqrt(static_cast<double>(_params.greedyMoveCells) *
                            layout.getArea() / static_cast<double>(numNodes));
  }

  double origOverfill = 0.;
  if (optimizeMode != HPWL) {
    unsigned xTiles = theMap->getNumHorizTiles();
    unsigned yTiles = theMap->getNumVertTiles();
    origOverfill = theMap->getTotalOverfill();
    uglyHackCongDiv = origOverfill;
    cout << "  Starting total overfill is " << origOverfill << endl;
    cout << "  Starting max overfill is " << theMap->getMaxOverfill() << endl;
    cout << "  Starting overfull tiles is " << theMap->getNumOverfullTiles()
         << " = "
         << 100. * static_cast<double>(theMap->getNumOverfullTiles()) /
                static_cast<double>(xTiles * yTiles)
         << "%" << endl;
    //    if(optimizeMode == ISPDContest)
    //    {
    //      double penalty =
    //      static_cast<ISPD06DensityMap*>(theMap)->getScaledOverflow();
    //      if(lessThanDouble(penalty,0.01))
    //      {
    //        if(_params.verb.getForActions() > 0)
    //        {
    //          cout << " Scaled overflow penalty is small, skipping
    //          optimization" << endl;
    //        }
    //        return;
    //      }
    //    }
  }

  for (unsigned i = 0; i < getNumCells(); ++i) {
    if (_params.verb.getForMajStats() > 3) {
      if (i % 100 == 0 || (i + 1) == getNumCells()) {
        cout << "Finding best greedy move for cell " << i << endl;
      }
    }

    if (!isCoreCell(i) || isFixed(i) || !_isInSubRow[i] ||
        _hgWDims->isNodeMacro(i))
      continue;

    allcells.push_back(i);

    //    if(optimizeMode == HPWL)
    //    {
    //      improvement[i] = make_pair(0.,0.);
    //    }
    //    else
    //    {
    //      improvement[i] = make_pair(0.,-0.005*startWL);
    //    }
    improvement[i] = make_pair(-DBL_MAX, -DBL_MAX);

    findBestGreedyMove(theMap, distCutoff, layout, i, newPos[i],
                       improvement[i]);

    if (newPos[i] != _placement[i] &&
        lessThanImprove(bestImprove, improvement[i])) {
      bestImprove = improvement[i];
      bestCell = i;
    }
  }

  std::vector<std::vector<unsigned> > affectedCells(getNumCells());

  for (unsigned i = 0; i < allcells.size(); ++i) {
    fillInAffectedCells(allcells[i], affectedCells[allcells[i]], allcells);
  }

  double improveCongCutoff = 0.;
  double improveWLCutoff = 0.;
  double maxHPWL = DBL_MAX;
  if (optimizeMode != HPWL) {
    maxHPWL = 1.005 * startWL;  // don't allow to make HPWL too much worse, in
                                // hopes to make congestion much better

    if (_params.verb.getForMajStats() > 3) {
      cout << "  Improvement Overfill cutoff " << improveCongCutoff << endl;
    }
  } else {
    if (1.e-6 * startWL > improveWLCutoff) improveWLCutoff = 1.e-6;
    if (_params.verb.getForMajStats() > 3) {
      cout << "  Improvement cutoff " << improveWLCutoff << endl;
    }
  }

  unsigned moves = 0;

  bool startedSecondPhase = (optimizeMode == HPWL);

  while (bestCell != UINT_MAX) {
    unsigned movedCell = bestCell;

    if (optimizeMode == Cong) {
      std::vector<unsigned> theCell(1, movedCell);
      static_cast<ISPD04CongMap *>(theMap)->removeDemandForAttachedNets(
          getHGraph(), theCell);
    } else if (optimizeMode == ISPDContest) {
      static_cast<ISPD06DensityMap *>(theMap)->removeCellUsage(movedCell);
    }

    Point origPos = _placement[movedCell];
    setLocation(movedCell, newPos[movedCell]);

    if (optimizeMode == Cong) {
      std::vector<unsigned> theCell(1, movedCell);
      static_cast<ISPD04CongMap *>(theMap)->addDemandForAttachedNets(
          getHGraph(), theCell);
    } else if (optimizeMode == ISPDContest) {
      static_cast<ISPD06DensityMap *>(theMap)->addCellUsage(movedCell);
    }

    ++moves;

    moved[movedCell] = true;

    if (_params.verb.getForMajStats() > 3) {
      cout << "  Just made move number " << moves << endl;
      cout << "     Cell id " << movedCell << endl;
      cout << "     From " << origPos << endl;
      cout << "     To   " << _placement[movedCell] << endl;
      if (optimizeMode != HPWL)
        cout << "     Overfill Improvement " << bestImprove.first << endl;
      cout << "     WL Improvement " << bestImprove.second << endl;
    }

    double currWL = evalHPWL(true);

    if (startedSecondPhase) {
      if (optimizeMode == ISPDContest) {
        //        double penalty =
        //        static_cast<ISPD06DensityMap*>(theMap)->getScaledOverflow();
        //        if(lessThanDouble(penalty,0.01))
        //        {
        //          break;
        //        }
        if (lessOrEqualDouble(bestImprove.first, improveCongCutoff)) break;
      } else if (optimizeMode == Cong) {
        if (lessOrEqualDouble(bestImprove.first, improveCongCutoff)) break;
      } else  // optimizeMode == HPWL
      {
        if (lessOrEqualDouble(bestImprove.second, improveWLCutoff)) break;
      }
    }

    bestImprove = make_pair(-DBL_MAX, -DBL_MAX);
    bestCell = UINT_MAX;
    unsigned majUpdates = 0, minUpdates = 0;
    for (unsigned i = 0; i < allcells.size(); ++i) {
      if (!moved[allcells[i]] &&
          // could allcells[i] candidate position overlap with movedCells new
          // position ?
          newPos[allcells[i]].y == _placement[movedCell].y &&
          ((lessOrEqualDouble(newPos[allcells[i]].x, _placement[movedCell].x) &&
            lessThanDouble(
                _placement[movedCell].x,
                newPos[allcells[i]].x + _hgWDims->getNodeWidth(allcells[i]))) ||
           (lessOrEqualDouble(_placement[movedCell].x, newPos[allcells[i]].x) &&
            lessThanDouble(newPos[allcells[i]].x,
                           _placement[movedCell].x +
                               _hgWDims->getNodeWidth(movedCell))))) {
        //        if(optimizeMode == HPWL)
        //        {
        //          improvement[allcells[i]] = make_pair(0.,0.);
        //        }
        //        else
        //        {
        //          improvement[allcells[i]] = make_pair(0.,currWL - maxHPWL);
        //        }
        improvement[allcells[i]] = make_pair(-DBL_MAX, -DBL_MAX);

        findBestGreedyMove(theMap, distCutoff, layout, allcells[i],
                           newPos[allcells[i]], improvement[allcells[i]]);
        ++majUpdates;
      }
      // is movedCell attached to allcells[i] ?
      else if (!moved[allcells[i]] &&
               binary_search(affectedCells[movedCell].begin(),
                             affectedCells[movedCell].end(), allcells[i])) {
        pair<double, double> oldImprove = improvement[allcells[i]];
        improvement[allcells[i]] = calculateGreedyMoveImprove(
            theMap, allcells[i], newPos[allcells[i]]);
        if ((greaterThanDouble(oldImprove.first, 0.) &&
             lessOrEqualDouble(improvement[allcells[i]].first, 0.)) ||
            (greaterThanDouble(oldImprove.second, 0.) &&
             lessOrEqualDouble(improvement[allcells[i]].second, 0.))) {
          //          if(optimizeMode == HPWL)
          //          {
          //            improvement[allcells[i]] = make_pair(0.,0.);
          //          }
          //          else
          //          {
          //            improvement[allcells[i]] = make_pair(0.,currWL -
          //            maxHPWL);
          //          }
          improvement[allcells[i]] = make_pair(-DBL_MAX, -DBL_MAX);

          findBestGreedyMove(theMap, distCutoff, layout, allcells[i],
                             newPos[allcells[i]], improvement[allcells[i]]);
          ++majUpdates;
        } else {
          ++minUpdates;
        }
      }

      if (!moved[allcells[i]] &&
          newPos[allcells[i]] != _placement[allcells[i]] &&
          lessThanImprove(bestImprove, improvement[allcells[i]])) {
        if (lessOrEqualDouble(currWL - improvement[allcells[i]].second,
                              maxHPWL)) {
          bestImprove = improvement[allcells[i]];
          bestCell = allcells[i];
        }
      }
    }

    if (_params.verb.getForMajStats() > 3) {
      cout << "     Major Updates needed " << majUpdates << endl;
      cout << "     Minor Updates needed " << minUpdates << endl;
    }

    if (lessOrEqualDouble(bestImprove.first, 0.) &&
        lessOrEqualDouble(bestImprove.second, 0.)) {
      bestImprove = make_pair(-DBL_MAX, -DBL_MAX);
      bestCell = UINT_MAX;
    }

    if (bestCell == UINT_MAX && !startedSecondPhase) {
      startedSecondPhase = true;
      uglyHackHPWLCut = -DBL_MAX;
      bestImprove = make_pair(-DBL_MAX, -DBL_MAX);

      if (_params.verb.getForMajStats() > 3) {
        cout << "    Starting the 2nd phase, more greedy" << endl;
      }

      for (unsigned i = 0; i < allcells.size(); ++i) {
        if (_params.verb.getForMajStats() > 3) {
          if (i % 100 == 0 || (i + 1) == allcells.size()) {
            cout << "Finding best greedy move for cell " << i << endl;
          }
        }

        moved[allcells[i]] = false;

        improvement[allcells[i]] = make_pair(-DBL_MAX, -DBL_MAX);

        findBestGreedyMove(theMap, distCutoff, layout, allcells[i],
                           newPos[allcells[i]], improvement[allcells[i]]);

        if (newPos[allcells[i]] != _placement[allcells[i]] &&
            lessThanImprove(bestImprove, improvement[allcells[i]])) {
          if (lessOrEqualDouble(currWL - improvement[allcells[i]].second,
                                maxHPWL)) {
            bestImprove = improvement[allcells[i]];
            bestCell = allcells[i];
          }
        }
      }
    }

    if (lessOrEqualDouble(bestImprove.first, 0.) &&
        lessOrEqualDouble(bestImprove.second, 0.)) {
      bestImprove = make_pair(-DBL_MAX, -DBL_MAX);
      bestCell = UINT_MAX;
    }
  }

  double endWL = evalHPWL(true);

  if (optimizeMode != HPWL) {
    unsigned xTiles = theMap->getNumHorizTiles();
    unsigned yTiles = theMap->getNumVertTiles();
    cout << "  Ending total overfill is " << theMap->getTotalOverfill() << endl;
    cout << "  Ending max overfill is " << theMap->getMaxOverfill() << endl;
    cout << "  Ending overfull tiles is " << theMap->getNumOverfullTiles()
         << " = "
         << 100. * static_cast<double>(theMap->getNumOverfullTiles()) /
                static_cast<double>(xTiles * yTiles)
         << "%" << endl;
  }

  tm.stop();

  if (_params.verb.getForMajStats()) {
    cout << " Greedy movement took " << tm << endl;
    cout << "  Moves performed: " << moves << endl;
    cout << "  HPWL after greedy movement is " << endWL << endl;
    double impr = (startWL == 0.) ? 0. : (startWL - endWL) / startWL;
    if (impr < 0.)
      cout << "  % increase in HPWL is " << floor(-impr * 100000. + 0.5) / 1000.
           << endl;
    else
      cout << "  % decrease in HPWL is " << floor(impr * 100000. + 0.5) / 1000.
           << endl;
  }
}
