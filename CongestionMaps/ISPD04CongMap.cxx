/**************************************************************************
***
*** Copyright (c) 2026 IonQ, Inc.
*** Copyright (c) 1995-2000 Regents of the University of California,
***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
*** Copyright (c) 2000-2006 Regents of the University of Michigan,
***               Saurabh N. Adya, Jarrod A. Roy, David A. Papa and
***               Igor L. Markov
***
***  Contact author(s): abk@cs.ucsd.edu, imarkov@umich.edu
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

#include "ISPD04CongMap.h"

#include <algorithm>
#include <fstream>
#include <vector>

#include "ABKCommon/abkassert.h"
#include "GeomTrees/FastSteiner/mst2.h"
#include "Placement/placeWOri.h"
#ifdef USEFLUTE
#include "Flute/flute.h"
#endif

using std::cout;
using std::endl;
using std::make_pair;
using std::max;
using std::min;
using std::ofstream;
using std::pair;
using std::sort;
using std::string;
using std::unique;
using std::vector;

ISPD04CongMap::ISPD04CongMap(const RBPlacement &rbplace,
                             const PlacementWOrient &pl, unsigned numHorizTiles,
                             unsigned numVertTiles, bool useSteiner)
    : _netbitboard(rbplace.getHGraph().getNumEdges()),
      _alpha(0.6),
      _pl(pl),
      _useSteiner(useSteiner) {
  abkfatal(numHorizTiles > 0, "0 horizontal tiles make no sense");
  abkfatal(numVertTiles > 0, "0 vertical tiles make no sense");

  setNumTiles(numHorizTiles, numVertTiles);
  BBox layoutBBox = rbplace.getBBox(false);
  setLayoutBBox(layoutBBox);
  setTileSize(layoutBBox.getWidth() / static_cast<double>(numHorizTiles),
              layoutBBox.getHeight() / static_cast<double>(numVertTiles));

  _tempMap = vector<vector<ResourceType> >(
      numHorizTiles, vector<ResourceType>(numVertTiles, ResourceZeroType()()));

  buildSupplyMap(rbplace.routingResources);

#ifdef USEFLUTE
  if (_useSteiner) readLUT();
#endif
}

void ISPD04CongMap::resetDemand(void) {
  for (unsigned i = 0; i < _numHorizTiles; ++i)
    for (unsigned j = 0; j < _numVertTiles; ++j)
      _demandMap[i][j] = make_pair(0., 0.);
}

double ISPD04CongMap::getXDemand(const BBox &bbox) const {
  double total = 0.;
  unsigned minXCell, minYCell, maxXCell, maxYCell;
  getGridCoords(bbox, minXCell, minYCell, maxXCell, maxYCell);

  BBox adjustedBox = _layoutBBox.intersect(bbox);
  Point trans(-_layoutBBox.xMin, -_layoutBBox.yMin);
  adjustedBox.TranslateBy(trans);

  double invGridTileArea = 1. / (_tileWidth * _tileHeight);
  BBox gridCell;

  gridCell.xMin = static_cast<double>(minXCell) * _tileWidth;
  gridCell.xMax = static_cast<double>(minXCell + 1) * _tileWidth;
  for (unsigned i = minXCell; i <= maxXCell; ++i) {
    gridCell.yMin = static_cast<double>(minYCell) * _tileHeight;
    gridCell.yMax = static_cast<double>(minYCell + 1) * _tileHeight;
    for (unsigned j = minYCell; j <= maxYCell; ++j) {
      BBox intersection = adjustedBox.intersect(gridCell);

      total +=
          _demandMap[i][j].first * intersection.getArea() * invGridTileArea;

      gridCell.yMin += _tileHeight;
      gridCell.yMax += _tileHeight;
    }
    gridCell.xMin += _tileWidth;
    gridCell.xMax += _tileWidth;
  }

  total *= _tileWidth;
  return total;
}

double ISPD04CongMap::getMaxXDemand(const BBox &bbox) const {
  double maxDemand = 0.;
  unsigned minXCell, minYCell, maxXCell, maxYCell;
  getGridCoords(bbox, minXCell, minYCell, maxXCell, maxYCell);

  for (unsigned i = minXCell; i <= maxXCell; ++i) {
    for (unsigned j = minYCell; j <= maxYCell; ++j) {
      maxDemand = std::max(maxDemand, _demandMap[i][j].first);
    }
  }

  return (maxDemand * _tileWidth);
}

double ISPD04CongMap::getXSupply(const BBox &bbox) const {
  double total = 0.;
  unsigned minXCell, minYCell, maxXCell, maxYCell;
  getGridCoords(bbox, minXCell, minYCell, maxXCell, maxYCell);

  BBox adjustedBox = _layoutBBox.intersect(bbox);
  Point trans(-_layoutBBox.xMin, -_layoutBBox.yMin);
  adjustedBox.TranslateBy(trans);

  double invGridTileArea = 1. / (_tileWidth * _tileHeight);
  BBox gridCell;

  gridCell.xMin = static_cast<double>(minXCell) * _tileWidth;
  gridCell.xMax = static_cast<double>(minXCell + 1) * _tileWidth;
  for (unsigned i = minXCell; i <= maxXCell; ++i) {
    gridCell.yMin = static_cast<double>(minYCell) * _tileHeight;
    gridCell.yMax = static_cast<double>(minYCell + 1) * _tileHeight;
    for (unsigned j = minYCell; j <= maxYCell; ++j) {
      BBox intersection = adjustedBox.intersect(gridCell);

      total +=
          _supplyMap[i][j].first * intersection.getArea() * invGridTileArea;

      gridCell.yMin += _tileHeight;
      gridCell.yMax += _tileHeight;
    }
    gridCell.xMin += _tileWidth;
    gridCell.xMax += _tileWidth;
  }

  total *= _tileWidth;
  return total;
}

double ISPD04CongMap::getYDemand(const BBox &bbox) const {
  double total = 0.;
  unsigned minXCell, minYCell, maxXCell, maxYCell;
  getGridCoords(bbox, minXCell, minYCell, maxXCell, maxYCell);

  BBox adjustedBox = _layoutBBox.intersect(bbox);
  Point trans(-_layoutBBox.xMin, -_layoutBBox.yMin);
  adjustedBox.TranslateBy(trans);

  double invGridTileArea = 1. / (_tileWidth * _tileHeight);
  BBox gridCell;

  gridCell.xMin = static_cast<double>(minXCell) * _tileWidth;
  gridCell.xMax = static_cast<double>(minXCell + 1) * _tileWidth;
  for (unsigned i = minXCell; i <= maxXCell; ++i) {
    gridCell.yMin = static_cast<double>(minYCell) * _tileHeight;
    gridCell.yMax = static_cast<double>(minYCell + 1) * _tileHeight;
    for (unsigned j = minYCell; j <= maxYCell; ++j) {
      BBox intersection = adjustedBox.intersect(gridCell);

      total +=
          _demandMap[i][j].second * intersection.getArea() * invGridTileArea;

      gridCell.yMin += _tileHeight;
      gridCell.yMax += _tileHeight;
    }
    gridCell.xMin += _tileWidth;
    gridCell.xMax += _tileWidth;
  }

  total *= _tileHeight;
  return total;
}

double ISPD04CongMap::getYSupply(const BBox &bbox) const {
  double total = 0.;
  unsigned minXCell, minYCell, maxXCell, maxYCell;
  getGridCoords(bbox, minXCell, minYCell, maxXCell, maxYCell);

  BBox adjustedBox = _layoutBBox.intersect(bbox);
  Point trans(-_layoutBBox.xMin, -_layoutBBox.yMin);
  adjustedBox.TranslateBy(trans);

  double invGridTileArea = 1. / (_tileWidth * _tileHeight);
  BBox gridCell;

  gridCell.xMin = static_cast<double>(minXCell) * _tileWidth;
  gridCell.xMax = static_cast<double>(minXCell + 1) * _tileWidth;
  for (unsigned i = minXCell; i <= maxXCell; ++i) {
    gridCell.yMin = static_cast<double>(minYCell) * _tileHeight;
    gridCell.yMax = static_cast<double>(minYCell + 1) * _tileHeight;
    for (unsigned j = minYCell; j <= maxYCell; ++j) {
      BBox intersection = adjustedBox.intersect(gridCell);

      total +=
          _supplyMap[i][j].second * intersection.getArea() * invGridTileArea;

      gridCell.yMin += _tileHeight;
      gridCell.yMax += _tileHeight;
    }
    gridCell.xMin += _tileWidth;
    gridCell.xMax += _tileWidth;
  }

  total *= _tileHeight;
  return total;
}

pair<double, double> ISPD04CongMap::getDemand(const BBox &bbox) const {
  double xTotal = 0.;
  double yTotal = 0.;
  unsigned minXCell, minYCell, maxXCell, maxYCell;
  getGridCoords(bbox, minXCell, minYCell, maxXCell, maxYCell);

  BBox adjustedBox = _layoutBBox.intersect(bbox);
  Point trans(-_layoutBBox.xMin, -_layoutBBox.yMin);
  adjustedBox.TranslateBy(trans);

  double invGridTileArea = 1. / (_tileWidth * _tileHeight);
  BBox gridCell;

  gridCell.xMin = static_cast<double>(minXCell) * _tileWidth;
  gridCell.xMax = static_cast<double>(minXCell + 1) * _tileWidth;
  for (unsigned i = minXCell; i <= maxXCell; ++i) {
    gridCell.yMin = static_cast<double>(minYCell) * _tileHeight;
    gridCell.yMax = static_cast<double>(minYCell + 1) * _tileHeight;
    for (unsigned j = minYCell; j <= maxYCell; ++j) {
      BBox intersection = adjustedBox.intersect(gridCell);

      xTotal +=
          _demandMap[i][j].first * intersection.getArea() * invGridTileArea;
      yTotal +=
          _demandMap[i][j].second * intersection.getArea() * invGridTileArea;

      gridCell.yMin += _tileHeight;
      gridCell.yMax += _tileHeight;
    }
    gridCell.xMin += _tileWidth;
    gridCell.xMax += _tileWidth;
  }

  xTotal *= _tileWidth;
  yTotal *= _tileHeight;
  return make_pair(xTotal, yTotal);
}

pair<double, double> ISPD04CongMap::getDemand(unsigned xTile,
                                              unsigned yTile) const {
  pair<double, double> cong = _demandMap[xTile][yTile];

  cong.first *= _tileWidth;
  cong.second *= _tileHeight;

  return cong;
}

pair<double, double> ISPD04CongMap::getSupply(unsigned xTile,
                                              unsigned yTile) const {
  pair<double, double> supply = _supplyMap[xTile][yTile];

  supply.first *= _tileWidth;
  supply.second *= _tileHeight;

  return supply;
}

pair<double, double> ISPD04CongMap::getSupply(const BBox &bbox) const {
  double xTotal = 0.;
  double yTotal = 0.;
  unsigned minXCell, minYCell, maxXCell, maxYCell;
  getGridCoords(bbox, minXCell, minYCell, maxXCell, maxYCell);

  BBox adjustedBox = _layoutBBox.intersect(bbox);
  Point trans(-_layoutBBox.xMin, -_layoutBBox.yMin);
  adjustedBox.TranslateBy(trans);

  double invGridTileArea = 1. / (_tileWidth * _tileHeight);
  BBox gridCell;

  gridCell.xMin = static_cast<double>(minXCell) * _tileWidth;
  gridCell.xMax = static_cast<double>(minXCell + 1) * _tileWidth;
  for (unsigned i = minXCell; i <= maxXCell; ++i) {
    gridCell.yMin = static_cast<double>(minYCell) * _tileHeight;
    gridCell.yMax = static_cast<double>(minYCell + 1) * _tileHeight;
    for (unsigned j = minYCell; j <= maxYCell; ++j) {
      BBox intersection = adjustedBox.intersect(gridCell);

      xTotal +=
          _supplyMap[i][j].first * intersection.getArea() * invGridTileArea;
      yTotal +=
          _supplyMap[i][j].second * intersection.getArea() * invGridTileArea;

      gridCell.yMin += _tileHeight;
      gridCell.yMax += _tileHeight;
    }
    gridCell.xMin += _tileWidth;
    gridCell.xMax += _tileWidth;
  }

  xTotal *= _tileWidth;
  yTotal *= _tileHeight;
  return make_pair(xTotal, yTotal);
}

void ISPD04CongMap::calculateDemand(const HGraphWDimensions &hg) {
  resetDemand();
  double alpha = 1.;
  for (unsigned n = 0; n < hg.getNumEdges(); ++n) {
    addNet(hg, n, alpha);
  }
}

void ISPD04CongMap::buildSupplyMap(const RoutingResources &rr) {
  vector<pair<double, unsigned> > xPosOfLongitudinalTracks,
      yPosOfLatitudinalTracks;
  unsigned lowestRoutingLayer = rr.getLowestRoutingLayer();
  BBox coreRegion = rr.getCoreRegion();

  for (itTrackPattern t = rr.trackPatternsBegin(); t != rr.trackPatternsEnd();
       ++t) {
    if (!rr.isRoutingLayer(t->getLayerNum())) continue;
    if (t->getLayerNum() == lowestRoutingLayer) continue;
    if (t->isLongitudinal()) {
      if (rr.getLayerPreference(t->getLayerNum()) == RoutingResources::HORIZ)
        continue;
      double xcoord = t->getStart();
      for (unsigned i = 0; i < t->getNumTracks(); ++i) {
        if (coreRegion.xMin <= xcoord && xcoord <= coreRegion.xMax &&
            _layoutBBox.xMin <= xcoord && xcoord <= _layoutBBox.xMax) {
          xPosOfLongitudinalTracks.push_back(
              make_pair(xcoord, t->getLayerNum()));
        }
        xcoord += t->getSpacing();
      }
    } else {
      if (rr.getLayerPreference(t->getLayerNum()) == RoutingResources::VERT)
        continue;
      double ycoord = t->getStart();
      for (unsigned i = 0; i < t->getNumTracks(); ++i) {
        if (coreRegion.yMin <= ycoord && ycoord <= coreRegion.yMax &&
            _layoutBBox.yMin <= ycoord && ycoord <= _layoutBBox.yMax) {
          yPosOfLatitudinalTracks.push_back(
              make_pair(ycoord, t->getLayerNum()));
        }
        ycoord += t->getSpacing();
      }
    }
  }

  sort(xPosOfLongitudinalTracks.begin(), xPosOfLongitudinalTracks.end());
  vector<pair<double, unsigned> >::iterator new_end =
      unique(xPosOfLongitudinalTracks.begin(), xPosOfLongitudinalTracks.end());
  xPosOfLongitudinalTracks.erase(new_end, xPosOfLongitudinalTracks.end());
  sort(yPosOfLatitudinalTracks.begin(), yPosOfLatitudinalTracks.end());
  new_end =
      unique(yPosOfLatitudinalTracks.begin(), yPosOfLatitudinalTracks.end());
  yPosOfLatitudinalTracks.erase(new_end, yPosOfLatitudinalTracks.end());

  unsigned minXCell, minYCell, maxXCell, maxYCell;
  const double alpha = 1.;
  for (unsigned i = 0; i < xPosOfLongitudinalTracks.size(); ++i) {
    BBox track = coreRegion;
    track.xMin = track.xMax = xPosOfLongitudinalTracks[i].first;
    getGridCoords(track, minXCell, minYCell, maxXCell, maxYCell);
    addFlatNet(_supplyMap, track, alpha, minXCell, minYCell, maxXCell,
               maxYCell);
  }
  for (unsigned i = 0; i < yPosOfLatitudinalTracks.size(); ++i) {
    BBox track = coreRegion;
    track.yMin = track.yMax = yPosOfLatitudinalTracks[i].first;
    getGridCoords(track, minXCell, minYCell, maxXCell, maxYCell);
    addFlatNet(_supplyMap, track, alpha, minXCell, minYCell, maxXCell,
               maxYCell);
  }
}

void ISPD04CongMap::getGridCoords(const BBox &box, unsigned &minXCell,
                                  unsigned &minYCell, unsigned &maxXCell,
                                  unsigned &maxYCell) const {
  BBox adjustedBox = _layoutBBox.intersect(box);
  Point trans(-_layoutBBox.xMin, -_layoutBBox.yMin);
  adjustedBox.TranslateBy(trans);

  if (adjustedBox.xMin <= 0.) {
    minXCell = 0;
  } else {
    minXCell = static_cast<unsigned>(floor(adjustedBox.xMin / _tileWidth));
  }
  if (minXCell == _numHorizTiles) --minXCell;

  if (adjustedBox.xMax <= 0.) {
    maxXCell = 0;
  } else {
    maxXCell = static_cast<unsigned>(floor(adjustedBox.xMax / _tileWidth));
  }
  if (maxXCell == _numHorizTiles) --maxXCell;

  if (adjustedBox.yMin <= 0.) {
    minYCell = 0;
  } else {
    minYCell = static_cast<unsigned>(floor(adjustedBox.yMin / _tileHeight));
  }
  if (minYCell == _numVertTiles) --minYCell;

  if (adjustedBox.yMax <= 0.) {
    maxYCell = 0;
  } else {
    maxYCell = static_cast<unsigned>(floor(adjustedBox.yMax / _tileHeight));
  }
  if (maxYCell == _numVertTiles) --maxYCell;
}

vector<Point> *ISPD04CongMap::getPointsOnNet(const HGraphWDimensions &hg,
                                             unsigned netId) const {
  const HGFEdge &edge = hg.getEdgeByIdx(netId);

  vector<Point> *pointsOnNet = new vector<Point>();

  unsigned nodeOffset = 0;
  for (itHGFNodeLocal nodeIt = edge.nodesBegin(); nodeIt != edge.nodesEnd();
       ++nodeIt, ++nodeOffset) {
    const HGFNode &node = (**nodeIt);
    Point pinOffset = hg.nodesOnEdgePinOffset(nodeOffset, netId,
                                              _pl.getOrient(node.getIndex()));

    unsigned cellId = node.getIndex();

    pointsOnNet->push_back(_pl[cellId] + pinOffset);
  }

  sort(pointsOnNet->begin(), pointsOnNet->end());
  vector<Point>::iterator new_end =
      unique(pointsOnNet->begin(), pointsOnNet->end());
  pointsOnNet->erase(new_end, pointsOnNet->end());

  return pointsOnNet;
}

pair<double, double> calculateIllegality(
    const pair<Point, Point> &pinPair, const std::vector<BBox> &obs) {
  double xFirstIntersectionLen = 0.;
  double yFirstIntersectionLen = 0.;

  double xMin = min(pinPair.first.x, pinPair.second.x);
  double xMax = max(pinPair.first.x, pinPair.second.x);
  double yMin = min(pinPair.first.y, pinPair.second.y);
  double yMax = max(pinPair.first.y, pinPair.second.y);

  for (unsigned j = 0; j < obs.size(); ++j) {
    // test the x first path
    if (lessThanDouble(pinPair.first.y, obs[j].yMax) &&
        greaterThanDouble(pinPair.first.y, obs[j].yMin) &&
        lessThanDouble(xMin, obs[j].xMax) &&
        greaterThanDouble(xMax, obs[j].xMin)) {
      // x segment of x first path fails
      xFirstIntersectionLen += min(xMax, obs[j].xMax) - max(xMin, obs[j].xMin);
    }
    if (lessThanDouble(pinPair.second.x, obs[j].xMax) &&
        greaterThanDouble(pinPair.second.x, obs[j].xMin) &&
        lessThanDouble(yMin, obs[j].yMax) &&
        greaterThanDouble(yMax, obs[j].yMin)) {
      // y segment of x first path fails
      xFirstIntersectionLen += min(yMax, obs[j].yMax) - max(yMin, obs[j].yMin);
    }

    // test the y first path
    if (lessThanDouble(pinPair.second.y, obs[j].yMax) &&
        greaterThanDouble(pinPair.second.y, obs[j].yMin) &&
        lessThanDouble(xMin, obs[j].xMax) &&
        greaterThanDouble(xMax, obs[j].xMin)) {
      // x segment of y first path fails
      yFirstIntersectionLen += min(xMax, obs[j].xMax) - max(xMin, obs[j].xMin);
    } else if (lessThanDouble(pinPair.first.x, obs[j].xMax) &&
               greaterThanDouble(pinPair.first.x, obs[j].xMin) &&
               lessThanDouble(yMin, obs[j].yMax) &&
               greaterThanDouble(yMax, obs[j].yMin)) {
      // y segment of y first path fails
      yFirstIntersectionLen += min(yMax, obs[j].yMax) - max(yMin, obs[j].yMin);
    }
  }

  return make_pair(xFirstIntersectionLen, yFirstIntersectionLen);
}

bool insideObs(const Point &pin, const std::vector<BBox> &obs) {
  for (unsigned i = 0; i < obs.size(); ++i) {
    if (obs[i].hasInside(pin)) return true;
  }
  return false;
}

vector<pair<Point, Point> > *ISPD04CongMap::getPinPairs(
    const HGraphWDimensions &hg, unsigned netId,
    vector<Point> *suppliedPoints) const {
  vector<Point> pointsOnNet;

  if (suppliedPoints == NULL) {
    vector<Point> *temp = getPointsOnNet(hg, netId);
    pointsOnNet.swap(*temp);
    delete temp;
  } else {
    pointsOnNet.swap(*suppliedPoints);
  }

  if (pointsOnNet.size() <= 1) return new vector<pair<Point, Point> >(0);

  vector<pair<Point, Point> > *pinPairs = new vector<pair<Point, Point> >();
  pinPairs->reserve(pointsOnNet.size() - 1);

#ifdef USEFLUTE
  if (pointsOnNet.size() > MAXD || !_useSteiner) {
#endif
    if (pointsOnNet.size() == 2) {
      pinPairs->push_back(make_pair(pointsOnNet[0], pointsOnNet[1]));
    } else {
      // call MST to get pin pairs
      vector<fastSteiner::Point> points(pointsOnNet.size());
      vector<long> parent(pointsOnNet.size());

      for (unsigned i = 0; i < pointsOnNet.size(); ++i) {
        points[i].x = pointsOnNet[i].x;
        points[i].y = pointsOnNet[i].y;
      }

      fastSteiner::mst2(pointsOnNet.size(), &points[0], &parent[0]);

      for (unsigned long i = 0; i < pointsOnNet.size(); ++i) {
        if (static_cast<unsigned long>(parent[i]) == i) continue;

        pinPairs->push_back(make_pair(pointsOnNet[i], pointsOnNet[parent[i]]));
      }

      abkfatal(pinPairs->size() == pointsOnNet.size() - 1,
               "Bad number of pin pairs from MST");

      fastSteiner::mst2_package_done();
    }
#ifdef USEFLUTE
  } else {
    double *xcoords = new double[pointsOnNet.size()];
    double *ycoords = new double[pointsOnNet.size()];

    for (unsigned i = 0; i < pointsOnNet.size(); ++i) {
      xcoords[i] = pointsOnNet[i].x;
      ycoords[i] = pointsOnNet[i].y;
    }

    Tree fltree;
    unsigned legal = 1;
    bool impossible = false;
    if (_obstacles == NULL) {
      fltree = flute(pointsOnNet.size(), xcoords, ycoords, ACCURACY);
    } else {
      vector<unsigned> relevantobs(_obstacles->size(), 0);
      for (unsigned i = 0; i < _obstacles->size(); ++i) {
        relevantobs[i] = i;
      }
      fltree = flautist(pointsOnNet.size(), xcoords, ycoords, ACCURACY,
                        *_obstacles, relevantobs, legal);

      for (unsigned i = 0; i < pointsOnNet.size(); ++i) {
        if (insideObs(pointsOnNet[i], *_obstacles)) {
          impossible = true;
          break;
        }
      }
    }

    /*
         if(impossible)
         {
           cout << "got an impossible net" << endl;
           cout << "net ID" << netId << endl;
           cout << "obstacles: " << endl;
           for(unsigned i = 0; i < _obstacles->size(); ++i)
           {
             cout << "\t" << (*_obstacles)[i] << endl;
           }
           cout << "pin points: " << endl;
           const HGFEdge &edge = hg.getEdgeByIdx(netId);

           unsigned nodeOffset = 0;
           for(itHGFNodeLocal nodeIt = edge.nodesBegin();
                              nodeIt != edge.nodesEnd(); ++nodeIt, ++nodeOffset)
           {
             const HGFNode& node=(**nodeIt);
             Point pinOffset = hg.nodesOnEdgePinOffset(nodeOffset, netId,
                               _pl.getOrient(node.getIndex()));

             unsigned cellId = node.getIndex();
             cout << "\t" << _pl[cellId] << " + " << pinOffset << " = " <<
       _pl[cellId] + pinOffset << endl;
           }
           cout << endl;
         }
    */

    if (legal || impossible || pointsOnNet.size() > 20) {
      for (int i = 0; i < 2 * fltree.deg - 2; ++i) {
        int n = fltree.branch[i].n;
        if (i == n) continue;
        Point first(fltree.branch[i].x, fltree.branch[i].y);
        Point second(fltree.branch[n].x, fltree.branch[n].y);

        if (first != second) {
          pinPairs->push_back(make_pair(first, second));
        }
      }

      free(fltree.branch);
      delete[] xcoords;
      delete[] ycoords;
    } else {
      free(fltree.branch);
      delete[] xcoords;
      delete[] ycoords;

      mst_simple steinerbuilder;
      vector<int> relevantobs(_obstacles->size(), 0);
      for (int i = 0; i < (int)_obstacles->size(); ++i) {
        relevantobs[i] = i;
      }
      steinerbuilder.getStePoints(*_obstacles, relevantobs, pointsOnNet);
      vector<int> Tree;
      double length = 0.;
      vector<double> allLength;
      steinerbuilder.calMST(*_obstacles, pointsOnNet, Tree, length, allLength);
      for (unsigned i = 0; i < Tree.size(); ++i) {
        if (Tree[i] == -1) continue;
        Point first = pointsOnNet[i];
        Point second = pointsOnNet[Tree[i]];
        pinPairs->push_back(make_pair(first, second));
      }
    }
  }
#endif

  return pinPairs;
}

void ISPD04CongMap::addNet(const HGraphWDimensions &hg, unsigned netId,
                           double alpha) {
  vector<pair<Point, Point> > pinPairs, *temp = getPinPairs(hg, netId);

  pinPairs.swap(*temp);
  delete temp;

  if (pinPairs.size() == 0) return;

  unsigned minXCell, minYCell, maxXCell, maxYCell;
  for (unsigned i = 0; i < pinPairs.size(); ++i) {
    BBox pinPairBBox;
    pinPairBBox += pinPairs[i].first;
    pinPairBBox += pinPairs[i].second;

    if (lessThanDouble(_layoutBBox.xMax, pinPairBBox.xMin) ||
        lessThanDouble(_layoutBBox.yMax, pinPairBBox.yMin) ||
        greaterThanDouble(_layoutBBox.xMin, pinPairBBox.xMax) ||
        greaterThanDouble(_layoutBBox.yMin, pinPairBBox.yMax)) {
      continue;
    }

    getGridCoords(pinPairBBox, minXCell, minYCell, maxXCell, maxYCell);

    bool flipped = (pinPairs[i].first != pinPairBBox.getBottomLeft() &&
                    pinPairs[i].second != pinPairBBox.getBottomLeft());

    pair<double, double> illegality = make_pair(0., 0.);
    if (_obstacles != NULL) {
      illegality = calculateIllegality(pinPairs[i], *_obstacles);
    }

    if (minXCell == maxXCell && minYCell == maxYCell) {
      // a short net
      addShortNet(_demandMap, pinPairBBox, alpha, minXCell, minYCell);
    } else if (minXCell == maxXCell || minYCell == maxYCell) {
      // a flat net
      addFlatNet(_demandMap, pinPairBBox, alpha, minXCell, minYCell, maxXCell,
                 maxYCell);
    } else {
      if (_obstacles == NULL) {
        if (maxXCell - minXCell == 1 && maxYCell - minYCell == 1) {
          // a net with only L shapes
          addLUsage(_demandMap, pinPairBBox, alpha, minXCell, minYCell,
                    maxXCell, maxYCell);
        } else {
          // a regular net
          addLUsage(_demandMap, pinPairBBox, alpha * _alpha, minXCell, minYCell,
                    maxXCell, maxYCell);
          addZUsage(_demandMap, pinPairBBox, alpha * (1. - _alpha), flipped,
                    minXCell, minYCell, maxXCell, maxYCell);
        }
      } else {
        if (lessThanDouble(illegality.first, illegality.second)) {
          // must be two flat nets, x direction first
          if (flipped) {
            // x direction first, flipped
            BBox temp = pinPairBBox;
            temp.yMin = temp.yMax;
            addFlatNet(_demandMap, temp, alpha, minXCell, maxYCell, maxXCell,
                       maxYCell);

            temp = pinPairBBox;
            temp.xMin = temp.xMax;
            addFlatNet(_demandMap, temp, alpha, maxXCell, minYCell, maxXCell,
                       maxYCell);
          } else {
            // x direction first, unflipped
            BBox temp = pinPairBBox;
            temp.yMax = temp.yMin;
            addFlatNet(_demandMap, temp, alpha, minXCell, minYCell, maxXCell,
                       minYCell);

            temp = pinPairBBox;
            temp.xMin = temp.xMax;
            addFlatNet(_demandMap, temp, alpha, maxXCell, minYCell, maxXCell,
                       maxYCell);
          }
        } else if (lessThanDouble(illegality.second, illegality.first)) {
          // must be two flat nets, y direction first
          if (flipped) {
            // y direction first, flipped
            BBox temp = pinPairBBox;
            temp.xMax = temp.xMin;
            addFlatNet(_demandMap, temp, alpha, minXCell, minYCell, minXCell,
                       maxYCell);

            temp = pinPairBBox;
            temp.yMax = temp.yMin;
            addFlatNet(_demandMap, temp, alpha, minXCell, minYCell, maxXCell,
                       minYCell);
          } else {
            // y direction first, unflipped
            BBox temp = pinPairBBox;
            temp.xMax = temp.xMin;
            addFlatNet(_demandMap, temp, alpha, minXCell, minYCell, minXCell,
                       maxYCell);

            temp = pinPairBBox;
            temp.yMin = temp.yMax;
            addFlatNet(_demandMap, temp, alpha, minXCell, maxYCell, maxXCell,
                       maxYCell);
          }
        } else  // illegality.first == illegality.second
        {
          if (maxXCell - minXCell == 1 && maxYCell - minYCell == 1) {
            // a net with only L shapes
            addLUsage(_demandMap, pinPairBBox, alpha, minXCell, minYCell,
                      maxXCell, maxYCell);
          } else {
            // a regular net
            addLUsage(_demandMap, pinPairBBox, alpha * _alpha, minXCell,
                      minYCell, maxXCell, maxYCell);
            addZUsage(_demandMap, pinPairBBox, alpha * (1. - _alpha), flipped,
                      minXCell, minYCell, maxXCell, maxYCell);
          }
        }
      }
    }
  }
}

void ISPD04CongMap::addNetToTemp(const HGraphWDimensions &hg, unsigned netId,
                                 double alpha, unsigned movingCell,
                                 const Point &newPos) const {
  // get the new points on the net
  const HGFEdge &edge = hg.getEdgeByIdx(netId);

  vector<Point> pointsOnNet;

  unsigned nodeOffset = 0;
  for (itHGFNodeLocal nodeIt = edge.nodesBegin(); nodeIt != edge.nodesEnd();
       ++nodeIt, ++nodeOffset) {
    const HGFNode &node = (**nodeIt);
    Point pinOffset = hg.nodesOnEdgePinOffset(nodeOffset, netId,
                                              _pl.getOrient(node.getIndex()));

    unsigned cellId = node.getIndex();

    if (cellId == movingCell) {
      pointsOnNet.push_back(newPos + pinOffset);
    } else {
      pointsOnNet.push_back(_pl[cellId] + pinOffset);
    }
  }

  sort(pointsOnNet.begin(), pointsOnNet.end());
  vector<Point>::iterator new_end =
      unique(pointsOnNet.begin(), pointsOnNet.end());
  pointsOnNet.erase(new_end, pointsOnNet.end());

  vector<pair<Point, Point> > pinPairs,
      *temp = getPinPairs(hg, netId, &pointsOnNet);

  pinPairs.swap(*temp);
  delete temp;

  if (pinPairs.size() == 0) return;

  unsigned minXCell, minYCell, maxXCell, maxYCell;
  for (unsigned i = 0; i < pinPairs.size(); ++i) {
    BBox pinPairBBox;
    pinPairBBox += pinPairs[i].first;
    pinPairBBox += pinPairs[i].second;

    if (lessThanDouble(_layoutBBox.xMax, pinPairBBox.xMin) ||
        lessThanDouble(_layoutBBox.yMax, pinPairBBox.yMin) ||
        greaterThanDouble(_layoutBBox.xMin, pinPairBBox.xMax) ||
        greaterThanDouble(_layoutBBox.yMin, pinPairBBox.yMax)) {
      continue;
    }

    getGridCoords(pinPairBBox, minXCell, minYCell, maxXCell, maxYCell);

    bool flipped = (pinPairs[i].first != pinPairBBox.getBottomLeft() &&
                    pinPairs[i].second != pinPairBBox.getBottomLeft());

    if (minXCell == maxXCell && minYCell == maxYCell) {
      // a short net
      addShortNet(_tempMap, pinPairBBox, alpha, minXCell, minYCell);
    } else if (minXCell == maxXCell || minYCell == maxYCell) {
      // a flat net
      addFlatNet(_tempMap, pinPairBBox, alpha, minXCell, minYCell, maxXCell,
                 maxYCell);
    } else if (maxXCell - minXCell == 1 && maxYCell - minYCell == 1) {
      // a net with only L shapes
      addLUsage(_tempMap, pinPairBBox, alpha, minXCell, minYCell, maxXCell,
                maxYCell);
    } else {
      // a regular net
      addLUsage(_tempMap, pinPairBBox, alpha * _alpha, minXCell, minYCell,
                maxXCell, maxYCell);
      addZUsage(_tempMap, pinPairBBox, alpha * (1. - _alpha), flipped, minXCell,
                minYCell, maxXCell, maxYCell);
    }
  }
}

void ISPD04CongMap::addShortNet(
    vector<vector<std::pair<double, double> > > &theMap, const BBox &twoPinNet,
    double alpha, unsigned xCell, unsigned yCell) const {
  theMap[xCell][yCell].first += alpha * twoPinNet.getWidth() / _tileWidth;
  theMap[xCell][yCell].second += alpha * twoPinNet.getHeight() / _tileHeight;
}

void ISPD04CongMap::addFlatNet(vector<vector<pair<double, double> > > &theMap,
                               const BBox &twoPinNet, double alpha,
                               unsigned minXCell, unsigned minYCell,
                               unsigned maxXCell, unsigned maxYCell) const {
  BBox adjustedNet = _layoutBBox.intersect(twoPinNet);
  Point trans(-_layoutBBox.xMin, -_layoutBBox.yMin);
  adjustedNet.TranslateBy(trans);

  if (minXCell == maxXCell) {
    double xUsage = alpha * adjustedNet.getWidth() /
                    (_tileWidth * static_cast<double>(maxYCell - minYCell + 1));
    theMap[minXCell][minYCell].first += xUsage;
    theMap[minXCell][minYCell].second +=
        alpha *
        max(0., static_cast<double>(minYCell + 1) * _tileHeight -
                    adjustedNet.yMin) /
        _tileHeight;
    for (unsigned i = minYCell + 1; i < maxYCell; ++i) {
      theMap[minXCell][i].first += xUsage;
      theMap[minXCell][i].second += alpha;
    }
    theMap[minXCell][maxYCell].first += xUsage;
    theMap[minXCell][maxYCell].second +=
        alpha *
        max(0.,
            adjustedNet.yMax - static_cast<double>(maxYCell) * _tileHeight) /
        _tileHeight;
  } else {
    // minYCell == maxYCell
    double yUsage =
        alpha * adjustedNet.getHeight() /
        (_tileHeight * static_cast<double>(maxXCell - minXCell + 1));
    theMap[minXCell][minYCell].first +=
        alpha *
        max(0.,
            static_cast<double>(minXCell + 1) * _tileWidth - adjustedNet.xMin) /
        _tileWidth;
    theMap[minXCell][minYCell].second += yUsage;
    for (unsigned i = minXCell + 1; i < maxXCell; ++i) {
      theMap[i][minYCell].first += alpha;
      theMap[i][minYCell].second += yUsage;
    }
    theMap[maxXCell][minYCell].first +=
        alpha *
        max(0., adjustedNet.xMax - static_cast<double>(maxXCell) * _tileWidth) /
        _tileHeight;
    theMap[maxXCell][minYCell].second += yUsage;
  }
}

void ISPD04CongMap::addLUsage(
    vector<vector<std::pair<double, double> > > &theMap, const BBox &twoPinNet,
    double alpha, unsigned minXCell, unsigned minYCell, unsigned maxXCell,
    unsigned maxYCell) const {
  BBox adjustedNet = _layoutBBox.intersect(twoPinNet);
  double scale_factor = 0.5 * alpha;

  // L usage is just 4 flat nets with appropriate weighting

  // lower horiz net
  BBox temp = adjustedNet;
  temp.yMax = temp.yMin;
  addFlatNet(theMap, temp, scale_factor, minXCell, minYCell, maxXCell,
             minYCell);

  // upper horiz net
  temp = adjustedNet;
  temp.yMin = temp.yMax;
  addFlatNet(theMap, temp, scale_factor, minXCell, maxYCell, maxXCell,
             maxYCell);

  // left vert net
  temp = adjustedNet;
  temp.xMax = temp.xMin;
  addFlatNet(theMap, temp, scale_factor, minXCell, minYCell, minXCell,
             maxYCell);

  // right vert net
  temp = adjustedNet;
  temp.xMin = temp.xMax;
  addFlatNet(theMap, temp, scale_factor, maxXCell, minYCell, maxXCell,
             maxYCell);
}

void ISPD04CongMap::addZUsage(
    vector<vector<std::pair<double, double> > > &theMap, const BBox &twoPinNet,
    double alpha, bool flipped, unsigned minXCell, unsigned minYCell,
    unsigned maxXCell, unsigned maxYCell) const {
  BBox adjustedNet = _layoutBBox.intersect(twoPinNet);
  Point trans(-_layoutBBox.xMin, -_layoutBBox.yMin);
  adjustedNet.TranslateBy(trans);

  unsigned gridCellsWideM2 = maxXCell - minXCell - 1;
  unsigned gridCellsHighM2 = maxYCell - minYCell - 1;
  double gridCellsWideM2D = static_cast<double>(gridCellsWideM2);
  double gridCellsHighM2D = static_cast<double>(gridCellsHighM2);

  // process vertical z shapes first
  if (gridCellsWideM2 > 0) {
    double scale_factor =
        (alpha * gridCellsWideM2D) / (gridCellsWideM2D + gridCellsHighM2D);
    double vertUsageBot =
        scale_factor *
        max(0., static_cast<double>(minYCell + 1) * _tileHeight -
                    adjustedNet.yMin) /
        (_tileHeight * gridCellsWideM2D);
    double vertUsageMid = scale_factor / gridCellsWideM2D;
    double vertUsageTop =
        scale_factor *
        max(0.,
            adjustedNet.yMax - static_cast<double>(maxYCell) * _tileHeight) /
        (_tileHeight * gridCellsWideM2D);
    if (flipped) {
      theMap[minXCell][maxYCell].first +=
          scale_factor *
          max(0., static_cast<double>(minXCell + 1) * _tileWidth -
                      adjustedNet.xMin) /
          _tileWidth;
      for (unsigned i = 1; i <= gridCellsWideM2; ++i) {
        double iD = static_cast<double>(i);
        theMap[minXCell + i][maxYCell].first +=
            scale_factor * (gridCellsWideM2D - iD + 0.5) / gridCellsWideM2D;
        theMap[minXCell + i][maxYCell].second += vertUsageBot;
        for (unsigned j = 1; j <= gridCellsHighM2; ++j) {
          theMap[minXCell + i][maxYCell - j].second += vertUsageMid;
        }
        theMap[minXCell + i][minYCell].first +=
            scale_factor * (iD - 0.5) / gridCellsWideM2D;
        theMap[minXCell + i][minYCell].second += vertUsageTop;
      }
      theMap[maxXCell][minYCell].first +=
          scale_factor *
          max(0.,
              adjustedNet.xMax - static_cast<double>(maxXCell) * _tileWidth) /
          _tileWidth;
    } else {
      theMap[minXCell][minYCell].first +=
          scale_factor *
          max(0., static_cast<double>(minXCell + 1) * _tileWidth -
                      adjustedNet.xMin) /
          _tileWidth;
      for (unsigned i = 1; i <= gridCellsWideM2; ++i) {
        double iD = static_cast<double>(i);
        theMap[minXCell + i][minYCell].first +=
            scale_factor * (gridCellsWideM2D - iD + 0.5) / gridCellsWideM2D;
        theMap[minXCell + i][minYCell].second += vertUsageBot;
        for (unsigned j = 1; j <= gridCellsHighM2; ++j) {
          theMap[minXCell + i][minYCell + j].second += vertUsageMid;
        }
        theMap[minXCell + i][maxYCell].first +=
            scale_factor * (iD - 0.5) / gridCellsWideM2D;
        theMap[minXCell + i][maxYCell].second += vertUsageTop;
      }
      theMap[maxXCell][maxYCell].first +=
          scale_factor *
          max(0.,
              adjustedNet.xMax - static_cast<double>(maxXCell) * _tileWidth) /
          _tileWidth;
    }
  }

  // now process the horizontal z shapes
  if (gridCellsHighM2 > 0) {
    double scale_factor =
        (alpha * gridCellsHighM2D) / (gridCellsWideM2D + gridCellsHighM2D);

    double horizUsageLeft =
        scale_factor *
        max(0.,
            static_cast<double>(minXCell + 1) * _tileWidth - adjustedNet.xMin) /
        (_tileWidth * gridCellsHighM2D);
    double horizUsageMid = scale_factor / gridCellsHighM2D;
    double horizUsageRight =
        scale_factor *
        max(0., adjustedNet.xMax - static_cast<double>(maxXCell) * _tileWidth) /
        (_tileWidth * gridCellsHighM2D);
    if (flipped) {
      theMap[maxXCell][minYCell].second +=
          scale_factor *
          max(0., static_cast<double>(minYCell + 1) * _tileHeight -
                      adjustedNet.yMin) /
          _tileHeight;
      for (unsigned i = 1; i <= gridCellsHighM2; ++i) {
        double iD = static_cast<double>(i);
        theMap[maxXCell][minYCell + i].first += horizUsageLeft;
        theMap[maxXCell][minYCell + i].second +=
            scale_factor * (gridCellsHighM2D - iD + 0.5) / gridCellsHighM2D;
        for (unsigned j = 1; j <= gridCellsWideM2; ++j) {
          theMap[maxXCell - j][minYCell + i].first += horizUsageMid;
        }
        theMap[minXCell][minYCell + i].first += horizUsageRight;
        theMap[minXCell][minYCell + i].second +=
            scale_factor * (iD - 0.5) / gridCellsHighM2D;
      }
      theMap[minXCell][maxYCell].second +=
          scale_factor *
          max(0.,
              adjustedNet.yMax - static_cast<double>(maxYCell) * _tileHeight) /
          _tileHeight;
    } else {
      theMap[minXCell][minYCell].second +=
          scale_factor *
          max(0., static_cast<double>(minYCell + 1) * _tileHeight -
                      adjustedNet.yMin) /
          _tileHeight;
      for (unsigned i = 1; i <= gridCellsHighM2; ++i) {
        double iD = static_cast<double>(i);
        theMap[minXCell][minYCell + i].first += horizUsageLeft;
        theMap[minXCell][minYCell + i].second +=
            scale_factor * (gridCellsHighM2D - iD + 0.5) / gridCellsHighM2D;
        for (unsigned j = 1; j <= gridCellsWideM2; ++j) {
          theMap[minXCell + j][minYCell + i].first += horizUsageMid;
        }
        theMap[maxXCell][minYCell + i].first += horizUsageRight;
        theMap[maxXCell][minYCell + i].second +=
            scale_factor * (iD - 0.5) / gridCellsHighM2D;
      }
      theMap[maxXCell][maxYCell].second +=
          scale_factor *
          max(0.,
              adjustedNet.yMax - static_cast<double>(maxYCell) * _tileHeight) /
          _tileHeight;
    }
  }
}

// plot congestion map in gnuplot format
void ISPD04CongMap::plot(const std::string& baseFileName) const {
  std::string fnameX = baseFileName.c_str();
  fnameX += "X.plt";
  std::string fnameY = baseFileName.c_str();
  fnameY += "Y.plt";
  std::string fname = baseFileName.c_str();
  fname += ".plt";

  cout << "Saving " << fnameX << ", " << fnameY << " and " << fname << endl;

  ofstream gpFileX(fnameX.c_str());
  ofstream gpFileY(fnameY.c_str());
  ofstream gpFile(fname.c_str());

  gpFileX << "#Use this file as a script for gnuplot" << endl
          << "#(See http://www.gnuplot.info/ for details)" << endl
          << "#" << TimeStamp() << User() << Platform() << endl
          << endl
          << "set title ' " << baseFileName << "X , #Cells= " << _pl.size()
          << ", #Nets= " << _netbitboard.getSize() << endl
          << endl
          << "set nokey" << endl
          << "#   Uncomment these two lines starting with \"set\"" << endl
          << "#   to save an EPS file for inclusion into a latex document"
          << endl
          << "# set terminal postscript eps color solid 10" << endl
          << "# set output \"" << baseFileName << "X.eps\"" << endl
          << endl
          << "#   Uncomment these two lines starting with \"set\"" << endl
          << "#   to save a PS file for printing" << endl
          << "# set terminal postscript portrait color solid 8" << endl
          << "# set output \"" << baseFileName << "X.ps\"" << endl
          << endl
          << endl
          << "set view 55,45 " << endl
          << endl
          << endl
          << "splot '-' w l " << endl
          << endl;

  gpFileY << "#Use this file as a script for gnuplot" << endl
          << "#(See http://www.gnuplot.info/ for details)" << endl
          << "#" << TimeStamp() << User() << Platform() << endl
          << endl
          << "set title ' " << baseFileName << "Y , #Cells= " << _pl.size()
          << ", #Nets= " << _netbitboard.getSize() << endl
          << endl
          << "set nokey" << endl
          << "#   Uncomment these two lines starting with \"set\"" << endl
          << "#   to save an EPS file for inclusion into a latex document"
          << endl
          << "# set terminal postscript eps color solid 10" << endl
          << "# set output \"" << baseFileName << "Y.eps\"" << endl
          << endl
          << "#   Uncomment these two lines starting with \"set\"" << endl
          << "#   to save a PS file for printing" << endl
          << "# set terminal postscript portrait color solid 8" << endl
          << "# set output \"" << baseFileName << "Y.ps\"" << endl
          << endl
          << endl
          << "set view 55,45 " << endl
          << endl
          << endl
          << "splot '-' w l " << endl
          << endl;

  gpFile << "#Use this file as a script for gnuplot" << endl
         << "#(See http://www.gnuplot.info/ for details)" << endl
         << "#" << TimeStamp() << User() << Platform() << endl
         << endl
         << "set title ' " << baseFileName << " , #Cells= " << _pl.size()
         << ", #Nets= " << _netbitboard.getSize() << endl
         << endl
         << "set nokey" << endl
         << "#   Uncomment these two lines starting with \"set\"" << endl
         << "#   to save an EPS file for inclusion into a latex document"
         << endl
         << "# set terminal postscript eps color solid 10" << endl
         << "# set output \"" << baseFileName << ".eps\"" << endl
         << endl
         << "#   Uncomment these two lines starting with \"set\"" << endl
         << "#   to save a PS file for printing" << endl
         << "# set terminal postscript portrait color solid 8" << endl
         << "# set output \"" << baseFileName << ".ps\"" << endl
         << endl
         << endl
         << "set view 55,45 " << endl
         << endl
         << endl
         << "splot '-' w l " << endl
         << endl;

  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    for (unsigned j = 0; j < _numVertTiles; ++j) {
      double regXMin = _layoutBBox.xMin + static_cast<double>(i) * _tileWidth;
      double regXMax =
          _layoutBBox.xMin + static_cast<double>(i + 1) * _tileWidth;
      double regYMin = _layoutBBox.yMin + static_cast<double>(j) * _tileHeight;
      double regYMax =
          _layoutBBox.yMin + static_cast<double>(j + 1) * _tileHeight;

      gpFileX << regXMin << " " << regYMin << " "
              << _demandMap[i][j].first * _tileWidth << endl;
      gpFileX << regXMax << " " << regYMin << " "
              << _demandMap[i][j].first * _tileWidth << endl;
      gpFileX << regXMax << " " << regYMax << " "
              << _demandMap[i][j].first * _tileWidth << endl;
      gpFileX << regXMin << " " << regYMax << " "
              << _demandMap[i][j].first * _tileWidth << endl;
      gpFileX << regXMin << " " << regYMin << " "
              << _demandMap[i][j].first * _tileWidth << endl
              << endl
              << endl;

      gpFileY << regXMin << " " << regYMin << " "
              << _demandMap[i][j].second * _tileHeight << endl;
      gpFileY << regXMax << " " << regYMin << " "
              << _demandMap[i][j].second * _tileHeight << endl;
      gpFileY << regXMax << " " << regYMax << " "
              << _demandMap[i][j].second * _tileHeight << endl;
      gpFileY << regXMin << " " << regYMax << " "
              << _demandMap[i][j].second * _tileHeight << endl;
      gpFileY << regXMin << " " << regYMin << " "
              << _demandMap[i][j].second * _tileHeight << endl
              << endl
              << endl;

      gpFile << regXMin << " " << regYMin << " "
             << _demandMap[i][j].first * _tileWidth +
                    _demandMap[i][j].second * _tileHeight
             << endl;
      gpFile << regXMax << " " << regYMin << " "
             << _demandMap[i][j].first * _tileWidth +
                    _demandMap[i][j].second * _tileHeight
             << endl;
      gpFile << regXMax << " " << regYMax << " "
             << _demandMap[i][j].first * _tileWidth +
                    _demandMap[i][j].second * _tileHeight
             << endl;
      gpFile << regXMin << " " << regYMax << " "
             << _demandMap[i][j].first * _tileWidth +
                    _demandMap[i][j].second * _tileHeight
             << endl;
      gpFile << regXMin << " " << regYMin << " "
             << _demandMap[i][j].first * _tileWidth +
                    _demandMap[i][j].second * _tileHeight
             << endl
             << endl
             << endl;
    }
  }
  gpFileX << "EOF" << endl;
  gpFileX << "pause -1 'Press any key' " << endl;
  gpFileX.close();

  gpFileY << "EOF" << endl;
  gpFileY << "pause -1 'Press any key' " << endl;
  gpFileY.close();

  gpFile << "EOF" << endl;
  gpFile << "pause -1 'Press any key' " << endl;
  gpFile.close();
}

// plot congestion map in matlab format
void ISPD04CongMap::plotMatlab(const std::string& baseFileName) const {
  std::string fnameX = baseFileName.c_str();
  fnameX += "X.m";
  std::string fnameY = baseFileName.c_str();
  fnameY += "Y.m";
  std::string fname = baseFileName.c_str();
  fname += ".m";

  cout << "Saving " << fname << ", " << fnameX << " and " << fnameY << endl;

  ofstream gpFileX(fnameX.c_str());
  ofstream gpFileY(fnameY.c_str());
  ofstream gpFile(fname.c_str());

  gpFile << "X = [ " << endl;
  gpFileX << "X = [ " << endl;
  gpFileY << "X = [ " << endl;
  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    double regXMin = _layoutBBox.xMin + static_cast<double>(i) * _tileWidth;
    double regXMax = _layoutBBox.xMin + static_cast<double>(i + 1) * _tileWidth;
    gpFile << 0.5 * (regXMin + regXMax) << " ; " << endl;
    gpFileX << 0.5 * (regXMin + regXMax) << " ; " << endl;
    gpFileY << 0.5 * (regXMin + regXMax) << " ; " << endl;
  }
  gpFile << "] ;" << endl << endl;
  gpFileX << "] ;" << endl << endl;
  gpFileY << "] ;" << endl << endl;

  gpFile << "Y = [ " << endl;
  gpFileX << "Y = [ " << endl;
  gpFileY << "Y = [ " << endl;
  for (unsigned i = 0; i < _numVertTiles; ++i) {
    double regYMin = _layoutBBox.yMin + static_cast<double>(i) * _tileHeight;
    double regYMax = _layoutBBox.xMin + static_cast<double>(i + 1) * _tileWidth;
    gpFile << 0.5 * (regYMin + regYMax) << " ; " << endl;
    gpFileX << 0.5 * (regYMin + regYMax) << " ; " << endl;
    gpFileY << 0.5 * (regYMin + regYMax) << " ; " << endl;
  }
  gpFile << "] ;" << endl << endl;
  gpFileX << "] ;" << endl << endl;
  gpFileY << "] ;" << endl << endl;

  gpFile << "Z = [" << endl << endl;
  gpFileX << "Z = [" << endl << endl;
  gpFileY << "Z = [" << endl << endl;
  for (unsigned j = 0; j < _numVertTiles; ++j) {
    for (unsigned i = 0; i < _numHorizTiles; ++i) {
      gpFile << _demandMap[i][j].first * _tileWidth +
                    _demandMap[i][j].second * _tileHeight
             << " ";
      gpFileX << _demandMap[i][j].first * _tileWidth << " ";
      gpFileY << _demandMap[i][j].second * _tileHeight << " ";
    }
    gpFile << " ; " << endl;
    gpFileX << " ; " << endl;
    gpFileY << " ; " << endl;
  }
  gpFile << "] ;" << endl << endl;
  gpFileX << "] ;" << endl << endl;
  gpFileY << "] ;" << endl << endl;

  gpFile << "surf(X,Y,Z)" << endl << endl;
  gpFileX << "surf(X,Y,Z)" << endl << endl;
  gpFileY << "surf(X,Y,Z)" << endl << endl;
  gpFile << "colormap hot" << endl;
  gpFileX << "colormap hot" << endl;
  gpFileY << "colormap hot" << endl;

  gpFile << "xlabel('X Axis')" << endl;
  gpFile << "ylabel('Y Axis')" << endl;
  gpFile << "zlabel('Congestion')" << endl;
  gpFile << "title('Total Congestion')" << endl;

  gpFileX << "xlabel('X Axis')" << endl;
  gpFileX << "ylabel('Y Axis')" << endl;
  gpFileX << "zlabel('XCongestion')" << endl;
  gpFileX << "title('Horizontal Congestion')" << endl;

  gpFileY << "xlabel('X Axis')" << endl;
  gpFileY << "ylabel('Y Axis')" << endl;
  gpFileY << "zlabel('YCongestion')" << endl;
  gpFileY << "title('Vertical Congestion')" << endl;

  gpFile << "view(0,90)" << endl;
  gpFileX << "view(0,90)" << endl;
  gpFileY << "view(0,90)" << endl;

  gpFile.close();
  gpFileX.close();
  gpFileY.close();
}

void ISPD04CongMap::plotXPM(const std::string& baseFileName) const {
  plotXPM(baseFileName, 0.0);
}

// plot congestion map in XPM format
void ISPD04CongMap::plotXPM(const std::string& baseFileName,
                            double maxCongestion) const {
  unsigned numColors = 64;
  Palette palette = getPalette();
  std::vector<string> ColorMap = palette.first;
  std::vector<string> Colors = palette.second;

  std::string fnameX = baseFileName.c_str();
  fnameX += "X.xpm";
  std::string fnameY = baseFileName.c_str();
  fnameY += "Y.xpm";
  std::string fname = baseFileName.c_str();
  fname += ".xpm";

  cout << "Saving " << fname << ", " << fnameX << " and " << fnameY << endl;

  ofstream xpmFileX(fnameX.c_str());
  ofstream xpmFileY(fnameY.c_str());
  ofstream xpmFile(fname.c_str());

  unsigned xBlowFactor =
      unsigned(ceil(400.0 / static_cast<double>(_numHorizTiles)));
  unsigned yBlowFactor =
      unsigned(ceil(400.0 / static_cast<double>(_numVertTiles)));

  xpmFile << "/* XPM */" << endl;
  xpmFile << "static char *congestion[] = {" << endl;
  xpmFile << "/* columns rows colors chars-per-pixel */" << endl;
  xpmFile << "\"" << _numHorizTiles * xBlowFactor << " "
          << _numVertTiles * yBlowFactor << " " << numColors << " 1\"," << endl;
  for (unsigned i = 0; i < numColors; ++i)
    xpmFile << "\"" << ColorMap[i] << "\"," << endl;
  xpmFile << "/* pixels */" << endl;

  xpmFileX << "/* XPM */" << endl;
  xpmFileX << "static char *congestionX[] = {" << endl;
  xpmFileX << "/* columns rows colors chars-per-pixel */" << endl;
  xpmFileX << "\"" << _numHorizTiles * xBlowFactor << " "
           << _numVertTiles * yBlowFactor << " " << numColors << " 1\","
           << endl;
  for (unsigned i = 0; i < numColors; ++i)
    xpmFileX << "\"" << ColorMap[i] << "\"," << endl;
  xpmFileX << "/* pixels */" << endl;

  xpmFileY << "/* XPM */" << endl;
  xpmFileY << "static char *congestionY[] = {" << endl;
  xpmFileY << "/* columns rows colors chars-per-pixel */" << endl;
  xpmFileY << "\"" << _numHorizTiles * xBlowFactor << " "
           << _numVertTiles * yBlowFactor << " " << numColors << " 1\","
           << endl;
  for (unsigned i = 0; i < numColors; ++i)
    xpmFileY << "\"" << ColorMap[i] << "\"," << endl;
  xpmFileY << "/* pixels */" << endl;

  vector<vector<double> > image;
  vector<vector<double> > imageX;
  vector<vector<double> > imageY;
  vector<pair<unsigned, unsigned> > maxLocs;
  for (unsigned rj = 1; rj <= _numVertTiles; ++rj) {
    unsigned j = _numVertTiles - rj;
    vector<double> horizLine;
    vector<double> horizLineX;
    vector<double> horizLineY;
    for (unsigned i = 0; i < _numHorizTiles; ++i) {
      double congestion = _demandMap[i][j].first * _tileWidth +
                          _demandMap[i][j].second * _tileHeight;
      if (congestion > maxCongestion) {
        maxLocs.clear();
        maxLocs.push_back(make_pair(i, rj - 1));
        maxCongestion = congestion;
      } else if (congestion == maxCongestion) {
        maxLocs.push_back(make_pair(i, rj - 1));
      }
      horizLine.push_back(congestion);
      horizLineX.push_back(_demandMap[i][j].first * _tileWidth);
      horizLineY.push_back(_demandMap[i][j].second * _tileHeight);
    }
    image.push_back(horizLine);
    imageX.push_back(horizLineX);
    imageY.push_back(horizLineY);
  }

  /*
    for(unsigned radius = 30; radius <= 34; ++radius)
    {
      for(unsigned i = 0; i < maxLocs.size(); ++i)
      {
        unsigned x = maxLocs[i].first;
        unsigned y = maxLocs[i].second;

        // draw left side
        if(x >= radius)
        {
          unsigned j = y;
          for(unsigned k = 0; k <= radius; ++k)
          {
            image[j][x-radius] = maxCongestion;
            if(j == 0) break;
            --j;
          }
          j = y;
          for(unsigned k = 0; k <= radius; ++k)
          {
            image[j][x-radius] = maxCongestion;
            ++j;
            if(j == _numVertTiles) break;
          }
        }

        // draw right side
        if(_numHorizTiles-x > radius)
        {
          unsigned j = y;
          for(unsigned k = 0; k <= radius; ++k)
          {
            image[j][x+radius] = maxCongestion;
            if(j == 0) break;
            --j;
          }
          j = y;
          for(unsigned k = 0; k <= radius; ++k)
          {
            image[j][x+radius] = maxCongestion;
            ++j;
            if(j == _numVertTiles) break;
          }
        }

        // draw bottom
        if(y >= radius)
        {
          unsigned j = x;
          for(unsigned k = 0; k <= radius; ++k)
          {
            image[y-radius][j] = maxCongestion;
            if(j == 0) break;
            --j;
          }
          j = x;
          for(unsigned k = 0; k <= radius; ++k)
          {
            image[y-radius][j] = maxCongestion;
            ++j;
            if(j == _numHorizTiles) break;
          }
        }

        // draw top
        if(_numVertTiles-y > radius)
        {
          unsigned j = x;
          for(unsigned k = 0; k <= radius; ++k)
          {
            image[y+radius][j] = maxCongestion;
            if(j == 0) break;
            --j;
          }
          j = x;
          for(unsigned k = 0; k <= radius; ++k)
          {
            image[y+radius][j] = maxCongestion;
            ++j;
            if(j == _numHorizTiles) break;
          }
        }
      }
    }
  */

  cout << "The maximum congestion is " << maxCongestion << endl;

  double conversion = (numColors - 1) / maxCongestion;
  for (unsigned j = 0; j < _numVertTiles; ++j) {
    for (unsigned k = 0; k < yBlowFactor; ++k) {
      xpmFile << "\"";
      xpmFileX << "\"";
      xpmFileY << "\"";
      for (unsigned i = 0; i < _numHorizTiles; ++i) {
        unsigned index = unsigned(floor(image[j][i] * conversion));
        unsigned indexX = unsigned(floor(imageX[j][i] * conversion));
        unsigned indexY = unsigned(floor(imageY[j][i] * conversion));

        for (unsigned l = 0; l < xBlowFactor; ++l) {
          xpmFile << Colors[index];
          xpmFileX << Colors[indexX];
          xpmFileY << Colors[indexY];
        }
      }
      if (j == _numVertTiles - 1 && k == yBlowFactor - 1) {
        xpmFile << "\"" << endl;
        xpmFileX << "\"" << endl;
        xpmFileY << "\"" << endl;
      } else {
        xpmFile << "\"," << endl;
        xpmFileX << "\"," << endl;
        xpmFileY << "\"," << endl;
      }
    }
  }
  xpmFile << "};" << endl;
  xpmFileX << "};" << endl;
  xpmFileY << "};" << endl;

  xpmFile.close();
  xpmFileX.close();
  xpmFileY.close();
}

void ISPD04CongMap::printDemandMap(void) {
  for (unsigned rj = 1; rj <= _numVertTiles; ++rj) {
    unsigned j = _numVertTiles - rj;
    for (unsigned i = 0; i < _numHorizTiles; ++i) {
      cout << "(" << _demandMap[i][j].first << "," << _demandMap[i][j].second
           << ") ";
    }
    cout << endl;
  }
}

void ISPD04CongMap::printSupplyMap(void) {
  for (unsigned rj = 1; rj <= _numVertTiles; ++rj) {
    unsigned j = _numVertTiles - rj;
    for (unsigned i = 0; i < _numHorizTiles; ++i) {
      cout << "(" << _supplyMap[i][j].first << "," << _supplyMap[i][j].second
           << ") ";
    }
    cout << endl;
  }
}

template <class CellVector>
void ISPD04CongMap::updateDemandForAttachedNets(const HGraphWDimensions &hg,
                                                const CellVector &cells,
                                                double alpha) {
  _netbitboard.clear();

  for (unsigned i = 0; i < cells.size(); ++i) {
    const HGFNode &node = hg.getNodeByIdx(cells[i]);

    for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); ++e) {
      if (_netbitboard.isBitSet((*e)->getIndex())) continue;

      _netbitboard.setBit((*e)->getIndex());

      addNet(hg, (*e)->getIndex(), alpha);
    }
  }
}

void ISPD04CongMap::addDemandForAttachedNets(const HGraphWDimensions &hg,
                                             const vector<unsigned> &cells) {
  updateDemandForAttachedNets(hg, cells, 1.);
}

void ISPD04CongMap::removeDemandForAttachedNets(const HGraphWDimensions &hg,
                                                const vector<unsigned> &cells) {
  updateDemandForAttachedNets(hg, cells, -1.);
}

double ISPD04CongMap::getTotalOverfill(void) const {
  double overfill = 0.;

  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    for (unsigned j = 0; j < _numVertTiles; ++j) {
      overfill +=
          max(0., _demandMap[i][j].first - _supplyMap[i][j].first) * _tileWidth;
      overfill += max(0., _demandMap[i][j].second - _supplyMap[i][j].second) *
                  _tileHeight;
    }
  }

  return overfill;
}

double ISPD04CongMap::getMaxOverfill(void) const {
  double overfill = 0.;

  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    for (unsigned j = 0; j < _numVertTiles; ++j) {
      overfill =
          max(overfill,
              max(0., _demandMap[i][j].first - _supplyMap[i][j].first) *
                      _tileWidth +
                  max(0., _demandMap[i][j].second - _supplyMap[i][j].second) *
                      _tileHeight);
    }
  }

  return overfill;
}

unsigned ISPD04CongMap::getNumOverfullTiles(void) const {
  unsigned num = 0;

  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    for (unsigned j = 0; j < _numVertTiles; ++j) {
      if ((_demandMap[i][j].first > _supplyMap[i][j].first) ||
          (_demandMap[i][j].second > _supplyMap[i][j].second))
        ++num;
    }
  }

  return num;
}

double ISPD04CongMap::estimateOverfullnessChangeFromCellMove(
    const HGraphWDimensions &hg, unsigned cellId, const Point &newPos) const {
  // add the usage for the cell's new position and remove for the old
  _netbitboard.clear();

  const HGFNode &node = hg.getNodeByIdx(cellId);
  Point oldPos = _pl[cellId];

  for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); ++e) {
    if (_netbitboard.isBitSet((*e)->getIndex())) continue;

    _netbitboard.setBit((*e)->getIndex());

    // remove the nets old contribution from temp
    addNetToTemp(hg, (*e)->getIndex(), -1., cellId, oldPos);
    // add the nets new contribution to temp
    addNetToTemp(hg, (*e)->getIndex(), 1., cellId, newPos);
  }

  // calculate the change in overfullness
  double change = 0.;

  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    for (unsigned j = 0; j < _numVertTiles; ++j) {
      double oldOverfill =
          max(0., _demandMap[i][j].first - _supplyMap[i][j].first) *
              _tileWidth +
          max(0., _demandMap[i][j].second - _supplyMap[i][j].second) *
              _tileHeight;
      double newOverfill =
          max(0., _demandMap[i][j].first + _tempMap[i][j].first -
                      _supplyMap[i][j].first) *
              _tileWidth +
          max(0., _demandMap[i][j].second + _tempMap[i][j].second -
                      _supplyMap[i][j].second) *
              _tileHeight;

      change += (newOverfill - oldOverfill);

      // clear the temp map
      _tempMap[i][j] = make_pair(0., 0.);
    }
  }

  return change;
}
