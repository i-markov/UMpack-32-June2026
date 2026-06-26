/**************************************************************************
***
*** Copyright (c) 2026 IonQ, Inc.
*** Copyright (c) 1995-2000 Regents of the University of California,
***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
*** Copyright (c) 2000-2007 Regents of the University of Michigan,
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

#include "ISPD06DensityMap.h"

#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

#include "ABKCommon/abkassert.h"
#include "RBPlace/RBPlacement.h"

using std::cout;
using std::endl;
using std::map;
using std::max;
using std::min;
using std::numeric_limits;
using std::ofstream;
using std::ostream;
using std::string;
using std::vector;

ISPD06DensityMap::ISPD06DensityMap(const RBPlacement& rbplace,
                                   bool makeMacrosFixed)
    : _rbplace(rbplace),
      _makeMacrosFixed(makeMacrosFixed),
      _total_overflow_amount(0),
      _numTiles(0),
      _total_movable_area(0) {
  setLayoutBBox(rbplace.getBBox(false));
}

void ISPD06DensityMap::compute(void) {
  _total_movable_area = 0;
  double row_height = _makeMacrosFixed
                          ? _rbplace.getCoreRow(0).getSite().height
                          : _rbplace.getBackupCoreRow(0).getSite().height;
  setTileSize(row_height * 10, row_height * 10);
  unsigned numHoriz = static_cast<unsigned>(
               ceil(_layoutBBox.getWidth() / (row_height * 10))),
           numVert = static_cast<unsigned>(
               ceil(_layoutBBox.getHeight() / (row_height * 10)));
  setNumTiles(numHoriz, numVert);
  _numTiles = _numHorizTiles * _numVertTiles;

  for (unsigned i = 0; i < numHoriz; ++i) {
    for (unsigned j = 0; j < numVert; ++j) {
      setSupplyValue(i, j, 0.);
      setDemandValue(i, j, 0.);
    }
  }

  const HGraphWDimensions& hgraph = _rbplace.getHGraph();

  // go over all the rows to compute supply
  if (_makeMacrosFixed) {
    for (unsigned cr = 0; cr < _rbplace.getNumCoreRows(); ++cr) {
      // go over all the sub rows
      for (itRBPSubRow srIt = _rbplace.getCoreRow(cr).subRowsBegin();
           srIt != _rbplace.getCoreRow(cr).subRowsEnd(); ++srIt) {
        double start_x = srIt->getXStart();
        double end_x = srIt->getXEnd();
        double y_coord = srIt->getYCoord();
        BBox subRowBBox(start_x, y_coord, end_x, y_coord + row_height);
        IBBox tiles = getTileCoords(subRowBBox);
        for (int tilei = tiles.xMin; tilei <= tiles.xMax; ++tilei)
          for (int tilej = tiles.yMin; tilej <= tiles.yMax; ++tilej) {
            BBox tileBBox = getTileBBox(tilei, tilej);
            double additional_area = subRowBBox.intersect(tileBBox).getArea();
            double currSupply = getSupply(tilei, tilej);
            setSupplyValue(tilei, tilej, currSupply + additional_area);
          }
      }
    }
  } else {
    for (unsigned cr = 0; cr < _rbplace.getNumBackupCoreRows(); ++cr) {
      // go over all the sub rows
      for (itRBPSubRow srIt = _rbplace.getBackupCoreRow(cr).subRowsBegin();
           srIt != _rbplace.getBackupCoreRow(cr).subRowsEnd(); ++srIt) {
        double start_x = srIt->getXStart();
        double end_x = srIt->getXEnd();
        double y_coord = srIt->getYCoord();
        BBox subRowBBox(start_x, y_coord, end_x, y_coord + row_height);
        IBBox tiles = getTileCoords(subRowBBox);
        for (int tilei = tiles.xMin; tilei <= tiles.xMax; ++tilei)
          for (int tilej = tiles.yMin; tilej <= tiles.yMax; ++tilej) {
            BBox tileBBox = getTileBBox(tilei, tilej);
            double additional_area = subRowBBox.intersect(tileBBox).getArea();
            double currSupply = getSupply(tilei, tilej);
            setSupplyValue(tilei, tilej, currSupply + additional_area);
          }
      }
    }
  }

  // go over all the movable nodes to compute demand
  for (itHGFNodeGlobal nodeIt = hgraph.nodesBegin();
       nodeIt != hgraph.nodesEnd(); ++nodeIt) {
    unsigned nodeId = (*nodeIt)->getIndex();

    if (!(_rbplace.isFixed(nodeId) ||
          (_makeMacrosFixed && hgraph.isNodeMacro(nodeId)))) {
      Point nodeLoc = _rbplace[nodeId];
      BBox nodeBBox(nodeLoc.x, nodeLoc.y,
                    nodeLoc.x + hgraph.getNodeWidth(nodeId),
                    nodeLoc.y + hgraph.getNodeHeight(nodeId));
      IBBox tiles = getTileCoords(nodeBBox);
      for (int tilei = tiles.xMin; tilei <= tiles.xMax; ++tilei) {
        for (int tilej = tiles.yMin; tilej <= tiles.yMax; ++tilej) {
          BBox tileBBox = getTileBBox(tilei, tilej);
          double utilized_area = nodeBBox.intersect(tileBBox).getArea();
          double currDemand = getDemand(tilei, tilej);
          setDemandValue(tilei, tilej, currDemand + utilized_area);
          _total_movable_area += utilized_area;
        }
      }
    } else if (_makeMacrosFixed && hgraph.isNodeMacro(nodeId)) {
      Point nodeLoc = _rbplace[nodeId];
      BBox nodeBBox(nodeLoc.x, nodeLoc.y,
                    nodeLoc.x + hgraph.getNodeWidth(nodeId),
                    nodeLoc.y + hgraph.getNodeHeight(nodeId));
      _total_movable_area += nodeBBox.getArea();
    }
  }

  _total_overflow_amount = 0.;
  unsigned violatingBins = 0;
  // loop over the tiles to get overflow amount
  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    for (unsigned j = 0; j < _numVertTiles; ++j) {
      double overflow_amount =
          getDemand(i, j) - getSupply(i, j) * _targetDensity;
      if (overflow_amount > 0) {
        _total_overflow_amount += overflow_amount;
        ++violatingBins;
      }
    }
  }
}

void ISPD06DensityMap::printHistogramInfo(ostream& os) const {
  // loop over the tiles to get average density
  map<double, unsigned> theMap;
  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    for (unsigned j = 0; j < _numVertTiles; ++j) {
      if (getSupply(i, j) > 0) {
        double tileWS = 1. - getDemand(i, j) / getSupply(i, j);
        if (theMap.find(tileWS) == theMap.end()) {
          theMap[tileWS] = 1;
        } else {
          theMap[tileWS] += 1;
        }
      }
    }
  }
  os << "Histogram Info" << endl;
  unsigned counter = 0;
  const double capInc = 0.025;
  double cap = capInc;
  for (map<double, unsigned>::iterator i = theMap.begin(); i != theMap.end();
       ++i) {
    if (lessThanDouble(i->first, cap)) {
      counter += i->second;
    } else {
      if (counter > 0) os << cap - capInc << " " << counter << endl;
      while (lessOrEqualDouble(cap, i->first)) {
        cap += capInc;
      }
      counter = i->second;
    }
  }
  if (counter > 0) os << cap - capInc << " " << counter << endl;
}

void ISPD06DensityMap::printInfo(ostream& os) const {
  os << "ISPD 2006 contest target density summary: " << endl;

  os << "\tTarget density: " << _targetDensity << endl;
  // printf "\tOverflow per bin: %f\tTotal overflow amount: %f\n",
  // $total_overflow_amount/$total_bin_count, $total_overflow_amount;
  os << "\tOverflow per bin: " << _total_overflow_amount / _numTiles
     << "\tTotal overflow amount: " << _total_overflow_amount << " ("
     << 100. * _total_overflow_amount / _total_movable_area << "%)" << endl;

  os << "\tScaled Overflow per bin: " << getScaledOverflow() << "%" << endl;
}

double ISPD06DensityMap::getScaledOverflow(void) const {
  //$scaled_overflow_per_bin =
  //($total_overflow_amount*$XUNIT*$YUNIT*$target_density)/($total_movable_area
  //* 400);
  //     double scaled_overflow_per_bin = (_total_overflow_amount * _tileWidth *
  //     _tileHeight * _targetDensity)/(_total_movable_area * 400);
  double scaled_overflow_per_bin =
      (_total_overflow_amount * _targetDensity * 36.) / _total_movable_area;
  // printf "\tScaled Overflow per bin: %f\n", $scaled_overflow_per_bin *
  // $scaled_overflow_per_bin;
  return (scaled_overflow_per_bin * scaled_overflow_per_bin);
}

void ISPD06DensityMap::setTargetDensity(double targetDensity) {
  abkassert(targetDensity >= 0 && targetDensity <= 1,
            "Target density must be between 0 and 1");
  _targetDensity = targetDensity;
}

Palette ISPD06DensityMap::getPalette(void) const {
  return PaletteBox::Greyscale();
}

IBBox ISPD06DensityMap::getTileCoords(const BBox& box) const {
  BBox adjustedBox = _layoutBBox.intersect(box);
  Point trans(-_layoutBBox.xMin, -_layoutBBox.yMin);
  adjustedBox.TranslateBy(trans);

  unsigned minXCell = 0, minYCell = 0, maxXCell = 0, maxYCell = 0;
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

  IBBox tileBBox(minXCell, minYCell, maxXCell, maxYCell);
  return tileBBox;
}

BBox ISPD06DensityMap::getGridCell(const Point& p) const {
  Point trans(-_layoutBBox.xMin, -_layoutBBox.yMin);
  Point adjustedPoint = p + trans;

  unsigned xIndex = 0;
  if (adjustedPoint.x <= 0.) {
    xIndex = 0;
  } else {
    xIndex = static_cast<unsigned>(floor(adjustedPoint.x / _tileWidth));
  }
  if (xIndex == _numHorizTiles) --xIndex;

  unsigned yIndex = 0;
  if (adjustedPoint.y <= 0.) {
    yIndex = 0;
  } else {
    yIndex = static_cast<unsigned>(floor(adjustedPoint.y / _tileHeight));
  }
  if (yIndex == _numHorizTiles) --yIndex;

  BBox theBox(_layoutBBox.xMin + static_cast<double>(xIndex) * _tileWidth,
              _layoutBBox.yMin + static_cast<double>(yIndex) * _tileHeight,
              _layoutBBox.xMin + static_cast<double>(xIndex + 1) * _tileWidth,
              _layoutBBox.yMin + static_cast<double>(yIndex + 1) * _tileHeight);

  theBox.xMax = min(theBox.xMax, _layoutBBox.xMax);
  theBox.yMax = min(theBox.yMax, _layoutBBox.yMax);

  return theBox;
}

void ISPD06DensityMap::plotXPM(const std::string& baseFileName) const {
  plotXPM(baseFileName, -numeric_limits<double>::max());
}

void ISPD06DensityMap::plotXPM(const std::string& baseFileName,
                               double maxResourceVal) const {
  unsigned numColors = 64;
  Palette palette = getPalette();
  vector<string> ColorMap = palette.first;
  vector<string> Colors = palette.second;

  std::string fname = baseFileName.c_str();
  fname += ".xpm";

  cout << "Saving " << fname << endl;

  ofstream xpmFile(fname.c_str());

  double maxResourceValue = maxResourceVal;

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

  // ensure that the max resource value is at least as big as actual values
  for (unsigned j = 0; j < _numVertTiles; ++j) {
    for (unsigned i = 0; i < _numHorizTiles; ++i) {
      double consumption = getDemand(i, j);
      if (consumption - maxResourceValue > 0.) {
        maxResourceValue = consumption;
      }
    }
  }

  // compute the overflow image
  vector<vector<unsigned> > image;
  for (unsigned rj = 1; rj <= _numVertTiles; ++rj) {
    unsigned j = _numVertTiles - rj;
    vector<unsigned> horizLine;
    for (unsigned i = 0; i < _numHorizTiles; ++i) {
      double consumption = getDemand(i, j);
      double target = getSupply(i, j) * _targetDensity;
      double conversionLow = (numColors / 2 - 1) / target;
      double conversionHigh = (numColors / 2 - 1) / (maxResourceValue - target);
      unsigned index = 0;
      if (equalDouble(target,
                      0.0))  // to specialize the case where demand==supply==0
        index = numColors / 2;
      else if (consumption <= target)
        index = unsigned(floor(consumption * conversionLow));
      else
        index = unsigned(floor((consumption - target) * conversionHigh)) +
                numColors / 2;

      horizLine.push_back(index);
    }
    image.push_back(horizLine);
  }

  // write the xpm file from image
  for (unsigned j = 0; j < _numVertTiles; ++j) {
    for (unsigned k = 0; k < yBlowFactor; ++k) {
      xpmFile << "\"";
      for (unsigned i = 0; i < _numHorizTiles; ++i) {
        unsigned index = image[j][i];

        for (unsigned l = 0; l < xBlowFactor; ++l) {
          xpmFile << Colors[index];
        }
      }
      if (j == _numVertTiles - 1 && k == yBlowFactor - 1) {
        xpmFile << "\"" << endl;
      } else {
        xpmFile << "\"," << endl;
      }
    }
  }
  xpmFile << "};" << endl;

  xpmFile.close();
}

void ISPD06DensityMap::removeCellUsage(unsigned cellId) {
  if (_rbplace.isFixed(cellId)) return;

  Point nodeLoc = _rbplace.getPlacement()[cellId];
  BBox nodeBBox(nodeLoc.x, nodeLoc.y,
                nodeLoc.x + _rbplace.getHGraph().getNodeWidth(cellId),
                nodeLoc.y + _rbplace.getHGraph().getNodeHeight(cellId));
  IBBox tiles = getTileCoords(nodeBBox);
  for (int tilei = tiles.xMin; tilei <= tiles.xMax; ++tilei) {
    for (int tilej = tiles.yMin; tilej <= tiles.yMax; ++tilej) {
      double currDemand = getDemand(tilei, tilej);
      double currSupplyTarget = getSupply(tilei, tilej) * _targetDensity;
      double old_overflow_amount = max(0., currDemand - currSupplyTarget);
      BBox tileBBox = getTileBBox(tilei, tilej);
      double utilized_area = nodeBBox.intersect(tileBBox).getArea();
      setDemandValue(tilei, tilej, currDemand - utilized_area);
      double new_overflow_amount =
          max(0., currDemand - utilized_area - currSupplyTarget);
      _total_movable_area -= utilized_area;
      _total_overflow_amount += (new_overflow_amount - old_overflow_amount);
    }
  }
}

void ISPD06DensityMap::addCellUsage(unsigned cellId) {
  if (_rbplace.isFixed(cellId)) return;

  Point nodeLoc = _rbplace.getPlacement()[cellId];
  BBox nodeBBox(nodeLoc.x, nodeLoc.y,
                nodeLoc.x + _rbplace.getHGraph().getNodeWidth(cellId),
                nodeLoc.y + _rbplace.getHGraph().getNodeHeight(cellId));
  IBBox tiles = getTileCoords(nodeBBox);
  for (int tilei = tiles.xMin; tilei <= tiles.xMax; ++tilei) {
    for (int tilej = tiles.yMin; tilej <= tiles.yMax; ++tilej) {
      double currDemand = getDemand(tilei, tilej);
      double currSupplyTarget = getSupply(tilei, tilej) * _targetDensity;
      double old_overflow_amount = max(0., currDemand - currSupplyTarget);
      BBox tileBBox = getTileBBox(tilei, tilej);
      double utilized_area = nodeBBox.intersect(tileBBox).getArea();
      setDemandValue(tilei, tilej, currDemand + utilized_area);
      double new_overflow_amount =
          max(0., currDemand + utilized_area - currSupplyTarget);
      _total_movable_area += utilized_area;
      _total_overflow_amount += (new_overflow_amount - old_overflow_amount);
    }
  }
}

double ISPD06DensityMap::estimateOverfullnessChangeFromCellMove(
    unsigned cellId, const Point& newPos) const {
  if (_rbplace.isFixed(cellId)) return 0.;

  Point oldLoc = _rbplace.getPlacement()[cellId];
  BBox oldNodeBBox(oldLoc.x, oldLoc.y,
                   oldLoc.x + _rbplace.getHGraph().getNodeWidth(cellId),
                   oldLoc.y + _rbplace.getHGraph().getNodeHeight(cellId));
  IBBox oldTiles = getTileCoords(oldNodeBBox);
  BBox newNodeBBox(newPos.x, newPos.y,
                   newPos.x + _rbplace.getHGraph().getNodeWidth(cellId),
                   newPos.y + _rbplace.getHGraph().getNodeHeight(cellId));
  IBBox newTiles = getTileCoords(newNodeBBox);

  double change = 0;

  for (int tilei = oldTiles.xMin; tilei <= oldTiles.xMax; ++tilei) {
    for (int tilej = oldTiles.yMin; tilej <= oldTiles.yMax; ++tilej) {
      BBox tileBBox = getTileBBox(tilei, tilej);
      double old_utilized_area = oldNodeBBox.intersect(tileBBox).getArea();
      double new_utilized_area = newNodeBBox.intersect(tileBBox).getArea();
      double currDemand = getDemand(tilei, tilej);
      double supply = getSupply(tilei, tilej) * _targetDensity;
      double newDemand = currDemand - old_utilized_area + new_utilized_area;
      double oldOverfill = max(0., currDemand - supply);
      double newOverfill = max(0., newDemand - supply);
      change += (newOverfill - oldOverfill);
    }
  }

  for (int tilei = newTiles.xMin; tilei <= newTiles.xMax; ++tilei) {
    for (int tilej = newTiles.yMin; tilej <= newTiles.yMax; ++tilej) {
      // don't go over the same tile twice
      if (oldTiles.xMin <= tilei && tilei <= oldTiles.xMax &&
          oldTiles.yMin <= tilej && tilej <= oldTiles.yMax)
        continue;

      BBox tileBBox = getTileBBox(tilei, tilej);
      double old_utilized_area = oldNodeBBox.intersect(tileBBox).getArea();
      double new_utilized_area = newNodeBBox.intersect(tileBBox).getArea();
      double currDemand = getDemand(tilei, tilej);
      double supply = getSupply(tilei, tilej) * _targetDensity;
      double newDemand = currDemand - old_utilized_area + new_utilized_area;
      double oldOverfill = max(0., currDemand - supply);
      double newOverfill = max(0., newDemand - supply);
      change += (newOverfill - oldOverfill);
    }
  }

  return change;
}

double ISPD06DensityMap::getTotalOverfill(void) const {
  double overfill = 0.;

  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    for (unsigned j = 0; j < _numVertTiles; ++j) {
      overfill += max(0., getDemand(i, j) - getSupply(i, j) * _targetDensity);
    }
  }

  return overfill;
}

double ISPD06DensityMap::getMaxOverfill(void) const {
  double maxOverfill = 0.;

  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    for (unsigned j = 0; j < _numVertTiles; ++j) {
      maxOverfill =
          max(maxOverfill, getDemand(i, j) - getSupply(i, j) * _targetDensity);
    }
  }

  return maxOverfill;
}

unsigned ISPD06DensityMap::getNumOverfullTiles(void) const {
  unsigned num = 0;

  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    for (unsigned j = 0; j < _numVertTiles; ++j) {
      if (greaterThanDouble(getDemand(i, j),
                            getSupply(i, j) * _targetDensity)) {
        ++num;
      }
    }
  }

  return num;
}
