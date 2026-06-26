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

#include "CongestionMaps.h"

#include "GeomTrees/FastSteiner/mst2.h"

#include <string>

#define NumNetMaxCost 51

using std::cout;
using std::endl;
using std::make_pair;
using std::ofstream;
using std::pair;
using std::string;
using std::vector;

double NetCostFactorTbl[NumNetMaxCost] = {
    0.0,      1.000000, 1.000000, 1.000000, 1.082797, 1.153598, 1.220592,
    1.282322, 1.338501, 1.399094, 1.449260, 1.497380, 1.545500, 1.593619,
    1.641739, 1.689859, 1.730376, 1.770893, 1.811410, 1.851927, 1.892444,
    1.928814, 1.965184, 2.001553, 2.037923, 2.074293, 2.106117, 2.137941,
    2.169766, 2.201590, 2.233414, 2.264622, 2.295830, 2.327038, 2.358246,
    2.389454, 2.418676, 2.447898, 2.477119, 2.506341, 2.535563, 2.560954,
    2.586346, 2.611737, 2.637129, 2.662520, 2.688684, 2.714848, 2.741013,
    2.767177, 2.793341};

static string congestionFileName(const char* baseFileName, const char* suffix) {
  return string(baseFileName) + suffix;
}

vector<pair<Point, Point> >* getMSTEdges(const vector<Point>& pointsOnNet) {
  // call MST to get pin pairs
  vector<fastSteiner::Point> points(pointsOnNet.size());
  vector<long> parent(pointsOnNet.size());

  for (unsigned i = 0; i < pointsOnNet.size(); ++i) {
    points[i].x = pointsOnNet[i].x;
    points[i].y = pointsOnNet[i].y;
  }

  fastSteiner::mst2(pointsOnNet.size(), &points[0], &parent[0]);

  vector<pair<Point, Point> >* pinPairs = new vector<pair<Point, Point> >();
  pinPairs->reserve(pointsOnNet.size() - 1);

  for (unsigned long i = 0; i < pointsOnNet.size(); ++i) {
    if (static_cast<unsigned long>(parent[i]) == i) continue;

    pinPairs->push_back(make_pair(pointsOnNet[i], pointsOnNet[parent[i]]));
  }

  abkfatal(pinPairs->size() == pointsOnNet.size() - 1,
           "Bad number of pin pairs from MST");

  fastSteiner::mst2_package_done();

  return pinPairs;
}

CongestionMaps::CongestionMaps(RBPlacement* rbplace, double gridSize)
    : _rbplace(rbplace), _initialized(false) {
  unsigned numCells = _rbplace->getNumCells();
  _layoutBBox = _rbplace->getBBox();
  double xSize = _layoutBBox.xMax - _layoutBBox.xMin;
  double ySize = _layoutBBox.yMax - _layoutBBox.yMin;
  /*
  double layoutAR = xSize/ySize;
  double numXCells = sqrt(numCells*layoutAR);
  double numYCells = numCells/numXCells;
  _horizGridSize = 2*xSize/numXCells;    //try to have 4 cells in each grid box
  _vertGridSize = 2*ySize/numYCells;
  */
  if (gridSize == -9999) {
    _horizGridSize =
        2 * sqrt(xSize * ySize / numCells);  // try to have 4 cells in
    _vertGridSize = _horizGridSize;          // each grid box
  } else {
    _horizGridSize = gridSize;
    _vertGridSize = gridSize;
  }

  if (_horizGridSize > xSize) _horizGridSize = xSize / 2;
  if (_vertGridSize > ySize) _vertGridSize = ySize / 2;

  cout << "Congestion map grid size is " << _horizGridSize << " X "
       << _vertGridSize << endl;

  _maxHorizIdx = unsigned(ceil(xSize / _horizGridSize));
  _maxVertIdx = unsigned(ceil(ySize / _vertGridSize));

  _congMapX.resize(_maxVertIdx);
  _congMapY.resize(_maxVertIdx);
  for (unsigned i = 0; i < _maxVertIdx; ++i) {
    _congMapX[i].resize(_maxHorizIdx, 0);
    _congMapY[i].resize(_maxHorizIdx, 0);
  }
  sortCongBins.resize(_maxVertIdx * _maxHorizIdx);
}

CongestionMaps::~CongestionMaps() {
  _congMapX.clear();
  _congMapY.clear();
  _F.clear();
}

void CongestionMaps::clearCongestionMap() {
  for (unsigned i = 0; i < _maxVertIdx; ++i) {
    for (unsigned j = 0; j < _maxHorizIdx; ++j) {
      _congMapX[i][j] = 0;
      _congMapY[i][j] = 0;
    }
  }
}

void CongestionMaps::createCongestionMap() {
  clearCongestionMap();
  const HGraphWDimensions& hgWDims = _rbplace->getHGraph();
  itHGFEdgeGlobal e;
  for (e = hgWDims.edgesBegin(); e != hgWDims.edgesEnd(); ++e)
    update1NetCong(e, 1);
}

void CongestionMaps::plotCongMap(const char* baseFileName) {
  const string fnameX = congestionFileName(baseFileName, "X.plt");
  const string fnameY = congestionFileName(baseFileName, "Y.plt");
  const string fname = congestionFileName(baseFileName, ".plt");
  cout << "Saving " << fnameX << ", " << fnameY << " and " << fname << endl;

  ofstream gpFileX(fnameX.c_str());
  ofstream gpFileY(fnameY.c_str());
  ofstream gpFile(fname.c_str());

  gpFileX << "#Use this file as a script for gnuplot" << endl;
  gpFileX << "#(See http://www.gnuplot.info/ for details)" << endl;
  gpFileX << "#" << TimeStamp() << User() << Platform() << endl << endl;
  gpFileX << "set title ' " << baseFileName
          << "X HPWL= " << _rbplace->evalHPWL(true)
          << ", #Cells= " << _rbplace->getNumCells()
          << ", #Nets= " << _rbplace->getHGraph().getNumEdges() << endl
          << endl;
  gpFileX << "set nokey" << endl;
  gpFileX << "#   Uncomment these two lines starting with \"set\"" << endl;
  gpFileX << "#   to save an EPS file for inclusion into a latex document"
          << endl;
  gpFileX << "# set terminal postscript eps color solid 10" << endl;
  gpFileX << "# set output \"" << baseFileName << ".eps\"" << endl << endl;
  gpFileX << "#   Uncomment these two lines starting with \"set\"" << endl;
  gpFileX << "#   to save a PS file for printing" << endl;
  gpFileX << "# set terminal postscript portrait color solid 8" << endl;
  gpFileX << "# set output \"" << baseFileName << ".ps\"" << endl
          << endl
          << endl;
  gpFileX << " set view 55,45 " << endl;
  gpFileX << endl << endl;
  gpFileX << "splot '-' w l " << endl << endl;

  gpFileY << "#Use this file as a script for gnuplot" << endl;
  gpFileY << "#(See http://www.gnuplot.info/ for details)" << endl;
  gpFileY << "#" << TimeStamp() << User() << Platform() << endl << endl;
  gpFileY << "set title ' " << baseFileName
          << "Y HPWL= " << _rbplace->evalHPWL(true)
          << ", #Cells= " << _rbplace->getNumCells()
          << ", #Nets= " << _rbplace->getHGraph().getNumEdges() << endl
          << endl;
  gpFileY << "set nokey" << endl;
  gpFileY << "#   Uncomment these two lines starting with \"set\"" << endl;
  gpFileY << "#   to save an EPS file for inclusion into a latex document"
          << endl;
  gpFileY << "# set terminal postscript eps color solid 10" << endl;
  gpFileY << "# set output \"" << baseFileName << ".eps\"" << endl << endl;
  gpFileY << "#   Uncomment these two lines starting with \"set\"" << endl;
  gpFileY << "#   to save a PS file for printing" << endl;
  gpFileY << "# set terminal postscript portrait color solid 8" << endl;
  gpFileY << "# set output \"" << baseFileName << ".ps\"" << endl
          << endl
          << endl;
  gpFileY << " set view 55,45 " << endl;
  gpFileY << endl << endl;
  gpFileY << "splot '-' w l " << endl << endl;

  gpFile << "#Use this file as a script for gnuplot" << endl;
  gpFile << "#(See http://www.gnuplot.info/ for details)" << endl;
  gpFile << "#" << TimeStamp() << User() << Platform() << endl << endl;
  gpFile << "set title ' " << baseFileName
         << "HPWL= " << _rbplace->evalHPWL(true)
         << ", #Cells= " << _rbplace->getNumCells()
         << ", #Nets= " << _rbplace->getHGraph().getNumEdges() << endl
         << endl;
  gpFile << "set nokey" << endl;
  gpFile << "#   Uncomment these two lines starting with \"set\"" << endl;
  gpFile << "#   to save an EPS file for inclusion into a latex document"
         << endl;
  gpFile << "# set terminal postscript eps color solid 10" << endl;
  gpFile << "# set output \"" << baseFileName << ".eps\"" << endl << endl;
  gpFile << "#   Uncomment these two lines starting with \"set\"" << endl;
  gpFile << "#   to save a PS file for printing" << endl;
  gpFile << "# set terminal postscript portrait color solid 8" << endl;
  gpFile << "# set output \"" << baseFileName << ".ps\"" << endl
         << endl
         << endl;
  gpFile << " set view 55,45 " << endl;
  gpFile << endl << endl;
  gpFile << "splot '-' w l " << endl << endl;

  for (unsigned i = 0; i < _maxVertIdx; ++i) {
    for (unsigned j = 0; j < _maxHorizIdx; ++j) {
      double regXMin = j * _horizGridSize + _layoutBBox.xMin;
      double regXMax = (j + 1) * _horizGridSize + _layoutBBox.xMin;
      double regYMin = i * _vertGridSize + _layoutBBox.yMin;
      double regYMax = (i + 1) * _vertGridSize + _layoutBBox.yMin;

      if (regXMax > _layoutBBox.xMax) regXMax = _layoutBBox.xMax;
      if (regYMax > _layoutBBox.yMax) regYMax = _layoutBBox.yMax;

      gpFileX << regXMin << " " << regYMin << " " << _congMapX[i][j] << endl;
      gpFileX << regXMax << " " << regYMin << " " << _congMapX[i][j] << endl;
      gpFileX << regXMax << " " << regYMax << " " << _congMapX[i][j] << endl;
      gpFileX << regXMin << " " << regYMax << " " << _congMapX[i][j] << endl;
      gpFileX << regXMin << " " << regYMin << " " << _congMapX[i][j] << endl
              << endl
              << endl;
      // gpFileX<<(regXMin+regXMax)/2<<" "<<(regYMin+regYMax)/2<<"
      // "<<_congMapX[i][j]<<endl;

      gpFileY << regXMin << " " << regYMin << " " << _congMapY[i][j] << endl;
      gpFileY << regXMax << " " << regYMin << " " << _congMapY[i][j] << endl;
      gpFileY << regXMax << " " << regYMax << " " << _congMapY[i][j] << endl;
      gpFileY << regXMin << " " << regYMax << " " << _congMapY[i][j] << endl;
      gpFileY << regXMin << " " << regYMin << " " << _congMapY[i][j] << endl
              << endl
              << endl;
      // gpFileY<<(regXMin+regXMax)/2<<" "<<(regYMin+regYMax)/2<<"
      // "<<_congMapY[i][j]<<endl;

      gpFile << regXMin << " " << regYMin << " "
             << _congMapX[i][j] + _congMapY[i][j] << endl;
      gpFile << regXMax << " " << regYMin << " "
             << _congMapX[i][j] + _congMapY[i][j] << endl;
      gpFile << regXMax << " " << regYMax << " "
             << _congMapX[i][j] + _congMapY[i][j] << endl;
      gpFile << regXMin << " " << regYMax << " "
             << _congMapX[i][j] + _congMapY[i][j] << endl;
      gpFile << regXMin << " " << regYMin << " "
             << _congMapX[i][j] + _congMapY[i][j] << endl
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

void CongestionMaps::plotMatCongMap(const char* baseFileName) {
  const string fnameX = congestionFileName(baseFileName, "X.m");
  const string fnameY = congestionFileName(baseFileName, "Y.m");
  const string fname = congestionFileName(baseFileName, ".m");
  cout << "Saving " << fname << " , " << fnameX << " and " << fnameY << endl;

  ofstream gpFileX(fnameX.c_str());
  ofstream gpFileY(fnameY.c_str());
  ofstream gpFile(fname.c_str());

  gpFile << "X = [ " << endl;
  gpFileX << "X = [ " << endl;
  gpFileY << "X = [ " << endl;
  for (unsigned j = 0; j < _maxHorizIdx; ++j) {
    double regXMin = j * _horizGridSize + _layoutBBox.xMin;
    double regXMax = (j + 1) * _horizGridSize + _layoutBBox.xMin;
    gpFile << (regXMin + regXMax) / 2 << " ; " << endl;
    gpFileX << (regXMin + regXMax) / 2 << " ; " << endl;
    gpFileY << (regXMin + regXMax) / 2 << " ; " << endl;
  }

  gpFile << "] ;" << endl << endl;
  gpFileX << "] ;" << endl << endl;
  gpFileY << "] ;" << endl << endl;

  gpFile << "Y = [ " << endl;
  gpFileX << "Y = [ " << endl;
  gpFileY << "Y = [ " << endl;
  for (unsigned i = 0; i < _maxVertIdx; ++i) {
    double regYMin = i * _vertGridSize + _layoutBBox.yMin;
    double regYMax = (i + 1) * _vertGridSize + _layoutBBox.yMin;
    gpFile << (regYMin + regYMax) / 2 << " ; " << endl;
    gpFileX << (regYMin + regYMax) / 2 << " ; " << endl;
    gpFileY << (regYMin + regYMax) / 2 << " ; " << endl;
  }
  gpFile << "] ;" << endl << endl;
  gpFileX << "] ;" << endl << endl;
  gpFileY << "] ;" << endl << endl;

  gpFile << "Z = [" << endl << endl;
  gpFileX << "Z = [" << endl << endl;
  gpFileY << "Z = [" << endl << endl;
  for (unsigned i = 0; i < _maxVertIdx; ++i) {
    for (unsigned j = 0; j < _maxHorizIdx; ++j) {
      gpFile << _congMapX[i][j] + _congMapY[i][j] << " ";
      gpFileX << _congMapX[i][j] << " ";
      gpFileY << _congMapY[i][j] << " ";
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

double CongestionMaps::getHorizCongestion(Point Pt) {
  int gridXIndex = int(floor((Pt.x - _layoutBBox.xMin) / _horizGridSize));
  int gridYIndex = int(floor((Pt.y - _layoutBBox.yMin) / _vertGridSize));
  if (gridXIndex < 0) {
    abkwarn(0, "Point outside layout ");
    gridXIndex = 0;
  } else if (gridXIndex >= int(_maxHorizIdx)) {
    abkwarn(0, "Point outside layout ");
    gridXIndex = _maxHorizIdx - 1;
  }

  if (gridYIndex < 0) {
    abkwarn(0, "Point outside layout ");
    gridYIndex = 0;
  } else if (gridYIndex >= int(_maxVertIdx)) {
    abkwarn(0, "Point outside layout ");
    gridYIndex = _maxVertIdx - 1;
  }
  return (_congMapX[gridYIndex][gridXIndex]);
}

double CongestionMaps::getVertCongestion(Point Pt) {
  int gridXIndex = int(floor((Pt.x - _layoutBBox.xMin) / _horizGridSize));
  int gridYIndex = int(floor((Pt.y - _layoutBBox.yMin) / _vertGridSize));
  if (gridXIndex < 0) {
    abkwarn(0, "Point outside layout ");
    gridXIndex = 0;
  } else if (gridXIndex >= int(_maxHorizIdx)) {
    abkwarn(0, "Point outside layout ");
    gridXIndex = _maxHorizIdx - 1;
  }

  if (gridYIndex < 0) {
    abkwarn(0, "Point outside layout ");
    gridYIndex = 0;
  } else if (gridYIndex >= int(_maxVertIdx)) {
    abkwarn(0, "Point outside layout ");
    gridYIndex = _maxVertIdx - 1;
  }
  return (_congMapY[gridYIndex][gridXIndex]);
}

double CongestionMaps::getCongestion(Point Pt) {
  double totalCongestion = getHorizCongestion(Pt) + getVertCongestion(Pt);
  return (totalCongestion);
}

template <class CellIds>
void CongestionMaps::_remCongByCells(const CellIds& cellIds, bool useProb) {
  const HGraphWDimensions& hgWDims = _rbplace->getHGraph();
  bit_vector seenNets(hgWDims.getNumEdges(), 0);

  for (unsigned i = 0; i < cellIds.size(); ++i) {
    unsigned cellId = cellIds[i];
    const HGFNode& node = hgWDims.getNodeByIdx(cellId);
    for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); e++) {
      const HGFEdge& edge = (**e);
      unsigned netId = edge.getIndex();
      if (!seenNets[netId]) {
        seenNets[netId] = true;
        if (!useProb)
          update1NetCong(e, 0);
        else
          update1NetCongProb(e, 0);
      }
    }
  }
}

void CongestionMaps::remCongByCells(const vector<unsigned>& cellIds,
                                    bool useProb) {
  _remCongByCells(cellIds, useProb);
}

template <class CellIds>
void CongestionMaps::_addCongByCells(const CellIds& cellIds, bool useProb) {
  const HGraphWDimensions& hgWDims = _rbplace->getHGraph();
  bit_vector seenNets(hgWDims.getNumEdges(), 0);

  for (unsigned i = 0; i < cellIds.size(); ++i) {
    unsigned cellId = cellIds[i];
    const HGFNode& node = hgWDims.getNodeByIdx(cellId);
    for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); e++) {
      const HGFEdge& edge = (**e);
      unsigned netId = edge.getIndex();
      if (!seenNets[netId]) {
        seenNets[netId] = true;
        if (!useProb)
          update1NetCong(e, 1);
        else
          update1NetCongProb(e, 1);
      }
    }
  }
}

void CongestionMaps::addCongByCells(const vector<unsigned>& cellIds,
                                    bool useProb) {
  _addCongByCells(cellIds, useProb);
}

void CongestionMaps::update1NetCong(itHGFEdgeGlobal e, bool addRemove) {
  const HGraphWDimensions& hgWDims = _rbplace->getHGraph();
  const PlacementWOrient& placement(*_rbplace);
  unsigned netDegree;
  BBox netBBox;
  const HGFEdge& edge = (**e);
  unsigned netId = edge.getIndex();
  netDegree = (*e)->getDegree();
  itHGFNodeLocal v;
  unsigned nodeOffset = 0;
  for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); v++, nodeOffset++) {
    unsigned cellId = (*v)->getIndex();
    Point pinOffset = hgWDims.nodesOnEdgePinOffset(nodeOffset, netId,
                                                   placement.getOrient(cellId));

    netBBox += (placement[cellId] + pinOffset);
  }
  // calculate the grid bins intersecting with netBBox and increment
  // the routing resource demand

  int minHorizIdx =
      int(floor((netBBox.xMin - _layoutBBox.xMin) / _horizGridSize));
  int maxHorizIdx =
      int(ceil((netBBox.xMax - _layoutBBox.xMin) / _horizGridSize));
  int minVertIdx =
      int(floor((netBBox.yMin - _layoutBBox.yMin) / _vertGridSize));
  int maxVertIdx = int(ceil((netBBox.yMax - _layoutBBox.yMin) / _vertGridSize));

  if (minHorizIdx < 0) minHorizIdx = 0;
  if (minVertIdx < 0) minVertIdx = 0;
  if (maxHorizIdx > int(_maxHorizIdx)) maxHorizIdx = _maxHorizIdx;
  if (maxVertIdx > int(_maxVertIdx)) maxVertIdx = _maxVertIdx;

  double Y = netBBox.getHeight();
  double X = netBBox.getWidth();

  double q;  // net weighting factor
  if (netDegree < NumNetMaxCost)
    q = NetCostFactorTbl[netDegree];
  else
    q = NetCostFactorTbl[NumNetMaxCost - 1];

  double deltaDX;
  double deltaDY;

  for (int i = minVertIdx; i < maxVertIdx; ++i) {
    for (int j = minHorizIdx; j < maxHorizIdx; ++j) {
      double regXMin = j * _horizGridSize + _layoutBBox.xMin;
      double regXMax = (j + 1) * _horizGridSize + _layoutBBox.xMin;
      double regYMin = i * _vertGridSize + _layoutBBox.yMin;
      double regYMax = (i + 1) * _vertGridSize + _layoutBBox.yMin;

      if (regXMax > _layoutBBox.xMax) regXMax = _layoutBBox.xMax;
      if (regYMax > _layoutBBox.yMax) regYMax = _layoutBBox.yMax;

      BBox region;
      region.add(regXMin, regYMin);
      region.add(regXMax, regYMax);

      BBox intersect = netBBox.intersect(region);
      double w = intersect.getHeight();
      double l = intersect.getWidth();
      double W = region.getHeight();
      double L = region.getWidth();

      if (Y == 0 && X != 0)  // degenerate cases first
      {
        deltaDX = q;
        deltaDY = 0;
      } else if (X == 0 && Y != 0) {
        deltaDX = 0;
        deltaDY = q;
      } else if (X == 0 && Y == 0) {
        deltaDX = 0;
        deltaDY = 0;
      } else {
        deltaDX = q * w * l / (Y * L);
        deltaDY = q * w * l / (X * W);
      }

      if (addRemove) {
        _congMapX[i][j] += deltaDX;
        _congMapY[i][j] += deltaDY;
      } else {
        _congMapX[i][j] -= deltaDX;
        _congMapY[i][j] -= deltaDY;
      }
    }
  }
}

void CongestionMaps::initializeFMatrix() {
  _F.resize(_maxVertIdx);
  for (unsigned i = 0; i < _maxVertIdx; ++i) _F[i].resize(_maxHorizIdx, 0);

  for (unsigned i = 0; i < _maxVertIdx; ++i) _F[i][0] = 1;
  for (unsigned i = 0; i < _maxHorizIdx; ++i) _F[0][i] = 1;

  for (unsigned i = 1; i < _maxVertIdx; ++i)
    for (unsigned j = 1; j < _maxHorizIdx; ++j)
      _F[i][j] = _F[i][j - 1] + _F[i - 1][j];
}

void CongestionMaps::createCongestionMapProb() {
  clearCongestionMap();
  const HGraphWDimensions& hgWDims = _rbplace->getHGraph();
  itHGFEdgeGlobal e;
  for (e = hgWDims.edgesBegin(); e != hgWDims.edgesEnd(); ++e) {
    update1NetCongProb(e, 1);
  }
}

void CongestionMaps::update1NetCongProb(itHGFEdgeGlobal e, bool addRemove) {
  const HGraphWDimensions& hgWDims = _rbplace->getHGraph();
  const PlacementWOrient& placement(*_rbplace);
  unsigned netDegree;

  const HGFEdge& edge = (**e);
  unsigned netId = edge.getIndex();
  netDegree = (*e)->getDegree();
  itHGFNodeLocal v;
  unsigned nodeOffset = 0;

  vector<Point> points;
  for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); v++, nodeOffset++) {
    unsigned cellId = (*v)->getIndex();
    Point pinOffset = hgWDims.nodesOnEdgePinOffset(nodeOffset, netId,
                                                   placement.getOrient(cellId));

    points.push_back((placement[cellId] + pinOffset));
  }

  // calculate the mst of the net and create 2 pin nets
  vector<pair<Point, Point> >* MSTEdges = getMSTEdges(points);

  for (unsigned net = 0; net != MSTEdges->size(); ++net) {
    Point point1 = (*MSTEdges)[net].first;
    Point point2 = (*MSTEdges)[net].second;
    Point ptLeft, ptRight;
    if (point1.x < point2.x) {
      ptLeft = point1;
      ptRight = point2;
    } else {
      ptLeft = point2;
      ptRight = point1;
    }

    if (ptLeft.x > _layoutBBox.xMax)
      ptLeft.x = _layoutBBox.xMax;
    else if (ptLeft.x < _layoutBBox.xMin)
      ptLeft.x = _layoutBBox.xMin;

    if (ptLeft.y > _layoutBBox.yMax)
      ptLeft.y = _layoutBBox.yMax;
    else if (ptLeft.y < _layoutBBox.yMin)
      ptLeft.y = _layoutBBox.yMin;

    if (ptRight.x > _layoutBBox.xMax)
      ptRight.x = _layoutBBox.xMax;
    else if (ptRight.x < _layoutBBox.xMin)
      ptRight.x = _layoutBBox.xMin;

    if (ptRight.y > _layoutBBox.yMax)
      ptRight.y = _layoutBBox.yMax;
    else if (ptRight.y < _layoutBBox.yMin)
      ptRight.y = _layoutBBox.yMin;

    ptLeft.x -= _layoutBBox.xMin;
    ptLeft.y -= _layoutBBox.yMin;
    ptRight.x -= _layoutBBox.xMin;
    ptRight.y -= _layoutBBox.yMin;

    int leftGridX = int(floor(ptLeft.x / _horizGridSize));
    int leftGridY = int(floor(ptLeft.y / _vertGridSize));
    int rightGridX = int(floor(ptRight.x / _horizGridSize));
    int rightGridY = int(floor(ptRight.y / _vertGridSize));

    if (leftGridY > rightGridY) {
      int temp = leftGridY;
      leftGridY = rightGridY;
      rightGridY = temp;
    }

    vector<vector<double> > Px;
    vector<vector<double> > Py;

    computePMatrix(Px, Py, ptLeft, ptRight);

    for (int i = leftGridY; i <= rightGridY; ++i)
      for (int j = leftGridX; j <= rightGridX; ++j) {
        if (addRemove) {
          _congMapX[i][j] += Px[i - leftGridY][j - leftGridX];
          _congMapY[i][j] += Py[i - leftGridY][j - leftGridX];
        } else {
          _congMapX[i][j] -= Px[i - leftGridY][j - leftGridX];
          _congMapY[i][j] -= Py[i - leftGridY][j - leftGridX];
        }
      }
  }

  delete MSTEdges;
}

void CongestionMaps::computePMatrix(vector<vector<double> >& Px,
                                    vector<vector<double> >& Py, Point& ptLeft,
                                    Point& ptRight) {
  int leftGridX = int(floor(ptLeft.x / _horizGridSize));
  int leftGridY = int(floor(ptLeft.y / _vertGridSize));
  int rightGridX = int(floor(ptRight.x / _horizGridSize));
  int rightGridY = int(floor(ptRight.y / _vertGridSize));

  bool flipVert = false;

  int XSize = abs(rightGridX - leftGridX) + 1;
  int YSize = abs(rightGridY - leftGridY) + 1;

  Px.resize(YSize);
  Py.resize(YSize);
  for (int i = 0; i < YSize; ++i) {
    Px[i].resize(XSize, 0);
    Py[i].resize(XSize, 0);
  }

  double dx1, dx2, dy1, dy2;
  if (leftGridY > rightGridY) {
    dx1 = _horizGridSize * (leftGridX + 1) - ptLeft.x;
    dy1 = ptLeft.y - _vertGridSize * (leftGridY);
    dx2 = ptRight.x - _horizGridSize * (rightGridX);
    dy2 = _vertGridSize * (rightGridY + 1) - ptRight.y;
    flipVert = true;
    int temp = leftGridY;
    leftGridY = rightGridY;
    rightGridY = temp;
  } else {
    dx1 = _horizGridSize * (leftGridX + 1) - ptLeft.x;
    dy1 = _vertGridSize * (leftGridY + 1) - ptLeft.y;
    dx2 = ptRight.x - _horizGridSize * (rightGridX);
    dy2 = ptRight.y - _vertGridSize * (rightGridY);
  }
  /*
  dx1 = _horizGridSize;
  dy1 = _vertGridSize;
  dx2 = _horizGridSize;
  dy2 = _vertGridSize;
  */

  int m = YSize - 1;
  int n = XSize - 1;

  if (leftGridX == rightGridX && leftGridY == rightGridY) {
    Px[0][0] = fabs(ptRight.x - ptLeft.x) / _horizGridSize;
    Py[0][0] = fabs(ptRight.y - ptLeft.y) / _vertGridSize;
  } else if (leftGridX == rightGridX) {
    for (int i = 0; i <= m; ++i) {
      Px[i][0] = fabs(dx1 + dx2 - _horizGridSize) / (2 * _horizGridSize);
      if (i == 0)
        Py[i][0] = dy1 / _vertGridSize;
      else if (i == m)
        Py[i][0] = dy2 / _vertGridSize;
      else
        Py[i][0] = 1;
    }
  } else if (leftGridY == rightGridY) {
    for (int i = 0; i <= n; ++i) {
      Py[0][i] = fabs(dy1 + dy2 - _vertGridSize) / (2 * _vertGridSize);
      if (i == 0)
        Px[0][i] = dx1 / _horizGridSize;
      else if (i == n)
        Px[0][i] = dx2 / _horizGridSize;
      else
        Px[0][i] = 1;
    }
  } else  // general case
  {
    for (int i = 0; i <= m; ++i) {
      for (int j = 0; j <= n; ++j) {
        if (i == 0 && j == 0) {
          Px[i][j] = _F[m][n - 1] * dx1 / _horizGridSize;
          Py[i][j] = _F[m - 1][n] * dy1 / _vertGridSize;
        } else if (i == m && j == n) {
          Px[i][j] = _F[m][n - 1] * dx2 / _horizGridSize;
          Py[i][j] = _F[m - 1][n] * dy2 / _vertGridSize;
        } else if (i == 0 && j == n) {
          Px[i][j] = dx2 / _horizGridSize;
          Py[i][j] = dy1 / _vertGridSize;
        } else if (i == m && j == 0) {
          Px[i][j] = dx1 / _horizGridSize;
          Py[i][j] = dy2 / _vertGridSize;
        } else if (j == 0) {
          Px[i][j] = _F[m - i][n - 1] * dx1 / _horizGridSize;
          Py[i][j] = (_F[m - i][n] + _F[m - i - 1][n]) / 2;
        } else if (j == n) {
          Px[i][j] = _F[i][n - 1] * dx2 / _horizGridSize;
          Py[i][j] = (_F[i][n] + _F[i - 1][n]) / 2;
        } else if (i == 0) {
          Px[i][j] = (_F[m][n - j] + _F[m][n - j - 1]) / 2;
          Py[i][j] = _F[m - 1][n - j] * dy1 / _vertGridSize;
        } else if (i == m) {
          Px[i][j] = (_F[m][j] + _F[m][j - 1]) / 2;
          Py[i][j] = _F[m - 1][j] * dy2 / _vertGridSize;
        } else {
          Px[i][j] = (_F[i][j] * _F[m - i][n - j - 1] +
                      _F[i][j - 1] * _F[m - i][n - j]) /
                     2;
          Py[i][j] = (_F[i][j] * _F[m - i - 1][n - j] +
                      _F[i - 1][j] * _F[m - i][n - j]) /
                     2;
        }
        Px[i][j] /= _F[m][n];
        Py[i][j] /= _F[m][n];
      }
    }
    if (flipVert == true) {
      double temp;
      for (int i = 0; i <= int(floor(m / 2.0)); ++i)
        for (int j = 0; j <= n; ++j) {
          temp = Px[i][j];
          Px[i][j] = Px[m - i][j];
          Px[m - i][j] = temp;

          temp = Py[i][j];
          Py[i][j] = Py[m - i][j];
          Py[m - i][j] = temp;
        }
    }
  }
}

void CongestionMaps::initializeSortBins() {
  for (unsigned i = 0; i < _maxVertIdx; ++i) {
    for (unsigned j = 0; j < _maxHorizIdx; ++j) {
      unsigned idx = i * _maxHorizIdx + j;
      sortCongBins[idx].horizIdx = j;
      sortCongBins[idx].vertIdx = i;
      sortCongBins[idx].totalCongestion = _congMapX[i][j] + _congMapY[i][j];
    }
  }
}

void CongestionMaps::sortCongestionBins() {
  initializeSortBins();
  sort(sortCongBins.begin(), sortCongBins.end(), congMap_less_mag());
}

double CongestionMaps::getCongestion(const BBox& bbox) {
  int xMin = int(floor((bbox.xMin - _layoutBBox.xMin) / _horizGridSize));
  int xMax = int(ceil((bbox.xMax - _layoutBBox.xMin) / _horizGridSize));
  int yMin = int(floor((bbox.yMin - _layoutBBox.yMin) / _vertGridSize));
  int yMax = int(ceil((bbox.yMax - _layoutBBox.yMin) / _vertGridSize));

  if (xMin < 0) {
    abkwarn(0, "BBox outside layout ");
    xMin = 0;
  }
  if (xMax > int(_maxHorizIdx)) {
    abkwarn(0, "BBox outside layout ");
    xMax = _maxHorizIdx;
  }
  if (yMin < 0) {
    abkwarn(0, "BBox outside layout ");
    yMin = 0;
  }
  if (yMax > int(_maxVertIdx)) {
    abkwarn(0, "BBox outside layout ");
    yMax = _maxVertIdx;
  }

  double totalCongestion = 0;
  double areaCongBin = _horizGridSize * _vertGridSize;

  for (int i = yMin; i < yMax; ++i) {
    for (int j = xMin; j < xMax; ++j) {
      BBox congBin;
      double x = j * _horizGridSize + _layoutBBox.xMin;
      double y = i * _vertGridSize + _layoutBBox.yMin;
      congBin.add(x, y);
      x = (j + 1) * _horizGridSize + _layoutBBox.xMin;
      y = (i + 1) * _vertGridSize + _layoutBBox.yMin;
      congBin.add(x, y);

      BBox intersectBBox = congBin.intersect(bbox);
      double areaIntersect = intersectBBox.getArea();
      totalCongestion +=
          (_congMapX[i][j] + _congMapY[i][j]) * areaIntersect / areaCongBin;
    }
  }
  return totalCongestion;
}

unsigned CongestionMaps::getAvgNumCellsPerGrid() {
  unsigned totalCells = _rbplace->getNumCells();
  double layoutArea = _layoutBBox.getArea();
  double gridRgnArea = _horizGridSize * _vertGridSize;
  unsigned avgCells = unsigned(ceil(totalCells * gridRgnArea / layoutArea));
  return avgCells;
}

const char* colorMap[] = {
    /* columns rows colors chars-per-pixel */
    "  c #030000", ". c #0b0000", "X c #140000", "o c #1c0000", "O c #220000",
    "+ c #2c0001", "@ c #340000", "# c #3b0000", "$ c #430000", "% c #4a0000",
    "& c #540000", "* c #5c0000", "= c #630001", "- c #6a0000", "; c #720000",
    ": c #7c0000", "> c #820001", ", c #8b0000", "< c #930000", "1 c #9c0000",
    "2 c #a50000", "3 c #ab0000", "4 c #b30000", "5 c #bb0000", "6 c #c30000",
    "7 c #cb0001", "8 c #d20002", "9 c #db0001", "0 c #e30000", "q c #eb0000",
    "w c #f00001", "e c #fc0101", "r c #ff0b00", "t c #ff1300", "y c #fe1a01",
    "u c #fd2300", "i c #fe2b00", "p c #fe3400", "a c #fe3c00", "s c #fe4201",
    "d c #fd4c00", "f c #fe5401", "g c #fe5c01", "h c #fe6200", "j c #ff6a01",
    "k c #fd7400", "l c #ff7c00", "z c #fe8301", "x c #fe8a00", "c c #ff9201",
    "v c #fd9d00", "b c #ffa401", "n c #feac01", "m c #ffb400", "M c #ffbb03",
    "N c #fdc400", "B c #fdcb00", "V c #fed300", "C c #fedb00", "Z c #fee300",
    "A c #feeb01", "S c #fef400", "D c #fefc00"};

const char* colors[] = {" ", ".", "X", "o", "O", "+", "@", "#", "$", "%", "&",
                        "*", "=", "-", ";", ":", ">", ",", "<", "1", "2", "3",
                        "4", "5", "6", "7", "8", "9", "0", "q", "w", "e", "r",
                        "t", "y", "u", "i", "p", "a", "s", "d", "f", "g", "h",
                        "j", "k", "l", "z", "x", "c", "v", "b", "n", "m", "M",
                        "N", "B", "V", "C", "Z", "A", "S", "D"};

void CongestionMaps::plotXPMCongMap(const char* baseFileName) {
  unsigned numColors = 63;

  /*
  cout<<"static char *colors[] = {"<<endl;
  for(unsigned i=0; i<numColors; ++i)
    {
      cout<<"\""<<colorMap[i][0]<<"\","<<endl;
    }
  cout<<"};"<<endl;
  exit(0);
  */

  const string fnameX = congestionFileName(baseFileName, "X.xpm");
  const string fnameY = congestionFileName(baseFileName, "Y.xpm");
  const string fname = congestionFileName(baseFileName, ".xpm");
  cout << "Saving " << fname << " , " << fnameX << " and " << fnameY << endl;

  ofstream xpmFileX(fnameX.c_str());
  ofstream xpmFileY(fnameY.c_str());
  ofstream xpmFile(fname.c_str());

  unsigned xBlowFactor = unsigned(ceil(400.0 / _maxHorizIdx));
  unsigned yBlowFactor = unsigned(ceil(400.0 / _maxVertIdx));

  xpmFile << "/* XPM */" << endl;
  xpmFile << "static char *congestion[] = {" << endl;
  xpmFile << "/* columns rows colors chars-per-pixel */" << endl;
  xpmFile << "\"" << _maxHorizIdx * xBlowFactor << " "
          << _maxVertIdx * yBlowFactor << " " << numColors << " 1\"," << endl;
  for (unsigned i = 0; i < numColors; ++i)
    xpmFile << "\"" << colorMap[i] << "\"," << endl;
  xpmFile << "/* pixels */" << endl;

  xpmFileX << "/* XPM */" << endl;
  xpmFileX << "static char *congestionX[] = {" << endl;
  xpmFileX << "/* columns rows colors chars-per-pixel */" << endl;
  xpmFileX << "\"" << _maxHorizIdx * xBlowFactor << " "
           << _maxVertIdx * yBlowFactor << " " << numColors << " 1\"," << endl;
  for (unsigned i = 0; i < numColors; ++i)
    xpmFileX << "\"" << colorMap[i] << "\"," << endl;
  xpmFileX << "/* pixels */" << endl;

  xpmFileY << "/* XPM */" << endl;
  xpmFileY << "static char *congestionY[] = {" << endl;
  xpmFileY << "/* columns rows colors chars-per-pixel */" << endl;
  xpmFileY << "\"" << _maxHorizIdx * xBlowFactor << " "
           << _maxVertIdx * yBlowFactor << " " << numColors << " 1\"," << endl;
  for (unsigned i = 0; i < numColors; ++i)
    xpmFileY << "\"" << colorMap[i] << "\"," << endl;
  xpmFileY << "/* pixels */" << endl;

  vector<vector<double> > image;
  vector<vector<double> > imageX;
  vector<vector<double> > imageY;
  double maxCongestion = -DBL_MAX;
  for (int i = _maxVertIdx - 1; i >= 0; --i) {
    vector<double> horizLine;
    vector<double> horizLineX;
    vector<double> horizLineY;
    for (unsigned j = 0; j < _maxHorizIdx; ++j) {
      double congestion = _congMapX[i][j] + _congMapY[i][j];
      if (congestion > maxCongestion) maxCongestion = congestion;
      horizLine.push_back(congestion);
      horizLineX.push_back(_congMapX[i][j]);
      horizLineY.push_back(_congMapY[i][j]);
    }
    image.push_back(horizLine);
    imageX.push_back(horizLineX);
    imageY.push_back(horizLineY);
  }

  double distance = maxCongestion / (numColors - 1);
  for (unsigned i = 0; i < _maxVertIdx; ++i) {
    for (unsigned k = 0; k < yBlowFactor; ++k) {
      xpmFile << "\"";
      xpmFileX << "\"";
      xpmFileY << "\"";
      for (unsigned j = 0; j < _maxHorizIdx; ++j) {
        unsigned index = unsigned(floor(image[i][j] / distance));
        unsigned indexX = unsigned(floor(imageX[i][j] / distance));
        unsigned indexY = unsigned(floor(imageY[i][j] / distance));

        for (unsigned l = 0; l < xBlowFactor; ++l) {
          xpmFile << colors[index];
          xpmFileX << colors[indexX];
          xpmFileY << colors[indexY];
        }
      }
      if (i == _maxVertIdx - 1 && k == yBlowFactor - 1) {
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

  xpmFileX.close();
  xpmFileY.close();
  xpmFile.close();
}

void CongestionMaps::printCongestionStats() {
  sortCongestionBins();
  int numCongBins = _maxVertIdx * _maxHorizIdx;

  cout << "Congestion Map Statistics " << endl;
  cout << "Congestion grid is " << _maxVertIdx << " x " << _maxHorizIdx << " = "
       << numCongBins << " bins" << endl;
  double peakCongestion = sortCongBins[numCongBins - 1].totalCongestion;
  double avgCongestion = 0;
  double medianCongestion = sortCongBins[int(numCongBins / 2)].totalCongestion;

  int i;
  for (i = 0; i < numCongBins; ++i)
    avgCongestion += sortCongBins[i].totalCongestion;

  avgCongestion /= numCongBins;

  int numLevels = 10;
  double minCongestion = sortCongBins[0].totalCongestion;
  double maxCongestion = peakCongestion;
  double rangeLevel = (maxCongestion - minCongestion) / numLevels;
  vector<int> binsPerLevel(numLevels, 0);

  double currLevel = minCongestion + rangeLevel;
  int currLevelIdx = 0;
  int currIdx = 0;
  while (currIdx < numCongBins) {
    while (currIdx < numCongBins) {
      if (sortCongBins[currIdx].totalCongestion <= currLevel) {
        ++binsPerLevel[currLevelIdx];
        ++currIdx;
      } else
        break;
    }
    if (currLevelIdx < numLevels - 1)  // boundary case because of numerical pbm
      ++currLevelIdx;
    currLevel += rangeLevel;
  }

  for (int i = 0; i < numLevels; ++i) {
    double minLevel = i * 100 / numLevels;
    double maxLevel = (i + 1) * 100 / numLevels;

    double percentBinsInLevel = 100 * binsPerLevel[i] / double(numCongBins);
    cout << "% bins with " << minLevel << "% - " << maxLevel
         << "% congestion : " << percentBinsInLevel << endl;
  }
  cout << endl;
  cout << "Median Congestion : " << medianCongestion << endl;
  cout << "Average Congestion : " << avgCongestion << endl;
  cout << "Peak Congestion : " << peakCongestion << endl;
}
