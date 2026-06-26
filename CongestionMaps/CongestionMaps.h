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
#ifndef _CONGESTIONMAPS_
#define _CONGESTIONMAPS_

#include <vector>

#include "RBPlace/RBPlacement.h"

struct ctrSortCongMap {
  unsigned horizIdx;
  unsigned vertIdx;
  double totalCongestion;
};

struct congMap_less_mag {
  bool operator()(const ctrSortCongMap& elem1, const ctrSortCongMap& elem2) {
    return (elem1.totalCongestion < elem2.totalCongestion);
  }
};

class CongestionMaps {
 private:
  const RBPlacement* _rbplace;
  bool _initialized;
  double _horizGridSize;
  double _vertGridSize;
  unsigned _maxHorizIdx;
  unsigned _maxVertIdx;
  BBox _layoutBBox;

  std::vector<std::vector<double> > _F;  // matrix for probabilistic analysis
                                        // needs to be initialised explicity
  std::vector<std::vector<double> > _congMapX;
  std::vector<std::vector<double> > _congMapY;

  template <class CellIds>
  void _remCongByCells(const CellIds& cellIds, bool useProb);

  template <class CellIds>
  void _addCongByCells(const CellIds& cellIds, bool useProb);

 public:
  std::vector<ctrSortCongMap> sortCongBins;

  CongestionMaps()
      : _rbplace(NULL),
        _initialized(false),
        _horizGridSize(0),
        _vertGridSize(0),
        _maxHorizIdx(0),
        _maxVertIdx(0) {}

  CongestionMaps(RBPlacement* rbplace, double gridSize = -9999);

  ~CongestionMaps();

  void clearCongestionMap();
  // statistical based congestion map
  void createCongestionMap();
  void update1NetCong(itHGFEdgeGlobal e, bool addRemove);
  // add or remove congestion contributed by one net

  // following functions use more accurate probabilistic based congestion
  // estimation. from "Estimating Routing Congestion using Probabilistic
  // analysis", J. Lou
  void initializeFMatrix();
  void computePMatrix(std::vector<std::vector<double> >& Px,
                      std::vector<std::vector<double> >& Py, Point& ptLeft,
                      Point& ptRight);
  void createCongestionMapProb();
  void update1NetCongProb(itHGFEdgeGlobal e, bool addRemove);

  // incremental updates to congestion map. Use first remove and then add
  void remCongByCells(const std::vector<unsigned>& cellIds,
                      bool useProb = false);
  // remove the congestion contributed by these cells
  // call this before updating the RBPlace
  void addCongByCells(const std::vector<unsigned>& cellIds,
                      bool useProb = false);
  // add congestion contributed by these cells
  // call this after updating RBPlace and after calling remCongByCells
  // sort the bins according to total congestion
  void initializeSortBins();
  void sortCongestionBins();
  const std::vector<ctrSortCongMap>& getSortedCongestionBins() {
    return sortCongBins;
  }

  void printCongestionStats();

  // plotting utilities
  void plotCongMap(const char* baseFileName);
  // plot congestion map in gnuplot format

  void plotMatCongMap(const char* baseFileName);
  // plot congestion map in matlab format

  void plotXPMCongMap(const char* baseFileName);
  // plot congestion map in XPM format

  double getHorizGridSize() { return _horizGridSize; }
  double getVertGridSize() { return _vertGridSize; }

  double getHorizCongestion(Point Pt);
  double getVertCongestion(Point Pt);
  double getCongestion(Point Pt);
  double getCongestion(const BBox& bbox);

  unsigned getAvgNumCellsPerGrid();
};

#endif
