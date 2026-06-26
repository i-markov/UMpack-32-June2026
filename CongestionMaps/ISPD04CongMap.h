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
#ifndef _ISPD04CONGMAP_H_
#define _ISPD04CONGMAP_H_

#include <string>
#include <utility>
#include <vector>

#include "Ctainers/bitBoard.h"
#include "Placement/placeWOri.h"
#include "RBPlace/RBPlacement.h"
#include "baseGeneric2DResourceMap.h"

class ISPD04CongMap
    : public Generic2DResourceMap<std::pair<double, double>, DoublePairDistT,
                                  DoublePairZeroT> {
 public:
  ISPD04CongMap(const RBPlacement &rbplace, const PlacementWOrient &pl,
                unsigned numHorizTiles, unsigned numVertTiles,
                bool useSteiner = false);

  double getXDemand(const BBox &bbox) const;
  double getYDemand(const BBox &bbox) const;
  double getMaxXDemand(const BBox &bbox) const;
  std::pair<double, double> getDemand(const BBox &bbox) const;
  std::pair<double, double> getDemand(unsigned xTile, unsigned yTile) const;

  double getXSupply(const BBox &bbox) const;
  double getYSupply(const BBox &bbox) const;
  std::pair<double, double> getSupply(const BBox &bbox) const;
  std::pair<double, double> getSupply(unsigned xTile, unsigned yTile) const;

  // the following functions use congestion prediction techniques from
  // J.Westra, C.Bartels and P.Groeneveld, "Probabilistic Congestion
  // Prediction", ISPD 2004.
  void calculateDemand(const HGraphWDimensions &hg);

  void compute(void) {}

  void printDemandMap(void);
  void printSupplyMap(void);

  double getTotalOverfill(void) const;
  double getMaxOverfill(void) const;
  unsigned getNumOverfullTiles(void) const;
  double estimateOverfullnessChangeFromCellMove(const HGraphWDimensions &hg,
                                                unsigned cellId,
                                                const Point &newPos) const;

  void addDemandForAttachedNets(const HGraphWDimensions &hg,
                                const std::vector<unsigned> &cells);
  void removeDemandForAttachedNets(const HGraphWDimensions &hg,
                                   const std::vector<unsigned> &cells);

  // plot congestion map in various formats
  void plot(const std::string& baseFileName) const;
  void plotMatlab(const std::string& baseFileName) const;
  void plotXPM(const std::string& baseFileName) const;
  void plotXPM(const std::string& baseFileName, double maxCongestion) const;

  void setObstacles(std::vector<BBox> *obs) { _obstacles = obs; }

 private:  // functions
  std::vector<Point> *getPointsOnNet(const HGraphWDimensions &hg,
                                      unsigned netId) const;
  std::vector<std::pair<Point, Point> > *getPinPairs(
      const HGraphWDimensions &hg, unsigned netId,
      std::vector<Point> *suppliedPoints = NULL) const;
  void buildSupplyMap(const RoutingResources &rr);
  void resetDemand(void);
  void getGridCoords(const BBox &box, unsigned &minXCell, unsigned &minYCell,
                     unsigned &maxXCell, unsigned &maxYCell) const;
  void addNet(const HGraphWDimensions &hg, unsigned netId, double alpha);
  void addNetToTemp(const HGraphWDimensions &hg, unsigned netId, double alpha,
                    unsigned movingCell, const Point &newPos) const;
  template <class CellVector>
  void updateDemandForAttachedNets(const HGraphWDimensions &hg,
                                   const CellVector &cells, double alpha);
  void addShortNet(
      std::vector<std::vector<std::pair<double, double> > > &theMap,
      const BBox &twoPinNet, double alpha, unsigned xCell,
      unsigned yCell) const;
  void addFlatNet(
      std::vector<std::vector<std::pair<double, double> > > &theMap,
      const BBox &twoPinNet, double alpha, unsigned minXCell, unsigned minYCell,
      unsigned maxXCell, unsigned maxYCell) const;
  void addLUsage(
      std::vector<std::vector<std::pair<double, double> > > &theMap,
      const BBox &twoPinNet, double alpha, unsigned minXCell, unsigned minYCell,
      unsigned maxXCell, unsigned maxYCell) const;
  void addZUsage(
      std::vector<std::vector<std::pair<double, double> > > &theMap,
      const BBox &twoPinNet, double alpha, bool flipped, unsigned minXCell,
      unsigned minYCell, unsigned maxXCell, unsigned maxYCell) const;

 private:  // data
  mutable ResourceMap _tempMap;
  mutable BitBoard _netbitboard;
  const double _alpha;
  const PlacementWOrient &_pl;
  const std::vector<BBox> *_obstacles;
  bool _useSteiner;
};

std::ostream &operator<<(std::ostream &os, std::pair<double, double> d);

#endif
