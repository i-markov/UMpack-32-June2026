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
// created by David Papa 11/02/06

#ifndef _BASEGENERIC2DRESOURCEMAP_H_
#define _BASEGENERIC2DRESOURCEMAP_H_

#include <map>
#include <string>
#include <vector>

#include "Geoms/bbox.h"
#include "palette.h"

class DoubleDistT {
 public:
  double operator()(const double& a, const double& b) { return a - b; }
};

class DoubleZeroT {
 public:
  double operator()(void) { return 0.0; }
};

class DoublePairDistT {
 public:
  double operator()(const std::pair<double, double>& a,
                    const std::pair<double, double>& b) {
    return sqrt((a.first - b.first) * (a.first - b.first) +
                (a.second - b.second) * (a.second - b.second));
  }
};

class DoublePairZeroT {
 public:
  std::pair<double, double> operator()(void) {
    return std::pair<double, double>(0.0, 0.0);
  }
};

std::ostream& operator<<(std::ostream& os, std::pair<double, double> d);

class BaseGeneric2DResourceMap {
 public:
  virtual ~BaseGeneric2DResourceMap() {}

  // input functions

  // will overwrite current map with a new default constructed one of the
  // specified size
  virtual void setNumTiles(unsigned numHorizTiles, unsigned numVertTiles) = 0;
  virtual void setTileSize(double tileWidth, double tileHeight) = 0;
  virtual void setLayoutBBox(const BBox& layoutBBox) = 0;
  virtual void compute(void) = 0;

  // inspectors
  virtual unsigned getNumHorizTiles(void) const = 0;
  virtual unsigned getNumVertTiles(void) const = 0;
  virtual double getTileWidth(void) const = 0;
  virtual double getTileHeight(void) const = 0;
  virtual BBox getLayoutBBox(void) const = 0;
  virtual BBox getTileBBox(unsigned i, unsigned j) const = 0;

  // Ploting functions in various formats
  virtual void plotGP(const std::string& baseFileName) const = 0;
  virtual void plotMatlab(const std::string& baseFileName) const = 0;
  virtual void plotXPM(const std::string& baseFileName) const = 0;

  // get the color palette to use in plotting
  // a virtual function so that derived classes
  // can change the palette without changing plotting code
  virtual Palette getPalette(void) const = 0;

  // print whatever information sumarizes this resouce map
  virtual void printInfo(std::ostream& os) const {}

  // estimate the change in overfullness caused by moving a cell
  virtual double getTotalOverfill(void) const = 0;
  virtual double getMaxOverfill(void) const = 0;
  virtual unsigned getNumOverfullTiles(void) const = 0;
  virtual double estimateOverfullnessChangeFromCellMove(
      unsigned cellId, const Point& newPos) const = 0;
};

template <typename ResourceT, typename ResourceDistT, typename ResourceZeroT>
class Generic2DResourceMap : public BaseGeneric2DResourceMap {
 public:
  typedef ResourceT ResourceType;
  typedef ResourceDistT ResourceDistType;
  typedef ResourceZeroT ResourceZeroType;
  typedef std::vector<std::vector<ResourceT> > ResourceMap;

  Generic2DResourceMap();
  virtual ~Generic2DResourceMap();

  // input functions

  // will overwrite current map with a new default constructed one of the
  // specified size
  void setNumTiles(unsigned numHorizTiles, unsigned numVertTiles);
  void setTileSize(double tileWidth, double tileHeight);
  void setLayoutBBox(const BBox& layoutBBox);

  void setDemandMap(const ResourceMap& demandMap);
  void setDemandValue(unsigned x, unsigned y, const ResourceT& value);

  void setSupplyMap(const ResourceMap& supplyMap);
  void setSupplyValue(unsigned x, unsigned y, const ResourceT& value);

  void compute(void) {}

  // inspectors
  unsigned getNumHorizTiles(void) const { return _numHorizTiles; }
  unsigned getNumVertTiles(void) const { return _numVertTiles; }
  double getTileWidth(void) const { return _tileWidth; }
  double getTileHeight(void) const { return _tileHeight; }
  BBox getLayoutBBox(void) const { return _layoutBBox; }
  BBox getTileBBox(unsigned i, unsigned j) const;

  ResourceT getDemand(const BBox& bbox) const;
  ResourceT getDemand(unsigned x, unsigned y) const;

  ResourceT getSupply(const BBox& bbox) const;
  ResourceT getSupply(unsigned x, unsigned y) const;

  double getTotalOverfill(void) const;
  double getMaxOverfill(void) const;
  unsigned getNumOverfullTiles(void) const;

  void printDemandMap(void) const;
  void printSupplyMap(void) const;

  // Ploting functions in various formats
  virtual void plotGP(const std::string& baseFileName) const;
  virtual void plotMatlab(const std::string& baseFileName) const;
  virtual void plotXPM(const std::string& baseFileName,
                       const ResourceT& maxResourceValue) const;

  // get the color palette to use in plotting
  // a virtual function so that derived classes
  // can change the palette without changing plotting code
  virtual Palette getPalette(void) const;

  // estimate the change in overfullness caused by moving a cell
  virtual double estimateOverfullnessChangeFromCellMove(
      unsigned cellId, const Point& newPos) const;

 protected:
  ResourceMap _demandMap;
  ResourceMap _supplyMap;

  unsigned _numHorizTiles, _numVertTiles;
  double _tileWidth, _tileHeight;
  BBox _layoutBBox;

 private:
  double getMagnitude(const ResourceT& r) const;
};

#include "baseGeneric2DResourceMap.cxx"

#endif
