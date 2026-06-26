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

// Created 970804 by Igor Markov and Paul Tucker, VLSI CAD ABKGROUP, UCLA/CS
// 970821 pat  added spacing to subRows
// 970823 ilm  multiple minor corrections (removing redundancies, adding refs)
// 970823 ilm  added orientations to RBPCoreRow and RBPIORow (Maogang's request)
//                  they are initialized now
// 970902 ilm  RBPlacement::RBPlacement(const DB&) does not populate rows
//                  with cell Id numbers anymore. This is done with
//                  RBPlacement::populateRowsWith(const Placement&)
//                  (the code has to be exteneded to put IOPads intoIORows,
//                  and a lot more)
// 970918 pat  add DB argument to populateRowsWith
// 970919 pat  added missing field to RBPSubRow copier
// 970930 ilm  added RBPCoreSubrow::getSymmetry() const;
// 970930 pat  added findSubRow and findCoreRow signatures
// 971001 pat  moved comparison functors from *.cxx
// 980222 pat  handle additional wierd cases in subrow construction
// 980222  dv  changed the interface of populateRowsWith, added two more
//             unsigned Data Members, added support for isPopulated.
// 980603  dv  removed orientedRBPlacement class
// 980309 pat  more const declarations; added siteRoundup method
// 980427 pat  added RBPCoreRow::findSubRow
// 990120 aec  complete revision. Merged with RBPwCheckups, changed
//      permissions, modified CoreRows and SubRows, added placement
//      functions, etc.
// 990312 aec  added parameters object
// 990527 aec  split RBPl into two packages. The base does not use DB
// 990608 s.m  added spaceCellsWithPinDimensionsAlg1
// 990911 aec  added file I/O.  changed Placement to PlacementWOrient
// 000803 ilm  fixed savePlacement() to save cell names
// 001128 s.a  added a plotter function
// 001128 s.a  added a function for overlap detection and removal

#ifndef _RBPLACEMENT_H_
#define _RBPLACEMENT_H_

#include <cfloat>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abktypes.h"
#include "ABKCommon/infolines.h"
#include <string>
#include <vector>
#include "ABKCommon/verbosity.h"
#include "CongestionMaps/baseGeneric2DResourceMap.h"
#include "Constraints/constraints.h"
#include "HGraphWDims/hgWDims.h"
#include "Placement/placeWOri.h"
#include "RBCellCoord.h"
#include "RBPlaceAux.h"
#include "RBPlacePlot.h"
#include "RBRows.h"
#include "overlapRem.h"
#include "routingResources.h"

#ifdef _MSC_VER
#pragma warning(disable : 4355)
#else
#include <unistd.h>
#endif

typedef std::vector<RBPCoreRow>::const_iterator itRBPCoreRow;

// This is used as the return type for a lookup member function
typedef std::pair<const RBPCoreRow *, const RBPSubRow *> RBRowPtrs;

enum clusterNodeOrient { N, FN, S, FS, W, FW, E, FE, I };
enum overlapData { All, Movable, OnlyMacro, OnlyFixed };
enum plotMode { plotAll, plotMovable, plotFixed };
enum greedyMoveMode { HPWL, Cong, ISPDContest };

class ISPD06DensityMap;

class RBPlacement {
  friend class HGraphFixed;

 public:
  class Parameters : public HGraphParameters {
   public:
    enum spaceCellsAlgType { EQUAL_SPACE, WITH_PIN_ALG1, WITH_PIN_ALG2 };

    Verbosity verb;
    unsigned numRowsToRemove;  // rbplace will remove this
    // many coreRows from the bottom of
    // the design
    spaceCellsAlgType spaceCellsAlg;  // spaceCells algorithm to use
    // added by sadya
    bool remCongestion;  // used to remove congestion by blowing
    // up cell widths. to be used only on
    //  routed placements
    bool adaptivePlotting;
    bool compressedPlotting;
    unsigned xpixels, ypixels;
    unsigned overlapCells;
    unsigned greedyMoveCells;
    bool macroOverlapDecreasing;
#ifdef USEFLOORIST
    bool enableFloorist;
#endif

    // Parameter added by DAP, to be used in conjunction with
    // Freeshape floorplanning.  It will compute the HPWL of nets
    // not added during shredding with locations set to the center
    // of gravity of the locations of the shreds
    bool centerOfGravityFSFPEval;
    bool pinToPinFSFPEval;

    std::string obstacleFile;
    bool allowSoftBlocks;

    Parameters();
    Parameters(int argc, const char *argv[]);
    Parameters(const Parameters &orig);
    friend std::ostream &operator<<(std::ostream &, const Parameters &);
  };

  Constraints constraints;  // defined in the Constraints package
  // each individual constraint within
  // this object will reference the
  // member "locations" defined above
  // These constraints copied directly from DB.spatial.constraints
  // copied only by RBPlFromDB, else remains empty

  RoutingResources routingResources;
  // These routing resources are populated by RBPlFromDB, else remain empty

 protected:
  HGraphWDimensions *_hgWDims;  // owned by RBPl, but
                                // a pointer so it can be populated by
                                // ctors of derived classes

  std::string _hgWDimsAuxName;  // name of an aux file where we dump the hgraph
                                 // to save memory in a pinch

  bool _ownHGWDimsAuxName;  // do we own the aux file where we dumped
                            // the hgraph (as opposed to reusing
                            // the aux file passed to us in the beginning)
                            // if we own it, we need to unlink it
                            // during destruction

  bool _populated;
  unsigned _cellsNotInRows;  // the #of cells not in a subrow

  std::vector<RBPCoreRow> _coreRows;

  // backUpCore rows is a copy of above. during mixed-size placement a copy
  // of initial is made so that one can revert when needed
  std::vector<RBPCoreRow> _backUpCoreRows;

  std::map<unsigned, std::string> _logicValues;

  PlacementWOrient _placement;  // location of each cell
  bit_vector _isInSubRow;
  std::vector<RBCellCoord> _cellCoords;  // for cells in subRows,
  // cellCoords allows for O(1)
  // time lookups
  bit_vector _isFixed;
  bit_vector _isCoreCell;
  bit_vector _storElt;

  Parameters _params;

  // both of following are initialized by initBBox()
  BBox _coreBBox;  // bbox of core region
  BBox _fullBBox;  // bbox of core region + terminals

  std::string origFileName;
  std::string dirName;
  std::vector<std::string> allFileNamesInAUX;

  void setCellsInRows();
  void sortSubRowsByX();
  void populateCC();
  void matchOrientToRow(unsigned id);
  void matchOrientsToRows();

  double oneNetHPWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                    unsigned netId, bool usePinLocs) const;
  double oneNetHPWL(unsigned netId, bool usePinLocs) const {
    return oneNetHPWL(getPlacement(), getHGraph(), netId, usePinLocs);
  }
  double oneNetMSTWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                     unsigned netId, bool usePinLocs) const;

  std::pair<double, double> oneNetXYHPWL(const PlacementWOrient &pl,
                                         const HGraphWDimensions &hg,
                                         unsigned netId, bool usePinLocs);

  void spaceCellsEqually(RBPSubRow &sr, unsigned firstCellOffset,
                         unsigned numCells);

  void spaceCellsWithPinDimensionsAlg1(RBPSubRow &sr, unsigned firstCellOffset,
                                       unsigned numCells);

  void spaceCellsWithPinDimensionsAlg2(RBPSubRow &sr, unsigned firstCellOffset,
                                       unsigned numCells);

  void constructRows(const char *sclFile);

  void coarsenRows(double newHeight);
  void saveAsSCLcoarsenRows(double newHeight,
                            const std::string &baseFileName) const;

  void parseLogic(const char *logicFile);

  RBPlacement(unsigned numCells, const Placement &pl,
              const std::vector<Orient> &ori, const Parameters &params);

  RBPlacement(unsigned numCells, const PlacementWOrient &pl,
              const Parameters &params);

 public:
  RBPlacement();

  RBPlacement(const char *auxFileName, const Parameters &params);

  // this ctor reads only the netlist from the auxfile.
  // it will create a layout region with the target aspect ratio
  // and whitespace.  It also assumes the following:
  // core node widths are integer
  // core node heights are all the same
  // pad node dimensions are all the same
  // whitespace is a percent (i.e. 10 == 10%)
  //  if aspect ratio is 0, there will be one row
  RBPlacement(const char *auxFile, double aspectRatio, double whiteSpace,
              const Parameters &params);

  ~RBPlacement();

  // added by dp for one step construction
  // Note:  The RBPlacement takes ownership of the HGraphWDimensions
  RBPlacement(HGraphWDimensions *hg, const PlacementWOrient &pl,
              const std::vector<RBPCoreRow> *cr,
              const Parameters &params = Parameters());
  RBPlacement(HGraphWDimensions *hg, const PlacementWOrient &pl,
              const std::vector<RBPCoreRow> &cr,
              const Parameters &params = Parameters());

  // added by royj just for fun!
  RBPlacement(const RBPlacement &orig);

  std::string getOrigFileName() const { return origFileName; }
  std::string getAuxName() const;  // derived from origFileName
  std::string getDirName() const { return dirName; }
  const std::vector<std::string> &getFileNames() const {
    return allFileNamesInAUX;
  }
  const Parameters &getParameters() const { return _params; }

  void resetPlacementWOri(const PlacementWOrient &pl);
  void resetPlacement(const Placement &pl);
  void resetPlacement();

  unsigned getNumCells() const { return _placement.getSize(); }
  bool isPopulated() const { return _populated; }
  bool allCellsPlaced() const { return !_cellsNotInRows; }
  unsigned numCellsNotPlaced() const { return _cellsNotInRows; }

  bool isCoreCell(unsigned id) const { return _isCoreCell[id]; }
  bool isFixed(unsigned id) const { return _isFixed[id]; }
  bool isStorElt(unsigned id) const { return _storElt[id]; }
  bool isInSubRow(unsigned id) const { return _isInSubRow[id]; }

  const bit_vector &getFixed() const { return _isFixed; }
  const bit_vector &getStorElts() const { return _storElt; }

  operator const PlacementWOrient &() const { return _placement; }
  const PlacementWOrient &getPlacement() const { return _placement; }
  // operator const Placement&()    const {return _placement;}

  const Point &operator[](unsigned id) const { return _placement[id]; }
  void setLocation(unsigned id, const Point &pt);
  Orient &getOrient(unsigned id) { return _placement.getOrient(id); }
  const Orient &getOrient(unsigned id) const {
    return _placement.getOrient(id);
  }
  void setOrient(unsigned id, Orient ori) { _placement.setOrient(id, ori); }

  // added by MRG to find the cell coord of a random point
  unsigned findCellIdx(const Point &point);
  // added by MRG
  const RBPCoreRow &getCoreRow(unsigned id) const { return _coreRows[id]; }
  // added by MRG to calculate total overlap
  double calcOverlap();
  double calcInstOverlap(std::vector<unsigned> &movables);

  // by sadya to calculate site area in given BBox
  double getContainedSiteAreaInBBox(const BBox &bbox) const;

  // by sadya to calculate overlaps regardless of rowstructure
  double calcOverlapGeneric(overlapData mode = All, bool print = false) const;
  double calcOverlapInBBox(const PlacementWOrient &pl, const BBox &bbox,
                           overlapData mode) const;
  double calcOverlapSpecific(const PlacementWOrient &pl, const BBox &bbox,
                             const std::vector<unsigned> &cellIds,
                             overlapData mode, bool print) const;
  // by sadya to align cells to closest core row y coordinates
  void alignCellsToRows();

  // check if all core cells have valid CellCoords
  bool checkCC();

  // added by MRG to calculate mean of set of pins
  Point calcMeanLoc(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                    unsigned cellId, bool usePinLocs = false);
  // added by MRG to update placement of an arbitrary set of cells
  void updateCells(const std::vector<unsigned> &movables,
                   const Placement &soln);

  // These return the row and subrow where a point would lie
  // These can't be const..it messes up the equal_range function
  const RBPCoreRow *findCoreRow(const Point &point) const;
  const unsigned findCoreRowIdx(const Point &point) const;
  const RBPSubRow *findSubRow(const Point &point) const;
  RBRowPtrs findBothRows(const Point &point);

  itRBPCoreRow coreRowsBegin() const { return _coreRows.begin(); }
  itRBPCoreRow coreRowsEnd() const { return _coreRows.end(); }
  unsigned getNumCoreRows() const { return _coreRows.size(); }

  bool isCellInRow(unsigned id) const { return _isInSubRow[id]; }

  const RBCellCoord &getCellCoord(unsigned id) const {
    abkassert(_isInSubRow[id], "Requested cellCorrd for cell not in a row");
    return _cellCoords[id];
  }

  // the following functions (un-)place cells from core or sub rows.

  //'extract' functions remove the cells from rows,
  //'embed' functions place the cells in rows.

  // see documentation in the DOCS subdirectory for more information

  void extractCellsFromSR(RBPSubRow &sr, std::vector<unsigned> &cellIds,
                          std::vector<double> &offsets, double beginOffset = 0,
                          double endOffset = DBL_MAX);

  void extractCellsFromSR(RBPSubRow &sr, std::vector<unsigned> &cellIds,
                          std::vector<unsigned> &siteOffsets,
                          unsigned beginSite = 0, unsigned endSite = UINT_MAX);

  void extractCellsFromSR(RBPSubRow &sr, std::vector<unsigned> &cellIds,
                          double beginOffset = 0, double endOffset = DBL_MAX);

  void embedCellsInSR(RBPSubRow &sr, const std::vector<unsigned> &cellIds,
                      const std::vector<double> &offsets);

  void embedCellsInSR(RBPSubRow &sr, const std::vector<unsigned> &cellIds,
                      const std::vector<unsigned> &siteOffsets);

  void embedCellsInSR(RBPSubRow &sr, const std::vector<unsigned> &cellIds,
                      unsigned beginSite = 0, unsigned endSite = UINT_MAX);

  void permuteCellsInSR(RBPSubRow &sr, unsigned firstCellOffset,
                        const std::vector<unsigned> &newCellOrder);

  void placeTerms();  // equally spaces the terminals around the edge of
  // the layout region.

  const HGraphWDimensions &getHGraph() const {
    abkfatal(_hgWDims != NULL, "HGraph not in memory.");
    return *_hgWDims;
  }
  const HGraphWDimensions *getHGraphPtr() const { return _hgWDims; }

  bool dumpHGraph() {
    if (_ownHGWDimsAuxName && _hgWDimsAuxName.empty()) {
#ifdef _MSC_VER
      std::string tempname = std::string(tmpnam(NULL));
      bool success = _hgWDims->saveAsNodesNetsWts(tempname.c_str());
      if (!success) return false;
      _hgWDimsAuxName = tempname + std::string(".aux");
#else
      char templ[] = "tempHGraphXXXXXX";
      const int fd = mkstemp(templ);
      if (fd == -1) return false;
      close(fd);
      remove(templ);
      bool success = _hgWDims->saveAsNodesNetsWts(templ);
      if (!success) return false;
      _hgWDimsAuxName = std::string(templ) + std::string(".aux");
#endif
    }
    return true;
  }

  bool destroyHGraph() {
    if (_hgWDimsAuxName.empty()) {
      return false;
    } else {
      delete _hgWDims;
      _hgWDims = NULL;
      return true;
    }
  }

  bool rebuildHGraph() {
    if (_hgWDims == NULL) {
      Parameters newParams(_params);
      newParams.makeAllSrcSnk = false;
      _hgWDims =
          new HGraphWDimensions(_hgWDimsAuxName.c_str(), NULL, NULL, newParams);
      markMacros();
      _hgWDims->expandStdCellWidthToFitSiteWidth(_coreRows[0].getSpacing());
      _hgWDims->setNodeWeightsToAreas(_coreRows[0].getSpacing(),
                                      _coreRows[0].getHeight());
      return true;
    } else
      return false;
  }

  void clearNamesFromHGraph() {
    if (_hgWDims != NULL) {
      _hgWDims->discardNames();
    }
  }

  double evalHPWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                  bool usePinLocs = false) const;
  double evalHPWL(bool usePinLocs = false) const {
    return evalHPWL(getPlacement(), getHGraph(), usePinLocs);
  }
  double evalTraditionalHPWL(const PlacementWOrient &pl,
                             const HGraphWDimensions &hg,
                             bool usePinLocs) const;
  double evalTraditionalWeightedHPWL(const PlacementWOrient &pl,
                                     const HGraphWDimensions &hg,
                                     bool usePinLocs) const;
  double evalTraditionalWeightedHPWL(bool usePinLocs = false) const {
    return evalTraditionalWeightedHPWL(getPlacement(), getHGraph(), usePinLocs);
  }
  double evalFSFPCenterOfGravityHPWL(const PlacementWOrient &pl,
                                     const HGraphWDimensions &hg,
                                     bool usePinLocs) const;
  double evalFSFPPinToPinHPWL(const PlacementWOrient &pl,
                              const HGraphWDimensions &hg,
                              bool usePinLocs) const;

  double evalMSTWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                   bool usePinLocs = false) const;
  double evalMSTWL(bool usePinLocs = false) const {
    return evalMSTWL(getPlacement(), getHGraph(), usePinLocs);
  }
  double evalWeightedWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                        bool usePinLocs = false) const;
  double evalWeightedWL(bool usePinLocs = false) const {
    return evalWeightedWL(getPlacement(), getHGraph(), usePinLocs);
  }

  double evalNetsHPWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                      const std::vector<unsigned> &netsToEval,
                      bool usePinLocs = false) const;
  double evalNetsHPWL(const std::vector<unsigned> &netsToEval,
                      bool usePinLocs = false) const {
    return evalNetsHPWL(getPlacement(), getHGraph(), netsToEval, usePinLocs);
  }

  // by sadya to evaluate HPWL in X and Y contributions
  std::pair<double, double> evalXYHPWL(const PlacementWOrient &pl,
                                       const HGraphWDimensions &hg,
                                       bool usePinLocs = false) const;
  std::pair<double, double> evalXYHPWL(bool usePinLocs = false) const {
    return evalXYHPWL(getPlacement(), getHGraph(), usePinLocs);
  }

  bool isConsistent();
  friend std::ostream &operator<<(std::ostream &, const RBPlacement &);

  // save a gnu-plot usable file drawing the sub-rows
  void printRows(const char *filename);

  // by sadya get total site area of the layout
  double getSitesArea();
  void saveAsSCL(const char *baseFileName) const;
  void saveAsLogic(const char *baseFileName) const;
  // save SCL with sites removed from congested regions
  void saveAsSCLCong(const char *baseFileName);
  void savePlNoDummy(const char *baseFileName);
  void removeSitesFromCongestedRgn(double whitespaceRatio = 40);

  void saveAsNodesNetsWts(const char *baseFileName, bool fixMacros = false,
                          bool saveNodeWts = true,
                          bool saveNetWts = true) const;

  // save with fractWS % of whitespace as dummy cells
  void saveAsNodesNetsWtsWDummy(const char *baseFileName, double fractWS = 0,
                                bool fixMacros = false);

  // save with fract % of cells tethered
  void saveAsNodesNetsWtsWTether(const char *baseFileName, double fract = 0,
                                 double rgnSizePercent = 0.005,
                                 bool takeNodeConstrFrmFile = false,
                                 bool takeNetConstrFrmFile = false,
                                 double tetherNewAR = -1,
                                 double tetherNewWS = -1);

  // save with new nodes added to congested nodes
  void saveAsNodesNetsWtsWCongInfo(const char *baseFileName,
                                   double whitespaceRatio = 20);

  // added by sadya to shred tall cells
  double getMaxHeightCoreRow(void) const;
  double getMinWidthNode(void) const;
  double getMaxWidthNode(void) const;
  double getAvgWidthNode(void) const;

  void saveAsNodesNetsWtsShred(const char *baseFileName) const;
  void saveAsNodesNetsWtsUnShred(const char *baseFileName) const;

  void saveAsNodesNetsWtsShredHW(const char *baseFileName, double maxHeight = 0,
                                 double maxWidth = 0,
                                 unsigned singleNetWt = 0) const;
  // shred in both H & V direction
  void savePlacementUnShredHW(const char *plFileName) const;
  clusterNodeOrient getOrientSubNode(unsigned index, unsigned indexRight,
                                     unsigned indexTop) const;
  const char *toChar(clusterNodeOrient orient) const;
  int getQuad(Point &pt) const;

  // save in floorplan format
  void saveAsNodesNetsPlFloorplan(const char *baseFileName) const;
  void saveMacrosAsNodesNetsPlFloorplan(const char *baseFileName) const;

  // added by sadya to remove sites occupied by macros
  void removeSitesBelowFixed(void);
  void removeSitesBelowMacros(void);

  // added by DAP to add a single obstacle, fairly slow i think
  // WARNING: I dont even want to know what happens if you try
  // to place an obstacle on top of a cell, do so at your own risk
  void addObstacle(const BBox &obstacle);

  // added by sadya to get number of macros in a design
  unsigned getNumMacros(void) const;
  unsigned getNumFixedCore(void) const;
  unsigned getNumFixedMacros(void) const;
  double getMovableArea(void) const;

  void saveAsNetDAre(const char *baseFileName) const;
  void savePlacement(const char *plFileName, bool labelFixed = false) const;
  //  { _placement.save(plFileName, _hgWDims->getNumTerminals());}
  void fillInDEF(const char *surrogate, const char *DEFName) const;

  void saveAsSpatialConstraints(const char *constFileName) const;

  // added by sadya to output design in cplace format
  void saveAsCplace(const char *baseFileName);
  // by sadya. output design in plato format
  void saveAsPlato(const char *baseFileName);

  // added by sadya. plots the nodes in gnuplot readable format
  void saveAsPlot(const Plotters::RBPlacePlotter::Parameters &plotterParams,
                  const char *baseFileName) const;

  // plot nodes with prefix dummy in their name
  void plotDummy(const Plotters::RBPlacePlotter::Parameters &plotterParams,
                 const char *baseFileName) const;

  // by sadya. plot only the requested nodes and nets connected to them
  void plot(const Plotters::RBPlacePlotter::Parameters &plotterParams,
            const std::vector<unsigned> &nodes,
            const char *baseFileName) const;

  //  void printGnuplotSiteMap(std::ostream& gpFile) const;
  //  void printGnuplotRows(std::ostream& gpFile) const;

  // plot nodes in the group constraints and their nets
  void plotGrpConstr(const Plotters::RBPlacePlotter::Parameters &plotterParams,
                     const char *baseFileName) const;

  // by sadya to plot the local whitespace histogram based on bins
  void plotWSHist(const char *baseFileName);

  // by sadya to pack all cells to the left corner of the layout
  void shiftCellsLeft();

  // by sadya to allocate whitespace equally in all rows
  void spaceCellsEquallyInRows(void);

  // by sadya to allocate whitespace in all rows according to congestion info
  void spaceCellsWCongInfoInRows();

  // by sadya to allocate whitespace in all rows according to congestion info
  void spaceCellsWCongInfoInRows1();

  // added by sadya to remove overlaps in a given placement by using the
  // whitespace in a subrow
  void remOverlapsSR(void);

  // added by sadya to fix cells in sites if possible
  void snapCellsInSites(void);

  // added by sadya to remove overlaps in a given placement by using the
  // whitespace in a subrow and adjoining subrows
  void remOverlapsVert(void);
  void shuffleSR(itRBPSubRow its);
  unsigned getNumOverlaps();

  // added by sadya to remove overlaps using BFS search strategy
  void remOverlaps(void);

  // added by royj for a different overlap strategy
  void swapCells(RBPSubRow *sr1, unsigned sr1Cell, const Point &sr1CellPos,
                 RBPSubRow *sr2, unsigned sr2Cell, const Point &sr2CellPos);
  double CalcSwapCost(const unsigned &cell1, const Point &cell1NewPos,
                      const bool skipLargeNets) const;
  double CalcSwapCost(const unsigned &cell1, const unsigned &cell2,
                      const Point &cell1NewPos, const Point &cell2NewPos,
                      const bool skipLargeNets) const;
  double CalcSwapCost(const unsigned &cell1, const unsigned &cell2,
                      const unsigned &cell3, const Point &cell1NewPos,
                      const Point &cell2NewPos, const Point &cell3NewPos,
                      const bool skipLargeNets) const;
  double calcSwapCost(const unsigned &cell1, const unsigned cell2,
                      const unsigned cell3, const Point &cell1NewPos,
                      const Point &cell2NewPos, const Point &cell3NewPos,
                      const bool skipLargeNets = false) const;
  double calcSwapCost(const itRBPSubRow &sr1, unsigned sr1Cell,
                      Point &sr1CellPos, const itRBPSubRow &sr2,
                      unsigned sr2Cell, Point &sr2CellPos) const;
  unsigned findBestSwap(
      const std::vector<unsigned> &ordering,
      const std::vector<std::pair<itRBPSubRow, double> > &haveSpace,
      const itRBPSubRow &sourceRow, unsigned sourceCell, Point &newSourcePos,
      unsigned &swapCell, Point &newSwapPos, double &bestImprovement,
      double &bestCost, const double &distCutOff);
  void royjRemOverlaps(MaxMem *maxMem);

  // added by royj for macro overlap removal
  bool findBestMacroMove(
      const PlacementWOrient &pl,
      const std::vector<std::vector<std::vector<unsigned> > > &grid,
      const double &gridXSize, const double &gridYSize,
      const unsigned &numGridRows, const unsigned &numGridCol,
      const unsigned &macro, const std::vector<macroData> &macros,
      Point &trans, double &improve, bool snapMacrosX, bool snapMacrosY);
  void moveMacro(const PlacementWOrient &pl,
                 std::vector<std::vector<std::vector<unsigned> > > &grid,
                 const double &gridXSize, const double &gridYSize,
                 const unsigned &numGridRows, const unsigned &numGridCol,
                 const unsigned &macro_to_move, const Point newPos,
                 const std::vector<macroData> &macros);
  bool doIOverlapHere(
      const PlacementWOrient &pl,
      const std::vector<std::vector<std::vector<unsigned> > > &grid,
      const double &gridXSize, const double &gridYSize,
      const unsigned &numGridRows, const unsigned &numGridCol,
      const unsigned &macro, const Point newPos,
      const std::vector<macroData> &macros);
  void removeMacroOverlaps(MaxMem *maxMem, bool snapMacrosX = false,
                           bool snapMacrosY = false, bool putCellsInRows = true,
                           bool print = true, unsigned badMoveCap = 1);
  void removeMacroOverlaps(MaxMem *maxMem, PlacementWOrient &pl,
                           const std::vector<unsigned> &cellIds,
                           bool snapMacrosX, bool snapMacrosY,
                           bool putCellsInRows, bool print,
                           unsigned badMoveCap);

  unsigned snapMacrosX(PlacementWOrient &pl,
                       const std::vector<unsigned> &cellIds,
                       const BBox &coreArea, double rowSpacing,
                       double siteSpacing);

  unsigned snapMacrosY(PlacementWOrient &pl,
                       const std::vector<unsigned> &cellIds,
                       const BBox &coreArea, double rowSpacing,
                       double siteSpacing);

#ifdef USEFLOORIST
  void doFloorist(PlacementWOrient &pl, const std::vector<unsigned> &cellIds);
#endif

  // added by royj for greedy detailed placement
  void findBestGreedySwap(greedySwapData &cells);
  void fillInAffectedCells(unsigned cell, std::vector<unsigned> &affectedCells,
                           const std::vector<unsigned> &allcells);
  double getAvailSpace(unsigned cell) const;
  bool movableCheck(const std::vector<unsigned> &possible_move) const;
  void greedySwapping(MaxMem *maxMem);

  std::pair<double, double> getSurroundingSpaceSpan(unsigned cell) const;
  std::pair<double, double> calculateGreedyMoveImprove(
      BaseGeneric2DResourceMap *theMap, unsigned cell, const Point &newPos);
  void findBestGreedyMove(BaseGeneric2DResourceMap *theMap, double cutoffRadius,
                          const BBox &layout, unsigned cell, Point &newPos,
                          std::pair<double, double> &improve);
  void findBestGreedyMoveInBox(BaseGeneric2DResourceMap *theMap, unsigned cell,
                               Point &newPos,
                               std::pair<double, double> &improve,
                               const BBox &box, const Point &weberCenter);
  void greedyMovement(greedyMoveMode optimizeMode,
                      BaseGeneric2DResourceMap *theMap);

  // added by sadya. checks for out of core cells after placement
  bool checkOutOfCoreCells(void);

  // by sadya to find closest whitespace to a given point. for annealing
  bool findClosestWS(Point &loc, Point &WSLoc, double width);

  // added by sadya to improve row ironing
  void permuteCellsAndWSInSR(RBPSubRow &sr, unsigned firstCellOffset,
                             const std::vector<unsigned> &newCellOrder);

  // added by sadya to update placement of 2 dim ironed cells
  void updateIronedCells(const std::vector<unsigned> &movables,
                         const Placement &soln, RBPSubRow &subrow1,
                         RBPSubRow &subrow2);

  // added by sadya to update placement of 2 dim ironed cells with clustering
  // L->R
  void updateIronedCellsLR(
      const std::vector<std::vector<unsigned> > &movables,
      const Placement &soln, RBPSubRow &subrow1, RBPSubRow &subrow2);

  // added by sadya to update placement of 2 dim ironed cells with clustering
  // R->L
  void updateIronedCellsRL(
      const std::vector<std::vector<unsigned> > &movables,
      const Placement &soln, RBPSubRow &subrow1, RBPSubRow &subrow2);

  // added by sadya to calculate the HPWL of a group of cells
  //  will be used in overlap removal and rowIroning
  double calcInstHPWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                      std::vector<unsigned> cellIds);

  // added by sadya to save LEF/DEF pair from current RB Placement
  // Generates own macro cells and routing information
  void saveLEFDEF(const char *baseFileName, bool markMacrosAsBlocks = false);

  // added by sadya to get the BBox of the layout core region + terminals
  void initBBox(void);
  BBox getBBox(bool addTerminals = true) const;

  // added by sadya to greedily assign pins to cells. Will change hgraph
  void assignPinsToCells(bool randomizeMacroDims = false);

  // only updates placement of cell. do not use if need to update row info as
  // well. instead use setLocation
  void updatePlacement(unsigned id, const Point &pt);
  void updatePlacementWOri(const PlacementWOrient &pl);

  // by sadya. gets ids of cells in given BBox
  void getCellsInBBox(BBox &bbox, std::vector<unsigned> &cellIds);

  std::vector<std::pair<BBox, std::vector<unsigned> > >
      groupRegionConstr;  // in SYNP
  // by sadya for groups in region constraints. populated from constraints
  std::vector<unsigned> cellToGrpMapping;
  void initGroupRegionConstr(void);
  std::vector<std::pair<BBox, double> >
      regionUtilization;  // by sadya for region
  // utilization constraints. added while at SYNP
  void initRegionUtilization(void);

  void reshapePlacement(double tetherNewAR, double tetherNewWS,
                        const char *newPlFileName);

  void splitCoreRowsByCell(const PlacementWOrient &pl, unsigned cellId);
  void _splitRowByObstacle(RBPCoreRow &cr, const BBox &obstacle);
  void addToCoreRows(const PlacementWOrient &pl, unsigned cellId);
  void removeFromCoreRows(const PlacementWOrient &pl, unsigned cellId);

  // adds rows from backup into _coreRows
  // needs backUpCoreRows ot be initialized earlier
  void addRowsByObstacleFrmBackup(const BBox &obstacle);

  // by sadya assumes the isMacro field in HGraphWDims is updated
  void markMacrosAsFixed(void);
  // void unmarkMacrosAsFixed(void);
  void markMacros(void);

  // by hhchan, assumes the isMacro field in HGraphWDims is updated
  void markMacrosAsSoft(double minAR = HGraphWDimensions::DEFAULT_MIN_AR,
                        double maxAR = HGraphWDimensions::DEFAULT_MAX_AR);
  unsigned getNumSoft() const;

  void backUpCoreRows(void);
  void reinstateCoreRows(void);
  void reinstateCoreRowsAndRepopulate(void);
  void clearStdCellsFromCoreRows(void);
  const RBPCoreRow &getBackupCoreRow(unsigned id) const;
  unsigned getNumBackupCoreRows() const { return _backUpCoreRows.size(); }

#ifdef USEFLOW
  void doFlow(ISPD06DensityMap *map, bool Xflow, bool Yflow, MaxMem *maxMem);
  void doFlowHelper(ISPD06DensityMap *map, bool Xflow, double Xdist, bool Yflow,
                    double Ydist, PlacementWOrient &pl,
                    std::vector<unsigned> *allowedMovables, MaxMem *maxMem);
  void doXFlowInBBox(ISPD06DensityMap *map,
                     std::vector<std::vector<unsigned> > &coreCols,
                     const BBox &theBox,
                     std::vector<unsigned> *allowedMovables, double Xdist,
                     MaxMem *maxMem);
  void doYFlowInBBox(ISPD06DensityMap *map,
                     std::vector<std::vector<unsigned> > &coreCols,
                     const BBox &theBox,
                     std::vector<unsigned> *allowedMovables, double Ydist,
                     MaxMem *maxMem);
  void doUnconstrainedXFlow(const PlacementWOrient &origPl,
                            const HGraphWDimensions &hg,
                            const std::vector<unsigned> &cellIds,
                            const std::vector<BBox> &cellBBoxes,
                            PlacementWOrient &changedPl, MaxMem *maxMem) const;
  void doUnconstrainedYFlow(const PlacementWOrient &origPl,
                            const HGraphWDimensions &hg,
                            const std::vector<unsigned> &cellIds,
                            const std::vector<BBox> &cellBBoxes,
                            PlacementWOrient &changedPl, MaxMem *maxMem) const;
#endif

  std::string getLogicValue(unsigned nodeIdx) const;

  double getRowHeight() const { return _coreRows[0].getHeight(); }
  double getRowSpacing() const { return _coreRows[0].getSpacing(); }

  void doCongestionMoves(void);
};

typedef RBPlacement::Parameters RBPlaceParams;

RBPlacement &operator<<(RBPlacement &db, const std::vector<Orient> &ori);

#endif
