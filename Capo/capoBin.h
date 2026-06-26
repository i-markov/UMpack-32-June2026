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

// created by Andy Caldwell on 09/14/97  caldwell@cs.ucla.edu

// CHANGES:
// mm  12/10/97 added const netsBegin(), netsEnd() mem fns
// ilm 980304   added const getNets(), netsEnd() method
// ilm 980316   added unsigned getNumSubrootChildren()
// ilm 980324   split clSubTreeH.h from this file (sim. for the .cxx)

// aec removed subRoot.  Bin is now just a vector of leaf clusters.
// aec 981118   re-worked CapoBins.  All Cluster classes are removed.
// aec 990203   moved to RBPlace Rows
// aec 990223   added getWhitespace and getRelativeWhitespace
// ilm 990319   bins now maintain a reflexive 4-directed neighborhood relation

#ifndef __CAPOBIN_H__
#define __CAPOBIN_H__

#include <algorithm>
#include <iostream>
#include <utility>

#include "ABKCommon/abkassert.h"
#include <string>
#include <vector>
#include "HGraphWDims/hgWDims.h"
#include "Partitioning/partitionData.h"
#include "RBPlace/RBRows.h"

// see Capo/DOCS/CapoBin.doc for a description of class CapoBin

typedef std::pair<unsigned, unsigned> CROffset;
class junkForCompilerCROffset : public CROffset {};

typedef std::vector<const RBPCoreRow*>::const_iterator itCBCoreRow;
typedef std::vector<CROffset>::const_iterator itCBRowOffset;

// forward declaration
class CapoPlacer;

class CapoBin {
  friend class CapoPlacer;
  friend class BaseBinSplitter;
  friend class CapoBinMatchRec;
  friend class LookAheadSplitter;

  std::vector<const RBPCoreRow*> _coreRows;  // core rows (at least partly)
                                              // contained in this bin
                                              // the first coreRow is the lowest
                                              //(smallest yCoord)

  std::vector<CROffset> _startOffsets;
  // for each CoreRow, the # of
  // subRows which are left of the
  // bin (same as index of 1st
  // subRow intersecting the bin),
  // and # of sites (in left-most
  // intersecting subRow) left of the
  // bin.
  // A start offset <i,j> means
  // that the jth cell of the ith
  // sub row of that CoreRow is the
  // first site actually in this bin
  std::vector<CROffset> _endOffsets;
  // for each CoreRow, the
  // index of the last subRow
  // intersecting the bin, and
  // of the last site in the bin.
  // An end offset of <i,j> means
  // that the jth cell of the ith
  // subRow is the last site actually
  // in the bin.
  // IMPORTANT: note that unlike STL begin/end
  // iterators where 'end' iterators are outside
  // the collection, the subrows/sites referenced
  // by 'endOffsets' ARE actually in the bin.
  // iteration is difficult otherwise...
  //(this means that start/endOffsets are
  // completely symetric.

  unsigned _numSites;

  double _capacity;
  // sum of the sites on the rows, within
  // the interval this bin covers

  BBox _bbox;
  // defines the layout boundries of
  // this bin.  The top and bottom of
  // the bbox corespond to the top of the
  // top-most row, and the bottom of the
  // bottom-most row, respecively

  std::vector<unsigned> _cellIds;
  // the Ids of the cells (movable,
  // and otherwise) contained in this
  // bin.  Capo has a single netlist-
  // level HGraph, to which these Ids
  // correspond

  std::vector<unsigned> _oracleCellIds;
  bool _ecoAllowed;

  // <aaronnn> any non-movable/placed cells in this bin get tracked here.
  // vector is pruned upon bin split and grown upon bin merges:
  std::vector<unsigned> _obstacleCellIds;

  double _totalCellArea;  // total area of all cells
                          // in this bin.
  double _totalMacroArea;  // total area of macros in this bin
  double _maxRecCellArea;  // whats the max amount
                           // of cell area that should be in this bin
  double _avgRowSpacing;   // vertical spacing between rows
  double _avgCellWidth;
  double _avgCellHeight;

  double _largestCellArea;  // the area of the largest
                            // cell in this bin.

  double _medianCellArea;  // area of median cell area

  std::string _name;  // the name of this bin.  This is
                       // used when saving problems

  unsigned _index;  // CapoBins are uniquely identified

  CapoPlacer& _capo;  // this is the capo placer object
                      // associated with this bin, this
                      // bin does not own the _capo object.

  void resetCellIds(const std::vector<unsigned>& newCellIds);
  // rips existing cells out of the bin
  // and repopulates the bin with newCellsIds

  bool _canSplitBin;
  bool _isNewBin;

  // by DAP to store COG of sites,
  // for using Center-of-Gravity as default location of cells in the bin,
  // and stop using the geometric center
  Point _centerOfGravity;

  //    std::vector<unsigned>* _savedSoln;

 protected:
  bool _needFP;        // this is set when a floorplanning attempt fails
  bool _floorplanned;  // this is set when a floorplanning attempt is
                       // successfull
  bool _ecoFP;
#ifdef USEPATOMA
  bool _patomaEndCase;  // "true" when the current bin is a Patoma end case
  bool _patomaMacrosCommitted;  // "true" if the macros are fixed
  bool _patomaPacked;  // "true" if this bin is checked to be able to fit
#endif
  std::vector<CapoBin*> _neighborsAbove;
  std::vector<CapoBin*> _neighborsBelow;
  std::vector<CapoBin*> _leftNeighbors;
  std::vector<CapoBin*> _rightNeighbors;

  CapoBin* _sibling;

  // by DAP to store COG of sites,
  // for using Center-of-Gravity as default location of cells in the bin,
  // and stop using the geometric center
  void computeCenterOfGravity();

  bool isSibling(CapoBin* candidate) const { return _sibling == candidate; }

  void _addNeighborAbove(CapoBin* neighbor) {
    _neighborsAbove.push_back(neighbor);
  }
  void _addNeighborBelow(CapoBin* neighbor) {
    _neighborsBelow.push_back(neighbor);
  }
  void _addLeftNeighbor(CapoBin* neighbor) {
    _leftNeighbors.push_back(neighbor);
  }
  void _addRightNeighbor(CapoBin* neighbor) {
    _rightNeighbors.push_back(neighbor);
  }

  void __removeFrom(std::vector<CapoBin*>& v, CapoBin* ptr) {
    unsigned sizeBefore = v.size();
    if (sizeBefore == 0) return;

    std::vector<CapoBin*>::iterator newEnd = remove(v.begin(), v.end(), ptr);
    v.erase(newEnd, v.end());
    if (sizeBefore <= v.size()) {
      std::cerr << std::endl
                << "  Bin [" << getIndex() << "] does not know about"
                << " its neighboring bin [" << ptr->getIndex() << "] ";
      abkfatal(0, "");
    }
  }

  void _removeNeighborAbove(CapoBin* neighbor) {
    __removeFrom(_neighborsAbove, neighbor);
  }
  void _removeNeighborBelow(CapoBin* neighbor) {
    __removeFrom(_neighborsBelow, neighbor);
  }
  void _removeLeftNeighbor(CapoBin* neighbor) {
    __removeFrom(_leftNeighbors, neighbor);
  }
  void _removeRightNeighbor(CapoBin* neighbor) {
    __removeFrom(_rightNeighbors, neighbor);
  }

  void setIndexToNextIndex();

  // populates coreRows and start/endOffsets with the intersection
  // of candidateRows and boundary
  void compStartAndEndOffsets(
      const BBox& boundary,
      const std::vector<const RBPCoreRow*>& candidateRows,
      std::vector<const RBPCoreRow*>& coreRows,
      std::vector<CROffset>& startOffsets, std::vector<CROffset>& endOffsets);

 public:
  CapoBin(const std::vector<unsigned>& cellIds, itCBCoreRow rowsBegin,
          itCBCoreRow rowsEnd, const std::vector<CROffset>& startOffsets,
          const std::vector<CROffset>& endOffsets,
          const std::vector<CapoBin*>& neighbAbove,
          const std::vector<CapoBin*>& neighbBelow,
          const std::vector<CapoBin*>& leftNeighb,
          const std::vector<CapoBin*>& rightNeighb, CapoPlacer& capo,
          const std::string& name = "");

  // this ctor calculates the included rows & their start/end offsets.
  // it includs all sites within the intersection of candidateRows
  // and boundary
  CapoBin(const std::vector<unsigned>& cellIds,
          const std::vector<const RBPCoreRow*>& candidateRows,
          const BBox& boundary, const std::vector<CapoBin*>& neighbAbove,
          const std::vector<CapoBin*>& neighbBelow,
          const std::vector<CapoBin*>& leftNeighb,
          const std::vector<CapoBin*>& rightNeighb, CapoPlacer& capo,
          const std::string& name = "");

  CapoBin(std::vector<CapoBin*>& src);

  // by sadya
  CapoBin(CapoBin* srcBin);
  CapoBin(std::vector<CapoBin*>& src, int);

  ~CapoBin();

  unsigned getNumCells() const { return _cellIds.size(); }

  std::vector<unsigned>::const_iterator cellIdsBegin() const {
    return _cellIds.begin();
  }
  std::vector<unsigned>::const_iterator cellIdsEnd() const {
    return _cellIds.end();
  }
  const std::vector<unsigned>& getCellIds() const { return _cellIds; }
  const std::vector<unsigned>& getOracleCellIds() const {
    return _oracleCellIds;
  }
  void setOracleCellIds(const std::vector<unsigned>& newIds) {
    _oracleCellIds = newIds;
  }
  void clearOracleCellIds() { _oracleCellIds.clear(); }
  bool isEcoAllowed() const { return _ecoAllowed; }
  void setEcoAllowed(bool allowed) { _ecoAllowed = allowed; }

  // <aaronnn>: bin's obstacle fns
  void setObstacleCellIds(const std::vector<unsigned>& newIds) {
    _obstacleCellIds = newIds;
  }
  const std::vector<unsigned>& getObstacleCellIds() const {
    return _obstacleCellIds;
  }
  unsigned getNumObstacles() const { return _obstacleCellIds.size(); }

  const std::vector<const RBPCoreRow*>& getRows() const { return _coreRows; }
  unsigned getNumRows() const { return _coreRows.size(); }
  itCBCoreRow rowsBegin() const { return _coreRows.begin(); }
  itCBCoreRow rowsEnd() const { return _coreRows.end(); }
  const RBPCoreRow& getRow(unsigned r) const { return *_coreRows[r]; }

  std::pair<unsigned, unsigned> getStartOffset(unsigned r) const {
    return _startOffsets[r];
  }
  std::pair<unsigned, unsigned> getEndOffset(unsigned r) const {
    return _endOffsets[r];
  }

  double getCapacity() const { return _capacity; }
  double getNumSites() const { return _numSites; }
  double getWhitespace() const { return _capacity - _totalCellArea; }
  double getRelativeWhitespace() const { return getWhitespace() / _capacity; }
  double getOverfill() const { return std::max(0.0, -getWhitespace()); }
  double getRelativeOverfill() const { return getOverfill() / _capacity; }

  double getContainedAreaInCoreRow(unsigned r) const;
  unsigned getContainedSitesInCoreRow(unsigned r) const;
  double getContainedLeftEdge(unsigned r) const;
  double getContainedRightEdge(unsigned r) const;
  double getContainedAreaInBBox(const BBox& rect) const;  // by sadya in SYNP

  const BBox& getBBox() const { return _bbox; }
  double getHeight() const { return _bbox.yMax - _bbox.yMin; }
  double getWidth() const { return _bbox.xMax - _bbox.xMin; }
  Point getCenter() const;

  double getTotalCellArea() const { return _totalCellArea; }
  // jflu: setTotalCellArea(double) is used for cell bloating
  void setTotalCellArea(double cellArea) { _totalCellArea = cellArea; }
  double getTotalMacroArea() const { return _totalMacroArea; }
  double getMaxRecommendedCellArea() const { return _maxRecCellArea; }
  void setMaxRecCellArea(double area) { _maxRecCellArea = area; }
  double getAvgCellWidth() const { return _avgCellWidth; }
  double getAvgCellHeight() const { return _avgCellHeight; }
  double getAvgRowSpacing() const { return _avgRowSpacing; }
  double getFactor() const {
    if (_avgCellWidth + _avgCellHeight < 0.1) return 1.0;
    return _avgRowSpacing / (_avgCellWidth + (1.0 / 3) * _avgCellHeight);
  }
  double getLargestCellArea() const { return _largestCellArea; }
  double getWidestCellWidth() const;

  double getMedianCellArea() const { return _medianCellArea; }

  unsigned getIndex() const { return _index; }

  unsigned getNumAdjacentCells() const;
  unsigned getNumAdjacentCellsWithDuplicates() const;

  // returns the x-line s.t. as close as possible to p0Cap
  // site area is left of the x-line.
  double findXSplit(double c0, double c1, double roundPct, double p0WSWeight,
                    double p1WSWeight, std::vector<double>& partCaps) const;
  double findXSplitUsingSkyline(double c0, double c1, double p0WSWeight,
                                double p1WSWeight,
                                std::vector<double>& partCaps) const;

  // shiftXSplitAndSnap(...) is a shell function for calling
  // shiftXSplit(...) followed by a call to skylineSnappingForOracle(...)
  // if useCutOracle is not enabled, then shiftXSplitAndSnap will accept
  // the result of shifting and not snap
  double shiftXSplitAndSnap(double xsplit, double localWS, double c0, double c1,
                            double roundPct,
                            std::vector<double>& partCaps) const;
  double shiftXSplit(double xsplit, double minLocalWS, double c0, double c1,
                     std::vector<double>& partCaps) const;
  double skylineSnapingForOracle(double xsplit, double localWS, double c0,
                                 double c1, double roundPct,
                                 std::vector<double>& partCaps) const;

  // in SYNP
  // finds split direction and point. tries to align the cut line to
  // boundary of a region. if no regions in the given bin or one region
  // occupies entire bin, then returns false. proceed with normal flow
  // returns desired splitPt
  bool findSplitPtDir(bool desiredDir,
                      std::vector<std::pair<BBox, std::vector<unsigned> > >&
                          groupRegionConstr,
                      bool& splitDir, double& splitPt);

  void computePartAreas(double xSplit, std::vector<double>& partCaps) const;

  // by sadya in SYNP to align cutlines to big obstacle edges
  bool shouldAlignCL2Obs(double& splitPt, bool& splitHoriz);

  // checks if macros are in the bin and satisfy the constraints of legality
  bool fpCondMet(bool desiredDir, double& splitPt);
  // same as fpCondMet, but only includes those conditions that imply
  // no legal solution is possible.  Currently this includes cond. I and II
  // for the current bin as is, as opposed to projected child bins).
  bool fpCondMet_Strict(void);
#ifdef USEPATOMA
  // hhchan: Patoma data-structures, setters and getters
  bool patomaCondMet() const;
  bool isPatomaEndCase() const {
    abkwarn((!_patomaMacrosCommitted) || _patomaEndCase,
            "PATOMA: macros are fixed but the bin isn't an end case");
    return _patomaEndCase;
  }
  bool isPatomaMacrosCommitted() const { return _patomaMacrosCommitted; }
  bool isPatomaPacked() const { return _patomaPacked; }

  void setPatomaEndCase(bool isEndCase) { _patomaEndCase = isEndCase; }
  void setPatomaMacrosCommitted(bool areCommitted) {
    _patomaMacrosCommitted = areCommitted;
  }
  void setPatomaPacked(bool isPacked) { _patomaPacked = isPacked; }
#endif
  double largestMacroArea();

  void linkNeighbors();
  void unLinkNeighbors();

  const std::string& getName() const { return _name; }

  bool canSplitBin() const { return _canSplitBin; }
  bool isNewBin() const { return _isNewBin; }

  //    std::vector<unsigned>* getSavedSoln() const
  //                      { return _savedSoln; }
  //    void setSavedSoln(std::vector<unsigned> *newSoln)
  //         { _savedSoln = newSoln; }

  friend std::ostream& operator<<(std::ostream& os, const CapoBin& ct);

  void decrementIdx();

  // by sadya. removes macros from bin. updates relevant info about bin.
  void removeMacros();

  // <aaronnn> removes only macros that were floorplanned
  void removeFloorplannedMacros();

  // resets row structure of each bin,
  // assumes the rows below macros have been split allready if that is the
  // reason for calling this routine
  void resetRows();

  class hCutLineInfo {
   public:
    unsigned splitRow;
    double cutCost, p0Area, p1Area;
  };

  class vCutLineInfo {
   public:
    double xSplit, cutCost, p0Area, p1Area;
  };

  hCutLineInfo scanCutLinesH(const SubHGraph& probSubHGraph, double p0Min,
                             double p0Max, double p1Min, double p1Max) const;
  vCutLineInfo scanCutLinesV(const SubHGraph& probSubHGraph, double p0Min,
                             double p0Max, double p1Min, double p1Max) const;
};

std::ostream& operator<<(std::ostream& os, const CapoBin& ct);

class CompareCellIdsByArea {
  const HGraphFixed& _hgraph;

 public:
  CompareCellIdsByArea(const HGraphFixed& hgraph) : _hgraph(hgraph) {}

  bool operator()(unsigned cellId1, unsigned cellId2) {
    if (_hgraph.getWeight(cellId1) < _hgraph.getWeight(cellId2)) return true;
    return false;
  }
};

struct CompareSrcBinsByX {
  CompareSrcBinsByX() {}
  bool operator()(const CapoBin* b1, const CapoBin* b2) {
    return (b1->getBBox().xMin < b2->getBBox().xMin);
  }
};

struct CompareSrcBinsByY {
  CompareSrcBinsByY() {}
  bool operator()(const CapoBin* b1, const CapoBin* b2) {
    return (b1->getBBox().yMin < b2->getBBox().yMin);
  }
};

struct CompareBinsByIndex {
  CompareBinsByIndex() {}
  bool operator()(const CapoBin* b1, const CapoBin* b2) {
    return (b1->getIndex() < b2->getIndex());
  }
};

#endif
