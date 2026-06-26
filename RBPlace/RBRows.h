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
// 970930 pat  added findSubRow and findCoreRow signatures
// 971001 pat  moved comparison functors from *.cxx
// 971007 pat  added getSpacing to subRow
// 980222 pat  handle additional wierd cases in subrow construction
// 980222  dv  changed the interface of populateRowsWith, added two more
//             unsigned Data Members, added support for isPopulated.
// 980603  dv  removed orientedRBPlacement class
// 980309 pat  more const declarations; added siteRoundup method
// 980427 pat  added RBPCoreRow::findSubRow

//======================================================================
//                           Row-based Placement
//======================================================================

#ifndef _RBROWS_H_
#define _RBROWS_H_

#include <iostream>
#include <utility>
#include <vector>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abktypes.h"
#include <string>
#include <vector>
#include "HGraphWDims/hgWDims.h"
#include "Placement/placeWOri.h"
#include "Placement/symOri.h"
#include "RBCellCoord.h"

//========================== CLASS INTERFACES =============================

class RBPSubRow;
struct RBPSite;
class RBPCoreRow;

typedef std::vector<unsigned>::const_iterator itRBPCellIds;

class RBPSubRow {
  friend class RBPlacement;
  friend class RBPlaceFromDB;
  friend class RBCPlacement;
  friend class RBPCoreRow;

  std::vector<unsigned> _cellIds;  // ids of cells in this subRow (in x-order)
  double _xStart, _xEnd;
  unsigned _numSites;
  bool _needsSorting;
  RBPCoreRow* _coreRow;  // coreRow containing this subRow

  void clearAllCells() { _cellIds.clear(); }
  void removeCell(unsigned id);  // NOTE: uses the location of <id> to find it
  void addCell(unsigned id);     // NOTE: uses the location of <id> to
                              // put the cell in the right order
  void addCell(unsigned id, const RBCellCoord& coord);
  void removeCell(unsigned id, const RBCellCoord& coord);

  void addCellToEnd(unsigned id);  // append the cell id. they will not be
                                   // in the correct order within this row,
                                   // and should be sorted before use.

  void sortCellsByLoc();
  void resetNumSites();

 public:
  RBPSubRow()
      : _xStart(0),
        _xEnd(0),
        _numSites(0),
        _needsSorting(false),
        _coreRow(NULL) {}

  RBPSubRow(double xStart, double xEnd)
      : _xStart(xStart),
        _xEnd(xEnd),
        _numSites(0),
        _needsSorting(false),
        _coreRow(NULL) {}

  RBPSubRow(double xStart, double xEnd, RBPCoreRow& coreRow);

  RBPSubRow(const RBPSubRow& orig)
      : _xStart(orig._xStart),
        _xEnd(orig._xEnd),
        _numSites(orig._numSites),
        _needsSorting(false),
        _coreRow(orig._coreRow) {
    abkassert(orig.getNumCells() == 0, "Can't copy constr. populated subRow");
  }

  double getXStart() const { return _xStart; }
  double getXEnd() const { return _xEnd; }
  double getLength() const { return _xEnd - _xStart; }
  unsigned getNumSites() const { return _numSites; }

  Orient getOrientation() const;
  const RBPSite& getSite() const;
  double getSpacing() const;
  double getSiteWidth() const;
  double getYCoord() const;
  double getHeight() const;
  const RBPCoreRow& getCoreRow() const { return *_coreRow; }

  unsigned getNumCells() const { return _cellIds.size(); }
  unsigned operator[](unsigned i) const;
  itRBPCellIds cellIdsBegin() const { return _cellIds.begin(); }
  itRBPCellIds cellIdsEnd() const { return _cellIds.end(); }

  bool operator<(const RBPSubRow& sr2) const {
    return this->_xEnd <= sr2._xStart;
  }

  friend std::ostream& operator<<(std::ostream&, const RBPSubRow&);
};

struct RBPSite {
  double height;
  double width;
  Symmetry symmetry;

  RBPSite() : height(1), width(1) {}

  RBPSite(double ht, double wth, Symmetry sym)
      : height(ht), width(wth), symmetry(sym) {}
  ~RBPSite() {}
  RBPSite(const RBPSite& orig)
      : height(orig.height), width(orig.width), symmetry(orig.symmetry) {}
  RBPSite& operator=(const RBPSite& orig) {
    height = orig.height;
    width = orig.width;
    symmetry = orig.symmetry;
    return *this;
  }
  bool operator==(const RBPSite& s2) const {
    return (height == s2.height && width == s2.width);
  }
  void setHeight(double ht) { height = ht; }
};

typedef std::vector<RBPSubRow>::const_iterator itRBPSubRow;

// RBPCoreRow is a set of rows of core cells with same y-coord
class RBPCoreRow {
  friend class RBPlacement;
  friend class RBPlaceFromDB;
  friend class RBCPlacement;
  friend class RBPSubRow;

  Orient _orient;  // information about allowed cell orientations
  RBPSite _site;
  double _spacing;  // distance between site beginings
  double _y;
  std::vector<unsigned> _allCells;
  mutable std::vector<std::pair<double, double> > _whitespaceLocs;
  mutable bool _whitespaceLocsPopulated;

  std::vector<RBPSubRow> _subRows;

  // this is a READ-ONLY reference!
  // a lot of implementation details can be simplified if this
  // reference is removed --DAP
  const PlacementWOrient& _placement;

  void save(std::ostream& os) const;
  void buildWhitespaceLocs(const HGraphWDimensions& hg) const;

 public:
  RBPCoreRow(double y, Orient orient, const RBPSite& site,
             const PlacementWOrient& pl, double spacing)
      : _orient(orient),
        _site(site),
        _spacing(spacing),
        _y(y),
        _whitespaceLocsPopulated(false),
        _placement(pl) {
    abkfatal(_spacing >= _site.width,
             "Error in input: site spacing"
             " is smaller than site width");
  }

  RBPCoreRow(double y)
      : _y(y),
        _whitespaceLocsPopulated(false),
        _placement(PlacementWOrient(0)) {}

  RBPCoreRow()
      : _y(0),
        _whitespaceLocsPopulated(false),
        _placement(PlacementWOrient(0)) {}

  RBPCoreRow(const RBPCoreRow& orig)
      : _orient(orig._orient),
        _site(orig._site),
        _spacing(orig._spacing),
        _y(orig._y),
        _allCells(orig._allCells),
        _whitespaceLocs(orig._whitespaceLocs),
        _whitespaceLocsPopulated(orig._whitespaceLocsPopulated),
        _subRows(orig._subRows),
        _placement(orig._placement) {
    abkfatal(_spacing >= _site.width,
             "Error in input: site spacing is smaller"
             " than site width");
    for (std::vector<RBPSubRow>::iterator it = _subRows.begin();
         it != _subRows.end(); ++it)
      it->_coreRow = this;
  }

  RBPCoreRow(const RBPCoreRow& orig, const PlacementWOrient& pl)
      : _orient(orig._orient),
        _site(orig._site),
        _spacing(orig._spacing),
        _y(orig._y),
        _allCells(orig._allCells),
        _whitespaceLocs(orig._whitespaceLocs),
        _whitespaceLocsPopulated(orig._whitespaceLocsPopulated),
        _subRows(orig._subRows),
        _placement(pl) {
    abkfatal(_spacing >= _site.width,
             "Error in input: site spacing is smaller"
             " than site width");
    for (std::vector<RBPSubRow>::iterator it = _subRows.begin();
         it != _subRows.end(); ++it)
      it->_coreRow = this;
  }

  RBPCoreRow& operator=(const RBPCoreRow& orig) {
    _orient = orig._orient;
    _site = orig._site;
    _spacing = orig._spacing;
    _y = orig._y;
    _allCells = orig._allCells;
    _whitespaceLocs = orig._whitespaceLocs;
    _whitespaceLocsPopulated = orig._whitespaceLocsPopulated;
    _subRows = orig._subRows;

    abkfatal(_spacing >= _site.width,
             "Error in input: site spacing is smaller"
             " than site width");
    // this is a huge ugly hack.
    // This is here because of the const PlacementWOrient&
    // this can only be changed at the time of construction,
    // so anytime I use the equal operator, I have to cast
    // away constness.   The elegant solution is to
    // redesign these classes so that they do not
    // depend on a const PlacementWOrient&  --DAP
    const_cast<PlacementWOrient&>(_placement) = orig._placement;

    for (std::vector<RBPSubRow>::iterator it = _subRows.begin();
         it != _subRows.end(); ++it)
      it->_coreRow = this;

    return *this;
  }

  void appendNewSubRow(double startX, double endX) {
    RBPSubRow r(startX, endX, *this);
    r.resetNumSites();
    _subRows.push_back(r);
    abkfatal(_subRows.back()._coreRow == this, "subRow init problem");
  }

  Orient getOrientation() const { return _orient; }
  const RBPSite& getSite() const { return _site; }
  RBPSite& getSite() { return _site; }
  double getSpacing() const { return _spacing; }
  double getSiteWidth() const { return _site.width; }
  double getSiteHeight() const { return _site.height; }
  double getSiteArea() const { return _site.width * _site.height; }
  double getYCoord() const { return _y; }
  double getHeight() const { return _site.height; }
  Symmetry getSymmetry() const { return _site.symmetry; }

  unsigned getNumSubRows() const { return _subRows.size(); }
  itRBPSubRow findSubRow(const double x) const;
  itRBPSubRow subRowsBegin() const { return _subRows.begin(); }
  itRBPSubRow subRowsEnd() const { return _subRows.end(); }

  std::vector<RBPSubRow>::iterator NCsubRowsBegin() {
    return _subRows.begin();
  }
  std::vector<RBPSubRow>::iterator NCsubRowsEnd() { return _subRows.end(); }
  std::vector<RBPSubRow>::iterator eraseSubRow(
      std::vector<RBPSubRow>::iterator sRowIt,
      std::vector<RBPSubRow>::iterator sRowIt1) {
    return _subRows.erase(sRowIt, sRowIt1);
  }

  double getXStart() const { return _subRows[0].getXStart(); }
  double getXEnd() const { return _subRows.back().getXEnd(); }

  double siteRoundUp(double length) const {
    return ceil(length / _spacing) * _spacing;
  }
  double siteRoundDown(double length) const {
    return floor(length / _spacing) * _spacing;
  }

  unsigned partSites(double length) const {
    return static_cast<unsigned>(ceil(length / _spacing));
  }
  unsigned wholeSites(double length) const {
    return static_cast<unsigned>(floor(length / _spacing));
  }

  const RBPSubRow& operator[](unsigned idx) const { return _subRows[idx]; }
  RBPSubRow& operator[](unsigned idx) { return _subRows[idx]; }

  friend std::ostream& operator<<(std::ostream&, const RBPCoreRow&);

  const std::vector<std::pair<double, double> >& getWhitespaceLocs(
      const HGraphWDimensions& hg) const {
    if (!_whitespaceLocsPopulated) {
      buildWhitespaceLocs(hg);
    }
    return _whitespaceLocs;
  }

  void clearWhitespaceLocs() const {
    _whitespaceLocs.clear();
    _whitespaceLocsPopulated = false;
  }

  void addCell(const PlacementWOrient& pl, unsigned cellId) {
    std::pair<std::vector<unsigned>::iterator,
              std::vector<unsigned>::iterator>
        range = std::equal_range(_allCells.begin(), _allCells.end(), cellId,
                                 CompareXWOri(pl));

    //     for(std::vector<unsigned>::iterator i = range.first; i !=
    //     range.second; ++i)
    //     {
    //       abkfatal(*i != cellId, "tried to add a cell and it was already
    //       there");
    //     }
    _allCells.insert(range.first, cellId);
  }

  void removeCell(const PlacementWOrient& pl, unsigned cellId) {
    std::pair<std::vector<unsigned>::iterator,
              std::vector<unsigned>::iterator>
        range = std::equal_range(_allCells.begin(), _allCells.end(), cellId,
                                 CompareXWOri(pl));

    for (std::vector<unsigned>::iterator i = range.first; i != range.second;
         ++i) {
      if (*i == cellId) {
        _allCells.erase(i);
        return;
      }
    }
  }

  void clearAllCells() { _allCells.clear(); }

  void clearMovableStdCells(const bit_vector& isFixed,
                            const HGraphWDimensions& hg) {
    std::vector<unsigned> newAllCells;

    for (unsigned i = 0; i < _allCells.size(); ++i) {
      unsigned cellId = _allCells[i];
      if (isFixed[cellId] || hg.isNodeMacro(cellId)) {
        newAllCells.push_back(cellId);
      }
    }

    _allCells = newAllCells;
  }

  const std::vector<unsigned>& getAllCells(void) const { return _allCells; }
};

inline Orient RBPSubRow::getOrientation() const {
  return _coreRow->getOrientation();
}
inline const RBPSite& RBPSubRow::getSite() const { return _coreRow->getSite(); }
inline double RBPSubRow::getSpacing() const { return _coreRow->getSpacing(); }
inline double RBPSubRow::getSiteWidth() const {
  return _coreRow->getSiteWidth();
}
inline double RBPSubRow::getYCoord() const { return _coreRow->getYCoord(); }
inline double RBPSubRow::getHeight() const { return _coreRow->getHeight(); }

inline unsigned RBPSubRow::operator[](unsigned i) const {
  abkassert(!_needsSorting, "addToEnd called w/o sortCellsByLoc");
  return _cellIds[i];
}

struct compareCoreRowY {
  bool operator()(const RBPCoreRow& cr1, const RBPCoreRow& cr2) const {
    return cr1.getYCoord() < cr2.getYCoord();
  }
};

struct compareSubRowX {
  bool operator()(const RBPSubRow& sr1, const RBPSubRow& sr2) const {
    return sr1.getXEnd() <= sr2.getXStart();
  }
};

struct srowXSlt {
  bool operator()(const RBPSubRow* sr1, const RBPSubRow* sr2) const {
    return sr1->getXStart() < sr2->getXStart();
  }
};

class srowXElt {
 public:
  srowXElt() {}
  bool operator()(const RBPSubRow* sr1, const RBPSubRow* sr2) const {
    return sr1->getXEnd() < sr2->getXEnd();
  }
};

class srowYlt {
 public:
  srowYlt() {}
  bool operator()(const RBPSubRow* sr1, const RBPSubRow* sr2) const {
    return sr1->getYCoord() < sr2->getYCoord();
  }
};

#endif
