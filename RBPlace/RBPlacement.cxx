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

// 97.8.8  Paul Tucker
// 970820 pat  added x-coord sort to subRows in RBPlacement
// 970821 pat  initialize subrow spacing
// 970823 ilm  added orientations to RBPCoreRow and RBPIORow (Maogang's request)
//                  they are initialized now
// 970902 ilm  RBPlacement::RBPlacement(const DB&) does not populate rows
//                  with cell Id numbers anymore. This is done with
//                  RBPlacement::populateRowsWith(const Placement&)
//                  (the code has to be exteneded to put IOPads intoIORows,
//                  and a lot more)
// 970917 pat  fix bug in RBCPlacement ctor for test creation
// 970918 pat  modify populateRowsWith to only put core cells in
//             core subrows, and use binary search.
//             improve operator<< for subRows.
// 970919 pat  more elaborate operator<< for RBCSubRowCol
// 970921 pat  upgraded constructors to split core subrows
//             according to special net wiring;
//             added real RBCPlacement ctor
// 970923 pat  finish RBCPlacement; delete dead code; obey object interfaces
// 970930 pat  erase subRow index vectors at start of populateRowsWith,
//             so legalizer can move cells around
// 970930 pat  add findCoreRow and findSubRow search functions
// 971001 pat  moved comparison functors into .h file
// 971129 pat  comment out old diagnostic message about # of
//             cells not placed in core rows
// 980211 pat  add handling for missing cases to construction of
//             RBPSubRows based on special nets
// 982202 pat  corrected error of creating subRows
// 982202  dv  changed the interface of populateRowsWith, added two more
//             unsigned Data Members, added support for isPopulated.
// 980306  dv  removed orientedRBPlacement and added its functionality to
//             RBPlacement
// 980309 pat  more const declarations
// 980427 pat  added findSubRow as a coreRow method.
// 990527 aec  split RBPlacement into a base class + derived class
//			base class does not use DB
// 990609 s.m  added spaceCellsWithPinDimensionsAlg1()
// 990610 s.m  added spaceCellsWithPinDimensionsAlg1()
// 000803 ilm  fixed savePlacement() to save cell names

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "RBPlacement.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

#include "ABKCommon/abkio.h"
#include "ABKCommon/abkfunc.h"
#include "ABKCommon/abktimer.h"
#include "Ctainers/bitBoard.h"
#include "Placement/transf.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::max;
using std::min;
using std::ofstream;
using std::ostringstream;
using std::ostream;
using std::pair;
using std::setw;
using std::stable_sort;
using std::string;
using std::vector;

#ifdef _MSC_VER
#ifndef rint
#define rint(a) floor((a) + 0.5)
#endif
#endif

void RBPlacement::splitCoreRowsByCell(const PlacementWOrient& pl,
                                      unsigned cellId) {
  unsigned angle = pl.getOrient(cellId).getAngle();
  bool notRotated = angle == 0 || angle == 180;
  double nodeWidth = notRotated ? _hgWDims->getNodeWidth(cellId)
                                : _hgWDims->getNodeHeight(cellId);
  double nodeHeight = notRotated ? _hgWDims->getNodeHeight(cellId)
                                 : _hgWDims->getNodeWidth(cellId);

  BBox obstacle;
  obstacle.add(pl[cellId].x, pl[cellId].y);
  obstacle.add(pl[cellId].x + nodeWidth, pl[cellId].y + nodeHeight);

  for (unsigned i = 0; i < _coreRows.size(); i++) {
    RBPCoreRow& cr = _coreRows[i];
    double coreRowYMin = cr.getYCoord();
    double coreRowYMax = cr.getYCoord() + cr.getHeight();

    if (greaterOrEqualDouble(obstacle.yMin, coreRowYMax)) continue;
    if (greaterOrEqualDouble(coreRowYMin, obstacle.yMax)) break;

    _splitRowByObstacle(cr, obstacle);
  }
}

void RBPlacement::removeFromCoreRows(const PlacementWOrient& pl,
                                     unsigned cellId) {
  unsigned angle = pl.getOrient(cellId).getAngle();
  bool notRotated = angle == 0 || angle == 180;
  double nodeWidth = notRotated ? _hgWDims->getNodeWidth(cellId)
                                : _hgWDims->getNodeHeight(cellId);
  double nodeHeight = notRotated ? _hgWDims->getNodeHeight(cellId)
                                 : _hgWDims->getNodeWidth(cellId);

  if (equalDouble(nodeWidth, 0.) || equalDouble(nodeHeight, 0.)) return;

  BBox obstacle;
  obstacle.add(pl[cellId].x, pl[cellId].y);
  obstacle.add(pl[cellId].x + nodeWidth, pl[cellId].y + nodeHeight);

  for (unsigned i = 0; i < _coreRows.size(); ++i) {
    RBPCoreRow& cr = _coreRows[i];
    double coreRowYMin = cr.getYCoord();
    double coreRowYMax = cr.getYCoord() + cr.getHeight();

    if (greaterOrEqualDouble(obstacle.yMin, coreRowYMax)) continue;
    if (greaterOrEqualDouble(coreRowYMin, obstacle.yMax)) break;

    cr.removeCell(pl, cellId);
  }
}

void RBPlacement::addToCoreRows(const PlacementWOrient& pl, unsigned cellId) {
  unsigned angle = pl.getOrient(cellId).getAngle();
  bool notRotated = angle == 0 || angle == 180;
  double nodeWidth = notRotated ? _hgWDims->getNodeWidth(cellId)
                                : _hgWDims->getNodeHeight(cellId);
  double nodeHeight = notRotated ? _hgWDims->getNodeHeight(cellId)
                                 : _hgWDims->getNodeWidth(cellId);

  if (equalDouble(nodeWidth, 0.) || equalDouble(nodeHeight, 0.)) return;

  BBox obstacle;
  obstacle.add(pl[cellId].x, pl[cellId].y);
  obstacle.add(pl[cellId].x + nodeWidth, pl[cellId].y + nodeHeight);

  for (unsigned i = 0; i < _coreRows.size(); ++i) {
    RBPCoreRow& cr = _coreRows[i];
    double coreRowYMin = cr.getYCoord();
    double coreRowYMax = cr.getYCoord() + cr.getHeight();

    if (greaterOrEqualDouble(obstacle.yMin, coreRowYMax)) continue;
    if (greaterOrEqualDouble(coreRowYMin, obstacle.yMax)) break;

    cr.addCell(pl, cellId);
  }
}

void RBPlacement::_splitRowByObstacle(RBPCoreRow& cr, const BBox& obstacle) {
  /*by sadya to make the halo sites parametrizable. by default 0*/
  const int numHaloSitesFront = 0;
  const int numHaloSitesEnd = 0;

  for (unsigned srIdx = 0; srIdx < cr.getNumSubRows(); ++srIdx) {
    RBPSubRow& sr = cr._subRows[srIdx];
    sr.clearAllCells();
  }

  for (unsigned srIdx = 0; srIdx < cr.getNumSubRows(); ++srIdx) {
    RBPSubRow& sr = cr._subRows[srIdx];
    double srXStart = sr.getXStart();
    double srXEnd = sr.getXEnd();

    if (obstacle.xMin >= srXEnd || obstacle.xMax <= srXStart) continue;

    // There are four ways the obstacle can intersect the subrow in x-dir
    // |----------| 		Obstacle
    //    |---|     		SubRow
    //
    //      |---|     	Obstacle
    //   |----------| 	SubRow
    //
    // |---|          	Obstacle
    //   |----------| 	SubRow
    //
    //            |---| 	Obstacle
    //   |----------|   	SubRow

    if ((obstacle.xMin <= srXStart) &&
        (obstacle.xMax >= srXEnd)) {  // Case1 : delete the subrow
      cr._subRows.erase(cr._subRows.begin() + srIdx);
      --srIdx;
    } else if ((obstacle.xMin <= srXStart) && (obstacle.xMax > srXStart) &&
               (obstacle.xMax <
                srXEnd)) {  // Case3: trim the front of the subrow.

      double newXStart = obstacle.xMax - srXStart;
      unsigned removedSites =
          static_cast<unsigned>(ceil(newXStart / sr.getSpacing())) +
          numHaloSitesEnd;
      newXStart = srXStart + removedSites * sr.getSpacing();

      sr._xStart = newXStart;
      sr.resetNumSites();
      if (sr.getNumSites() < 1) {
        cr._subRows.erase(cr._subRows.begin() + srIdx);
        --srIdx;
      }
    } else if ((obstacle.xMin > srXStart) && (obstacle.xMin < srXEnd) &&
               (obstacle.xMax >=
                srXEnd)) {  // Case4: trim the end of the subrow

      double newXEnd = srXEnd - obstacle.xMin;
      unsigned removedSites =
          static_cast<unsigned>(ceil(newXEnd / sr.getSpacing())) +
          numHaloSitesFront;
      newXEnd = srXEnd - removedSites * sr.getSpacing();

      sr._xEnd = newXEnd;
      sr.resetNumSites();
      if (sr.getNumSites() < 1) {
        cr._subRows.erase(cr._subRows.begin() + srIdx);
        --srIdx;
      }
    } else if ((obstacle.xMin > srXStart) &&
               (obstacle.xMax <
                srXEnd)) {  // Case2: split the subrow in two
                            // The new subRow will be the second one

      double newXStart = obstacle.xMax - srXStart;
      unsigned removedSites = static_cast<unsigned>(
          ceil(newXStart / sr.getSpacing()) + numHaloSitesEnd);
      newXStart = srXStart + removedSites * sr.getSpacing();
      RBPSubRow newSR(newXStart, sr.getXEnd(), cr);

      double newXEnd = srXEnd - obstacle.xMin;
      removedSites = static_cast<unsigned>(ceil(newXEnd / sr.getSpacing())) +
                     numHaloSitesFront;
      newXEnd = srXEnd - removedSites * sr.getSpacing();

      sr._xEnd = newXEnd;
      sr.resetNumSites();

      // must make this check before the insert as the
      // the insert will invalidate the reference "sr"
      bool deleteCurrSR = (sr.getNumSites() < 1);

      if (newSR.getNumSites() >= 1) {
        cr._subRows.insert(cr._subRows.begin() + srIdx + 1, newSR);
      }

      if (deleteCurrSR) {
        cr._subRows.erase(cr._subRows.begin() + srIdx);
        --srIdx;
      }
    }
  }
}

void RBPlacement::resetPlacementWOri(const PlacementWOrient& pl) {
  _placement = pl;
  setCellsInRows();
  populateCC();
}

void RBPlacement::resetPlacement(const Placement& pl) {
  _placement = PlacementWOrient(pl, vector<Orient>(pl.size()));
  setCellsInRows();
  populateCC();
}

void RBPlacement::resetPlacement() {
  setCellsInRows();
  populateCC();
}

void RBPlacement::updatePlacementWOri(const PlacementWOrient& pl) {
  _placement = pl;
}

void RBPlacement::setCellsInRows()
// Puts all the cellIds onto their appropriate subRows.
{
  _populated = true;
  _cellsNotInRows = 0;

  // clear all rows
  std::fill(_isInSubRow.begin(), _isInSubRow.end(), false);

  for (unsigned i = 0; i < _coreRows.size(); i++) {
    for (unsigned j = 0; j < _coreRows[i].getNumSubRows(); j++) {
      RBPSubRow& sRow = _coreRows[i][j];
      sRow.clearAllCells();
    }
    _coreRows[i].clearMovableStdCells(_isFixed, getHGraph());
  }

  for (unsigned i = 0; i < getNumCells(); i++) {
    if (_isCoreCell[i] && !_isFixed[i] && !_hgWDims->isNodeMacro(i)) {
      RBRowPtrs rowPtrs = findBothRows(_placement[i]);

      if (rowPtrs.first != NULL && rowPtrs.second != NULL) {
        const_cast<RBPSubRow*>(rowPtrs.second)->addCellToEnd(i);
        _isInSubRow[i] = true;
      } else {
        //          cout<<"Cell "<<_hgWDims->getNodeNameByIndex(i)<<" (" << i <<
        //          ") not placed "<<_placement[i]<<endl;
        _cellsNotInRows++;
      }
    }
  }
  if (_cellsNotInRows)
    cout << " RBPlacement: " << _cellsNotInRows
         << " core cell(s) not in rows or overlap with obstacles" << endl;

  sortSubRowsByX();
}

void RBPlacement::populateCC() {
  for (unsigned cellId = 0; cellId < _cellCoords.size(); ++cellId) {
    _cellCoords[cellId] = RBCellCoord();
  }

  for (unsigned crIdx = 0; crIdx < _coreRows.size(); ++crIdx) {
    const RBPCoreRow& cr = _coreRows[crIdx];
    for (unsigned srIdx = 0; srIdx < cr.getNumSubRows(); ++srIdx) {
      const RBPSubRow& sr = cr[srIdx];
      for (unsigned cOffset = 0; cOffset < sr.getNumCells(); ++cOffset) {
        _cellCoords[sr[cOffset]] = RBCellCoord(crIdx, srIdx, cOffset);
      }
    }
  }
  matchOrientsToRows();
}

void RBPlacement::setLocation(unsigned id, const Point& pt) {
  // puts cell <id> at location <pt>. updates the row structures
  if (_isInSubRow[id])  // need to remove the cell from the sub-row it was in
  {
    RBCellCoord& coord = _cellCoords[id];
    RBPSubRow& sr = _coreRows[coord.coreRowId][coord.subRowId];
    sr.removeCell(id);

    for (unsigned idx = 0; idx < sr.getNumCells(); ++idx) {
      _cellCoords[sr[idx]].cellOffset = idx;
    }

    _isInSubRow[id] = false;
    ++_cellsNotInRows;
  }

  _placement[id] = pt;
  RBRowPtrs newRows = findBothRows(pt);

  //    abkassert(newRows.first != NULL, "Invalid y-coordinate in setLocation");
  //    abkassert(newRows.second != NULL,"Invalid x-coordinate in setLocation");

  if (newRows.first != NULL && newRows.second != NULL) {
    _cellCoords[id].coreRowId = newRows.first - &(*coreRowsBegin());
    _cellCoords[id].subRowId =
        newRows.second - &(*newRows.first->subRowsBegin());
    RBPSubRow* newSR = const_cast<RBPSubRow*>(newRows.second);

    newSR->addCell(id);

    for (unsigned idx = 0; idx < newSR->getNumCells(); ++idx) {
      _cellCoords[(*newSR)[idx]].cellOffset = idx;
    }

    _isInSubRow[id] = true;
    --_cellsNotInRows;
  }
}

void RBPlacement::matchOrientToRow(unsigned id) {
  unsigned coreRowId = _cellCoords[id].coreRowId;
  Orient newCellOri(_coreRows[coreRowId].getOrientation().getAngle(),
                    _placement.getOrient(id).isFaceUp());
  _placement.setOrient(id, newCellOri);
}

void RBPlacement::matchOrientsToRows(void) {
  for (unsigned i = 0; i < getNumCells(); i++) {
    if (_isCoreCell[i] && !_isFixed[i] && _isInSubRow[i]) {
      matchOrientToRow(i);
    }
  }
}

void RBPlacement::sortSubRowsByX() {
  unsigned i, j;
  for (i = 0; i < _coreRows.size(); i++)
    for (j = 0; j < _coreRows[i].getNumSubRows(); j++)
      _coreRows[i][j].sortCellsByLoc();
}

bool RBPlacement::isConsistent() {
  bool looksOK = true;

  for (unsigned crId = 0; crId != _coreRows.size(); crId++) {
    RBPCoreRow& cr = _coreRows[crId];
    if (cr.getNumSubRows() == 0) {
      looksOK = false;
      cerr << "CR: " << crId << endl;
      cerr << "CoreRow does not have any SubRows" << endl;
      continue;
    }

    double prevXEnd = cr[0].getXStart();
    for (unsigned srId = 0; srId != cr.getNumSubRows(); srId++) {
      RBPSubRow& sr = cr[srId];
      if (sr._coreRow != &cr) {
        looksOK = false;
        cerr << "CR: " << crId << " SR: " << srId << endl;
        cerr << "SubRow does not point to the correct row" << endl;
        if (sr._coreRow == NULL)
          cerr << "SubRow has NULL coreRow pointer" << endl;
      }
      if (sr.getOrientation() != cr.getOrientation()) {
        looksOK = false;
        cerr << "Subrow orientation does not match its coreRow's" << endl;
      }
      if (sr.getXStart() < prevXEnd) {
        looksOK = false;
        cerr << "Mis-sorted x starts: " << sr.getXStart() << " " << prevXEnd
             << endl;
      }
      if (sr.getNumSites() <= 0) {
        looksOK = false;
        cerr << "Empty SubRow" << endl;
        cerr << sr << endl;
        cerr << "Site width " << sr.getSite().width << endl;
      }

      for (unsigned cOffset = 0; cOffset != sr.getNumCells(); cOffset++) {
        unsigned cellId = sr[cOffset];
        if (!_isInSubRow[cellId]) {
          looksOK = false;
          cerr << "subRow contains a cell marked as not in a SR" << endl;
        } else {
          const RBCellCoord& cc = _cellCoords[cellId];
          if ((cc.coreRowId != crId) || (cc.subRowId != srId) ||
              (cc.cellOffset != cOffset)) {
            looksOK = false;
            cerr << "CellCoord does not match cells actual location";
            cerr << endl;
          }
          if (_placement[cellId].y != sr.getYCoord()) {
            looksOK = false;
            cerr << "Cell in subRow with mismatched y" << endl;
          }
          if (_placement[cellId].x < sr.getXStart() ||
              _placement[cellId].x > sr.getXEnd()) {
            looksOK = false;
            cerr << "Cell in subRow with mismatched x" << endl;
          }
        }
      }
    }
  }

  return looksOK;
}

RBRowPtrs RBPlacement::findBothRows(const Point& point)
//   We assume that the RBPlacement contains _coreRows sorted by
//   y coordinate, and subrows sorted by x.  Return pointers
//   to the corerow and subrow that should contain the argument point.
{
  const RBPCoreRow* coreRowIdx = findCoreRow(point);

  if (coreRowIdx != NULL) {
    itRBPSubRow sr;
    for (sr = coreRowIdx->subRowsBegin(); sr != coreRowIdx->subRowsEnd();
         sr++) {
      if (point.x >= sr->getXStart() && point.x <= sr->getXEnd())
        return RBRowPtrs(
            std::pair<const RBPCoreRow*, const RBPSubRow*>(coreRowIdx, &(*sr)));
    }

    return RBRowPtrs(coreRowIdx, static_cast<const RBPSubRow*>(NULL));
  } else
    return RBRowPtrs(static_cast<const RBPCoreRow*>(NULL),
                     static_cast<const RBPSubRow*>(NULL));
}

const RBPSubRow* RBPlacement::findSubRow(const Point& point) const
//   Return a pointer to the subrow containing <point>
{
  const RBPCoreRow* coreRowIdx = findCoreRow(point);

  if (coreRowIdx == NULL) return NULL;

  itRBPSubRow sr;
  for (sr = coreRowIdx->subRowsBegin(); sr != coreRowIdx->subRowsEnd(); sr++) {
    if (point.x >= sr->getXStart() && point.x <= sr->getXEnd()) return &(*sr);
  }

  return NULL;
}

const RBPCoreRow* RBPlacement::findCoreRow(const Point& point) const
// Return a pointer to the _coreRow that should
// contain the argument point.
{
  RBPCoreRow cr;
  cr._y = point.y;

  pair<itRBPCoreRow, itRBPCoreRow> coreRowIts = std::equal_range(
      _coreRows.begin(), _coreRows.end(), cr, compareCoreRowY());

  if (coreRowIts.first == coreRowIts.second) {
    // not found
    return NULL;
  } else {
    return &(*coreRowIts.first);
  }
}

const unsigned RBPlacement::findCoreRowIdx(const Point& point) const {
  RBPCoreRow cr;
  cr._y = point.y;

  pair<itRBPCoreRow, itRBPCoreRow> coreRowIts = std::equal_range(
      _coreRows.begin(), _coreRows.end(), cr, compareCoreRowY());

  if (coreRowIts.first == coreRowIts.second) {
    // not found
    return UINT_MAX;
  } else {
    return (coreRowIts.first - _coreRows.begin());
  }
}

void RBPlacement::extractCellsFromSR(RBPSubRow& sr,
                                     std::vector<unsigned>& cellIds,
                                     std::vector<double>& offsets,
                                     double beginOffset, double endOffset) {
  cellIds.erase(cellIds.begin(), cellIds.end());
  offsets.erase(offsets.begin(), offsets.end());
  cellIds.reserve(sr.getNumCells());
  offsets.reserve(sr.getNumCells());
  bool foundEraseBegin = false;
  bool foundEraseEnd = false;
  std::vector<unsigned>::iterator eraseBegin = sr._cellIds.begin();
  std::vector<unsigned>::iterator eraseEnd = sr._cellIds.end();

  double epsilon = 0.01;

  beginOffset = max(0.0, beginOffset);
  endOffset = min(sr.getLength() + epsilon, endOffset);

  double rangeBegin = sr.getXStart() + beginOffset;
  double rangeEnd = sr.getXStart() + endOffset;

  for (unsigned i = 0; i < sr.getNumCells(); i++) {
    unsigned cellId = sr._cellIds[i];
    abkassert(_isInSubRow[cellId], "cell not marked as being in a sr");
    double placedLoc = _placement[cellId].x;

    if (placedLoc < rangeBegin) {
      abkassert(!foundEraseBegin, "cell left of range after start found");
      continue;
    } else if (placedLoc >= rangeEnd) {
      if (!foundEraseBegin)  // all are right of the range
        return;
      else if (!foundEraseEnd) {  // first one right of the range
        eraseEnd = sr._cellIds.begin() + i;
        foundEraseEnd = true;
      }

      break;
    } else  // in the range..remove it
    {
      if (!foundEraseBegin) {
        eraseBegin = sr._cellIds.begin() + i;
        foundEraseBegin = true;
      }
      abkassert(!foundEraseEnd, "sorting order of CellIds is incorrect");

      cellIds.push_back(cellId);
      offsets.push_back(placedLoc - sr.getXStart());

      _cellsNotInRows++;
      _isInSubRow[cellId] = false;
      // maintain _CellCoords??
    }
  }
  if (!foundEraseBegin) return;
  if (!foundEraseEnd) eraseEnd = sr._cellIds.end();

  sr._cellIds.erase(eraseBegin, eraseEnd);
}

void RBPlacement::extractCellsFromSR(RBPSubRow& sr,
                                     std::vector<unsigned>& cellIds,
                                     std::vector<unsigned>& siteOffsets,
                                     unsigned beginSite, unsigned endSite) {
  std::vector<double> offsets;
  extractCellsFromSR(sr, cellIds, offsets, beginSite * sr.getSpacing(),
                     endSite * sr.getSpacing());

  siteOffsets.erase(siteOffsets.begin(), siteOffsets.end());
  siteOffsets.insert(siteOffsets.begin(), offsets.size(), 0);
  for (unsigned i = 0; i < offsets.size(); i++)
    siteOffsets[i] = static_cast<unsigned>(rint(offsets[i] / sr.getSpacing()));
}

void RBPlacement::extractCellsFromSR(RBPSubRow& sr,
                                     std::vector<unsigned>& cellIds,
                                     double beginOffset, double endOffset) {
  cellIds.erase(cellIds.begin(), cellIds.end());
  cellIds.reserve(sr.getNumCells());
  bool foundEraseBegin = false;
  bool foundEraseEnd = false;
  std::vector<unsigned>::iterator eraseBegin = sr._cellIds.begin();
  std::vector<unsigned>::iterator eraseEnd = sr._cellIds.end();

  double epsilon = 0.01;

  beginOffset = max(0.0, beginOffset);
  endOffset = min(sr.getLength() + epsilon, endOffset);

  double rangeBegin = sr.getXStart() + beginOffset;
  double rangeEnd = sr.getXStart() + endOffset;

  for (unsigned i = 0; i < sr.getNumCells(); i++) {
    unsigned cellId = sr._cellIds[i];
    abkassert(_isInSubRow[cellId], "cell not marked as being in a sr");
    double placedLoc = _placement[cellId].x;

    if (placedLoc < rangeBegin) {
      abkfatal(!foundEraseBegin, "cell left of range after start found");
      continue;
    } else if (placedLoc >= rangeEnd) {
      if (!foundEraseBegin)  // all are right of the range
        return;
      else if (!foundEraseEnd) {  // first one right of the range
        eraseEnd = sr._cellIds.begin() + i;
        foundEraseEnd = true;
      }

      break;
    } else  // in the range..remove it
    {
      if (!foundEraseBegin) {
        eraseBegin = sr._cellIds.begin() + i;
        foundEraseBegin = true;
      }
      abkassert(!foundEraseEnd, "sorting order of CellIds is incorrect");

      cellIds.push_back(cellId);

      _cellsNotInRows++;
      _isInSubRow[cellId] = false;
      // maintain _CellCoords??
    }
  }
  if (!foundEraseBegin) return;
  if (!foundEraseEnd) eraseEnd = sr._cellIds.end();

  sr._cellIds.erase(eraseBegin, eraseEnd);
}

void RBPlacement::embedCellsInSR(RBPSubRow& sr,
                                 const std::vector<unsigned>& cellIds,
                                 const std::vector<double>& spacing) {
  abkfatal(spacing.size() == cellIds.size(),
           "cellIds and spacing do not have the same number of entries");

  unsigned i = 0;
  for (unsigned newCell = 0; newCell < cellIds.size(); newCell++) {
    unsigned cellId = cellIds[newCell];
    double newCellOffset = spacing[newCell];
    // find where cellId goes...

    for (; i < sr._cellIds.size(); i++)
      if (_placement[sr._cellIds[i]].x - sr.getXStart() > newCellOffset) break;

    sr._cellIds.insert(sr._cellIds.begin() + i, cellId);

    double loc = sr.getXStart() + newCellOffset;
    _placement[cellId].x = loc;
    _placement[cellId].y = sr.getYCoord();

    _cellsNotInRows--;
    _isInSubRow[cellId] = true;
    // maintain _CellCoords??
  }
}

void RBPlacement::embedCellsInSR(RBPSubRow& sr,
                                 const std::vector<unsigned>& cellIds,
                                 const std::vector<unsigned>& siteSpacing) {
  abkfatal(siteSpacing.size() == cellIds.size(),
           "cellIds and spacing do not have the same number of entries");

  unsigned i = 0;
  for (unsigned newCell = 0; newCell < cellIds.size(); newCell++) {
    unsigned cellId = sr._cellIds[newCell];
    double newCellOffset = siteSpacing[newCell] * sr.getSpacing();
    // find where cellId goes...

    for (; i < sr._cellIds.size(); i++)
      if (_placement[sr._cellIds[i]].x - sr.getXStart() > newCellOffset) break;

    sr._cellIds.insert(sr._cellIds.begin() + i, cellId);

    double loc = sr.getXStart() + newCellOffset;
    _placement[cellId].x = loc;
    _placement[cellId].y = sr.getYCoord();

    _cellsNotInRows--;
    _isInSubRow[cellId] = true;
    // maintain _CellCoords??
  }
}

void RBPlacement::embedCellsInSR(RBPSubRow& sr,
                                 const std::vector<unsigned>& cellIds,
                                 unsigned beginSite, unsigned endSite) {
  if (cellIds.size() == 0) return;

  sr._cellIds.reserve(sr.getNumCells() + cellIds.size());

  beginSite = max((unsigned)0, beginSite);
  endSite = min(sr.getNumSites(), endSite);

  double rangeBegin = sr.getXStart() + sr.getSpacing() * beginSite;
  double rangeEnd = sr.getXStart() + sr.getSpacing() * endSite;

  cout << "embedding range [" << rangeBegin << "," << rangeEnd << ")" << endl;

  double totalCellWidth = 0;
  unsigned i;
  for (i = 0; i < cellIds.size(); i++) {
    double cellWidth = _hgWDims->getNodeWidth(cellIds[i]);
    totalCellWidth += sr.getCoreRow().siteRoundUp(cellWidth);
  }

  double spaceBetween =
      ((rangeEnd - rangeBegin) - totalCellWidth) / cellIds.size();

  cout << "Total Site Width:      " << rangeEnd - rangeBegin << endl;
  cout << "Total Cell Width:      " << totalCellWidth << endl;
  cout << " ->Space Between Each: " << spaceBetween << endl << endl;

  double leftEdge = rangeBegin;
  double cumExtraSpace = 0;  // extra space is accumulated as cells are
                             // placed. when the extra space is > k*sitewidth,
                             // k sites are skipped (for k any non-zero
                             // integer)
  unsigned insPoint = 0;

  for (i = 0; i < cellIds.size(); i++) {
    unsigned newCellId = cellIds[i];

    for (; insPoint < sr._cellIds.size(); insPoint++) {
      if (_placement[sr._cellIds[insPoint]].x > leftEdge) break;
    }

    sr._cellIds.insert(sr._cellIds.begin() + insPoint, newCellId);
    _placement[newCellId].x = leftEdge;
    _placement[newCellId].y = sr.getYCoord();

    _cellsNotInRows--;
    _isInSubRow[newCellId] = true;
    // maintain _CellCoords??

    double cellWidth = _hgWDims->getNodeWidth(newCellId);
    leftEdge += sr.getCoreRow().siteRoundUp(cellWidth);

    cumExtraSpace += spaceBetween;

    int numSitesOverlap =
        static_cast<unsigned>(floor(cumExtraSpace / sr.getSpacing()));
    leftEdge += sr.getSpacing() * numSitesOverlap;
    cumExtraSpace -= sr.getSpacing() * numSitesOverlap;
  }
}

void RBPlacement::permuteCellsInSR(
    RBPSubRow& sr, unsigned firstCellOffset,
    const std::vector<unsigned>& cellPermutation) {
  abkassert(sr.getNumCells() - firstCellOffset >= cellPermutation.size(),
            "permuteCellsInSR has too many cells in the new ordering");

  std::vector<unsigned>::iterator permBegin =
      sr._cellIds.begin() + firstCellOffset;
  unsigned numToPermute = cellPermutation.size();

  double deltaSpace = sr.getSiteWidth();  // meaningless distance,
                                          // as a spacing function will be
                                          // called after permuting

  double leftEdge = _placement[*permBegin].x;

  vector<unsigned> newCellOrder(numToPermute);
  unsigned k;
  for (k = 0; k < numToPermute; k++)
    newCellOrder[k] = *(permBegin + cellPermutation[k]);

  for (k = 0; k < numToPermute; k++) {
    unsigned newCellId = newCellOrder[k];
    *(permBegin + k) = newCellId;
    _placement[newCellId].x = leftEdge;

    leftEdge += deltaSpace;
  }

  switch (_params.spaceCellsAlg) {
    case RBPlacement::Parameters::EQUAL_SPACE:
      spaceCellsEqually(sr, firstCellOffset, numToPermute);
      break;
    case RBPlacement::Parameters::WITH_PIN_ALG1:
      spaceCellsWithPinDimensionsAlg1(sr, firstCellOffset, numToPermute);
      break;
    case RBPlacement::Parameters::WITH_PIN_ALG2:
      spaceCellsWithPinDimensionsAlg2(sr, firstCellOffset, numToPermute);
      break;
    default:
      abkfatal(0, "unknown spacecells value");
  }
}

ostream& operator<<(ostream& out, const RBPlacement& arg) {
  unsigned i, j;
  unsigned cellCount = 0;

  out << "RBPlacement object\n";
  out << static_cast<const Placement>(static_cast<const PlacementWOrient>(arg));
  out << "  " << arg._coreRows.size() << " core rows: ";
  for (i = 0; i < arg._coreRows.size(); i++) {
    out << endl << "   " << arg._coreRows[i];
    for (j = 0; j < arg._coreRows[i].getNumSubRows(); j++)
      cellCount += arg._coreRows[i][j].getNumCells();
  }
  out << endl;
  out << "   containing " << cellCount << " cell Id(s)" << endl;
  out << endl;
  return out;
}

void RBPlacement::savePlacement(const char* plFileName, bool labelFixed) const {
  ofstream out(plFileName);

  out << "UCLA pl 1.0" << endl;
  out << TimeStamp() << User() << Platform() << endl << endl;

  for (unsigned i = 0; i < _placement.getSize(); i++) {
    const std::string namePtr = _hgWDims->getNodeNameByIndex(i).c_str();

    out << setw(8) << (!namePtr.empty() ? namePtr : "name unavailable") << "  "
        << setw(10) << _placement[i].x << "  " << setw(10)
        << (_placement[i].y / _params.yScale) << " : "
        << _placement.getOrient(i);
    if (labelFixed && _isFixed[i]) {
      out << "/FIXED";
    }
    out << endl;
  }

  out.close();
}

void RBPlacement::fillInDEF(const char* surrogate, const char* DEFName) const {
  ifstream olddef(surrogate);
  if (!olddef) {
    cerr << "Could not open surrogate DEF file";
    return;
  }

  std::string header;
  if (!std::getline(olddef, header) || header != "NOT A REAL DEF FILE!") {
    cerr << "Surrogate DEF corruption detected" << endl;
    return;
  }

  ofstream newdef(DEFName);
  if (!newdef) {
    cerr << "Could not open file for writing";
    return;
  }

  std::string line;
  while (std::getline(olddef, line)) {
    size_t locPos = line.find("UNPLACED");
    if (locPos != std::string::npos) {
      size_t cellStart = line.find_first_not_of(' ', 4);
      size_t cellEnd = line.find(' ', cellStart);
      std::string cellnameStr = line.substr(cellStart, cellEnd - cellStart);
      unsigned cellIdx = _hgWDims->getNodeByName(cellnameStr.c_str()).getIndex();
      ostringstream repl;
      repl << "  PLACED ( " << static_cast<long>(_placement[cellIdx].x) << "  "
           << static_cast<long>(_placement[cellIdx].y / _params.yScale) << " ) "
           << _placement.getOrient(cellIdx) << " ;\n";
      line.erase(locPos);
      line += repl.str();
    } else {
      locPos = line.find("PLACED");
      if (locPos != std::string::npos) {
        size_t cellStart = line.find_first_not_of(' ', 4);
        size_t cellEnd = line.find(' ', cellStart);
        std::string cellnameStr = line.substr(cellStart, cellEnd - cellStart);
        unsigned cellIdx =
            _hgWDims->getNodeByName(cellnameStr.c_str()).getIndex();
        ostringstream repl;
        repl << "PLACED ( " << static_cast<long>(_placement[cellIdx].x) << "  "
             << static_cast<long>(_placement[cellIdx].y / _params.yScale)
             << " ) " << _placement.getOrient(cellIdx) << " ;\n";
        line.erase(locPos);
        line += repl.str();
      }
    }

    newdef << line << '\n';
  }
}

void RBPlacement::savePlacementUnShredHW(const char* plFileName) const {
  ofstream out(plFileName);

  out << "UCLA pl 1.0" << endl;
  out << TimeStamp() << User() << Platform() << endl << endl;

  unsigned i = 0;
  for (i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (_placement[i].x > 1e308 || _placement[i].x > 1e308) continue;

    // HGFNode& node   = _hgWDims->getNodeByIdx(i);
    std::string nodeName = _hgWDims->getNodeNameByIndex(i).c_str();
    std::string actualName = nodeName;
    unsigned long pos = actualName.find('@'), pos1 = actualName.find('#');

    if (pos != std::string::npos && pos1 != std::string::npos)  // found
    {
      actualName.erase(pos);
      std::string subNodeName = nodeName;

      double avgYLoc = 0;
      double avgXLoc = 0;
      unsigned numSubNodes = 0;
      unsigned subVertIndex = 0;
      unsigned subHorizIndex = 0;
      std::string indices = nodeName;

      double nodeHeight = 0;
      double nodeWidth = 0;
      // Orientation nodeOrient = _placement.getOrient(i);

      vector<vector<unsigned> > subNodes;
      vector<unsigned> horizSubNodes;
      unsigned lastVertIdx = 0;
      unsigned currVertIdx = 0;

      while (i < _hgWDims->getNumNodes()) {
        // HGFNode& subNode   = _hgWDims->getNodeByIdx(i);
        nodeName = _hgWDims->getNodeNameByIndex(i).c_str();
        subNodeName = nodeName;
        pos = subNodeName.find('@');
        if (pos == std::string::npos) {
          subNodes.push_back(horizSubNodes);
          horizSubNodes.clear();
          break;
        }
        subNodeName.erase(pos);
        if (actualName != subNodeName) {
          subNodes.push_back(horizSubNodes);
          horizSubNodes.clear();
          break;
        }

        subNodeName = nodeName;
        pos = subNodeName.find('@');

        std::string rest(subNodeName, pos + 1);
        subNodeName.erase(pos);

        indices = rest;
        pos = indices.find('#');
        std::string rest2(indices, pos + 1);
        indices.erase(pos);
        istringstream ss1(indices.c_str());
        ss1 >> subVertIndex;
        istringstream ss2(rest2.c_str());
        ss2 >> subHorizIndex;

        currVertIdx = subVertIndex;
        if (currVertIdx == lastVertIdx)
          horizSubNodes.push_back(i);
        else {
          subNodes.push_back(horizSubNodes);
          horizSubNodes.clear();
          horizSubNodes.push_back(i);
          lastVertIdx = currVertIdx;
        }

        if (subHorizIndex == 0) nodeHeight += _hgWDims->getNodeHeight(i);
        if (subVertIndex == 0) nodeWidth += _hgWDims->getNodeWidth(i);

        avgYLoc += (_placement[i].y + _hgWDims->getNodeHeight(i) / 2);
        avgXLoc += (_placement[i].x + _hgWDims->getNodeWidth(i) / 2);
        numSubNodes++;
        i++;
      }
      if (horizSubNodes.size() != 0) {
        subNodes.push_back(horizSubNodes);
        horizSubNodes.clear();
      }
      i--;
      avgXLoc /= numSubNodes;
      avgYLoc /= numSubNodes;

      vector<int> orientScore(9, 0);
      unsigned k;
      for (k = 0; k < subNodes.size() - 1; ++k) {
        for (unsigned l = 0; l < subNodes[k].size() - 1; ++l) {
          clusterNodeOrient orient = getOrientSubNode(
              subNodes[k][l], subNodes[k][l + 1], subNodes[k + 1][l]);

          switch (orient) {
            case N:
              ++orientScore[N];
              break;
            case FN:
              ++orientScore[FN];
              break;
            case S:
              ++orientScore[S];
              break;
            case FS:
              ++orientScore[FS];
              break;
            case W:
              ++orientScore[W];
              break;
            case FW:
              ++orientScore[FW];
              break;
            case E:
              ++orientScore[E];
              break;
            case FE:
              ++orientScore[FE];
              break;
            case I:
              ++orientScore[I];
              break;
            default:
              cout << "Error in computing orientation " << orient << endl;
          };
        }
      }

      cout << actualName << endl;
      for (k = 0; k < 9; ++k) cout << setw(5) << toChar(clusterNodeOrient(k));
      cout << endl;

      int maxScore = -1;
      int maxScoreOrient = 0;
      for (unsigned k = 0; k < 8; ++k) {
        if (orientScore[k] > maxScore) {
          maxScore = orientScore[k];
          maxScoreOrient = k;
        }
        cout << setw(5) << orientScore[k];
      }

      // maxScoreOrient = N;

      double nodeYLoc;
      double nodeXLoc;
      if (maxScoreOrient == N || maxScoreOrient == S || maxScoreOrient == FN ||
          maxScoreOrient == FS) {
        nodeXLoc = avgXLoc - nodeWidth / 2;
        nodeYLoc = avgYLoc - nodeHeight / 2;
      } else {
        nodeXLoc = avgXLoc - nodeHeight / 2;
        nodeYLoc = avgYLoc - nodeWidth / 2;
      }

      out << actualName << "  ";
      out << setw(10) << nodeXLoc << "  ";
      out << setw(10) << nodeYLoc << " :  ";

      cout << setw(5) << orientScore[8];
      cout << endl;
      out << toChar(clusterNodeOrient(maxScoreOrient)) << endl;
      // out<<" N"<<endl;
    } else {
      out << _hgWDims->getNodeNameByIndex(i) << "  " << setw(10)
          << _placement[i].x << "  " << setw(10) << _placement[i].y << " : "
          << _placement.getOrient(i) << endl;
    }
  }
  out.close();
}

RBPlacement::~RBPlacement() {
  if (_hgWDims) delete _hgWDims;
  //  unsigned k;
  //  for(k=0; k!=allFileNamesInAUX.size(); k++) free(allFileNamesInAUX[k]);
}

void RBPlacement::permuteCellsAndWSInSR(
    RBPSubRow& sr, unsigned firstCellOffset,
    const std::vector<unsigned>& cellPermutation) {
  std::vector<unsigned>::iterator permBegin =
      sr._cellIds.begin() + firstCellOffset;
  unsigned numToPermute = cellPermutation.size();
  vector<unsigned> movables;
  vector<unsigned> WSpositions;
  vector<double> whiteSpaceWidths;
  unsigned totalMovables = 0;
  unsigned nodeCtr = firstCellOffset;
  unsigned totalCells = sr.getNumCells();

  int numSitesWS = 0;
  int sitesPerWS = 0;
  unsigned i;
  for (i = firstCellOffset; i < firstCellOffset + numToPermute; ++i) {
    if (nodeCtr == totalCells - 1)  // last cell of the subRow
    {
      movables.push_back(*(sr.cellIdsBegin() + nodeCtr));
      ++totalMovables;
      ++nodeCtr;
      if (totalMovables == numToPermute) break;

      // if whitespace present after last cell
      double currWSWidth =
          sr.getXEnd() - _placement[*(sr.cellIdsBegin() + nodeCtr - 1)].x -
          _hgWDims->getNodeWidth(*(sr.cellIdsBegin() + nodeCtr - 1));

      if (currWSWidth >= 2 * sr.getSpacing()) {
        numSitesWS = int(currWSWidth / sr.getSpacing());

        movables.push_back(UINT_MAX);
        sitesPerWS = int(floor(numSitesWS / 2 + 0.5));
        numSitesWS -= sitesPerWS;
        whiteSpaceWidths.push_back(sitesPerWS * sr.getSpacing());
        WSpositions.push_back(totalMovables);
        ++totalMovables;
        if (totalMovables == numToPermute) break;

        movables.push_back(UINT_MAX);
        sitesPerWS = numSitesWS;
        numSitesWS -= sitesPerWS;
        whiteSpaceWidths.push_back(sitesPerWS * sr.getSpacing());
        WSpositions.push_back(totalMovables);
        ++totalMovables;
        if (totalMovables == numToPermute) break;
      } else {
        movables.push_back(UINT_MAX);
        whiteSpaceWidths.push_back(currWSWidth);
        WSpositions.push_back(totalMovables);
        ++totalMovables;
        if (totalMovables == numToPermute) break;
      }
    } else {
      // if whitespace found
      if (_placement[*(sr.cellIdsBegin() + nodeCtr)].x +
              _hgWDims->getNodeWidth(*(sr.cellIdsBegin() + nodeCtr)) <
          _placement[*(sr.cellIdsBegin() + nodeCtr + 1)].x) {
        movables.push_back(*(sr.cellIdsBegin() + nodeCtr));
        ++totalMovables;
        ++nodeCtr;
        if (totalMovables == numToPermute) break;

        // better white space handling
        double currWSWidth =
            _placement[*(sr.cellIdsBegin() + nodeCtr)].x -
            _placement[*(sr.cellIdsBegin() + nodeCtr - 1)].x -
            _hgWDims->getNodeWidth(*(sr.cellIdsBegin() + nodeCtr - 1));

        if (currWSWidth >= 3 * sr.getSpacing()) {
          numSitesWS = int(currWSWidth / sr.getSpacing());

          movables.push_back(UINT_MAX);
          sitesPerWS = int(floor(numSitesWS / 3 + 0.5));
          numSitesWS -= sitesPerWS;
          whiteSpaceWidths.push_back(sitesPerWS * sr.getSpacing());
          WSpositions.push_back(totalMovables);
          ++totalMovables;
          if (totalMovables == numToPermute) break;

          movables.push_back(UINT_MAX);
          sitesPerWS = int(floor(numSitesWS / 2 + 0.5));
          numSitesWS -= sitesPerWS;
          whiteSpaceWidths.push_back(sitesPerWS * sr.getSpacing());
          WSpositions.push_back(totalMovables);
          ++totalMovables;
          if (totalMovables == numToPermute) break;

          movables.push_back(UINT_MAX);
          sitesPerWS = numSitesWS;
          numSitesWS -= sitesPerWS;
          whiteSpaceWidths.push_back(sitesPerWS * sr.getSpacing());
          WSpositions.push_back(totalMovables);
          ++totalMovables;
          if (totalMovables == numToPermute) break;
        } else if (currWSWidth < 3 * sr.getSpacing() &&
                   currWSWidth >= 2 * sr.getSpacing())

        // if(currWSWidth >=2*sr.getSpacing())
        {
          movables.push_back(UINT_MAX);
          whiteSpaceWidths.push_back(currWSWidth / 2);
          WSpositions.push_back(totalMovables);
          ++totalMovables;
          if (totalMovables == numToPermute) break;
          movables.push_back(UINT_MAX);
          whiteSpaceWidths.push_back(currWSWidth / 2);
          WSpositions.push_back(totalMovables);
          ++totalMovables;
          if (totalMovables == numToPermute) break;
        } else

        {
          movables.push_back(UINT_MAX);
          whiteSpaceWidths.push_back(currWSWidth);
          WSpositions.push_back(totalMovables);
          ++totalMovables;
          if (totalMovables == numToPermute) break;
        }
      } else {
        movables.push_back(*(sr.cellIdsBegin() + nodeCtr));
        ++totalMovables;
        ++nodeCtr;
        if (totalMovables == numToPermute) break;
      }
    }
  }

  abkassert(movables.size() == cellPermutation.size(),
            "permuteCellsInSR has too many cells in the new ordering");

  for (i = 0; i < WSpositions.size(); ++i) {
    for (unsigned j = 0; j < numToPermute; ++j) {
      if (WSpositions[i] == cellPermutation[j]) {
        WSpositions[i] = j;
        break;
      }
    }
  }

  double leftEdge = _placement[movables[0]].x;
  nodeCtr = 0;
  for (unsigned i = 0; i < numToPermute; ++i) {
    if (movables[cellPermutation[i]] != UINT_MAX) {
      _placement[movables[cellPermutation[i]]].x = leftEdge;
      leftEdge += _hgWDims->getNodeWidth(movables[cellPermutation[i]]);
      *(permBegin + nodeCtr) = movables[cellPermutation[i]];
      ++nodeCtr;
    } else {
      for (unsigned int j = 0; j < WSpositions.size(); ++j) {
        if (i == WSpositions[j]) {
          leftEdge += whiteSpaceWidths[j];
          break;
        }
      }
    }
  }
}

// added by sadya to update placement of ironed cells

void RBPlacement::updateIronedCells(const std::vector<unsigned>& movables,
                                    const Placement& soln, RBPSubRow& subrow1,
                                    RBPSubRow& subrow2) {
  double yCoord1 = subrow1.getYCoord();
  double yCoord2 = subrow2.getYCoord();

  unsigned i;
  for (i = 0; i < movables.size(); ++i) {
    if (movables[i] != UINT_MAX) {
      if (_placement[movables[i]].y == yCoord1) {
        subrow1.removeCell(movables[i]);
        _isInSubRow[movables[i]] = false;
      } else if (_placement[movables[i]].y == yCoord2) {
        subrow2.removeCell(movables[i]);
        _isInSubRow[movables[i]] = false;
      } else
        abkfatal(0, "invalid small placement problem\n");
    }
  }
  for (i = 0; i < movables.size(); i++) {
    if (movables[i] != UINT_MAX) {
      _placement[movables[i]] = Point(soln[i].x, soln[i].y);

      if (soln[i].y == yCoord1) {
        subrow1.addCell(movables[i]);
        _isInSubRow[movables[i]] = true;
      } else if (soln[i].y == yCoord2) {
        subrow2.addCell(movables[i]);
        _isInSubRow[movables[i]] = true;
      } else
        abkfatal(0, "invalid small placement solution\n");
    }
  }
}

// added by sadya to update placement of ironed cells with clustering

void RBPlacement::updateIronedCellsLR(
    const std::vector<std::vector<unsigned> >& movables,
    const Placement& soln, RBPSubRow& subrow1, RBPSubRow& subrow2) {
  double yCoord1 = subrow1.getYCoord();
  double yCoord2 = subrow2.getYCoord();

  unsigned i;
  for (i = 0; i < movables.size(); ++i) {
    if (movables[i][0] != UINT_MAX) {
      for (unsigned k = 0; k < movables[i].size(); ++k) {
        if (_placement[movables[i][k]].y == yCoord1) {
          subrow1.removeCell(movables[i][k]);
          _isInSubRow[movables[i][k]] = false;
        } else if (_placement[movables[i][k]].y == yCoord2) {
          subrow2.removeCell(movables[i][k]);
          _isInSubRow[movables[i][k]] = false;
        } else
          abkfatal(0, "invalid small placement problem\n");
      }
    }
  }
  for (i = 0; i < movables.size(); i++) {
    if (movables[i][0] != UINT_MAX) {
      double leftEdge = soln[i].x;
      for (unsigned k = 0; k < movables[i].size(); ++k) {
        _placement[movables[i][k]] = Point(leftEdge, soln[i].y);
        leftEdge += (*_hgWDims).getNodeWidth(movables[i][k]);

        if (soln[i].y == yCoord1) {
          if (_placement[movables[i][k]].x < subrow1.getXStart())
            _placement[movables[i][k]].x = subrow1.getXStart();
          if (_placement[movables[i][k]].x > subrow1.getXEnd())
            _placement[movables[i][k]].x = subrow1.getXEnd();
          subrow1.addCell(movables[i][k]);
          _isInSubRow[movables[i][k]] = true;
        } else if (soln[i].y == yCoord2) {
          if (_placement[movables[i][k]].x < subrow2.getXStart())
            _placement[movables[i][k]].x = subrow2.getXStart();
          if (_placement[movables[i][k]].x > subrow2.getXEnd())
            _placement[movables[i][k]].x = subrow2.getXEnd();
          subrow2.addCell(movables[i][k]);
          _isInSubRow[movables[i][k]] = true;
        } else
          abkfatal(0, "invalid small placement solution\n");
      }
    }
  }
}

// added by sadya to update placement of ironed cells with clustering R->L

void RBPlacement::updateIronedCellsRL(
    const std::vector<std::vector<unsigned> >& movables,
    const Placement& soln, RBPSubRow& subrow1, RBPSubRow& subrow2) {
  double yCoord1 = subrow1.getYCoord();
  double yCoord2 = subrow2.getYCoord();

  unsigned i;
  for (i = 0; i < movables.size(); ++i) {
    if (movables[i][0] != UINT_MAX) {
      for (unsigned k = 0; k < movables[i].size(); ++k) {
        if (_placement[movables[i][k]].y == yCoord1) {
          subrow1.removeCell(movables[i][k]);
          _isInSubRow[movables[i][k]] = false;
        } else if (_placement[movables[i][k]].y == yCoord2) {
          subrow2.removeCell(movables[i][k]);
          _isInSubRow[movables[i][k]] = false;
        } else
          abkfatal(0, "invalid small placement problem\n");
      }
    }
  }

  for (i = 0; i < movables.size(); i++) {
    if (movables[i][0] != UINT_MAX) {
      double leftEdge = soln[i].x;
      for (unsigned k = 0; k < movables[i].size(); ++k) {
        _placement[movables[i][k]] = Point(leftEdge, soln[i].y);

        leftEdge += (*_hgWDims).getNodeWidth(movables[i][k]);

        if (soln[i].y == yCoord1) {
          if (_placement[movables[i][k]].x < subrow1.getXStart())
            _placement[movables[i][k]].x = subrow1.getXStart();
          if (_placement[movables[i][k]].x > subrow1.getXEnd())
            _placement[movables[i][k]].x = subrow1.getXEnd();
          subrow1.addCell(movables[i][k]);
          _isInSubRow[movables[i][k]] = true;
        } else if (soln[i].y == yCoord2) {
          if (_placement[movables[i][k]].x < subrow2.getXStart())
            _placement[movables[i][k]].x = subrow2.getXStart();
          if (_placement[movables[i][k]].x > subrow2.getXEnd())
            _placement[movables[i][k]].x = subrow2.getXEnd();
          subrow2.addCell(movables[i][k]);
          _isInSubRow[movables[i][k]] = true;
        } else
          abkfatal(0, "invalid small placement solution\n");
      }
    }
  }
}

// added by sadya to get the BBox of the core layout
void RBPlacement::initBBox(void) {
  _coreBBox.clear();
  _fullBBox.clear();
  itRBPCoreRow itc;
  for (itc = coreRowsBegin(); itc != coreRowsEnd(); ++itc) {
    _coreBBox.add(itc->getXStart(), itc->getYCoord());
    _coreBBox.add(itc->getXEnd(), itc->getYCoord() + itc->getHeight());
    _fullBBox.add(itc->getXStart(), itc->getYCoord());
    _fullBBox.add(itc->getXEnd(), itc->getYCoord() + itc->getHeight());
  }

  // add to bbox all the fixed cells
  for (unsigned i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (isFixed(i)) {
      _fullBBox.add(_placement[i].x, _placement[i].y);
      _fullBBox.add(_placement[i].x + _hgWDims->getNodeWidth(i),
                    _placement[i].y + _hgWDims->getNodeHeight(i));
    }
  }
}

// added by sadya to get the BBox of the core layout
BBox RBPlacement::getBBox(bool addTerminals) const {
  BBox bbox;
  if (addTerminals)
    bbox = _fullBBox;
  else
    bbox = _coreBBox;
  return bbox;
}

clusterNodeOrient RBPlacement::getOrientSubNode(unsigned index,
                                                unsigned indexRight,
                                                unsigned indexTop) const {
  clusterNodeOrient orient = I;
  Point node(_placement[index]);
  Point nodeRight(_placement[indexRight]);
  Point nodeAbove(_placement[indexTop]);

  Point V2(nodeRight.x - node.x, nodeRight.y - node.y);
  Point V1(nodeAbove.x - node.x, nodeAbove.y - node.y);

  int quad1 = getQuad(V1);
  int quad2 = getQuad(V2);

  // cout<<"V1: "<<V1<<" quad1: "<<quad1<<"V2: "<<V2<<" quad2: "<<quad2<<endl;

  switch (quad2) {
    case 0:
      switch (quad1) {
        case 0:
        case 2:
          break;
        case 1:
          orient = N;
          break;
        case 3:
          orient = FS;
          break;
      };
      break;

    case 1:
      switch (quad1) {
        case 1:
        case 3:
          break;
        case 0:
          orient = FW;
          break;
        case 2:
          orient = W;
          break;
      };
      break;

    case 2:
      switch (quad1) {
        case 2:
        case 0:
          break;
        case 1:
          orient = FN;
          break;
        case 3:
          orient = S;
          break;
      };
      break;

    case 3:
      switch (quad1) {
        case 3:
        case 1:
          break;
        case 0:
          orient = E;
          break;
        case 2:
          orient = FE;
          break;
      };
      break;
  };

  // cout<<"quad1: "<<quad1<<" quad2: "<<quad2<<" orient
  // "<<toChar(orient)<<endl;
  /*
  if(_placement[index].x < _placement[indexRight].x &&
     _placement[index].y <= _placement[indexTop].y)
    orient = N;
  else if(_placement[index].x > _placement[indexRight].x &&
          _placement[index].y <= _placement[indexTop].y)
    orient = FN;
  else if(_placement[index].x <= _placement[indexRight].x &&
          _placement[index].y > _placement[indexTop].y)
    orient = S;
  else if(_placement[index].x >= _placement[indexRight].x &&
          _placement[index].y > _placement[indexTop].y)
    orient = FS;
  else if(_placement[index].x >= _placement[indexTop].x &&
          _placement[index].y <= _placement[indexRight].y)
    orient = W;
  else if(_placement[index].x <= _placement[indexTop].x &&
          _placement[index].y <= _placement[indexRight].y)
    orient = FW;
  else if(_placement[index].x <= _placement[indexTop].x &&
          _placement[index].y >= _placement[indexRight].y)
    orient = E;
  else if(_placement[index].x >= _placement[indexTop].x &&
          _placement[index].y >= _placement[indexRight].y)
    orient = FE;
  else
    cout<<"ERROR"<<endl;
  */
  // orient = N;

  return orient;
}

int RBPlacement::getQuad(Point& pt) const {
  if (pt.x > pt.y) {
    if (pt.x > -pt.y)
      return 0;
    else
      return 3;
  } else {
    if (pt.x > -pt.y)
      return 1;
    else
      return 2;
  }
}

const char* RBPlacement::toChar(clusterNodeOrient orient) const {
  switch (orient) {
    case N:
      return ("N");
      break;
    case FN:
      return ("FN");
      break;
    case S:
      return ("S");
      break;
    case FS:
      return ("FS");
      break;
    case W:
      return ("W");
      break;
    case FW:
      return ("FW");
      break;
    case E:
      return ("E");
      break;
    case FE:
      return ("FE");
      break;
    default:
      return ("I");
      break;
  };
  cout << orient << " ERROR" << endl;
  return ("N");
}

void RBPlacement::addObstacle(const BBox& obstacle) {
  double coreRowYMin, coreRowYMax;

  for (unsigned i = 0; i < _coreRows.size(); i++) {
    RBPCoreRow& cr = _coreRows[i];
    coreRowYMin = cr.getYCoord();
    coreRowYMax = cr.getYCoord() + cr.getHeight();

    if ((obstacle.yMin >= coreRowYMax) || (obstacle.yMax <= coreRowYMin))
      continue;

    _splitRowByObstacle(cr, obstacle);
  }

  // doing this screws up the coreRows pointers inside the sub rows...
  unsigned crId;
  for (crId = 0; crId < _coreRows.size(); crId++)
    if (_coreRows[crId].getNumSubRows() == 0) {
      _coreRows.erase(_coreRows.begin() + crId);
      crId--;
    }
  //...so we have to reset them manually
  for (crId = 0; crId < _coreRows.size(); crId++) {
    RBPCoreRow& cr = _coreRows[crId];
    unsigned srId;
    for (srId = 0; srId < cr.getNumSubRows(); srId++)
      cr._subRows[srId]._coreRow = &cr;
  }

  // finally, put the cells back
  setCellsInRows();
}

void RBPlacement::removeSitesBelowFixed(void) {
  for (unsigned i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (!isFixed(i)) continue;

    splitCoreRowsByCell(_placement, i);
    addToCoreRows(_placement, i);
  }

  // reading manually specified external obstacle file
  if (!_params.obstacleFile.empty()) {
    ifstream obstacles;
    obstacles.open(_params.obstacleFile.c_str());
    while (obstacles.good()) {
      double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
      obstacles >> x1 >> y1 >> x2 >> y2;
      if (!obstacles.good()) break;
      BBox obstacle(x1, y1, x2, y2);

      for (unsigned i = 0; i < _coreRows.size(); i++) {
        RBPCoreRow& cr = _coreRows[i];
        double coreRowYMin = cr.getYCoord();
        double coreRowYMax = cr.getYCoord() + cr.getHeight();

        if ((obstacle.yMin >= coreRowYMax) || (obstacle.yMax <= coreRowYMin))
          continue;

        _splitRowByObstacle(cr, obstacle);
      }
    }
  }

  // doing this screws up the coreRows pointers inside the sub rows...
  unsigned crId;
  for (crId = 0; crId < _coreRows.size(); crId++)
    if (_coreRows[crId].getNumSubRows() == 0) {
      _coreRows.erase(_coreRows.begin() + crId);
      crId--;
    }
  //...so we have to reset them manually
  for (crId = 0; crId < _coreRows.size(); crId++) {
    RBPCoreRow& cr = _coreRows[crId];
    unsigned srId;
    for (srId = 0; srId < cr.getNumSubRows(); srId++)
      cr._subRows[srId]._coreRow = &cr;
  }
}

void RBPlacement::removeSitesBelowMacros(void) {
  for (unsigned i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (!(isCoreCell(i) && _hgWDims->isNodeMacro(i) && !isFixed(i))) continue;

    splitCoreRowsByCell(_placement, i);
  }

  setCellsInRows();
  populateCC();
}

unsigned RBPlacement::getNumMacros(void) const {
  unsigned numMacros = 0;

  for (unsigned i = 0; i < _hgWDims->getNumNodes(); ++i) {
    if (isCoreCell(i) && _hgWDims->isNodeMacro(i) && !isFixed(i)) ++numMacros;
  }

  return numMacros;
}

unsigned RBPlacement::getNumFixedCore(void) const {
  unsigned numFixedCore = 0;

  for (unsigned i = 0; i < _hgWDims->getNumNodes(); ++i) {
    if (isCoreCell(i) && isFixed(i)) ++numFixedCore;
  }

  return numFixedCore;
}

unsigned RBPlacement::getNumFixedMacros(void) const {
  unsigned numFixedMacros = 0;

  for (unsigned i = 0; i < _hgWDims->getNumNodes(); ++i) {
    if (isCoreCell(i) && _hgWDims->isNodeMacro(i) && isFixed(i))
      ++numFixedMacros;
  }

  return numFixedMacros;
}

double RBPlacement::getMovableArea(void) const {
  double totalarea = 0;

  for (unsigned i = 0; i < _hgWDims->getNumNodes(); ++i) {
    if (isCoreCell(i) && !isFixed(i))
      totalarea += _hgWDims->getNodeWidth(i) * _hgWDims->getNodeHeight(i);
  }

  return totalarea;
}

void RBPlacement::updatePlacement(unsigned id, const Point& pt) {
  _placement[id] = pt;
}

// added by MRG to update placement of arbitrary list of cells
void RBPlacement::updateCells(const std::vector<unsigned>& movables,
                              const Placement& soln) {
  abkassert(movables.size() == soln.size(),
            "misuse of update cells, # of movables doesn't equal # of points");

  for (unsigned i = 0; i < movables.size(); ++i) {
    setLocation(movables[i], soln[i]);
  }
}

// added by MRG
// finds the current cellcoord at the point or makes one
unsigned RBPlacement::findCellIdx(const Point& point) {
  unsigned coreidx, subidx, i;
  for (i = 0; i < _coreRows.size(); i++) {
    if (_coreRows[i].getYCoord() == point.y) {
      coreidx = i;
      break;
    }
  }
  abkassert(point.y == _coreRows[coreidx].getYCoord(),
            "Point isn't in this row");
  for (i = 0; i < _coreRows[coreidx].getNumSubRows(); i++) {
    if (point.x >= _coreRows[coreidx][i].getXStart() &&
        point.x < _coreRows[coreidx][i].getXEnd()) {
      subidx = i;
      break;
    }
  }
  abkassert(point.x >= _coreRows[coreidx][subidx].getXStart() &&
                point.x < _coreRows[coreidx][subidx].getXEnd(),
            "Point isn't in this subrow");
  for (i = 0; i < _coreRows[coreidx][subidx].getNumCells(); i++) {
    // if there's a cell, return the cellcoord
    if (point.x == _placement[_coreRows[coreidx][subidx][i]].x)
      return (_coreRows[coreidx][subidx][i]);
    else if (point.x < _placement[_coreRows[coreidx][subidx][i]].x)
      break;
  }
  return (UINT_MAX);
}

bool RBPlacement::checkCC() {
  bool isOK = true;
  for (unsigned i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (isCoreCell(i)) {
      if (!_isInSubRow[i]) {
        cout << endl << "Cell not in valid row" << endl;
        cout << _hgWDims->getNodeNameByIndex(i);
        cout << _placement[i] << endl;
        isOK = false;
      }
    }
  }
  return (isOK);
}

void RBPlacement::getCellsInBBox(BBox& bbox, std::vector<unsigned>& cellIds) {
  Point bottomLeft(bbox.xMin, bbox.yMin);
  Point topRight(bbox.xMax, bbox.yMax);

  const RBPCoreRow* bottomCR = findCoreRow(bottomLeft);
  const RBPCoreRow* topCR = findCoreRow(topRight);

  if (bottomCR == NULL || topCR == NULL) return;

  unsigned crIdx;
  for (crIdx = 0; crIdx < _coreRows.size(); crIdx++) {
    const RBPCoreRow& cr = _coreRows[crIdx];
    if (bottomCR == &cr) break;
  }
  if (crIdx == _coreRows.size()) return;

  while (crIdx < _coreRows.size()) {
    const RBPCoreRow& cr = _coreRows[crIdx];

    if (cr.getYCoord() >= topRight.y) break;

    for (unsigned srIdx = 0; srIdx < cr.getNumSubRows(); ++srIdx) {
      const RBPSubRow& sr = cr[srIdx];
      double srXStart = sr.getXStart();
      double srXEnd = sr.getXEnd();
      if (srXStart < bottomLeft.x && srXEnd <= bottomLeft.x)
        continue;
      else if (srXStart >= topRight.x)
        continue;

      for (unsigned cOffset = 0; cOffset < sr.getNumCells(); cOffset++) {
        unsigned cellId = sr[cOffset];
        Point pl = _placement[cellId];
        double width = _hgWDims->getNodeWidth(cellId);
        if (pl.x >= bottomLeft.x && pl.x <= topRight.x)
          cellIds.push_back(cellId);
        else if (pl.x <= bottomLeft.x && (pl.x + width) >= bottomLeft.x)
          cellIds.push_back(cellId);
        if (pl.x >= topRight.x) break;
      }
    }
    if (topCR == &cr) break;
    ++crIdx;
  }
}

struct pair_less_mag
//: public binary_function< pair<unsigned, double>,
//		       pair<unsigned, double> , bool >
{
  bool operator()(pair<unsigned, double> x, pair<unsigned, double> y) {
    return (x.second < y.second);
  }
};

// To Do
void RBPlacement::assignPinsToCells(bool randomizeMacroDims) {
  itHGFEdgeLocal edge;

  double maxCoreHeight = getMaxHeightCoreRow();
  for (unsigned i = _hgWDims->getNumTerminals(); i < _hgWDims->getNumNodes();
       i++) {
    const HGFNode& node = _hgWDims->getNodeByIdx(i);

    if (randomizeMacroDims) {
      double origNodeHeight = _hgWDims->getNodeHeight(i);
      double origNodeWidth = _hgWDims->getNodeWidth(i);
      double nodeArea = origNodeHeight * origNodeWidth;
      if (origNodeHeight > maxCoreHeight) {
        double newAR = RandomDouble(0.5, 2.0);
        double newHeight = sqrt(nodeArea / newAR);
        if (newHeight <= maxCoreHeight) newHeight = maxCoreHeight * 2.0;
        newHeight = ceil(newHeight / maxCoreHeight) * maxCoreHeight;
        double newWidth = rint(nodeArea / newHeight);
        if (newWidth < 1) newWidth = 1;

        _hgWDims->setNodeWidth(newWidth, i);
        _hgWDims->setNodeHeight(newHeight, i);
      }
    }

    double nodeHeight = _hgWDims->getNodeHeight(i);
    double nodeWidth = _hgWDims->getNodeWidth(i);

    unsigned N = node.getDegree();
    vector<pair<unsigned, double> > pinAngles;

    Point center(_placement[i].x + nodeWidth * 0.5,
                 _placement[i].y + nodeHeight * 0.5);

    unsigned e = 0;  // edge offset
    for (edge = node.edgesBegin(); edge != node.edgesEnd(); edge++, e++) {
      Point optLoc(0, 0);
      unsigned edgeIdx = (*edge)->getIndex();
      unsigned nodeOffset = 0;
      for (itHGFNodeLocal currNode = (*edge)->nodesBegin();
           currNode != (*edge)->nodesEnd(); ++currNode, ++nodeOffset) {
        unsigned currNodeIdx = (*currNode)->getIndex();
        Point currNodePOffset = _hgWDims->nodesOnEdgePinOffset(
            nodeOffset, edgeIdx, _placement.getOrient(currNodeIdx));
        Point currNodePlacement = _placement[currNodeIdx] + currNodePOffset;
        optLoc += currNodePlacement;
      }

      Point pinLoc = (_placement[i] + _hgWDims->edgesOnNodePinOffset(
                                          e, i, _placement.getOrient(i)));

      optLoc -= pinLoc;
      optLoc.x /= (nodeOffset - 1);
      optLoc.y /= (nodeOffset - 1);

      Point netLoc = optLoc - center;
      double angle = atan2(netLoc.y, netLoc.x);
      pinAngles.push_back(pair<unsigned, double>(e, angle));
    }

    stable_sort(pinAngles.begin(), pinAngles.end(), pair_less_mag());

    double avgDistance = 2.0 * (nodeHeight + nodeWidth) / N;
    Point newNodeOffset(0, 0);
    vector<Point> newNodeOffsets(N, newNodeOffset);
    for (unsigned j = 0; j < pinAngles.size(); ++j) {
      unsigned index = pinAngles[j].first;
      double stLineDist = j * avgDistance;

      if (stLineDist <= nodeHeight * 0.5) {
        newNodeOffsets[index].x = 0;
        newNodeOffsets[index].y = stLineDist;
      } else if (stLineDist >= nodeHeight * 0.5 &&
                 stLineDist <= (nodeHeight * 0.5 + nodeWidth)) {
        newNodeOffsets[index].x = stLineDist - nodeHeight * 0.5;
        newNodeOffsets[index].y = 0;
      } else if (stLineDist >= (nodeHeight * 0.5 + nodeWidth) &&
                 stLineDist <= (1.5 * nodeHeight + nodeWidth)) {
        newNodeOffsets[index].x = nodeWidth;
        newNodeOffsets[index].y = stLineDist - (nodeHeight * 0.5 + nodeWidth);
      } else if (stLineDist >= (1.5 * nodeHeight + nodeWidth) &&
                 stLineDist <= (1.5 * nodeHeight + 2.0 * nodeWidth)) {
        newNodeOffsets[index].x =
            nodeWidth - (stLineDist - 1.5 * nodeHeight - nodeWidth);
        newNodeOffsets[index].y = nodeHeight;
      } else {
        newNodeOffsets[index].x = 0;
        newNodeOffsets[index].y =
            nodeHeight - (stLineDist - 1.5 * nodeHeight - 2.0 * nodeWidth);

        if (stLineDist > (2.0 * nodeHeight + 2.0 * nodeWidth))
          cout << "Error in assigning pins for Node "
               << _hgWDims->getNodeNameByIndex(node.getIndex()) << " "
               << stLineDist << " greater than perimeter "
               << (2.0 * nodeHeight + 2.0 * nodeWidth) << endl;
      }
    }

    e = 0;  // edge offset
    for (edge = node.edgesBegin(); edge != node.edgesEnd(); edge++, e++) {
      _hgWDims->updateEONPinOffsets(e, i, newNodeOffsets[e]);
    }
  }
}

double RBPlacement::getSitesArea() {
  double sitesArea = 0;
  for (unsigned i = 0; i < _coreRows.size(); i++) {
    for (unsigned j = 0; j < _coreRows[i].getNumSubRows(); j++) {
      RBPSubRow& sRow = _coreRows[i][j];
      sitesArea += sRow.getNumSites() * sRow.getHeight() * sRow.getSiteWidth();
    }
  }
  return (sitesArea);
}

void RBPlacement::alignCellsToRows() {
  for (unsigned i = 0; i < getNumCells(); i++) {
    if (_isCoreCell[i] && !_isFixed[i] && !_hgWDims->isNodeMacro(i)) {
      double cellWidth = _hgWDims->getNodeWidth(i);
      double cellHeight = _hgWDims->getNodeHeight(i);
      if (_placement[i].y < _coreBBox.yMin) _placement[i].y = _coreBBox.yMin;
      if (_placement[i].y > (_coreBBox.yMax - cellHeight))
        _placement[i].y = _coreBBox.yMax - cellHeight;
      if (_placement[i].x < _coreBBox.xMin) _placement[i].x = _coreBBox.xMin;
      if (_placement[i].x > (_coreBBox.xMax - cellWidth))
        _placement[i].x = _coreBBox.xMax - cellWidth;

      RBPCoreRow cr;
      const RBPCoreRow* lowerBound;

      cr._y = _placement[i].y;
      /*
      lowerBound = std::lower_bound(&(*_coreRows.begin()),
                                   &(*_coreRows.end()), cr,
                                   compareCoreRowY());
      */
      pair<const RBPCoreRow*, const RBPCoreRow*> coreRowIdx;
      coreRowIdx = std::equal_range(&(*_coreRows.begin()), &(*_coreRows.end()),
                                    cr, compareCoreRowY());
      if (coreRowIdx.first < coreRowIdx.second)
        lowerBound = coreRowIdx.first;
      else
        lowerBound = NULL;

      if (lowerBound != NULL) {
        _placement[i].y = lowerBound->getYCoord();

        itRBPSubRow sr;
        bool foundSR = false;
        for (sr = lowerBound->subRowsBegin(); sr != lowerBound->subRowsEnd();
             sr++) {
          if (_placement[i].x >= sr->getXStart() &&
              _placement[i].x <= sr->getXEnd()) {
            foundSR = true;
            break;
          }
        }

        if (!foundSR)  // not found. find the closest sub-row and fix it
        {
          int srId = -1;
          int minSrId = -1;
          bool beginEnd = false;  // closest to beginning or end of subrow
          double minDiff = DBL_MAX;
          double diff;
          for (sr = lowerBound->subRowsBegin(); sr != lowerBound->subRowsEnd();
               sr++) {
            ++srId;
            diff = fabs(_placement[i].x - sr->getXStart());
            if (diff < minDiff) {
              minSrId = srId;
              beginEnd = 0;
              minDiff = diff;
            }
            diff = fabs(_placement[i].x - (sr->getXEnd() - cellWidth));
            if (diff < minDiff) {
              minSrId = srId;
              beginEnd = 1;
              minDiff = diff;
            }
          }
          if (minSrId != -1) {
            foundSR = true;
            sr = lowerBound->subRowsBegin() + minSrId;
            if (beginEnd == false)
              _placement[i].x = sr->getXStart();
            else
              _placement[i].x = sr->getXEnd() - cellWidth;
          }
        }

        if (foundSR) {
          if (_placement[i].x > (sr->getXEnd() - cellWidth))
            _placement[i].x = sr->getXEnd() - cellWidth;
          if (_placement[i].x < sr->getXStart())
            _placement[i].x = sr->getXStart();
        }
      }
    }
  }
  setCellsInRows();
  populateCC();
}

double RBPlacement::getContainedSiteAreaInBBox(const BBox& bbox) const {
  BBox layoutBBox = getBBox(false);
  double containedArea = 0;
  BBox intersect = layoutBBox.intersect(bbox);
  if (intersect.getArea() == 0) return (0);

  itRBPCoreRow CR;
  itRBPSubRow SR;
  for (CR = coreRowsBegin(); CR != coreRowsEnd(); ++CR) {
    const RBPCoreRow& curRow = *CR;
    double yCoord = curRow.getYCoord();
    double rowHeight = curRow.getHeight();
    double thisRowHeight = rowHeight;

    if ((yCoord < intersect.yMin && (yCoord + rowHeight) <= intersect.yMin) ||
        (yCoord >= intersect.yMax))
      continue;

    if ((yCoord < intersect.yMin && (yCoord + rowHeight) >= intersect.yMin))
      thisRowHeight = (yCoord + rowHeight) - intersect.yMin;
    else if (yCoord < intersect.yMax && (yCoord + rowHeight) > intersect.yMax)
      thisRowHeight = intersect.yMax - yCoord;

    // this allows for fractional sites in the region
    double numSites = 0;
    double siteSpacing = curRow.getSpacing();
    double siteWidth = curRow.getSiteWidth();

    for (SR = (CR)->subRowsBegin(); SR != (CR)->subRowsEnd(); ++SR) {
      double xStart = SR->getXStart();
      double xEnd = SR->getXEnd();
      if (xStart < intersect.xMin) xStart = intersect.xMin;
      if (xEnd > intersect.xMax) xEnd = intersect.xMax;

      if (xEnd <= xStart) continue;
      double sitesInThisSR = (xEnd - xStart) / siteSpacing;
      numSites += sitesInThisSR;
    }
    containedArea += numSites * thisRowHeight * siteWidth;
  }
  return (containedArea);
}

void RBPlacement::initGroupRegionConstr(void) {
  for (unsigned i = 0; i < constraints.getSize(); ++i) {
    const Constraint& constraint = constraints[i];
    if (!constraint.getType().isRectRgnType()) continue;
    Constraint* nonConstConstraint = const_cast<Constraint*>(&constraint);
    RectRegionConstraint* rectConstr =
        static_cast<RectRegionConstraint*>(nonConstConstraint);
    if (rectConstr == NULL) {
      cout << "Not RectRegionConstraint. Actual Constraint Type"
           << (const char*)constraint.getType() << " " << endl;
      abkfatal(0,
               "Cannot convert type Constraint to type RectRegionConstraint\n")
    }
    BBox rgn = rectConstr->getBBox();

    const Mapping& mapping = constraint.getMapping();
    vector<unsigned> groupIds;
    for (unsigned j = 0; j < constraint.getSize(); ++j) {
      unsigned idx = mapping[j];
      groupIds.push_back(idx);
    }
    pair<BBox, vector<unsigned> > temp(rgn, groupIds);
    groupRegionConstr.push_back(temp);
    groupIds.clear();
  }

  // TO DO REMOVE THE FOLLOWING FAKE GROUP CONSTRUCT
  /*
  BBox layoutBBox = getBBox();
  BBox regionBBox;

  regionBBox.xMin = layoutBBox.xMin+layoutBBox.getWidth()/4;
  regionBBox.xMax = layoutBBox.xMax-layoutBBox.getWidth()/4;
  regionBBox.yMin = layoutBBox.yMin+layoutBBox.getHeight()/4;
  regionBBox.yMax = layoutBBox.yMax-layoutBBox.getHeight()/4;
  vector<unsigned> groupIds;
  for(unsigned i=0; i<rint(_hgWDims->getNumNodes()/10); i++)
    {
      unsigned ru = RandomUnsigned(_hgWDims->getNumTerminals(),
                                   _hgWDims->getNumNodes());
      groupIds.push_back(ru);
    }
  pair<BBox, vector<unsigned> > temp(regionBBox, groupIds);
  groupRegionConstr.push_back(temp);
  */

  cellToGrpMapping.resize(_hgWDims->getNumNodes(), UINT_MAX);
  for (unsigned j = 0; j < groupRegionConstr.size(); ++j) {
    for (unsigned k = 0; k < groupRegionConstr[j].second.size(); ++k)
      cellToGrpMapping[groupRegionConstr[j].second[k]] = j;
  }

  for (unsigned j = 0; j < groupRegionConstr.size(); ++j) {
    cout << "Group Region Constraint " << groupRegionConstr[j].first << endl;
    double cellArea = 0;
    for (unsigned k = 0; k < groupRegionConstr[j].second.size(); ++k) {
      // cout<<_hgWDims->getNodeByIdx(groupRegionConstr[j].second[k]).getName()<<"
      // ";
      unsigned id = groupRegionConstr[j].second[k];
      cellArea += _hgWDims->getNodeWidth(id) * _hgWDims->getNodeHeight(id);
    }
    double rgnArea = getContainedSiteAreaInBBox(groupRegionConstr[j].first);
    cout << "Region Utilization " << 100 * cellArea / rgnArea
         << " % : " << " Cell Area " << cellArea << ": Placeable sites area "
         << rgnArea << endl;
  }
}

void RBPlacement::initRegionUtilization(void) {
  for (unsigned i = 0; i < constraints.getSize(); ++i) {
    const Constraint& constraint = constraints[i];
    if (!constraint.getType().isUtilRectRgnType()) continue;
    Constraint* nonConstConstraint = const_cast<Constraint*>(&constraint);
    UtilRectRegionConstraint* utilRectConstr =
        static_cast<UtilRectRegionConstraint*>(nonConstConstraint);
    if (utilRectConstr == NULL) {
      cout << "Not UtilRectRegionConstraint. Actual Constraint Type"
           << (const char*)constraint.getType() << " " << endl;
      abkfatal(
          0,
          "Cannot convert type Constraint to type UtilRectRegionConstraint\n")
    }
    BBox rgn = utilRectConstr->getBBox();
    double utilization = utilRectConstr->getUtilization();
    pair<BBox, double> temp(rgn, utilization);
    regionUtilization.push_back(temp);
  }

  /* TODO REMOVE THE FOLLOWING FAKE CONSTRAINT
  BBox layoutBBox = getBBox();
  BBox regionBBox;

  regionBBox.xMin = layoutBBox.xMin+layoutBBox.getWidth()/4;
  regionBBox.xMax = layoutBBox.xMax-layoutBBox.getWidth()/4;
  regionBBox.yMin = layoutBBox.yMin+layoutBBox.getHeight()/4;
  regionBBox.yMax = layoutBBox.yMax-layoutBBox.getHeight()/4;
  pair <BBox , double> rgn1(regionBBox, 0.0);
  regionUtilization.push_back(rgn1);

  regionBBox.xMin = layoutBBox.xMin;
  regionBBox.xMax = layoutBBox.xMin+layoutBBox.getWidth()/3;
  regionBBox.yMin = layoutBBox.yMin;
  regionBBox.yMax = layoutBBox.yMin+layoutBBox.getHeight()/3;
  pair <BBox , double> rgn1(regionBBox, 0.05);
  regionUtilization.push_back(rgn1);

  regionBBox.xMin = layoutBBox.xMax-layoutBBox.getWidth()/3;
  regionBBox.xMax = layoutBBox.xMax;
  regionBBox.yMin = layoutBBox.yMax-layoutBBox.getHeight()/3;
  regionBBox.yMax = layoutBBox.yMax;
  pair <BBox , double> rgn2(regionBBox, 0.05);
  regionUtilization.push_back(rgn2);
  */

  for (unsigned i = 0; i < regionUtilization.size(); ++i) {
    cout << "Utilization Region Constraint " << regionUtilization[i].first
         << " : " << regionUtilization[i].second << endl;
  }
}

void RBPlacement::reshapePlacement(double tetherNewAR, double tetherNewWS,
                                   const char* newPlFileName) {
  BBox layoutBBox = getBBox(false);
  BBox newLayoutBBox = layoutBBox;
  double currWidth = layoutBBox.getWidth();
  double currHeight = layoutBBox.getHeight();
  double currAR = currWidth / currHeight;
  double layoutArea = getSitesArea();
  double totalNodesArea = _hgWDims->getNodesArea();
  if (tetherNewWS != -1) {
    double newLayoutArea = totalNodesArea / (1 - 0.01 * tetherNewWS);
    double newLayoutBBoxArea =
        layoutBBox.getArea() * newLayoutArea / layoutArea;
    double newWidth = sqrt(newLayoutBBoxArea * currAR);
    double newHeight = newWidth / currAR;
    newLayoutBBox.xMax = layoutBBox.xMin + newWidth;
    newLayoutBBox.yMax = layoutBBox.yMin + newHeight;
  }
  if (tetherNewAR != -1) {
    double newLayoutBBoxArea = newLayoutBBox.getArea();
    double newWidth = sqrt(newLayoutBBoxArea * tetherNewAR);
    double newHeight = newWidth / tetherNewAR;
    newLayoutBBox.xMax = newLayoutBBox.xMin + newWidth;
    newLayoutBBox.yMax = newLayoutBBox.yMin + newHeight;
  }

  if (tetherNewWS != -1 || tetherNewAR != -1) {
    BBoxToBBox reshapeBBox(layoutBBox, newLayoutBBox);
    for (unsigned i = 0; i < _hgWDims->getNumNodes(); i++) {
      if (_placement[i].x > 1e308 || _placement[i].x > 1e308) continue;
      Point cellLoc = _placement[i];
      Point newCellLoc = reshapeBBox.apply(cellLoc);
      _placement[i] = newCellLoc;
    }
  }
  if (tetherNewAR != -1 && newPlFileName != NULL) {
    // have to parse the placement file here, as it may have
    // node names, which only the HGraph knows.
    PlacementWOrient newPlacement =
        PlacementWOrient(_hgWDims->getNumNodes(), Point(DBL_MAX, DBL_MAX));

    ifstream plFile(newPlFileName);
    abkfatal2(plFile, " Could not open ", newPlFileName);
    cout << " Reading " << newPlFileName << " ... " << endl;

    string nodeName;
    Point pt;
    string oriStr;
    int lineno = 0;

    plFile >> needcaseword("UCLA") >> needcaseword("pl") >> needword("1.0") >>
        skiptoeol;
    plFile >> eathash(lineno);

    while (!plFile.eof()) {
      plFile >> eathash(lineno) >> nodeName >> pt.x >> pt.y >> eatblank;
      unsigned nodeId = _hgWDims->getNodeByName(nodeName.c_str()).getIndex();
      newPlacement[nodeId] = pt;

      if (plFile.peek() == ':') {
        plFile >> needword(":") >> oriStr;
        newPlacement.setOrient(nodeId, Orientation(oriStr.c_str()));
      } else if (plFile.peek() == '/') {
        plFile.get();
        plFile >> oriStr;
        newPlacement.setOrient(nodeId, Orientation(oriStr.c_str()));
      }
      plFile >> skiptoeol >> eathash(lineno);
    }
    plFile.close();

    // now that we have new terminal locations we can rotate each node
    // according to new terminal locations
    Point center = newLayoutBBox.getGeomCenter();
    vector<pair<unsigned, double> > initTermAngles;
    vector<double> newTermAngles;

    for (unsigned i = 0; i < _hgWDims->getNumTerminals(); i++) {
      Point cellLoc = _placement[i] - center;
      double angle = atan2(cellLoc.y, cellLoc.x);
      initTermAngles.push_back(pair<unsigned, double>(i, angle));

      cellLoc = newPlacement[i] - center;
      angle = atan2(cellLoc.y, cellLoc.x);
      newTermAngles.push_back(angle);
    }

    stable_sort(initTermAngles.begin(), initTermAngles.end(), pair_less_mag());
    for (unsigned i = 0; i < initTermAngles.size(); ++i) {
      // cout<<_hgWDims->getNodeByIdx(initTermAngles[i].first).getName()<<"
      // "<<initTermAngles[i].second<<"
      // "<<newTermAngles[initTermAngles[i].first]<<endl;
    }

    BBox centeredBBox = newLayoutBBox;
    centeredBBox.xMin -= center.x;
    centeredBBox.xMax -= center.x;
    centeredBBox.yMin -= center.y;
    centeredBBox.yMax -= center.y;
    double Q1 = atan2(centeredBBox.yMin, centeredBBox.xMax);
    double Q2 = atan2(centeredBBox.yMax, centeredBBox.xMax);
    double Q3 = atan2(centeredBBox.yMax, centeredBBox.xMin);
    double Q4 = atan2(centeredBBox.yMin, centeredBBox.xMin);

    for (unsigned i = 0; i < _hgWDims->getNumNodes(); i++) {
      if (_placement[i].x > 1e308 || _placement[i].x > 1e308) continue;
      Point cellLoc = _placement[i] - center;
      if (cellLoc.x == 0 && cellLoc.y == 0) continue;

      double initAngle = atan2(cellLoc.y, cellLoc.x);

      // now find between which terminals this point rests
      unsigned j;
      for (j = 0; j < _hgWDims->getNumTerminals() - 1; j++) {
        if (initTermAngles[j].second <= initAngle &&
            initTermAngles[j + 1].second >= initAngle)
          break;
      }
      unsigned term1 = j;
      unsigned term2 = j + 1;
      if (j == _hgWDims->getNumTerminals() - 1) term2 = 0;

      double alpha1 = initTermAngles[term1].second;
      double alpha2 = initTermAngles[term2].second;
      double beta1 = newTermAngles[initTermAngles[term1].first];
      double beta2 = newTermAngles[initTermAngles[term2].first];

      double initDiff = alpha2 - alpha1;
      double newDiff = beta2 - beta1;
      double angleDiff = initAngle - alpha1;

      if (initDiff <= -Pi) initDiff += 2 * Pi;
      if (newDiff <= -Pi) newDiff += 2 * Pi;
      if (angleDiff <= -Pi) angleDiff += 2 * Pi;

      double newAngle = beta1 + (angleDiff / initDiff) * newDiff;

      if (newAngle > Pi)
        newAngle -= 2 * Pi;
      else if (newAngle <= -Pi)
        newAngle += 2 * Pi;

      double blkWidth = centeredBBox.getWidth();
      double blkHeight = centeredBBox.getHeight();

      double initR = sqrt(cellLoc.x * cellLoc.x + cellLoc.y * cellLoc.y);

      double length;

      // cout<<_hgWDims->getNodeByIdx(i).getName()<<" "<<initAngle<<"
      // "<<newAngle<<endl;

      if (Q1 <= initAngle && initAngle <= Q2)
        length = (blkWidth / 2.0) / cos(initAngle);
      else if (Q2 <= initAngle && initAngle <= Q3)
        length = (blkHeight / 2.0) / cos(Pi / 2.0 - initAngle);
      else if (Q3 <= initAngle || initAngle <= Q4)
        length = (blkWidth / 2.0) / cos(Pi - initAngle);
      else
        length = (blkHeight / 2.0) / cos(Pi / 2.0 + initAngle);

      initR /= length;

      if (Q1 <= newAngle && newAngle <= Q2)
        length = (blkWidth / 2.0) / cos(newAngle);
      else if (Q2 <= newAngle && newAngle <= Q3)
        length = (blkHeight / 2.0) / cos(Pi / 2.0 - newAngle);
      else if (Q3 <= newAngle || newAngle <= Q4)
        length = (blkWidth / 2.0) / cos(Pi - newAngle);
      else
        length = (blkHeight / 2.0) / cos(Pi / 2.0 + newAngle);

      double newR = initR * length;

      Point newCellLoc;
      newCellLoc.x = newR * cos(newAngle);
      newCellLoc.y = newR * sin(newAngle);

      _placement[i] = newCellLoc + center;
    }
  }
}

void RBPlacement::markMacrosAsFixed() {
  unsigned j;
  for (j = 0; j < _hgWDims->getNumNodes(); j++) {
    if (_hgWDims->isNodeMacro(j)) _isFixed[j] = true;
  }
}

// void RBPlacement::unmarkMacrosAsFixed()
//{
//   unsigned j;
//   for(j=0; j<_hgWDims->getNumNodes(); j++)
//     {
//       if(_hgWDims->isNodeMacro(j))
//	_isFixed[j] = false;
//     }
// }

void RBPlacement::markMacrosAsSoft(double minAR, double maxAR) {
  _hgWDims->markMacrosAsSoft(minAR, maxAR);
}

unsigned RBPlacement::getNumSoft() const {
  unsigned numSoft = 0;
  for (unsigned n = 0; n < _hgWDims->getNumNodes(); n++)
    if (_hgWDims->isNodeSoft(n)) numSoft++;

  return numSoft;
}

void RBPlacement::backUpCoreRows() {
  unsigned numRows = _coreRows.size();
  _backUpCoreRows.clear();
  _backUpCoreRows.reserve(numRows);

  for (unsigned rNum = 0; rNum < numRows; rNum++) {
    _backUpCoreRows.push_back(RBPCoreRow(
        _coreRows[rNum].getYCoord(), _coreRows[rNum].getOrientation(),
        _coreRows[rNum].getSite(), _placement, _coreRows[rNum].getSpacing()));

    RBPCoreRow& newRow = _backUpCoreRows.back();

    newRow._allCells = _coreRows[rNum]._allCells;
    newRow._whitespaceLocs = _coreRows[rNum]._whitespaceLocs;
    newRow._whitespaceLocsPopulated = _coreRows[rNum]._whitespaceLocsPopulated;

    for (itRBPSubRow sRowIt = _coreRows[rNum].subRowsBegin();
         sRowIt != _coreRows[rNum].subRowsEnd(); ++sRowIt)
      newRow.appendNewSubRow(sRowIt->getXStart(), sRowIt->getXEnd());
  }
}

void RBPlacement::reinstateCoreRows() {
  unsigned numRows = _backUpCoreRows.size();
  _coreRows.clear();
  _coreRows.reserve(numRows);

  for (unsigned rNum = 0; rNum < numRows; rNum++) {
    _coreRows.push_back(RBPCoreRow(_backUpCoreRows[rNum].getYCoord(),
                                   _backUpCoreRows[rNum].getOrientation(),
                                   _backUpCoreRows[rNum].getSite(), _placement,
                                   _backUpCoreRows[rNum].getSpacing()));

    RBPCoreRow& newRow = _coreRows.back();

    newRow._allCells = _backUpCoreRows[rNum]._allCells;
    newRow._whitespaceLocs = _backUpCoreRows[rNum]._whitespaceLocs;
    newRow._whitespaceLocsPopulated =
        _backUpCoreRows[rNum]._whitespaceLocsPopulated;

    for (itRBPSubRow sRowIt = _backUpCoreRows[rNum].subRowsBegin();
         sRowIt != _backUpCoreRows[rNum].subRowsEnd(); ++sRowIt)
      newRow.appendNewSubRow(sRowIt->getXStart(), sRowIt->getXEnd());
  }
}

void RBPlacement::reinstateCoreRowsAndRepopulate() {
  reinstateCoreRows();
  setCellsInRows();
  populateCC();
}

void RBPlacement::addRowsByObstacleFrmBackup(const BBox& obstacle) {
  unsigned numRows = _coreRows.size();
  for (unsigned rNum = 0; rNum < numRows; rNum++) {
    RBPCoreRow& curRow = _coreRows[rNum];
    if (curRow.getYCoord() >= obstacle.yMin &&
        curRow.getYCoord() <= obstacle.yMax) {
      RBPCoreRow& backUpRow = _backUpCoreRows[rNum];
      for (itRBPSubRow sRowIt = backUpRow.subRowsBegin();
           sRowIt != backUpRow.subRowsEnd(); ++sRowIt) {
        double xStart = sRowIt->getXStart();
        double xEnd = sRowIt->getXEnd();
        if (xStart <= obstacle.xMin && xEnd >= obstacle.xMax)
          curRow.appendNewSubRow(obstacle.xMin, obstacle.xMax);
        else if (xStart <= obstacle.xMin && xEnd <= obstacle.xMax &&
                 xEnd >= obstacle.xMin)
          curRow.appendNewSubRow(obstacle.xMin, xEnd);
        else if (xStart >= obstacle.xMin && xStart <= obstacle.xMax &&
                 xEnd >= obstacle.xMax)
          curRow.appendNewSubRow(xStart, obstacle.xMax);
      }
      std::vector<RBPSubRow>::iterator sRowBegin = curRow.NCsubRowsBegin();
      std::vector<RBPSubRow>::iterator sRowEnd = curRow.NCsubRowsEnd();
      std::sort(sRowBegin, sRowEnd);

      // now check for consistency and merge subrows if possible
      for (std::vector<RBPSubRow>::iterator sRowIt = curRow.NCsubRowsBegin();
           sRowIt != (curRow.NCsubRowsEnd() - 1); ++sRowIt) {
        std::vector<RBPSubRow>::iterator sRowIt1 = sRowIt + 1;
        if (sRowIt->getXEnd() > sRowIt1->getXStart()) {
          cout << *sRowIt << " " << *sRowIt1 << endl;
          abkfatal(0, "Sub-rows Inconsistent while merging from backup rows\n");
        }
        if (sRowIt->getXEnd() == sRowIt1->getXStart()) {
          // need to merge these sub-rows
          double xStart = sRowIt->getXStart();
          double xEnd = sRowIt1->getXEnd();
          sRowIt = curRow.eraseSubRow(sRowIt, sRowIt1);

          curRow.appendNewSubRow(xStart, xEnd);
          sRowBegin = curRow.NCsubRowsBegin();
          sRowEnd = curRow.NCsubRowsEnd();
          std::sort(sRowBegin, sRowEnd);
          sRowIt = curRow.NCsubRowsBegin();
        }
        if (sRowIt == (curRow.NCsubRowsEnd() - 1)) break;
      }
    }
  }
}

RBPlacement& operator<<(RBPlacement& rbplace, const vector<Orient>& ori) {
  for (unsigned i = 0; i < ori.size(); i++) {
    rbplace.setOrient(i, ori[i]);
  }
  return rbplace;
}

void RBPlacement::markMacros(void) {
  _hgWDims->markNodesAsMacro(FLT_MAX, getMaxHeightCoreRow());
}

RBPlacement::RBPlacement(HGraphWDimensions* hg, const PlacementWOrient& pl,
                         const std::vector<RBPCoreRow>* cr,
                         const Parameters& params)
    : RBPlacement(hg, pl, std::vector<RBPCoreRow>(cr->begin(), cr->end()),
                  params) {}

RBPlacement::RBPlacement(HGraphWDimensions* hg, const PlacementWOrient& pl,
                         const std::vector<RBPCoreRow>& cr,
                         const Parameters& params)
    : _hgWDims(hg),
      _ownHGWDimsAuxName(true),
      _populated(true),
      _placement(pl),
      _isInSubRow(hg->getNumNodes(), false),
      _cellCoords(hg->getNumNodes()),
      _isFixed(hg->getNumNodes(), false),
      _isCoreCell(hg->getNumNodes(), true),
      _storElt(hg->getNumNodes(), false),
      _params(params) {
  for (itHGFNodeGlobal it = hg->terminalsBegin(); it != hg->terminalsEnd();
       ++it) {
    unsigned term_idx = (*it)->getIndex();
    _isFixed[term_idx] = true;
    _isCoreCell[term_idx] = false;
    _isInSubRow[term_idx] = true;
  }

  // the following code is used to make the core rows passed
  // in refer to the local PlacementWOrient which is a
  // member of the RBPlacement class, rather than whatever
  // external PlacementWOrient was used to construct the rows
  // vector.
  // note that this is ugly, but it is necessary because
  // of the tight coupling between rows, placements, and rbplaces
  // this can be fixed by removing the dependency of
  // rows on placements.
  _coreRows.reserve(cr.size());
  for (std::vector<RBPCoreRow>::const_iterator it = cr.begin();
       it != cr.end(); ++it) {
    RBPCoreRow r(*it, _placement);
    _coreRows.push_back(r);
  }

  markMacros();
  // the following two lines are to ensure that node weights relate to area
  _hgWDims->expandStdCellWidthToFitSiteWidth(_coreRows[0].getSpacing());
  _hgWDims->setNodeWeightsToAreas(_coreRows[0].getSpacing(),
                                  _coreRows[0].getHeight());
  removeSitesBelowFixed();
  backUpCoreRows();
  setCellsInRows();
  populateCC();
  initBBox();
}

void RBPlacement::findBestGreedySwap(greedySwapData& cells) {
  cells.pos1 = _placement[cells.c1];
  cells.pos2 = _placement[cells.c2];
  cells.improvement = -DBL_MAX;

  double improvement = 0.;

  double space1 = getAvailSpace(cells.c1);
  double space2 = getAvailSpace(cells.c2);
  double len1 = _hgWDims->getNodeWidth(cells.c1);
  double len2 = _hgWDims->getNodeWidth(cells.c2);

  // swap 1 and 2, valid with 2 and 3 cells
  if (len1 <= space2 && len2 <= space1) {
    improvement =
        -calcSwapCost(cells.c1, cells.c2, UINT_MAX, _placement[cells.c2],
                      _placement[cells.c1], cells.pos3);
    if (improvement > cells.improvement) {
      cells.improvement = improvement;
      cells.pos1 = _placement[cells.c2];
      cells.pos2 = _placement[cells.c1];
    }
  }
  // the rest are only valid with 3 cells
  if (cells.c3 == UINT_MAX) return;
  cells.pos3 = _placement[cells.c3];
  double space3 = getAvailSpace(cells.c3);
  double len3 = _hgWDims->getNodeWidth(cells.c3);

  // swap 1 and 3
  if (len1 <= space3 && len3 <= space1) {
    improvement =
        -calcSwapCost(cells.c1, cells.c3, UINT_MAX, _placement[cells.c3],
                      _placement[cells.c1], cells.pos2);
    if (improvement > cells.improvement) {
      cells.improvement = improvement;
      cells.pos1 = _placement[cells.c3];
      cells.pos2 = _placement[cells.c2];
      cells.pos3 = _placement[cells.c1];
    }
  }

  // swap 2 and 3
  if (len2 <= space3 && len3 <= space2) {
    improvement =
        -calcSwapCost(cells.c2, cells.c3, UINT_MAX, _placement[cells.c3],
                      _placement[cells.c2], cells.pos1);
    if (improvement > cells.improvement) {
      cells.improvement = improvement;
      cells.pos1 = _placement[cells.c1];
      cells.pos2 = _placement[cells.c3];
      cells.pos3 = _placement[cells.c2];
    }
  }

  // rotate 1 backward
  if (len1 <= space3 && len2 <= space1 && len3 <= space2) {
    improvement =
        -calcSwapCost(cells.c1, cells.c2, cells.c3, _placement[cells.c3],
                      _placement[cells.c1], _placement[cells.c2]);
    if (improvement > cells.improvement) {
      cells.improvement = improvement;
      cells.pos1 = _placement[cells.c3];
      cells.pos2 = _placement[cells.c1];
      cells.pos3 = _placement[cells.c2];
    }
  }

  // rotate 1 forward
  if (len1 <= space2 && len2 <= space3 && len3 <= space1) {
    improvement =
        -calcSwapCost(cells.c1, cells.c2, cells.c3, _placement[cells.c2],
                      _placement[cells.c3], _placement[cells.c1]);
    if (improvement > cells.improvement) {
      cells.improvement = improvement;
      cells.pos1 = _placement[cells.c2];
      cells.pos2 = _placement[cells.c3];
      cells.pos3 = _placement[cells.c1];
    }
  }
}

void RBPlacement::fillInAffectedCells(unsigned cell,
                                      std::vector<unsigned>& affectedCells,
                                      const std::vector<unsigned>& allcells) {
  static BitBoard temp(0);
  if (temp.getSize() < _placement.size()) temp = BitBoard(_placement.size());
  temp.clear();

  itHGFEdgeLocal e;
  itHGFNodeLocal n;
  const HGFNode& node = _hgWDims->getNodeByIdx(cell);
  for (e = node.edgesBegin(); e != node.edgesEnd(); ++e) {
    for (n = (*e)->nodesBegin(); n != (*e)->nodesEnd(); ++n) {
      const unsigned& idx = (*n)->getIndex();
      if (std::binary_search(allcells.begin(), allcells.end(), idx)) {
        temp.setBit(idx);
      }
    }
  }

  const auto& setBits = temp.getIndicesOfSetBits();
  affectedCells.assign(setBits.begin(), setBits.end());
  std::sort(affectedCells.begin(), affectedCells.end());
}

std::pair<double, double> RBPlacement::getSurroundingSpaceSpan(
    unsigned cell) const {
  std::pair<double, double> span(0., 0.);

  abkassert(_isInSubRow[cell],
            "getSurroundingSpaceSpan() should not be called on a cell not in a "
            "subrow.");

  const RBCellCoord& cellCoord = _cellCoords[cell];

  abkassert(
      _coreRows[cellCoord.coreRowId][cellCoord.subRowId].getNumCells() > 0,
      "Cell is in a subrow that contains no cells?");

  if (cellCoord.cellOffset == 0) {
    span.first = _coreRows[cellCoord.coreRowId][cellCoord.subRowId].getXStart();
  } else {
    unsigned frontneighbor = _coreRows[cellCoord.coreRowId][cellCoord.subRowId]
                                      [cellCoord.cellOffset - 1];
    span.first =
        _placement[frontneighbor].x + _hgWDims->getNodeWidth(frontneighbor);
  }

  if (cellCoord.cellOffset >=
      _coreRows[cellCoord.coreRowId][cellCoord.subRowId].getNumCells() - 1) {
    span.second = _coreRows[cellCoord.coreRowId][cellCoord.subRowId].getXEnd();
  } else {
    unsigned backneighbor = _coreRows[cellCoord.coreRowId][cellCoord.subRowId]
                                     [cellCoord.cellOffset + 1];
    span.second = _placement[backneighbor].x;
  }

  return span;
}

double RBPlacement::getAvailSpace(unsigned cell) const {
  double space = 0.;

  abkassert(_isInSubRow[cell],
            "getAvailSpace() should not be called on a cell not in a subrow.");

  const RBCellCoord& cellCoord = _cellCoords[cell];

  abkassert(
      _coreRows[cellCoord.coreRowId][cellCoord.subRowId].getNumCells() > 0,
      "Cell is in a subrow that contains no cells?");

  if (cellCoord.cellOffset >=
      _coreRows[cellCoord.coreRowId][cellCoord.subRowId].getNumCells() - 1) {
    space = _coreRows[cellCoord.coreRowId][cellCoord.subRowId].getXEnd() -
            _placement[cell].x;
  } else {
    unsigned neighbor = _coreRows[cellCoord.coreRowId][cellCoord.subRowId]
                                 [cellCoord.cellOffset + 1];
    space = _placement[neighbor].x - _placement[cell].x;
  }

  return space;
}

bool RBPlacement::movableCheck(
    const std::vector<unsigned>& possible_move) const {
  if (possible_move.size() == 2) {
    unsigned cell1 = possible_move[0];
    unsigned cell2 = possible_move[1];

    if (!isCoreCell(cell1) || isFixed(cell1) || !_isInSubRow[cell1] ||
        _hgWDims->isNodeMacro(cell1))
      return false;
    if (!isCoreCell(cell2) || isFixed(cell2) || !_isInSubRow[cell2] ||
        _hgWDims->isNodeMacro(cell2))
      return false;

    return true;
  } else {
    unsigned cell1 = possible_move[0];
    unsigned cell2 = possible_move[1];
    unsigned cell3 = possible_move[2];

    bool movable1 = isCoreCell(cell1) && !isFixed(cell1) &&
                    _isInSubRow[cell1] && !_hgWDims->isNodeMacro(cell1);
    bool movable2 = isCoreCell(cell2) && !isFixed(cell2) &&
                    _isInSubRow[cell2] && !_hgWDims->isNodeMacro(cell2);
    bool movable3 = isCoreCell(cell3) && !isFixed(cell3) &&
                    _isInSubRow[cell3] && !_hgWDims->isNodeMacro(cell3);

    if (!movable1 && !movable2) return false;
    if (!movable1 && !movable3) return false;
    if (!movable2 && !movable3) return false;

    return true;
  }
}

void RBPlacement::greedySwapping(MaxMem* maxMem) {
  if (_params.verb.getForActions())
    cout << " Running greedy swapping ..." << endl;

  Timer tm;

  double startWL = evalHPWL(true);
  if (_params.verb.getForMajStats()) {
    cout << "  HPWL before greedy swapping is " << startWL << endl;
  }

  std::set<std::vector<unsigned> > potential_moves;
  std::vector<unsigned> base, move;
  base.reserve(3);
  move.reserve(3);
  unsigned cell1, cell2, cell3;
  bool extended;
  itHGFEdgeLocal e;
  itHGFNodeLocal c;

  for (itHGFEdgeGlobal n = getHGraph().edgesBegin();
       n != getHGraph().edgesEnd(); ++n) {
    if ((*n)->getDegree() == 2) {
      extended = false;
      cell1 = (*((*n)->nodesBegin()))->getIndex();
      cell2 = (*((*n)->nodesBegin() + 1))->getIndex();

      base.clear();
      base.push_back(cell1);
      base.push_back(cell2);

      const HGFNode& node1 = _hgWDims->getNodeByIdx(cell1);
      for (e = node1.edgesBegin(); e != node1.edgesEnd(); ++e) {
        if ((*e)->getDegree() != 2) continue;
        for (c = (*e)->nodesBegin(); c != (*e)->nodesEnd(); ++c) {
          if ((*c)->getIndex() != cell1 && (*c)->getIndex() != cell2) {
            move = base;
            move.push_back((*c)->getIndex());
            std::sort(move.begin(), move.end());
            std::vector<unsigned>::iterator new_end =
                std::unique(move.begin(), move.end());
            move.erase(new_end, move.end());
            std::set<std::vector<unsigned> >::iterator new_pos =
                potential_moves.lower_bound(move);
            if ((new_pos == potential_moves.end() || *new_pos != move) &&
                movableCheck(move)) {
              potential_moves.insert(new_pos, move);
            }
            extended = true;
          }
        }
      }

      const HGFNode& node2 = _hgWDims->getNodeByIdx(cell2);
      for (e = node2.edgesBegin(); e != node2.edgesEnd(); ++e) {
        if ((*e)->getDegree() != 2) continue;
        for (c = (*e)->nodesBegin(); c != (*e)->nodesEnd(); ++c) {
          if ((*c)->getIndex() != cell1 && (*c)->getIndex() != cell2) {
            move = base;
            move.push_back((*c)->getIndex());
            std::sort(move.begin(), move.end());
            std::vector<unsigned>::iterator new_end =
                std::unique(move.begin(), move.end());
            move.erase(new_end, move.end());
            std::set<std::vector<unsigned> >::iterator new_pos =
                potential_moves.lower_bound(move);
            if ((new_pos == potential_moves.end() || *new_pos != move) &&
                movableCheck(move)) {
              potential_moves.insert(new_pos, move);
            }
            extended = true;
          }
        }
      }

      if (!extended) {
        move = base;
        std::sort(move.begin(), move.end());
        std::vector<unsigned>::iterator new_end =
            std::unique(move.begin(), move.end());
        move.erase(new_end, move.end());
        std::set<std::vector<unsigned> >::iterator new_pos =
            potential_moves.lower_bound(move);
        if ((new_pos == potential_moves.end() || *new_pos != move) &&
            movableCheck(move)) {
          potential_moves.insert(new_pos, move);
        }
      }
    } else if ((*n)->getDegree() == 3) {
      cell1 = (*((*n)->nodesBegin()))->getIndex();
      cell2 = (*((*n)->nodesBegin() + 1))->getIndex();
      cell3 = (*((*n)->nodesBegin() + 2))->getIndex();

      move.clear();
      move.push_back(cell1);
      move.push_back(cell2);
      move.push_back(cell3);

      std::sort(move.begin(), move.end());
      std::vector<unsigned>::iterator new_end =
          std::unique(move.begin(), move.end());
      move.erase(new_end, move.end());
      std::set<std::vector<unsigned> >::iterator new_pos =
          potential_moves.lower_bound(move);
      if ((new_pos == potential_moves.end() || *new_pos != move) &&
          movableCheck(move)) {
        potential_moves.insert(new_pos, move);
      }
    }
  }
  base.clear();
  move.clear();

  maxMem->update("Greedy Swapping: After building potential moves");

  if (_params.verb.getForMajStats() > 3) {
    cout << "Number of potential moves " << potential_moves.size() << endl;
  }

  bool movable1, movable2, movable3;
  double len1, len2, len3, space1, space2, space3;
  std::vector<greedySwapData> themoves;
  themoves.reserve(potential_moves.size());

  for (std::set<std::vector<unsigned> >::iterator s = potential_moves.begin();
       s != potential_moves.end(); ++s) {
    if (s->size() == 2) {
      cell1 = (*s)[0];
      cell2 = (*s)[1];

      if (!isCoreCell(cell1) || isFixed(cell1) || !_isInSubRow[cell1] ||
          _hgWDims->isNodeMacro(cell1))
        continue;
      if (!isCoreCell(cell2) || isFixed(cell2) || !_isInSubRow[cell2] ||
          _hgWDims->isNodeMacro(cell2))
        continue;

      len1 = _hgWDims->getNodeWidth(cell1);
      len2 = _hgWDims->getNodeWidth(cell2);
      space1 = getAvailSpace(cell1);
      space2 = getAvailSpace(cell2);

      if (len1 > space2 || len2 > space1) continue;
      themoves.push_back(greedySwapData(cell1, cell2, UINT_MAX));
    } else if (s->size() == 3) {
      cell1 = (*s)[0];
      cell2 = (*s)[1];
      cell3 = (*s)[2];

      movable1 = isCoreCell(cell1) && !isFixed(cell1) && _isInSubRow[cell1] &&
                 !_hgWDims->isNodeMacro(cell1);
      movable2 = isCoreCell(cell2) && !isFixed(cell2) && _isInSubRow[cell2] &&
                 !_hgWDims->isNodeMacro(cell2);
      movable3 = isCoreCell(cell3) && !isFixed(cell3) && _isInSubRow[cell3] &&
                 !_hgWDims->isNodeMacro(cell3);

      if (!movable1 && !movable2) continue;
      if (!movable1 && !movable3) continue;
      if (!movable2 && !movable3) continue;

      len1 = _hgWDims->getNodeWidth(cell1);
      len2 = _hgWDims->getNodeWidth(cell2);
      len3 = _hgWDims->getNodeWidth(cell3);

      if (!movable1) {
        space2 = getAvailSpace(cell2);
        space3 = getAvailSpace(cell3);
        if (len2 <= space3 && len3 <= space2)
          themoves.push_back(greedySwapData(cell2, cell3, UINT_MAX));
      } else if (!movable2) {
        space1 = getAvailSpace(cell1);
        space3 = getAvailSpace(cell3);
        if (len1 <= space3 && len3 <= space1)
          themoves.push_back(greedySwapData(cell1, cell3, UINT_MAX));
      } else if (!movable3) {
        space1 = getAvailSpace(cell1);
        space2 = getAvailSpace(cell2);
        if (len1 <= space2 && len2 <= space1)
          themoves.push_back(greedySwapData(cell1, cell2, UINT_MAX));
      } else {
        themoves.push_back(greedySwapData(cell1, cell2, cell3));
      }
    }
  }

  maxMem->update("Greedy Swapping: After weeding out impossible moves");

  if (_params.verb.getForMajStats() > 3) {
    cout << "Number of initial possible moves " << themoves.size() << endl;
  }

  double bestImprovement = 0.;
  double improveCutoff = 1.e-6;
  unsigned i, bestMove = UINT_MAX;
  std::vector<greedySwapData> tempmoves;
  tempmoves.swap(themoves);

  for (i = 0; i < tempmoves.size(); ++i) {
    findBestGreedySwap(tempmoves[i]);
    if (tempmoves[i].improvement != -DBL_MAX &&
        tempmoves[i].improvement >= improveCutoff) {
      if (tempmoves[i].improvement > bestImprovement) {
        bestMove = themoves.size();
        bestImprovement = tempmoves[i].improvement;
      }
      themoves.push_back(tempmoves[i]);
    }
  }

  maxMem->update(
      "Greedy Swapping: After weeding out moves with not enough improvement");

  if (_params.verb.getForMajStats() > 3) {
    cout << "Number of final possible moves " << themoves.size() << endl;
  }

  std::vector<std::vector<unsigned> > affectedCells(getNumCells());
  {
    std::vector<unsigned> allcells;
    allcells.reserve(3 * themoves.size());
    for (i = 0; i < themoves.size(); ++i) {
      allcells.push_back(themoves[i].c1);
      allcells.push_back(themoves[i].c2);
      if (themoves[i].c3 != UINT_MAX) allcells.push_back(themoves[i].c3);
    }

    std::sort(allcells.begin(), allcells.end());
    std::vector<unsigned>::iterator new_end =
        std::unique(allcells.begin(), allcells.end());
    allcells.erase(new_end, allcells.end());

    for (i = 0; i < allcells.size(); ++i) {
      fillInAffectedCells(allcells[i], affectedCells[allcells[i]], allcells);
    }

    maxMem->update("Greedy Swapping: After filling in affected cells");
  }

  if (_params.verb.getForMajStats() > 3) {
    cout << "Finished filling in affected cells" << endl;
  }

  if (0.01 * bestImprovement > improveCutoff)
    improveCutoff = 0.01 * bestImprovement;
  if (1.e-6 * startWL > improveCutoff) improveCutoff = 1.e-6 * startWL;
  if (_params.verb.getForMajStats() > 3) {
    cout << "Improvement cutoff " << improveCutoff << endl
         << " which is "
         << floor(improveCutoff * 1000000. / startWL + 0.5) / 10000.
         << "% of starting HPWL" << endl;
  }

  unsigned oldMove, twocellmoves = 0, threecellmoves = 0;
  while (bestMove != UINT_MAX) {
    setLocation(themoves[bestMove].c1, themoves[bestMove].pos1);
    setLocation(themoves[bestMove].c2, themoves[bestMove].pos2);
    if (themoves[bestMove].c3 == UINT_MAX)
      ++twocellmoves;
    else {
      setLocation(themoves[bestMove].c3, themoves[bestMove].pos3);
      ++threecellmoves;
    }

    if (_params.verb.getForMajStats() > 3) {
      cout << "Just made move number " << twocellmoves + threecellmoves << endl;
      cout << "   Improvement " << bestImprovement << endl;
      cout << "   Two-cell moves so far " << twocellmoves << endl;
      cout << "   Three-cell moves so far " << threecellmoves << endl;
    }
    if (bestImprovement < improveCutoff) break;

    oldMove = bestMove;
    bestMove = UINT_MAX;
    bestImprovement = 1.e-6;
    unsigned updates = 0;
    for (i = 0; i < themoves.size(); ++i) {
      if (binary_search(affectedCells[themoves[oldMove].c1].begin(),
                        affectedCells[themoves[oldMove].c1].end(),
                        themoves[i].c1) ||
          binary_search(affectedCells[themoves[oldMove].c1].begin(),
                        affectedCells[themoves[oldMove].c1].end(),
                        themoves[i].c2) ||
          binary_search(affectedCells[themoves[oldMove].c1].begin(),
                        affectedCells[themoves[oldMove].c1].end(),
                        themoves[i].c3) ||
          binary_search(affectedCells[themoves[oldMove].c2].begin(),
                        affectedCells[themoves[oldMove].c2].end(),
                        themoves[i].c1) ||
          binary_search(affectedCells[themoves[oldMove].c2].begin(),
                        affectedCells[themoves[oldMove].c2].end(),
                        themoves[i].c2) ||
          binary_search(affectedCells[themoves[oldMove].c2].begin(),
                        affectedCells[themoves[oldMove].c2].end(),
                        themoves[i].c3) ||
          (themoves[oldMove].c3 != UINT_MAX &&
           (binary_search(affectedCells[themoves[oldMove].c3].begin(),
                          affectedCells[themoves[oldMove].c3].end(),
                          themoves[i].c1) ||
            binary_search(affectedCells[themoves[oldMove].c3].begin(),
                          affectedCells[themoves[oldMove].c3].end(),
                          themoves[i].c2) ||
            binary_search(affectedCells[themoves[oldMove].c3].begin(),
                          affectedCells[themoves[oldMove].c3].end(),
                          themoves[i].c3)))) {
        findBestGreedySwap(themoves[i]);
        ++updates;
      }
      if (themoves[i].improvement > bestImprovement) {
        bestMove = i;
        bestImprovement = themoves[i].improvement;
      }
    }
    if (_params.verb.getForMajStats() > 3) {
      cout << "   Updates needed " << updates << endl;
    }
  }

  double endWL = evalHPWL(true);

  tm.stop();

  if (_params.verb.getForMajStats()) {
    cout << " Greedy swapping took " << tm << endl;
    cout << "  Two-cell moves performed: " << twocellmoves << endl;
    cout << "  Three-cell moves performed: " << threecellmoves << endl;
    cout << "  HPWL after greedy swapping is " << endWL << endl;
    double impr = (startWL == 0.) ? 0. : (startWL - endWL) / startWL;
    if (impr < 0.)
      cout << "  % increase in HPWL is " << floor(-impr * 100000. + 0.5) / 1000.
           << endl;
    else
      cout << "  % decrease in HPWL is " << floor(impr * 100000. + 0.5) / 1000.
           << endl;
  }

  matchOrientsToRows();
}

RBPlacement::RBPlacement(unsigned numCells, const Placement& pl,
                         const vector<Orient>& ori, const Parameters& params)
    : _hgWDims(NULL),
      _ownHGWDimsAuxName(true),
      _populated(false),
      _placement(pl, ori),
      _isInSubRow(numCells, false),
      _cellCoords(numCells),
      _isFixed(numCells, false),
      _isCoreCell(numCells, false),
      _storElt(numCells, false),
      _params(params) {}

RBPlacement::RBPlacement(unsigned numCells, const PlacementWOrient& pl,
                         const Parameters& params)
    : _hgWDims(NULL),
      _ownHGWDimsAuxName(true),
      _populated(false),
      _placement(pl),
      _isInSubRow(numCells, false),
      _cellCoords(numCells),
      _isFixed(numCells, false),
      _isCoreCell(numCells, false),
      _storElt(numCells, false),
      _params(params) {}

RBPlacement::RBPlacement()
    : _hgWDims(NULL),
      _ownHGWDimsAuxName(true),
      _populated(false),
      _cellsNotInRows(0),
      _placement(0),
      _params() {
  initBBox();
}

// added by royj just for fun!
RBPlacement::RBPlacement(const RBPlacement& orig)
    : _ownHGWDimsAuxName(orig._ownHGWDimsAuxName),
      _populated(orig._populated),
      _cellsNotInRows(orig._cellsNotInRows),
      _placement(orig._placement),
      _isInSubRow(orig._isInSubRow),
      _cellCoords(orig._cellCoords),
      _isFixed(orig._isFixed),
      _isCoreCell(orig._isCoreCell),
      _storElt(orig._storElt),
      _params(orig._params),
      _coreBBox(orig._coreBBox),
      _fullBBox(orig._fullBBox),
      origFileName(orig.origFileName),
      dirName(orig.dirName),
      allFileNamesInAUX(orig.allFileNamesInAUX) {
  if (orig._hgWDims == NULL) {
    _hgWDims = NULL;
  } else {
    _hgWDims = new HGraphWDimensions(*orig._hgWDims);
  }
  _coreRows.clear();
  _coreRows.reserve(orig._coreRows.size());
  for (unsigned i = 0; i < orig._coreRows.size(); ++i)
    _coreRows.push_back(RBPCoreRow(orig._coreRows[i], _placement));

  _backUpCoreRows.clear();
  _backUpCoreRows.reserve(orig._backUpCoreRows.size());
  for (unsigned i = 0; i < orig._backUpCoreRows.size(); ++i)
    _backUpCoreRows.push_back(RBPCoreRow(orig._backUpCoreRows[i], _placement));
}

void RBPlacement::clearStdCellsFromCoreRows(void) {
  for (unsigned i = 0; i < _coreRows.size(); ++i) {
    _coreRows[i].clearMovableStdCells(_isFixed, getHGraph());
  }
}

std::string RBPlacement::getLogicValue(unsigned nodeIdx) const {
  if (_logicValues.find(nodeIdx) == _logicValues.end())
    return "NODE DOES NOT EXIST";
  else
    return _logicValues.find(nodeIdx)->second;
}

const RBPCoreRow& RBPlacement::getBackupCoreRow(unsigned id) const {
  abkassert(id < _backUpCoreRows.size(),
            "Row index out of range in getBackupCoreRow");
  return _backUpCoreRows[id];
}
