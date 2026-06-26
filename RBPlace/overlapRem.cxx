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

// created by Saurabh Adya, Dec 3 2000
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <cmath>
#include <fstream>
#include <algorithm>
#include <list>
#include <set>
#include <unordered_map>

#include "ABKCommon/abkfunc.h"
#include "ABKCommon/abktimer.h"
#ifdef USEFLOORIST
#include "Floorist/Floorist.h"
#endif
#include "RBPlacement.h"

using std::cout;
using std::endl;
using std::list;
using std::make_pair;
using std::map;
using std::max;
using std::min;
using std::ofstream;
using std::pair;
using std::set;
using std::vector;

#ifdef _MSC_VER
#ifndef rint
#define rint(a) floor((a) + 0.5)
#endif
#endif

// added by sadya to remove overlaps in each subrow
void RBPlacement::remOverlapsSR(void) {
  if (_params.verb.getForActions())
    cout << " Running overlap removal ..." << endl;
  double initHPWL = evalHPWL();
  if (_params.verb.getForMajStats())
    cout << "  HPWL before overlap remover is " << initHPWL << endl;
  Timer tm;

  itRBPCoreRow itc;
  itRBPSubRow its;

  unsigned initOverlaps = getNumOverlaps();

  // for each core row, for each subrow scan left to right if overlap found then
  // take action
  for (itc = _coreRows.begin(); itc != _coreRows.end(); ++itc) {
    for (its = itc->subRowsBegin(); its != itc->subRowsEnd(); ++its)
      shuffleSR(its);
  }
  tm.stop();
  double HPWLafter = evalHPWL();
  double HPWLwPins = evalHPWL(true);
  unsigned finalOverlaps = getNumOverlaps();

  if (_params.verb.getForMajStats()) {
    cout << " Removing overlaps took " << tm << endl;
    cout << "  Discovered " << initOverlaps << " overlaps " << endl;
    cout << "  Overlaps remaining: " << finalOverlaps << endl;
    cout << "  HPWL after  overlap removal is " << HPWLafter << endl;
    cout << "  HPWL w pins                 is " << HPWLwPins << endl;
    double change = (HPWLafter == 0) ? 0 : (HPWLafter - initHPWL) / HPWLafter;
    if (change < 0)
      cout << "  % decrease in HPWL is " << floor(-change * 100000 + 0.5) / 1000
           << endl;
    else
      cout << "  % increase in HPWL is " << floor(change * 100000 + 0.5) / 1000
           << endl;
  }
}

// added by sadya to remove overlaps in each subrow
void RBPlacement::remOverlapsVert(void) {
  if (_params.verb.getForActions())
    cout << " Running overlap removal with Vertical Moves..." << endl;
  double initHPWL = evalHPWL();
  if (_params.verb.getForMajStats())
    cout << "  HPWL before overlap remover is " << initHPWL << endl;
  Timer tm;

  itRBPCoreRow itc;
  itRBPSubRow its;
  itRBPCellIds itn;
  itRBPCellIds itntemp;
  itRBPCoreRow itcNext;       // used to scan immediate next row
  itRBPSubRow itsNext;        // used to scan immediate next sub row
  itRBPCellIds itnNext;       // used to scan immediate next sub row cells
  itRBPCellIds itnright;      // used to scan immediate next sub row cells
  itRBPCellIds itnNextLeft;   // used to scan
  itRBPCellIds itnNextRight;  // used to scan

  RBPSubRow tempSR;

  bool leftRight = false;  // false means move left else move right
  double WS = 0;
  unsigned numOverlaps = 0;
  unsigned totalOverlaps = 0;
  double rowHeight = 0;
  double xDiff = 0;
  // for each core row, for each subrow scan left to right if overlap found then
  // take action
  for (itc = _coreRows.begin(); itc != _coreRows.end(); ++itc) {
    rowHeight = itc->getHeight();
    for (its = itc->subRowsBegin(); its != itc->subRowsEnd(); ++its) {
      for (itn = its->cellIdsBegin(); itn != its->cellIdsEnd(); ++itn) {
        // check for overlaps in the same row
        xDiff = 1;
        if (itn != its->cellIdsEnd() - 1) {
          itntemp = itn + 1;
          xDiff = _placement[*itn].x + _hgWDims->getNodeWidth(*itn) -
                  _placement[*itntemp].x;
          if (xDiff > 1e-10) {
            // overlap found, now deal with it
            numOverlaps++;
            totalOverlaps++;

            itnright = itntemp;

            // check for WS in upper rows
            if (itc != _coreRows.end() - 1) {
              itcNext = itc + 1;
              for (itsNext = itcNext->subRowsBegin();
                   itsNext != itcNext->subRowsEnd(); ++itsNext) {
                for (itnNext = itsNext->cellIdsBegin();
                     itnNext != itsNext->cellIdsEnd() - 1; ++itnNext) {
                  // check for WS in the upper adjacent  row
                  if ((_placement[*itnNext].x <= _placement[*itn].x) &&
                      (_placement[*(itnNext + 1)].x >= _placement[*itn].x))
                    break;  // closest cell found in upper
                            // row
                }
                itnNextLeft = itnNext;
                itnNextRight = itnNext;

                while (itnNextRight != itsNext->cellIdsEnd() - 1 ||
                       itnNextLeft != itsNext->cellIdsBegin()) {
                  if (leftRight == false &&
                      itnNextRight != itsNext->cellIdsEnd() - 1)
                    leftRight = true;
                  else if (leftRight == true &&
                           itnNextLeft != itsNext->cellIdsBegin())
                    leftRight = false;

                  if (leftRight == true)  // search right
                  {
                    itntemp = itnNextRight + 1;
                    WS = _placement[*itntemp].x - _placement[*itnNextRight].x -
                         _hgWDims->getNodeWidth(*itnNextRight);
                    if (WS > 0)  // white space found
                    {
                      if (WS >= _hgWDims->getNodeWidth(*itn)) {
                        tempSR = *its;
                        tempSR.removeCell(*itn);
                        _placement[*itn].x = _placement[*itnNextRight].x;
                        _placement[*itn].y = itsNext->getYCoord();
                        tempSR = *itsNext;
                        tempSR.addCell(*itn);
                        --numOverlaps;
                        break;
                      } else if (WS >= _hgWDims->getNodeWidth(*itnright)) {
                        tempSR = *its;
                        tempSR.removeCell(*itnright);
                        _placement[*itnright].x = _placement[*itnNextRight].x;
                        _placement[*itnright].y = itsNext->getYCoord();
                        tempSR = *itsNext;
                        tempSR.addCell(*itnright);
                        --numOverlaps;
                        break;
                      } else {
                      }
                    }
                    ++itnNextRight;
                  } else  // search left for WS
                  {
                    itntemp = itnNextLeft - 1;
                    WS = _placement[*itnNextLeft].x - _placement[*itntemp].x -
                         _hgWDims->getNodeWidth(*itntemp);
                    if (WS > 0)  // white space found
                    {
                      if (WS >= _hgWDims->getNodeWidth(*itn)) {
                        tempSR = *its;
                        tempSR.removeCell(*itn);
                        _placement[*itn].x = _placement[*itntemp].x;
                        _placement[*itn].y = itsNext->getYCoord();
                        tempSR = *itsNext;
                        tempSR.addCell(*itn);
                        --numOverlaps;
                        break;
                      } else if (WS >= _hgWDims->getNodeWidth(*itnright)) {
                        tempSR = *its;
                        tempSR.removeCell(*itnright);
                        _placement[*itnright].x = _placement[*itntemp].x;
                        _placement[*itnright].y = itsNext->getYCoord();
                        tempSR = *itsNext;
                        tempSR.addCell(*itnright);
                        --numOverlaps;
                        break;
                      }
                    }
                    itnNextLeft--;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  tm.stop();
  double HPWLafter = evalHPWL();
  double HPWLwPins = evalHPWL(true);
  if (_params.verb.getForMajStats()) {
    cout << " Removing overlaps took " << tm << endl;
    cout << "  Discovered " << totalOverlaps << " overlaps " << endl;
    cout << "  Overlaps remaining: " << numOverlaps << endl;
    cout << "  HPWL after  overlap removal is " << HPWLafter << endl;
    cout << "  HPWL w pins                 is " << HPWLwPins << endl;
    double change = (HPWLafter == 0) ? 0 : (HPWLafter - initHPWL) / HPWLafter;
    if (change < 0)
      cout << "  % decrease in HPWL is " << floor(-change * 100000 + 0.5) / 1000
           << endl;
    else
      cout << "  % increase in HPWL is " << floor(change * 100000 + 0.5) / 1000
           << endl;
  }
}

void RBPlacement::snapCellsInSites(void) {
  itRBPCoreRow itc;
  itRBPSubRow its;
  itRBPCellIds itn;
  // itRBPCellIds itntemp;
  double siteSpacing;
  int index;

  for (itc = _coreRows.begin(); itc != _coreRows.end(); ++itc) {
    for (its = itc->subRowsBegin(); its != itc->subRowsEnd(); ++its) {
      siteSpacing = its->getSpacing();
      for (itn = its->cellIdsBegin(); itn != its->cellIdsEnd(); ++itn) {
        if (int(_placement[*itn].x) % int(siteSpacing) != 0 ||
            _placement[*itn].x != int(_placement[*itn].x)) {
          // cout<<"cell not in site "<<endl;
          index = *itn;
          int newPosX = int(_placement[*itn].x) -
                        int(_placement[*itn].x) % int(siteSpacing);

          _placement[*itn].x = newPosX;
        }
      }
    }
  }
}

void RBPlacement::remOverlaps(void) {
  if (_params.verb.getForActions())
    cout << " Running overlap removal ..." << endl;
  double initHPWL = evalHPWL();
  if (_params.verb.getForMajStats())
    cout << "  HPWL before overlap remover is " << initHPWL << endl;
  Timer tm;

  itRBPCoreRow itc;
  itRBPCoreRow itcNext;  // used to scan immediate next row
  itRBPCoreRow itcCurr;  // used to scan immediate next row
  itRBPSubRow its;
  itRBPSubRow itsCurr;  // used to track current sub row in bfs
  itRBPCellIds itn;
  itRBPCellIds itntemp;
  itRBPCellIds itnToMove;
  itRBPSubRow itsNext;  // used to scan next sub row
  double rowHeight = 0;
  double nodeHeight = 0;
  double xDiff = 0;
  const unsigned threshold = 50;

  unsigned initOverlaps = getNumOverlaps();

  // for each core row, for each subrow scan left to right if overlap found then
  // take action
  for (itc = _coreRows.begin(); itc != _coreRows.end(); ++itc) {
    rowHeight = itc->getHeight();
    for (its = itc->subRowsBegin(); its != itc->subRowsEnd(); ++its) {
      shuffleSR(its);
      for (itn = its->cellIdsBegin(); itn != its->cellIdsEnd(); ++itn) {
        nodeHeight = _hgWDims->getNodeHeight(*itn);
        // check for overlaps in the same row
        xDiff = 1;
        if (itn != its->cellIdsEnd() - 1) {
          itntemp = itn + 1;
          xDiff = _placement[*itn].x + _hgWDims->getNodeWidth(*itn) -
                  _placement[*itntemp].x;

          if (xDiff > 1e-10) {
            // overlap found, now deal with it

            double srOrigStart = its->getXStart();
            double srOrigEnd =
                srOrigStart + its->getNumSites() * its->getSiteWidth();
            list<itRBPSubRow> bfsList;
            // seenList is a pure visited-set (membership test + insert only,
            // never iterated), so a std::set replaces the former O(n) linear
            // findInList scan with O(log n) lookups.
            std::set<itRBPSubRow> seenList;
            bfsList.push_back(its);
            seenList.insert(its);

            // try to move the smallest node
            if (_hgWDims->getNodeWidth(*itn) >
                _hgWDims->getNodeWidth(*(itn + 1)))
              itnToMove = itn + 1;
            else
              itnToMove = itn;

            unsigned cutoffCtr = 0;
            while (cutoffCtr < threshold && !bfsList.empty()) {
              ++cutoffCtr;
              itsCurr = *(bfsList.begin());
              bfsList.pop_front();
              // get an iterator to the core row for the parent of this subrow
              itcCurr = coreRowsBegin() +
                        (&itsCurr->getCoreRow() - &(*coreRowsBegin()));

              double srStart = itsCurr->getXStart();
              double srEnd =
                  srStart + itsCurr->getNumSites() * itsCurr->getSiteWidth();

              double srWidth = srEnd - srStart;
              double cellWidths = 0;
              for (itntemp = itsCurr->cellIdsBegin();
                   itntemp != itsCurr->cellIdsEnd(); ++itntemp)
                cellWidths += _hgWDims->getNodeWidth(*itntemp);

              double widthToMove = _hgWDims->getNodeWidth(*itnToMove);

              // cout<<" overlap x "<<_placement[*itn].x<<" YCoord:
              // "<<itsCurr->getYCoord()<<" XCoords: "<<srStart<<"  "<<srEnd<<"
              // list Size: "<<bfsList.size()<<" widthToMove "<<widthToMove<<"
              // space "<<(srWidth-cellWidths)<<endl; found sr to move
              if (widthToMove <= (srWidth - cellWidths)) {
                // cout<<"found subrow "<<widthToMove<<"
                // "<<(srWidth-cellWidths)<<"
                // "<<_hgWDims->getNodeByIdx(*itnToMove).getName()<<"
                // "<<*itnToMove<<endl;
                double xLocToMoveTo = 0;

                if (srEnd <= srOrigStart)  // move left most cell
                {
                  double leftNodeWidth =
                      _hgWDims->getNodeWidth(*its->cellIdsBegin());
                  if (leftNodeWidth <= widthToMove) {
                    itnToMove = its->cellIdsBegin();
                    widthToMove = leftNodeWidth;
                  }
                  xLocToMoveTo = srEnd - widthToMove;
                } else if (srStart >= srOrigEnd)  // move right cell
                {
                  double rightNodeWidth =
                      _hgWDims->getNodeWidth(*(its->cellIdsEnd() - 1));
                  if (rightNodeWidth <= widthToMove) {
                    itnToMove = its->cellIdsEnd() - 1;
                    widthToMove = rightNodeWidth;
                  }
                  xLocToMoveTo = srStart;
                } else  // move the current cell
                {
                  if (srEnd <= _placement[*itnToMove].x)
                    xLocToMoveTo = srEnd - widthToMove;
                  else if (srStart >= _placement[*itnToMove].x)
                    xLocToMoveTo = srStart;
                  else
                    xLocToMoveTo = _placement[*itnToMove].x;
                }
                unsigned indexToMove = *itnToMove;

                RBPSubRow &tempSR1 = const_cast<RBPSubRow &>(*its);
                tempSR1.removeCell(indexToMove);
                _placement[indexToMove].x = xLocToMoveTo;
                _placement[indexToMove].y = itsCurr->getYCoord();
                RBPSubRow &tempSR2 = const_cast<RBPSubRow &>(*itsCurr);
                tempSR2.addCell(indexToMove);
                shuffleSR(itsCurr);
                shuffleSR(its);
                break;
              }

              // add new elements to bfsList
              if (itsCurr != itcCurr->subRowsEnd() - 1) {
                itsNext = itsCurr + 1;
                if (seenList.find(itsNext) == seenList.end()) {
                  bfsList.push_back(itsNext);
                  seenList.insert(itsNext);
                }
              }
              if (itsCurr != itcCurr->subRowsBegin()) {
                itsNext = itsCurr - 1;
                if (seenList.find(itsNext) == seenList.end()) {
                  bfsList.push_back(itsNext);
                  seenList.insert(itsNext);
                }
              }
              if (itcCurr != _coreRows.end() - 1) {
                itcNext = itcCurr + 1;
                for (itsNext = itcNext->subRowsBegin();
                     itsNext != itcNext->subRowsEnd(); ++itsNext) {
                  double nextsrStart = itsNext->getXStart();
                  double nextsrEnd = nextsrStart + itsNext->getNumSites() *
                                                       itsNext->getSiteWidth();

                  if ((nextsrStart <= srStart && nextsrEnd >= srStart) ||
                      (nextsrStart >= srStart && nextsrEnd <= srEnd) ||
                      (nextsrStart <= srEnd &&
                       nextsrEnd >= srEnd))  // overlapping SR
                  {
                    if (seenList.find(itsNext) == seenList.end()) {
                      bfsList.push_back(itsNext);
                      seenList.insert(itsNext);
                    }
                  }
                }
              }
              if (itcCurr != _coreRows.begin()) {
                itcNext = itcCurr - 1;
                for (itsNext = itcNext->subRowsBegin();
                     itsNext != itcNext->subRowsEnd(); ++itsNext) {
                  double nextsrStart = itsNext->getXStart();
                  double nextsrEnd = nextsrStart + itsNext->getNumSites() *
                                                       itsNext->getSiteWidth();

                  if ((nextsrStart <= srStart && nextsrEnd >= srStart) ||
                      (nextsrStart >= srStart && nextsrEnd <= srEnd) ||
                      (nextsrStart <= srEnd &&
                       nextsrEnd >= srEnd))  // overlapping SR
                  {
                    if (seenList.find(itsNext) == seenList.end()) {
                      bfsList.push_back(itsNext);
                      seenList.insert(itsNext);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  tm.stop();
  double HPWLafter = evalHPWL();
  double HPWLwPins = evalHPWL(true);
  unsigned finalOverlaps = getNumOverlaps();

  if (1 /*_params.verb.getForMajStats()*/) {
    cout << " Removing overlaps took " << tm << endl;
    cout << "  Discovered " << initOverlaps << " overlaps " << endl;
    cout << "  Overlaps remaining: " << finalOverlaps << endl;
    cout << "  HPWL after  overlap removal is " << HPWLafter << endl;
    cout << "  HPWL w pins                 is " << HPWLwPins << endl;
    double change = (HPWLafter == 0) ? 0 : (HPWLafter - initHPWL) / HPWLafter;
    if (change < 0)
      cout << "  % decrease in HPWL is " << floor(-change * 100000 + 0.5) / 1000
           << endl;
    else
      cout << "  % increase in HPWL is " << floor(change * 100000 + 0.5) / 1000
           << endl;
  }
}

void RBPlacement::shuffleSR(itRBPSubRow its) {
  itRBPCellIds itn;
  itRBPCellIds itntemp;
  itRBPCellIds itnleft;
  itRBPCellIds itnright;

  if (its->getNumCells() == 0) return;

  // itRBPSubRow itsNext;    //used to scan immediate next sub row
  // itRBPCellIds itnNext;   //used to scan immediate next sub row cells

  bool leftRight = false;  // false means move left else move right
  double WS = 0;

  // snap all cells outside row to the row
  bool snapped = false;
  for (itn = its->cellIdsBegin(); itn != its->cellIdsEnd(); ++itn) {
    if (its->getXEnd() < (_placement[*itn].x + _hgWDims->getNodeWidth(*itn))) {
      _placement[*itn].x = its->getXEnd() - _hgWDims->getNodeWidth(*itn);
      snapped = true;
    }
    if (_placement[*itn].x < its->getXStart()) {
      _placement[*itn].x = its->getXStart();
      snapped = true;
    }
  }

  if (snapped) {
    const_cast<RBPSubRow *>(&(*its))->sortCellsByLoc();
    for (unsigned idx = 0; idx < its->getNumCells(); ++idx) {
      _cellCoords[(*its)[idx]].cellOffset = idx;
    }
  }

  if (its->getNumCells() == 1) return;

  double xDiff = 0;
  double totalCellWidths = 0;
  double totalOverlap = 0.;

  for (itn = its->cellIdsBegin(); itn != its->cellIdsEnd(); ++itn) {
    // check for overlaps in the same row
    xDiff = 0.;
    totalCellWidths += _hgWDims->getNodeWidth(*itn);
    itntemp = itn + 1;
    if (itntemp != its->cellIdsEnd()) {
      xDiff = _placement[*itn].x + _hgWDims->getNodeWidth(*itn) -
              _placement[*itntemp].x;

      if (xDiff > 1e-10) {
        // overlap found, now deal with it
        totalOverlap += xDiff;
      }
    }
  }

  if (totalOverlap <= 1e-10) return;

  double totalSiteWidths = its->getXEnd() - its->getXStart();
  itn = its->cellIdsBegin();
  double startWS = _placement[*itn].x - its->getXStart();
  itn = its->cellIdsEnd() - 1;
  double endWS =
      its->getXEnd() - _placement[*itn].x - _hgWDims->getNodeWidth(*itn);
  double WSNeeded = totalCellWidths + startWS + endWS - totalSiteWidths;
  if (WSNeeded > 0) {
    itn = its->cellIdsBegin();
    if (WSNeeded <= startWS) {
      _placement[*itn].x = _placement[*itn].x - WSNeeded;
      WSNeeded = 0;
    } else if (startWS > 0) {
      _placement[*itn].x = _placement[*itn].x - startWS;
      WSNeeded = WSNeeded - startWS;
    }
    if (_placement[*itn].x < its->getXStart())
      _placement[*itn].x = its->getXStart();

    itn = its->cellIdsEnd() - 1;
    if (WSNeeded <= endWS) {
      _placement[*itn].x = _placement[*itn].x + WSNeeded;
      WSNeeded = 0;
    } else if (endWS > 0) {
      _placement[*itn].x = _placement[*itn].x + endWS;
      WSNeeded = WSNeeded - endWS;
    }
    if (its->getXEnd() < (_placement[*itn].x + _hgWDims->getNodeWidth(*itn)))
      _placement[*itn].x = its->getXEnd() - _hgWDims->getNodeWidth(*itn);
  }

  WS = 0;
  xDiff = 0;

  for (itn = its->cellIdsBegin(); itn != its->cellIdsEnd(); ++itn) {
    // check for overlaps in the same row
    xDiff = 0.;
    if (itn != its->cellIdsEnd() - 1) {
      itntemp = itn + 1;
      xDiff = _placement[*itn].x + _hgWDims->getNodeWidth(*itn) -
              _placement[*itntemp].x;

      if (xDiff > 1e-10) {
        // overlap found, now deal with it

        itnleft = itn;
        itnright = itn;
        leftRight = false;

        while (xDiff > 1e-10 && (itnleft != (its->cellIdsBegin()) ||
                                 itnright != (its->cellIdsEnd() - 2))) {
          // choose the direction of search
          if (leftRight == false && itnright != its->cellIdsEnd() - 2)
            leftRight = true;
          else if (leftRight == true && itnleft != its->cellIdsBegin())
            leftRight = false;

          if (leftRight == true)  // search right for WS
          {
            ++itnright;
            itntemp = itnright + 1;
            WS = _placement[*itntemp].x - _placement[*itnright].x -
                 _hgWDims->getNodeWidth(*itnright);
            if (WS > 0)  // WS found
            {
              if (WS > xDiff) WS = xDiff;
              xDiff -= WS;
              for (itntemp = itn + 1; itntemp != itnright + 1; ++itntemp)
                _placement[*itntemp].x += WS;
            }
          } else  // search left for WS
          {
            itntemp = itnleft;
            --itnleft;
            WS = _placement[*itntemp].x - _placement[*itnleft].x -
                 _hgWDims->getNodeWidth(*itnleft);
            if (WS > 0)  // WS found
            {
              if (WS > xDiff) WS = xDiff;
              xDiff -= WS;
              for (itntemp = itn; itntemp != itnleft; --itntemp) {
                _placement[*itntemp].x -= WS;
              }
            }
          }
        }
      }
    }
  }

  // again snap all cells outside row to the row
  snapped = false;
  for (itn = its->cellIdsBegin(); itn != its->cellIdsEnd(); ++itn) {
    if (its->getXEnd() < (_placement[*itn].x + _hgWDims->getNodeWidth(*itn))) {
      _placement[*itn].x = its->getXEnd() - _hgWDims->getNodeWidth(*itn);
      snapped = true;
    }
    if (_placement[*itn].x < its->getXStart()) {
      _placement[*itn].x = its->getXStart();
      snapped = true;
    }
  }

  if (snapped) {
    const_cast<RBPSubRow *>(&(*its))->sortCellsByLoc();
    for (unsigned idx = 0; idx < its->getNumCells(); ++idx) {
      _cellCoords[(*its)[idx]].cellOffset = idx;
    }
  }

  return;
}

unsigned RBPlacement::getNumOverlaps() {
  unsigned overlaps = 0;
  itRBPCellIds itn;
  itRBPCellIds itntemp;
  itRBPCoreRow itc;
  itRBPSubRow its;

  double xDiff = 0;
  for (itc = _coreRows.begin(); itc != _coreRows.end(); ++itc) {
    for (its = itc->subRowsBegin(); its != itc->subRowsEnd(); ++its) {
      for (itn = its->cellIdsBegin(); itn != its->cellIdsEnd(); ++itn) {
        // check for overlaps in the same row
        xDiff = 1;
        if (itn != its->cellIdsEnd() - 1) {
          itntemp = itn + 1;
          xDiff = _placement[*itn].x + _hgWDims->getNodeWidth(*itn) -
                  _placement[*itntemp].x;

          if (xDiff > 1e-10) overlaps++;
        }
      }
    }
  }
  return (overlaps);
}

bool RBPlacement::checkOutOfCoreCells(void) {
  BBox layout = getBBox(false);
  bool consistent = true;
  double rowHeight = coreRowsBegin()->getHeight();

  for (unsigned i = 0; i < getNumCells(); i++) {
    if (_isCoreCell[i] && !_isFixed[i]) {
      double cellWidth = 0;
      double cellHeight = 0;
      Orient nodeOrient = _placement.getOrient(i);
      if (nodeOrient == Orient("N") || nodeOrient == Orient("S") ||
          nodeOrient == Orient("FN") || nodeOrient == Orient("FS")) {
        cellWidth = _hgWDims->getNodeWidth(i);
        cellHeight = _hgWDims->getNodeHeight(i) * 0.9;
      } else {
        cellWidth = _hgWDims->getNodeHeight(i);
        cellHeight = _hgWDims->getNodeWidth(i) * 0.9;
      }

      if (_placement[i].x < layout.xMin ||
          (_placement[i].x + cellWidth) > layout.xMax ||
          _placement[i].y < layout.yMin ||
          (_placement[i].y + cellHeight) > layout.yMax) {
        consistent = false;
        if (cellHeight > rowHeight)
          cout << "Macro (" << i << ") ";
        else
          cout << "Cell (" << i << ") ";
        cout << _hgWDims->getNodeNameByIndex(i) << " (" << _placement[i].x
             << " " << _placement[i].y << "->" << _placement[i].x + cellWidth
             << " " << _placement[i].y + cellHeight << ")" << " out of core "
             << layout << endl;
      }
    }
  }
  return consistent;
}

double RBPlacement::calcOverlapGeneric(overlapData mode, bool print) const {
  vector<unsigned> cellIds(getNumCells());

  for (unsigned i = 0; i < getNumCells(); ++i) {
    cellIds[i] = i;
  }

  const PlacementWOrient &pl = getPlacement();
  BBox bbox = getBBox(true);
  double overlap = calcOverlapSpecific(pl, bbox, cellIds, mode, print);

  double siteArea = getContainedSiteAreaInBBox(bbox);
  double percentOverlap = floor(100000. * (overlap / siteArea) + 0.5) / 1000.;
  if (mode == OnlyFixed) {
    cout << "Total overlap between fixed objects is " << overlap << " OR "
         << percentOverlap << " % of layout area" << endl;
    if (percentOverlap > 0.01)
      cout << "  (most likely this overlap was not created by the placer)"
           << endl;
  } else if (mode == OnlyMacro) {
    cout << "Total overlap between macros is " << overlap << " OR "
         << percentOverlap << " % of layout area" << endl;
    if (percentOverlap > 0 && percentOverlap < 0.01)
      cout << "  (one more run should get rid of this overlap)" << endl;
  } else if (mode == Movable) {
    cout << "Total overlap of all movable objects is " << overlap << " OR "
         << percentOverlap << " % of layout area" << endl;
    if (percentOverlap > 0 && percentOverlap < 0.01)
      cout << "  (one more run should get rid of this overlap)" << endl;
  } else {
    cout << "Total overlap is " << overlap << " OR " << percentOverlap
         << " % of layout area" << endl;
    if (percentOverlap > 0 && percentOverlap < 0.01)
      cout << "  (one more run should get rid of this overlap)" << endl;
  }

  return overlap;
}

double RBPlacement::calcOverlapInBBox(const PlacementWOrient &pl,
                                      const BBox &bbox,
                                      overlapData mode) const {
  vector<unsigned> cellIds;
  for (unsigned i = 0; i < getNumCells(); ++i) {
    unsigned angle = pl.getOrient(i).getAngle();
    bool notRotated = angle == 0 || angle == 180;
    double cellWidth =
        notRotated ? _hgWDims->getNodeWidth(i) : _hgWDims->getNodeHeight(i);
    double cellHeight =
        notRotated ? _hgWDims->getNodeHeight(i) : _hgWDims->getNodeWidth(i);

    BBox cell;
    cell += pl[i];
    cell.add(pl[i].x + cellWidth, pl[i].y + cellHeight);

    BBox intersection = cell.intersect(bbox);
    if (!equalDouble(intersection.getArea(), 0.)) cellIds.push_back(i);
  }

  bool print = false;
  double overlap = calcOverlapSpecific(pl, bbox, cellIds, mode, print);

  return overlap;
}

double RBPlacement::calcOverlapSpecific(const PlacementWOrient &pl,
                                        const BBox &bbox,
                                        const vector<unsigned> &cellIds,
                                        overlapData mode, bool print) const {
  const double P = 1.0;
  unsigned numNodes = cellIds.size();
  unsigned numGridRows = static_cast<unsigned>(sqrt(numNodes / P));
  unsigned numGridCol = numGridRows;

  double gridXSize = bbox.getWidth() / static_cast<double>(numGridCol);
  double gridYSize = bbox.getHeight() / static_cast<double>(numGridRows);

  map<pair<unsigned, unsigned>, BBox> overlaps;

  //  set< BBox > overlaps;

  //  ofstream plot;

  BBox totalBBox = getBBox(false);

  /*
    if(print)
    {
      plot.open("plot.tex");

      BBox coreBBox = getBBox(false);

      for(unsigned i = 0; i < cellIds.size(); ++i)
      {
        unsigned angle = pl.getOrient(i).getAngle();
        bool notRotated = angle == 0 || angle == 180;
        double cellWidth =  notRotated ? _hgWDims->getNodeWidth(i) :
    _hgWDims->getNodeHeight(i); double cellHeight = notRotated ?
    _hgWDims->getNodeHeight(i) : _hgWDims->getNodeWidth(i);

        BBox cell;
        cell += pl[i];
        cell.add(pl[i].x+cellWidth, pl[i].y+cellHeight);

        totalBBox.expandToInclude(cell);
      }

      double unitLen = 10.0/max(totalBBox.getWidth(),totalBBox.getHeight());

      plot << "\\documentclass{article}" << endl
           << "\\usepackage[papersize={" << totalBBox.getWidth()*unitLen <<
    "in," << totalBBox.getHeight()*unitLen
           << "in},nohead,nofoot,margin=0cm]{geometry}" << endl
           << "\\usepackage{graphics}" << endl
           << "\\usepackage{epsfig}" << endl
           << "\\usepackage{color}" << endl
           << "\\usepackage{graphicx,epic,eepicemu}" << endl
           << "\\usepackage{verbatim}" << endl
           << "\\usepackage{nopageno}" << endl
           << "\\begin{document}" << endl
           << "{\\setlength{\\fboxsep}{0in}" << endl
           << "{\\setlength{\\fboxrule}{0.03in}" << endl
           << "\\begin{figure}[t]" << endl
           << "\\setlength{\\unitlength}{" << unitLen << "in}" << endl
           << "\\centering" << endl
           << "\\begin{picture}(" << totalBBox.getWidth() << "," <<
    totalBBox.getHeight() << ")" << endl
           << "\\thicklines" << endl
           << "\\put(-1,-1){\\textcolor{red}{\\drawline[50](" <<
    coreBBox.xMin-totalBBox.xMin << "," << coreBBox.yMin-totalBBox.yMin << ")("
                                                              <<
    coreBBox.xMin-totalBBox.xMin << "," << coreBBox.yMax-totalBBox.yMin << ")}}"
    << endl
           << "\\put(-1,-1){\\textcolor{red}{\\drawline[50](" <<
    coreBBox.xMin-totalBBox.xMin << "," << coreBBox.yMax-totalBBox.yMin << ")("
                                                              <<
    coreBBox.xMax-totalBBox.xMin << "," << coreBBox.yMax-totalBBox.yMin << ")}}"
    << endl
           << "\\put(-1,-1){\\textcolor{red}{\\drawline[50](" <<
    coreBBox.xMax-totalBBox.xMin << "," << coreBBox.yMax-totalBBox.yMin << ")("
                                                              <<
    coreBBox.xMax-totalBBox.xMin << "," << coreBBox.yMin-totalBBox.yMin << ")}}"
    << endl
           << "\\put(-1,-1){\\textcolor{red}{\\drawline[50](" <<
    coreBBox.xMax-totalBBox.xMin << "," << coreBBox.yMin-totalBBox.yMin << ")("
                                                              <<
    coreBBox.xMin-totalBBox.xMin << "," << coreBBox.yMin-totalBBox.yMin << ")}}"
    << endl
           << "\\thinlines" << endl;

      RandomRawDouble rng;

      for(unsigned i = 0; i < cellIds.size(); ++i)
      {
        unsigned angle = pl.getOrient(i).getAngle();
        bool notRotated = angle == 0 || angle == 180;
        double cellWidth =  notRotated ? _hgWDims->getNodeWidth(i) :
    _hgWDims->getNodeHeight(i); double cellHeight = notRotated ?
    _hgWDims->getNodeHeight(i) : _hgWDims->getNodeWidth(i);

        BBox cell;
        cell += pl[i];
        cell.add(pl[i].x+cellWidth, pl[i].y+cellHeight);

        BBox constrained = coreBBox.intersect(cell);

        double grey = rng * 0.15 + 0.75;

        if(constrained != cell)
        {
          plot << "\\put(" << cell.xMin-totalBBox.xMin << "," <<
    cell.yMin-totalBBox.yMin << "){\\colorbox{red}"
               << "{\\makebox(" << cell.getWidth() << "," << cell.getHeight() <<
    "){}}}" << endl;

          plot << "\\put(" << constrained.xMin-totalBBox.xMin << "," <<
    constrained.yMin-totalBBox.yMin << "){\\colorbox[rgb]{"
               << grey << "," << grey << "," << grey << "}{\\makebox(" <<
    constrained.getWidth() << "," << constrained.getHeight() << "){}}}" << endl;
        }
        else
        {
          plot << "\\put(" << constrained.xMin-totalBBox.xMin << "," <<
    constrained.yMin-totalBBox.yMin << "){\\colorbox[rgb]{"
               << grey << "," << grey << "," << grey << "}{\\makebox(" <<
    constrained.getWidth() << "," << constrained.getHeight() << "){}}}" << endl;
        }
      }
    }
  */

  vector<vector<vector<unsigned> > > grid(
      numGridRows, vector<vector<unsigned> >(numGridCol));

  for (unsigned k = 0; k < cellIds.size(); ++k) {
    const unsigned &i = cellIds[k];
    if (_isCoreCell[i] || _hgWDims->isTerminal(i)) {
      if (mode == OnlyFixed) {
        if (!isFixed(i)) continue;
      } else if (mode == OnlyMacro) {
        if (!(isFixed(i) || _hgWDims->isNodeMacro(i))) continue;
      }
      // unfixed terms can overlap with anything
      if (_hgWDims->isTerminal(i) && !isFixed(i)) continue;

      unsigned angle = pl.getOrient(i).getAngle();
      bool notRotated = angle == 0 || angle == 180;
      double cellWidth =
          notRotated ? _hgWDims->getNodeWidth(i) : _hgWDims->getNodeHeight(i);
      double cellHeight =
          notRotated ? _hgWDims->getNodeHeight(i) : _hgWDims->getNodeWidth(i);

      BBox cell;
      cell += pl[i];
      cell.add(pl[i].x + cellWidth, pl[i].y + cellHeight);

      cell = bbox.intersect(cell);
      unsigned startGridCol =
          unsigned(floor((cell.xMin - bbox.xMin) / gridXSize));
      unsigned startGridRow =
          unsigned(floor((cell.yMin - bbox.yMin) / gridYSize));
      unsigned endGridCol = unsigned(ceil((cell.xMax - bbox.xMin) / gridXSize));
      unsigned endGridRow = unsigned(ceil((cell.yMax - bbox.yMin) / gridYSize));
      if (startGridRow > numGridRows) startGridRow = numGridRows;
      if (endGridRow > numGridRows) endGridRow = numGridRows;
      if (startGridCol > numGridCol) startGridCol = numGridCol;
      if (endGridCol > numGridCol) endGridCol = numGridCol;

      for (unsigned rIdx = startGridRow; rIdx < endGridRow; ++rIdx)
        for (unsigned cIdx = startGridCol; cIdx < endGridCol; ++cIdx)
          grid[rIdx][cIdx].push_back(i);
    }
  }

  for (unsigned k = 0; k < cellIds.size(); ++k) {
    const unsigned &i = cellIds[k];
    if (_isCoreCell[i] || _hgWDims->isTerminal(i)) {
      if (mode == OnlyFixed) {
        if (!isFixed(i)) continue;
      } else if (mode == OnlyMacro) {
        if (!_hgWDims->isNodeMacro(i)) continue;
      } else if (mode == Movable) {
        if (isFixed(i)) continue;
      }
      // unfixed terms can overlap with anything
      if (_hgWDims->isTerminal(i) && !isFixed(i)) continue;

      unsigned angle = pl.getOrient(i).getAngle();
      bool notRotated = angle == 0 || angle == 180;
      double cellWidth =
          notRotated ? _hgWDims->getNodeWidth(i) : _hgWDims->getNodeHeight(i);
      double cellHeight =
          notRotated ? _hgWDims->getNodeHeight(i) : _hgWDims->getNodeWidth(i);

      BBox cell;
      cell += pl[i];
      cell.add(pl[i].x + cellWidth, pl[i].y + cellHeight);

      cell = bbox.intersect(cell);
      unsigned startGridCol =
          static_cast<unsigned>(floor((cell.xMin - bbox.xMin) / gridXSize));
      unsigned startGridRow =
          static_cast<unsigned>(floor((cell.yMin - bbox.yMin) / gridYSize));
      unsigned endGridCol =
          static_cast<unsigned>(ceil((cell.xMax - bbox.xMin) / gridXSize));
      unsigned endGridRow =
          static_cast<unsigned>(ceil((cell.yMax - bbox.yMin) / gridYSize));
      if (startGridRow > numGridRows) startGridRow = numGridRows;
      if (endGridRow > numGridRows) endGridRow = numGridRows;
      if (startGridCol > numGridCol) startGridCol = numGridCol;
      if (endGridCol > numGridCol) endGridCol = numGridCol;

      set<unsigned> potentialIntersections;
      for (unsigned rIdx = startGridRow; rIdx < endGridRow; ++rIdx) {
        for (unsigned cIdx = startGridCol; cIdx < endGridCol; ++cIdx) {
          for (unsigned idx = 0; idx < grid[rIdx][cIdx].size(); ++idx) {
            unsigned j = grid[rIdx][cIdx][idx];
            if (i != j) potentialIntersections.insert(j);
          }
        }
      }

      for (set<unsigned>::const_iterator j = potentialIntersections.begin();
           j != potentialIntersections.end(); ++j) {
        unsigned angle2 = pl.getOrient(*j).getAngle();
        bool notRotated2 = angle2 == 0 || angle2 == 180;
        double cellWidth2 = notRotated2 ? _hgWDims->getNodeWidth(*j)
                                        : _hgWDims->getNodeHeight(*j);
        double cellHeight2 = notRotated2 ? _hgWDims->getNodeHeight(*j)
                                         : _hgWDims->getNodeWidth(*j);

        BBox cell2;
        cell2 += pl[*j];
        cell2.add(pl[*j].x + cellWidth2, pl[*j].y + cellHeight2);

        // capture the amount of intersection within the bbox
        BBox intersection = cell.intersect(cell2).intersect(bbox);
        if (greaterThanDouble(intersection.getArea(), 0.)) {
          unsigned smaller = min(i, *j);
          unsigned larger = max(i, *j);
          overlaps[make_pair(smaller, larger)] = intersection;
          //              overlaps.insert(intersection);
          /*
          cout<<_hgWDims->getNodeByIdx(i).getName()<<"
          "<<_hgWDims->getNodeByIdx(*j).getName()
              <<" "<<intersection.getArea()<<endl;
          */
        }
      }
    }
  }

  double overlap = 0.;

  //  set< Point > seenCenters;

  for (map<pair<unsigned, unsigned>, BBox>::const_iterator b = overlaps.begin();
       b != overlaps.end(); ++b)
  //  for(set< BBox >::const_iterator b = overlaps.begin(); b != overlaps.end();
  //  ++b)
  {
    double area = b->second.getArea();
    //    double area = b->getArea();
    overlap += area;
    if (print)
      cout << "Overlap found at " << b->second << " between " << b->first.first
           << " and " << b->first.second << endl;
    //    if(print) cout << "Overlap found at " << *b << endl;
    /*
        if(print)
        {
          double xspan = 300;
          double yspan = 300;

          Point center = b->second.getGeomCenter();

          if(seenCenters.find(center) == seenCenters.end())
          {
            seenCenters.insert(center);

            cout << center.x - 0.5*xspan << " " << center.y + 0.5*yspan << endl
                 << center.x + 0.5*xspan << " " << center.y - 0.5*yspan << endl
       << endl
                 << center.x + 0.5*xspan << " " << center.y + 0.5*yspan << endl
                 << center.x - 0.5*xspan << " " << center.y - 0.5*yspan << endl
       << endl;
          }
        }
    */

    /*
        if(print)
        {

          plot << "\\put(" << b->xMin-totalBBox.xMin << "," <<
       b->yMin-totalBBox.yMin << "){\\colorbox{red}{\\makebox("
               << b->getWidth() << "," << b->getHeight() << "){}}}" << endl;
        }
    */
  }

  /*
    if(print)
    {
      plot << "\\end{picture}" << endl
           << "\\end{figure}" << endl
           << "\\end{document}" << endl;

      plot.close();
    }
  */

  return overlap;
}

// added by MRG to calculate total linear overlap
double RBPlacement::calcOverlap() {
  itRBPCoreRow itc;
  itRBPSubRow its;
  itRBPCellIds itn;
  itRBPCellIds itntemp;

  unsigned numOverlaps = 0;
  unsigned totalOverlaps = 0;
  double xDiff = 0;
  double totalOverlapLength = 0;

  // for each core row, for each subrow scan left to right if overlap found
  // then take action
  for (itc = _coreRows.begin(); itc != _coreRows.end(); ++itc) {
    for (its = itc->subRowsBegin(); its != itc->subRowsEnd(); ++its) {
      for (itn = its->cellIdsBegin(); itn != its->cellIdsEnd(); ++itn) {
        // check for overlaps in the same row
        xDiff = 1;
        if (itn != its->cellIdsEnd() - 1) {
          itntemp = itn + 1;
          if (_placement[*itntemp].x <
              _placement[*itn].x + _hgWDims->getNodeWidth(*itn)) {
            // overlap found, now deal with it
            numOverlaps++;
            totalOverlaps++;

            xDiff = _placement[*itn].x + _hgWDims->getNodeWidth(*itn) -
                    _placement[*itntemp].x;
            totalOverlapLength += (xDiff * itc->getHeight());
          }
        }
      }
    }
  }
  return (totalOverlapLength);
}

// added by MRG to incrementally calculate total linear overlap
double RBPlacement::calcInstOverlap(std::vector<unsigned> &movables) {
  double xDiff = 0;
  double totalOverlapLength = 0;
  itRBPCellIds itn;
  for (unsigned i = 0; i < movables.size(); i++) {
    abkassert(_isInSubRow[movables[i]],
              "Requested cellCoord for cell not in a row");

    RBCellCoord &coord = _cellCoords[movables[i]];
    RBPSubRow &sr = _coreRows[coord.coreRowId][coord.subRowId];
    RBPCoreRow &cr = _coreRows[coord.coreRowId];
    double rowHeight = cr.getHeight();
    double myLeftEdge = _placement[movables[i]].x;
    double myRightEdge =
        _placement[movables[i]].x + _hgWDims->getNodeWidth(movables[i]);
    for (itn = sr.cellIdsBegin(); itn != sr.cellIdsEnd(); ++itn) {
      double theirLeftEdge = _placement[*itn].x;
      double theirRightEdge = _placement[*itn].x + _hgWDims->getNodeWidth(*itn);

      // we're done when we've past the right of our cell
      if (theirLeftEdge > myRightEdge) break;
      if (movables[i] == *itn) continue;

      // i'm contained within another
      if ((myLeftEdge >= theirLeftEdge) && (myRightEdge <= theirRightEdge)) {
        totalOverlapLength += _hgWDims->getNodeWidth(movables[i]) * rowHeight;
      }
      // another is contained within me
      else if ((myLeftEdge <= theirLeftEdge) &&
               (myRightEdge >= theirRightEdge)) {
        totalOverlapLength += _hgWDims->getNodeWidth(*itn) * rowHeight;
      }
      // left side overlap
      else if ((myLeftEdge > theirLeftEdge) && (myLeftEdge < theirRightEdge)) {
        xDiff = theirRightEdge - myLeftEdge;
        totalOverlapLength += (xDiff * rowHeight);
      }
      // right side overlap
      else if ((myRightEdge > theirLeftEdge) &&
               (myRightEdge < theirRightEdge)) {
        xDiff = myRightEdge - theirLeftEdge;
        totalOverlapLength += (xDiff * rowHeight);
      }
    }
  }
  return (totalOverlapLength);
}

bool RBPlacement::findClosestWS(Point &loc, Point &WSLoc, double width) {
  RBRowPtrs rowPtrs = findBothRows(loc);
  RBPSubRow *its = const_cast<RBPSubRow *>(rowPtrs.second);

  itRBPCellIds itn;
  itRBPCellIds itntemp;
  itRBPCellIds itnleft;
  itRBPCellIds itnright;

  bool leftRight = false;  // false means move left else move right
  double WS = 0;
  double xDiff = 0;

  itn = its->cellIdsBegin();
  if (itn == its->cellIdsEnd())  // no cells in subrow
  {
    return (true);
  }
  for (itn = its->cellIdsBegin(); itn != (its->cellIdsEnd() - 1); ++itn) {
    xDiff = 1;
    itntemp = itn + 1;
    if (_placement[*itn].x <= loc.x && _placement[*itntemp].x >= loc.x) break;
  }
  if (itn != (its->cellIdsEnd() - 1))
    xDiff = _placement[*itntemp].x -
            (_placement[*itn].x + _hgWDims->getNodeWidth(*itn));
  else
    xDiff =
        its->getXEnd() - (_placement[*itn].x + _hgWDims->getNodeWidth(*itn));

  if (xDiff >= width) {
    WSLoc.x = _placement[*itn].x + _hgWDims->getNodeWidth(*itn);
    return (true);
  }

  // closest node found, now try to find WS

  itnleft = itn;
  itnright = itn;
  leftRight = false;

  while ((itnleft != (its->cellIdsBegin()) ||
          itnright != (its->cellIdsEnd() - 1))) {
    // chose the direction of search
    if (leftRight == false && itnright != its->cellIdsEnd() - 1)
      leftRight = true;
    else if (leftRight == true && itnleft != its->cellIdsBegin())
      leftRight = false;

    if (leftRight == true)  // search right for WS
    {
      itntemp = itnright + 1;
      WS = _placement[*itntemp].x - _placement[*itnright].x -
           _hgWDims->getNodeWidth(*itnright);
      if (WS >= width)  // WS found
      {
        WSLoc.x = _placement[*itnright].x + _hgWDims->getNodeWidth(*itnright);
        return (true);
      }
      ++itnright;
    } else  // search left for WS
    {
      itntemp = itnleft - 1;
      WS = _placement[*itnleft].x - _placement[*itntemp].x -
           _hgWDims->getNodeWidth(*itntemp);
      if (WS >= width)  // WS found
      {
        WSLoc.x = _placement[*itnleft].x - width;
        return (true);
      }
      --itnleft;
    }
  }
  return (false);
}

void RBPlacement::swapCells(RBPSubRow *sr1, unsigned sr1Cell,
                            const Point &sr1CellPos, RBPSubRow *sr2,
                            unsigned sr2Cell, const Point &sr2CellPos) {
  //  sr1->removeCell(sr1Cell);
  //  _placement[sr1Cell] = sr1CellPos;
  //  sr2->addCell(sr1Cell);
  setLocation(sr1Cell, sr1CellPos);

  if (sr2Cell != UINT_MAX) {
    //    sr2->removeCell(sr2Cell);
    //    _placement[sr2Cell] = sr2CellPos;
    //    sr1->addCell(sr2Cell);
    setLocation(sr2Cell, sr2CellPos);
  }
}

inline double rectPointDist(const Point &p1, const Point &p2) {
  return fabs(p1.x - p2.x) + fabs(p1.y - p2.y);
}

double RBPlacement::CalcSwapCost(const unsigned &cell1,
                                 const Point &cell1NewPos,
                                 const bool skipLargeNets) const {
  bool contains1, endearly;
  BBox temp, beforetemp, aftertemp, c1cur, c1new, pinOffset;
  double diff = 0.;
  itHGFNodeLocal c;
  itHGFEdgeLocal e;
  const Orient &cell1O = getOrient(cell1);
  unsigned netoffset = 0, offset = 0;
  const unsigned maxNetDegree = 32;

  const HGFNode &node1 = _hgWDims->getNodeByIdx(cell1);
  for (e = node1.edgesBegin(); e != node1.edgesEnd(); ++e, ++netoffset) {
    if ((*e)->getDegree() > 1 &&
        (!skipLargeNets || (*e)->getDegree() <= maxNetDegree)) {
      const unsigned &n = (*e)->getIndex();
      endearly = false;
      c1cur = c1new =
          _hgWDims->edgesOnNodePinOffsetBBox(netoffset, cell1, cell1O);
      c1cur.TranslateBy(_placement[cell1]);
      c1new.TranslateBy(cell1NewPos);
      const HGFEdge &edge = _hgWDims->getEdgeByIdx(n);
      temp.clear();
      offset = 0;
      for (c = edge.nodesBegin(); c != edge.nodesEnd(); ++c, ++offset) {
        if ((*c)->getIndex() != cell1) {
          pinOffset = _hgWDims->nodesOnEdgePinOffsetBBox(
              offset, n, getOrient((*c)->getIndex()));
          pinOffset.TranslateBy(_placement[(*c)->getIndex()]);
          temp.expandToInclude(pinOffset);
          contains1 = temp.contains(c1cur) && temp.contains(c1new);
          if (contains1) {
            endearly = true;
            break;
          }
        }
      }
      if (endearly) continue;
      beforetemp = aftertemp = temp;
      beforetemp.expandToInclude(c1cur);
      aftertemp.expandToInclude(c1new);
      diff += (*e)->getWeight() *
              (aftertemp.getHalfPerim() - beforetemp.getHalfPerim());
    }
  }

  return diff;
}

double RBPlacement::CalcSwapCost(const unsigned &cell1, const unsigned &cell2,
                                 const Point &cell1NewPos,
                                 const Point &cell2NewPos,
                                 const bool skipLargeNets) const {
  vector<unsigned> thenets;
  std::unordered_map<unsigned, unsigned, std::hash<unsigned>, std::equal_to<unsigned> > c1offsets,
      c2offsets;
  const unsigned maxNetDegree = 32;

  const HGFNode &node1 = _hgWDims->getNodeByIdx(cell1);
  unsigned offset = 0;
  for (itHGFEdgeLocal e = node1.edgesBegin(); e != node1.edgesEnd();
       ++e, ++offset) {
    if ((*e)->getDegree() > 1 &&
        (!skipLargeNets || (*e)->getDegree() <= maxNetDegree)) {
      thenets.push_back((*e)->getIndex());
      c1offsets[(*e)->getIndex()] = offset;
    }
  }
  offset = 0;
  const HGFNode &node2 = _hgWDims->getNodeByIdx(cell2);
  for (itHGFEdgeLocal e = node2.edgesBegin(); e != node2.edgesEnd();
       ++e, ++offset) {
    if ((*e)->getDegree() > 1 &&
        (!skipLargeNets || (*e)->getDegree() <= maxNetDegree)) {
      thenets.push_back((*e)->getIndex());
      c2offsets[(*e)->getIndex()] = offset;
    }
  }

  sort(thenets.begin(), thenets.end());
  thenets.erase(unique(thenets.begin(), thenets.end()), thenets.end());

  bool has1, has2, contains1, contains2, endearly;
  BBox temp, beforetemp, aftertemp;
  BBox c1cur, c1new, c2cur, c2new, pinOffset;
  double diff = 0.;
  itHGFNodeLocal c;
  const Orient &cell1O = getOrient(cell1);
  const Orient &cell2O = getOrient(cell2);

  for (vector<unsigned>::iterator n = thenets.begin(); n != thenets.end();
       ++n) {
    has1 = has2 = endearly = false;
    if (c1offsets.find(*n) != c1offsets.end()) {
      has1 = true;
      c1cur = c1new =
          _hgWDims->edgesOnNodePinOffsetBBox(c1offsets[*n], cell1, cell1O);
      c1cur.TranslateBy(_placement[cell1]);
      c1new.TranslateBy(cell1NewPos);
    }
    if (c2offsets.find(*n) != c2offsets.end()) {
      has2 = true;
      c2cur = c2new =
          _hgWDims->edgesOnNodePinOffsetBBox(c2offsets[*n], cell2, cell2O);
      c2cur.TranslateBy(_placement[cell2]);
      c2new.TranslateBy(cell2NewPos);
    }
    const HGFEdge &edge = _hgWDims->getEdgeByIdx(*n);
    temp.clear();
    contains1 = !has1;
    contains2 = !has2;
    offset = 0;
    for (c = edge.nodesBegin(); c != edge.nodesEnd(); ++c, ++offset) {
      if ((*c)->getIndex() != cell1 && (*c)->getIndex() != cell2) {
        pinOffset = _hgWDims->nodesOnEdgePinOffsetBBox(
            offset, *n, getOrient((*c)->getIndex()));
        pinOffset.TranslateBy(_placement[(*c)->getIndex()]);
        temp.expandToInclude(pinOffset);
        contains1 = contains1 || (temp.contains(c1cur) && temp.contains(c1new));
        contains2 = contains2 || (temp.contains(c2cur) && temp.contains(c2new));
        if (contains1 && contains2) {
          endearly = true;
          break;
        }
      }
    }
    if (endearly) continue;
    beforetemp = aftertemp = temp;
    if (has1) {
      beforetemp.expandToInclude(c1cur);
      aftertemp.expandToInclude(c1new);
    }
    if (has2) {
      beforetemp.expandToInclude(c2cur);
      aftertemp.expandToInclude(c2new);
    }
    diff += edge.getWeight() *
            (aftertemp.getHalfPerim() - beforetemp.getHalfPerim());
  }

  return diff;
}

double RBPlacement::CalcSwapCost(const unsigned &cell1, const unsigned &cell2,
                                 const unsigned &cell3,
                                 const Point &cell1NewPos,
                                 const Point &cell2NewPos,
                                 const Point &cell3NewPos,
                                 const bool skipLargeNets) const {
  vector<unsigned> thenets;
  std::unordered_map<unsigned, unsigned, std::hash<unsigned>, std::equal_to<unsigned> > c1offsets,
      c2offsets, c3offsets;
  const unsigned maxNetDegree = 32;

  const HGFNode &node1 = _hgWDims->getNodeByIdx(cell1);
  unsigned offset = 0;
  for (itHGFEdgeLocal e = node1.edgesBegin(); e != node1.edgesEnd();
       ++e, ++offset) {
    if ((*e)->getDegree() > 1 &&
        (!skipLargeNets || (*e)->getDegree() <= maxNetDegree)) {
      thenets.push_back((*e)->getIndex());
      c1offsets[(*e)->getIndex()] = offset;
    }
  }
  offset = 0;
  const HGFNode &node2 = _hgWDims->getNodeByIdx(cell2);
  for (itHGFEdgeLocal e = node2.edgesBegin(); e != node2.edgesEnd();
       ++e, ++offset) {
    if ((*e)->getDegree() > 1 &&
        (!skipLargeNets || (*e)->getDegree() <= maxNetDegree)) {
      thenets.push_back((*e)->getIndex());
      c2offsets[(*e)->getIndex()] = offset;
    }
  }
  offset = 0;
  const HGFNode &node3 = _hgWDims->getNodeByIdx(cell3);
  for (itHGFEdgeLocal e = node3.edgesBegin(); e != node3.edgesEnd();
       ++e, ++offset) {
    if ((*e)->getDegree() > 1 &&
        (!skipLargeNets || (*e)->getDegree() <= maxNetDegree)) {
      thenets.push_back((*e)->getIndex());
      c3offsets[(*e)->getIndex()] = offset;
    }
  }

  sort(thenets.begin(), thenets.end());
  thenets.erase(unique(thenets.begin(), thenets.end()), thenets.end());

  bool has1, has2, has3, contains1, contains2, contains3, endearly;
  BBox temp, beforetemp, aftertemp;
  BBox c1cur, c1new, c2cur, c2new, c3cur, c3new, pinOffset;
  double diff = 0.;
  itHGFNodeLocal c;
  const Orient &cell1O = getOrient(cell1);
  const Orient &cell2O = getOrient(cell2);
  const Orient &cell3O = getOrient(cell3);

  for (vector<unsigned>::iterator n = thenets.begin(); n != thenets.end();
       ++n) {
    has1 = has2 = has3 = endearly = false;
    if (c1offsets.find(*n) != c1offsets.end()) {
      has1 = true;
      c1cur = c1new =
          _hgWDims->edgesOnNodePinOffsetBBox(c1offsets[*n], cell1, cell1O);
      c1cur.TranslateBy(_placement[cell1]);
      c1new.TranslateBy(cell1NewPos);
    }
    if (c2offsets.find(*n) != c2offsets.end()) {
      has2 = true;
      c2cur = c2new =
          _hgWDims->edgesOnNodePinOffsetBBox(c2offsets[*n], cell2, cell2O);
      c2cur.TranslateBy(_placement[cell2]);
      c2new.TranslateBy(cell2NewPos);
    }
    if (c3offsets.find(*n) != c3offsets.end()) {
      has3 = true;
      c3cur = c3new =
          _hgWDims->edgesOnNodePinOffsetBBox(c3offsets[*n], cell3, cell3O);
      c3cur.TranslateBy(_placement[cell3]);
      c3new.TranslateBy(cell3NewPos);
    }
    const HGFEdge &edge = _hgWDims->getEdgeByIdx(*n);
    temp.clear();
    contains1 = !has1;
    contains2 = !has2;
    contains3 = !has3;
    offset = 0;
    for (c = edge.nodesBegin(); c != edge.nodesEnd(); ++c, ++offset) {
      if ((*c)->getIndex() != cell1 && (*c)->getIndex() != cell2 &&
          (*c)->getIndex() != cell3) {
        pinOffset = _hgWDims->nodesOnEdgePinOffsetBBox(
            offset, *n, getOrient((*c)->getIndex()));
        pinOffset.TranslateBy(_placement[(*c)->getIndex()]);
        temp.expandToInclude(pinOffset);
        contains1 = contains1 || (temp.contains(c1cur) && temp.contains(c1new));
        contains2 = contains2 || (temp.contains(c2cur) && temp.contains(c2new));
        contains3 = contains3 || (temp.contains(c3cur) && temp.contains(c3new));
        if (contains1 && contains2 && contains3) {
          endearly = true;
          break;
        }
      }
    }
    if (endearly) continue;
    beforetemp = aftertemp = temp;
    if (has1) {
      beforetemp.expandToInclude(c1cur);
      aftertemp.expandToInclude(c1new);
    }
    if (has2) {
      beforetemp.expandToInclude(c2cur);
      aftertemp.expandToInclude(c2new);
    }
    if (has3) {
      beforetemp.expandToInclude(c3cur);
      aftertemp.expandToInclude(c3new);
    }
    diff += edge.getWeight() *
            (aftertemp.getHalfPerim() - beforetemp.getHalfPerim());
  }

  return diff;
}

double RBPlacement::calcSwapCost(const unsigned &cell1, unsigned cell2,
                                 unsigned cell3, const Point &cell1NewPos,
                                 const Point &cell2NewPos,
                                 const Point &cell3NewPos,
                                 const bool skipLargeNets) const {
  if (cell2 == UINT_MAX)
    return CalcSwapCost(cell1, cell1NewPos, skipLargeNets);
  else if (cell3 == UINT_MAX)
    return CalcSwapCost(cell1, cell2, cell1NewPos, cell2NewPos, skipLargeNets);
  else
    return CalcSwapCost(cell1, cell2, cell3, cell1NewPos, cell2NewPos,
                        cell3NewPos, skipLargeNets);
}

double RBPlacement::calcSwapCost(const itRBPSubRow &sr1, unsigned sr1Cell,
                                 Point &sr1CellPos, const itRBPSubRow &sr2,
                                 unsigned sr2Cell, Point &sr2CellPos) const {
  sr1CellPos.x = _placement[sr1Cell].x;
  sr1CellPos.y = sr2->getYCoord();
  if (sr1CellPos.x > (sr2->getXEnd() - _hgWDims->getNodeWidth(sr1Cell)))
    sr1CellPos.x = sr2->getXEnd() - _hgWDims->getNodeWidth(sr1Cell);
  if (sr1CellPos.x < sr2->getXStart()) sr1CellPos.x = sr2->getXStart();

  if (sr2Cell != UINT_MAX) {
    sr2CellPos.x = _placement[sr2Cell].x;
    sr2CellPos.y = sr1->getYCoord();
    if (sr2CellPos.x > (sr1->getXEnd() - _hgWDims->getNodeWidth(sr2Cell)))
      sr2CellPos.x = sr1->getXEnd() - _hgWDims->getNodeWidth(sr2Cell);
    if (sr2CellPos.x < sr1->getXStart()) sr2CellPos.x = sr1->getXStart();
  }

  if (sr2Cell == UINT_MAX)
    return CalcSwapCost(sr1Cell, sr1CellPos, true);
  else
    return CalcSwapCost(sr1Cell, sr2Cell, sr1CellPos, sr2CellPos, true);
}

double distToSR(const itRBPSubRow &its, const Point &pt) {
  if (pt.x > its->getXEnd()) {
    return (pt.x - its->getXEnd()) + fabs(pt.y - its->getYCoord());
  } else if (pt.x < its->getXStart()) {
    return (its->getXStart() - pt.x) + fabs(pt.y - its->getYCoord());
  } else {
    return fabs(pt.y - its->getYCoord());
  }
}

unsigned RBPlacement::findBestSwap(
    const vector<unsigned> &ordering,
    const vector<pair<itRBPSubRow, double> > &haveSpace,
    const itRBPSubRow &sourceRow, unsigned sourceCell, Point &newSourcePos,
    unsigned &swapCell, Point &newSwapPos, double &bestImprovement,
    double &bestCost, const double &distCutOff) {
  const double lenNeeded = _hgWDims->getNodeWidth(sourceCell);
  double lenTarget, cost, improvement, dist;
  Point sourcePos, swapPos;
  unsigned bestSR = UINT_MAX;
  itRBPCellIds c;
  bestImprovement = -DBL_MAX;
  bestCost = DBL_MAX;

  // examine in the order of decreasing available space
  for (unsigned i = 0; i < ordering.size(); ++i) {
    const unsigned &rIdx = ordering[i];
    if (haveSpace[rIdx].second < 1.e-6)
      break;  // since available space will only get smaller
    if (distCutOff != DBL_MAX) {
      dist = distToSR(haveSpace[rIdx].first, _placement[sourceCell]);
      if (dist > distCutOff) continue;
    }
    if (lenNeeded > haveSpace[rIdx].second) {
      if (haveSpace[rIdx].second < bestImprovement)
        break;  // since available space will only get smaller
      lenTarget = lenNeeded - haveSpace[rIdx].second;
      if (lenTarget >= lenNeeded) continue;
      for (c = haveSpace[rIdx].first->cellIdsBegin();
           c != haveSpace[rIdx].first->cellIdsEnd(); ++c) {
        const double &curLen = _hgWDims->getNodeWidth(*c);
        if (curLen < lenTarget || curLen >= lenNeeded) continue;
        improvement = lenNeeded - curLen;
        if (improvement < bestImprovement) continue;
        cost = calcSwapCost(sourceRow, sourceCell, sourcePos,
                            haveSpace[rIdx].first, *c, swapPos);
        if (improvement == bestImprovement && cost >= bestCost) continue;
        bestImprovement = improvement;
        bestCost = cost;
        bestSR = rIdx;
        swapCell = *c;
        newSourcePos = sourcePos;
        newSwapPos = swapPos;
      }
    } else {
      improvement = lenNeeded;
      cost = calcSwapCost(sourceRow, sourceCell, sourcePos,
                          haveSpace[rIdx].first, UINT_MAX, swapPos);
      if (improvement == bestImprovement && cost >= bestCost) continue;
      bestImprovement = improvement;
      bestCost = cost;
      bestSR = rIdx;
      swapCell = UINT_MAX;
      newSourcePos = sourcePos;
      newSwapPos = swapPos;
    }
  }
  return bestSR;
}

void RBPlacement::royjRemOverlaps(MaxMem *maxMem) {
  if (_params.verb.getForActions())
    cout << " Running overlap removal ..." << endl;

  Timer tm;

  double startWL = 0.;
  if (_params.verb.getForMajStats()) {
    startWL = evalHPWL(true);
    cout << "  HPWL before overlap remover is " << startWL << endl;
  }

  itRBPCoreRow itc;
  itRBPSubRow its;

  // find any cells not in a row and put them in the closest one
  double dist, bestDist;
  Point newPos, bestNewPos;
  unsigned oormoves = 0;
  for (unsigned i = 0; i < getNumCells(); i++) {
    if (_isCoreCell[i] && !_isFixed[i] && !_hgWDims->isNodeMacro(i)) {
      if (_isInSubRow[i]) continue;

      bestDist = DBL_MAX;
      for (itc = _coreRows.begin(); itc != _coreRows.end(); ++itc) {
        for (its = itc->subRowsBegin(); its != itc->subRowsEnd(); ++its) {
          newPos.x = _placement[i].x;
          newPos.y = its->getYCoord();
          if (newPos.x > (its->getXEnd() - _hgWDims->getNodeWidth(i)))
            newPos.x = its->getXEnd() - _hgWDims->getNodeWidth(i);
          if (newPos.x < its->getXStart()) newPos.x = its->getXStart();

          dist = rectPointDist(_placement[i], newPos);
          if (dist >= bestDist) continue;
          bestDist = dist;
          bestNewPos = newPos;
        }
      }

      if (bestDist == DBL_MAX)
        continue;  // this should not happen, but necessary for safety

      setLocation(i, bestNewPos);
      ++oormoves;
    }
  }

  if (_params.verb.getForMajStats() > 1) {
    cout << "  Number of out of row moves " << oormoves << endl;
    cout << "  HPWL after moving cells to rows is " << evalHPWL(true) << endl;
  }

  unsigned initOverlaps = 0, moves = 0;
  if (_params.verb.getForMajStats()) initOverlaps = getNumOverlaps();

  if (_params.verb.getForMajStats() > 1) {
    cout << "  Number of row overlaps " << initOverlaps << endl;
  }

  vector<pair<itRBPSubRow, double> > haveSpace, needSpace;
  vector<pair<itRBPSubRow, double> >::iterator r;
  vector<unsigned> ordering(0);
  double totalCellLen, surplus;
  itRBPCellIds c;
  double sitewidth =
      (_coreRows.size() > 0) ? _coreRows[0].getSiteWidth() : 1.e-6;
  for (itc = _coreRows.begin(); itc != _coreRows.end(); ++itc) {
    for (its = itc->subRowsBegin(); its != itc->subRowsEnd(); ++its) {
      totalCellLen = 0.;
      for (c = its->cellIdsBegin(); c != its->cellIdsEnd(); ++c) {
        totalCellLen += _hgWDims->getNodeWidth(*c);
      }
      surplus = its->getLength() - totalCellLen;
      if (surplus >= 0.) {
        shuffleSR(its);
        if (surplus >= sitewidth) {
          haveSpace.push_back(pair<itRBPSubRow, double>(its, surplus));
          ordering.push_back(ordering.size());
        }
      } else if (surplus < -1.e-6)
        needSpace.push_back(pair<itRBPSubRow, double>(its, -surplus));
    }
  }

  resetPlacement();

  unsigned initoverfull = needSpace.size();
  unsigned curoverfull = initoverfull;

  if (_params.verb.getForMajStats() > 1) {
    cout << "  Number of overfull rows " << needSpace.size() << endl;
    cout << "  Number of underfull rows " << haveSpace.size() << endl;
  }

  sort(ordering.begin(), ordering.end(), OrderBySpace(haveSpace));

  vector<vector<cellMoveData> > themoves(needSpace.size());
  vector<unsigned> rowbestmove(needSpace.size());
  unsigned overallbestmove = UINT_MAX, rIdx, SR, swapCell;
  double cost, bestCost, rowBestCost, improvement, bestImprovement,
      rowBestImprovement;
  Point newSwapPos;
  const unsigned numNodes = getNumCells() - _hgWDims->getNumTerminals();
  double distCutOff = DBL_MAX;
  if (numNodes >= _params.overlapCells) {
    distCutOff =
        sqrt(static_cast<double>(_params.overlapCells) *
             getBBox(false).getArea() / static_cast<double>(2 * numNodes));
  }

  // figure out all the possible moves, as well as the best per row and best
  // overall
  bestImprovement = -DBL_MAX;
  bestCost = DBL_MAX;
  unsigned possibleswaps = 0;
  for (rIdx = 0, r = needSpace.begin(); r != needSpace.end(); ++rIdx, ++r) {
    rowBestImprovement = -DBL_MAX;
    rowBestCost = DBL_MAX;
    rowbestmove[rIdx] = UINT_MAX;
    for (c = r->first->cellIdsBegin(); c != r->first->cellIdsEnd(); ++c) {
      SR = findBestSwap(ordering, haveSpace, r->first, *c, newPos, swapCell,
                        newSwapPos, improvement, cost, distCutOff);
      if (SR == UINT_MAX) continue;
      abkassert(SR < haveSpace.size(), "bad dest row");
      themoves[rIdx].push_back(cellMoveData(rIdx, *c, newPos, SR, swapCell,
                                            newSwapPos, improvement, cost));
      ++possibleswaps;
      if (improvement < rowBestImprovement) continue;
      if (improvement == rowBestImprovement && cost >= rowBestCost) continue;
      rowbestmove[rIdx] = themoves[rIdx].size() - 1;
      rowBestImprovement = improvement;
      rowBestCost = cost;
    }
    if (_params.verb.getForMajStats() > 3) {
      cout << "Finished finding swaps for overfull row " << rIdx + 1 << endl;
      cout << possibleswaps << " possible swaps so far" << endl;
    }
    if (rowbestmove[rIdx] == UINT_MAX) continue;
    if (themoves[rIdx][rowbestmove[rIdx]].WLcost > bestCost) continue;
    if (themoves[rIdx][rowbestmove[rIdx]].WLcost == bestCost &&
        themoves[rIdx][rowbestmove[rIdx]].improvement <= bestImprovement)
      continue;
    overallbestmove = rIdx;
    bestImprovement = themoves[rIdx][rowbestmove[rIdx]].improvement;
    bestCost = themoves[rIdx][rowbestmove[rIdx]].WLcost;
  }

  if (_params.verb.getForMajStats() > 1) {
    cout << "  Number of possible swaps " << possibleswaps << endl;
    tm.split();
    cout << "  Determining possible swaps took " << tm << endl;
    tm.resume();
  }

  maxMem->update("Standard cell overlap remover");

  cellMoveData bestmove;
  while (overallbestmove != UINT_MAX) {
    bestmove.sourceSRIdx =
        themoves[overallbestmove][rowbestmove[overallbestmove]].sourceSRIdx;
    bestmove.sourceCell =
        themoves[overallbestmove][rowbestmove[overallbestmove]].sourceCell;
    bestmove.sourceNewPos =
        themoves[overallbestmove][rowbestmove[overallbestmove]].sourceNewPos;
    bestmove.destSRIdx =
        themoves[overallbestmove][rowbestmove[overallbestmove]].destSRIdx;
    bestmove.destCell =
        themoves[overallbestmove][rowbestmove[overallbestmove]].destCell;
    bestmove.destNewPos =
        themoves[overallbestmove][rowbestmove[overallbestmove]].destNewPos;
    bestmove.improvement =
        themoves[overallbestmove][rowbestmove[overallbestmove]].improvement;
    bestmove.WLcost =
        themoves[overallbestmove][rowbestmove[overallbestmove]].WLcost;

    // make the swap and record it
    swapCells(
        const_cast<RBPSubRow *>(&(*needSpace[bestmove.sourceSRIdx].first)),
        bestmove.sourceCell, bestmove.sourceNewPos,
        const_cast<RBPSubRow *>(&(*haveSpace[bestmove.destSRIdx].first)),
        bestmove.destCell, bestmove.destNewPos);
    ++moves;

    // update the amount of space everything has now
    needSpace[bestmove.sourceSRIdx].second -=
        _hgWDims->getNodeWidth(bestmove.sourceCell);
    haveSpace[bestmove.destSRIdx].second -=
        _hgWDims->getNodeWidth(bestmove.sourceCell);
    if (bestmove.destCell != UINT_MAX) {
      needSpace[bestmove.sourceSRIdx].second +=
          _hgWDims->getNodeWidth(bestmove.destCell);
      haveSpace[bestmove.destSRIdx].second +=
          _hgWDims->getNodeWidth(bestmove.destCell);
    }
    sort(ordering.begin(), ordering.end(), OrderBySpace(haveSpace));

    unsigned updates = 0;
    if (needSpace[bestmove.sourceSRIdx].second < 1.e-6) {
      --curoverfull;
      // if the row has space now, put it on the haveSpace list
      if (needSpace[bestmove.sourceSRIdx].second <= -sitewidth) {
        haveSpace.push_back(
            pair<itRBPSubRow, double>(needSpace[bestmove.sourceSRIdx].first,
                                      -needSpace[bestmove.sourceSRIdx].second));
        ordering.push_back(ordering.size());
        sort(ordering.begin(), ordering.end(), OrderBySpace(haveSpace));
      }
      // get rid of all the moves for this row since its fixed
      themoves[overallbestmove].clear();
      // shuffle this row for better HPWL estimates
      shuffleSR(needSpace[bestmove.sourceSRIdx].first);
    } else {
      // erase the move we just made
      themoves[overallbestmove].erase(themoves[overallbestmove].begin() +
                                      rowbestmove[overallbestmove]);
      // add a new possible move for the cell we just added to this row, if we
      // can
      if (bestmove.destCell != UINT_MAX) {
        ++updates;
        SR = findBestSwap(ordering, haveSpace,
                          needSpace[bestmove.sourceSRIdx].first,
                          bestmove.destCell, newPos, swapCell, newSwapPos,
                          improvement, cost, distCutOff);
        if (SR != UINT_MAX)
          themoves[overallbestmove].push_back(
              cellMoveData(bestmove.sourceSRIdx, bestmove.destCell, newPos, SR,
                           swapCell, newSwapPos, improvement, cost));
      }
    }
    // signify that this rows best move need to be updated
    rowbestmove[overallbestmove] = UINT_MAX;
    // pick up a new best move and update other moves that may not be valid
    overallbestmove = UINT_MAX;
    bestImprovement = -DBL_MAX;
    bestCost = DBL_MAX;
    for (unsigned i = 0; i < themoves.size(); ++i) {
      // check all the moves and see if they should be updated
      bool madeachange = false;
      for (unsigned j = 0; j < themoves[i].size(); ++j) {
        if (themoves[i][j].destSRIdx == bestmove.destSRIdx) {
          // do a quick check first
          double space = haveSpace[bestmove.destSRIdx].second -
                         _hgWDims->getNodeWidth(themoves[i][j].sourceCell);
          if (themoves[i][j].destCell != UINT_MAX)
            space += _hgWDims->getNodeWidth(themoves[i][j].destCell);
          if (space >= 0.) continue;
          // curses! we need to update
          madeachange = true;
          ++updates;
          SR = findBestSwap(ordering, haveSpace,
                            needSpace[themoves[i][j].sourceSRIdx].first,
                            themoves[i][j].sourceCell, newPos, swapCell,
                            newSwapPos, improvement, cost, distCutOff);
          if (SR == UINT_MAX) {
            themoves[i].erase(themoves[i].begin() + j);
            --j;
          } else {
            themoves[i][j].sourceNewPos = newPos;
            themoves[i][j].destSRIdx = SR;
            themoves[i][j].destCell = swapCell;
            themoves[i][j].destNewPos = newSwapPos;
            themoves[i][j].improvement = improvement;
            themoves[i][j].WLcost = cost;
          }
        }
      }

      if (madeachange || rowbestmove[i] == UINT_MAX) {
        rowBestImprovement = -DBL_MAX;
        rowBestCost = DBL_MAX;
        rowbestmove[i] = UINT_MAX;
        for (unsigned j = 0; j < themoves[i].size(); ++j) {
          if (themoves[i][j].improvement < rowBestImprovement) continue;
          if (themoves[i][j].improvement == rowBestImprovement &&
              themoves[i][j].WLcost >= rowBestCost)
            continue;
          rowbestmove[i] = j;
          rowBestImprovement = improvement;
          rowBestCost = cost;
        }
      }

      if (rowbestmove[i] == UINT_MAX) continue;
      if (themoves[i][rowbestmove[i]].WLcost > bestCost) continue;
      if (themoves[i][rowbestmove[i]].WLcost == bestCost &&
          themoves[i][rowbestmove[i]].improvement <= bestImprovement)
        continue;
      overallbestmove = i;
      bestImprovement = themoves[i][rowbestmove[i]].improvement;
      bestCost = themoves[i][rowbestmove[i]].WLcost;
    }
    if (_params.verb.getForMajStats() > 3) {
      cout << "Just made move " << moves << endl;
      cout << "  Space still needed by source row: "
           << needSpace[bestmove.sourceSRIdx].second << endl;
      cout << "  Overfull rows left: " << curoverfull << endl;
      cout << "  Updates needed: " << updates << endl;
    }
  }

  if (_params.verb.getForMajStats() > 1)
    cout << "  HPWL after cell movement is " << evalHPWL(true) << endl;

  // remove the overlaps from all the subrows, since they should have room now
  for (itc = _coreRows.begin(); itc != _coreRows.end(); ++itc)
    for (its = itc->subRowsBegin(); its != itc->subRowsEnd(); ++its)
      shuffleSR(its);

  resetPlacement();

  double endWL = 0.;
  unsigned overfullrows = 0, finalOverlaps = 0;
  if (_params.verb.getForMajStats()) {
    endWL = evalHPWL(true);
    finalOverlaps = getNumOverlaps();
    for (r = needSpace.begin(); r != needSpace.end(); ++r)
      if (r->second > 1.e-6) ++overfullrows;
  }

  if (_params.verb.getForMajStats() > 1)
    cout << "  HPWL after shuffling is " << endWL << endl;

  tm.stop();

  if (_params.verb.getForMajStats()) {
    cout << " Removing overlaps took " << tm << endl;
    if (oormoves > 0)
      cout << "  Out of row moves performed: " << oormoves << endl;
    if (initOverlaps > 0 || finalOverlaps > 0) {
      cout << "  Discovered " << initOverlaps << " overlaps " << endl;
      cout << "  Overlaps remaining: " << finalOverlaps << endl;
    }
    if (initoverfull > 0 || overfullrows > 0) {
      cout << "  Discovered " << initoverfull << " overfull subrows " << endl;
      cout << "  Overfull subrows remaining: " << overfullrows << endl;
      cout << "  Overfull subrow moves performed: " << moves << endl;
    }
    if (endWL != startWL) {
      cout << "  HPWL after overlap removal is " << endWL << endl;
      double change = (startWL == 0) ? 0 : (endWL - startWL) / startWL;
      if (change < 0)
        cout << "  % decrease in HPWL is "
             << floor(-change * 100000 + 0.5) / 1000 << endl;
      else
        cout << "  % increase in HPWL is "
             << floor(change * 100000 + 0.5) / 1000 << endl;
    }
  }
}

bool RBPlacement::findBestMacroMove(
    const PlacementWOrient &pl, const vector<vector<vector<unsigned> > > &grid,
    const double &gridXSize, const double &gridYSize,
    const unsigned &numGridRows, const unsigned &numGridCol,
    const unsigned &macro, const vector<macroData> &macros, Point &trans,
    double &improve, bool snapMacrosX, bool snapMacrosY) {
  BBox layout = getBBox(true), core = getBBox(false), cell1, cell2,
       intersection;

  unsigned startGridCol, startGridRow, endGridCol, endGridRow;
  unsigned rIdx, cIdx, idx;
  double area;

  double siteSpacing = coreRowsBegin()->getSpacing();
  double rowSpacing =
      (coreRowsBegin() + 1)->getYCoord() - coreRowsBegin()->getYCoord();

  set<unsigned> potential_overlaps;

  cell1 += pl[macros[macro].id];
  cell1.add(pl[macros[macro].id].x + macros[macro].width,
            pl[macros[macro].id].y + macros[macro].height);

  startGridCol = unsigned(floor((cell1.xMin - layout.xMin) / gridXSize));
  startGridRow = unsigned(floor((cell1.yMin - layout.yMin) / gridYSize));
  endGridCol = unsigned(ceil((cell1.xMax - layout.xMin) / gridXSize));
  endGridRow = unsigned(ceil((cell1.yMax - layout.yMin) / gridYSize));
  if (startGridRow > numGridRows) startGridRow = numGridRows;
  if (endGridRow > numGridRows) endGridRow = numGridRows;
  if (startGridCol > numGridCol) startGridCol = numGridCol;
  if (endGridCol > numGridCol) endGridCol = numGridCol;

  for (rIdx = startGridRow; rIdx < endGridRow; ++rIdx)
    for (cIdx = startGridCol; cIdx < endGridCol; ++cIdx)
      for (idx = 0; idx < grid[rIdx][cIdx].size(); ++idx) {
        const unsigned &potential = grid[rIdx][cIdx][idx];
        if (potential != macro) potential_overlaps.insert(potential);
      }

  if (potential_overlaps.empty()) return false;

  vector<Point> possible_moves;
  BBox newpos;
  Point translation;

  double currentoverlap = 0.;

  unsigned overlaps = 0;
  BBox totaloverlap;
  for (set<unsigned>::iterator p = potential_overlaps.begin();
       p != potential_overlaps.end(); ++p) {
    cell2.clear();
    cell2 += pl[macros[*p].id];
    cell2.add(pl[macros[*p].id].x + macros[*p].width,
              pl[macros[*p].id].y + macros[*p].height);

    intersection = cell1.intersect(cell2);
    area = intersection.getArea();
    if (lessOrEqualDouble(area, 0.)) continue;
    totaloverlap.expandToInclude(intersection);
    currentoverlap += area;
    ++overlaps;

    // try 8 possibilities for removing the overlap
    if (snapMacrosX)
      translation.x = ceil(intersection.getWidth() / siteSpacing) * siteSpacing;
    else
      translation.x = intersection.getWidth();
    translation.y = 0.;
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);

    if (snapMacrosX)
      translation.x =
          -ceil(intersection.getWidth() / siteSpacing) * siteSpacing;
    else
      translation.x = -intersection.getWidth();
    translation.y = 0.;
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);

    translation.x = 0.;
    if (snapMacrosY)
      translation.y = ceil(intersection.getHeight() / rowSpacing) * rowSpacing;
    else
      translation.y = intersection.getHeight();
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);

    translation.x = 0.;
    if (snapMacrosY)
      translation.y = -ceil(intersection.getHeight() / rowSpacing) * rowSpacing;
    else
      translation.y = -intersection.getHeight();
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);

    if (snapMacrosX)
      translation.x = ceil(intersection.getWidth() / siteSpacing) * siteSpacing;
    else
      translation.x = intersection.getWidth();
    if (snapMacrosY)
      translation.y = ceil(intersection.getHeight() / rowSpacing) * rowSpacing;
    else
      translation.y = intersection.getHeight();
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);

    if (snapMacrosX)
      translation.x = ceil(intersection.getWidth() / siteSpacing) * siteSpacing;
    else
      translation.x = intersection.getWidth();
    if (snapMacrosY)
      translation.y = -ceil(intersection.getHeight() / rowSpacing) * rowSpacing;
    else
      translation.y = -intersection.getHeight();
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);

    if (snapMacrosX)
      translation.x =
          -ceil(intersection.getWidth() / siteSpacing) * siteSpacing;
    else
      translation.x = -intersection.getWidth();
    if (snapMacrosY)
      translation.y = ceil(intersection.getHeight() / rowSpacing) * rowSpacing;
    else
      translation.y = intersection.getHeight();
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);

    if (snapMacrosX)
      translation.x =
          -ceil(intersection.getWidth() / siteSpacing) * siteSpacing;
    else
      translation.x = -intersection.getWidth();
    if (snapMacrosY)
      translation.y = -ceil(intersection.getHeight() / rowSpacing) * rowSpacing;
    else
      translation.y = -intersection.getHeight();
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);
  }

  if (overlaps > 1) {
    // try 8 possibilities for removing the total overlap
    if (snapMacrosX)
      translation.x = ceil(totaloverlap.getWidth() / siteSpacing) * siteSpacing;
    else
      translation.x = totaloverlap.getWidth();
    translation.y = 0.;
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);

    if (snapMacrosX)
      translation.x =
          -ceil(totaloverlap.getWidth() / siteSpacing) * siteSpacing;
    else
      translation.x = -totaloverlap.getWidth();
    translation.y = 0.;
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);

    translation.x = 0.;
    if (snapMacrosY)
      translation.y = ceil(totaloverlap.getHeight() / rowSpacing) * rowSpacing;
    else
      translation.y = totaloverlap.getHeight();
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);

    translation.x = 0.;
    if (snapMacrosY)
      translation.y = -ceil(totaloverlap.getHeight() / rowSpacing) * rowSpacing;
    else
      translation.y = -totaloverlap.getHeight();
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);

    if (snapMacrosX)
      translation.x = ceil(totaloverlap.getWidth() / siteSpacing) * siteSpacing;
    else
      translation.x = totaloverlap.getWidth();
    if (snapMacrosY)
      translation.y = ceil(totaloverlap.getHeight() / rowSpacing) * rowSpacing;
    else
      translation.y = totaloverlap.getHeight();
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);

    if (snapMacrosX)
      translation.x = ceil(totaloverlap.getWidth() / siteSpacing) * siteSpacing;
    else
      translation.x = totaloverlap.getWidth();
    if (snapMacrosY)
      translation.y = -ceil(totaloverlap.getHeight() / rowSpacing) * rowSpacing;
    else
      translation.y = -totaloverlap.getHeight();
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);

    if (snapMacrosX)
      translation.x =
          -ceil(totaloverlap.getWidth() / siteSpacing) * siteSpacing;
    else
      translation.x = -totaloverlap.getWidth();
    if (snapMacrosY)
      translation.y = ceil(totaloverlap.getHeight() / rowSpacing) * rowSpacing;
    else
      translation.y = totaloverlap.getHeight();
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);

    if (snapMacrosX)
      translation.x =
          -ceil(totaloverlap.getWidth() / siteSpacing) * siteSpacing;
    else
      translation.x = -totaloverlap.getWidth();
    if (snapMacrosY)
      translation.y = -ceil(totaloverlap.getHeight() / rowSpacing) * rowSpacing;
    else
      translation.y = -totaloverlap.getHeight();
    newpos = cell1;
    newpos.TranslateBy(translation);
    if (core.contains(newpos)) possible_moves.push_back(translation);
  }

  if (possible_moves.size() == 0) return false;

  double bestoverlap = DBL_MAX;
  double bestmovedist = DBL_MAX, overlap, movedist;
  for (vector<Point>::iterator t = possible_moves.begin();
       t != possible_moves.end(); ++t) {
    movedist = fabs(t->x) + fabs(t->y);
    newpos = cell1;
    newpos.TranslateBy(*t);

    startGridCol = unsigned(floor((newpos.xMin - layout.xMin) / gridXSize));
    startGridRow = unsigned(floor((newpos.yMin - layout.yMin) / gridYSize));
    endGridCol = unsigned(ceil((newpos.xMax - layout.xMin) / gridXSize));
    endGridRow = unsigned(ceil((newpos.yMax - layout.yMin) / gridYSize));
    if (startGridRow > numGridRows) startGridRow = numGridRows;
    if (endGridRow > numGridRows) endGridRow = numGridRows;
    if (startGridCol > numGridCol) startGridCol = numGridCol;
    if (endGridCol > numGridCol) endGridCol = numGridCol;

    potential_overlaps.clear();
    for (rIdx = startGridRow; rIdx < endGridRow; ++rIdx)
      for (cIdx = startGridCol; cIdx < endGridCol; ++cIdx)
        for (idx = 0; idx < grid[rIdx][cIdx].size(); ++idx) {
          const unsigned &potential = grid[rIdx][cIdx][idx];
          if (potential != macro) potential_overlaps.insert(potential);
        }

    overlap = 0.;
    for (set<unsigned>::const_iterator p = potential_overlaps.begin();
         p != potential_overlaps.end(); ++p) {
      cell2.clear();
      cell2 += pl[macros[*p].id];
      cell2.add(pl[macros[*p].id].x + macros[*p].width,
                pl[macros[*p].id].y + macros[*p].height);

      intersection = newpos.intersect(cell2);
      area = intersection.getArea();
      if (greaterThanDouble(area, 0.)) {
        overlap += area;
      }
    }

    if (overlap > bestoverlap) continue;
    if (overlap == bestoverlap && movedist >= bestmovedist) continue;
    bestmovedist = movedist;
    bestoverlap = overlap;
    trans = *t;
  }

  if (bestoverlap == DBL_MAX) return false;  // this shouldn't happen

  improve = currentoverlap - bestoverlap;

  return true;
}

static void getRowBoundaries(const RBPCoreRow &cr, double macroXLoc,
                             double macroWidth, double movement,
                             double &minLeftBoundary,
                             double &maxRightBoundary) {
  minLeftBoundary = DBL_MAX;
  maxRightBoundary = -DBL_MAX;

  for (itRBPSubRow sr = cr.subRowsBegin(); sr != cr.subRowsEnd(); ++sr) {
    double leftBoundary = max(sr->getXStart(), macroXLoc - movement);
    double rightBoundary =
        min(sr->getXEnd() - macroWidth, macroXLoc + movement);
    minLeftBoundary = min(leftBoundary, minLeftBoundary);
    maxRightBoundary = max(rightBoundary, maxRightBoundary);
  }
}

static void roundLeftBoundary(double xMin, double siteSpacing,
                              double &leftBoundary) {
  leftBoundary = ceil((leftBoundary - xMin) / siteSpacing) * siteSpacing + xMin;
}

static void roundRightBoundary(double xMin, double siteSpacing,
                               double &rightBoundary) {
  rightBoundary =
      floor((rightBoundary - xMin) / siteSpacing) * siteSpacing + xMin;
}

void RBPlacement::removeMacroOverlaps(MaxMem *maxMem, bool snapMacrosX,
                                      bool snapMacrosY, bool putCellsInRows,
                                      bool print, unsigned badMoveCap) {
  PlacementWOrient &pl = _placement;
  vector<unsigned> cellIds(getNumCells());

  for (unsigned i = 0; i < getNumCells(); ++i) {
    cellIds[i] = i;
  }

  removeMacroOverlaps(maxMem, pl, cellIds, snapMacrosX, snapMacrosY,
                      putCellsInRows, print, badMoveCap);
}

unsigned RBPlacement::snapMacrosX(PlacementWOrient &pl,
                                  const vector<unsigned> &cellIds,
                                  const BBox &coreArea, double rowSpacing,
                                  double siteSpacing) {
  unsigned snaps = 0;

  vector<set<pair<double, unsigned> > > rightMostSnappedX(_coreRows.size());

  vector<unsigned> ordering(cellIds);
  sort(ordering.begin(), ordering.end(), CompareXWOri(pl));

  for (unsigned i = 0; i < ordering.size(); ++i) {
    unsigned cellId = ordering[i];

    if (_hgWDims->isNodeMacro(cellId) && !isFixed(cellId)) {
      unsigned angle = pl.getOrient(cellId).getAngle();
      bool notRotated = angle == 0 || angle == 180;
      double nodeHeight = notRotated ? _hgWDims->getNodeHeight(cellId)
                                     : _hgWDims->getNodeWidth(cellId);
      double nodeWidth = notRotated ? _hgWDims->getNodeWidth(cellId)
                                    : _hgWDims->getNodeHeight(cellId);
      nodeWidth = ceil(nodeWidth / siteSpacing) * siteSpacing;

      unsigned rowsHigh = static_cast<unsigned>(ceil(nodeHeight / rowSpacing));
      std::size_t yMinIdx = 0, yMaxIdx = 0;
      if (pl[cellId].y < coreArea.yMin) {
        yMinIdx = 0;
        unsigned rowsBelow = static_cast<unsigned>(
            ceil((coreArea.yMin - pl[cellId].y) / rowSpacing));
        yMaxIdx = (rowsBelow >= rowsHigh) ? 0 : rowsHigh - rowsBelow;
      } else {
        yMinIdx = static_cast<unsigned>(
            floor((pl[cellId].y - coreArea.yMin) / rowSpacing));
        yMaxIdx = yMinIdx + rowsHigh;
        unsigned yMaxIdx2 = static_cast<unsigned>(
            floor((pl[cellId].y + nodeHeight - coreArea.yMin) / rowSpacing));
        if (yMaxIdx <= yMaxIdx2) {
          ++yMaxIdx;
        }
      }
      yMinIdx = std::min(yMinIdx, rightMostSnappedX.size());
      yMaxIdx = std::min(yMaxIdx, rightMostSnappedX.size());

      double snappedX =
          floor((pl[cellId].x - coreArea.xMin) / siteSpacing) * siteSpacing +
          coreArea.xMin;

      for (std::size_t j = yMinIdx; j < yMaxIdx; ++j) {
        for (set<pair<double, unsigned> >::reverse_iterator k =
                 rightMostSnappedX[j].rbegin();
             k != rightMostSnappedX[j].rend(); ++k) {
          if (greaterOrEqualDouble(snappedX, k->first)) break;

          unsigned angle2 = pl.getOrient(k->second).getAngle();
          bool notRotated2 = angle2 == 0 || angle2 == 180;
          double nodeHeight2 = notRotated2 ? _hgWDims->getNodeHeight(k->second)
                                           : _hgWDims->getNodeWidth(k->second);

          if (lessThanDouble(pl[cellId].y, pl[k->second].y + nodeHeight2) &&
              lessThanDouble(pl[k->second].y, pl[cellId].y + nodeHeight)) {
            snappedX = k->first;
            break;
          }
        }
      }

      if (!equalDouble(pl[cellId].x, snappedX)) ++snaps;

      pl[cellId].x = snappedX;

      double rightSide = pl[cellId].x + nodeWidth;

      for (std::size_t j = yMinIdx; j < yMaxIdx; ++j) {
        rightMostSnappedX[j].insert(make_pair(rightSide, cellId));
      }
    }
  }

  return snaps;
}

unsigned RBPlacement::snapMacrosY(PlacementWOrient &pl,
                                  const vector<unsigned> &cellIds,
                                  const BBox &coreArea, double rowSpacing,
                                  double siteSpacing) {
  unsigned snaps = 0;

  unsigned numCols =
      static_cast<unsigned>((coreArea.xMax - coreArea.xMin) / siteSpacing);

  vector<set<pair<double, unsigned> > > topMostSnappedY(numCols);

  vector<unsigned> ordering(cellIds);
  sort(ordering.begin(), ordering.end(), CompareYWOri(pl));

  for (unsigned i = 0; i < ordering.size(); ++i) {
    unsigned cellId = ordering[i];

    if (_hgWDims->isNodeMacro(cellId) && !isFixed(cellId)) {
      unsigned angle = pl.getOrient(cellId).getAngle();
      bool notRotated = angle == 0 || angle == 180;
      double nodeHeight = notRotated ? _hgWDims->getNodeHeight(cellId)
                                     : _hgWDims->getNodeWidth(cellId);
      double nodeWidth = notRotated ? _hgWDims->getNodeWidth(cellId)
                                    : _hgWDims->getNodeHeight(cellId);
      nodeHeight = ceil(nodeHeight / rowSpacing) * rowSpacing;

      unsigned colsWide = static_cast<unsigned>(ceil(nodeWidth / siteSpacing));
      std::size_t xMinIdx = 0, xMaxIdx = 0;
      if (pl[cellId].x < coreArea.xMin) {
        xMinIdx = 0;
        unsigned colsLeft = static_cast<unsigned>(
            ceil((coreArea.xMin - pl[cellId].x) / siteSpacing));
        xMaxIdx = (colsLeft >= colsWide) ? 0 : colsWide - colsLeft;
      } else {
        xMinIdx = static_cast<unsigned>(
            floor((pl[cellId].x - coreArea.xMin) / siteSpacing));
        xMaxIdx = xMinIdx + colsWide;
        unsigned xMaxIdx2 = static_cast<unsigned>(
            floor((pl[cellId].x + nodeWidth - coreArea.xMin) / siteSpacing));
        if (xMaxIdx <= xMaxIdx2) {
          ++xMaxIdx;
        }
      }
      xMinIdx = std::min(xMinIdx, topMostSnappedY.size());
      xMaxIdx = std::min(xMaxIdx, topMostSnappedY.size());

      double snappedY =
          floor((pl[cellId].y - coreArea.yMin) / rowSpacing) * rowSpacing +
          coreArea.yMin;

      for (std::size_t j = xMinIdx; j < xMaxIdx; ++j) {
        for (set<pair<double, unsigned> >::reverse_iterator k =
                 topMostSnappedY[j].rbegin();
             k != topMostSnappedY[j].rend(); ++k) {
          if (greaterOrEqualDouble(snappedY, k->first)) break;

          unsigned angle2 = pl.getOrient(k->second).getAngle();
          bool notRotated2 = angle2 == 0 || angle2 == 180;
          double nodeWidth2 = notRotated2 ? _hgWDims->getNodeWidth(k->second)
                                          : _hgWDims->getNodeHeight(k->second);

          if (lessThanDouble(pl[cellId].x, pl[k->second].x + nodeWidth2) &&
              lessThanDouble(pl[k->second].x, pl[cellId].x + nodeWidth)) {
            snappedY = k->first;
            break;
          }
        }
      }

      if (!equalDouble(pl[cellId].y, snappedY)) ++snaps;

      pl[cellId].y = snappedY;

      double topSide = pl[cellId].y + nodeHeight;

      for (std::size_t j = xMinIdx; j < xMaxIdx; ++j) {
        topMostSnappedY[j].insert(make_pair(topSide, cellId));
      }
    }
  }

  return snaps;
}

void RBPlacement::removeMacroOverlaps(MaxMem *maxMem, PlacementWOrient &pl,
                                      const vector<unsigned> &cellIds,
                                      bool snapMacrosX, bool snapMacrosY,
                                      bool putCellsInRows, bool print,
                                      unsigned badMoveCap) {
  if (print && _params.verb.getForActions())
    cout << " Running macro overlap removal ..." << endl;

  Timer tm;

  double startWL = 0.;
  double currentOverlap = 0.;

  if (print && _params.verb.getForMajStats()) {
    startWL = evalHPWL(pl, getHGraph(), true);
  }

#ifdef USEFLOORIST
  if (_params.enableFloorist) {
    currentOverlap =
        calcOverlapSpecific(pl, getBBox(true), cellIds, OnlyMacro, false);

    if (greaterThanDouble(currentOverlap, 0.)) {
      if (print && _params.verb.getForMajStats()) {
        cout << "  HPWL before Floorist overlap remover is " << startWL << endl;
        cout << "  Current macro overlap is " << currentOverlap << endl;
      }

      doFloorist(pl, cellIds);
    }
  }
#endif

  currentOverlap =
      calcOverlapSpecific(pl, getBBox(true), cellIds, OnlyMacro, false);

  if (print && _params.verb.getForMajStats()) {
    cout << "  HPWL before macro overlap remover is "
         << evalHPWL(pl, getHGraph(), true) << endl;
    cout << "  Current macro overlap is " << currentOverlap << endl;
  }

  unsigned numNodes = cellIds.size();
  unsigned numGridRows =
      static_cast<unsigned>(sqrt(static_cast<double>(numNodes)));
  unsigned numGridCol = numGridRows;
  unsigned startGridCol, startGridRow, endGridCol, endGridRow;
  unsigned rIdx, cIdx;
  double siteSpacing = coreRowsBegin()->getSpacing();
  double rowSpacing =
      (coreRowsBegin() + 1)->getYCoord() - coreRowsBegin()->getYCoord();

  BBox layout = getBBox(true), core = getBBox(false);
  double gridXSize = layout.getWidth() / numGridCol;
  double gridYSize = layout.getHeight() / numGridRows;
  double cellWidth, cellHeight;
  bool movedsomething = false, localmoved, ooc = false;

  vector<vector<vector<unsigned> > > grid(
      numGridRows, vector<vector<unsigned> >(numGridCol));
  set<unsigned> potential_overlaps;
  vector<macroData> macros;
  map<unsigned, Point> savedFixedLocations;

  unsigned snapmoves = 0;
  if (snapMacrosX) {
    snapmoves += this->snapMacrosX(pl, cellIds, core, rowSpacing, siteSpacing);
  }

  if (snapMacrosY) {
    snapmoves += this->snapMacrosY(pl, cellIds, core, rowSpacing, siteSpacing);
  }

  unsigned oocmoves = 0;
  for (unsigned j = 0; j < cellIds.size(); ++j) {
    const unsigned &i = cellIds[j];

    if (!(isFixed(i) || _hgWDims->isNodeMacro(i))) continue;

    unsigned angle = pl.getOrient(i).getAngle();
    bool notRotated = angle == 0 || angle == 180;
    cellWidth =
        notRotated ? _hgWDims->getNodeWidth(i) : _hgWDims->getNodeHeight(i);
    cellHeight =
        notRotated ? _hgWDims->getNodeHeight(i) : _hgWDims->getNodeWidth(i);

    if (isFixed(i)) {
      if (snapMacrosX || snapMacrosY) {
        // increase size to ensure row and site alignment
        //  snap to sites and rows (moving left and down)
        double xpos = pl[i].x;
        double ypos = pl[i].y;
        if (snapMacrosX) {
          xpos = floor((pl[i].x - core.xMin) / siteSpacing) * siteSpacing +
                 core.xMin;
        }
        if (snapMacrosY) {
          ypos = floor((pl[i].y - core.yMin) / rowSpacing) * rowSpacing +
                 core.yMin;
        }

        if (pl[i].x != xpos || pl[i].y != ypos) {
          savedFixedLocations[i] = pl[i];
          pl[i].x = xpos;
          pl[i].y = ypos;
        }
      }
    } else {
      double snappedWidth = cellWidth;
      if (snapMacrosX)
        snappedWidth = ceil(snappedWidth / siteSpacing) * siteSpacing;

      double snappedHeight = cellHeight;
      if (snapMacrosY)
        snappedHeight = ceil(snappedHeight / rowSpacing) * rowSpacing;

      // pull out of core macros into the core, no doubt causing overlaps
      if (pl[i].x < core.xMin) {
        pl[i].x = core.xMin;
        ooc = movedsomething = true;
      }
      if (pl[i].x > (core.xMax - snappedWidth)) {
        pl[i].x = core.xMax - snappedWidth;
        ooc = movedsomething = true;
      }
      if (pl[i].y < core.yMin) {
        pl[i].y = core.yMin;
        ooc = movedsomething = true;
      }
      if (pl[i].y > (core.yMax - snappedHeight)) {
        pl[i].y = core.yMax - snappedHeight;
        ooc = movedsomething = true;
      }
      if (ooc) {
        ++oocmoves;
        ooc = false;
      }
    }

    BBox cellOutline(pl[i].x, pl[i].y, pl[i].x + cellWidth,
                     pl[i].y + cellHeight);

    if (greaterThanDouble(cellOutline.intersect(core).getArea(), 0.)) {
      // no need to consider anything that doesn't intersect the core
      startGridCol =
          unsigned(floor((cellOutline.xMin - layout.xMin) / gridXSize));
      startGridRow =
          unsigned(floor((cellOutline.yMin - layout.yMin) / gridYSize));
      endGridCol = unsigned(ceil((cellOutline.xMax - layout.xMin) / gridXSize));
      endGridRow = unsigned(ceil((cellOutline.yMax - layout.yMin) / gridYSize));
      if (startGridRow > numGridRows) startGridRow = numGridRows;
      if (endGridRow > numGridRows) endGridRow = numGridRows;
      if (startGridCol > numGridCol) startGridCol = numGridCol;
      if (endGridCol > numGridCol) endGridCol = numGridCol;

      for (rIdx = startGridRow; rIdx < endGridRow; ++rIdx)
        for (cIdx = startGridCol; cIdx < endGridCol; ++cIdx)
          grid[rIdx][cIdx].push_back(macros.size());

      macros.push_back(macroData(i, isFixed(i), cellWidth, cellHeight));
    }
  }

  if (oocmoves > 0 || snapmoves > 0) {
    movedsomething = true;
    currentOverlap = calcOverlapSpecific(pl, layout, cellIds, OnlyMacro, false);
  }

  if (print && _params.verb.getForMajStats() &&
      (oocmoves > 0 || snapmoves > 0)) {
    cout << "  HPWL after out of core moves and snaps is "
         << evalHPWL(pl, getHGraph(), true) << endl;
    cout << "  Current macro overlap is " << currentOverlap << endl;
  }

  unsigned moves = 0, lastmove = macros.size();
  vector<unsigned> nonimprovingmovesleft(macros.size(), badMoveCap);

  maxMem->update("Macro Overlap Removal");

  if (greaterThanDouble(currentOverlap, 0.)) {
    if (print && _params.verb.getForMajStats() > 1) {
      cout << "  Number of bad moves allowed per movable macro: " << badMoveCap
           << endl;
    }

    do {
      localmoved = false;
      double bestimprove = -DBL_MAX;
      double bestmovedist = DBL_MAX;
      unsigned macro_to_move;
      Point bestTrans;
      for (unsigned i = 0; i < macros.size(); ++i) {
        // if the last move did not improve overlap,
        // disallow that macro from moving this round
        if (macros[i].fixed || i == lastmove) continue;
        Point trans;
        double improve;
        bool needstobemoved = findBestMacroMove(
            pl, grid, gridXSize, gridYSize, numGridRows, numGridCol, i, macros,
            trans, improve, snapMacrosX, snapMacrosY);
        if (!needstobemoved) continue;
        if (improve < bestimprove) continue;
        if (lessOrEqualDouble(improve, 0.) && nonimprovingmovesleft[i] == 0)
          continue;
        double movedist = fabs(trans.x) + fabs(trans.y);
        if (improve == bestimprove && movedist >= bestmovedist) continue;

        bestTrans = trans;
        macro_to_move = i;
        bestimprove = improve;
        bestmovedist = movedist;
      }

      if (bestimprove == -DBL_MAX) {
        if (lastmove != macros.size()) {
          // so we can move back this last bad move
          lastmove = macros.size();
          localmoved = true;
        }
        continue;
      }

      if (lessOrEqualDouble(bestimprove, 0.)) {
        --nonimprovingmovesleft[macro_to_move];
        lastmove = macro_to_move;
      } else
        lastmove = macros.size();

      unsigned macro = macros[macro_to_move].id;

      Point newPos = pl[macro] + bestTrans;

      moveMacro(pl, grid, gridXSize, gridYSize, numGridRows, numGridCol,
                macro_to_move, newPos, macros);

      pl[macro] = newPos;
      localmoved = movedsomething = true;
      ++moves;
    } while (localmoved);

    currentOverlap = calcOverlapSpecific(pl, layout, cellIds, OnlyMacro, false);

    if (print && _params.verb.getForMajStats() > 1) {
      cout << "  HPWL after phase one is " << evalHPWL(pl, getHGraph(), true)
           << endl;
      cout << "  Current macro overlap is " << currentOverlap << endl;
    }
  }

  if (movedsomething) {
    // reset the core rows
    reinstateCoreRows();
    // remove the sites below the movable macros
    for (unsigned i = 0; i < macros.size(); ++i) {
      if (macros[i].fixed) continue;

      splitCoreRowsByCell(pl, macros[i].id);
    }
    if (putCellsInRows) {
      resetPlacement();
    }
  }

  unsigned phasetwomoves = 0;

  if (!equalDouble(currentOverlap, 0.)) {
    movedsomething = false;

    // now for part 2 ....
    vector<unsigned> macrosThatStillOverlap;

    for (unsigned i = 0; i < macros.size(); ++i) {
      if (macros[i].fixed) continue;
      if (doIOverlapHere(pl, grid, gridXSize, gridYSize, numGridRows,
                         numGridCol, i, pl[macros[i].id], macros)) {
        macrosThatStillOverlap.push_back(i);
      }
    }

    maxMem->update("Macro Overlap Removal");

    if (macrosThatStillOverlap.size() > 0) {
      if (print && _params.verb.getForMajStats() > 1) {
        cout << "Going into phase two with " << macrosThatStillOverlap.size()
             << " macros that still overlap" << endl;
      }

      if (_params.macroOverlapDecreasing) {
        sort(macrosThatStillOverlap.begin(), macrosThatStillOverlap.end(),
             OrderByDecreasingArea(macros));
      } else {
        sort(macrosThatStillOverlap.begin(), macrosThatStillOverlap.end(),
             OrderByIncreasingArea(macros));
      }

      for (unsigned i = 0; i < macrosThatStillOverlap.size(); ++i) {
        // check if it still overlaps in its position
        const unsigned &themacro = macrosThatStillOverlap[i];
        if (!doIOverlapHere(pl, grid, gridXSize, gridYSize, numGridRows,
                            numGridCol, themacro, pl[macros[themacro].id],
                            macros)) {
          continue;
        }

        Point pos = pl[macros[themacro].id];
        pos.y = floor((pl[macros[themacro].id].y - core.yMin) / rowSpacing) *
                    rowSpacing +
                core.yMin;
        unsigned baseCrId = findCoreRowIdx(pos);

        double maxmovement =
            2 * (macros[themacro].width + macros[themacro].height);
        double lastMaxMovement = 0.;
        bool done = false;

        while (true) {
          bool goUp = false;
          unsigned above = baseCrId + 1;
          unsigned below = baseCrId;
          while (true) {
            double movementleftup = -DBL_MAX;
            double movementleftdown = -DBL_MAX;
            if (goUp && above < getNumCoreRows()) {
              const RBPCoreRow &cr = getCoreRow(above);
              movementleftup =
                  maxmovement - cr.getYCoord() + pl[macros[themacro].id].y;
              bool exceedsTopOfCore =
                  (cr.getYCoord() + macros[themacro].height > core.yMax)
                      ? true
                      : false;  // <aaronnn> check to exclude rows that would
                                // place macro too high

              if (movementleftup > 0 && !exceedsTopOfCore) {
                double bestDist = DBL_MAX;
                Point bestPos;

                double lastMovementLeftUp = lastMaxMovement - cr.getYCoord() +
                                            _placement[macros[themacro].id].y;
                double previousLeftBoundary = DBL_MAX;
                double previousRightBoundary = -DBL_MAX;
                if (lastMovementLeftUp > 0.) {
                  // get previous boundaries of this row, wrt the
                  // lastMaxMovement
                  getRowBoundaries(cr, _placement[macros[themacro].id].x,
                                   macros[themacro].width, lastMovementLeftUp,
                                   previousLeftBoundary, previousRightBoundary);
                  roundLeftBoundary(core.xMin, siteSpacing,
                                    previousLeftBoundary);
                  roundRightBoundary(core.xMin, siteSpacing,
                                     previousRightBoundary);
                }

                // expand the search space, but exclude search strips we've seen
                // before
                for (itRBPSubRow sr = cr.subRowsBegin(); sr != cr.subRowsEnd();
                     ++sr) {
                  if (sr->getLength() < macros[themacro].width) continue;

                  double leftboundary =
                      max(sr->getXStart(),
                          pl[macros[themacro].id].x - movementleftup);
                  roundLeftBoundary(core.xMin, siteSpacing, leftboundary);

                  double rightboundary =
                      min(sr->getXEnd() - macros[themacro].width,
                          pl[macros[themacro].id].x + movementleftup);
                  roundRightBoundary(core.xMin, siteSpacing, rightboundary);

                  for (double xcoord = leftboundary; xcoord <= rightboundary;
                       xcoord += siteSpacing) {
                    if (xcoord >= previousLeftBoundary &&
                        xcoord <= previousRightBoundary) {
                      // this horizontal strip has been searched earlier, skip
                      // it
                      continue;
                    }

                    Point newPlace(xcoord, cr.getYCoord());

                    if (doIOverlapHere(pl, grid, gridXSize, gridYSize,
                                       numGridRows, numGridCol, themacro,
                                       newPlace, macros)) {
                      continue;
                    }

                    double movedist =
                        fabs(newPlace.x - pl[macros[themacro].id].x) +
                        fabs(newPlace.y - pl[macros[themacro].id].y);
                    if (movedist >= bestDist) {
                      continue;
                    }

                    bestDist = movedist;
                    bestPos = newPlace;
                  }
                }
                if (bestDist != DBL_MAX) {
                  moveMacro(pl, grid, gridXSize, gridYSize, numGridRows,
                            numGridCol, themacro, bestPos, macros);
                  movedsomething = true;
                  pl[macros[themacro].id] = bestPos;
                  ++phasetwomoves;
                  done = true;
                  break;
                }
                ++above;
              }
            } else if (!goUp && below != UINT_MAX) {
              const RBPCoreRow &cr = getCoreRow(below);
              movementleftdown =
                  maxmovement + cr.getYCoord() - pl[macros[themacro].id].y;
              if (movementleftdown > 0) {
                double bestDist = DBL_MAX;
                Point bestPos;

                double lastMovementLeftDown = lastMaxMovement + cr.getYCoord() -
                                              _placement[macros[themacro].id].y;
                double previousLeftBoundary = DBL_MAX;
                double previousRightBoundary = -DBL_MAX;
                if (lastMovementLeftDown > 0.) {
                  // get previous boundaries of this row, wrt to lastMaxMovement
                  getRowBoundaries(cr, _placement[macros[themacro].id].x,
                                   macros[themacro].width, lastMovementLeftDown,
                                   previousLeftBoundary, previousRightBoundary);
                  roundLeftBoundary(core.xMin, siteSpacing,
                                    previousLeftBoundary);
                  roundRightBoundary(core.xMin, siteSpacing,
                                     previousRightBoundary);
                }

                for (itRBPSubRow sr = cr.subRowsBegin(); sr != cr.subRowsEnd();
                     ++sr) {
                  if (sr->getLength() < macros[themacro].width) continue;

                  double leftboundary =
                      max(sr->getXStart(),
                          pl[macros[themacro].id].x - movementleftdown);
                  roundLeftBoundary(core.xMin, siteSpacing, leftboundary);

                  double rightboundary =
                      min(sr->getXEnd() - macros[themacro].width,
                          pl[macros[themacro].id].x + movementleftdown);
                  roundRightBoundary(core.xMin, siteSpacing, rightboundary);

                  for (double xcoord = leftboundary; xcoord <= rightboundary;
                       xcoord += siteSpacing) {
                    if (xcoord >= previousLeftBoundary &&
                        xcoord <= previousRightBoundary) {
                      // this horizontal strip has been considered earlier, skip
                      // it
                      continue;
                    }

                    Point newPlace(xcoord, cr.getYCoord());
                    if (doIOverlapHere(pl, grid, gridXSize, gridYSize,
                                       numGridRows, numGridCol, themacro,
                                       newPlace, macros)) {
                      continue;
                    }

                    double movedist =
                        fabs(newPlace.x - pl[macros[themacro].id].x) +
                        fabs(newPlace.y - pl[macros[themacro].id].y);

                    if (movedist >= bestDist) {
                      continue;
                    }

                    bestDist = movedist;
                    bestPos = newPlace;
                  }
                }
                if (bestDist != DBL_MAX) {
                  moveMacro(pl, grid, gridXSize, gridYSize, numGridRows,
                            numGridCol, themacro, bestPos, macros);
                  movedsomething = true;
                  pl[macros[themacro].id] = bestPos;
                  ++phasetwomoves;
                  done = true;
                  break;
                }
                --below;
              }
            }

            if (movementleftup <= 0. && movementleftdown <= 0.) break;
            goUp = !goUp;
          }

          if (done) break;

          // didn't find anything, expand movement if it isn't already too big
          if (maxmovement >= core.getHalfPerim()) {
            break;
          }

          lastMaxMovement = maxmovement;
          maxmovement += macros[themacro].width + macros[themacro].height;
        }
      }
    }

    currentOverlap = calcOverlapSpecific(pl, layout, cellIds, OnlyMacro, false);

    if (print && _params.verb.getForMajStats() > 1) {
      cout << "  HPWL after phase two is " << evalHPWL(pl, getHGraph(), true)
           << endl;
      cout << "  Current macro overlap is " << currentOverlap << endl;
    }

    if (movedsomething) {
      // reset the core rows
      reinstateCoreRows();
      // remove the sites below the movable macros
      for (unsigned i = 0; i < macros.size(); ++i) {
        if (macros[i].fixed) continue;

        splitCoreRowsByCell(pl, macros[i].id);
      }
      if (putCellsInRows) {
        resetPlacement();
      }
    }
  }

  // reinstate fixed locations if any
  for (map<unsigned, Point>::const_iterator i = savedFixedLocations.begin();
       i != savedFixedLocations.end(); ++i) {
    pl[i->first] = i->second;
  }

  double endWL = 0.;
  if (print && _params.verb.getForMajStats()) {
    endWL = evalHPWL(pl, getHGraph(), true);
  }

  tm.stop();

  if (print && _params.verb.getForMajStats()) {
    cout << " Removing macro overlaps took " << tm << endl;
    if (oocmoves > 0)
      cout << "  Out of core moves performed: " << oocmoves << endl;
    if (snapmoves > 0)
      cout << "  Macro snap moves performed: " << snapmoves << endl;
    if (savedFixedLocations.size() > 0)
      cout << "  Fixed object snap moves performed: "
           << savedFixedLocations.size() << endl;
    if (moves > 0 || phasetwomoves > 0) {
      cout << "  Macro moves performed: " << moves << endl;
      cout << "  Phase two moves performed: " << phasetwomoves << endl;
      cout << "  HPWL after macro overlap removal is " << endWL << endl;
      cout << "  Current macro overlap is " << currentOverlap << endl;
    }
    if (oocmoves > 0 || snapmoves > 0 || moves > 0 || phasetwomoves > 0) {
      double change = (startWL == 0) ? 0 : (endWL - startWL) / startWL;
      if (change < 0)
        cout << "  % decrease in HPWL is "
             << floor(-change * 100000 + 0.5) / 1000 << endl;
      else
        cout << "  % increase in HPWL is "
             << floor(change * 100000 + 0.5) / 1000 << endl;
    }
  }
}

void RBPlacement::moveMacro(const PlacementWOrient &pl,
                            vector<vector<vector<unsigned> > > &grid,
                            const double &gridXSize, const double &gridYSize,
                            const unsigned &numGridRows,
                            const unsigned &numGridCol,
                            const unsigned &macro_to_move, const Point newPos,
                            const vector<macroData> &macros) {
  BBox cell, layout = getBBox(true);

  const unsigned &macro = macros[macro_to_move].id;
  unsigned rIdx, cIdx, startGridCol, startGridRow, endGridCol, endGridRow;

  cell += pl[macro];
  cell.add(pl[macro].x + macros[macro_to_move].width,
           pl[macro].y + macros[macro_to_move].height);

  // remove old grid listings
  startGridCol = unsigned(floor((cell.xMin - layout.xMin) / gridXSize));
  startGridRow = unsigned(floor((cell.yMin - layout.yMin) / gridYSize));
  endGridCol = unsigned(ceil((cell.xMax - layout.xMin) / gridXSize));
  endGridRow = unsigned(ceil((cell.yMax - layout.yMin) / gridYSize));
  if (startGridRow > numGridRows) startGridRow = numGridRows;
  if (endGridRow > numGridRows) endGridRow = numGridRows;
  if (startGridCol > numGridCol) startGridCol = numGridCol;
  if (endGridCol > numGridCol) endGridCol = numGridCol;

  for (rIdx = startGridRow; rIdx < endGridRow; ++rIdx) {
    for (cIdx = startGridCol; cIdx < endGridCol; ++cIdx) {
      vector<unsigned>::iterator tmp =
          find(grid[rIdx][cIdx].begin(), grid[rIdx][cIdx].end(), macro_to_move);
      grid[rIdx][cIdx].erase(tmp);
    }
  }

  // add new grid listings
  cell.clear();
  cell += newPos;
  cell.add(newPos.x + macros[macro_to_move].width,
           newPos.y + macros[macro_to_move].height);

  startGridCol = unsigned(floor((cell.xMin - layout.xMin) / gridXSize));
  startGridRow = unsigned(floor((cell.yMin - layout.yMin) / gridYSize));
  endGridCol = unsigned(ceil((cell.xMax - layout.xMin) / gridXSize));
  endGridRow = unsigned(ceil((cell.yMax - layout.yMin) / gridYSize));
  if (startGridRow > numGridRows) startGridRow = numGridRows;
  if (endGridRow > numGridRows) endGridRow = numGridRows;
  if (startGridCol > numGridCol) startGridCol = numGridCol;
  if (endGridCol > numGridCol) endGridCol = numGridCol;

  for (rIdx = startGridRow; rIdx < endGridRow; ++rIdx) {
    for (cIdx = startGridCol; cIdx < endGridCol; ++cIdx) {
      grid[rIdx][cIdx].push_back(macro_to_move);
    }
  }
}

bool RBPlacement::doIOverlapHere(const PlacementWOrient &pl,
                                 const vector<vector<vector<unsigned> > > &grid,
                                 const double &gridXSize,
                                 const double &gridYSize,
                                 const unsigned &numGridRows,
                                 const unsigned &numGridCol,
                                 const unsigned &macro, const Point newPos,
                                 const vector<macroData> &macros) {
  BBox layout = getBBox(true), cell1, cell2;

  unsigned startGridCol, startGridRow, endGridCol, endGridRow;
  unsigned rIdx, cIdx, idx;

  vector<bool> touched(macros.size(), false);

  cell1 += newPos;
  cell1.add(newPos.x + macros[macro].width, newPos.y + macros[macro].height);

  startGridCol = unsigned(floor((cell1.xMin - layout.xMin) / gridXSize));
  startGridRow = unsigned(floor((cell1.yMin - layout.yMin) / gridYSize));
  endGridCol = unsigned(ceil((cell1.xMax - layout.xMin) / gridXSize));
  endGridRow = unsigned(ceil((cell1.yMax - layout.yMin) / gridYSize));
  if (startGridRow > numGridRows) startGridRow = numGridRows;
  if (endGridRow > numGridRows) endGridRow = numGridRows;
  if (startGridCol > numGridCol) startGridCol = numGridCol;
  if (endGridCol > numGridCol) endGridCol = numGridCol;

  for (rIdx = startGridRow; rIdx < endGridRow; ++rIdx)
    for (cIdx = startGridCol; cIdx < endGridCol; ++cIdx)
      for (idx = 0; idx < grid[rIdx][cIdx].size(); ++idx) {
        const unsigned &potential = grid[rIdx][cIdx][idx];
        if (potential == macro || touched[potential]) continue;
        cell2.clear();
        cell2 += pl[macros[potential].id];
        cell2.add(pl[macros[potential].id].x + macros[potential].width,
                  pl[macros[potential].id].y + macros[potential].height);
        if (cell1.intersects(cell2)) return true;
      }

  return false;
}

#ifdef USEFLOORIST
void RBPlacement::doFloorist(PlacementWOrient &pl,
                             const vector<unsigned> &cellIds) {
  Floorist::Parameters flooristParams;
  flooristParams.verb = _params.verb;

  Floorist solver(flooristParams);

  BBox core = getBBox(false);
  BBox layout = getBBox(true);
  solver.setOutline(core.getWidth(), core.getHeight());
  double originX = layout.xMin;
  double originY = layout.yMin;

  // package up inputs for Floorist
  vector<double> locX, locY, widths, heights;
  vector<bool> fixed;
  vector<unsigned> flooristBlockIDs(cellIds.size(), 0);
  for (unsigned i = 0; i < cellIds.size(); i++) {
    unsigned cellId = cellIds[i];

    if (!_hgWDims->isNodeMacro(cellId) && !isFixed(cellId)) {
      continue;
    }

    unsigned angle = pl.getOrient(cellId).getAngle();
    bool notRotated = angle == 0 || angle == 180;
    double cellWidth = notRotated ? _hgWDims->getNodeWidth(cellId)
                                  : _hgWDims->getNodeHeight(cellId);
    double cellHeight = notRotated ? _hgWDims->getNodeHeight(cellId)
                                   : _hgWDims->getNodeWidth(cellId);
    double cellX = pl[cellId].x;
    double cellY = pl[cellId].y;

    if (isFixed(cellId) && pl[cellId].x > core.xMax) continue;
    if (isFixed(cellId) && pl[cellId].x + cellWidth < core.xMin) continue;
    if (isFixed(cellId) && pl[cellId].y > core.yMax) continue;
    if (isFixed(cellId) && pl[cellId].y + cellHeight < core.yMin) continue;

    // cout << "ANDBG RBPlace pushing back locX: " << pl[cellId].x - originX <<
    // endl; cout << "ANDBG RBPlace pushing back locY: " << pl[cellId].y -
    // originY << endl;

    // trim fixed cells that are partially in the core
    if (isFixed(cellId) && pl[cellId].x + cellWidth > core.xMax)
      cellWidth = core.xMax - pl[cellId].x;
    if (isFixed(cellId) && pl[cellId].x < core.xMin) {
      cellWidth = cellWidth - (core.xMin - pl[cellId].x);
      cellX += (core.xMin - pl[cellId].x);
    }
    if (isFixed(cellId) && pl[cellId].y + cellHeight > core.yMax)
      cellHeight = core.yMax - pl[cellId].y;
    if (isFixed(cellId) && pl[cellId].y < core.yMin) {
      cellHeight = cellHeight - (core.yMin - pl[cellId].y);
      cellY += (core.yMin - pl[cellId].y);
    }

    fixed.push_back(isFixed(cellId));

    locX.push_back(cellX - originX);
    locY.push_back(cellY - originY);

    widths.push_back(cellWidth);
    heights.push_back(cellHeight);

    flooristBlockIDs[i] = locX.size() - 1;
  }

  if (locX.empty())  // nothing to be done
    return;

  solver.setBlocks(locX, locY, widths, heights, fixed);
  solver.run();

  // get result and update placement
  solver.getSolutionPts(locX, locY);

  for (unsigned i = 0; i < cellIds.size(); i++) {
    unsigned cellId = cellIds[i];
    if (!_hgWDims->isNodeMacro(cellId) && !isFixed(cellId)) continue;

    unsigned srcIdx = flooristBlockIDs[i];

    if (!isFixed(cellId)) {
      pl[cellId].x = locX[srcIdx] + originX;
      pl[cellId].y = locY[srcIdx] + originY;
    }
  }
}
#endif
