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

// 990917 aec  added to RBPlacement
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <cctype>
#include <cmath>
#include <fstream>

#include "ABKCommon/abkassert.h"
#include "RBPlacement.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::vector;

#ifdef _MSC_VER
#ifndef rint
#define rint(a) floor((a) + 0.5)
#endif
#endif

bool nocase_compare(char c1, char c2) { return toupper(c1) == toupper(c2); }

void RBPlacement::placeTerms()
// place the terminals around the edge of the layout region.
{
  unsigned k, numNodes = _hgWDims->getNumNodes();
  for (k = 0; k < numNodes; k++) _placement[k] = Point(0, 0);

  unsigned numRows = getNumCoreRows();
  abkfatal(numRows, " Zero core rows");
  unsigned numTerms = _hgWDims->getNumTerminals();
  if (numRows == 1) {
    Orient nOri("N");
    bool onLeft = true;
    for (k = 0; k != numTerms; k++) {
      std::string fotl = "Fixed__on__the__left";
      std::string fotr = "Fixed__on__the__right";
      if (getHGraph().getNodeNames()[k].size() == fotl.size() &&
          equal(getHGraph().getNodeNames()[k].begin(),
                getHGraph().getNodeNames()[k].end(), fotl.begin(),
                nocase_compare)) {
        _placement[k] = Point(-1.0, 0.0);
        _placement.setOrient(k, nOri);
      } else if (getHGraph().getNodeNames()[k].size() == fotl.size() &&
                 equal(getHGraph().getNodeNames()[k].begin(),
                       getHGraph().getNodeNames()[k].end(), fotr.begin(),
                       nocase_compare)) {
        _placement[k] = Point(numNodes - numTerms, 0);
        _placement.setOrient(k, nOri);
      } else {
        abkwarn(
            0, "In linear placement terms need to be fixed on left or right\n");
        if (onLeft == true)
          _placement[k] = Point(-1.0, 0.0);
        else
          _placement[k] = Point(numNodes - numTerms, 0);
        _placement.setOrient(k, nOri);
        onLeft = !onLeft;
      }
    }
    return;
  }

  // first, find the BBox of the layout region.
  BBox layoutBBox;
  for (itRBPCoreRow cr = coreRowsBegin(); cr != coreRowsEnd(); cr++) {
    Point crOrigin(cr->subRowsBegin()->getXStart(), cr->getYCoord());
    Point crTopLeft(cr->_subRows.back().getXEnd(),
                    cr->getYCoord() + cr->getHeight());
    layoutBBox += crOrigin;
    layoutBBox += crTopLeft;
  }

  if (_params.verb.getForMajStats() > 4)
    cout << "Layout BBox is: " << layoutBBox << endl;

  double layoutWidth = layoutBBox.getWidth();
  double layoutHeight = layoutBBox.getHeight();
  double layoutPerim = 2 * layoutWidth + 2 * layoutHeight;

  unsigned termsOnTop, termsOnBottom, termsOnLeft, termsOnRight;

  termsOnTop = termsOnBottom =
      static_cast<unsigned>(rint(numTerms * (layoutWidth / layoutPerim)));

  termsOnLeft = numTerms - (termsOnTop + termsOnBottom);
  termsOnRight = static_cast<unsigned>(ceil(termsOnLeft / 2.0));
  termsOnLeft -= termsOnRight;

  abkfatal(termsOnTop + termsOnBottom + termsOnLeft + termsOnRight == numTerms,
           "Mismatch in number of terminals.\n total != sum on each side");

  cout << "Terminals: " << numTerms << endl;
  cout << "  Top:     " << termsOnTop << endl;
  cout << "  Bottom:  " << termsOnBottom << endl;
  cout << "  Left:    " << termsOnLeft << endl;
  cout << "  Right:   " << termsOnRight << endl;

  // note: uses the assumption that terminals are numbered
  // consecutively from 0

  Point nextLoc;
  double delta;
  double crHeight = coreRowsBegin()->getHeight();
  // assumes that all terminals have the same height & width
  double termHeight = _hgWDims->getNodeHeight(0);
  double termWidth = _hgWDims->getNodeWidth(0);
  unsigned t;

  nextLoc = Point(floor(layoutBBox.xMin), ceil(layoutBBox.yMax + 2 * crHeight));
  delta = layoutWidth / termsOnTop;

  Orient fsOri("FS");
  for (t = 0; t < termsOnTop; t++) {
    _placement.setOrient(t, fsOri);
    _placement[t] = nextLoc;
    _placement[t].x = rint(nextLoc.x);
    if (t && _placement[t].x < _placement[t - 1].x + termWidth) {
      _placement[t].x++;  // can't be off by more than one
      nextLoc.x++;
    }

    nextLoc.x += delta;
  }

  nextLoc = Point(floor(layoutBBox.xMin),
                  floor(layoutBBox.yMin - 2 * crHeight - termHeight));
  delta = layoutWidth / termsOnBottom;
  Orient nOri("N");
  for (; t < termsOnTop + termsOnBottom; t++) {
    _placement.setOrient(t, nOri);
    _placement[t] = nextLoc;
    _placement[t].x = rint(nextLoc.x);
    if (t != termsOnTop && _placement[t].x < _placement[t - 1].x + termWidth) {
      _placement[t].x++;  // can't be off by more than one
      nextLoc.x++;
    }

    nextLoc.x += delta;
  }

  nextLoc = Point(floor(layoutBBox.xMin - 2 * crHeight - termWidth),
                  floor(layoutBBox.yMin));
  delta = layoutHeight / termsOnLeft;
  Orient eOri("E");
  for (; t < termsOnTop + termsOnBottom + termsOnLeft; t++) {
    _placement.setOrient(t, eOri);
    _placement[t] = nextLoc;
    _placement[t].y = rint(nextLoc.y);
    if (_placement[t].y < _placement[t - 1].y + termHeight) {
      _placement[t].y++;  // can't be off by more than one
      nextLoc.y++;
    }

    nextLoc.y += delta;
  }

  nextLoc = Point(ceil(layoutBBox.xMax + 2 * crHeight), floor(layoutBBox.yMin));
  delta = layoutHeight / termsOnRight;
  Orient fwOri("FW");
  for (; t < termsOnTop + termsOnBottom + termsOnLeft + termsOnRight; t++) {
    _placement.setOrient(t, fwOri);
    _placement[t] = nextLoc;
    _placement[t].y = rint(nextLoc.y);
    if (_placement[t].y < _placement[t - 1].y + termHeight) {
      _placement[t].y++;  // can't be off by more than one
      nextLoc.y++;
    }

    nextLoc.y += delta;
  }
}

RBPlacement::RBPlacement(const char* auxFileName, double aspectRatio,
                         double whiteSpace, const Parameters& params)
    : _hgWDims(NULL),
      _populated(false),
      _cellsNotInRows(0),
      _placement(0),
      _params(params) {
  if (aspectRatio != 0.0)
    abkfatal(aspectRatio >= 0.1 && aspectRatio <= 10,
             "AspectRatio should be between 0.0 or 1/10th and 10");
  abkfatal(whiteSpace >= 0 && whiteSpace <= 100,
           "Whitespace should be a percentage, between 0 and 100%");

  Parameters newParams(params);
  newParams.makeAllSrcSnk = false;
  _hgWDims = new HGraphWDimensions(auxFileName, NULL, NULL, newParams);
  _hgWDimsAuxName = auxFileName;
  _ownHGWDimsAuxName = false;

  unsigned numCoreCells = _hgWDims->getNumNodes() - _hgWDims->getNumTerminals();

  _placement = PlacementWOrient(_hgWDims->getNumNodes());
  _isCoreCell = bit_vector(_hgWDims->getNumNodes(), true);
  _isInSubRow = bit_vector(_hgWDims->getNumNodes(), true);
  _isFixed = bit_vector(_hgWDims->getNumNodes(), false);
  _storElt = bit_vector(_hgWDims->getNumNodes(), false);
  _cellCoords = std::vector<RBCellCoord>(_hgWDims->getNumNodes());

  unsigned t;
  for (t = 0; t != _hgWDims->getNumTerminals(); t++) {
    _isFixed[t] = true;
    _isCoreCell[t] = false;
  }

  if (_params.verb.getForActions()) {
    cout << "Loaded netlist only from auxFile" << endl;
    cout << " Netlist contains:" << endl;
    cout << "    " << _hgWDims->getNumTerminals() << " terminals" << endl;
    cout << "    " << numCoreCells << " core cells" << endl;
  }

  // check that the netlist meets the dimension requirments
  double padHeight = _hgWDims->getNodeHeight(0);
  double padWidth = _hgWDims->getNodeWidth(0);
  double coreHeight = _hgWDims->getNodeHeight(_hgWDims->getNumTerminals());
  bool dimsAreOK = true;
  unsigned n;

  for (n = 0; n < _hgWDims->getNumTerminals(); n++) {
    if (_hgWDims->getNodeHeight(n) != padHeight ||
        _hgWDims->getNodeWidth(n) != padWidth) {
      cout << "Pad node " << n << "(" << _hgWDims->getNodeNameByIndex(n)
           << ") is not the same size as other pads" << endl;
      dimsAreOK = false;
    }
  }

  double totalWidth;
  totalWidth = 0.0;
  double coreCellArea = 0;
  for (n = _hgWDims->getNumTerminals(); n < _hgWDims->getNumNodes(); n++) {
    double height = _hgWDims->getNodeHeight(n);
    if (height != coreHeight) {
      cout << "Core node " << n << "(" << _hgWDims->getNodeNameByIndex(n)
           << ") is not the same height as other core nodes" << endl;
      dimsAreOK = false;
    }
    double width = _hgWDims->getNodeWidth(n);
    if ((rint(width) - width) != 0) {
      cout << "Core node " << n << "(" << _hgWDims->getNodeNameByIndex(n)
           << ") has non-integer width" << endl;
      dimsAreOK = false;
    }
    totalWidth += width;

    coreCellArea += width * height;
  }
  dimsAreOK = true;
  abkfatal(dimsAreOK, "input netlist failed dimension test");

  // compute the layout dimensions

  // double coreCellArea  = totalWidth * coreHeight;
  bool areaOK = (coreCellArea != 0.);
  abkfatal(areaOK, "No core cells to be placed. Exiting.. ");
  double targetLayoutArea = coreCellArea / (1.0 - (whiteSpace / 100.0));
  unsigned numCoreRows = 1;
  if (aspectRatio > 0.0)
    numCoreRows = static_cast<unsigned>(
        rint(sqrt(targetLayoutArea / aspectRatio) / coreHeight));

  double layoutHeight = numCoreRows * coreHeight;
  double layoutWidth = rint(targetLayoutArea / layoutHeight);
  double actualLayoutArea = layoutWidth * layoutHeight;
  // layoutWidth must be integer as the sites wll have width 1
  if (_params.verb.getForMajStats()) {
    cout << "All nodes meet dimension requirments" << endl;
    cout << "Core cell height(==site height) is " << coreHeight << endl;
    cout << "Total core cell width is " << totalWidth << " giving an"
         << " ave width of " << totalWidth / numCoreCells << endl
         << endl;
    cout << "Target Dimensions:" << endl;
    cout << "  Area:   " << targetLayoutArea << endl;
    cout << "  WS%:    " << whiteSpace << endl;
    cout << "  AR:     " << aspectRatio << endl << endl;
    cout << "Actual Dimensions:" << endl;
    cout << "  Width:  " << layoutWidth << endl;
    cout << "  Height: " << layoutHeight << "(" << numCoreRows << " rows)"
         << endl;
    cout << "  Area:   " << actualLayoutArea << endl;
    cout << "  WS%:    "
         << 100 * (actualLayoutArea - coreCellArea) / actualLayoutArea << endl;
    cout << "  AR:     " << layoutWidth / layoutHeight << endl << endl;
  }

  // construct the coreRows
  Point origin(0, 0);
  Orient rowOrients[2] = {"N", "FS"};
  RBPSite crSite(coreHeight, 1.0, Symmetry("Y"));
  double spacing = 1.0;

  _coreRows.reserve(numCoreRows);

  for (unsigned rowId = 0; rowId < numCoreRows; rowId++) {
    _coreRows.push_back(RBPCoreRow(origin.y, rowOrients[rowId % 2], crSite,
                                   _placement, spacing));
    origin.y += coreHeight;

    RBPCoreRow& curCoreRow = _coreRows.back();
    double xStart = origin.x;
    double xEnd = xStart + layoutWidth;
    curCoreRow.appendNewSubRow(xStart, xEnd);
  }

  _hgWDims->expandStdCellWidthToFitSiteWidth(_coreRows[0].getSpacing());
  _hgWDims->setNodeWeightsToAreas(_coreRows[0].getSpacing(),
                                  _coreRows[0].getHeight());

  if (_params.verb.getForActions())
    cout << "Done constructing layout region" << endl;

  // layout construction is done...add in terminal locations
  // if (aspectRatio>0.0)
  placeTerms();
  backUpCoreRows();
  initBBox();

  // check if a .logic file exists, and if so parse it
  ifstream auxFile(auxFileName);
  abkfatal(auxFile, "no aux file found");
  char oneLine[1024];
  auxFile.getline(oneLine, 1023);
  std::string auxLine(oneLine);
  std::string::size_type loc = auxLine.find(".logic");
  if (loc != std::string::npos) {
    std::string::size_type start_loc = auxLine.rfind(" ", loc);
    ++start_loc;
    std::string::size_type end_loc = auxLine.find(" ", start_loc);
    std::string logicFName(auxLine, start_loc, end_loc);
    cout << "DPDEBUG, logic file name: \"" << logicFName << "\"" << endl;
    parseLogic(logicFName.c_str());
  }
  // end parse .logic
}
