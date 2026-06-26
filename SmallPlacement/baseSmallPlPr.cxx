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

// Small Placement Interface
//  - Small Placement Problem
//    created by Andrew Caldwell  981206
//  990429  ilm  added colons to the format
//               so that it looks similar to SmallPartProblem;
//               added indentation
//  990429  ilm  changed format version number from 1.0 to 1.1
//  990501  ilm  renamed into BaseSmallPlacementProblem

#include "baseSmallPlPr.h"

#include "ABKCommon/abkio.h"
#include "ABKCommon/infolines.h"

#include "ABKCommon/pathDelims.h"
#include "HGraphWDims/hgWDims.h"

#include <string>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::setw;
using std::string;
using std::vector;

static bool isAbsoluteSmallPlacementPath(const char* fileName) {
  return fileName != NULL &&
         (fileName[0] == pathDelimWindows || fileName[0] == pathDelimUnix);
}

static string smallPlacementFileName(const char* baseFileName,
                                     const char* suffix) {
  return string(baseFileName) + suffix;
}

static string makeSmallPlacementAuxPath(const string& dir,
                                        const char* fileName) {
  if (fileName == NULL) return string();
  if (isAbsoluteSmallPlacementPath(fileName)) return string(fileName);
  return dir + fileName;
}

void BaseSmallPlacementProblem::readAUX(const char* auxFileName) {
  ifstream auxFile(auxFileName);
  abkfatal(auxFile, "can't open auxfile for reading\n");

  bool success = false;
  char* newPln = NULL;
  char* newPla = NULL;
  char* newPlc = NULL;

  string dir;
  const char* auxFileDirEnd = strrchr(auxFileName, pathDelim);
  if (auxFileDirEnd) {
    dir.assign(auxFileName, auxFileDirEnd - auxFileName + 1);
  }

  bool found = false;
  int lineNo = 1;
  char oneLine[1023];
  string word1, word2;

  while (!found && !auxFile.eof()) {
    auxFile >> eathash(lineNo) >> word1 >> noeol(lineNo) >> word2 >>
        noeol(lineNo);
    abkfatal(word2 == ":",
             " Error in aux file: space-separated column expected\n");

    abkfatal(!strcasecmp(word1.c_str(), "SmallPlaceProb"),
             " Expecting \"SmallPlaceProb\" -- other formats not recognized\n");
    {
      found = true;
      auxFile.getline(oneLine, 1023);
      unsigned len = strlen(oneLine), fileNum = 0;
      bool space = true;
      char* fileNames[3] = {NULL, NULL, NULL};
      unsigned j;
      for (j = 0; j != len; j++) {
        if (isspace(oneLine[j])) {
          if (!space) oneLine[j] = '\0';
          space = true;
        } else if (space) {
          space = false;
          abkfatal(fileNum < 3, " Too many filenames in AUX file");
          fileNames[fileNum++] = oneLine + j;
        }
      }

      for (j = 0; j != fileNum; j++) {
        if (strstr(fileNames[j], ".pln")) newPln = fileNames[j];
        if (strstr(fileNames[j], ".PLN")) newPln = fileNames[j];
        if (strstr(fileNames[j], ".plc")) newPlc = fileNames[j];
        if (strstr(fileNames[j], ".PLC")) newPlc = fileNames[j];
        if (strstr(fileNames[j], ".pla")) newPla = fileNames[j];
        if (strstr(fileNames[j], ".PLA")) newPla = fileNames[j];
      }

      abkfatal(newPln != NULL && newPlc != NULL,
               "missing a file in smallPartProb");

      string tmpPln = makeSmallPlacementAuxPath(dir, newPln);
      string tmpPla = makeSmallPlacementAuxPath(dir, newPla);
      string tmpPlc = makeSmallPlacementAuxPath(dir, newPlc);

      cout << " Reading PLN File " << tmpPln << endl;
      ifstream plnStr(tmpPln.c_str());
      success = readPLN(plnStr);
      plnStr.close();

      cout << " Reading PLC File " << tmpPlc << endl;
      ifstream plcStr(tmpPlc.c_str());
      success &= readPLC(plcStr);
      plcStr.close();

      if (!tmpPla.empty()) {
        cout << " Reading PLA File " << tmpPla << endl;
        ifstream plaStr(tmpPla.c_str());
        success &= readPLA(plaStr);
        plaStr.close();
      } else
        _cellWidths =
            vector<double>(_netlist.getNumCells(), _rows[0].siteInterval);
    }
  }

  abkwarn(success,
          "Reading BaseSmallPlacementProblem from files may have failed\n");
}

bool BaseSmallPlacementProblem::readPLN(ifstream& plnStream) {
  abkfatal(plnStream, "can't open .pln file for reading\n");

  int lineNo = 1;
  unsigned numNets = 0, numCells = 0;
  plnStream >> eathash(lineNo) >> needcaseword("UCLA", lineNo) >>
      needcaseword("pln", lineNo) >> skiptoeol;
  lineNo++;
  plnStream >> eathash(lineNo) >> needcaseword("Cells", lineNo) >>
      needcaseword(":", lineNo) >> my_isnumber(lineNo) >> numCells >>
      needeol(lineNo);
  plnStream >> eathash(lineNo) >> needcaseword("Nets", lineNo) >>
      needcaseword(":", lineNo) >> my_isnumber(lineNo) >> numNets >>
      needeol(lineNo);

  _netlist.setup(numCells, numNets);
  _solnBuffer.push_back(SmallPlacementSolution(numCells));

  unsigned netIdx;
  unsigned cellIdx;
  Point newPt;
  BBox newBox;
  unsigned netDegree;

  while (!plnStream.eof() && plnStream.peek() != -1) {
    netDegree = 0;

    plnStream >> eathash(lineNo) >> needcaseword("Net", lineNo) >>
        my_isnumber(lineNo) >> netIdx >> needword(":", lineNo);

    plnStream >> eatblank;

    while (!plnStream.eof() && plnStream.peek() >= '0' &&
           plnStream.peek() <= '9') {
      plnStream >> my_isnumber(lineNo) >> cellIdx;
      plnStream >> needword("(", lineNo) >> my_isnumber(lineNo) >> newPt.x >>
          my_isnumber(lineNo) >> newPt.y;
      newBox = newPt;

      plnStream >> eatblank;
      if (plnStream.peek() != ')') {
        plnStream >> my_isnumber(lineNo) >> newPt.x >> my_isnumber(lineNo) >>
            newPt.y >> needword(")", lineNo);
        newBox += newPt;
      }
      plnStream >> skiptoeol >> eatblank;

      lineNo++;
      while (plnStream.peek() == ' ' || plnStream.peek() == '\t' ||
             plnStream.peek() == '\n' || plnStream.peek() == '\r')
        plnStream.get();

      _netlist.setPinOffset(cellIdx, netIdx, newBox);
      netDegree++;
    }

    abkassert(netDegree >= 1, "nets must have degree >= 1");
  }

  return true;
}

bool BaseSmallPlacementProblem::readPLA(ifstream& plaStream) {
  abkfatal(plaStream, "can't open .pla file for reading\n");

  int lineNo = 1;
  unsigned numCells = 0;
  string cName;
  plaStream >> eathash(lineNo) >> needcaseword("UCLA", lineNo) >>
      needcaseword("pla", lineNo) >> skiptoeol;
  lineNo++;

  plaStream >> eathash(lineNo) >> needcaseword("Cells", lineNo) >>
      needword(":", lineNo) >> my_isnumber(lineNo) >> numCells >>
      needeol(lineNo);

  abkfatal(numCells == _netlist.getNumCells(),
           "numCells in .pla and .pln do not match\n");

  _cellWidths = vector<double>(numCells, DBL_MAX);
  _fixed = bit_vector(numCells, 0);

  unsigned cellIdx;
  double cellWidth;
  bool isFixed = 0;

  for (unsigned c = 0; c < numCells; c++) {
    plaStream >> needcaseword("Cell") >> my_isnumber(lineNo) >> cellIdx >>
        needword(":", lineNo) >> my_isnumber(lineNo) >> cellWidth;
    while (!plaStream.eof() && plaStream.peek() != '\n' &&
           plaStream.peek() != '\r') {
      plaStream >> eatblank;
      if (plaStream.peek() != '\n' && plaStream.peek() != '\r') {
        plaStream >> cName;
        if (!strcasecmp(cName.c_str(), "FIXED")) {
          isFixed = 1;
          break;
        } else {
          abkfatal2(0, "Expecting a fixed cell at ", lineNo);
          break;
        }
      }
    }

    plaStream >> skiptoeol;

    lineNo++;
    abkassert(cellIdx < numCells, "cellIdx too large");

    if (isFixed == 1) _fixed[cellIdx] = 1;

    _cellWidths[cellIdx] = cellWidth;

    isFixed = 0;
  }

  return true;
}

bool BaseSmallPlacementProblem::readPLC(ifstream& plcStream) {
  abkfatal(plcStream, "can't open .plc file for reading");

  int lineNo = 1;

  plcStream >> eathash(lineNo) >> needcaseword("UCLA", lineNo) >>
      needcaseword("plc", lineNo) >> skiptoeol;
  lineNo++;

  unsigned numRows;
  unsigned totalNets;
  unsigned netBBoxes;

  plcStream >> eathash(lineNo) >> needcaseword("Rows", lineNo) >>
      needword(":", lineNo) >> my_isnumber(lineNo) >> numRows >>
      needeol(lineNo);

  plcStream >> eathash(lineNo) >> needcaseword("Total", lineNo) >>
      needcaseword("Nets") >> needword(":", lineNo) >> my_isnumber(lineNo) >>
      totalNets >> needeol(lineNo);

  plcStream >> eathash(lineNo) >> needcaseword("Net", lineNo) >>
      needcaseword("BBoxes") >> needword(":", lineNo) >> my_isnumber(lineNo) >>
      netBBoxes >> needeol(lineNo);

  plcStream >> eatblank;

  abkfatal(totalNets == _netlist.getNumNets(),
           "num nets from pln and plc files do not match\n");

  _rows = vector<SmallPlacementRow>(numRows);

  unsigned rowIdx;
  double xMin, yMin, xMax, yMax, siteInterval;
  unsigned ctr;
  for (ctr = 0; ctr < numRows; ctr++) {
    plcStream >> eathash(lineNo) >> needcaseword("Row", lineNo) >>
        my_isnumber(lineNo) >> rowIdx >> needword(":", lineNo) >>
        my_isnumber(lineNo) >> xMin >> my_isnumber(lineNo) >> yMin >>
        my_isnumber(lineNo) >> xMax >> my_isnumber(lineNo) >> yMax >>
        my_isnumber(lineNo) >> siteInterval >> needeol(lineNo);

    abkassert(rowIdx < numRows, "row index too large");

    BBox rowBBox(xMin, yMin, xMax, yMax);
    _rows[rowIdx] = SmallPlacementRow(rowBBox, siteInterval);
  }

  unsigned netIdx;
  _netTerminalBBoxes = vector<BBox>(totalNets);
  for (ctr = 0; ctr < netBBoxes; ctr++) {
    plcStream >> eathash(lineNo) >> needcaseword("Net", lineNo) >>
        my_isnumber(lineNo) >> netIdx >> needword(":", lineNo) >>
        my_isnumber(lineNo) >> xMin >> my_isnumber(lineNo) >> yMin >>
        my_isnumber(lineNo) >> xMax >> my_isnumber(lineNo) >> yMax >>
        needeol(lineNo);

    abkassert(netIdx < totalNets, "netIdx is too large");

    BBox netBBox(xMin, yMin, xMax, yMax);
    _netTerminalBBoxes[netIdx] = netBBox;
  }

  return true;
}

void BaseSmallPlacementProblem::save(const char* baseFileName) {
  const string newAUX = smallPlacementFileName(baseFileName, ".aux");
  const string newPLN = smallPlacementFileName(baseFileName, ".pln");
  const string newPLC = smallPlacementFileName(baseFileName, ".plc");
  const string newPLA = smallPlacementFileName(baseFileName, ".pla");

  bool success;

  ofstream plnStream(newPLN.c_str());
  abkfatal(plnStream, "can't open pln file for writing");
  cout << " Writing " << newPLN << endl;
  success = savePLN(plnStream);

  ofstream plcStream(newPLC.c_str());
  abkfatal(plcStream, "can't open plc file for writing");
  cout << " Writing " << newPLC << endl;
  success &= savePLC(plcStream);

  ofstream plaStream(newPLA.c_str());
  abkfatal(plaStream, "can't open pla file for writing");
  cout << " Writing " << newPLA << endl;
  success &= savePLA(plaStream);

  if (success) {
    ofstream auxFile(newAUX.c_str());
    abkfatal(auxFile, " Can't open AUX file for writing");
    cout << " Writing " << newAUX << endl;
    auxFile << "SmallPlaceProb : " << newPLN << " " << newPLC << " " << newPLA
            << endl;
  }
}

bool BaseSmallPlacementProblem::savePLN(ofstream& plnStream) {
  plnStream << "UCLA pln 1.1 " << endl
            << TimeStamp() << User() << Platform() << endl;

  plnStream << "Cells : " << _netlist.getNumCells() << endl;
  plnStream << "Nets : " << _netlist.getNumNets() << endl << endl;

  for (unsigned n = 0; n < _netlist.getNumNets(); n++) {
    plnStream << "Net " << setw(4) << n << " : ";
    for (unsigned c = 0; c < _netlist.getNumCells(); c++) {
      const BBox& offset = _netlist.getPinOffset(c, n);
      if (!offset.isEmpty()) {
        if (offset.getBottomLeft() == offset.getTopRight()) {
          plnStream << setw(4) << c << " ( " << setw(10) << offset.xMin << " "
                    << setw(10) << offset.yMin << " )" << endl;
        } else {
          plnStream << setw(4) << c << " ( " << setw(10) << offset.xMin << " "
                    << setw(10) << offset.yMin << " " << setw(10) << offset.xMax
                    << " " << setw(10) << offset.yMax << " )" << endl;
        }
      }
    }
    plnStream << endl;
  }

  return true;
}

bool BaseSmallPlacementProblem::savePLC(ofstream& plcStream) {
  plcStream << "UCLA plc 1.1 " << endl
            << TimeStamp() << User() << Platform() << endl;

  unsigned n, numBBoxes = 0;
  for (n = 0; n < _netTerminalBBoxes.size(); n++)
    if (!_netTerminalBBoxes[n].isEmpty()) numBBoxes++;

  plcStream << "Rows : " << _rows.size() << endl;
  plcStream << "Total Nets : " << _netlist.getNumNets() << endl;
  plcStream << "Net BBoxes : " << numBBoxes << endl << endl;

  for (unsigned r = 0; r < _rows.size(); r++) {
    plcStream << "Row " << r << " : " << _rows[r].xMin << " " << _rows[r].yMin
              << " " << _rows[r].xMax << " " << _rows[r].yMax << " "
              << _rows[r].siteInterval << endl;
  }

  for (n = 0; n < _netTerminalBBoxes.size(); n++) {
    const BBox& box = _netTerminalBBoxes[n];
    if (!box.isEmpty())
      plcStream << "Net " << setw(4) << n << " : " << setw(10) << box.xMin
                << " " << setw(10) << box.yMin << " " << setw(10) << box.xMax
                << " " << setw(10) << box.yMax << endl;
  }

  return true;
}

bool BaseSmallPlacementProblem::savePLA(ofstream& plaStream) {
  plaStream << "UCLA pla 1.1 " << endl
            << TimeStamp() << User() << Platform() << endl;

  plaStream << "Cells : " << _netlist.getNumCells() << endl << endl;

  for (unsigned c = 0; c < _cellWidths.size(); c++) {
    plaStream << "Cell " << setw(4) << c << " : " << setw(8) << _cellWidths[c];
    if (_fixed[c]) plaStream << setw(10) << "FIXED";

    plaStream << endl;
  }
  return true;
}

ostream& operator<<(ostream& os, const BaseSmallPlacementProblem& prob) {
  os << "Small Placement Problem" << endl << endl;
  os << "Netlist/Pin Offsets" << endl;
  os << prob.getNetlist() << endl;
  os << "Cell Widths" << endl;
  os << prob.getCellWidths() << endl;
  os << "Rows " << endl;
  for (unsigned r = 0; r < prob._rows.size(); r++)
    os << r << "  " << prob._rows[r];
  os << "Net Bounding Boxes" << endl;
  for (unsigned n = 0; n < prob._netTerminalBBoxes.size(); n++)
    if (!prob._netTerminalBBoxes[n].isEmpty())
      os << n << "  " << prob._netTerminalBBoxes[n];
    else
      os << n << " none" << endl;

  return os;
}

// added by sadya to plot the small placement solution
void BaseSmallPlacementProblem::plot(const char* baseFileName) {
  const string newGPL = smallPlacementFileName(baseFileName, ".plt");
  double cellWidth = 0;
  double rowHeight = 1;
  double xMin = DBL_MAX;
  double xMax = -DBL_MAX;
  double rowsXMin = DBL_MAX;
  double rowsXMax = -DBL_MAX;
  double rowsYMin = DBL_MAX;
  double rowsYMax = -DBL_MAX;

  Point netStarPoint;
  vector<Point> netPoints;
  unsigned c = 0;

  Placement& placement = getBestSoln().placement;
  abkfatal(getBestSoln().populated,
           "you can't very well plot an unpopulated placement, can you");
  const SmallPlacementNetlist& netlist = getNetlist();
  double WL = calculateWL(placement);
  rowHeight = _rows[0].yMax - _rows[0].yMin;
  double netStep = (1.5 * rowHeight) / netlist.getNumNets();

  ofstream gpFile(newGPL.c_str());

  gpFile << "#This file should be used as a script for gnuplot " << endl;
  gpFile << "#(See http://www.gnuplot.info/ for details)" << endl;
  gpFile << "#" << TimeStamp() << User() << Platform() << endl << endl;
  gpFile << "set title ' " << baseFileName << " HPWL= " << WL
         << " , #Cells= " << getNumCells() << " ,#Nets= " << getNumNets()
         << endl
         << endl;
  gpFile << "set nokey " << endl << endl;
  gpFile << "# uncomment these two lines starting with \"set\"" << endl;
  gpFile << "# to save an EPS file for inclusion into a latex document" << endl;
  gpFile << "# set terminal postscript eps color solid 10" << endl;
  gpFile << "# set output \"" << baseFileName << ".eps\"" << endl << endl;

  gpFile << "# uncomment these two lines starting with \"set\"" << endl;
  gpFile << "# to save an PS file for viewing and printing" << endl;
  gpFile << "# set terminal postscript portrait color solid 8" << endl;
  gpFile << "# set output \"" << baseFileName << ".ps\"" << endl
         << endl
         << endl;

  for (c = 0; c < getNumCells(); c++) {
    gpFile << "set label \"" << c << "\" at "
           << placement[c].x + getCellWidth(c) / 3 << " , "
           << placement[c].y + rowHeight / 3 << endl;
  }

  gpFile << " plot[:][:] '-' w l, '-' w l 3 " << endl;

  for (c = 0; c < getNumCells(); c++) {
    cellWidth = getCellWidth(c);
    gpFile << placement[c].x << " " << placement[c].y + 0.05 * rowHeight
           << endl;
    gpFile << placement[c].x + 0.05 * cellWidth << " " << placement[c].y
           << endl;
    gpFile << placement[c].x + 0.95 * cellWidth << " " << placement[c].y
           << endl;
    gpFile << placement[c].x + cellWidth << " "
           << placement[c].y + 0.05 * rowHeight << endl;
    gpFile << placement[c].x + cellWidth << " "
           << placement[c].y + 0.95 * rowHeight << endl;
    gpFile << placement[c].x + 0.95 * cellWidth << " "
           << placement[c].y + rowHeight << endl;
    gpFile << placement[c].x + 0.05 * cellWidth << " "
           << placement[c].y + rowHeight << endl;
    gpFile << placement[c].x << " " << placement[c].y + 0.95 * rowHeight
           << endl;
    gpFile << placement[c].x << " " << placement[c].y + 0.05 * rowHeight << endl
           << endl;

    if ((placement[c].y + rowHeight) > rowsYMax)
      rowsYMax = placement[c].y + rowHeight;
    if (placement[c].y < rowsYMin) rowsYMin = placement[c].y;
  }

  gpFile << "EOF" << endl;

  double instHeight = rowsYMax - rowsYMin;

  // get the span of the instance
  for (unsigned i = 0; i < _rows.size(); ++i) {
    if (_rows[i].xMin < rowsXMin) rowsXMin = _rows[i].xMin;
    if (_rows[i].xMax > rowsXMax) rowsXMax = _rows[i].xMax;
  }

  // plot the nets
  for (unsigned n = netlist.getNumNets(); n != 0;) {
    --n;
    xMin = DBL_MAX;
    xMax = -DBL_MAX;
    netStarPoint.x = 0;
    netStarPoint.y = 0;

    for (c = 0; c < getNumCells(); c++) {
      const BBox& offset = netlist.getPinOffset(c, n);
      if (!offset.isEmpty()) {
        cellWidth = getCellWidth(c);

        if (placement[c].x + netlist.getPinOffset(c, n).xMin < xMin)
          xMin = placement[c].x + netlist.getPinOffset(c, n).xMin;
        if (placement[c].x + netlist.getPinOffset(c, n).xMax > xMax)
          xMax = placement[c].x + netlist.getPinOffset(c, n).xMax;

        double xMid = 0.5 * (offset.xMin + offset.xMax);
        double yMid = 0.5 * (offset.yMin + offset.yMax);

        Point netPoint(placement[c].x + xMid, placement[c].y + yMid);
        netPoints.push_back(netPoint);
        netStarPoint.x += placement[c].x + xMid;
        netStarPoint.y += placement[c].y + yMid;

        // Point netPoint(placement[c].x+netlist.getPinOffset(c,n).x,
        //	       placement[c].y+netlist.getPinOffset(c,n).y);
        // netPoints.push_back(netPoint);
        // netStarPoint.x += placement[c].x+netlist.getPinOffset(c,n).x;
        // netStarPoint.y += placement[c].y+netlist.getPinOffset(c,n).y;
      }
    }

    netStarPoint.x /= netPoints.size();
    netStarPoint.y /= netPoints.size();

    if (netPoints.size() > 1) {
      double starPtWidth = (rowsXMax - rowsXMin) * 0.005;
      netStarPoint.y += instHeight + (n + 1) * netStep;
      for (unsigned i = 0; i < netPoints.size(); ++i) {
        gpFile << netPoints[i].x << " " << netPoints[i].y << endl;
        gpFile << netStarPoint.x << " " << netStarPoint.y << endl;
        gpFile << netPoints[i].x << " " << netPoints[i].y << endl << endl;
      }
      gpFile << netStarPoint.x - starPtWidth << " " << netStarPoint.y << endl;
      gpFile << netStarPoint.x << " " << netStarPoint.y - starPtWidth << endl;
      gpFile << netStarPoint.x + starPtWidth << " " << netStarPoint.y << endl;
      gpFile << netStarPoint.x << " " << netStarPoint.y + starPtWidth << endl;
      gpFile << netStarPoint.x - starPtWidth << " " << netStarPoint.y << endl
             << endl;
    }
    // Net BBox span

    BBox NetBBox = getNetTerminalBBox(n);
    if (fabs(NetBBox.xMin) < DBL_MAX - 10 &&
        fabs(NetBBox.xMax) < DBL_MAX - 10) {
      if (NetBBox.xMin < rowsXMin)
        NetBBox.xMin = rowsXMin;
      else if (NetBBox.xMin > rowsXMax)
        NetBBox.xMin = rowsXMax;

      if (NetBBox.xMax > rowsXMax)
        NetBBox.xMax = rowsXMax;
      else if (NetBBox.xMax < rowsXMin)
        NetBBox.xMax = rowsXMin;

      gpFile << netStarPoint.x << " " << netStarPoint.y << endl;
      gpFile << (NetBBox.xMin + NetBBox.xMax) / 2 << " "
             << netStarPoint.y + (n + 1) * netStep + instHeight << endl;
      gpFile << netStarPoint.x << " " << netStarPoint.y << endl << endl;

      gpFile << NetBBox.xMin << " "
             << netStarPoint.y + (n + 1) * netStep + instHeight << endl;
      gpFile << NetBBox.xMax << " "
             << netStarPoint.y + (n + 1) * netStep + instHeight << endl;
      gpFile << NetBBox.xMin << " "
             << netStarPoint.y + (n + 1) * netStep + instHeight << endl
             << endl;

      gpFile << NetBBox.xMin << " "
             << netStarPoint.y + (n + 1) * netStep + instHeight << endl;
      gpFile << (NetBBox.xMin + NetBBox.xMax) / 2 << " "
             << netStarPoint.y + (n + 1) * netStep + instHeight + 0.95 * netStep
             << endl;
      gpFile << NetBBox.xMin << " "
             << netStarPoint.y + (n + 1) * netStep + instHeight << endl
             << endl;

      gpFile << NetBBox.xMax << " "
             << netStarPoint.y + (n + 1) * netStep + instHeight << endl;
      gpFile << (NetBBox.xMin + NetBBox.xMax) / 2 << " "
             << netStarPoint.y + (n + 1) * netStep + instHeight + 0.95 * netStep
             << endl;
      gpFile << NetBBox.xMax << " "
             << netStarPoint.y + (n + 1) * netStep + instHeight << endl
             << endl;
    }
    netPoints.clear();
  }

  gpFile << "EOF" << endl;
  gpFile << "pause -1 'Press any key' " << endl;

  cout << "Writing " << newGPL << endl;
  gpFile.close();
}

struct less_mag  //: public binary_function<Point, Point, bool>
{
  bool operator()(Point pt1, Point pt2) { return pt1.x < pt2.x; }
};

// do not use this function. Instead use similar one in RBPlace package
void BaseSmallPlacementProblem::remOverlaps(void) {
  // this function checks if any overlaps present in the placement
  // if yes then attempts to remove them

  // the following code sorts the cells according to xlocation
  vector<Point> cellsInfo;
  unsigned c = 0;
  Placement& placement = getBestSoln().placement;
  for (c = 0; c < getNumCells(); c++) {
    Point cellInfo(placement[c].x, c);
    cellsInfo.push_back(cellInfo);
  }
  std::sort(cellsInfo.begin(), cellsInfo.end(), less_mag());

  //
  for (c = 0; c < getNumCells() - 1; c++) {
    unsigned idx = static_cast<unsigned>(cellsInfo[c].y);
    if (cellsInfo[c + 1].x < cellsInfo[c].x + getCellWidth(idx))
    // overlap found
    {
      double xDiff = cellsInfo[c].x + getCellWidth(idx) - cellsInfo[c + 1].x;
      for (unsigned i = c + 1; i < getNumCells(); ++i) {
        cellsInfo[i].x += xDiff;
        placement[static_cast<unsigned>(cellsInfo[i].y)].x += xDiff;
      }
    }
  }
}

// this function is added to calculate the WL of the placement
double BaseSmallPlacementProblem::calculateWL(const Placement& placement) {
  double wl = 0;

  BBox netBBox;

  for (unsigned n = _netlist.getNumNets(); n != 0;) {
    n--;
    netBBox = getNetTerminalBBox(n);
    for (unsigned c = 0; c < getNumCells(); c++) {
      const BBox& offset = _netlist.getPinOffset(c, n);
      if (!offset.isEmpty()) {
        netBBox += Point(placement[c].x + _netlist.getPinOffset(c, n).xMin,
                         placement[c].y + _netlist.getPinOffset(c, n).yMin);
        netBBox += Point(placement[c].x + _netlist.getPinOffset(c, n).xMax,
                         placement[c].y + _netlist.getPinOffset(c, n).yMax);
      }
    }
    wl += netBBox.getHalfPerim();
  }
  return wl;
}
