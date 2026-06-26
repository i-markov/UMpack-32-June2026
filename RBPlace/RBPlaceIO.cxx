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

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <strings.h>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkio.h"
#include "ABKCommon/infolines.h"
#include "ABKCommon/paramproc.h"
#include "ABKCommon/pathDelims.h"
#include "CongestionMaps/CongestionMaps.h"
#include "Geoms/bbox.h"
#include "HGraph/hgFNode.h"
#include "RBPlacement.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::ostringstream;
using std::map;
using std::max;
using std::min;
using std::ofstream;
using std::setw;
using std::string;
using std::vector;

namespace {
inline string makeName(const char* baseFileName, const char* suffix) {
  return string(baseFileName) + suffix;
}
}  // namespace

RBPlacement::RBPlacement(const char* auxFileName, const Parameters& params)
    : _hgWDims(NULL),
      _ownHGWDimsAuxName(true),
      _populated(false),
      _cellsNotInRows(0),
      _placement(0),
      _isFixed(0),
      _isCoreCell(0),
      _storElt(0),
      _params(params),
      origFileName(auxFileName),
      dirName() {
  ifstream auxFile(auxFileName);
  abkfatal(auxFile, "no aux file found");

  const char* auxFileDirEnd = strrchr(auxFileName, pathDelim);
  if (auxFileDirEnd) {
    dirName = string(auxFileName, 0, auxFileDirEnd - auxFileName + 1);
  }

  bool found = false;
  int i = 0, lineNo = 1;
  cout << "Reading " << auxFileName << endl;
  char oneLine[1023];
  string word1, word2;

  string newScl, newPl, newLGT, newBlocks, newLogic;

  while (!found && !auxFile.eof()) {
    auxFile >> eathash(lineNo) >> word1 >> noeol(lineNo) >> word2 >>
        noeol(lineNo);
    abkfatal(word2 == ":",
             " Error in aux file: space-separated column expected");
    if (!strcasecmp(word1.c_str(), "CD")) {
      auxFile >> word1;
      auxFile >> needeol(lineNo++);
      if (word1[0] == pathDelimWindows || word1[0] == pathDelimUnix)
        dirName = word1;
      else
        dirName += word1;
      char fDel[2];
      fDel[0] = pathDelim;
      fDel[1] = '\0';
      if (word1[word1.length() - 1] != pathDelimWindows ||
          word1[0] == pathDelimUnix)
        dirName += fDel;
    } else if (!strcasecmp(word1.c_str(), "LEFDEF")) {
      cerr << endl;
      cerr << "This utility (from the RBPlace package)";
      cerr << " cannot read LEFDEF files." << endl;
      cerr << "Please try LDToBookConvert.exe";
      cerr << " (a.k.a. RBPlFromDBTest0.exe)" << endl;
      cerr << "or MetaPlacer executables --- each of them can plot." << endl;

      exit(1);
    } else if (!strcasecmp(word1.c_str(), "RowBasedPlacement")) {
      found = true;
      auxFile.getline(oneLine, 1023);
      unsigned j, len = strlen(oneLine), fileNum = 0;
      bool space = true;
      char* fileNames[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
      for (j = 0; j != len; j++) {
        if (isspace(oneLine[j])) {
          if (!space) oneLine[j] = '\0';
          space = true;
        } else if (space) {
          abkfatal(fileNum < 6, " Too many filenames in AUX file");
          fileNames[fileNum++] = oneLine + j;
          space = false;
        }
      }

      for (j = 0; j != fileNum; j++) {
        if (strstr(fileNames[j], ".scl")) newScl = fileNames[j];
        if (strstr(fileNames[j], ".SCL")) newScl = fileNames[j];
        if (strstr(fileNames[j], ".pl")) newPl = fileNames[j];
        if (strstr(fileNames[j], ".PL")) newPl = fileNames[j];
        if (strstr(fileNames[j], ".lgt")) newLGT = fileNames[j];
        if (strstr(fileNames[j], ".LGT")) newLGT = fileNames[j];
        if (strstr(fileNames[j], ".blocks")) newBlocks = fileNames[j];
        if (strstr(fileNames[j], ".BLOCKS")) newBlocks = fileNames[j];
        if (strstr(fileNames[j], ".logic")) newLogic = fileNames[j];
        if (strstr(fileNames[j], ".LOGIC")) newLogic = fileNames[j];
      }

      if (!newScl.empty()) {
        if (newScl[0] != pathDelimWindows && newScl[0] != pathDelimUnix) {
          newScl = dirName + newScl;
        }
        allFileNamesInAUX.push_back(newScl);
      }

      if (!newPl.empty()) {
        if (newPl[0] != pathDelimWindows && newPl[0] != pathDelimUnix) {
          newPl = dirName + newPl;
        }
        allFileNamesInAUX.push_back(newPl);
      }

      if (!newLGT.empty()) {
        if (newLGT[0] != pathDelimWindows && newLGT[0] != pathDelimUnix) {
          newLGT = dirName + newLGT;
        }
        allFileNamesInAUX.push_back(newLGT);
      }

      if (!newLogic.empty()) {
        if (newLogic[0] != pathDelimWindows && newLogic[0] != pathDelimUnix) {
          newLogic = dirName + newLogic;
        }
        allFileNamesInAUX.push_back(newLogic);
      }

      if (!newBlocks.empty()) {
        if (newBlocks[0] != pathDelimWindows && newBlocks[0] != pathDelimUnix) {
          newBlocks = dirName + newBlocks;
        }
        allFileNamesInAUX.push_back(newBlocks);
      }
    }
    if (found)
      while (allFileNamesInAUX.size() < 20 && auxFile >> word1) {
        allFileNamesInAUX.push_back(word1);
      }
    if (i++ == 5) abkfatal(0, "Done");
  }

  Parameters newParams(params);
  newParams.makeAllSrcSnk = false;
  _hgWDims = new HGraphWDimensions(auxFileName, NULL, NULL, newParams);
  _hgWDimsAuxName = auxFileName;
  _ownHGWDimsAuxName = false;

  if (!newBlocks.empty() && _params.allowSoftBlocks) {
    _hgWDims->markSelectedAsSoftFromFile(newBlocks);
  }

  _isCoreCell = bit_vector(_hgWDims->getNumNodes(), true);
  fill_n(_isCoreCell.begin(), _hgWDims->getNumTerminals(), false);

  _isFixed = bit_vector(_hgWDims->getNumNodes(), false);
  fill_n(_isFixed.begin(), _hgWDims->getNumTerminals(), true);

  _isInSubRow = bit_vector(_hgWDims->getNumNodes(), false);
  fill_n(_isInSubRow.begin(), _hgWDims->getNumTerminals(), true);

  _storElt = bit_vector(_hgWDims->getNumNodes(), false);

  _cellCoords = std::vector<RBCellCoord>(_hgWDims->getNumNodes());

  constructRows(newScl.c_str());

  if (!newPl.empty()) {
    // have to parse the placement file here, as it may have
    // node names, which only the HGraph knows.
    _placement =
        PlacementWOrient(_hgWDims->getNumNodes(), Point(DBL_MAX, DBL_MAX));

    ifstream plFile(newPl.c_str());
    abkfatal2(plFile, " Could not open ", newPl);
    cout << " Reading " << newPl << " ... " << endl;

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

      _placement[nodeId] = pt;
      if (plFile.peek() == ':' || plFile.peek() == '/') {
        // plFile>>needword(":");
        plFile.get();
        eatblank(plFile);
        if (plFile.peek() == '\n' || (plFile.peek() == '\r')) {
        } else {
          plFile >> oriStr;
          _placement.setOrient(nodeId, Orientation(oriStr.c_str()));
        }
        eatblank(plFile);
        if (plFile.peek() == ':' || plFile.peek() == '/') {
          plFile.get();
          eatblank(plFile);
          if (plFile.peek() == '\n' || (plFile.peek() == '\r')) {
          } else {
            string fixedStr;
            plFile >> fixedStr;
            if (fixedStr.find("FIXED") != string::npos ||
                fixedStr.find("fixed") != string::npos) {
              _isFixed[nodeId] = true;
            }
          }
        }
      }
      plFile >> skiptoeol >> eathash(lineNo);
    }
  } else
    _placement =
        PlacementWOrient(_hgWDims->getNumNodes(), Point(DBL_MAX, DBL_MAX));

  if (!newLGT.empty()) {
    // have to parse the LGT file here, as it may have
    // node names, which only the HGraph knows.

    ifstream LGTFile(newLGT.c_str());
    abkfatal2(LGTFile, " Could not open ", newLGT);
    cout << " Reading " << newLGT << " ... " << endl;

    string nodeName;  //, buf[256];
    int lineno = 0;

    LGTFile >> needcaseword("UCLA") >> needcaseword("LGT") >> needword("1.0") >>
        skiptoeol;
    LGTFile >> eathash(lineno);

    unsigned numStor = 0;
    while (!LGTFile.eof()) {
      LGTFile >> eathash(lineno) >> nodeName >> needcaseword("STOR");
      unsigned nodeId = _hgWDims->getNodeByName(nodeName.c_str()).getIndex();
      _storElt[nodeId] = true;
      numStor++;
      LGTFile >> skiptoeol >> eathash(lineNo);
    }
    cout << " Read " << numStor << " elements " << endl;
  }
  cout << "Finished reading " << auxFileName << " and referenced files" << endl;

  if (!newLogic.empty()) {
    parseLogic(newLogic.c_str());
  }

  markMacros();
  _hgWDims->expandStdCellWidthToFitSiteWidth(_coreRows[0].getSpacing());
  _hgWDims->setNodeWeightsToAreas(_coreRows[0].getSpacing(),
                                  _coreRows[0].getHeight());
  removeSitesBelowFixed();
  backUpCoreRows();
  setCellsInRows();
  populateCC();
  initBBox();
}

struct RowHeightCmp {
  bool operator()(const RBPCoreRow& r1, const RBPCoreRow& r2) {
    return r1.getYCoord() < r2.getYCoord();
  }
};

void RBPlacement::constructRows(const char* sclFileName) {
  ifstream sclFile(sclFileName);
  abkfatal2(sclFile, " Could not open ", sclFileName);
  cout << " Reading " << sclFileName << " ... " << endl;

  int lineNo = 0;
  unsigned numRows = 0;
  string oneLine;

  sclFile >> eathash(lineNo) >> needcaseword("UCLA", lineNo) >>
      needcaseword("scl", lineNo) >> skiptoeol;
  lineNo++;
  sclFile >> eathash(lineNo);

  sclFile >> needcaseword("NumRows") >> needword(":") >> numRows >> skiptoeol;
  _coreRows.reserve(numRows * 2);

  double yCoord = -1.23456;
  double height = -1.23456;
  double width = -1.23456;
  double spacing = -1.23456;
  string oriStr;
  string symStr;
  double xStart = -1.23456;
  unsigned numSites = UINT_MAX;

  lineNo = 3;

  for (unsigned rNum = 0; rNum < numRows; rNum++) {
    oriStr = " ";
    symStr = " ";
    spacing = DBL_MAX;

    sclFile >> needcaseword("CoreRow") >> needcaseword("Horizontal") >>
        skiptoeol;
    lineNo++;

    sclFile >> needcaseword("Coordinate") >> needword(":", lineNo) >>
        my_isnumber(lineNo) >> yCoord >> skiptoeol;
    lineNo++;

    sclFile >> eathash(lineNo);

    if (sclFile.peek() == 'H' || sclFile.peek() == 'h') {
      sclFile >> needcaseword("Height") >> needword(":", lineNo) >>
          my_isnumber(lineNo) >> height >> needeol();
      lineNo++;
    }

    sclFile >> eathash(lineNo) >> oneLine;
    while (strcasecmp(oneLine.c_str(), "SubrowOrigin")) {
      if (!strcasecmp(oneLine.c_str(), "Sitewidth"))
        sclFile >> needword(":", lineNo) >> my_isnumber(lineNo) >> width >>
            skiptoeol;
      else if (!strcasecmp(oneLine.c_str(), "Sitespacing"))
        sclFile >> needword(":", lineNo) >> my_isnumber(lineNo) >> spacing >>
            skiptoeol;
      else if (!strcasecmp(oneLine.c_str(), "Siteorient"))
        sclFile >> needword(":", lineNo) >> oriStr >> skiptoeol;
      else if (!strcasecmp(oneLine.c_str(), "Sitesymmetry"))
        sclFile >> needword(":", lineNo) >> symStr >> skiptoeol;
      else
        abkfatal3(0, "Unexpected word : ", oneLine.c_str(),
                  "\n(expected one of SubRowOrigin, Sitewidth, Siteorient, "
                  "Sitesymmetry)\n");
      lineNo++;
      sclFile >> eathash(lineNo) >> oneLine;
    }

    spacing = (spacing == DBL_MAX) ? width : spacing;

    _coreRows.push_back(
        RBPCoreRow(yCoord * _params.yScale, Orient(oriStr.c_str()),
                   RBPSite(height * _params.yScale, width,
                           Symmetry(symStr.c_str())),
                   _placement, spacing));

    RBPCoreRow& curCoreRow = _coreRows.back();

    while (!strcasecmp(oneLine.c_str(), "SubrowOrigin")) {
      sclFile >> needword(":", lineNo) >> xStart >> needcaseword("Numsites") >>
          needword(":") >> numSites >> skiptoeol;

      double xEnd = (numSites - 1) * spacing + width + xStart;
      curCoreRow.appendNewSubRow(xStart, xEnd);

      sclFile >> oneLine;
      lineNo++;
    }
    abkfatal2(!strcasecmp(oneLine.c_str(), "End") || sclFile.eof(),
              "in SCL file, expected End, got ", oneLine.c_str());
    sclFile >> eatblank;

    lineNo++;
  }

  sort(_coreRows.begin(), _coreRows.end(), RowHeightCmp());

  itRBPCoreRow CR;
  if (_coreRows.size() > 1) {
    for (CR = (coreRowsBegin() + 1); CR != coreRowsEnd(); ++CR) {
      double rowYCoord = (CR)->getYCoord();
      double prevRowYCoord = (CR - 1)->getYCoord();
      if (rowYCoord == prevRowYCoord) {
        cerr << "Two core rows with same Y-coordinate " << rowYCoord;
        abkfatal(0, "Please use subrows.\n");
      }
    }
  }
}

void RBPlacement::saveAsNodesNetsWts(const char* baseFileName, bool fixMacros,
                                     bool saveNodeWts, bool saveNetWts) const {
  if (fixMacros) {
    float maxCoreHeight = getMaxHeightCoreRow();
    // double maxCellWidth = 10*getAvgWidthNode();
    float maxCellWidth = FLT_MAX;
    _hgWDims->markNodesAsMacro(maxCellWidth, maxCoreHeight);
  } else

    _hgWDims->saveAsNodesNetsWts(baseFileName, fixMacros, saveNodeWts,
                                 saveNetWts);

  saveAsSCL(baseFileName);
  saveAsLogic(baseFileName);

  string plName = makeName(baseFileName, ".pl");
  ofstream plOut(plName.c_str());
  plOut << "UCLA pl 1.0" << endl;
  plOut << TimeStamp() << User() << Platform() << endl << endl;

  unsigned i = 0;
  for (i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (_placement[i].x > 1e308 || _placement[i].x > 1e308) continue;

    plOut << _hgWDims->getNodeNameByIndex(i) << "  " << setw(10)
          << _placement[i].x << "  " << setw(10) << _placement[i].y << " : "
          << _placement.getOrient(i) << endl;
  }
  plOut.close();
  bool haveStorElts = false;
  for (i = 0; i < _storElt.size(); i++)
    if (_storElt[i]) {
      haveStorElts = true;
      break;
    }

  if (haveStorElts) {
    string lgtName = string(baseFileName) + ".lgt";

    ofstream lgtOut(lgtName.c_str());
    lgtOut << "UCLA lgt 1.0" << endl;
    lgtOut << TimeStamp() << User() << Platform() << endl << endl;

    unsigned i = 0;
    for (i = 0; i < _hgWDims->getNumNodes(); i++) {
      if (!_storElt[i]) continue;

      lgtOut << setw(20) << _hgWDims->getNodeNameByIndex(i) << "    STOR "
             << endl;
    }
    lgtOut.close();
  }

  string auxName = makeName(baseFileName, ".aux");
  ofstream auxFile(auxName.c_str());

  auxFile << "RowBasedPlacement : ";

  const char* tmp = NULL;
  tmp = strrchr(baseFileName, pathDelim);
  if (tmp++) {
    auxFile << tmp << ".nodes " << tmp << ".nets ";
    if (saveNodeWts || saveNetWts) auxFile << tmp << ".wts ";
    auxFile << tmp << ".pl " << tmp << ".scl ";
    if (haveStorElts) auxFile << tmp << ".lgt ";
    if (!_logicValues.empty()) auxFile << tmp << ".logic ";
  } else {
    auxFile << baseFileName << ".nodes " << baseFileName << ".nets ";
    if (saveNodeWts || saveNetWts) auxFile << baseFileName << ".wts ";
    auxFile << baseFileName << ".pl " << baseFileName << ".scl ";
    if (haveStorElts) auxFile << baseFileName << ".lgt ";
    if (!_logicValues.empty()) auxFile << baseFileName << ".logic ";
  }
  auxFile << endl;

  auxFile.close();
}

void RBPlacement::saveAsNodesNetsWtsWCongInfo(const char* baseFileName,
                                              double whitespaceRatio) {
  CongestionMaps congMap(this);
  congMap.initializeFMatrix();
  congMap.createCongestionMapProb();
  unsigned numNodes = _hgWDims->getNumNodes();
  vector<double> nodesXCongestion(numNodes, 0);
  vector<double> nodesYCongestion(numNodes, 0);
  unsigned i;
  for (i = 0; i < numNodes; ++i) {
    nodesXCongestion[i] = congMap.getHorizCongestion(_placement[i]);
    nodesYCongestion[i] = congMap.getVertCongestion(_placement[i]);
  }

  //  BBox layoutBBox = getBBox(false);
  double layoutArea = getSitesArea();
  double totalNodesArea = _hgWDims->getNodesArea();
  double whitespace = 0.01 * whitespaceRatio * (layoutArea - totalNodesArea);
  double rowHeight = getMaxHeightCoreRow();
  _hgWDims->saveAsNodesNetsWtsWCongInfo(baseFileName, nodesXCongestion,
                                        nodesYCongestion, whitespace,
                                        static_cast<float>(rowHeight));

  saveAsSCL(baseFileName);
  saveAsLogic(baseFileName);

  string plName = makeName(baseFileName, ".pl");
  ofstream plOut(plName.c_str());
  plOut << "UCLA pl 1.0" << endl;
  plOut << TimeStamp() << User() << Platform() << endl << endl;

  for (i = 0; i < _hgWDims->getNumNodes(); ++i) {
    if (_placement[i].x > 1e308 || _placement[i].x > 1e308) continue;

    plOut << _hgWDims->getNodeNameByIndex(i) << "  " << setw(10)
          << _placement[i].x << "  " << setw(10) << _placement[i].y << " : "
          << _placement.getOrient(i) << endl;
  }
  plOut.close();
  bool haveStorElts = false;
  for (i = 0; i < _storElt.size(); ++i)
    if (_storElt[i]) {
      haveStorElts = true;
      break;
    }

  if (haveStorElts) {
    string lgtName = makeName(baseFileName, ".lgt");
    ofstream lgtOut(lgtName.c_str());
    lgtOut << "UCLA lgt 1.0" << endl;
    lgtOut << TimeStamp() << User() << Platform() << endl << endl;

    for (i = 0; i < _hgWDims->getNumNodes(); ++i) {
      if (!_storElt[i]) continue;

      lgtOut << setw(20) << _hgWDims->getNodeNameByIndex(i) << "    STOR "
             << endl;
    }
    lgtOut.close();
  }

  string auxName = makeName(baseFileName, ".aux");
  ofstream auxFile(auxName.c_str());

  auxFile << "RowBasedPlacement : ";

  const char* tmp = NULL;
  tmp = strrchr(baseFileName, pathDelim);
  if (tmp++) {
    auxFile << tmp << ".nodes " << tmp << ".nets " << tmp << ".wts " << tmp
            << ".pl " << tmp << ".scl ";
    if (haveStorElts) auxFile << tmp << ".lgt ";
    if (!_logicValues.empty()) auxFile << tmp << ".logic ";
  } else {
    auxFile << baseFileName << ".nodes " << baseFileName << ".nets "
            << baseFileName << ".wts " << baseFileName << ".pl " << baseFileName
            << ".scl ";
    if (haveStorElts) auxFile << baseFileName << ".lgt ";
    if (!_logicValues.empty()) auxFile << baseFileName << ".logic ";
  }
  auxFile << endl;

  auxFile.close();
}

void RBPlacement::savePlNoDummy(const char* baseFileName) {
  string plName = string(baseFileName) + string(".pl");

  ofstream plOut(plName.c_str());
  plOut << "UCLA pl 1.0" << endl;
  plOut << TimeStamp() << User() << Platform() << endl << endl;

  unsigned i = 0;
  for (i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (_placement[i].x > 1e308 || _placement[i].x > 1e308) continue;

    // TODO: THIS CODE IS INCORRECT, PROBABLY DOESNT WORK
    // string temp_str = _hgWDims->getNodeByIdx(i).getName();
    // if(sscanf(_hgWDims->getNodeByIdx(i).getName().c_str(),
    // "dummy%s",temp_str.c_str())!=1) reverting to old code by sadya. If it
    // doesn't compile on windows then comment next 2 lines and uncomment above
    // 2 lines.
    if (_hgWDims->getNodeNameByIndex(i).compare(0, 5, "dummy") != 0) {
      plOut << _hgWDims->getNodeNameByIndex(i) << "  " << setw(10)
            << _placement[i].x << "  " << setw(10) << _placement[i].y << " : "
            << _placement.getOrient(i) << endl;
    }
  }
  plOut.close();
}

void RBPlacement::saveAsNetDAre(const char* baseFileName) const {
  _hgWDims->saveAsNetDAre(baseFileName);
  saveAsSCL(baseFileName);
  saveAsLogic(baseFileName);

  string plName = makeName(baseFileName, ".pl");
  _placement.save(plName.c_str(), _hgWDims->getNumTerminals());

  string auxName = makeName(baseFileName, ".aux");
  ofstream auxFile(auxName.c_str());

  auxFile << "RowBasedPlacement : ";
  auxFile << baseFileName << ".netD " << baseFileName << ".areM "
          << baseFileName << ".pl " << baseFileName << ".scl" << endl;

  auxFile.close();
}

void RBPlacement::saveAsSCL(const char* baseFileName) const {
  string fname = makeName(baseFileName, ".scl");
  ofstream sclFile(fname.c_str());

  sclFile << "UCLA scl 1.0" << endl;
  sclFile << TimeStamp() << User() << Platform() << endl << endl;

  sclFile << "Numrows : " << _coreRows.size() << endl << endl;

  unsigned r;
  for (r = 0; r < _coreRows.size(); r++) _coreRows[r].save(sclFile);
  sclFile.close();
}

void RBPlacement::saveAsLogic(const char* baseFileName) const {
  if (_logicValues.empty()) {
    return;
  }

  string fname = baseFileName;
  fname += ".logic";
  ofstream logicFile(fname.c_str());

  logicFile << "UCLA logic 1.0" << endl;
  logicFile << TimeStamp() << User() << Platform() << endl << endl;

  for (map<unsigned, string>::const_iterator logicIt = _logicValues.begin();
       logicIt != _logicValues.end(); logicIt++) {
    // get node name or make one up if there isnt one
    // note that at the time of writing May 08, 2006,
    // this name will be identical to what HGraphWDims makes up
    string nodeName = _hgWDims->getNodeNameByIndex(logicIt->first).c_str();
    if (nodeName.empty()) {
      if (_hgWDims->isTerminal(logicIt->first)) {
        ostringstream ss;
        ss << "p" << (logicIt->first + 1);
        nodeName = ss.str();
      } else {
        ostringstream ss;
        ss << "a" << (logicIt->first - _hgWDims->getNumTerminals());
        nodeName = ss.str();
      }
    }
    // end get node name or make one up if there isnt one

    logicFile << setw(10) << nodeName << "  " << (logicIt->second) << endl;
  }
  logicFile.close();
}

// TO DO
void RBPlacement::saveAsSCLCong(const char* baseFileName) {
  removeSitesFromCongestedRgn();
  saveAsSCL(baseFileName);
}

void RBPlacement::removeSitesFromCongestedRgn(double whitespaceRatio) {
  CongestionMaps congMap(this);
  congMap.initializeFMatrix();
  congMap.createCongestionMapProb();

  congMap.sortCongestionBins();
  const std::vector<ctrSortCongMap>& sortedCongBins =
      congMap.getSortedCongestionBins();

  BBox layoutBBox = getBBox(false);
  double layoutArea = getSitesArea();
  double totalNodesArea = _hgWDims->getNodesArea();
  double whitespace = 0.01 * whitespaceRatio * (layoutArea - totalNodesArea);

  //    cout<<whitespaceRatio<<"  "<<whitespace<<endl<<endl;
  // cout<<"layoutArea "<<layoutArea<<" nodesArea "<<totalNodesArea<<endl;
  // cout<<"whitespace "<<(layoutArea-totalNodesArea)/totalNodesArea<<endl;

  double horizGridSize = congMap.getHorizGridSize();
  double vertGridSize = congMap.getVertGridSize();
  for (int s = sortedCongBins.size() - 1; s >= 0; --s) {
    if (whitespace <= 0) break;
    ctrSortCongMap bin = sortedCongBins[s];
    double xlow = max(layoutBBox.xMin, bin.horizIdx * horizGridSize);
    double xhigh = min(layoutBBox.xMax, (bin.horizIdx + 1) * horizGridSize);
    double ylow = max(layoutBBox.yMin, bin.vertIdx * vertGridSize);
    double yhigh = min(layoutBBox.yMax, (bin.vertIdx + 1) * vertGridSize);
    BBox obstacle;
    obstacle.add(xlow, ylow);
    obstacle.add(xhigh, yhigh);
    // cout<<obstacle<<"gridSize "<<horizGridSize<<" "<<vertGridSize<<endl;
    for (unsigned i = 0; i < _coreRows.size(); i++) {
      RBPCoreRow& cr = _coreRows[i];
      double coreRowYMin = cr.getYCoord();
      double coreRowYMax = cr.getYCoord() + cr.getHeight();

      if ((obstacle.yMin > coreRowYMin) || (obstacle.yMax < coreRowYMax))
        continue;

      _splitRowByObstacle(cr, obstacle);
    }
    whitespace -= obstacle.getArea();
  }
  // cout<<"Final whitespace "<<whitespace<<endl;
}

void RBPlacement::saveAsSpatialConstraints(const char* constFileName) const {
  ofstream constFile(constFileName);

  vector<double> validYLocs;
  BBox layoutBBox;

  for (itRBPCoreRow cr = coreRowsBegin(); cr != coreRowsEnd(); cr++) {
    layoutBBox += Point((cr)->getXStart(), (cr)->getYCoord());
    layoutBBox += Point((cr)->getXEnd(), (cr)->getYCoord() + (cr)->getHeight());

    validYLocs.push_back((cr)->getYCoord());
  }

  constFile << "Layout BBox (" << layoutBBox.xMin << ", " << layoutBBox.yMin
            << ", " << layoutBBox.xMax << ", " << layoutBBox.yMax << ")"
            << endl;

  constFile << "validY      (" << validYLocs[0];
  for (unsigned y = 1; y < validYLocs.size(); y++)
    constFile << ", " << validYLocs[y];
  constFile << ")" << endl;

  constFile.close();
}

string RBPlacement::getAuxName() const {
  string::size_type left = origFileName.rfind(pathDelim);
  if (left == string::npos)
    left = 0;
  else
    left++;
  string temp = origFileName.substr(left, string::npos);

  string::size_type right = temp.rfind('.');
  if (right == 0)
    return string("noname");
  else
    return temp.substr(0, right);
}

double RBPlacement::getMaxHeightCoreRow(void) const {
  double maxHeight = -DBL_MAX;
  unsigned r;
  for (r = 0; r < _coreRows.size(); r++) {
    if (maxHeight < _coreRows[r].getHeight())
      maxHeight = _coreRows[r].getHeight();
  }
  return maxHeight;
}

double RBPlacement::getMinWidthNode(void) const {
  double minWidth = DBL_MAX;
  unsigned i;
  for (i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (minWidth > _hgWDims->getNodeWidth(i) && !_hgWDims->isTerminal(i))
      minWidth = _hgWDims->getNodeWidth(i);
  }
  return minWidth;
}

double RBPlacement::getMaxWidthNode(void) const {
  double maxWidth = -DBL_MAX;
  unsigned i;
  for (i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (maxWidth < _hgWDims->getNodeWidth(i) && !_hgWDims->isTerminal(i))
      maxWidth = _hgWDims->getNodeWidth(i);
  }
  return maxWidth;
}

double RBPlacement::getAvgWidthNode(void) const {
  double avgWidth = 0;
  unsigned i;
  for (i = 0; i < _hgWDims->getNumNodes(); i++) {
    avgWidth += _hgWDims->getNodeWidth(i);
  }
  if (_hgWDims->getNumNodes() != 0) avgWidth /= _hgWDims->getNumNodes();

  return avgWidth;
}

void RBPlacement::saveAsNodesNetsWtsUnShred(const char* baseFileName) const {
  double maxCoreRowHeight;
  maxCoreRowHeight = getMaxHeightCoreRow();
  _hgWDims->saveAsNodesNetsWtsUnShred(baseFileName, maxCoreRowHeight);
  saveAsSCL(baseFileName);
  saveAsLogic(baseFileName);

  string plName = makeName(baseFileName, ".pl");
  ofstream plOut(plName.c_str());
  plOut << "UCLA pl 1.0" << endl;
  plOut << TimeStamp() << User() << Platform() << endl << endl;

  unsigned i = 0;
  for (i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (_placement[i].x > 1e308 || _placement[i].x > 1e308) continue;

    // HGFNode& node   = _hgWDims->getNodeByIdx(i);
    string nodeName = _hgWDims->getNodeNameByIndex(i).c_str();
    string actualName = nodeName;
    unsigned long pos = actualName.find('@');

    if (pos != string::npos)  // found
    {
      actualName.erase(pos);
      string subNodeName;

      double avgYLoc = 0;
      double avgXLoc = 0;
      unsigned numSubNodes = 0;
      double nodeHeight = 0;
      // double nodeWidth = _hgWDims->getNodeWidth(i);
      while (i < _hgWDims->getNumNodes()) {
        // HGFNode& subNode   = _hgWDims->getNodeByIdx(i);
        nodeName = _hgWDims->getNodeNameByIndex(i);
        subNodeName = nodeName;
        pos = subNodeName.find('@');
        if (pos == string::npos) break;
        subNodeName.erase(pos);
        if (actualName == subNodeName) break;

        nodeHeight += _hgWDims->getNodeHeight(i);
        avgYLoc += _placement[i].y;
        avgXLoc += _placement[i].x;
        numSubNodes++;
        i++;
      }
      i--;
      avgXLoc /= numSubNodes;
      avgYLoc /= numSubNodes;
      double nodeYLoc = avgYLoc - nodeHeight / 2;
      double nodeXLoc = avgXLoc;

      plOut << actualName << "  ";
      plOut << setw(10) << nodeXLoc << "  ";
      plOut << setw(10) << nodeYLoc << " :  " << _placement.getOrient(i - 1)
            << endl;
    } else {
      plOut << _hgWDims->getNodeNameByIndex(i) << "  " << setw(10)
            << _placement[i].x << "  " << setw(10) << _placement[i].y << " : "
            << _placement.getOrient(i) << endl;
    }
  }
  plOut.close();

  bool haveStorElts = false;
  for (i = 0; i < _storElt.size(); i++)
    if (_storElt[i]) {
      haveStorElts = true;
      break;
    }

  if (haveStorElts) {
    string lgtName = makeName(baseFileName, ".lgt");
    ofstream lgtOut(lgtName.c_str());
    lgtOut << "UCLA lgt 1.0" << endl;
    lgtOut << TimeStamp() << User() << Platform() << endl << endl;

    unsigned i = 0;
    for (i = 0; i < _hgWDims->getNumNodes(); i++) {
      if (!_storElt[i]) continue;

      lgtOut << setw(20) << _hgWDims->getNodeNameByIndex(i) << "    STOR "
             << endl;
    }
    lgtOut.close();
  }

  string auxName = makeName(baseFileName, ".aux");
  ofstream auxFile(auxName.c_str());

  auxFile << "RowBasedPlacement : ";

  const char* tmp = NULL;
  tmp = strrchr(baseFileName, pathDelim);
  if (tmp++) {
    auxFile << tmp << ".nodes " << tmp << ".nets " << tmp << ".wts " << tmp
            << ".pl " << tmp << ".scl ";
    if (haveStorElts) auxFile << tmp << ".lgt ";
    if (!_logicValues.empty()) auxFile << tmp << ".logic ";
  } else {
    auxFile << baseFileName << ".nodes " << baseFileName << ".nets "
            << baseFileName << ".wts " << baseFileName << ".pl " << baseFileName
            << ".scl ";
    if (haveStorElts) auxFile << baseFileName << ".lgt ";
    if (!_logicValues.empty()) auxFile << baseFileName << ".logic ";
  }
  auxFile << endl;

  auxFile.close();
}

void RBPlacement::saveAsNodesNetsWtsShred(const char* baseFileName) const {
  double maxCoreRowHeight;
  maxCoreRowHeight = getMaxHeightCoreRow();
  _hgWDims->saveAsNodesNetsWtsShred(baseFileName, maxCoreRowHeight);
  saveAsSCL(baseFileName);
  saveAsLogic(baseFileName);

  string plName = makeName(baseFileName, ".pl");
  ofstream plOut(plName.c_str());
  plOut << "UCLA pl 1.0" << endl;
  plOut << TimeStamp() << User() << Platform() << endl << endl;

  unsigned i = 0;
  for (i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (_placement[i].x > 1e308 || _placement[i].x > 1e308) continue;

    HGFNode& node = _hgWDims->getNodeByIdx(i);

    if (_hgWDims->getNodeHeight(i) > maxCoreRowHeight &&
        !_hgWDims->isTerminal(i)) {
      double nodeHeight = _hgWDims->getNodeHeight(i);
      double nodeYLoc = _placement[i].y;
      string nodeName = _hgWDims->getNodeNameByIndex(i).c_str();
      unsigned numNewNodes = unsigned(ceil(nodeHeight / maxCoreRowHeight));
      for (unsigned j = 0; j < numNewNodes; ++j) {
        ostringstream ss;
        ss << _hgWDims->getNodeNameByIndex(node.getIndex()) << "@" << j;
        nodeName = ss.str();
        plOut << nodeName << "  ";
        plOut << setw(10) << _placement[i].x << "  ";
        plOut << setw(10) << nodeYLoc << " :  " << _placement.getOrient(i)
              << endl;

        nodeYLoc += maxCoreRowHeight;
      }
    } else {
      plOut << _hgWDims->getNodeNameByIndex(i) << "  " << setw(10)
            << _placement[i].x << "  " << setw(10) << _placement[i].y << " : "
            << _placement.getOrient(i) << endl;
    }
  }
  plOut.close();

  bool haveStorElts = false;
  for (i = 0; i < _storElt.size(); i++)
    if (_storElt[i]) {
      haveStorElts = true;
      break;
    }

  if (haveStorElts) {
    string lgtName = makeName(baseFileName, ".lgt");
    ofstream lgtOut(lgtName.c_str());
    lgtOut << "UCLA lgt 1.0" << endl;
    lgtOut << TimeStamp() << User() << Platform() << endl << endl;

    unsigned i = 0;
    for (i = 0; i < _hgWDims->getNumNodes(); i++) {
      if (!_storElt[i]) continue;

      lgtOut << setw(20) << _hgWDims->getNodeNameByIndex(i) << "    STOR "
             << endl;
    }
    lgtOut.close();
  }

  string auxName = makeName(baseFileName, ".aux");
  ofstream auxFile(auxName.c_str());

  auxFile << "RowBasedPlacement : ";

  const char* tmp = NULL;
  tmp = strrchr(baseFileName, pathDelim);
  if (tmp++) {
    auxFile << tmp << ".nodes " << tmp << ".nets " << tmp << ".wts " << tmp
            << ".pl " << tmp << ".scl ";
    if (haveStorElts) auxFile << tmp << ".lgt ";
    if (!_logicValues.empty()) auxFile << tmp << ".logic ";
  } else {
    auxFile << baseFileName << ".nodes " << baseFileName << ".nets "
            << baseFileName << ".wts " << baseFileName << ".pl " << baseFileName
            << ".scl ";
    if (haveStorElts) auxFile << baseFileName << ".lgt ";
    if (!_logicValues.empty()) auxFile << baseFileName << ".logic ";
  }
  auxFile << endl;

  auxFile.close();
}

void RBPlacement::coarsenRows(double newHeight) {
  std::vector<RBPCoreRow> newRows;
  std::vector<RBPCoreRow> sortedCoreRows(_coreRows);
  sort(sortedCoreRows.begin(), sortedCoreRows.end(), RowHeightCmp());
  std::vector<RBPCoreRow>::iterator rIt = sortedCoreRows.begin();
  double last_height = -1.0e10;
  while (rIt != sortedCoreRows.end()) {
    if (rIt->getHeight() < newHeight) {
      double y_coord = rIt->getYCoord();
      RBPCoreRow tempRow(*rIt);
      tempRow.getSite().setHeight(newHeight);
      newRows.push_back(tempRow);
      ++rIt;
      while ((rIt != sortedCoreRows.end()) &&
             (rIt->getYCoord() < y_coord + newHeight)) {
        last_height = rIt->getYCoord();
        ++rIt;
        if ((rIt != sortedCoreRows.end()) && last_height >= rIt->getYCoord()) {
          cout << "Values: " << last_height << " vs. " << rIt->getYCoord()
               << endl;
          abkassert(last_height < rIt->getYCoord(),
                    "This function assumes _coreRows is sorted");
        }
      }
    } else {
      newRows.push_back(RBPCoreRow(*rIt));
      ++rIt;
    }
  }
  _coreRows = newRows;
}

void RBPlacement::saveAsSCLcoarsenRows(double newHeight,
                                       const string& baseFileName) const {
  // Segment copied from coarsenRows
  std::vector<RBPCoreRow> newRows;
  std::vector<RBPCoreRow> sortedCoreRows(_coreRows);
  sort(sortedCoreRows.begin(), sortedCoreRows.end(), RowHeightCmp());
  std::vector<RBPCoreRow>::const_iterator rIt = sortedCoreRows.begin();
  double last_height = -1.0e10;
  while (rIt != sortedCoreRows.end()) {
    if (rIt->getHeight() < newHeight) {
      double y_coord = rIt->getYCoord();
      RBPCoreRow tempRow(*rIt);
      tempRow.getSite().setHeight(newHeight);
      newRows.push_back(tempRow);
      ++rIt;
      while ((rIt != sortedCoreRows.end()) &&
             (rIt->getYCoord() < y_coord + newHeight)) {
        last_height = rIt->getYCoord();
        ++rIt;
        if ((rIt != sortedCoreRows.end()) && last_height >= rIt->getYCoord()) {
          cout << "Values: " << last_height << " vs. " << rIt->getYCoord()
               << endl;
          abkassert(last_height < rIt->getYCoord(),
                    "This function assumes _coreRows is sorted");
        }
      }
    } else {
      newRows.push_back(RBPCoreRow(*rIt));
      ++rIt;
    }
  }
  // end copy from coarsenRows

  // Segment copied and partially modified from the saveAsSCL function
  string sclFileName(baseFileName);
  sclFileName += ".scl";
  ofstream sclFile(sclFileName.c_str());

  sclFile << "UCLA scl 1.0" << endl;
  sclFile << TimeStamp() << User() << Platform() << endl << endl;

  sclFile << "Numrows : " << newRows.size() << endl << endl;

  unsigned r;
  for (r = 0; r < newRows.size(); r++) newRows[r].save(sclFile);
  sclFile.close();
  // end saveAsSCL copy
}

void RBPlacement::saveAsNodesNetsWtsShredHW(const char* baseFileName,
                                            double maxHeight, double maxWidth,
                                            unsigned singleNetWt) const {
  double maxCoreRowHeight;
  if (maxHeight == 0)
    maxCoreRowHeight = getMaxHeightCoreRow();
  else
    maxCoreRowHeight = maxHeight;

  double minWidthNode = getMinWidthNode();
  // double minWidthNode = ceil(getAvgWidthNode());
  double maxWidthNode = getMaxWidthNode();
  // double avgWidthNode = 5*getAvgWidthNode();
  double avgWidthNode = 100 * getAvgWidthNode();

  if (maxWidthNode / minWidthNode > 50)  // limit the number of horizontal split
    minWidthNode = maxWidthNode / 50;

  if (maxWidth != 0) minWidthNode = maxWidth;

  // minWidthNode = 1;
  _hgWDims->markNodesAsMacro(maxCoreRowHeight, avgWidthNode);
  _hgWDims->saveAsNodesNetsWtsShredHW(baseFileName, minWidthNode,
                                      maxCoreRowHeight, singleNetWt);

  saveAsSCLcoarsenRows(maxCoreRowHeight, baseFileName);

  string plName = makeName(baseFileName, ".pl");
  ofstream plOut(plName.c_str());
  plOut << "UCLA pl 1.0" << endl;
  plOut << TimeStamp() << User() << Platform() << endl << endl;

  unsigned i = 0;
  for (i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (_placement[i].x > 1e308 || _placement[i].x > 1e308) continue;

    HGFNode& node = _hgWDims->getNodeByIdx(i);

    if (_hgWDims->isNodeMacro(i) && !_hgWDims->isTerminal(i)) {
      double nodeHeight = _hgWDims->getNodeHeight(i);
      double nodeWidth = _hgWDims->getNodeWidth(i);
      double nodeYLoc = _placement[i].y;
      double nodeXLoc = _placement[i].x;

      string nodeName = _hgWDims->getNodeNameByIndex(i).c_str();
      unsigned numNewVertNodes = unsigned(ceil(nodeHeight / maxCoreRowHeight));
      unsigned numNewHorizNodes = unsigned(ceil(nodeWidth / minWidthNode));

      for (unsigned j = 0; j < numNewVertNodes; ++j) {
        nodeXLoc = _placement[i].x;
        for (unsigned k = 0; k < numNewHorizNodes; ++k) {
          ostringstream ss;
          ss << _hgWDims->getNodeNameByIndex(node.getIndex()) << "@" << j << "#"
             << k;
          nodeName = ss.str();
          plOut << nodeName << "  ";
          plOut << setw(10) << nodeXLoc << "  ";
          plOut << setw(10) << nodeYLoc << " :  " << _placement.getOrient(i)
                << endl;
          nodeXLoc += minWidthNode;
        }
        nodeYLoc += maxCoreRowHeight;
      }
    } else {
      plOut << _hgWDims->getNodeNameByIndex(i) << "  " << setw(10)
            << _placement[i].x << "  " << setw(10) << _placement[i].y << " : "
            << _placement.getOrient(i) << endl;
    }
  }
  plOut.close();

  bool haveStorElts = false;
  for (i = 0; i < _storElt.size(); i++)
    if (_storElt[i]) {
      haveStorElts = true;
      break;
    }

  if (haveStorElts) {
    string lgtName = makeName(baseFileName, ".lgt");
    ofstream lgtOut(lgtName.c_str());
    lgtOut << "UCLA lgt 1.0" << endl;
    lgtOut << TimeStamp() << User() << Platform() << endl << endl;

    unsigned i = 0;
    for (i = 0; i < _hgWDims->getNumNodes(); i++) {
      if (!_storElt[i]) continue;

      lgtOut << setw(20) << _hgWDims->getNodeNameByIndex(i) << "    STOR "
             << endl;
    }
    lgtOut.close();
  }

  string auxName = makeName(baseFileName, ".aux");
  ofstream auxFile(auxName.c_str());

  auxFile << "RowBasedPlacement : ";

  const char* tmp = NULL;
  tmp = strrchr(baseFileName, pathDelim);
  if (tmp++) {
    auxFile << tmp << ".nodes " << tmp << ".nets " << tmp << ".wts " << tmp
            << ".pl " << tmp << ".scl ";
    if (haveStorElts) auxFile << tmp << ".lgt ";
    if (!_logicValues.empty()) auxFile << tmp << ".logic ";
  } else {
    auxFile << baseFileName << ".nodes " << baseFileName << ".nets "
            << baseFileName << ".wts " << baseFileName << ".pl " << baseFileName
            << ".scl ";
    if (haveStorElts) auxFile << baseFileName << ".lgt ";
    if (!_logicValues.empty()) auxFile << baseFileName << ".logic ";
  }
  auxFile << endl;

  auxFile.close();
}

void RBPlacement::saveAsNodesNetsPlFloorplan(const char* baseFileName) const {
  _hgWDims->saveAsNodesNetsFloorplan(baseFileName);

  string plName = makeName(baseFileName, ".pl");
  savePlacement(plName.c_str());
}

void RBPlacement::saveMacrosAsNodesNetsPlFloorplan(
    const char* baseFileName) const {
  double maxCoreRowHeight = getMaxHeightCoreRow();
  _hgWDims->markNodesAsMacro(FLT_MAX, maxCoreRowHeight);
  _hgWDims->saveMacrosAsNodesNetsFloorplan(baseFileName);

  string plName = makeName(baseFileName, ".pl");
  ofstream plOut(plName.c_str());
  plOut << "UCLA pl 1.0" << endl;
  plOut << TimeStamp() << User() << Platform() << endl << endl;

  unsigned i = 0;
  for (i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (_placement[i].x > 1e308 || _placement[i].x > 1e308) continue;

    if (_hgWDims->isTerminal(i)) {
      plOut << _hgWDims->getNodeNameByIndex(i) << "  " << setw(10)
            << _placement[i].x << "  " << setw(10) << _placement[i].y << " : "
            << _placement.getOrient(i) << endl;
    } else {
      if (_hgWDims->isNodeMacro(i)) {
        plOut << _hgWDims->getNodeNameByIndex(i) << "  " << setw(10)
              << _placement[i].x << "  " << setw(10) << _placement[i].y << " : "
              << _placement.getOrient(i) << endl;
      }
    }
  }
  plOut.close();
}

void RBPlacement::saveAsCplace(const char* baseFileName) {
  itHGFEdgeLocal edge;
  unsigned edgeId;
  string fname = makeName(baseFileName, "_cpin");
  cout << "Saving in Cplace format. Output File is " << fname << endl;
  ofstream cpFile(fname.c_str());

  BBox imageSize = getBBox();

  cpFile << "begin_control" << endl;
  cpFile << "  multilevel 4" << endl;
  cpFile << "  default_number_of_tries 4" << endl;
  cpFile << "  density_target .9" << endl;
  cpFile << "  start_window  " << imageSize.xMin << " " << imageSize.yMin << " "
         << imageSize.xMax << " " << imageSize.yMax << " " << endl;
  cpFile << "  number_of_cuts -9999" << endl;
  cpFile << "  auto_dp easy" << endl;
  cpFile << "  write_capo /tmp/trash\n";
  cpFile << "  write_cplaceout /tmp/cpout\n";
  cpFile << "end_control" << endl << endl;

  cpFile << "begin_orientation_spec" << endl;
  cpFile << "  mirror=no" << endl;
  cpFile << "  rotate=0" << endl;
  cpFile << "end_orientation_spec" << endl << endl;

  cpFile << "begin_prewire" << endl;
  cpFile << "end_prewire" << endl;
  cpFile << "begin_multi_bit_area" << endl;
  cpFile << "end_multi_bit_area" << endl;
  cpFile << "begin_blocked_areas" << endl;
  cpFile << "end_blocked_areas" << endl;
  cpFile << "begin_edge" << endl;
  cpFile << "end_edge" << endl;
  cpFile << "begin_io_level" << endl;
  cpFile << "end_io_level" << endl;
  cpFile << "begin_blocked_io" << endl;
  cpFile << "end_blocked_io" << endl;
  cpFile << "begin_bit_def" << endl;
  cpFile << "end_bit_def" << endl;
  cpFile << "begin_clk_regen" << endl;
  cpFile << "end_clk_regen" << endl << endl;

  cpFile << "begin_image" << endl;
  cpFile << "io_boundary " << imageSize.xMin << " " << imageSize.yMin << " "
         << imageSize.xMax << " " << imageSize.yMax << " " << endl;

  itRBPCoreRow itCore;
  for (itCore = coreRowsBegin(); itCore != coreRowsEnd(); ++itCore) {
    itRBPSubRow itSub;
    for (itSub = itCore->subRowsBegin(); itSub != itCore->subRowsEnd();
         ++itSub) {
      cpFile << "row " << itCore->getYCoord() << " " << itSub->getXStart()
             << " " << itCore->getYCoord() << " " << itSub->getXEnd() << " "
             << itCore->getYCoord() + itCore->getHeight() << " -9999" << endl;
    }
  }
  cpFile << "end_image" << endl;

  cpFile << "begin_net_weight" << endl;
  cpFile << "end_net_weight" << endl;

  cpFile << "begin_object_list" << endl;

  for (unsigned i = 0; i < _hgWDims->getNumNodes(); i++) {
    const HGFNode& node = _hgWDims->getNodeByIdx(i);
    cpFile << "object  " << _hgWDims->getNodeNameByIndex(i) << "  ";
    if (_hgWDims->isTerminal(i))
      cpFile << "fixedio  " << _placement[i].x << "  " << _placement[i].y
             << "  ";
    else
      cpFile << "movable  -9999  -9999  ";

    cpFile << _hgWDims->getNodeWidth(i) << "  " << _hgWDims->getNodeHeight(i)
           << " 0 0 0 0" << endl;
    int pinNum = 0;
    for (edge = node.edgesBegin(); edge != node.edgesEnd(); ++edge) {
      edgeId = (*edge)->getIndex();
      string netName = _hgWDims->getNetNameByIndex(edgeId).c_str();
      ostringstream pinName;
      pinName << "dummy_" << pinNum;
      cpFile << "pin " << pinName.str() << "  ";
      if (!netName.empty())
        cpFile << setw(10) << netName;
      else
        cpFile << " net" << edgeId;

#ifdef SIGNAL_DIRECTIONS
      if (node.isEdgeSrc(*edge))
        cpFile << " O ";
      else if (node.isEdgeSrcSnk(*edge))
        cpFile << " B ";
      else
        cpFile << " I ";
#else
      cpFile << " B ";
#endif

      cpFile << endl;
      Point pinOffset = _hgWDims->edgesOnNodePinOffset(pinNum, i);
      cpFile << "pinrpt " << pinOffset.x << " " << pinOffset.y << endl;
      pinNum++;
    }
  }
  cpFile << "end_object_list" << endl;
  cpFile.close();
}

void RBPlacement::saveAsPlato(const char* baseFileName) {
  string fname = makeName(baseFileName, ".spc");
  cout << "Saving in Plato format. Output File is " << fname << endl;
  ofstream spcFile(fname.c_str());

  BBox imageSize = getBBox(false);

  // write the layout
  spcFile << "add_region_rectangle CORE " << imageSize.xMin << " "
          << imageSize.yMin << " " << imageSize.getWidth() << " "
          << imageSize.getHeight() << " 1 0 0 " << endl;
  for (unsigned i = 0; i < _coreRows.size(); i++) {
    for (unsigned j = 0; j < _coreRows[i].getNumSubRows(); j++) {
      RBPSubRow& sRow = _coreRows[i][j];
      spcFile << "add_region_row CORE " << sRow.getXStart() << " "
              << sRow.getYCoord() << " " << sRow.getLength() << " "
              << sRow.getHeight() << endl;
    }
  }

  // write the cell library
  for (unsigned i = 0; i < _hgWDims->getNumNodes(); ++i) {
    HGFNode& node = _hgWDims->getNodeByIdx(i);
    spcFile << "add_module " << _hgWDims->getNodeNameByIndex(i) << "_proto "
            << _hgWDims->getNodeWidth(i) << " " << _hgWDims->getNodeHeight(i);
    if (_hgWDims->isTerminal(i))
      spcFile << " 0 ";
    else
      spcFile << " 2 ";
    spcFile << "0 " << endl;

    unsigned offset = 0;
    for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd();
         e++, offset++) {
      Point pOffset =
          _hgWDims->edgesOnNodePinOffset(offset, i, _placement.getOrient(i));
      spcFile << "add_pin_celllib_sequential pin" << offset << " 0 0 "
              << pOffset.x << " " << pOffset.y << " 3 " << endl;
    }
  }
  spcFile << "finish_read_celllib" << endl;

  // write the actual netlist
  for (unsigned i = 0; i < _hgWDims->getNumNodes(); ++i) {
    HGFNode& node = _hgWDims->getNodeByIdx(i);
    spcFile << "add_instance " << _hgWDims->getNodeNameByIndex(i) << " "
            << _hgWDims->getNodeNameByIndex(i) << "_proto " << _placement[i].x
            << " " << _placement[i].y << " 0 ";
    if (_hgWDims->isTerminal(i))
      spcFile << "0 0 ";
    else
      spcFile << "1 2 ";
    spcFile << endl;

    unsigned offset = 0;
    for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd();
         e++, offset++) {
      unsigned edgeId = (*e)->getIndex();
      string netName = getHGraph().getNetNames()[edgeId].c_str();
      spcFile << "add_pin pin" << offset << " "
              << _hgWDims->getNodeNameByIndex(node.getIndex());
      if (!netName.empty())
        spcFile << " " << netName;
      else
        spcFile << " edge" << edgeId;
      spcFile << " 3 " << endl;
    }
  }
  // write netlist->region specification
  for (unsigned i = 0; i < _hgWDims->getNumNodes(); ++i) {
    // HGFNode& node = _hgWDims->getNodeByIdx(i);
    if (!_hgWDims->isTerminal(i)) {
      spcFile << "set_instance_location "
              << ((*_hgWDims) << _hgWDims->getNodeByIdx(i)) << " CORE " << endl;
    }
  }
  spcFile << "finish_read_netlist" << endl;

  spcFile.close();
}

void RBPlacement::saveAsNodesNetsWtsWDummy(const char* baseFileName,
                                           double fractWS, bool fixMacros) {
  double layoutArea = getSitesArea();
  double totalNodesArea = _hgWDims->getNodesArea();
  double newWS = (0.01 * fractWS) * layoutArea;
  double oldWS = layoutArea - totalNodesArea;
  cout << "Old White Space is " << 100 * oldWS / layoutArea << " %" << endl;

  //    double areaDummy = 0.01*fractWS*(layoutArea - totalNodesArea);
  double areaDummy = oldWS - newWS;
  if (areaDummy < 0) areaDummy = 0;

  double siteHeight = getMaxHeightCoreRow();
  double siteWidth = _coreRows[0].getSiteWidth();
  if (fixMacros) {
    double maxCoreHeight = getMaxHeightCoreRow();
    _hgWDims->markNodesAsMacro(FLT_MAX, maxCoreHeight);
    _hgWDims->saveAsNodesNetsWtsWDummy(baseFileName, siteWidth, siteHeight,
                                       areaDummy, fixMacros);
  } else
    _hgWDims->saveAsNodesNetsWtsWDummy(baseFileName, siteWidth, siteHeight,
                                       areaDummy, 0);

  saveAsSCL(baseFileName);
  saveAsLogic(baseFileName);

  string plName = makeName(baseFileName, ".pl");
  ofstream plOut(plName.c_str());
  plOut << "UCLA pl 1.0" << endl;
  plOut << TimeStamp() << User() << Platform() << endl << endl;

  unsigned i = 0;
  for (i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (_placement[i].x > 1e308 || _placement[i].x > 1e308) continue;

    plOut << _hgWDims->getNodeNameByIndex(i) << "  " << setw(10)
          << _placement[i].x << "  " << setw(10) << _placement[i].y << " : "
          << _placement.getOrient(i) << endl;
  }
  plOut.close();

  bool haveStorElts = false;
  for (i = 0; i < _storElt.size(); i++)
    if (_storElt[i]) {
      haveStorElts = true;
      break;
    }

  if (haveStorElts) {
    string lgtName = makeName(baseFileName, ".lgt");
    ofstream lgtOut(lgtName.c_str());
    lgtOut << "UCLA lgt 1.0" << endl;
    lgtOut << TimeStamp() << User() << Platform() << endl << endl;

    unsigned i = 0;
    for (i = 0; i < _hgWDims->getNumNodes(); i++) {
      if (!_storElt[i]) continue;

      lgtOut << setw(20) << _hgWDims->getNodeNameByIndex(i) << "    STOR "
             << endl;
    }
    lgtOut.close();
  }

  string auxName = makeName(baseFileName, ".aux");
  ofstream auxFile(auxName.c_str());

  auxFile << "RowBasedPlacement : ";

  const char* tmp = NULL;
  tmp = strrchr(baseFileName, pathDelim);
  if (tmp++) {
    auxFile << tmp << ".nodes " << tmp << ".nets " << tmp << ".wts " << tmp
            << ".pl " << tmp << ".scl ";
    if (haveStorElts) auxFile << tmp << ".lgt ";
    if (!_logicValues.empty()) auxFile << tmp << ".logic ";
  } else {
    auxFile << baseFileName << ".nodes " << baseFileName << ".nets "
            << baseFileName << ".wts " << baseFileName << ".pl " << baseFileName
            << ".scl ";
    if (haveStorElts) auxFile << baseFileName << ".lgt ";
    if (!_logicValues.empty()) auxFile << baseFileName << ".logic ";
  }
  auxFile << endl;

  auxFile.close();
}

void RBPlacement::saveAsNodesNetsWtsWTether(const char* baseFileName,
                                            double fract, double rgnSizePercent,
                                            bool takeNodeConstrFrmFile,
                                            bool takeNetConstrFrmFile,
                                            double tetherNewAR,
                                            double tetherNewWS) {
  vector<unsigned> tetheredCells;
  _hgWDims->saveAsNodesNetsWtsWTether(baseFileName, fract, tetheredCells,
                                      takeNodeConstrFrmFile,
                                      takeNetConstrFrmFile);

  saveAsSCL(baseFileName);
  saveAsLogic(baseFileName);

  string plName = makeName(baseFileName, ".pl");
  ofstream plOut(plName.c_str());
  plOut << "UCLA pl 1.0" << endl;
  plOut << TimeStamp() << User() << Platform() << endl << endl;

  unsigned i = 0;

  BBox layoutBBox = getBBox(false);
  double rgnXMax = rgnSizePercent * layoutBBox.getWidth();
  double rgnYMax = rgnSizePercent * layoutBBox.getHeight();

  reshapePlacement(tetherNewAR, tetherNewWS, "newPl.pl");

  for (i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (_placement[i].x > 1e308 || _placement[i].x > 1e308) continue;

    plOut << _hgWDims->getNodeNameByIndex(i) << "  " << setw(10)
          << _placement[i].x << "  " << setw(10) << _placement[i].y << " : "
          << _placement.getOrient(i) << endl;
  }

  vector<Point> tetheredCellsPlacement;
  for (i = 0; i < tetheredCells.size(); ++i) {
    unsigned nodeIdx = tetheredCells[i];
    Point thisNodePl = _placement[nodeIdx];
    tetheredCellsPlacement.push_back(thisNodePl);
  }

  unsigned tetherTermIdx = 0;
  for (i = 0; i < tetheredCells.size(); ++i) {
    plOut << "dummyTerm" << tetherTermIdx << " "
          << tetheredCellsPlacement[i].x - rgnXMax / 2.0 << " "
          << tetheredCellsPlacement[i].y - rgnYMax / 2.0 << endl;
    ++tetherTermIdx;
    plOut << "dummyTerm" << tetherTermIdx << " "
          << tetheredCellsPlacement[i].x + rgnXMax / 2.0 << " "
          << tetheredCellsPlacement[i].y - rgnYMax / 2.0 << endl;
    ++tetherTermIdx;
    plOut << "dummyTerm" << tetherTermIdx << " "
          << tetheredCellsPlacement[i].x + rgnXMax / 2.0 << " "
          << tetheredCellsPlacement[i].y + rgnYMax / 2.0 << endl;
    ++tetherTermIdx;
    plOut << "dummyTerm" << tetherTermIdx << " "
          << tetheredCellsPlacement[i].x - rgnXMax / 2.0 << " "
          << tetheredCellsPlacement[i].y + rgnYMax / 2.0 << endl;
    ++tetherTermIdx;
  }

  plOut.close();

  bool haveStorElts = false;
  for (i = 0; i < _storElt.size(); i++)
    if (_storElt[i]) {
      haveStorElts = true;
      break;
    }

  if (haveStorElts) {
    string lgtName = makeName(baseFileName, ".lgt");
    ofstream lgtOut(lgtName.c_str());
    lgtOut << "UCLA lgt 1.0" << endl;
    lgtOut << TimeStamp() << User() << Platform() << endl << endl;

    unsigned i = 0;
    for (i = 0; i < _hgWDims->getNumNodes(); i++) {
      if (!_storElt[i]) continue;

      lgtOut << setw(20) << _hgWDims->getNodeNameByIndex(i) << "    STOR "
             << endl;
    }
    lgtOut.close();
  }

  string auxName = makeName(baseFileName, ".aux");
  ofstream auxFile(auxName.c_str());

  auxFile << "RowBasedPlacement : ";

  const char* tmp = NULL;
  tmp = strrchr(baseFileName, pathDelim);
  if (tmp++) {
    auxFile << tmp << ".nodes " << tmp << ".nets " << tmp << ".wts " << tmp
            << ".pl " << tmp << ".scl ";
    if (haveStorElts) auxFile << tmp << ".lgt ";
    if (!_logicValues.empty()) auxFile << tmp << ".logic ";
  } else {
    auxFile << baseFileName << ".nodes " << baseFileName << ".nets "
            << baseFileName << ".wts " << baseFileName << ".pl " << baseFileName
            << ".scl ";
    if (haveStorElts) auxFile << baseFileName << ".lgt ";
    if (!_logicValues.empty()) auxFile << baseFileName << ".logic ";
  }
  auxFile << endl;

  auxFile.close();
}

void RBPlacement::parseLogic(const char* logicFile) {
  ifstream logic(logicFile);
  abkfatal2(logic, " Could not open ", logicFile);
  cout << " Reading " << logicFile << " ... " << endl;

  int lineno = 1;
  logic >> needcaseword("UCLA") >> needcaseword("logic") >>
      needcaseword("1.0") >> skiptoeol;
  logic >> eathash(lineno);

  while (logic.good()) {
    string nodeName;
    logic >> nodeName;
    unsigned nodeId = _hgWDims->getNodeByName(nodeName.c_str()).getIndex();
    const unsigned maxSize = 256;
    char logicData[maxSize];
    logic.getline(logicData, maxSize);
    _logicValues[nodeId] = logicData;

    logic >> eathash(lineno);
    ++lineno;
  }
}
