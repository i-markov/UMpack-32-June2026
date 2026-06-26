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

// Created: Nov 30, 1999 by Andrew Caldwell

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>
#include <fstream>
#include <string>
#include <strings.h>
#include <unordered_map>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abk_hash_common.h"
#include "ABKCommon/abkio.h"
#include "ABKCommon/abktimer.h"
#include "PlaceEvals/placeEvals.h"
#include "Placement/xyDraw.h"
#include "RBPlace/RBPlacement.h"
#include "RBPlace/slicingTree.h"
#include "Stats/trivStats.h"
#include "capoBin.h"
#include "capoPlacer.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;

namespace {

string extractQuotedAttribute(const char* token) {
  abkfatal(token != NULL, "null token while parsing HCL");
  const char* begin = strchr(token, '\"');
  abkfatal(begin != NULL, "expected quoted attribute value in HCL");
  ++begin;
  const char* end = strchr(begin, '\"');
  abkfatal(end != NULL, "unterminated quoted attribute value in HCL");
  return string(begin, end - begin);
}

}  // namespace

CapoPlacer::CapoPlacer(RBPlacement& rbplace, const char* hclFileName,
                       const char* binPLFileName, const CapoParameters& params)
    : _params(params),
      _rbplace(rbplace),
      _coreBBox(/*BBoxFromRBPlace*/),
      _placement(static_cast<PlacementWOrient>(rbplace)),
      _oraclePlacement(
          (_params.splitterParams.WSOracle || _params.splitterParams.eco)
              ? static_cast<PlacementWOrient>(rbplace)
              : PlacementWOrient(0)),
      _ispd04CongMap(NULL),
      _ispd04CongMapPlacement(NULL),
      _placed(rbplace.getFixed()),
      _totalLayerNum(0),
      _mainLayerNum(0),
      _feedbackLayerNum(0),
      nextBinNum(0),
      _cellToBinMap(rbplace.getHGraph().getNumNodes(),
                    static_cast<const CapoBin*>(NULL)),
      _netUpperBounds(rbplace.getHGraph().getNumEdges()),
      _netLowerBounds(rbplace.getHGraph().getNumEdges()),
      _nodeSeen(rbplace.getHGraph().getNumNodes()),
      _edgeSeen(rbplace.getHGraph().getNumEdges()),
      _nodesNeedingFP(rbplace.getHGraph().getNumNodes(), false),
      partProbEdgesVisited(rbplace.getHGraph().getNumEdges()),
      thetoWeights(params.wtCut ? rbplace.getHGraph().getNumEdges() : 0,
                   vector<double>(3, 0.)),
      _numOrigNets(100, 0),
      _numEssentialNets(100, 0),
      _origNetPins(100, 0),
      _essentialNetPins(100, 0),
      _numProblemsOfSize(100, 0),
      _netHasBeenCut(rbplace.getHGraph().getNumEdges(), false),
      _numNotSolved(0),
      _numSmallPartProbs(0),
      //        _slicingTree(rbplace.getBBox(true)),
      constraints(_rbplace.constraints),
      regionUtilization(_rbplace.regionUtilization),
      groupRegionConstr(_rbplace.groupRegionConstr),
      cellToGrpMapping(_rbplace.cellToGrpMapping),
      //        _altCellNames(0)
      _fastStPt(NULL),
      _fastStParent(NULL) {
  // compilerCheck();

  Timer capoTimer;

  if (_params.verb.getForActions() > 0)
    cerr << endl
         << "Launching UCLA's CapoPlacer " << " from a Floorplanning" << endl;

  setupAndCheck();

  readBinsFile(binPLFileName, hclFileName);
  _mainLayerNum = static_cast<unsigned>(
      ceil(sqrt(static_cast<double>(_curLayout->size()))));
  _totalLayerNum = _mainLayerNum;
  cout << "Starting at layer: " << _mainLayerNum << "." << _feedbackLayerNum
       << endl;

  if (_curLayout == &_layout[0]) {
    _layout[1] = vector<CapoBin*>(_layout[0].size());
    std::copy(_layout[0].begin(), _layout[0].end(), _layout[1].begin());
  } else if (_curLayout == &_layout[1]) {
    _layout[0] = vector<CapoBin*>(_layout[1].size());
    std::copy(_layout[1].begin(), _layout[1].end(), _layout[0].begin());
  } else
    abkfatal(0, "_curLayout does not equal either layout");

  unsigned maxBins = _params.stopAtBins ? _params.stopAtBins : UINT_MAX;

  bool haveUnplacedBins = true;
  bool doVCuts = true;
  while (haveUnplacedBins && _curLayout->size() < maxBins) {
    AllowedPartDir partDir;
    if (_mainLayerNum < 3)
      partDir = HAndV;
    else if (_mainLayerNum == 3 || doVCuts) {
      partDir = VOnly;
      doVCuts = false;  // next time will be horizontal
    } else {
      partDir = HOnly;
      doVCuts = true;
    }

    // we have no idea at this point which 'layout' will be choosen..
    // and which will be cleared..so they are made to be both the same.

    bool dummy;
    haveUnplacedBins = doOneLayer(partDir, dummy);
  }

  if (_params.replaceSmallBins == AtTheEnd) {
    Timer replaceTimer;
    cout << " Replacing Small Bins" << endl;
    replaceSmallBins();
    replaceTimer.stop();
    cout << " Replacing Small Bins took: " << replaceTimer.getUserTime()
         << "sec" << endl;
  }

  capoTimer.stop();

  if (_params.verb.getForMajStats() > 0)
    cout << "   CapoPlacer took : " << capoTimer << endl;

  rbplace.resetPlacementWOri(_placement);
}

struct eqstr {
  bool operator()(const char* s1, const char* s2) const {
    return strcmp(s1, s2) == 0;
  }
};

void CapoPlacer::readBinsFile(const char* binPlName, const char* hclName) {
  ifstream binsFile(binPlName);
  abkfatal2(binsFile, "Could not open bin file ", binPlName);

  CapoBin* topLvlBin = (*_curLayout)[0];
  _curLayout->clear();

  std::unordered_map<string, unsigned> nameToBinIdMap;

  binsFile >> needword("BoundingBox") >> skiptoeol;
  string buffer;
  binsFile >> buffer;

  while (strcasecmp(buffer.c_str(), "Name")) {
    binsFile >> buffer;
  }

  while (!binsFile.eof()) {
    abkfatal(!strcasecmp(buffer.c_str(), "Name"),
             "expected 'Name' at the start of the line");

    binsFile >> buffer;
    Point origin;
    binsFile >> needword("location") >> needword("(") >> origin.x >>
        needword(",") >> origin.y >> needword(")");
    BBox binBBox(origin);

    binsFile >> needword("ori") >> needword("N");
    binsFile >> needword("width") >> origin.x >> needword("height") >>
        origin.y >> skiptoeol;

    binBBox += (Point(binBBox.xMin + origin.x, binBBox.yMin + origin.y));

    vector<unsigned> emptyCellIds;
    vector<CapoBin*> emptyNeighbors;

    _curLayout->push_back(new CapoBin(
        emptyCellIds, topLvlBin->_coreRows, binBBox, emptyNeighbors,
        emptyNeighbors, emptyNeighbors, emptyNeighbors, *this,
        buffer.c_str()));

    nameToBinIdMap[_curLayout->back()->getName()] = _curLayout->size() - 1;
    binsFile >> buffer;
  }
  binsFile.close();

  ifstream hclFile(hclName);
  abkfatal2(hclFile, "Could not open hcl file ", hclName);

  hclFile >> needword("<hcl>") >> skiptoeol;
  hclFile >> buffer;
  while (buffer != "<cluster") {
    hclFile >> buffer;
  }

  vector<vector<unsigned> > cellIds(_curLayout->size());

  while (buffer != "</body>") {
    hclFile >> buffer;
    const string nodeName = extractQuotedAttribute(buffer.c_str());

    hclFile >> buffer;
    const string binName = extractQuotedAttribute(buffer.c_str());

    if (!strcasecmp(nodeName.c_str(), binName.c_str())) {  // the root
      hclFile >> skiptoeol;
      hclFile >> needcaseword("<cluster");
      continue;
    }

    unsigned binId = nameToBinIdMap[binName];
    unsigned nId = getNetlistHGraph().getNodeByName(nodeName.c_str()).getIndex();
    cellIds[binId].push_back(nId);

    hclFile >> skiptoeol;
    hclFile >> buffer;
  }
  hclFile.close();

  double totalCap = 0;
  for (unsigned b = 0; b < _curLayout->size(); b++) {
    (*_curLayout)[b]->resetCellIds(cellIds[b]);

    cout << *(*_curLayout)[b] << endl;
    totalCap += (*_curLayout)[b]->getCapacity();
  }

  cout << "Total Bin Cap: " << totalCap << endl;
  cout << "Orig Bin  Cap: " << topLvlBin->getCapacity() << endl;

  delete topLvlBin;
}
