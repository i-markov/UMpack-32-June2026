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

#include <cctype>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkio.h"
#include "Ctainers/bitBoard.h"
#include "Geoms/bbox.h"
#include "hgWDims.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::string;
using std::vector;

void HGraphWDimensions::readNodes(const char* nodesFileName) {
  // the base-class parser just populates the cell names, etc.
  // this allows us to do lookups by cell name when populating
  // things like height, width, etc.
  HGraphFixed::readNodes(nodesFileName);

  _nodeSymmetries = vector<Symmetry>(getNumNodes());
  _nodeWidths = vector<float>(getNumNodes(), 1.0);
  _nodeHeights = vector<float>(getNumNodes(), 1.0);

  ifstream nodes(nodesFileName);
  abkfatal2(nodes, " Could not open ", nodesFileName);

  int lineno = 1;
  nodes >> needcaseword("UCLA") >> needcaseword("nodes") >> needword("1.0") >>
      skiptoeol >> eathash(lineno);

  unsigned junkUnsigned;
  nodes >> needcaseword("NumNodes") >> needword(":") >> my_isnumber(lineno) >>
      junkUnsigned >> skiptoeol;
  lineno++;

  nodes >> needcaseword("NumTerminals") >> needword(":") >>
      my_isnumber(lineno) >> junkUnsigned >> skiptoeol;
  lineno++;
  // this is a huge security hole:
  //  char tmpStr[1024];
  std::string tmpStr;
  nodes >> tmpStr;

  for (unsigned nodeCtr = 0; nodeCtr < _nodes.size(); nodeCtr++) {
    abkfatal(!tmpStr.empty(), "not enough nodeNames in names file");
    HGFNode& node = getNodeByName(tmpStr.c_str());
    unsigned nId = node.getIndex();

    nodes >> eatblank;
    char nxtChar;

    while (!nodes.eof() && nodes.peek() != '\n' && nodes.peek() != '\r') {
      nxtChar = nodes.peek();
      if (nxtChar == 't' || nxtChar == 'T') {
        nodes >> tmpStr >> eatblank;
        continue;  // base class processes this...
      }
      if (nxtChar == ':' || nxtChar == '/') {
        char c;
        //              nodes>>needword(":")>>eatblank;
        nodes >> c >> eatblank;
        // tmpStr[0] = 0;	//makes strcat start from the beginning
        // char symStr[10];
        std::string symStr;
        while (nodes.peek() == 'X' || nodes.peek() == 'Y' ||
               nodes.peek() == 'R' || nodes.peek() == 'x' ||
               nodes.peek() == 'y' || nodes.peek() == 'r') {
          nodes >> symStr >> eatblank;
          tmpStr += symStr;
        }
        _nodeSymmetries[nId] = Symmetry(tmpStr.c_str());
      } else if (isdigit(nxtChar))  // node dimensions
      {
        nodes >> _nodeWidths[nId] >> _nodeHeights[nId] >> eatblank;
        _nodeHeights[nId] *= _param.yScale;
      } else
        abkfatal2(0,
                  "Error in .nodes (missing 'Terminals'?), unexpected string: ",
                  tmpStr);
    }
    lineno++;

    if (!nodes.eof()) nodes >> tmpStr;
  }

  nodes.close();
}

void HGraphWDimensions::readNets(const char* netsFileName) {
  ifstream nets(netsFileName);
  abkfatal2(nets, " Could not open ", netsFileName);
  cout << " Reading " << netsFileName << " ... " << endl;

  int lineno = 1;
  nets >> needcaseword("UCLA") >> needcaseword("nets") >> needcaseword("1.0") >>
      skiptoeol;
  nets >> eathash(lineno);

  unsigned expectedNumPins = _numPins = 0;
  string cName;

  nets >> cName;
  if (cName == "NumNets") {
    unsigned expectedNumNets = 0;
    nets >> needword(":") >> my_isnumber(lineno) >> expectedNumNets >>
        skiptoeol;
    _edges.reserve(expectedNumNets);
    nets >> needcaseword("NumPins") >> needword(":") >> my_isnumber(lineno) >>
        expectedNumPins >> skiptoeol;
  } else if (cName == "NumPins")
    nets >> needword(":") >> my_isnumber(lineno) >> expectedNumPins >>
        skiptoeol;

#ifdef SIGNAL_DIRECTIONS
  _srcsWPins.reserve(1 + (expectedNumPins / 3));
  _snksWPins.reserve(1 + (expectedNumPins / 3));
  _srcSnksWPins.reserve(1 + (expectedNumPins / 3));
#else
  _srcSnksWPins.reserve(expectedNumPins);
#endif

  // bbox of a set of pins on a node
  std::unordered_map<unsigned, BBox, std::hash<unsigned>, std::equal_to<unsigned> > pinBBoxes;
  // direction of a set of pins on a node (B, I, O)
  std::unordered_map<unsigned, char, std::hash<unsigned>, std::equal_to<unsigned> > pinBBoxDirs;
  vector<unsigned> nodeOrder;

  unsigned netDegree;
  char tmpChar;
  string tmpNetName;

  while (_numPins < expectedNumPins) {
    nets >> eathash(lineno) >> needcaseword("NetDegree") >> needword(":") >>
        my_isnumber(lineno) >> netDegree >> eatblank;
    char c = nets.peek();
    if (c != '\n' && c != '\r')
      nets >> tmpNetName >> skiptoeol;
    else
      tmpNetName.clear();

    HGFEdge* edge;
    if (netDegree > 1) {
      edge = addEdge();
      if (!tmpNetName.empty()) {
        // char * newName=new char[strlen(tmpNetName)+1];
        // strcpy(newName, tmpNetName);
        getNetNames().resize(getNetNames().size() + 1);
        getNetNames().back() = tmpNetName;
        getNetNamesMap()[getNetNames().back()] = edge->getIndex();
      } else
        getNetNames().resize(getNetNames().size() + 1);
    } else {
      // cerr << " Net of degree 1 found ";
      //  MRG - to run golem3
      // if (strlen(tmpNetName)) cerr << " named " << tmpNetName;
      // abkfatal(0,"This parser does not allow nets of degree < 2");
      // makingEdge = false;
      edge = NULL;
    }

    char pinDir;
    Point offset;

    /* <aaronnn>
     * commented out because we want to consider repeated pins
     */
    // nodesSeen.clear();

    pinBBoxes.clear();
    pinBBoxDirs.clear();
    nodeOrder.clear();

    for (unsigned n = 0; n < netDegree; n++) {
      // MRG - to skip degree 1 nets
      if (!edge) {
        nets >> skiptoeol >> eathash(lineno) >> cName >> eatblank >> skiptoeol;
        _numPins++;
        continue;
      }

      nets >> eathash(lineno) >> cName >> eatblank;

      abkfatal(cName != "NetDegree",
               "did not find expected number of node names");
      abkfatal2(haveSuchNode(cName.c_str()), "no such node: ", cName.c_str());
      HGFNode& node = getNodeByName(cName.c_str());
      unsigned nId = node.getIndex();

      pinDir = 'B';
      offset.x = _nodeWidths[nId] / 2.0;
      offset.y = _nodeHeights[nId] / 2.0;

      tmpChar = nets.peek();

      while (tmpChar != '\n' && tmpChar != '\r') {
        nets >> tmpChar;

        tmpChar = toupper(tmpChar);
        if (tmpChar == 'B' || tmpChar == 'I' || tmpChar == 'O') {
          pinDir = tmpChar;
          if (_param.makeAllSrcSnk && (tmpChar == 'O' || tmpChar == 'I'))
            pinDir = 'B';
        } else if (tmpChar == ':') {
          nets >> offset.x >> offset.y;
          offset.y *= _param.yScale;
          offset.x += _nodeWidths[nId] / 2.0;
          offset.y += _nodeHeights[nId] / 2.0;
        } else if (tmpChar == '/') {
          std::string pinName;
          nets >> pinName;
          nets >> skiptoeol;
        } else
          abkfatal2(0, "unknown char in nets format: ", tmpChar);

        nets >> eatblank;
        tmpChar = nets.peek();
      }

      if (pinDir != 'B' && pinDir != 'I' && pinDir != 'O')
        abkfatal3(0, pinDir, ": unexpected symbol in line ", lineno);

      // expand the bbox of repeated pins on a node
      // (trust the bbox default constructor for proper init)
      pinBBoxes[nId] += offset;

      if (pinBBoxDirs.find(nId) == pinBBoxDirs.end()) {
        // first time seeing a pin for this node
        nodeOrder.push_back(nId);
        pinBBoxDirs[nId] = pinDir;
      } else {
        if (pinBBoxDirs[nId] != 'B') {
          if (pinBBoxDirs[nId] != pinDir) pinBBoxDirs[nId] = 'B';
        }
      }

      _numPins++;
      nets >> skiptoeol;
      nets >> eathash(lineno);
      lineno++;
    }

    for (unsigned i = 0; i < nodeOrder.size(); ++i) {
      unsigned nId = nodeOrder[i];
      HGFNode& node = getNodeByIdx(nId);
      char pinDir = pinBBoxDirs[nId];
      const Point& offset = pinBBoxes[nId].getGeomCenter();
      PinPoint pinPoint(static_cast<float>(offset.x),
                        static_cast<float>(offset.y),
                        static_cast<float>(pinBBoxes[nId].getWidth() / 2.0),
                        static_cast<float>(pinBBoxes[nId].getHeight() / 2.0));

      switch (pinDir) {
        case ('I'):
#ifdef SIGNAL_DIRECTIONS
          addSnkWPin(node, *edge, pinPoint);
          break;
#endif
        case ('O'):
#ifdef SIGNAL_DIRECTIONS
          addSrcWPin(node, *edge, pinPoint);
          break;
#endif
        case ('B'):
          addSrcSnkWPin(node, *edge, pinPoint);
          break;
      };
    }
  }
  nets.close();
}

// determine which nodes are soft from a .blocks file
void HGraphWDimensions::markSelectedAsSoftFromFile(const string& fileName) {
  ifstream infile(fileName.c_str());
  abkfatal2(infile.good(), " Could not open ", fileName);
  cout << " Reading " << fileName << " ... " << endl;

  std::string line;
  getline(infile, line);
  while (infile.good()) {
    std::string::size_type loc = line.find("softrectangular", 0);
    if (loc != std::string::npos) {
      std::stringstream ss;
      ss << line;
      std::string nodeName, dummy;
      ss >> nodeName >> dummy;

      abkfatal2(haveSuchNode(nodeName.c_str()), nodeName.c_str(),
                " does not exist.");

      const HGFNode& currNode = getNodeByName(nodeName.c_str());
      const unsigned nodeId = currNode.getIndex();

      double nodeArea, nodeMinAR, nodeMaxAR;
      ss >> nodeArea >> nodeMinAR >> nodeMaxAR;

      updateNodeSoftInfo(nodeId, nodeMinAR, nodeMaxAR);
      //          printf("[%d] name: %s area: %.2f minAR: %.2f maxAR: %.2f\n",
      //                 nodeId, nodeName.c_str(), nodeArea, nodeMinAR,
      //                 nodeMaxAR); //-----
    }
    getline(infile, line);
  }
}
// --------------------------------------------------------
