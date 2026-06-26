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

// created by Saurabh Adya, Mar 23 2000
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <vector>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/infolines.h"
#include "RBPlacement.h"

using std::cout;
using std::endl;
using std::ofstream;
using std::ostringstream;
using std::setw;
using std::vector;
using std::string;

#ifdef _MSC_VER
#ifndef rint
#define rint(a) floor((a) + 0.5)
#endif
#endif

// added by sadya to save the current placement as LEF/DEF pair.
// Macro cells,routing information are automatically generated

void RBPlacement::saveLEFDEF(const char* baseFileName,
                             bool markMacrosAsBlocks) {
  abkfatal3(strlen(baseFileName) < 60,
            "Base file name too long in saveLEFDEF: ", baseFileName, " \n");
  double rowHeight = coreRowsBegin()->getHeight();
  const unsigned units = 100;
  string lefFileName = string(baseFileName) + ".lef";

  cout << " Opening " << lefFileName << " ... " << endl;
  ofstream lef(lefFileName.c_str());
  abkfatal(lef, "Could not open Lef File for output");

  lef << endl << TimeStamp() << Platform() << User() << endl;

  vector<vector<double> > nodesInfo;
  vector<double> currNodeInfo;

  itHGFNodeLocal v;
  itHGFEdgeGlobal e;

  std::string cellName;
  // populate nodesInfo
  for (v = _hgWDims->nodesBegin(); v != _hgWDims->nodesEnd(); ++v) {
    cellName = _hgWDims->getNodeNameByIndex((*v)->getIndex());
    //      HGFNode& node = _hgWDims->getNodeByName(cellName.c_str());
    //      unsigned currNodeIdx=node.getIndex();
    unsigned currNodeIdx = (*v)->getIndex();
    currNodeInfo.push_back(_hgWDims->getNodeWidth(currNodeIdx));
    currNodeInfo.push_back(_hgWDims->getNodeHeight(currNodeIdx));
    unsigned numEdges = 0;
    for (e = (*v)->edgesBegin(), numEdges = 0; e != (*v)->edgesEnd();
         ++e, ++numEdges) {
    }

    currNodeInfo.push_back(numEdges);
    unsigned edgeOffset = 0;
    for (e = (*v)->edgesBegin(), edgeOffset = 0; e != (*v)->edgesEnd();
         ++e, ++edgeOffset) {
      Point currPinOffset = _hgWDims->edgesOnNodePinOffset(
          edgeOffset, currNodeIdx, _placement.getOrient(currNodeIdx));
      currNodeInfo.push_back(currPinOffset.x);
      currNodeInfo.push_back(currPinOffset.y);
    }
    nodesInfo.push_back(currNodeInfo);
    currNodeInfo.erase(currNodeInfo.begin(), currNodeInfo.end());
  }

  // sort the info according to sort_fn function
  std::sort(nodesInfo.begin(), nodesInfo.end(), sort_fn());

  vector<masterNodeInfo> masterNodes;

  bool newMasterNode = 0;
  unsigned masterNodeIdx = 0;
  char masterNodeName[20] = "MAS";
  unsigned i;
  for (i = 0; i < nodesInfo.size(); ++i) {
    if (i == nodesInfo.size() - 1)
      newMasterNode = 1;
    else if (nodesInfo[i] == nodesInfo[i + 1])
      newMasterNode = 0;
    else
      newMasterNode = 1;

    // one master node found
    if (newMasterNode == 1) {
      masterNodeInfo currMasterNode;
      currMasterNode.nodeInfo = nodesInfo[i];
      ostringstream masterName;
      masterName << masterNodeName << masterNodeIdx;
      strncpy(currMasterNode.name, masterName.str().c_str(),
              sizeof(currMasterNode.name) - 1);
      currMasterNode.name[sizeof(currMasterNode.name) - 1] = '\0';
      currMasterNode.index = masterNodeIdx;
      masterNodes.push_back(currMasterNode);
      ++masterNodeIdx;
    }
  }
  // start writing the lef file now
  lef << "NAMESCASESENSITIVE ON ;" << endl;
  lef << "UNITS" << endl;
  lef << "DATABASE MICRONS 100 ;" << endl;
  lef << "END UNITS" << endl << endl;

  // add routing info here
  lef << "LAYER poly" << endl;
  lef << "\tTYPE MASTERSLICE ;" << endl;
  lef << "END poly" << endl;

  lef << "LAYER cont" << endl;
  lef << "\tTYPE CUT ;" << endl;
  lef << "END cont" << endl;

  lef << "LAYER metal1" << endl;
  lef << "\t" << "TYPE ROUTING ;" << endl;
  lef << "\t" << "PITCH 1.7 ;" << endl;
  lef << "\t" << "WIDTH .8 ;" << endl;
  lef << "\t" << "SPACING .6 ;" << endl;
  lef << "\t" << "DIRECTION HORIZONTAL ;" << endl;
  lef << "END metal1" << endl << endl;

  lef << "LAYER via1" << endl;
  lef << "\t" << "TYPE CUT ;" << endl;
  lef << "END via1" << endl << endl;

  lef << "LAYER metal2" << endl;
  lef << "\t" << "TYPE ROUTING ;" << endl;
  lef << "\t" << "PITCH 2 ;" << endl;
  lef << "\t" << "WIDTH .8 ;" << endl;
  lef << "\t" << "SPACING .9 ;" << endl;
  lef << "\t" << "DIRECTION VERTICAL ;" << endl;
  lef << "END metal2" << endl << endl;

  lef << "LAYER via2" << endl;
  lef << "\t" << "TYPE CUT ;" << endl;
  lef << "END via2" << endl << endl;

  lef << "LAYER metal3" << endl;
  lef << "\t" << "TYPE ROUTING ;" << endl;
  lef << "\t" << "PITCH 1.7 ;" << endl;
  lef << "\t" << "WIDTH .8 ;" << endl;
  lef << "\t" << "SPACING .9 ;" << endl;
  lef << "\t" << "DIRECTION HORIZONTAL ;" << endl;
  lef << "END metal3" << endl << endl;

  lef << "LAYER OVERLAP" << endl;
  lef << "\tTYPE OVERLAP ;" << endl;
  lef << "END OVERLAP" << endl;

  // add via spec and SPACING info here to do
  lef << "VIA cont DEFAULT" << endl;
  lef << "\tLAYER poly ;" << endl;
  lef << "\t\tRECT -.6 -.6 .6 .6 ;" << endl;
  lef << "\tLAYER cont ;" << endl;
  lef << "\t\tRECT -.3 -.3 .3 .3 ;" << endl;
  lef << "\tLAYER metal1 ;" << endl;
  lef << "\t\tRECT -.4 -.4 .4 .4 ;" << endl;
  lef << "END cont" << endl << endl;

  lef << "VIA via1 DEFAULT" << endl;
  lef << "\tLAYER metal1 ;" << endl;
  lef << "\t\tRECT -.7 -.7 .7 .7 ;" << endl;
  lef << "\tLAYER via1 ;" << endl;
  lef << "\t\tRECT -.3 -.3 .3 .3 ;" << endl;
  lef << "\tLAYER metal2 ;" << endl;
  lef << "\t\tRECT -.4 -.4 .4 .4 ;" << endl;
  lef << "END via1" << endl << endl;

  lef << "VIA via2 DEFAULT" << endl;
  lef << "\tLAYER metal2 ;" << endl;
  lef << "\t\tRECT -.65 -.65 .65 .65 ;" << endl;
  lef << "\tLAYER via2 ;" << endl;
  lef << "\t\tRECT -.3 -.3 .3 .3 ;" << endl;
  lef << "\tLAYER metal3 ;" << endl;
  lef << "\t\tRECT -.4 -.4 .4 .4 ;" << endl;
  lef << "END via2" << endl << endl;

  lef << "SPACING" << endl;
  lef << "\tSAMENET metal1 metal1 .6 ;" << endl;
  lef << "\tSAMENET metal2 metal2 .9 ;" << endl;
  lef << "\tSAMENET metal3 metal3 .9 ;" << endl;
  lef << "\tSAMENET cont cont .8 ;" << endl;
  lef << "\tSAMENET cont via1 .65 ;" << endl;
  lef << "\tSAMENET via1 via1 1 ;" << endl;
  lef << "\tSAMENET via1 via2 .6 ;" << endl;
  lef << "\tSAMENET via2 via1 1 ;" << endl;
  lef << "END SPACING" << endl;

  // add site info here
  lef << "SITE core" << endl;
  lef << "\t" << "SIZE ";
  lef << coreRowsBegin()->getSiteWidth() / units << " BY "
      << coreRowsBegin()->getHeight() / units << " ;" << endl;
  lef << "\t" << "CLASS CORE ;" << endl;
  lef << "\t" << "SYMMETRY ";
  unsigned sym = coreRowsBegin()->getSymmetry();
  if (sym)
    lef << static_cast<const char*>(coreRowsBegin()->getSymmetry()) << " ;"
        << endl;
  else
    lef << "X ;" << endl;

  lef << "END core" << endl << endl;

  // to do what are correct offsets? what is correct symmetry
  // add macro info here
  for (i = 0; i < masterNodes.size(); ++i) {
    lef << "MACRO " << masterNodes[i].name << endl;
    if (masterNodes[i].nodeInfo[1] > rowHeight && markMacrosAsBlocks == true)
      lef << "\t" << "CLASS BLOCK ;" << endl;
    else
      lef << "\t" << "CLASS CORE ;" << endl;
    lef << "\t" << "SIZE " << masterNodes[i].nodeInfo[0] / units << " BY "
        << masterNodes[i].nodeInfo[1] / units << " ;" << endl;
    lef << "\t" << "ORIGIN 0 0 ;" << endl;
    lef << "\t" << "SYMMETRY X ;" << endl;
    lef << "\t" << "SITE core 0 0 N DO 1 BY 1 STEP 0 0 ;" << endl;  // to do
    // pin info
    double xOffset = 0.005 * masterNodes[i].nodeInfo[0];
    double yOffset = 0.005 * masterNodes[i].nodeInfo[1];
    for (unsigned j = 3; j < masterNodes[i].nodeInfo.size(); j += 2) {
      Point pinOffset(masterNodes[i].nodeInfo[j],
                      masterNodes[i].nodeInfo[j + 1]);
      char* pinName;
      pinName = new char[20];
      ostringstream pinNameStr;
      pinNameStr << "P" << ((j - 1) / 2);
      strncpy(pinName, pinNameStr.str().c_str(), 19);
      pinName[19] = '\0';
      masterNodes[i].pinNames.push_back(pinName);
      masterNodes[i].pinOffsets.push_back(pinOffset);

      lef << "\t" << "PIN " << pinName << endl;
      lef << "\t" << "DIRECTION INOUT ;" << endl;
      lef << "\t" << "USE SIGNAL ;" << endl;
      lef << "\t\t" << "PORT" << endl;
      lef << "\t\t\t" << "LAYER metal1 ;" << endl;
      lef << "\t\t\t\t" << "RECT ";
      lef << setw(4) << (masterNodes[i].nodeInfo[j] - xOffset) / units << " "
          << (masterNodes[i].nodeInfo[j + 1] - yOffset) / units << " "
          << (masterNodes[i].nodeInfo[j] + xOffset) / units << " "
          << (masterNodes[i].nodeInfo[j + 1] + yOffset) / units << " ;" << endl;
      lef << "\t\t" << "END" << endl;
      lef << "\t" << "END " << pinName << endl;
    }
    lef << "END " << masterNodes[i].name << endl << endl;
  }

  lef << "END LIBRARY" << endl << endl;
  lef.close();

  // now start def file
  const unsigned unitsDEF = 1;

  string defFileName = string(baseFileName) + ".def";
  cout << " Opening " << defFileName << " ... " << endl;
  ofstream def(defFileName.c_str());
  abkfatal(def, "Could not open Def File for output");

  // mapping between nodes and master nodes
  nodesMasInfo nodes(_hgWDims->getNumNodes());

  def << endl;
  def << TimeStamp();
  def << Platform();
  def << User() << endl;
  def << "DESIGN " << getOrigFileName() << " ;" << endl << endl;
  def << "UNITS DISTANCE MICRONS 100" << " ;" << endl << endl;

  BBox layoutBBox = getBBox();
  def << endl
      << "DIEAREA ( " << layoutBBox.xMin * unitsDEF << " "
      << layoutBBox.yMin * unitsDEF << " ) ( " << layoutBBox.xMax * unitsDEF
      << " " << layoutBBox.yMax * unitsDEF << " ) ;" << endl
      << endl;

  // to do SITES and TRACKS
  itRBPCoreRow itCore;
  for (itCore = coreRowsBegin(); itCore != coreRowsEnd(); ++itCore) {
    def << "SITE core ";
    def << itCore->getXStart() << " ";
    def << itCore->getYCoord() * unitsDEF << " ";
    def << itCore->getOrientation() << " ";
    def << "DO ";
    unsigned numSites = unsigned((itCore->getXEnd() - itCore->getXStart()) /
                                 itCore->getSiteWidth());
    def << numSites << " BY 1 STEP " << itCore->getSiteWidth() * unitsDEF;
    def << " " << itCore->getHeight() * unitsDEF << " ;" << endl;
  }

  /*TO DO*/
  // give 2 routes per sitewidth in both dimensions
  double siteWidth = coreRowsBegin()->getSiteWidth();
  int trackDist = int(siteWidth / 2.0);
  if (trackDist == 0) trackDist = 1;
  int numXTracks = int(layoutBBox.getWidth() * 1.5 / trackDist);
  int numYTracks = int(layoutBBox.getHeight() * 1.5 / trackDist);
  int xStart = int(layoutBBox.xMin);
  int yStart = int(layoutBBox.yMin);
  def << "TRACKS Y " << yStart << " DO " << numYTracks << " STEP " << trackDist
      << " LAYER metal1 metal2 metal3 ;" << endl;
  def << "TRACKS X " << xStart << " DO " << numXTracks << " STEP " << trackDist
      << " LAYER metal1 metal2 metal3 ;" << endl;

  def << "GCELLGRID Y " << yStart << " DO " << int(numYTracks / 10) << " STEP "
      << trackDist * 10 << " ;" << endl;
  def << "GCELLGRID X " << xStart << " DO " << int(numXTracks / 10) << " STEP "
      << trackDist * 10 << " ;" << endl;

  def << "VIAS 3 ;" << endl;
  def << "- cont + RECT metal1 ( -40 -40 ) ( 40 40 ) + RECT cont ( -30 -30 )"
      << endl
      << "( 30 30 ) + RECT poly ( -60 -60 ) ( 60 60 ) ;" << endl;
  def << "- via1 + RECT metal2 ( -40 -40 ) ( 40 40 ) + RECT via1 ( -30 -30 )"
      << endl
      << "( 30 30 ) + RECT metal1 ( -70 -70 ) ( 70 70 ) ;" << endl;
  def << "- via2 + RECT metal3 ( -40 -40 ) ( 40 40 ) + RECT via2 ( -30 -30 )"
      << endl
      << "( 30 30 ) + RECT metal2 ( -65 -65 ) ( 65 65 ) ;" << endl;
  def << "END VIAS" << endl << endl;

  // for each cell initialize the pinoffset counter to false
  //  make this counter true for every net connection. Used if multiple pins map
  //  to same location
  vector<vector<bool> > pinCounter;
  pinCounter.resize(getNumCells());

  def << "COMPONENTS " << getNumCells() << " ;" << endl;

  for (v = _hgWDims->nodesBegin(); v != _hgWDims->nodesEnd(); ++v) {
    unsigned cellIdx = (*v)->getIndex();
    cellName = _hgWDims->getNodeNameByIndex((*v)->getIndex());
    //      HGFNode& node = _hgWDims->getNodeByName(cellName.c_str());
    //      unsigned currNodeIdx=node.getIndex();
    unsigned currNodeIdx = (*v)->getIndex();
    currNodeInfo.push_back(_hgWDims->getNodeWidth(currNodeIdx));
    currNodeInfo.push_back(_hgWDims->getNodeHeight(currNodeIdx));
    unsigned numEdges = 0;
    for (e = (*v)->edgesBegin(), numEdges = 0; e != (*v)->edgesEnd();
         ++e, ++numEdges) {
    }

    currNodeInfo.push_back(numEdges);

    unsigned edgeOffset = 0;
    for (e = (*v)->edgesBegin(), edgeOffset = 0; e != (*v)->edgesEnd();
         ++e, ++edgeOffset) {
      Point currPinOffset = _hgWDims->edgesOnNodePinOffset(
          edgeOffset, currNodeIdx, _placement.getOrient(currNodeIdx));
      currNodeInfo.push_back(currPinOffset.x);
      currNodeInfo.push_back(currPinOffset.y);
    }

    // now match the currNodeInfo with matching master cell
    bool brokeFromLoop = 0;
    unsigned midx = 0;
    for (midx = 0; midx < masterNodes.size(); ++midx) {
      if (currNodeInfo == masterNodes[midx].nodeInfo) {
        pinCounter[cellIdx].resize(masterNodes[midx].pinOffsets.size());
        std::fill(pinCounter[cellIdx].begin(), pinCounter[cellIdx].end(), 0);
        brokeFromLoop = 1;
        break;
      }
    }
    abkfatal(brokeFromLoop,
             "Could Not match a node to master node in saveLEFDEF.cxx\n");
    def << "  - " << setw(10) << _hgWDims->getNodeNameByIndex((*v)->getIndex())
        << " " << setw(10) << masterNodes[midx].name << " ";
    nodes.putIndex(masterNodes[midx].index, cellIdx);

    if (isFixed(cellIdx)) {
      def << " + " << setw(11) << "FIXED ( "
          << static_cast<long int>(_placement[cellIdx].x * unitsDEF) << " "
          << static_cast<long int>(_placement[cellIdx].y * unitsDEF)
          << " ) "
          //<< getOrient(cellIdx) << " ";
          << " N " << " ";
    } else if (isInSubRow(cellIdx)) {
      def << " + " << setw(11) << "PLACED ( "
          << static_cast<long int>(_placement[cellIdx].x * unitsDEF) << " "
          << static_cast<long int>(_placement[cellIdx].y * unitsDEF) << " ) "
          << getOrient(cellIdx) << " ";
    } else {
      def << " + " << setw(11) << "UNPLACED ( "
          << static_cast<long int>(_placement[cellIdx].x * unitsDEF) << " "
          << static_cast<long int>(_placement[cellIdx].y * unitsDEF) << " ) "
          << getOrient(cellIdx) << " ";
    }
    def << ";" << endl;
    currNodeInfo.erase(currNodeInfo.begin(), currNodeInfo.end());
  }
  def << "END COMPONENTS" << endl << endl;

  unsigned numEdges = 0;
  for (e = _hgWDims->edgesBegin(); e != _hgWDims->edgesEnd(); ++e) ++numEdges;

  bool brokeFromLoop = 0;
  def << "NETS " << numEdges << " ;" << endl;
  for (e = _hgWDims->edgesBegin(); e != _hgWDims->edgesEnd(); ++e) {
    unsigned netId = (*e)->getIndex();
    ostringstream netName;
    netName << "net" << netId;
    def << "- " << netName.str() << " ";
    unsigned nodeOffset = 0;
    for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); ++v, ++nodeOffset) {
      unsigned cellId = (*v)->getIndex();
      Point pinOffset = _hgWDims->nodesOnEdgePinOffset(
          nodeOffset, netId, _placement.getOrient(cellId));
      unsigned masIndex = nodes.getMasterIndex(cellId);
      unsigned midx = 0;
      brokeFromLoop = 0;
      for (midx = 0; midx < masterNodes.size(); ++midx) {
        if (masIndex == masterNodes[midx].index) {
          brokeFromLoop = 1;
          break;
        }
      }

      abkfatal(brokeFromLoop, "Could not find masterNode in saveLEFDEF\n");
      unsigned pidx = 0;
      brokeFromLoop = 0;
      // search for the pin in masterNode
      for (pidx = 0; pidx < masterNodes[midx].pinOffsets.size(); ++pidx) {
        if (pinOffset == masterNodes[midx].pinOffsets[pidx]) {
          if (pinCounter[cellId][pidx] == 1) continue;
          pinCounter[cellId][pidx] = 1;
          brokeFromLoop = 1;
          break;
        }
      }
      abkfatal(brokeFromLoop,
               "Could not find reqd Pin in masterNode in saveLEFDEF\n");

      def << endl
          << "( " << _hgWDims->getNodeNameByIndex((*v)->getIndex()) << " "
          << masterNodes[midx].pinNames[pidx] << " ) ";
    }
    def << " ;" << endl;
  }
  def << "END NETS " << endl << endl;
  def << "END DESIGN" << endl;
  def.close();

  // now start def file
  string auxFileName = string(baseFileName) + ".aux";
  cout << " Opening " << auxFileName << " ... " << endl;
  ofstream auxF(auxFileName.c_str());
  abkfatal(auxF, "Could not open Aux File for output");
  auxF << " LEFDEF : " << " " << lefFileName << " " << defFileName << endl;
}
