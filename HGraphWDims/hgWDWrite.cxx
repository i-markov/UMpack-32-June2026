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

#include <algorithm>
#include <cstring>
#include <deque>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkio.h"
#include "ABKCommon/abkrand.h"
#include "ABKCommon/infolines.h"
#include "ABKCommon/pathDelims.h"
#include "hgWDims.h"

using std::cout;
using std::deque;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::setw;
using std::string;
using std::vector;

static string fileNameWithSuffix(const char* baseFileName, const char* suffix) {
  return string(baseFileName) + suffix;
}

bool HGraphWDimensions::saveAsNodesNetsWts(const char* baseFileName,
                                           bool fixMacros, bool saveNodeWts,
                                           bool saveNetWts) const {
  const string netsFileName = fileNameWithSuffix(baseFileName, ".nets");
  const string nodesFileName = fileNameWithSuffix(baseFileName, ".nodes");
  const string wtsFileName = fileNameWithSuffix(baseFileName, ".wts");
  const string auxFileName = fileNameWithSuffix(baseFileName, ".aux");

  ofstream nets(netsFileName.c_str());

  if (!nets.is_open()) return false;

  ofstream nodes(nodesFileName.c_str());

  if (!nodes.is_open()) return false;

  ofstream wts;
  if (saveNodeWts || saveNetWts) {
    wts.open(wtsFileName.c_str());
    if (!wts.is_open()) return false;
  }

  ofstream aux(auxFileName.c_str());

  if (!aux.is_open()) return false;

  itHGFNodeGlobal nIt;

  // get number of macros
  unsigned numMacros = 0;
  if (fixMacros) {
    for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
      HGFNode& node = **nIt;
      unsigned nodeId = node.getIndex();
      if (!isTerminal(nodeId)) {
        if (_isMacro[nodeId]) ++numMacros;
      }
    }
  }
  // write the nodes file
  nodes << "UCLA nodes 1.0" << endl;
  nodes << TimeStamp() << User() << Platform() << endl;

  nodes << "NumNodes :     " << setw(10) << getNumNodes() << endl;
  nodes << "NumTerminals : " << setw(10) << getNumTerminals() + numMacros
        << endl;

  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    Symmetry sym = _nodeSymmetries[nodeId];

    if (!getNodeNameByIndex(node.getIndex()).empty())
      nodes << setw(10) << getNodeNameByIndex(node.getIndex()) << "  ";
    else {
      if (isTerminal(nodeId))
        nodes << " p" << nodeId + 1 << " ";
      else
        nodes << " a" << nodeId - getNumTerminals() << " ";
    }

    if (getNodeWidth(nodeId) != 0 || getNodeHeight(nodeId) != 0)
      nodes << setw(10) << getNodeWidth(nodeId) << " " << setw(10)
            << getNodeHeight(nodeId) << " ";
    if (isTerminal(nodeId))
      nodes << "terminal ";
    else if (fixMacros && _isMacro[nodeId])
      nodes << "terminal ";

    if (sym.rot90 || sym.y || sym.x) nodes << " : " << sym;
    nodes << endl;
  }
  nodes.close();

  unsigned num1PinEdges = 0, num0PinEdges = 0;
  itHGFEdgeGlobal e;
  for (e = edgesBegin(); e != edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    if (edge.getDegree() > 1) continue;
    if (edge.getDegree() == 1)
      num1PinEdges++;
    else
      num0PinEdges++;
  }

  // write the nets file
  nets << "UCLA nets 1.0" << endl;
  nets << TimeStamp() << User() << Platform() << endl;

  nets << "NumNets : " << getNumEdges() - num1PinEdges << endl;
  nets << "NumPins : " << getNumPins() - num1PinEdges << endl;

  for (e = edgesBegin(); e != edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    if (edge.getDegree() < 2) continue;
    unsigned edgeId = edge.getIndex();
    std::string netName = getNetNameByIndex(edgeId).c_str();

    nets << "NetDegree : " << (*e)->getDegree() << " ";
    if (!netName.empty())
      nets << setw(10) << netName;
    else if (fabs(edge.getWeight() - 1.0) > 1e-10)
      nets << " net" << edgeId;
    else
      nets << " net" << edgeId;
    nets << endl;

    unsigned nodeOffset = 0;

    itHGFNodeLocal v;
#ifdef SIGNAL_DIRECTIONS
    for (v = (*e)->srcsBegin(); v != (*e)->srcsEnd(); v++, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " O ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;
      //	    if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;
      else
        nets << endl;
    }
    for (v = (*e)->srcSnksBegin(); v != (*e)->srcSnksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " B ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;
      //          if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;
      else
        nets << endl;
    }
    for (v = (*e)->snksBegin(); v != (*e)->snksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " I ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;

      //	    if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;

      else
        nets << endl;
    }
#else
    for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " B ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;

      //	    if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;

      else
        nets << endl;
    }
#endif
  }
  nets.close();

  // write the wts file

  if (saveNodeWts || saveNetWts) {
    wts << "UCLA wts 1.0" << endl;
    wts << TimeStamp() << User() << Platform() << endl;

    unsigned n;
    if (saveNodeWts && getNumWeights() >= 1) {
      wts << "# Node weights" << endl;
      for (n = 0; n != getNumNodes(); n++) {
        HGFNode& node = *_nodes[n];
        wts << setw(10) << getNodeNameByIndex(node.getIndex()) << "  ";

        for (unsigned i = 0; i != getNumWeights(); ++i)
          wts << " " << setw(6) << getWeight(node.getIndex(), i);

        wts << endl;
      }
    }

    if (saveNetWts) {
      wts << "# Net weights" << endl;
      for (n = 0; n != getNumEdges(); ++n) {
        const HGFEdge& net = getEdgeByIdx(n);
        if (fabs(net.getWeight() - 1.0) < 1e-10) continue;
        wts << " ";
        std::string netName = getNetNameByIndex(n).c_str();
        if (!netName.empty())
          wts << setw(10) << netName;
        else
          wts << "net" << n;
        wts << " " << setw(10) << net.getWeight() << endl;
      }
    }

    wts.close();
  }

  aux << "HGraphWPins : ";

  const char* tmp = NULL;
  tmp = strrchr(baseFileName, pathDelim);
  if (tmp++) {
    aux << tmp << ".nets " << tmp << ".nodes ";
    if (saveNodeWts || saveNetWts) aux << tmp << ".wts";
    aux << endl;
  } else {
    aux << baseFileName << ".nets " << baseFileName << ".nodes ";
    if (saveNodeWts || saveNetWts) aux << baseFileName << ".wts";
    aux << endl;
  }
  aux.close();

  return true;
}

void HGraphWDimensions::saveAsNodesNetsWtsWCongInfo(
    const char* baseFileName, std::vector<double>& nodesXCongestion,
    std::vector<double>& nodesYCongestion, double whitespace,
    float rowHeight) {
  unsigned numExtraNodes = 1;
  const string netsFileName = fileNameWithSuffix(baseFileName, ".nets");
  const string nodesFileName = fileNameWithSuffix(baseFileName, ".nodes");
  const string wtsFileName = fileNameWithSuffix(baseFileName, ".wts");
  const string auxFileName = fileNameWithSuffix(baseFileName, ".aux");

  ofstream nets(netsFileName.c_str());

  ofstream nodes(nodesFileName.c_str());

  ofstream wts(wtsFileName.c_str());

  ofstream aux(auxFileName.c_str());

  itHGFNodeGlobal nIt;

  unsigned numNodes = getNumNodes();
  vector<ctrSortDoublePair> sortedNodesByCongestion(numNodes);
  for (unsigned i = 0; i < numNodes; ++i) {
    sortedNodesByCongestion[i].index = double(i);
    sortedNodesByCongestion[i].value =
        nodesXCongestion[i] + nodesYCongestion[i];
  }
  sort(sortedNodesByCongestion.begin(), sortedNodesByCongestion.end(),
       doublePair_less_mag());
  //    float dummyNodeWidth = floor(getAvgNodeWidth()/numExtraNodes);
  float dummyNodeWidth = floor(2.0 * getAvgNodeWidth());

  if (dummyNodeWidth == 0) dummyNodeWidth = 1;
  float dummyNodeHeight = rowHeight;
  float dummyNodeArea = dummyNodeHeight * dummyNodeWidth;

  vector<unsigned> congestedNodes;
  int idx = numNodes - 1;
  double currWS = whitespace;
  while (currWS > 0 && idx >= 0) {
    if (!isTerminal(unsigned(sortedNodesByCongestion[idx].index))) {
      congestedNodes.push_back(unsigned(sortedNodesByCongestion[idx].index));
      currWS -= dummyNodeArea * numExtraNodes;
    }
    --idx;
  }

  // TODO
  // write the nodes file
  nodes << "UCLA nodes 1.0" << endl;
  nodes << TimeStamp() << User() << Platform() << endl;

  nodes << "NumNodes :     " << setw(10)
        << getNumNodes() + numExtraNodes * congestedNodes.size() << endl;
  nodes << "NumTerminals : " << setw(10) << getNumTerminals() << endl;

  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    Symmetry sym = _nodeSymmetries[nodeId];

    if (!getNodeNameByIndex(node.getIndex()).empty())
      nodes << setw(10) << getNodeNameByIndex(node.getIndex()) << "  ";
    else {
      if (isTerminal(nodeId))
        nodes << " p" << nodeId + 1 << " ";
      else
        nodes << " a" << nodeId - getNumTerminals() << " ";
    }

    if (getNodeWidth(nodeId) != 0 || getNodeHeight(nodeId) != 0)
      nodes << setw(10) << getNodeWidth(nodeId) << " " << setw(10)
            << getNodeHeight(nodeId) << " ";
    if (isTerminal(nodeId)) nodes << "terminal ";

    if (sym.rot90 || sym.y || sym.x) nodes << " : " << sym;
    nodes << endl;
  }

  for (unsigned i = 0; i < congestedNodes.size(); ++i) {
    unsigned nodeId = congestedNodes[i];
    HGFNode& node = getNodeByIdx(nodeId);

    for (unsigned k = 0; k < numExtraNodes; ++k) {
      if (!getNodeNameByIndex(node.getIndex()).empty())
        nodes << "dummy_" << k << "_" << getNodeNameByIndex(node.getIndex())
              << "  ";
      else
        nodes << "dummy_a" << k << "_" << nodeId - getNumTerminals() << " ";

      nodes << setw(10) << dummyNodeWidth << " " << setw(10) << dummyNodeHeight
            << " " << endl;
    }
  }
  nodes.close();

  unsigned num1PinEdges = 0, num0PinEdges = 0;
  itHGFEdgeGlobal e;
  for (e = edgesBegin(); e != edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    if (edge.getDegree() > 1) continue;
    if (edge.getDegree() == 1)
      num1PinEdges++;
    else
      num0PinEdges++;
  }

  // write the nets file
  nets << "UCLA nets 1.0" << endl;
  nets << TimeStamp() << User() << Platform() << endl;

  nets << "NumNets : "
       << getNumEdges() - num1PinEdges + numExtraNodes * congestedNodes.size()
       << endl;
  nets << "NumPins : "
       << getNumPins() - num1PinEdges +
              2 * (numExtraNodes * congestedNodes.size())
       << endl;

  for (e = edgesBegin(); e != edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    if (edge.getDegree() < 2) continue;
    unsigned edgeId = edge.getIndex();
    const char* netName = getNetNameByIndex(edgeId).c_str();

    nets << "NetDegree : " << (*e)->getDegree() << " ";
    if (netName)
      nets << setw(10) << netName;
    else if (fabs(edge.getWeight() - 1.0) > 1e-10)
      nets << " net" << edgeId;
    else
      nets << " net" << edgeId;
    nets << endl;

    unsigned nodeOffset = 0;

    itHGFNodeLocal v;
#ifdef SIGNAL_DIRECTIONS
    for (v = (*e)->srcsBegin(); v != (*e)->srcsEnd(); v++, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " O ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;
      //	    if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;
      else
        nets << endl;
    }
    for (v = (*e)->srcSnksBegin(); v != (*e)->srcSnksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " B ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;
      //          if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;
      else
        nets << endl;
    }
    for (v = (*e)->snksBegin(); v != (*e)->snksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " I ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;

      //	    if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;

      else
        nets << endl;
    }
#else
    for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " B ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;
      //          if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;
      else
        nets << endl;
    }
#endif
  }

  for (unsigned i = 0; i < congestedNodes.size(); ++i) {
    unsigned nodeId = congestedNodes[i];
    HGFNode& node = getNodeByIdx(nodeId);
    unsigned dummyNetDegree = 2;
    for (unsigned k = 0; k < numExtraNodes; ++k) {
      nets << "NetDegree : " << dummyNetDegree << " ";
      if (!getNodeNameByIndex(node.getIndex()).empty())
        nets << "dummy_net_" << k << "_" << getNodeNameByIndex(node.getIndex())
             << "  ";
      else
        nets << "dummy_net_a" << k << "_" << nodeId - getNumTerminals() << " ";
      nets << endl;

      if (!getNodeNameByIndex(node.getIndex()).empty())
        nets << setw(10) << getNodeNameByIndex(node.getIndex()) << "  ";
      else
        nets << setw(10) << " a" << nodeId - getNumTerminals() << " ";
      nets << "B " << endl;

      if (!getNodeNameByIndex(node.getIndex()).empty())
        nets << setw(10) << "dummy_" << k << "_"
             << getNodeNameByIndex(node.getIndex()) << "  ";
      else
        nets << setw(10) << "dummy_" << k << "_" << nodeId - getNumTerminals()
             << " ";
      nets << "B " << endl;
    }
  }
  nets.close();

  // write the wts file

  wts << "UCLA wts 1.0" << endl;
  wts << TimeStamp() << User() << Platform() << endl;

  unsigned n;
  if (getNumWeights() >= 1) {
    for (n = 0; n != getNumNodes(); n++) {
      HGFNode& node = *_nodes[n];
      wts << setw(10) << getNodeNameByIndex(node.getIndex()) << "  ";

      for (unsigned i = 0; i != getNumWeights(); ++i)
        wts << " " << setw(6) << getWeight(node.getIndex(), i);

      wts << endl;
    }
  }

  for (n = 0; n != getNumEdges(); ++n) {
    const HGFEdge& net = getEdgeByIdx(n);
    if (fabs(net.getWeight() - 1.0) < 1e-10) continue;
    wts << " ";
    std::string netName = getNetNameByIndex(n).c_str();
    if (!netName.empty())
      wts << setw(10) << netName;
    else
      wts << "net" << n;
    wts << " " << setw(10) << net.getWeight() << endl;
  }

  wts.close();

  aux << "HGraphWPins : ";

  const char* tmp = NULL;
  tmp = strrchr(baseFileName, pathDelim);
  if (tmp++)
    aux << tmp << ".nets " << tmp << ".nodes " << tmp << ".wts" << endl;
  else
    aux << baseFileName << ".nets " << baseFileName << ".nodes " << baseFileName
        << ".wts" << endl;
  aux.close();

}

void HGraphWDimensions::saveAsNodesNetsFloorplan(
    const char* baseFileName) const {
  const string netsFileName = fileNameWithSuffix(baseFileName, ".nets");
  const string blocksFileName = fileNameWithSuffix(baseFileName, ".blocks");

  ofstream nets(netsFileName.c_str());

  ofstream nodes(blocksFileName.c_str());

  // write the nodes file
  nodes << "UCLA blocks 1.0" << endl;
  nodes << "#" << TimeStamp() << User() << Platform() << endl;

  nodes << "NumSoftRectangularBlocks :     " << setw(10)
        << getNumNodes() - getNumTerminals() << endl;
  nodes << "NumHardRectilinearBlocks :     " << setw(10) << 0 << endl;
  nodes << "NumTerminals : " << setw(10) << getNumTerminals() << endl;

  itHGFNodeGlobal nIt;
  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    Symmetry sym = _nodeSymmetries[nodeId];

    if (!getNodeNameByIndex(node.getIndex()).empty()) {
      nodes << setw(10) << getNodeNameByIndex(node.getIndex()) << " ";
    } else {
      if (isTerminal(nodeId))
        nodes << " p" << nodeId + 1 << " ";
      else
        nodes << " a" << nodeId - getNumTerminals() << " ";
    }

    if (isTerminal(nodeId))
      nodes << "terminal ";
    else {
      float area = getNodeWidth(nodeId) * getNodeHeight(nodeId);
      float ar = getNodeWidth(nodeId) / getNodeHeight(nodeId);
      nodes << "softrectangular ";
      nodes << area << " ";
      nodes << ar << " " << ar;
      if (sym.rot90 || sym.y || sym.x) nodes << " : " << sym;
    }
    nodes << endl;
  }
  nodes.close();

  unsigned num1PinEdges = 0, num0PinEdges = 0;
  itHGFEdgeGlobal e;
  for (e = edgesBegin(); e != edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    if (edge.getDegree() > 1) continue;
    if (edge.getDegree() == 1)
      num1PinEdges++;
    else
      num0PinEdges++;
  }

  // write the nets file
  nets << "UCLA nets 1.0" << endl;
  nets << "#" << TimeStamp() << User() << Platform() << endl;

  nets << "NumNets : " << getNumEdges() - num1PinEdges << endl;
  nets << "NumPins : " << getNumPins() - num1PinEdges << endl;

  for (e = edgesBegin(); e != edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    if (edge.getDegree() < 2) continue;
    unsigned edgeId = edge.getIndex();
    // const char*     netName=getNetNameByIndex(edgeId);

    nets << "NetDegree : " << (*e)->getDegree() << " ";
    nets << endl;

    unsigned nodeOffset = 0;

    itHGFNodeLocal v;
#ifdef SIGNAL_DIRECTIONS
    for (v = (*e)->srcsBegin(); v != (*e)->srcsEnd(); v++, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " O ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;
      dx = (dx / _nodeWidths[idx]) * 100;
      dy = (dy / _nodeHeights[idx]) * 100;

      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : %" << dx << " %" << dy << endl;
      else
        nets << endl;
    }
    for (v = (*e)->srcSnksBegin(); v != (*e)->srcSnksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " B ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;
      dx = (dx / _nodeWidths[idx]) * 100;
      dy = (dy / _nodeHeights[idx]) * 100;

      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : %" << dx << " %" << dy << endl;
      else
        nets << endl;
    }
    for (v = (*e)->snksBegin(); v != (*e)->snksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " I ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;
      dx = (dx / _nodeWidths[idx]) * 100;
      dy = (dy / _nodeHeights[idx]) * 100;

      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : %" << dx << " %" << dy << endl;
      else
        nets << endl;
    }
#else
    for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " B ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;
      dx = (dx / _nodeWidths[idx]) * 100;
      dy = (dy / _nodeHeights[idx]) * 100;

      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : %" << dx << " %" << dy << endl;
      else
        nets << endl;
    }
#endif
  }
  nets.close();
}

void HGraphWDimensions::saveMacrosAsNodesNetsFloorplan(
    const char* baseFileName) const {
  const string netsFileName = fileNameWithSuffix(baseFileName, ".nets");
  const string blocksFileName = fileNameWithSuffix(baseFileName, ".blocks");

  ofstream nets(netsFileName.c_str());

  ofstream nodes(blocksFileName.c_str());

  unsigned numNodes = 0;
  itHGFNodeGlobal nIt;
  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    if (!isTerminal(nodeId)) {
      if (_isMacro[nodeId]) ++numNodes;
    }
  }

  // write the nodes file
  nodes << "UCLA blocks 1.0" << endl;
  nodes << "#" << TimeStamp() << User() << Platform() << endl;

  nodes << "NumSoftRectangularBlocks :     " << setw(10) << numNodes << endl;
  nodes << "NumHardRectilinearBlocks :     " << setw(10) << 0 << endl;
  nodes << "NumTerminals : " << setw(10) << getNumTerminals() << endl;

  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    Symmetry sym = _nodeSymmetries[nodeId];

    if (isTerminal(nodeId)) {
      if (!getNodeNameByIndex(node.getIndex()).empty())
        nodes << setw(10) << getNodeNameByIndex(node.getIndex()) << " ";
      else
        nodes << " p" << nodeId + 1 << " ";
      nodes << "terminal ";
      nodes << endl;
    } else {
      if (_isMacro[nodeId]) {
        if (!getNodeNameByIndex(node.getIndex()).empty())
          nodes << setw(10) << getNodeNameByIndex(node.getIndex()) << " ";
        else
          nodes << " a" << nodeId - getNumTerminals() << " ";

        nodes << "softrectangular ";
        nodes << getNodeWidth(nodeId) * getNodeHeight(nodeId) << " ";
        nodes << "1.0 1.0";
        if (sym.rot90 || sym.y || sym.x) nodes << " : " << sym;
        nodes << endl;
      }
    }
  }
  nodes.close();

  // write the nets file
  nets << "UCLA nets 1.0" << endl;
  nets << "#" << TimeStamp() << User() << Platform() << endl;

  nets << "NumNets : " << 0 << endl;
  nets << "NumPins : " << 0 << endl;
  nets.close();
}

void HGraphWDimensions::saveAsNodesNetsWtsWDummy(const char* baseFileName,
                                                 float siteWidth,
                                                 float siteHeight,
                                                 float areaDummy,
                                                 bool fixMacros) const {
  const string netsFileName = fileNameWithSuffix(baseFileName, ".nets");
  const string nodesFileName = fileNameWithSuffix(baseFileName, ".nodes");
  const string wtsFileName = fileNameWithSuffix(baseFileName, ".wts");
  const string auxFileName = fileNameWithSuffix(baseFileName, ".aux");

  ofstream nets(netsFileName.c_str());

  ofstream nodes(nodesFileName.c_str());

  ofstream wts(wtsFileName.c_str());

  ofstream aux(auxFileName.c_str());

  itHGFNodeGlobal nIt;

  // get number of macros
  unsigned numMacros = 0;
  if (fixMacros) {
    for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
      HGFNode& node = **nIt;
      unsigned nodeId = node.getIndex();
      if (!isTerminal(nodeId)) {
        if (_isMacro[nodeId]) ++numMacros;
      }
    }
  }
  // write the nodes file
  float avgWidth = 0;
  float numCore = 0;

  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    if (isTerminal(nodeId) || _isMacro[nodeId]) continue;
    numCore++;
    avgWidth += getNodeWidth(nodeId);
  }
  avgWidth /= numCore;
  float dummyWidth = 4 * floor(avgWidth / siteWidth) * siteWidth;
  float dummyHeight = siteHeight;
  float dummyArea = dummyWidth * dummyHeight;
  unsigned numDummy = unsigned(floor(areaDummy / dummyArea));

  nodes << "UCLA nodes 1.0" << endl;
  nodes << TimeStamp() << User() << Platform() << endl;

  nodes << "NumNodes :     " << setw(10) << getNumNodes() + numDummy << endl;
  nodes << "NumTerminals : " << setw(10) << getNumTerminals() + numMacros
        << endl;

  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    Symmetry sym = _nodeSymmetries[nodeId];

    if (!getNodeNameByIndex(node.getIndex()).empty())
      nodes << setw(10) << getNodeNameByIndex(node.getIndex()) << "  ";
    else {
      if (isTerminal(nodeId))
        nodes << " p" << nodeId + 1 << " ";
      else
        nodes << " a" << nodeId - getNumTerminals() << " ";
    }

    if (getNodeWidth(nodeId) != 0 || getNodeHeight(nodeId) != 0)
      nodes << setw(10) << getNodeWidth(nodeId) << " " << setw(10)
            << getNodeHeight(nodeId) << " ";
    if (isTerminal(nodeId))
      nodes << "terminal ";
    else if (fixMacros && _isMacro[nodeId])
      nodes << "terminal ";

    if (sym.rot90 || sym.y || sym.x) nodes << " : " << sym;
    nodes << endl;
  }
  for (unsigned i = 0; i < numDummy; ++i)
    nodes << setw(10) << "dummy" << i << " " << dummyWidth << "  "
          << dummyHeight << endl;

  nodes.close();

  unsigned num1PinEdges = 0, num0PinEdges = 0;
  itHGFEdgeGlobal e;
  for (e = edgesBegin(); e != edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    if (edge.getDegree() > 1) continue;
    if (edge.getDegree() == 1)
      num1PinEdges++;
    else
      num0PinEdges++;
  }

  // write the nets file
  nets << "UCLA nets 1.0" << endl;
  nets << TimeStamp() << User() << Platform() << endl;

  nets << "NumNets : " << getNumEdges() - num1PinEdges << endl;
  nets << "NumPins : " << getNumPins() - num1PinEdges << endl;

  for (e = edgesBegin(); e != edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    if (edge.getDegree() < 2) continue;
    unsigned edgeId = edge.getIndex();
    const char* netName = getNetNameByIndex(edgeId).c_str();

    nets << "NetDegree : " << (*e)->getDegree() << " ";
    if (netName)
      nets << setw(10) << netName;
    else if (fabs(edge.getWeight() - 1.0) > 1e-10)
      nets << " net" << edgeId;
    else
      nets << " net" << edgeId;
    nets << endl;

    unsigned nodeOffset = 0;

    itHGFNodeLocal v;
#ifdef SIGNAL_DIRECTIONS
    for (v = (*e)->srcsBegin(); v != (*e)->srcsEnd(); v++, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " O ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;
      //	    if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;
      else
        nets << endl;
    }
    for (v = (*e)->srcSnksBegin(); v != (*e)->srcSnksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " B ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;
      //          if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;
      else
        nets << endl;
    }
    for (v = (*e)->snksBegin(); v != (*e)->snksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " I ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;

      //	    if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;

      else
        nets << endl;
    }
#else
    for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " B ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;
      //          if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;
      else
        nets << endl;
    }
#endif
  }
  nets.close();

  // write the wts file

  wts << "UCLA wts 1.0" << endl;
  wts << TimeStamp() << User() << Platform() << endl;

  unsigned n;
  if (getNumWeights() >= 1) {
    for (n = 0; n != getNumNodes(); n++) {
      HGFNode& node = *_nodes[n];
      wts << setw(10) << getNodeNameByIndex(node.getIndex()) << "  ";

      for (unsigned i = 0; i != getNumWeights(); ++i)
        wts << " " << setw(6) << getWeight(node.getIndex(), i);

      wts << endl;
    }
  }

  for (n = 0; n != getNumEdges(); ++n) {
    const HGFEdge& net = getEdgeByIdx(n);
    if (fabs(net.getWeight() - 1.0) < 1e-10) continue;
    wts << " ";
    std::string netName = getNetNameByIndex(n).c_str();
    if (!netName.empty())
      wts << setw(10) << netName;
    else
      wts << "net" << n;
    wts << " " << setw(10) << net.getWeight() << endl;
  }

  wts.close();

  aux << "HGraphWPins : ";

  const char* tmp = NULL;
  tmp = strrchr(baseFileName, pathDelim);
  if (tmp++)
    aux << tmp << ".nets " << tmp << ".nodes " << tmp << ".wts" << endl;
  else
    aux << baseFileName << ".nets " << baseFileName << ".nodes " << baseFileName
        << ".wts" << endl;
  aux.close();

}

void HGraphWDimensions::getCellsToTetherBFS(
    std::vector<unsigned>& tetheredCells, unsigned numTetherCells,
    bool takeNodeConstrFrmFile, bool takeNetConstrFrmFile) const {
  unsigned numActTerms = getNumTerminals();
  unsigned numNodes = getNumNodes();
  vector<bool> iamTethered(numNodes, false);

  if (0) {
    for (unsigned i = 0; i < numTetherCells;) {
      unsigned randIdx = RandomUnsigned(numActTerms, numNodes);
      if (isTerminal(randIdx) || iamTethered[randIdx]) continue;
      tetheredCells.push_back(randIdx);
      iamTethered[randIdx] = true;
      ++i;
    }
  } else {
    vector<bool> seenCells(numNodes, false);
    while (1) {
      deque<unsigned> cells;
      for (unsigned i = numActTerms; i < numNodes; ++i) {
        if (!seenCells[i]) {
          cells.push_back(i);
          break;
        }
      }
      if (cells.size() == 0) break;
      // the BFS loop
      seenCells[cells[0]] = true;
      vector<unsigned> groupCells;
      groupCells.push_back(cells[0]);
      while (1) {
        if (cells.size() == 0) break;
        unsigned rootIdx = cells[0];
        const HGFNode& origNode = getNodeByIdx(rootIdx);
        itHGFEdgeLocal e;
        itHGFNodeLocal v;
        for (e = origNode.edgesBegin(); e != origNode.edgesEnd(); ++e) {
          for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); v++) {
            HGFNode& node = (**v);
            unsigned nodeIndex = node.getIndex();
            if (!seenCells[nodeIndex]) {
              seenCells[nodeIndex] = true;
              cells.push_back(nodeIndex);
              groupCells.push_back(nodeIndex);
            }
          }
        }
        // std::vector<unsigned>::iterator cellsBegin = cells.begin();
        // cout<<"cellsBegin "<<*cellsBegin<<endl;
        // cells.erase(cellsBegin);
        cells.pop_front();
      }
      for (unsigned i = 0; i < groupCells.size(); ++i) {
        unsigned randIdx = RandomUnsigned(0, groupCells.size());
        if (!iamTethered[groupCells[randIdx]] &&
            !isTerminal(groupCells[randIdx])) {
          tetheredCells.push_back(groupCells[randIdx]);
          iamTethered[groupCells[randIdx]] = true;
          break;
        }
      }
      cout << "Group Size " << groupCells.size() << endl;
      groupCells.clear();
      cells.clear();
    }
    unsigned numBFSTethered = tetheredCells.size();

    // now check if a node constraint file specified.
    if (takeNodeConstrFrmFile) {
      ifstream constrFile("constraints.nodes");
      abkfatal(constrFile, " Could not open constraints.nodes");
      cout << " Reading constraints.nodes ... " << endl;
      string nodeName;
      int lineno = 0;

      constrFile >> eathash(lineno);
      while (!constrFile.eof()) {
        constrFile >> eathash(lineno) >> nodeName >> eatblank;
        unsigned nodeId = getNodeByName(nodeName.c_str()).getIndex();
        if (!isTerminal(nodeId) && !iamTethered[nodeId]) {
          tetheredCells.push_back(nodeId);
          iamTethered[nodeId] = true;
        }
        constrFile >> skiptoeol >> eathash(lineno);
      }
      constrFile.close();
    }
    unsigned numFileTethered = tetheredCells.size() - numBFSTethered;

    // read the net constraint file. Find nodes on nets. do not tether these
    // nodes during random tethering
    vector<bool> nodesNotToTether(numNodes, false);
    unsigned numNotTethered = 0;
    if (takeNetConstrFrmFile) {
      ifstream constrFile("constraints.nets");
      abkfatal(constrFile, " Could not open constraints.nets");
      cout << " Reading constraints.nets ... " << endl;
      string netName;
      int lineno = 0;

      constrFile >> eathash(lineno);
      while (!constrFile.eof()) {
        constrFile >> eathash(lineno) >> netName >> eatblank;
        const HGFEdge& edge = getNetByName(netName.c_str());
        for (itHGFNodeLocal nodeOnEdge = edge.nodesBegin();
             nodeOnEdge != edge.nodesEnd(); ++nodeOnEdge) {
          nodesNotToTether[(*nodeOnEdge)->getIndex()] = true;
          numNotTethered++;
        }
        constrFile >> skiptoeol >> eathash(lineno);
      }
      constrFile.close();
    }

    // the remaining cells are random
    for (unsigned i = tetheredCells.size(); i < numTetherCells;) {
      unsigned randIdx = RandomUnsigned(numActTerms, numNodes);
      if (isTerminal(randIdx) || iamTethered[randIdx] ||
          nodesNotToTether[randIdx])
        continue;
      tetheredCells.push_back(randIdx);
      iamTethered[randIdx] = true;
      ++i;
    }
    unsigned numRandomTethered =
        tetheredCells.size() - numBFSTethered - numFileTethered;
    cout << "numBFSTethered " << numBFSTethered << " numFileTethered "
         << numFileTethered << " numRandomTethered " << numRandomTethered
         << " numNotTethered " << numNotTethered << endl;
  }
}

void HGraphWDimensions::saveAsNodesNetsWtsWTether(
    const char* baseFileName, double fract, std::vector<unsigned>& tetheredCells,
    bool takeNodeConstrFrmFile, bool takeNetConstrFrmFile) const {
  const string netsFileName = fileNameWithSuffix(baseFileName, ".nets");
  const string nodesFileName = fileNameWithSuffix(baseFileName, ".nodes");
  const string wtsFileName = fileNameWithSuffix(baseFileName, ".wts");
  const string auxFileName = fileNameWithSuffix(baseFileName, ".aux");

  ofstream nets(netsFileName.c_str());

  ofstream nodes(nodesFileName.c_str());

  ofstream wts(wtsFileName.c_str());

  ofstream aux(auxFileName.c_str());

  itHGFNodeGlobal nIt;

  // get number of macros
  unsigned numMacros = 0;
  bool fixMacros = 0;
  if (fixMacros) {
    for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
      HGFNode& node = **nIt;
      unsigned nodeId = node.getIndex();
      if (!isTerminal(nodeId)) {
        if (_isMacro[nodeId]) ++numMacros;
      }
    }
  }
  // write the nodes file
  //  float avgWidth=0;
  float numCore = getNumNodes() - getNumTerminals();

  unsigned numTetherCells = unsigned(floor(fract * numCore / 100));
  unsigned numTetherTerms = 4 * numTetherCells;
  //  unsigned numActTerms = getNumTerminals();
  //  unsigned numNodes = getNumNodes();

  getCellsToTetherBFS(tetheredCells, numTetherCells, takeNodeConstrFrmFile,
                      takeNetConstrFrmFile);
  if (numTetherCells != tetheredCells.size()) {
    cout << "numTetherCells " << numTetherCells << " vs. tetheredCells.size() "
         << tetheredCells.size() << endl;
    abkwarn(0, "numTetherCells != tetheredCells.size() \n");
  }
  numTetherCells = tetheredCells.size();
  numTetherTerms = 4 * numTetherCells;

  nodes << "UCLA nodes 1.0" << endl;
  nodes << TimeStamp() << User() << Platform() << endl;

  nodes << "NumNodes :     " << setw(10) << getNumNodes() + numTetherTerms
        << endl;
  nodes << "NumTerminals : " << setw(10)
        << getNumTerminals() + numMacros + numTetherTerms << endl;

  unsigned numTermsAdded = 0;
  for (unsigned i = 0; i < numTetherCells; ++i) {
    for (unsigned j = 0; j < 4; ++j) {
      nodes << setw(10) << "dummyTerm" << numTermsAdded << " 1 1 terminal "
            << endl;
      ++numTermsAdded;
    }
  }

  if (numTermsAdded != numTetherTerms) {
    cout << "numTermsAdded " << numTermsAdded << " numTermsExpected "
         << numTetherTerms << endl;
    abkfatal(0, "numTermsAdded != numTetherTerms\n");
  }

  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    Symmetry sym = _nodeSymmetries[nodeId];

    if (!getNodeNameByIndex(node.getIndex()).empty())
      nodes << setw(10) << getNodeNameByIndex(node.getIndex()) << "  ";
    else {
      if (isTerminal(nodeId))
        nodes << " p" << nodeId + 1 << " ";
      else
        nodes << " a" << nodeId - getNumTerminals() << " ";
    }

    if (getNodeWidth(nodeId) != 0 || getNodeHeight(nodeId) != 0)
      nodes << setw(10) << getNodeWidth(nodeId) << " " << setw(10)
            << getNodeHeight(nodeId) << " ";
    if (isTerminal(nodeId))
      nodes << "terminal ";
    else if (fixMacros && _isMacro[nodeId])
      nodes << "terminal ";

    if (sym.rot90 || sym.y || sym.x) nodes << " : " << sym;
    nodes << endl;
  }

  nodes.close();

  unsigned num1PinEdges = 0, num0PinEdges = 0;
  itHGFEdgeGlobal e;
  for (e = edgesBegin(); e != edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    if (edge.getDegree() > 1) continue;
    if (edge.getDegree() == 1)
      num1PinEdges++;
    else
      num0PinEdges++;
  }

  unsigned numTetherNets = numTetherTerms;

  // write the nets file
  nets << "UCLA nets 1.0" << endl;
  nets << TimeStamp() << User() << Platform() << endl;

  nets << "NumNets : " << getNumEdges() - num1PinEdges + numTetherNets << endl;
  nets << "NumPins : " << getNumPins() - num1PinEdges + 2 * numTetherNets
       << endl;

  unsigned tetherTermIdx = 0;
  for (unsigned i = 0; i < numTetherCells; ++i) {
    const HGFNode& node = getNodeByIdx(tetheredCells[i]);
    for (unsigned j = 0; j < 4; ++j) {
      nets << "NetDegree : 2 " << endl;
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " B " << endl;
      nets << setw(10) << "dummyTerm" << tetherTermIdx << " B " << endl;
      ++tetherTermIdx;
    }
  }

  for (e = edgesBegin(); e != edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    if (edge.getDegree() < 2) continue;
    unsigned edgeId = edge.getIndex();
    const char* netName = getNetNameByIndex(edgeId).c_str();

    nets << "NetDegree : " << (*e)->getDegree() << " ";
    if (netName)
      nets << setw(10) << netName;
    else if (fabs(edge.getWeight() - 1.0) > 1e-10)
      nets << " net" << edgeId;
    else
      nets << " net" << edgeId;
    nets << endl;

    unsigned nodeOffset = 0;

    itHGFNodeLocal v;
#ifdef SIGNAL_DIRECTION
    for (v = (*e)->srcsBegin(); v != (*e)->srcsEnd(); v++, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " O ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;
      //	    if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;
      else
        nets << endl;
    }
    for (v = (*e)->srcSnksBegin(); v != (*e)->srcSnksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " B ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;
      //          if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;
      else
        nets << endl;
    }
    for (v = (*e)->snksBegin(); v != (*e)->snksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " I ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;

      //	    if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;

      else
        nets << endl;
    }
#else
    for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << getNodeNameByIndex(node.getIndex()) << " B ";

      Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
      float dx = pOffset.x - _nodeWidths[idx] / 2,
            dy = pOffset.y - _nodeHeights[idx] / 2;
      //          if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;
      else
        nets << endl;
    }
#endif
  }
  nets.close();

  // write the wts file

  wts << "UCLA wts 1.0" << endl;
  wts << TimeStamp() << User() << Platform() << endl;

  unsigned n;
  if (getNumWeights() >= 1) {
    for (n = 0; n != getNumNodes(); n++) {
      HGFNode& node = *_nodes[n];
      wts << setw(10) << getNodeNameByIndex(node.getIndex()) << "  ";

      for (unsigned i = 0; i != getNumWeights(); ++i)
        wts << " " << setw(6) << getWeight(node.getIndex(), i);

      wts << endl;
    }
  }

  // increase the weight of these nets by certain factor
  vector<bool> netsToIncreaseWt(getNumEdges(), false);
  if (takeNetConstrFrmFile) {
    ifstream constrFile("constraints.nets");
    abkfatal(constrFile, " Could not open constraints.nets");
    string netName;
    int lineno = 0;
    constrFile >> eathash(lineno);
    while (!constrFile.eof()) {
      constrFile >> eathash(lineno) >> netName >> eatblank;
      netsToIncreaseWt[getNetByName(netName.c_str()).getIndex()] = true;
      constrFile >> skiptoeol >> eathash(lineno);
    }
    constrFile.close();
  }

  int boostFactor = 2;
  for (n = 0; n != getNumEdges(); ++n) {
    const HGFEdge& net = getEdgeByIdx(n);
    if (fabs(net.getWeight() - 1.0) < 1e-10 && !netsToIncreaseWt[n]) continue;
    wts << " ";
    const char* netName = getNetNameByIndex(n).c_str();
    if (netName)
      wts << setw(10) << netName;
    else
      wts << "net" << n;
    if (netsToIncreaseWt[n])
      wts << " " << setw(10) << boostFactor * net.getWeight() << endl;
    else
      wts << " " << setw(10) << net.getWeight() << endl;
  }

  wts.close();

  aux << "HGraphWPins : ";

  const char* tmp = NULL;
  tmp = strrchr(baseFileName, pathDelim);
  if (tmp++)
    aux << tmp << ".nets " << tmp << ".nodes " << tmp << ".wts" << endl;
  else
    aux << baseFileName << ".nets " << baseFileName << ".nodes " << baseFileName
        << ".wts" << endl;
  aux.close();

}
