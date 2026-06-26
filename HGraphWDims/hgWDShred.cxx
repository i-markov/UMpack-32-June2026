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

#include <cstdlib>
#include <cerrno>
#include <climits>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/infolines.h"
#include "ABKCommon/pathDelims.h"
#include "hgWDims.h"

using std::endl;
using std::ofstream;
using std::ostringstream;
using std::setw;
using std::string;

#define NUM_NEW_VERT_NETS 1   // vertical nets added per node during shredding
#define NUM_NEW_HORIZ_NETS 1  // horizontal nets added per node during shredding

static string shredFileNameWithSuffix(const char* baseFileName,
                                      const char* suffix) {
  return string(baseFileName) + suffix;
}

struct ShreddedNodeName {
  string baseName;
  unsigned subIndex;
  bool hasSubIndex;
};

static unsigned parseUnsignedSuffix(const char* token, const char* context) {
  abkfatal(token != NULL && *token != '\0', context);
  errno = 0;
  char* end = NULL;
  const unsigned long parsed = strtoul(token, &end, 10);
  abkfatal(end != token && *end == '\0' && errno != ERANGE &&
               parsed <= UINT_MAX,
           context);
  return static_cast<unsigned>(parsed);
}

static ShreddedNodeName parseShreddedNodeName(const string& nodeName) {
  ShreddedNodeName parsed;
  parsed.baseName = nodeName;
  parsed.subIndex = 0;
  parsed.hasSubIndex = false;

  const string::size_type pos = nodeName.find('@');
  if (pos == string::npos) return parsed;

  parsed.baseName = nodeName.substr(0, pos);
  parsed.subIndex = parseUnsignedSuffix(nodeName.c_str() + pos + 1,
                                        "Malformed shredded node suffix");
  parsed.hasSubIndex = true;
  return parsed;
}

void HGraphWDimensions::saveAsNodesNetsWtsUnShred(
    const char* baseFileName, float maxCoreRowHeight) const {
  const string netsFileName = shredFileNameWithSuffix(baseFileName, ".nets");
  const string nodesFileName = shredFileNameWithSuffix(baseFileName, ".nodes");
  const string wtsFileName = shredFileNameWithSuffix(baseFileName, ".wts");
  const string auxFileName = shredFileNameWithSuffix(baseFileName, ".aux");
  std::map<string, float> widths;
  std::map<string, float> heights;

  ofstream nets(netsFileName.c_str());

  ofstream nodes(nodesFileName.c_str());

  ofstream wts(wtsFileName.c_str());

  ofstream aux(auxFileName.c_str());

  unsigned numExtraNodes = 0;
  unsigned numExtraNets = 0;
  itHGFNodeGlobal nIt;

  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    string nodeName = getNodeNameByIndex(node.getIndex()).c_str();
    string actualName = nodeName;
    unsigned long pos = actualName.find('@');

    if (pos != string::npos)  // found
    {
      string rest(actualName, pos);
      actualName.resize(pos);
      string subNodeName;
      while (nIt != nodesEnd()) {
        HGFNode& subNode = **nIt;
        nodeName = getNodeNameByIndex(subNode.getIndex()).c_str();
        subNodeName = nodeName;
        pos = subNodeName.find('@');
        if (pos == string::npos) break;
        string rest(subNodeName, pos);
        subNodeName.resize(pos);
        if (actualName != subNodeName) break;

        ++numExtraNodes;
        numExtraNets += NUM_NEW_VERT_NETS;
        ++nIt;
      }
      --nIt;
      --numExtraNodes;
      numExtraNets -= NUM_NEW_VERT_NETS;
    }
  }

  // write the nodes file
  nodes << "UCLA nodes 1.0" << endl;
  nodes << TimeStamp() << User() << Platform() << endl;

  nodes << "NumNodes :     " << setw(10) << getNumNodes() - numExtraNodes
        << endl;
  nodes << "NumTerminals : " << setw(10) << getNumTerminals() << endl;

  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    Symmetry sym = _nodeSymmetries[nodeId];

    string nodeName = getNodeNameByIndex(node.getIndex()).c_str();
    string actualName = nodeName;
    unsigned long pos = actualName.find('@');

    string tempName;

    if (pos != string::npos)  // found
    {
      string rest(actualName, pos);
      string subNodeName;
      float nodeHeight = 0;
      float nodeWidth = getNodeWidth(nodeId);

      while (nIt != nodesEnd()) {
        HGFNode& subNode = **nIt;
        unsigned subNodeId = subNode.getIndex();
        nodeName = getNodeNameByIndex(subNode.getIndex()).c_str();
        subNodeName = nodeName;
        pos = subNodeName.find('@');
        if (pos == string::npos) break;
        string rest2(subNodeName, pos);
        subNodeName.resize(pos);
        if (actualName != subNodeName) break;

        nodeHeight += getNodeHeight(subNodeId);
        nIt++;
      }
      --nIt;
      nodes << setw(10) << actualName << "  ";
      nodes << setw(10) << nodeWidth << " ";
      nodes << setw(10) << nodeHeight << " ";
      if (sym.rot90 || sym.y || sym.x) nodes << " : " << sym;
      nodes << endl;

      tempName = actualName;
      widths[tempName] = nodeWidth;
      heights[tempName] = nodeHeight;
    } else {
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

  nets << "NumNets : " << getNumEdges() - num1PinEdges - numExtraNets << endl;
  nets << "NumPins : " << getNumPins() - num1PinEdges - numExtraNets * 2
       << endl;

  for (e = edgesBegin(); e != edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    if (edge.getDegree() < 2) continue;
    unsigned edgeId = edge.getIndex();
    const char* netName = getNetNameByIndex(edgeId).c_str();
    const char* pos = strchr(netName, '@');
    if (pos != NULL) continue;

    nets << "NetDegree : " << (*e)->getDegree() << " ";
    if (netName)
      nets << setw(10) << netName;
    else if (fabs(edge.getWeight() - 1.0) > 1e-10)
      nets << " net" << edgeId;
    nets << endl;

    unsigned nodeOffset = 0;

    itHGFNodeLocal v;
    float nodeHeight;
    float nodeWidth;
    unsigned subIdx;
#ifdef SIGNAL_DIRECTIONS
    for (v = (*e)->srcsBegin(); v != (*e)->srcsEnd(); v++, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      const ShreddedNodeName nodeName =
          parseShreddedNodeName(getNodeNameByIndex(node.getIndex()));
      if (nodeName.hasSubIndex) {
        subIdx = nodeName.subIndex;
        nodeWidth = widths[nodeName.baseName];
        nodeHeight = heights[nodeName.baseName];

        Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
        float dx = pOffset.x - nodeWidth / 2, dy = pOffset.y;
        dy = dy + subIdx * maxCoreRowHeight - nodeHeight / 2;

        nets << setw(10) << nodeName.baseName << " O ";
        if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
          nets << " : " << dx << " " << dy << endl;
        else
          nets << endl;
      } else {
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
    }
    for (v = (*e)->srcSnksBegin(); v != (*e)->srcSnksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      const ShreddedNodeName nodeName =
          parseShreddedNodeName(getNodeNameByIndex(node.getIndex()));

      if (nodeName.hasSubIndex) {
        subIdx = nodeName.subIndex;
        nodeWidth = widths[nodeName.baseName];
        nodeHeight = heights[nodeName.baseName];

        //		nodeWidth = widths[actualName];
        // nodeHeight = heights[actualName];
        Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
        float dx = pOffset.x - nodeWidth / 2, dy = pOffset.y;
        dy = dy + subIdx * maxCoreRowHeight - nodeHeight / 2;

        // cout<<actualName<<nodeWidth<<"\t"<<nodeHeight<<endl;

        nets << setw(10) << nodeName.baseName << " B ";
        if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
          nets << " : " << dx << " " << dy << endl;
        else
          nets << endl;
      } else {
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
    }
    for (v = (*e)->snksBegin(); v != (*e)->snksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      const ShreddedNodeName nodeName =
          parseShreddedNodeName(getNodeNameByIndex(node.getIndex()));

      if (nodeName.hasSubIndex) {
        subIdx = nodeName.subIndex;
        nodeWidth = widths[nodeName.baseName];
        nodeHeight = heights[nodeName.baseName];
        Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
        float dx = pOffset.x - nodeWidth / 2, dy = pOffset.y;
        dy = dy + subIdx * maxCoreRowHeight - nodeHeight / 2;

        nets << setw(10) << nodeName.baseName << " I ";
        if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
          nets << " : " << dx << " " << dy << endl;
        else
          nets << endl;
      } else {
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
    }
#else
    for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      const ShreddedNodeName nodeName =
          parseShreddedNodeName(getNodeNameByIndex(node.getIndex()));

      if (nodeName.hasSubIndex) {
        subIdx = nodeName.subIndex;
        nodeWidth = widths[nodeName.baseName];
        nodeHeight = heights[nodeName.baseName];

        //		nodeWidth = widths[actualName];
        // nodeHeight = heights[actualName];
        Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
        float dx = pOffset.x - nodeWidth / 2, dy = pOffset.y;
        dy = dy + subIdx * maxCoreRowHeight - nodeHeight / 2;

        // cout<<actualName<<nodeWidth<<"\t"<<nodeHeight<<endl;

        nets << setw(10) << nodeName.baseName << " B ";
        if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
          nets << " : " << dx << " " << dy << endl;
        else
          nets << endl;
      } else {
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

  float weight;

  for (n = 0; n != getNumEdges(); ++n) {
    const HGFEdge& net = getEdgeByIdx(n);
    const char* netName = getNetNameByIndex(n).c_str();
    const char* pos = strchr(netName, '@');
    if (pos != NULL) continue;

    weight = net.getWeight();

    if (fabs(weight - 1.0) < 1e-10) continue;
    wts << " ";
    if (netName)
      wts << setw(10) << netName;
    else
      wts << "net" << n;

    wts << " " << setw(10) << weight << endl;
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
void HGraphWDimensions::saveAsNodesNetsWtsShred(const char* baseFileName,
                                                float maxCoreRowHeight) const {
  const string netsFileName = shredFileNameWithSuffix(baseFileName, ".nets");
  const string nodesFileName = shredFileNameWithSuffix(baseFileName, ".nodes");
  const string wtsFileName = shredFileNameWithSuffix(baseFileName, ".wts");
  const string auxFileName = shredFileNameWithSuffix(baseFileName, ".aux");

  ofstream nets(netsFileName.c_str());

  ofstream nodes(nodesFileName.c_str());

  ofstream wts(wtsFileName.c_str());

  ofstream aux(auxFileName.c_str());

  unsigned numExtraNodes = 0;
  unsigned numExtraNets = 0;
  itHGFNodeGlobal nIt;

  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    float nodeHeight = getNodeHeight(nodeId);
    if (nodeHeight > maxCoreRowHeight && !isTerminal(nodeId)) {
      unsigned numNewNodes = unsigned(ceil(nodeHeight / maxCoreRowHeight));
      numExtraNodes += numNewNodes - 1;
      numExtraNets += (numNewNodes - 1) * NUM_NEW_VERT_NETS;
    }
  }

  // write the nodes file
  nodes << "UCLA nodes 1.0" << endl;
  nodes << TimeStamp() << User() << Platform() << endl;

  nodes << "NumNodes :     " << setw(10) << getNumNodes() + numExtraNodes
        << endl;
  nodes << "NumTerminals : " << setw(10) << getNumTerminals() << endl;

  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    Symmetry sym = _nodeSymmetries[nodeId];

    if (getNodeHeight(nodeId) > maxCoreRowHeight && !isTerminal(nodeId))
    // shred these nodes
    {
      float nodeHeight = getNodeHeight(nodeId);
      float remainingHeight = nodeHeight;
      float thisNodeHeight = maxCoreRowHeight;
      string nodeName;
      unsigned numNewNodes = unsigned(ceil(nodeHeight / maxCoreRowHeight));
      for (unsigned i = 0; i < numNewNodes; ++i) {
        ostringstream ss;
        ss << getNodeNameByIndex(node.getIndex()) << "@" << i;
        nodeName = ss.str().c_str();
        nodes << setw(10) << nodeName << "  ";
        nodes << setw(10) << getNodeWidth(nodeId) << " ";
        if (remainingHeight < maxCoreRowHeight)
          thisNodeHeight = remainingHeight;
        else
          thisNodeHeight = maxCoreRowHeight;

        nodes << setw(10) << thisNodeHeight << " ";
        remainingHeight -= thisNodeHeight;

        if (sym.rot90 || sym.y || sym.x) nodes << " : " << sym;
        nodes << endl;
      }
    } else {
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

  nets << "NumNets : " << getNumEdges() - num1PinEdges + numExtraNets << endl;
  nets << "NumPins : " << getNumPins() - num1PinEdges + numExtraNets * 2
       << endl;

  for (e = edgesBegin(); e != edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    if (edge.getDegree() < 2) continue;
    unsigned edgeId = edge.getIndex();
    string netName = getNetNameByIndex(edgeId).c_str();

    nets << "NetDegree : " << (*e)->getDegree() << " ";
    if (!netName.empty())
      nets << setw(10) << netName;
    else if (fabs(edge.getWeight() - 1.0) > 1e-10)
      nets << " net" << edgeId;
    nets << endl;

    unsigned nodeOffset = 0;

    itHGFNodeLocal v;
    float nodeHeight;
    unsigned subIdx;
    unsigned maxSubIdx;
    float thisNodeHeight;
#ifdef SIGNAL_DIRECTIONS
    for (v = (*e)->srcsBegin(); v != (*e)->srcsEnd(); v++, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();

      string nodeName;
      nodeHeight = getNodeHeight(idx);
      if (nodeHeight > maxCoreRowHeight) {
        maxSubIdx = unsigned(ceil(nodeHeight / maxCoreRowHeight) - 1);
        Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
        float dx = pOffset.x - _nodeWidths[idx] / 2, dy = pOffset.y;
        if (dy > maxCoreRowHeight)
          subIdx = unsigned(ceil(dy / maxCoreRowHeight) - 1);
        else
          subIdx = 0;

        if (subIdx > maxSubIdx) subIdx = maxSubIdx;

        ostringstream ss;
        ss << getNodeNameByIndex(node.getIndex()) << "@" << subIdx;
        nodeName = ss.str().c_str();
        nets << setw(10) << nodeName << " O ";

        if (subIdx == maxSubIdx)
          thisNodeHeight = nodeHeight - subIdx * maxCoreRowHeight;
        else
          thisNodeHeight = maxCoreRowHeight;

        dy -= subIdx * maxCoreRowHeight;
        dy -= thisNodeHeight / 2;

        if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
          nets << " : " << dx << " " << dy << endl;
        else
          nets << endl;
      } else {
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
    }
    for (v = (*e)->srcSnksBegin(); v != (*e)->srcSnksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();

      string nodeName;
      nodeHeight = getNodeHeight(idx);
      if (nodeHeight > maxCoreRowHeight) {
        maxSubIdx = unsigned(ceil(nodeHeight / maxCoreRowHeight) - 1);
        Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
        float dx = pOffset.x - _nodeWidths[idx] / 2, dy = pOffset.y;
        if (dy > maxCoreRowHeight)
          subIdx = unsigned(ceil(dy / maxCoreRowHeight) - 1);
        else
          subIdx = 0;

        if (subIdx > maxSubIdx) subIdx = maxSubIdx;

        ostringstream ss;
        ss << getNodeNameByIndex(node.getIndex()) << "@" << subIdx;
        nodeName = ss.str().c_str();

        nets << setw(10) << nodeName << " B ";

        if (subIdx == maxSubIdx)
          thisNodeHeight = nodeHeight - subIdx * maxCoreRowHeight;
        else
          thisNodeHeight = maxCoreRowHeight;

        dy -= subIdx * maxCoreRowHeight;
        dy -= thisNodeHeight / 2;

        if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
          nets << " : " << dx << " " << dy << endl;
        else
          nets << endl;
      } else {
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
    }
    for (v = (*e)->snksBegin(); v != (*e)->snksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();

      string nodeName;
      nodeHeight = getNodeHeight(idx);
      if (nodeHeight > maxCoreRowHeight) {
        maxSubIdx = unsigned(ceil(nodeHeight / maxCoreRowHeight) - 1);
        Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
        float dx = pOffset.x - _nodeWidths[idx] / 2, dy = pOffset.y;
        if (dy > maxCoreRowHeight)
          subIdx = unsigned(ceil(dy / maxCoreRowHeight) - 1);
        else
          subIdx = 0;

        if (subIdx > maxSubIdx) subIdx = maxSubIdx;

        ostringstream ss;
        ss << getNodeNameByIndex(node.getIndex()) << "@" << subIdx;
        nodeName = ss.str().c_str();
        nets << setw(10) << nodeName << " I ";

        if (subIdx == maxSubIdx)
          thisNodeHeight = nodeHeight - subIdx * maxCoreRowHeight;
        else
          thisNodeHeight = maxCoreRowHeight;

        dy -= subIdx * maxCoreRowHeight;
        dy -= thisNodeHeight / 2;

        if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
          nets << " : " << dx << " " << dy << endl;
        else
          nets << endl;
      } else {
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
    }
#else
    for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();

      string nodeName;
      nodeHeight = getNodeHeight(idx);
      if (nodeHeight > maxCoreRowHeight) {
        maxSubIdx = unsigned(ceil(nodeHeight / maxCoreRowHeight) - 1);
        Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
        float dx = pOffset.x - _nodeWidths[idx] / 2, dy = pOffset.y;
        if (dy > maxCoreRowHeight)
          subIdx = unsigned(ceil(dy / maxCoreRowHeight) - 1);
        else
          subIdx = 0;

        if (subIdx > maxSubIdx) subIdx = maxSubIdx;

        ostringstream ss;
        ss << getNodeNameByIndex(node.getIndex()) << "@" << subIdx;
        nodeName = ss.str().c_str();

        nets << setw(10) << nodeName << " B ";

        if (subIdx == maxSubIdx)
          thisNodeHeight = nodeHeight - subIdx * maxCoreRowHeight;
        else
          thisNodeHeight = maxCoreRowHeight;

        dy -= subIdx * maxCoreRowHeight;
        dy -= thisNodeHeight / 2;

        if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
          nets << " : " << dx << " " << dy << endl;
        else
          nets << endl;
      } else {
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
    }
#endif
  }

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
  float maxWeight = 0;
  float weight;

  for (n = 0; n != getNumEdges(); ++n) {
    const HGFEdge& net = getEdgeByIdx(n);
    weight = net.getWeight();

    if (fabs(weight - 1.0) < 1e-10) continue;
    wts << " ";
    string netName = getNetNameByIndex(n).c_str();
    if (!netName.empty())
      wts << setw(10) << netName;
    else
      wts << "net" << n;

    wts << " " << setw(10) << weight << endl;
    if (maxWeight > weight) maxWeight = weight;
  }

  if (maxWeight == 0)
    maxWeight = 1;
  else
    maxWeight *= 1;

  // add edges to shredded cells
  unsigned netId = 0;
  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    if (getNodeHeight(nodeId) > maxCoreRowHeight && !isTerminal(nodeId))
    // shred these nodes
    {
      float nodeHeight = getNodeHeight(nodeId);
      //          float nodeWidth = getNodeWidth(nodeId);
      //          float dx;
      //          float dyAbove = -maxCoreRowHeight/2;
      //          float dyBelow = maxCoreRowHeight/2;

      string nodeName;
      string nodeNameNext;
      string netName;

      unsigned numNewNodes = unsigned(ceil(nodeHeight / maxCoreRowHeight));
      for (unsigned i = 0; i < numNewNodes - 1; ++i) {
        ostringstream ss1;
        ss1 << getNodeNameByIndex(node.getIndex()) << "@" << i;
        nodeName = ss1.str().c_str();
        ostringstream ss2;
        ss2 << getNodeNameByIndex(node.getIndex()) << "@" << i + 1;
        nodeNameNext = ss2.str().c_str();

        for (unsigned j = 0; j < NUM_NEW_VERT_NETS; ++j) {
          ostringstream ss3;
          ss3 << "temp@" << netId;
          netName = ss3.str().c_str();

          netId++;

          nets << "NetDegree : 2 ";
          nets << setw(10) << netName << endl;
          /*
             dx = j*(nodeWidth/(NUM_NEW_VERT_NETS-1));
             dx -= nodeWidth/2;

             nets<<setw(10)<<nodeName<<" B ";
             nets<<" : "<<dx <<" "  <<dyBelow<<endl;
             nets<<setw(10)<<nodeNameNext<<" B ";
             nets<<" : "<<dx <<" "  <<dyAbove<<endl;
           */
          nets << setw(10) << nodeName << " B " << endl;
          ;
          nets << setw(10) << nodeNameNext << " B " << endl;

          wts << setw(10) << netName;
          wts << " " << setw(10) << maxWeight << endl;
        }
      }
    }
  }

  nets.close();
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

// split in both dimensions
void HGraphWDimensions::saveAsNodesNetsWtsShredHW(const char* baseFileName,
                                                  float minWidthNode,
                                                  float maxCoreRowHeight,
                                                  unsigned singleNetWt) const {
  unsigned numNewSingleNets = 20;
  bool useSingleNet = false;
  if (singleNetWt > 0) {
    useSingleNet = true;
    numNewSingleNets = singleNetWt;
  }
  const string netsFileName = shredFileNameWithSuffix(baseFileName, ".nets");
  const string nodesFileName = shredFileNameWithSuffix(baseFileName, ".nodes");
  const string wtsFileName = shredFileNameWithSuffix(baseFileName, ".wts");
  const string auxFileName = shredFileNameWithSuffix(baseFileName, ".aux");

  ofstream nets(netsFileName.c_str());

  ofstream nodes(nodesFileName.c_str());

  ofstream wts(wtsFileName.c_str());

  ofstream aux(auxFileName.c_str());

  unsigned numExtraNodes = 0;
  unsigned numExtraNets = 0;
  unsigned numExtraPins = 0;

  itHGFNodeGlobal nIt;

  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    float nodeHeight = getNodeHeight(nodeId);
    float nodeWidth = getNodeWidth(nodeId);
    if (_isMacro[nodeId] == true && !isTerminal(nodeId)) {
      unsigned numNewVertNodes = unsigned(ceil(nodeHeight / maxCoreRowHeight));
      unsigned numNewHorizNodes = unsigned(ceil(nodeWidth / minWidthNode));
      numExtraNodes += (numNewVertNodes * numNewHorizNodes) - 1;
      if (_param.shredTopology == HGraphParameters::SINGLE_NET) {
        numExtraNets += numNewSingleNets;
        numExtraPins += numNewSingleNets * numNewVertNodes * numNewHorizNodes;
      } else if (_param.shredTopology == HGraphParameters::GRID) {
        unsigned numExtraNetsHere =
            (numNewVertNodes * (numNewHorizNodes - 1) * NUM_NEW_HORIZ_NETS +
             numNewHorizNodes * (numNewVertNodes - 1) * NUM_NEW_VERT_NETS);

        numExtraNets += numExtraNetsHere;
        numExtraPins += numExtraNetsHere * 2;
      } else if (_param.shredTopology == HGraphParameters::STAR) {
        unsigned numExtraNetsHere = numNewVertNodes * numNewHorizNodes;
        numExtraNets += numExtraNetsHere;
        numExtraPins += numExtraNetsHere * 2;
      }
    }
  }

  // write the nodes file and part of wts file

  wts << "UCLA wts 1.0" << endl;
  wts << TimeStamp() << User() << Platform() << endl;

  nodes << "UCLA nodes 1.0" << endl;
  nodes << TimeStamp() << User() << Platform() << endl;

  nodes << "NumNodes :     " << setw(10) << getNumNodes() + numExtraNodes
        << endl;
  nodes << "NumTerminals : " << setw(10) << getNumTerminals() << endl;

  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    Symmetry sym = _nodeSymmetries[nodeId];

    if (_isMacro[nodeId] == true && !isTerminal(nodeId))
    // shred these nodes
    {
      float nodeHeight = getNodeHeight(nodeId);
      float nodeWidth = getNodeWidth(nodeId);
      float remainingHeight = nodeHeight;
      float remainingWidth = nodeWidth;

      float thisNodeHeight = maxCoreRowHeight;
      float thisNodeWidth = minWidthNode;

      string nodeName;
      unsigned numNewVertNodes = unsigned(ceil(nodeHeight / maxCoreRowHeight));
      unsigned numNewHorizNodes = unsigned(ceil(nodeWidth / minWidthNode));

      for (unsigned i = 0; i < numNewVertNodes; ++i) {
        if (remainingHeight < maxCoreRowHeight)
          thisNodeHeight = remainingHeight;
        else
          thisNodeHeight = maxCoreRowHeight;

        remainingWidth = nodeWidth;
        for (unsigned j = 0; j < numNewHorizNodes; ++j) {
          ostringstream ss;
          ss << getNodeNameByIndex(node.getIndex()) << "@" << i << "#" << j;
          nodeName = ss.str().c_str();
          nodes << setw(10) << nodeName << "  ";

          if (remainingWidth < minWidthNode)
            thisNodeWidth = remainingWidth;
          else
            thisNodeWidth = minWidthNode;

          if (thisNodeWidth < 1e-5)  // float roundoff error
            thisNodeWidth = minWidthNode / 100;

          nodes << setw(10) << thisNodeWidth << " ";
          nodes << setw(10) << thisNodeHeight << " ";
          if (sym.rot90 || sym.y || sym.x) nodes << " : " << sym;
          nodes << endl;

          remainingWidth -= thisNodeWidth;

          wts << setw(10) << nodeName << "  ";
          wts << setw(10) << thisNodeWidth * thisNodeHeight << " " << endl;
        }
        remainingHeight -= thisNodeHeight;
      }
    } else {
      if (!getNodeNameByIndex(node.getIndex()).empty()) {
        nodes << setw(10) << getNodeNameByIndex(node.getIndex()) << "  ";
        wts << setw(10) << getNodeNameByIndex(node.getIndex()) << "  ";
      } else {
        if (isTerminal(nodeId)) {
          nodes << " p" << nodeId + 1 << " ";
          wts << setw(10) << " p" << nodeId + 1;
        } else {
          nodes << " a" << nodeId - getNumTerminals() << " ";
          wts << setw(10) << " a" << nodeId - getNumTerminals();
        }
      }

      if (getNodeWidth(nodeId) != 0 || getNodeHeight(nodeId) != 0)
        nodes << setw(10) << getNodeWidth(nodeId) << " " << setw(10)
              << getNodeHeight(nodeId) << " ";
      if (isTerminal(nodeId)) nodes << "terminal ";
      if (sym.rot90 || sym.y || sym.x) nodes << " : " << sym;
      nodes << endl;

      for (unsigned i = 0; i != getNumWeights(); ++i)
        wts << " " << setw(6) << getWeight(node.getIndex(), i);
      wts << endl;
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

  nets << "NumNets : " << getNumEdges() - num1PinEdges + numExtraNets << endl;
  nets << "NumPins : " << getNumPins() - num1PinEdges + numExtraPins << endl;

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
    nets << endl;

    unsigned nodeOffset = 0;

    itHGFNodeLocal v;
    float nodeHeight;
    float nodeWidth;
    unsigned subVertIdx;
    unsigned subHorizIdx;
    unsigned maxSubVertIdx;
    unsigned maxSubHorizIdx;

    float thisNodeHeight;
    float thisNodeWidth;
#ifdef SIGNAL_DIRECTIONS
    for (v = (*e)->srcsBegin(); v != (*e)->srcsEnd(); v++, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();

      string nodeName;
      nodeHeight = getNodeHeight(idx);
      nodeWidth = getNodeWidth(idx);
      if (_isMacro[idx] == true && !isTerminal(idx)) {
        maxSubVertIdx = unsigned(ceil(nodeHeight / maxCoreRowHeight) - 1);
        maxSubHorizIdx = unsigned(ceil(nodeWidth / minWidthNode) - 1);

        Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
        float dx = pOffset.x, dy = pOffset.y;

        if (dx < 0) dx = 0;
        if (dy < 0) dy = 0;

        if (fabs(dy) < 0.0001)
          subVertIdx = 0;
        else
          subVertIdx = unsigned(ceil(dy / maxCoreRowHeight) - 1);

        if (fabs(dx) < 0.0001)
          subHorizIdx = 0;
        else
          subHorizIdx = unsigned(ceil(dx / minWidthNode) - 1);

        if (subHorizIdx > maxSubHorizIdx) subHorizIdx = maxSubHorizIdx;

        if (subVertIdx > maxSubVertIdx) subVertIdx = maxSubVertIdx;

        ostringstream ss2;
        ss2 << getNodeNameByIndex(node.getIndex()) << "@" << subVertIdx << "#"
            << subHorizIdx;
        nodeName = ss2.str().c_str();
        nets << setw(10) << nodeName << " O ";

        if (subVertIdx == maxSubVertIdx)
          thisNodeHeight = nodeHeight - subVertIdx * maxCoreRowHeight;
        else
          thisNodeHeight = maxCoreRowHeight;

        if (subHorizIdx == maxSubHorizIdx)
          thisNodeWidth = nodeWidth - subHorizIdx * minWidthNode;
        else
          thisNodeWidth = minWidthNode;

        dy -= subVertIdx * maxCoreRowHeight;
        dy -= thisNodeHeight / 2;
        dx -= subHorizIdx * minWidthNode;
        dx -= thisNodeWidth / 2;

        if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
          nets << " : " << dx << " " << dy << endl;
        else
          nets << endl;
      } else {
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
    }
    for (v = (*e)->srcSnksBegin(); v != (*e)->srcSnksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();

      string nodeName;
      nodeHeight = getNodeHeight(idx);
      nodeWidth = getNodeWidth(idx);
      if (_isMacro[idx] == true && !isTerminal(idx)) {
        maxSubVertIdx = unsigned(ceil(nodeHeight / maxCoreRowHeight) - 1);
        maxSubHorizIdx = unsigned(ceil(nodeWidth / minWidthNode) - 1);

        Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
        float dx = pOffset.x, dy = pOffset.y;

        if (dx < 0) dx = 0;
        if (dy < 0) dy = 0;

        if (fabs(dy) < 0.0001)
          subVertIdx = 0;
        else
          subVertIdx = unsigned(ceil(dy / maxCoreRowHeight) - 1);

        if (fabs(dx) < 0.0001)
          subHorizIdx = 0;
        else
          subHorizIdx = unsigned(ceil(dx / minWidthNode) - 1);

        if (subHorizIdx > maxSubHorizIdx) subHorizIdx = maxSubHorizIdx;

        if (subVertIdx > maxSubVertIdx) subVertIdx = maxSubVertIdx;

        ostringstream ss3;
        ss3 << getNodeNameByIndex(node.getIndex()) << "@" << subVertIdx << "#"
            << subHorizIdx;
        nodeName = ss3.str().c_str();
        nets << setw(10) << nodeName << " B ";

        if (subVertIdx == maxSubVertIdx)
          thisNodeHeight = nodeHeight - subVertIdx * maxCoreRowHeight;
        else
          thisNodeHeight = maxCoreRowHeight;

        if (subHorizIdx == maxSubHorizIdx)
          thisNodeWidth = nodeWidth - subHorizIdx * minWidthNode;
        else
          thisNodeWidth = minWidthNode;

        dy -= subVertIdx * maxCoreRowHeight;
        dy -= thisNodeHeight / 2;
        dx -= subHorizIdx * minWidthNode;
        dx -= thisNodeWidth / 2;

        if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
          nets << " : " << dx << " " << dy << endl;
        else
          nets << endl;
      } else {
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
    }
    for (v = (*e)->snksBegin(); v != (*e)->snksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();

      string nodeName;
      nodeHeight = getNodeHeight(idx);
      nodeWidth = getNodeWidth(idx);
      if (_isMacro[idx] == true && !isTerminal(idx)) {
        maxSubVertIdx = unsigned(ceil(nodeHeight / maxCoreRowHeight) - 1);
        maxSubHorizIdx = unsigned(ceil(nodeWidth / minWidthNode) - 1);

        Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
        float dx = pOffset.x, dy = pOffset.y;

        if (dx < 0) dx = 0;
        if (dy < 0) dy = 0;

        if (fabs(dy) < 0.0001)
          subVertIdx = 0;
        else
          subVertIdx = unsigned(ceil(dy / maxCoreRowHeight) - 1);

        if (fabs(dx) < 0.0001)
          subHorizIdx = 0;
        else
          subHorizIdx = unsigned(ceil(dx / minWidthNode) - 1);

        if (subHorizIdx > maxSubHorizIdx) subHorizIdx = maxSubHorizIdx;

        if (subVertIdx > maxSubVertIdx) subVertIdx = maxSubVertIdx;

        ostringstream ss4;
        ss4 << getNodeNameByIndex(node.getIndex()) << "@" << subVertIdx << "#"
            << subHorizIdx;
        nodeName = ss4.str().c_str();
        nets << setw(10) << nodeName << " I ";

        if (subVertIdx == maxSubVertIdx)
          thisNodeHeight = nodeHeight - subVertIdx * maxCoreRowHeight;
        else
          thisNodeHeight = maxCoreRowHeight;

        if (subHorizIdx == maxSubHorizIdx)
          thisNodeWidth = nodeWidth - subHorizIdx * minWidthNode;
        else
          thisNodeWidth = minWidthNode;

        dy -= subVertIdx * maxCoreRowHeight;
        dy -= thisNodeHeight / 2;
        dx -= subHorizIdx * minWidthNode;
        dx -= thisNodeWidth / 2;

        if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
          nets << " : " << dx << " " << dy << endl;
        else
          nets << endl;
      } else {
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
    }
#else
    for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();

      string nodeName;
      nodeHeight = getNodeHeight(idx);
      nodeWidth = getNodeWidth(idx);
      if (_isMacro[idx] == true && !isTerminal(idx)) {
        maxSubVertIdx = unsigned(ceil(nodeHeight / maxCoreRowHeight) - 1);
        maxSubHorizIdx = unsigned(ceil(nodeWidth / minWidthNode) - 1);

        Point pOffset = nodesOnEdgePinOffset(nodeOffset, edgeId);
        float dx = pOffset.x, dy = pOffset.y;

        if (dx < 0) dx = 0;
        if (dy < 0) dy = 0;

        if (fabs(dy) < 0.0001)
          subVertIdx = 0;
        else
          subVertIdx = unsigned(ceil(dy / maxCoreRowHeight) - 1);

        if (fabs(dx) < 0.0001)
          subHorizIdx = 0;
        else
          subHorizIdx = unsigned(ceil(dx / minWidthNode) - 1);

        if (subHorizIdx > maxSubHorizIdx) subHorizIdx = maxSubHorizIdx;

        if (subVertIdx > maxSubVertIdx) subVertIdx = maxSubVertIdx;

        ostringstream ss3;
        ss3 << getNodeNameByIndex(node.getIndex()) << "@" << subVertIdx << "#"
            << subHorizIdx;
        nodeName = ss3.str().c_str();
        nets << setw(10) << nodeName << " B ";

        if (subVertIdx == maxSubVertIdx)
          thisNodeHeight = nodeHeight - subVertIdx * maxCoreRowHeight;
        else
          thisNodeHeight = maxCoreRowHeight;

        if (subHorizIdx == maxSubHorizIdx)
          thisNodeWidth = nodeWidth - subHorizIdx * minWidthNode;
        else
          thisNodeWidth = minWidthNode;

        dy -= subVertIdx * maxCoreRowHeight;
        dy -= thisNodeHeight / 2;
        dx -= subHorizIdx * minWidthNode;
        dx -= thisNodeWidth / 2;

        if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
          nets << " : " << dx << " " << dy << endl;
        else
          nets << endl;
      } else {
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
    }
#endif
  }

  // write the wts for nets
  float maxWeight = 0;
  float weight;
  int n;
  // non fake nets should have 10 times weight of fake nets
  for (n = 0; n != int(getNumEdges()); ++n) {
    const HGFEdge& net = getEdgeByIdx(n);
    weight = net.getWeight();
    if (maxWeight < weight) maxWeight = weight;

    if (fabs(weight - 1.0) < 1e-10) continue;
    wts << " ";
    const char* netName = getNetNameByIndex(n).c_str();
    if (netName)
      wts << setw(10) << netName;
    else
      wts << "net" << n;

    // wts << " " << setw(10) << 10*weight << endl;
    wts << " " << setw(10) << weight << endl;
  }

  if (maxWeight == 0)
    maxWeight = 1;
  else
    maxWeight *= 100;

  // add edges to shredded cells
  if (_param.shredTopology == HGraphParameters::GRID)
    AddNetsToShredsGridTopology(nets, wts, maxCoreRowHeight, minWidthNode,
                                maxWeight);
  else if (_param.shredTopology == HGraphParameters::STAR)
    AddNetsToShredsStarTopology(nets, wts, maxCoreRowHeight, minWidthNode,
                                maxWeight);
  else if (_param.shredTopology == HGraphParameters::SINGLE_NET)
    AddNetsToShredsSingleNet(nets, wts, maxCoreRowHeight, minWidthNode,
                             maxWeight, numNewSingleNets);
  else
    abkfatal(0, "Unknown Shred Net Topology.\n");

  nets.close();
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

void HGraphWDimensions::AddNetsToShredsStarTopology(ofstream& nets,
                                                    ofstream& wts,
                                                    float maxCoreRowHeight,
                                                    float minWidthNode,
                                                    double maxWeight) const {
  unsigned netId = 0;
  itHGFNodeGlobal nIt;
  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    if (_isMacro[nodeId] == true && !isTerminal(nodeId))
    // then these nodes are shreds
    {
      float nodeHeight = getNodeHeight(nodeId);
      float nodeWidth = getNodeWidth(nodeId);

      string StarNodeName;
      string nodeName;
      string nodeNameNext;
      string netName;

      unsigned numNewVertNodes = unsigned(ceil(nodeHeight / maxCoreRowHeight));
      unsigned numNewHorizNodes = unsigned(ceil(nodeWidth / minWidthNode));

      ostringstream starPt_ss;
      starPt_ss << getNodeNameByIndex(node.getIndex()) << "@" << 0 << "#" << 0;
      StarNodeName = starPt_ss.str().c_str();

      for (unsigned i = 0; i < numNewVertNodes; ++i) {
        for (unsigned j = 0; j < numNewHorizNodes; ++j) {
          ostringstream ss2;
          ss2 << getNodeNameByIndex(node.getIndex()) << "@" << i << "#" << j;
          nodeNameNext = ss2.str().c_str();

          ostringstream ss3;
          ss3 << "temp@" << netId;
          netName = ss3.str().c_str();
          netId++;

          nets << "NetDegree : 2 ";
          nets << setw(10) << netName << endl;

          nets << setw(10) << StarNodeName << " B " << endl;
          ;
          nets << setw(10) << nodeNameNext << " B " << endl;

          wts << setw(10) << netName;
          wts << " " << setw(10) << maxWeight << endl;
        }
      }
    }
  }
}

void HGraphWDimensions::AddNetsToShredsGridTopology(ofstream& nets,
                                                    ofstream& wts,
                                                    float maxCoreRowHeight,
                                                    float minWidthNode,
                                                    double maxWeight) const {
  unsigned netId = 0;
  itHGFNodeGlobal nIt;
  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    if (_isMacro[nodeId] == true && !isTerminal(nodeId))
    // then these nodes are shreds
    {
      float nodeHeight = getNodeHeight(nodeId);
      float nodeWidth = getNodeWidth(nodeId);

      string nodeName;
      string nodeNameNext;
      string netName;

      unsigned numNewVertNodes = unsigned(ceil(nodeHeight / maxCoreRowHeight));
      unsigned numNewHorizNodes = unsigned(ceil(nodeWidth / minWidthNode));

      for (unsigned i = 0; i < numNewVertNodes; ++i) {
        for (unsigned j = 0; j < numNewHorizNodes; ++j) {
          if (j != numNewHorizNodes - 1) {
            ostringstream ss;
            ss << getNodeNameByIndex(node.getIndex()) << "@" << i << "#" << j;
            nodeName = ss.str().c_str();
            ostringstream ss2;
            ss2 << getNodeNameByIndex(node.getIndex()) << "@" << i << "#"
                << j + 1;
            nodeNameNext = ss2.str().c_str();

            for (unsigned k = 0; k < NUM_NEW_HORIZ_NETS; ++k) {
              ostringstream ss3;
              ss3 << "temp@" << netId;
              netName = ss3.str().c_str();
              netId++;

              nets << "NetDegree : 2 ";
              nets << setw(10) << netName << endl;

              nets << setw(10) << nodeName << " B " << endl;
              ;
              nets << setw(10) << nodeNameNext << " B " << endl;

              wts << setw(10) << netName;
              wts << " " << setw(10) << maxWeight << endl;
            }
          }

          if (i != numNewVertNodes - 1) {
            ostringstream ss;
            ss << getNodeNameByIndex(node.getIndex()) << "@" << i << "#" << j;
            nodeName = ss.str().c_str();
            ostringstream ss2;
            ss2 << getNodeNameByIndex(node.getIndex()) << "@" << i + 1 << "#"
                << j;
            nodeNameNext = ss2.str().c_str();

            for (unsigned k = 0; k < NUM_NEW_VERT_NETS; ++k) {
              ostringstream ss45;
              ss45 << "temp@" << netId;

              netName = ss45.str().c_str();
              netId++;

              nets << "NetDegree : 2 ";
              nets << setw(10) << netName << endl;
              nets << setw(10) << nodeName << " B " << endl;
              ;
              nets << setw(10) << nodeNameNext << " B " << endl;

              wts << setw(10) << netName;
              wts << " " << setw(10) << maxWeight << endl;
            }
          }
        }
      }
    }
  }
}

// add single nets instead of distributed nets
void HGraphWDimensions::AddNetsToShredsSingleNet(
    ofstream& nets, ofstream& wts, float maxCoreRowHeight, float minWidthNode,
    double maxWeight, unsigned numNewSingleNets) const {
  unsigned netId = 0;
  itHGFNodeGlobal nIt;
  for (nIt = nodesBegin(); nIt != nodesEnd(); nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    if (_isMacro[nodeId] == true && !isTerminal(nodeId))
    // then these nodes are shreds
    {
      float nodeHeight = getNodeHeight(nodeId);
      float nodeWidth = getNodeWidth(nodeId);
      string nodeName;
      string netName;

      unsigned numNewVertNodes = unsigned(ceil(nodeHeight / maxCoreRowHeight));
      unsigned numNewHorizNodes = unsigned(ceil(nodeWidth / minWidthNode));
      for (unsigned l = 0; l < numNewSingleNets; ++l) {
        ostringstream ss3;
        ss3 << "temp@" << netId;
        netName = ss3.str().c_str();
        netId++;
        nets << "NetDegree :  " << (numNewVertNodes * numNewHorizNodes);
        nets << " " << setw(10) << netName << endl;
        wts << setw(10) << netName;
        wts << " " << setw(10) << maxWeight << endl;
        for (unsigned i = 0; i < numNewVertNodes; ++i) {
          for (unsigned j = 0; j < numNewHorizNodes; ++j) {
            ostringstream ss;
            ss << getNodeNameByIndex(node.getIndex()) << "@" << i << "#" << j;
            nodeName = ss.str().c_str();
            nets << setw(10) << nodeName << " B " << endl;
            ;
          }
        }
      }
    }
  }
}
