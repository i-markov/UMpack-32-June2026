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

// Mini-benchmark slicing driver.
// It marks cells that cross the cut window, removes fully fixed nets, and
// writes a reduced benchmark for placement experimentation.

/*
 * <aaronnn> create a new benchmark by cutting up a region of an existing
 * benchmark. NOTE: add obstacles by hand to mask excluded regions using capo's
 * command-line option.
 *
 * For a less hackish way to do this, construct RBPlacement and hgraph from
 * scratch using API.
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/infolines.h"
#include "ABKCommon/paramproc.h"
#include "Geoms/bbox.h"
#include "PlaceEvals/pinPlEval.h"
#include "RBPlacement.h"

void sliceAndSaveBenchmark(const RBPlacement& rbplace, const BBox& cutWindow);

int main(int argc, const char* argv[]) {
  RBPlacement::Parameters rbParams(argc, argv);

  StringParam auxFile("f", argc, argv);
  DoubleParam cutXMin_("cutXMin", argc, argv);
  DoubleParam cutXMax_("cutXMax", argc, argv);
  DoubleParam cutYMin_("cutYMin", argc, argv);
  DoubleParam cutYMax_("cutYMax", argc, argv);

  abkfatal(auxFile.found(), "must have -f <auxfilename>");

  if (!cutXMin_.found() && !cutXMax_.found())
    abkfatal(0, "no X cut(s); specify -cutXMin and -cutXMax");
  if (!cutYMin_.found() && !cutYMax_.found())
    abkfatal(0, "no Y cut(s); specify -cutYMin and -cutYMax");

  RBPlacement rbplace(auxFile, rbParams);

  unsigned getTerminals = true;
  BBox origCoreBBox = rbplace.getBBox(getTerminals);

  double cutXMin = origCoreBBox.xMin;
  double cutXMax = origCoreBBox.xMax;
  double cutYMin = origCoreBBox.yMin;
  double cutYMax = origCoreBBox.yMax;

  if (cutXMin_.found()) cutXMin = cutXMin_;
  if (cutXMax_.found()) cutXMax = cutXMax_;
  if (cutYMin_.found()) cutYMin = cutYMin_;
  if (cutYMax_.found()) cutYMax = cutYMax_;

  sliceAndSaveBenchmark(rbplace, BBox(cutXMin, cutYMin, cutXMax, cutYMax));

  return 0;
}

void sliceAndSaveBenchmark(const RBPlacement& rbplace, const BBox& cutWindow) {
  const HGraphWDimensions& hgraph = rbplace.getHGraph();

  cout << "Cutting up a benchmark with core region " << cutWindow << endl;

  vector<bool> nodesToFix(hgraph.getNumNodes(), false);
  vector<bool> edgesToRemove(hgraph.getNumEdges(), false);
  unsigned numExtraNodesFixed = 0;
  unsigned numEdgesToRemove = 0;
  unsigned numPinsToRemove = 0;

  // mark nodes to fix as terminals (but only those that aren't already fixed)
  for (unsigned nodeID = 0; nodeID < hgraph.getNumNodes(); nodeID++) {
    const HGFNode& node = hgraph.getNodeByIdx(nodeID);

    // if any pins are outside the window, fix it
    unsigned edgeOffset = 0;
    for (itHGFEdgeLocal edgeIt = node.edgesBegin(); edgeIt != node.edgesEnd();
         edgeIt++, edgeOffset++) {
      const HGFEdge& edge = **edgeIt;
      unsigned edgeID = edge.getIndex();

      if (hgraph.isTerminal(nodeID) || rbplace.isFixed(nodeID)) {
        // node is already a terminal or fixed, don't fix it
        continue;
      }

      // get relative pin offset
      Point pinPoint =
          hgraph.edgesOnNodePinOffsetBBox(edgeOffset, nodeID).getGeomCenter();

      // get absolute pin position
      pinPoint += rbplace.getPlacement()[nodeID];

      if (cutWindow.contains(pinPoint)) continue;

      // this cell has a pin outside the cut window
      // and isn't a terminal, mark this cell for fixing as a terminal
      nodesToFix[nodeID] = true;
      numExtraNodesFixed++;

      break;
    }
  }

  // mark edges to remove
  for (itHGFEdgeGlobal edgeIt = hgraph.edgesBegin();
       edgeIt != hgraph.edgesEnd(); edgeIt++) {
    const HGFEdge& edge = **edgeIt;
    unsigned edgeID = edge.getIndex();

    itHGFNodeLocal node;
    for (node = edge.nodesBegin(); node != edge.nodesEnd(); node++) {
      unsigned nodeID = (*node)->getIndex();

      if (!hgraph.isTerminal(nodeID) && !rbplace.isFixed(nodeID) &&
          !nodesToFix[nodeID])
        break;
    }

    if (node == edge.nodesEnd()) {
      // all the nodes are fixed
      // -- remove this edge
      edgesToRemove[edgeID];
      numEdgesToRemove++;
      numPinsToRemove += edge.getDegree();
    }
  }

  // write the nodes file
  cout << "Saving out-mini.nodes" << endl;
  ofstream nodes("out-mini.nodes");
  nodes << "UCLA nodes 1.0" << endl;
  nodes << TimeStamp() << User() << Platform() << endl;

  nodes << "NumNodes :     " << setw(10) << hgraph.getNumNodes() << endl;
  nodes << "NumTerminals : " << setw(10)
        << hgraph.getNumTerminals() + numExtraNodesFixed << endl;

  for (itHGFNodeGlobal nIt = hgraph.nodesBegin(); nIt != hgraph.nodesEnd();
       nIt++) {
    HGFNode& node = **nIt;
    unsigned nodeId = node.getIndex();
    Symmetry sym = hgraph.getNodeSymmetry(nodeId);

    if (!node.getName().empty())
      nodes << setw(10) << node.getName() << "  ";
    else {
      if (hgraph.isTerminal(nodeId))
        nodes << " p" << nodeId + 1 << " ";
      else
        nodes << " a" << nodeId - hgraph.getNumTerminals() << " ";
    }

    if (hgraph.getNodeWidth(nodeId) != 0 || hgraph.getNodeHeight(nodeId) != 0)
      nodes << setw(10) << hgraph.getNodeWidth(nodeId) << " " << setw(10)
            << hgraph.getNodeHeight(nodeId) << " ";
    if (hgraph.isTerminal(nodeId))
      nodes << "terminal ";
    else if (nodesToFix[nodeId])
      nodes << "terminal ";

    if (sym.rot90 || sym.y || sym.x) nodes << " : " << sym;
    nodes << endl;
  }
  nodes.close();

  unsigned num1PinEdges = 0, num0PinEdges = 0;
  for (itHGFEdgeGlobal e = hgraph.edgesBegin(); e != hgraph.edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    if (edge.getDegree() > 1) continue;
    if (edge.getDegree() == 1)
      num1PinEdges++;
    else
      num0PinEdges++;
  }

  // write the nets file
  cout << "Saving out-mini.nets" << endl;
  ofstream nets("out-mini.nets");
  nets << "UCLA nets 1.0" << endl;
  nets << TimeStamp() << User() << Platform() << endl;

  nets << "NumNets : " << hgraph.getNumEdges() - num1PinEdges - numEdgesToRemove
       << endl;
  nets << "NumPins : " << hgraph.getNumPins() - num1PinEdges - numPinsToRemove
       << endl;

  for (itHGFEdgeGlobal e = hgraph.edgesBegin(); e != hgraph.edgesEnd(); ++e) {
    HGFEdge& edge = **e;
    if (edge.getDegree() < 2) continue;
    unsigned edgeId = edge.getIndex();
    string netName = hgraph.getNetNameByIndex(edgeId);

    // remove net we marked earlier (the edge only spans fixed nodes)
    if (edgesToRemove[edgeId]) continue;

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
    for (v = (*e)->srcsBegin(); v != (*e)->srcsEnd(); v++, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << node.getName() << " O ";

      Point pOffset = hgraph.nodesOnEdgePinOffset(nodeOffset, edgeId);
      double dx = pOffset.x - hgraph.getNodeWidth(idx) / 2,
             dy = pOffset.y - hgraph.getNodeHeight(idx) / 2;
      //	    if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;
      else
        nets << endl;
    }
    for (v = (*e)->srcSnksBegin(); v != (*e)->srcSnksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << node.getName() << " B ";

      Point pOffset = hgraph.nodesOnEdgePinOffset(nodeOffset, edgeId);
      double dx = pOffset.x - hgraph.getNodeWidth(idx) / 2,
             dy = pOffset.y - hgraph.getNodeHeight(idx) / 2;
      //          if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;
      else
        nets << endl;
    }
    for (v = (*e)->snksBegin(); v != (*e)->snksEnd(); ++v, nodeOffset++) {
      const HGFNode& node = (**v);
      unsigned idx = node.getIndex();
      nets << setw(10) << node.getName() << " I ";

      Point pOffset = hgraph.nodesOnEdgePinOffset(nodeOffset, edgeId);
      double dx = pOffset.x - hgraph.getNodeWidth(idx) / 2,
             dy = pOffset.y - hgraph.getNodeHeight(idx) / 2;

      //	    if(pOffset.x != 0 || pOffset.y != 0)
      if (fabs(dx) > 1e-10 || fabs(dy) > 1e-10)
        nets << " : " << dx << " " << dy << endl;

      else
        nets << endl;
    }
  }
  nets.close();

  cout << "Done. Remember to specify an obstacle file for Metaplacer to remove "
          "sites"
       << endl;
}
