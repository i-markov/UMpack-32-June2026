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

// Created by David Papa on 06/28/04

#include <sstream>
#include <string>

#include "ABKCommon/abkassert.h"
#include "hgWDims.h"

using std::string;
using std::stringstream;
using std::vector;

// added by dp for easy manual construction
// call createNewManualHGraphWDimensions to get an object to work with
// Instructions for manually constructing a HGraphWDimensions:
// for( all edges )
//   if(not all nodes on this edge are created)
//     create missing nodes
//   Create an edge
//   Add all nodes on the edge
// finalize()

HGraphWDimensions* HGraphWDimensions::createNewManualHGraphWDimensions(
    unsigned numNodes) {
  unsigned numWeights = 1;
  HGraphParameters defaultParams;
  vector<float> heights(numNodes, 1.0);
  vector<float> widths(numNodes, 1.0);
  HGraphWDimensions* h = new HGraphWDimensions(numNodes, numWeights,
                                               defaultParams, heights, widths);
  h->_nodeSymmetries = vector<Symmetry>(numNodes);
  h->_numTerminals = 0;
  return h;
}

unsigned HGraphWDimensions::create_edge() { return addEdge()->getIndex(); }

unsigned HGraphWDimensions::create_edge(std::string name) {
  unsigned edgeRVal = addEdge()->getIndex();
  setNetNameManual(edgeRVal, name.c_str());
  return edgeRVal;
}

unsigned HGraphWDimensions::create_node(std::string name, size_t width,
                                        size_t height, bool term) {
  int newNodeIdx = -1;
  if (term) ++_numTerminals;
  if (!haveSuchNode(name.c_str())) {
    newNodeIdx = getNodeNames().size();
    getNodeNames().resize(getNodeNames().size() + 1);
    getNodeNames().back() = name.c_str();
    getNodeNamesMap()[getNodeNames().back()] = newNodeIdx;
    abkfatal(_nodes.size() > static_cast<unsigned>(newNodeIdx),
             "vector index too large");
    _nodes[newNodeIdx] = (new HGFNode(newNodeIdx));
    setWeight(newNodeIdx, 1.0);
    setNodeWidth(width, newNodeIdx);
    setNodeHeight(height, newNodeIdx);

    // nonTermNodeNames.push_back(name.c_str());
    // termNodeNames.push_back(name.c_str());
  }
  HGFNode& node = getNodeByName(name.c_str());
  abkfatal(newNodeIdx >= 0, "broken create_node");
  abkfatal(node.getIndex() == static_cast<unsigned>(newNodeIdx),
           "broken create_node");
  return node.getIndex();
}

unsigned HGraphWDimensions::create_node(size_t width, size_t height,
                                        bool term) {
  if (term) ++_numTerminals;
  stringstream ss;
  ss << "dp" << getNumNodes();
  string name = ss.str().c_str();
  unsigned newNodeIdx = getNodeNames().size();
  getNodeNames().resize(getNodeNames().size() + 1);
  getNodeNames().back() = name.c_str();
  getNodeNamesMap()[getNodeNames().back()] = newNodeIdx;
  HGFNode* nn = new HGFNode(newNodeIdx);
  setWeight(nn->getIndex(), 1.0);
  _nodes[newNodeIdx] = nn;
  setNodeWidth(width, newNodeIdx);
  setNodeHeight(height, newNodeIdx);
  return newNodeIdx;
}

void HGraphWDimensions::add_node_to_edge(unsigned node_idx, unsigned edge_idx,
                                         char direction) {
  Point offset;
  offset.x = _nodeWidths[node_idx] / 2.0;
  offset.y = _nodeHeights[node_idx] / 2.0;

  switch (direction) {
    case 'I':
#ifdef SIGNAL_DIRECTIONS
      addSnkWPin(getNodeByIdx(node_idx), getEdgeByIdx(edge_idx));
      break;
#endif
    case 'O':
#ifdef SIGNAL_DIRECTIONS
      addSrcWPin(getNodeByIdx(node_idx), getEdgeByIdx(edge_idx));
      break;
#endif
    case 'B':
      addSrcSnkWPin(getNodeByIdx(node_idx), getEdgeByIdx(edge_idx));
      break;
    default:
      abkfatal(0, "no such pin direction");
  }
  _numPins++;
}

void HGraphWDimensions::done(void) {
  finalize();
  setupUnitDimensions();
  setupCenterPinOffsets();
  setupEmptySymmetries();
}
