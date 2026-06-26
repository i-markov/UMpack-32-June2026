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

#include "hgForHierarchy.h"

#include <string>

using std::string;
using std::vector;

// oNodes is a vector of the nodes to make a subHGraph of.
// all other vectors (nodeWeights, adjacentEdges) are indexed
// BY THE VALUES IN oNodes.

BitBoard HGraphForHierarchy::_madeEdge(1);  // have to give it some size
vector<unsigned> HGraphForHierarchy::_clusteredDegree;

HGraphForHierarchy::HGraphForHierarchy(const vector<unsigned>& oNodes,
                                       unsigned numTerminals,
                                       const FillableHierarchy& fillH)
    : SubHGraph(numTerminals, oNodes.size() - numTerminals,
                fillH.getEdgeWeights().size()) {
  //    _nodeNames = vector<string>(_nodes.size());

  const vector<double>& edgeWeights = fillH.getEdgeWeights();
  const vector<vector<unsigned> >& adjacentEdges = fillH.getAdjacentEdges();
  const vector<double>& nodeWeights = fillH.getClusterAreas();

  origNodeToNew.reserve(oNodes.size());
  origEdgeToNew.reserve(edgeWeights.size());

  unsigned n;
  for (n = 0; n < oNodes.size(); n++) {  // similar to the addNode function
    unsigned oNodeId = oNodes[n];

    setWeight(n, nodeWeights[oNodeId]);

    origNodeToNew[oNodeId] = n;
    newNodeToOrig[n] = oNodeId;
  }

  // if the size is correct, this just calls clear. Otherwise,
  // it allocates the correct amount of space 1st
  _madeEdge.reset(edgeWeights.size());

  if (_clusteredDegree.size() < edgeWeights.size())
    _clusteredDegree = vector<unsigned>(edgeWeights.size(), 0);
  else
    std::fill(_clusteredDegree.begin(), _clusteredDegree.end(), 0);

  // count the edge degree, so we can exclude degree 1 edges
  for (n = 0; n < oNodes.size(); n++) {
    unsigned oNodeId = oNodes[n];
    for (unsigned e = 0; e < adjacentEdges[oNodeId].size(); e++)
      _clusteredDegree[adjacentEdges[oNodeId][e]]++;
  }

  for (n = 0; n < oNodes.size(); n++) {
    unsigned oNodeId = oNodes[n];
    const HGFNode& node = *_nodes[n];

    for (unsigned e = 0; e < adjacentEdges[oNodeId].size(); e++) {
      unsigned eId = adjacentEdges[oNodeId][e];
      if (_clusteredDegree[eId] < 2) continue;

      if (!_madeEdge.isBitSet(eId)) {
        _madeEdge.setBit(eId);
        HGFEdge* newEdge = addEdge(edgeWeights[eId]);

        origEdgeToNew[eId] = newEdge->getIndex();
        if (newEdgeToOrig)
          (*newEdgeToOrig)[newEdge->getIndex()] = eId;
        else
          abkfatal(0,
                   "Filled heirarchy must have a newEdgeToOrig, and hence a "
                   "net bound.");
      }

      addSrcSnk(node, getNewEdgeByOrigIdx(eId));
    }
  }

  finalize();
}
