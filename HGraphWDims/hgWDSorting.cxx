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

#include "ABKCommon/abkassert.h"
#include "Stats/stats.h"
#include "hgWDims.h"

using std::sort;
using std::vector;

struct HGWDimsNodeOnEdgeIdxComp {
  const ctainerHGFNodesLocal* localNds;

  bool operator()(unsigned nO1, unsigned nO2) {
    return (*localNds)[nO1]->getIndex() < (*localNds)[nO2]->getIndex();
  }
};

struct HGWDimsNodeOnEdgeAscDegreeComp {
  const ctainerHGFNodesLocal* localNds;

  bool operator()(unsigned nO1, unsigned nO2) {
    return (*localNds)[nO1]->getDegree() < (*localNds)[nO2]->getDegree();
  }
};

struct HGWDimsNodeOnEdgeDesDegreeComp {
  const ctainerHGFNodesLocal* localNds;

  bool operator()(unsigned nO1, unsigned nO2) {
    return (*localNds)[nO1]->getDegree() > (*localNds)[nO2]->getDegree();
  }
};

void HGraphWDimensions::sortNodes() {
  switch (_param.nodeSortMethod) {
    case (HGraphParameters::DONT_SORT_NODES):
      return;
    case (HGraphParameters::SORT_NODES_BY_INDEX):
      break;
    case (HGraphParameters::SORT_NODES_BY_ASCENDING_DEGREE):
    case (HGraphParameters::SORT_NODES_BY_DESCENDING_DEGREE):
      abkfatal(0, "all node sorts but index are disabled");
      break;
    default:
      abkfatal(0, "unknown nodeSortMethod value in Parameters");
  }

  HGWDimsNodeOnEdgeIdxComp compFun;

  for (itHGFEdgeGlobal edgeIt = edgesBegin(); edgeIt != edgesEnd(); ++edgeIt) {
    HGFEdge& edge = **edgeIt;
    compFun.localNds = &edge._nodes;

    vector<unsigned> nOffsets(edge._nodes.size());

    for (unsigned c = 0; c < nOffsets.size(); ++c) {
      nOffsets[c] = c;
    }

    vector<unsigned>::iterator rangeBgn, rangeEnd;
#ifdef SIGNAL_DIRECTIONS
    rangeBgn = rangeEnd = nOffsets.begin();
    rangeEnd += edge.getNumSrcs();
    sort(rangeBgn, rangeEnd, compFun);

    rangeBgn = rangeEnd;
    rangeEnd += edge.getNumSrcSnks();
    sort(rangeBgn, rangeEnd, compFun);

    rangeBgn = rangeEnd;
    rangeEnd = nOffsets.end();
    sort(rangeBgn, rangeEnd, compFun);
#else
    rangeBgn = rangeEnd = nOffsets.begin();
    rangeEnd += edge.getDegree();
    sort(rangeBgn, rangeEnd, compFun);
#endif

    // take nOffsets as a permutation, and apply it to
    // the edge's _node vector and the pin-offsets vector

    ctainerHGFNodesLocal nodesCpy(edge._nodes);
    vector<unsigned> pinOffsetCpy;
    unsigned pOffsetBgn = _noeBeginEnd[edge.getIndex()];
    unsigned pOffsetEnd = _noeBeginEnd[edge.getIndex() + 1];

    pinOffsetCpy.insert(pinOffsetCpy.begin(),
                        _nodesOnEdgePinIdx.begin() + pOffsetBgn,
                        _nodesOnEdgePinIdx.begin() + pOffsetEnd);

    for (unsigned n = 0; n < nOffsets.size(); n++) {
      edge._nodes[n] = nodesCpy[nOffsets[n]];
      _nodesOnEdgePinIdx[pOffsetBgn + n] = pinOffsetCpy[nOffsets[n]];
    }
  }
}

struct HGWDimsEdgeOnNodeIdxComp {
  const ctainerHGFEdgesLocal* localEgs;

  bool operator()(unsigned nO1, unsigned nO2) {
    return (*localEgs)[nO1]->getIndex() < (*localEgs)[nO2]->getIndex();
  }
};

struct HGWDimsEdgeOnNodeAscDegreeComp {
  const ctainerHGFEdgesLocal* localEgs;

  bool operator()(unsigned nO1, unsigned nO2) {
    return (*localEgs)[nO1]->getDegree() < (*localEgs)[nO2]->getDegree();
  }
};

struct HGWDimsEdgeOnNodeDesDegreeComp {
  const ctainerHGFEdgesLocal* localEgs;

  bool operator()(unsigned nO1, unsigned nO2) {
    return (*localEgs)[nO1]->getDegree() > (*localEgs)[nO2]->getDegree();
  }
};

void HGraphWDimensions::sortEdges() {
  HGWDimsEdgeOnNodeIdxComp compFun;

  switch (_param.edgeSortMethod) {
    case (HGraphParameters::SORT_EDGES_BY_INDEX):
      break;
    case (HGraphParameters::DONT_SORT_EDGES):
      return;
    case (HGraphParameters::SORT_EDGES_BY_ASCENDING_DEGREE):
      abkfatal(0, "SORT_EDGES_BY_ASCENDING_DEGREE is disabled");
    case (HGraphParameters::SORT_EDGES_BY_DESCENDING_DEGREE):
      abkfatal(0, "SORT_EDGES_BY_DESCENTING_DEGREE is disabled");
    default:
      abkfatal(0, "unknown edgeSortMethod value in Parameters");
  };

  vector<unsigned> edgesOnNodePinIdxCpy = _edgesOnNodePinIdx;

  for (itHGFNodeGlobal nodeIt = nodesBegin(); nodeIt != nodesEnd(); ++nodeIt) {
    HGFNode& node = **nodeIt;
    if (node.getDegree() == 0) continue;

    compFun.localEgs = &node._edges;

    vector<unsigned> eOffsets(node._edges.size());

    for (unsigned c = 0; c < eOffsets.size(); c++) {
      eOffsets[c] = c;
    }

    vector<unsigned>::iterator rangeBgn, rangeEnd;
#ifdef SIGNAL_DIRECTIONS
    rangeBgn = rangeEnd = eOffsets.begin();
    rangeEnd += node.getNumSrcEdges();
    sort(rangeBgn, rangeEnd, compFun);

    rangeBgn = rangeEnd;
    rangeEnd += node.getNumSrcSnkEdges();
    sort(rangeBgn, rangeEnd, compFun);

    rangeBgn = rangeEnd;
    rangeEnd = eOffsets.end();
    sort(rangeBgn, rangeEnd, compFun);
#else
    rangeBgn = rangeEnd = eOffsets.begin();
    rangeEnd += node.getDegree();
    sort(rangeBgn, rangeEnd, compFun);
#endif

    // take eOffsets as a permutation, and apply it to
    // the node's _edges vector and the pin-offsets vector

    ctainerHGFEdgesLocal edgesCpy(node._edges);
    unsigned pOffsetBgn = _eonBeginEnd[node.getIndex()];

    for (unsigned n = 0; n < eOffsets.size(); n++) {
      node._edges[n] = edgesCpy[eOffsets[n]];
      _edgesOnNodePinIdx[pOffsetBgn + n] =
          edgesOnNodePinIdxCpy[pOffsetBgn + eOffsets[n]];
    }
  }
}

void HGraphWDimensions::sortAsDB() {
  _param.edgeSortMethod = HGraphParameters::SORT_EDGES_BY_INDEX;
  sortEdges();
}
