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

#include "constrainingHG.h"

#include <string>
#include <vector>

using std::string;
using std::vector;

ConstrainingHGraph::ConstrainingHGraph(const vector<CapoBin*>& allBins,
                                       const HGraphFixed& hgraphOfNetlist,
                                       const Verbosity verb)
    : HGraphFixed(hgraphOfNetlist.getNumNodes() + allBins.size() * 4, 1),
      _allBins(allBins),
      _netlistHG(hgraphOfNetlist) {
  static_cast<void>(verb);

  _numTerminals = _netlistHG.getNumTerminals();
  //  unsigned numNonTerms  = getNumNodes() - _numTerminals;

  //    _nodeNames = vector<string>(getNumNodes());
  _edges.reserve(_netlistHG.getNumEdges() + getNumNodes() * 4);

  // all existing nodes are the same...add them in first.

  itHGFNodeGlobal n;
  for (n = _netlistHG.nodesBegin(); n != _netlistHG.nodesEnd(); n++) {
    unsigned nId = (*n)->getIndex();

    // note: names are in fact not very essential in this context..
    // though they will be helpful in debugging. This process
    // can be speedup by not having node names in the new HG
    //	string newName = (*n)->getName();
    //	_nodeNamesMap[newName] = nId;
    //	_nodeNames[nId] = newName;

    // in the interest of speed..this could prob. be done w/
    // a bin memory copy..but this would be much more fragile
    setWeight(nId, getWeight((*n)->getIndex()));
  }

  // setup (name and give 0 weight) to the 4*numBins new nodes

  unsigned nId;
  for (nId = _netlistHG.getNumNodes(); nId < _nodes.size(); nId++) {
    //	char* newName = new char[20];

    //	unsigned binNum = (nId - _netlistHG.getNumNodes())/4;
    //	unsigned subNum   = (nId - _netlistHG.getNumNodes()) % 4;

    //	sprintf(newName, "BinNode%d-%d", binNum, subNum);

    //	_nodeNamesMap[newName] = nId;
    setWeight(nId, 0);
    //	_nodeNames[nId] = newName;
  }

  // copy all existing edges

  itHGFEdgeGlobal e;
  for (e = _netlistHG.edgesBegin(); e != _netlistHG.edgesEnd(); e++) {
    const HGFEdge& origE = **e;
    HGFEdge& newEdge = *(addEdge(origE.getWeight()));

    itHGFNodeLocal nIt;
    for (nIt = origE.nodesBegin(); nIt != origE.nodesEnd(); nIt++)
      addSrcSnk(**nIt, newEdge);
  }

  // add in the new edges

  unsigned numOrigNodes = _netlistHG.getNumNodes();
  for (unsigned bin = 0; bin < _allBins.size(); bin++) {
    const CapoBin& Bin = *_allBins[bin];

    const HGFNode& subNode0 = *_nodes[numOrigNodes + bin * 4];
    const HGFNode& subNode1 = *_nodes[numOrigNodes + bin * 4 + 1];
    const HGFNode& subNode2 = *_nodes[numOrigNodes + bin * 4 + 2];
    const HGFNode& subNode3 = *_nodes[numOrigNodes + bin * 4 + 3];

    vector<unsigned>::const_iterator cellIds;
    for (cellIds = Bin.cellIdsBegin(); cellIds != Bin.cellIdsEnd(); cellIds++) {
      const HGFNode& node = _netlistHG.getNodeByIdx(*cellIds);
      unsigned origDegree = node.getDegree();
      unsigned numToAdd =
          static_cast<unsigned>(ceil((float)(origDegree) * 0.25));

      abkfatal3(origDegree == 0 || numToAdd >= 1, origDegree,
                " numToAdd is 0 for a node with degree > 0 ", numToAdd);

      unsigned k;
      for (k = 0; k < numToAdd; k++) {
        HGFEdge* newE = addEdge(1);
        addSrcSnk(subNode0, *newE);
        addSrcSnk(node, *newE);
      }
      for (k = 0; k < numToAdd; k++) {
        HGFEdge* newE = addEdge(1);
        addSrcSnk(subNode1, *newE);
        addSrcSnk(node, *newE);
      }
      for (k = 0; k < numToAdd; k++) {
        HGFEdge* newE = addEdge(1);
        addSrcSnk(subNode2, *newE);
        addSrcSnk(node, *newE);
      }
      for (k = 0; k < numToAdd; k++) {
        HGFEdge* newE = addEdge(1);
        addSrcSnk(subNode3, *newE);
        addSrcSnk(node, *newE);
      }
    }
  }

  finalize();
}
