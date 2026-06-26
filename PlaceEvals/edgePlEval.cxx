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

// created by Igor Markov on 02/25/2000
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "edgePlEval.h"

#include <iostream>
#include <limits>
#include "HGraph/hgFixed.h"
#include "Placement/placement.h"

using std::endl;
using std::ostream;

EdgePlEval::EdgePlEval(const HGraphFixed& hg, const Placement& loc,
                       unsigned maxNetSize) {
  double maxCost = 0.0;
  double minCost = std::numeric_limits<double>::max();
  double totalCost = 0.0;
  unsigned numEdges = 0;
  for (itHGFEdgeGlobal edgeIt = hg.edgesBegin(); edgeIt != hg.edgesEnd();
       edgeIt++) {
    HGFEdge& edge = (**edgeIt);
    if (edge.getDegree() > maxNetSize) continue;
    //    cout << " Edge idx: " << edge.getIndex() << endl;
    double edgeCost;
#ifdef SIGNAL_DIRECTIONS
    // first go over edges whose sources are pure sources
    for (itHGFNodeLocal srcIt = edge.srcsBegin(); srcIt != edge.srcsEnd();
         srcIt++) {
      unsigned srcIdx = (*srcIt)->getIndex();
      //        cout << "   source idx: " << srcIdx << endl;
      itHGFNodeLocal snkIt;
      for (snkIt = edge.srcSnksBegin(); snkIt != edge.srcSnksEnd(); snkIt++) {
        unsigned snkIdx = (*snkIt)->getIndex();
        //           cout << "      src-sink idx: " << snkIdx << endl;
        edgeCost = mDist(loc[srcIdx], loc[snkIdx]);
        totalCost += edgeCost;
        if (edgeCost > maxCost) maxCost = edgeCost;
        if (edgeCost < minCost) minCost = edgeCost;
        numEdges++;
      }
      for (snkIt = edge.snksBegin(); snkIt != edge.snksEnd(); snkIt++) {
        unsigned snkIdx = (*snkIt)->getIndex();
        //           cout << "      sink idx: " << snkIdx << endl;
        edgeCost = mDist(loc[srcIdx], loc[snkIdx]);
        totalCost += edgeCost;
        if (edgeCost > maxCost) maxCost = edgeCost;
        if (edgeCost < minCost) minCost = edgeCost;
        numEdges++;
      }
    }
#endif

    // then go over edges whose sources are sink-sources
#ifdef SIGNAL_DIRECTIONS
    for (itHGFNodeLocal srcSnkIt = edge.srcSnksBegin();
         srcSnkIt != edge.srcSnksEnd(); ++srcSnkIt)
#else
    for (itHGFNodeLocal srcSnkIt = edge.nodesBegin();
         srcSnkIt != edge.nodesEnd(); ++srcSnkIt)
#endif
    {
      unsigned srcIdx = (*srcSnkIt)->getIndex();
      //        cout << "   source-sink idx: " << srcIdx << endl;
      itHGFNodeLocal snkIt;
#ifdef SIGNAL_DIRECTIONS
      for (snkIt = edge.srcSnksBegin(); snkIt != edge.srcSnksEnd(); ++snkIt)
#else
      for (snkIt = edge.nodesBegin(); snkIt != edge.nodesEnd(); ++snkIt)
#endif
      {
        unsigned snkIdx = (*snkIt)->getIndex();
        if (snkIdx == srcIdx) continue;
        //           cout << "      src-sink idx: " << snkIdx << endl;
        edgeCost = mDist(loc[srcIdx], loc[snkIdx]);
        totalCost += edgeCost;
        if (edgeCost > maxCost) maxCost = edgeCost;
        if (edgeCost < minCost) minCost = edgeCost;
        numEdges++;
      }
#ifdef SIGNAL_DIRECTIONS
      for (snkIt = edge.snksBegin(); snkIt != edge.snksEnd(); snkIt++) {
        unsigned snkIdx = (*snkIt)->getIndex();
        //           cout << "      sink idx: " << snkIdx << endl;
        edgeCost = mDist(loc[srcIdx], loc[snkIdx]);
        totalCost += edgeCost;
        if (edgeCost > maxCost) maxCost = edgeCost;
        if (edgeCost < minCost) minCost = edgeCost;
        numEdges++;
      }
#endif
    }
  }
  _maxEdgeLen = maxCost;
  _minEdgeLen = minCost;
  _avgEdgeLen = totalCost / numEdges;
  _numEdges = numEdges;
}

ostream& operator<<(ostream& out, const EdgePlEval& eval) {
  out << " Total edges: " << eval.getNumEdges()
      << "  Edge lengths: " << eval.getAvgEdge() << "(avg) "
      << eval.getMaxEdge() << "(max) " << eval.getMinEdge() << "(min) " << endl;
  return out;
}
