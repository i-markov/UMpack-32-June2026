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

// Created 051103 by David Papa to move all auxiliary class code and
//       definitions out of RBPlacement.h

#ifndef _RBPLACE_AUX_H_
#define _RBPLACE_AUX_H_

#include <vector>

class CompareCellIdsByX {
  const Placement& _placement;

 public:
  CompareCellIdsByX(const Placement& place) : _placement(place) {}

  bool operator()(unsigned id1, unsigned id2) {
    return _placement[id1].x < _placement[id2].x;
  }
};

// added by sadya to generate LEF/DEF files from RBPlacement format
// used to sort the master nodes
struct sort_fn
{
  bool operator()(const std::vector<double>& node1,
                  const std::vector<double>& node2) {
    unsigned size1 = node1.size();
    unsigned size2 = node2.size();
    unsigned minSize = (size1 < size2) ? size1 : size2;
    if (node1[0] * node1[1] != node2[0] * node2[1])
      return (node1[0] * node1[1] < node2[0] * node2[1]);

    for (unsigned i = 0; i < minSize; ++i) {
      if (node1[i] != node2[i]) return (node1[i] < node2[i]);
    }
    return (node1[0] < node2[0]);
  }
};

struct masterNodeInfo {
  std::vector<double> nodeInfo;
  char name[20];
  unsigned index;
  std::vector<char*> pinNames;
  std::vector<Point> pinOffsets;
};

class nodesMasInfo  // this class maintains the mapping between node and master
                    // node
{
 public:
  std::vector<unsigned> masterNodeIndex;

  nodesMasInfo(unsigned numCells) { masterNodeIndex.resize(numCells); }

  void putIndex(unsigned masIndex, unsigned index) {
    masterNodeIndex[index] = masIndex;
  }

  unsigned getMasterIndex(unsigned index) { return masterNodeIndex[index]; }
};

#endif
