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

#include "primMST.h"

#include <ABKCommon/abkassert.h>

using std::vector;

PrimMST::PrimMST(const vector<Point>& points) : MSTBase(points) {
  int size = points.size();
  // vector<double> distToTree; //size, 0.0);
  double distToTree[1000];  // assuming the size is < 1000 pts
  abkfatal(size <= 1000, "size too large\n");
  int root = 0;
  int i, j, k;
  double d1, max_d, best_d;
  int best_i;

  max_d = 0;
  for (i = 0; i < size; i++) {
    _mstNodes[i].parent = root;
    distToTree[i] = mDist(_mstNodes[i].edge.first, _mstNodes[root].edge.first);
    if (distToTree[i] > max_d) max_d = distToTree[i];
  }
  _mstNodes[root].parent = root;
  distToTree[root] = 0;

  for (k = 1; k < size; k++)  // n-1 points to add
  {
    best_d = max_d + 1;
    best_i = -1;
    for (i = 0; i < size; i++) {
      if ((distToTree[i] != 0) && (distToTree[i] < best_d)) {
        best_d = distToTree[i];
        best_i = i;
      }
    }
    // pt[best_i] enters the tree
    distToTree[best_i] = 0;

    // update distToTree for remaining nodes
    for (j = 0; j < size; j++) {
      if (distToTree[j] != 0) {
        d1 = mDist(_mstNodes[j].edge.first, _mstNodes[best_i].edge.first);
        if (d1 < distToTree[j]) {
          distToTree[j] = d1;
          _mstNodes[j].parent = best_i;
        }
      }
    }
  }

  // correct the parent's coordinate and compute cost
  _mstCost = 0;
  for (i = 0; i < size; i++) {
    _mstNodes[i].edge.second = _mstNodes[_mstNodes[i].parent].edge.first;
    _mstCost += _mstNodes[i].edge.getMDist();
  }
  /*
  print(cout);
  for(i = 0; i < size; i++)
  {
      cout << "Node: " << i << endl
           << _mstNodes[i].edge.first << endl
           << _mstNodes[i].edge.second << endl << flush;
  }
  cout << "End primMST\n" << flush;
  cout << "Cost = " << _mstCost << endl;
  */
}
