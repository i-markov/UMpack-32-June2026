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

#ifndef _SLICING_TREE_H_
#define _SLICING_TREE_H_

#include <deque>
#include <map>

#include "Geoms/bbox.h"

class oaSlicingTree;

namespace RBPlace {
typedef unsigned slicingTreeIdxT;
typedef double slicingTreeCoordT;

// This is a binary tree, i.e. not a normalized slicing tree
class SlicingTree {
  friend class oaSlicingTree;

 public:
  SlicingTree(const BBox&, const std::string& name = std::string());
  SlicingTree(slicingTreeCoordT left, slicingTreeCoordT bottom,
              slicingTreeCoordT right, slicingTreeCoordT top,
              const std::string& name = std::string());
  SlicingTree(const SlicingTree& orig);
  // SlicingTree& operator=()  // you want it you write it.

  // Use this sparingly, its a linear algorithm == slow
  slicingTreeIdxT getIdxFromCoords(slicingTreeCoordT left,
                                   slicingTreeCoordT bottom,
                                   slicingTreeCoordT right,
                                   slicingTreeCoordT top);
  // similar for a BBox TODO

  BBox getBBoxFromIdx(slicingTreeIdxT);
  slicingTreeIdxT getIdxFromName(std::string);
  slicingTreeCoordT getCutLine(slicingTreeIdxT);
  std::pair<slicingTreeIdxT, slicingTreeIdxT> getChildren(slicingTreeIdxT);
  slicingTreeIdxT getParent(slicingTreeIdxT);
  void cutH(slicingTreeIdxT, slicingTreeCoordT);
  void cutV(slicingTreeIdxT, slicingTreeCoordT);
  void eraseDecentants(slicingTreeIdxT);
  void erase(slicingTreeIdxT);

  bool isRoot(slicingTreeIdxT);
  bool isLeaf(slicingTreeIdxT);
  static slicingTreeCoordT getCutBetweenSiblings(BBox b1, BBox b2);
  // getHeight(void);  //You want it, impelement it once :) //height of the tree
  // of coz size(); //you want it, implement it once :)
 private:
  std::map<std::string, slicingTreeIdxT> _nameToIdx;
  std::deque<BBox> _idxToRect;
  std::vector<slicingTreeIdxT> _parent;
  std::vector<slicingTreeIdxT> _leftOrBottom;
  std::vector<slicingTreeIdxT> _rightOrTop;
  std::vector<bool> _leftIsNull;
  std::vector<bool> _rightIsNull;
};
}  // namespace RBPlace

#endif
