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

#include "slicingTree.h"

namespace RBPlace {
typedef unsigned slicingTreeIdxT;
typedef double slicingTreeCoordT;

SlicingTree::SlicingTree(const BBox& _init, const std::string& name) {
  slicingTreeIdxT new_idx(0);
  _nameToIdx.clear();
  _nameToIdx[name] = new_idx;

  _idxToRect.clear();
  _idxToRect.push_back(_init);

  _parent.clear();
  _parent.push_back(new_idx);  // parent of root is itself

  _leftOrBottom.clear();
  _leftOrBottom.push_back(UINT_MAX);

  _rightOrTop.clear();
  _rightOrTop.push_back(UINT_MAX);

  _leftIsNull.clear();
  _leftIsNull.push_back(true);

  _rightIsNull.clear();
  _rightIsNull.push_back(true);
}

SlicingTree::SlicingTree(slicingTreeCoordT left, slicingTreeCoordT bottom,
                         slicingTreeCoordT right, slicingTreeCoordT top,
                         const std::string& name) {
  slicingTreeIdxT new_idx(0);
  _nameToIdx.clear();
  _nameToIdx[name] = new_idx;

  _idxToRect.clear();
  _idxToRect.push_back(BBox(left, bottom, right, top));

  _parent.clear();
  _parent.push_back(new_idx);  // parent of root is itself

  _leftOrBottom.clear();
  _leftOrBottom.push_back(UINT_MAX);

  _rightOrTop.clear();
  _rightOrTop.push_back(UINT_MAX);

  _leftIsNull.clear();
  _leftIsNull.push_back(true);

  _rightIsNull.clear();
  _rightIsNull.push_back(true);
}

SlicingTree::SlicingTree(const SlicingTree& orig)
    : _nameToIdx(orig._nameToIdx),
      _idxToRect(orig._idxToRect),
      _parent(orig._parent),
      _leftOrBottom(orig._leftOrBottom),
      _rightOrTop(orig._rightOrTop),
      _leftIsNull(orig._leftIsNull),
      _rightIsNull(orig._rightIsNull) {}

// Use this sparingly, its a linear algorithm == slow
slicingTreeIdxT SlicingTree::getIdxFromCoords(slicingTreeCoordT left,
                                              slicingTreeCoordT bottom,
                                              slicingTreeCoordT right,
                                              slicingTreeCoordT top) {
  return slicingTreeIdxT(0);
  // TODO:  this is just a stub, impelement if u want to use it ;)
}
// TODO: similar for a BBox, bridge to call above for convenience

BBox SlicingTree::getBBoxFromIdx(slicingTreeIdxT idx) {
  abkfatal(idx < _idxToRect.size(),
           "Index out of range in slicing tree lookup");
  return _idxToRect[idx];
}

slicingTreeIdxT SlicingTree::getIdxFromName(std::string name) {
  if (_nameToIdx.find(name) != _nameToIdx.end())
    return _nameToIdx[name];
  else
    abkfatal(0, "Get slicing tree rect by name failed");
}

slicingTreeCoordT SlicingTree::getCutLine(slicingTreeIdxT idx) {
  // In the direction where not both coordinates agree
  // return the middle coordinate, it should agree

  abkfatal(idx < _idxToRect.size(),
           "Index out of range in SlicingTree::getCutLine");
  abkfatal(!isLeaf(idx), "Cannot getCutLine of a Leaf node");
  if (_idxToRect[_leftOrBottom[idx]].yMax ==
          _idxToRect[_rightOrTop[idx]].yMax &&
      _idxToRect[_leftOrBottom[idx]].yMin ==
          _idxToRect[_rightOrTop[idx]].yMin) {
    abkfatal(_idxToRect[_leftOrBottom[idx]].xMax ==
                 _idxToRect[_rightOrTop[idx]].xMin,
             "Get cutline internal error");
    return _idxToRect[_leftOrBottom[idx]].xMax;
  } else if (_idxToRect[_leftOrBottom[idx]].xMin ==
                 _idxToRect[_rightOrTop[idx]].xMin &&
             _idxToRect[_leftOrBottom[idx]].xMax ==
                 _idxToRect[_rightOrTop[idx]].xMax) {
    abkfatal(_idxToRect[_leftOrBottom[idx]].yMax ==
                 _idxToRect[_rightOrTop[idx]].yMin,
             "Get cutline internal error");
    return _idxToRect[_leftOrBottom[idx]].yMax;
  } else {
    abkfatal(0, "Get cutline internal error");
  }
  abkfatal(0, "Get cutline internal error, control flow shouldnt reach here");
  return _idxToRect[_leftOrBottom[idx]].xMax;
}

// assumes b1 = left or bottom
//         b2 = right or top
slicingTreeCoordT SlicingTree::getCutBetweenSiblings(BBox b1, BBox b2) {
  // In the direction where not both coordinates agree
  // return the middle coordinate, it should agree
  if (b1.yMax == b2.yMax && b1.yMin == b2.yMin) {
    // these blocks are siblings of a horizontal cut
    abkfatal(b1.xMax == b2.xMin, "Get cutline internal error");
    return b1.xMax;
  } else if (b1.xMin == b2.xMin && b1.xMax == b2.xMax) {
    abkfatal(b2.yMax == b1.yMin, "Get cutline internal error");
    return b2.yMax;
  } else if ((b1.xMin == b2.xMin) ^ (b1.xMax == b2.xMax) &&
             (b1.yMin == b2.yMax) ^ (b1.yMax == b2.yMin)) {
    // exactly one equal in both domains,
    // will take the common coordinate from the
    // dimension with differing coordinates of
    // greatest percentage differance

    // TODO account for this in the slicing tree,
    // here we will just return the cutline, but
    // the slicing tree does not account for this "hole"

    bool diffx_min = b1.xMin != b2.xMin;
    bool diffy_minmax = b1.yMin != b2.yMax;
    double xpct_diff = 0.0, ypct_diff = 0.0;

    if (diffx_min)
      xpct_diff = abs(b1.xMin - b2.xMin) / b1.xMin;
    else
      xpct_diff = abs(b1.xMax - b2.xMax) / b1.xMax;

    if (diffy_minmax)
      ypct_diff = abs(b1.yMin - b2.yMax) / b1.yMin;
    else
      ypct_diff = abs(b1.yMax - b2.yMin) / b1.yMax;

    if (ypct_diff > xpct_diff)
      if (diffx_min)
        return b1.xMax;
      else
        return b1.xMin;
    else if (diffy_minmax)
      return b1.yMax;
    else
      return b1.yMin;

  } else {
    cerr << endl << endl;
    cerr << "B1: " << b1 << endl;
    cerr << "B2: " << b2 << endl;
    return -6;
    abkfatal(0, "getCutBetweenSiblings called on non-sibling \'BBox\'es");
  }
}

std::pair<slicingTreeIdxT, slicingTreeIdxT> SlicingTree::getChildren(
    slicingTreeIdxT idx) {
  abkfatal(_leftOrBottom.size() == _idxToRect.size(),
           "Slicing tree abused and/or broken");
  abkfatal(_rightOrTop.size() == _idxToRect.size(),
           "Slicing tree abused and/or broken");
  abkfatal(idx < _leftOrBottom.size(),
           "Index out of range in slicing tree lookup");
  return std::pair<slicingTreeIdxT, slicingTreeIdxT>(_leftOrBottom[idx],
                                                     _rightOrTop[idx]);
}

slicingTreeIdxT SlicingTree::getParent(slicingTreeIdxT idx) {
  abkfatal(_parent.size() == _idxToRect.size(),
           "Slicing tree abused and/or broken");
  abkfatal(idx < _parent.size(), "Index out of range in slicing tree lookup");
  return _parent[idx];
}

void SlicingTree::cutH(slicingTreeIdxT idx, slicingTreeCoordT cut) {
  abkfatal(idx < _idxToRect.size(),
           "Index out of range in slicing tree lookup");

  BBox boxBottom(_idxToRect[idx]);
  boxBottom.yMax = cut;  // Top of bottom box is cutline
  slicingTreeIdxT idxBottom = _idxToRect.size();
  _idxToRect.push_back(boxBottom);
  _parent.push_back(idx);
  _leftOrBottom.push_back(INT_MAX);
  _rightOrTop.push_back(INT_MAX);
  _leftIsNull.push_back(true);
  _rightIsNull.push_back(true);

  BBox boxTop(_idxToRect[idx]);
  boxTop.yMin = cut;  // bottom of the top box is cutline
  slicingTreeIdxT idxTop = _idxToRect.size();
  _idxToRect.push_back(boxTop);
  _parent.push_back(idx);
  _leftOrBottom.push_back(INT_MAX);
  _rightOrTop.push_back(INT_MAX);
  _leftIsNull.push_back(true);
  _rightIsNull.push_back(true);

  _leftOrBottom[idx] = idxBottom;
  _rightOrTop[idx] = idxTop;
  _leftIsNull[idx] = false;
  _rightIsNull[idx] = false;

  abkfatal(_idxToRect.size() == _parent.size(), "SlicingTree abused");
  abkfatal(_idxToRect.size() == _leftOrBottom.size(), "SlicingTree abused");
  abkfatal(_idxToRect.size() == _rightOrTop.size(), "SlicingTree abused");
  abkfatal(_idxToRect.size() == _leftIsNull.size(), "SlicingTree abused");
  abkfatal(_idxToRect.size() == _rightIsNull.size(), "SlicingTree abused");

  return;
}

void SlicingTree::cutV(slicingTreeIdxT idx, slicingTreeCoordT cut) {
  abkfatal(idx < _idxToRect.size(),
           "Index out of range in slicing tree lookup");

  BBox boxLeft(_idxToRect[idx]);
  boxLeft.xMax = cut;  // Right of left box is cutline
  slicingTreeIdxT idxLeft = _idxToRect.size();
  _idxToRect.push_back(boxLeft);
  _parent.push_back(idx);
  _leftOrBottom.push_back(INT_MAX);
  _rightOrTop.push_back(INT_MAX);
  _leftIsNull.push_back(true);
  _rightIsNull.push_back(true);

  BBox boxRight(_idxToRect[idx]);
  boxRight.xMin = cut;  // Left of the right box is cutline
  slicingTreeIdxT idxRight = _idxToRect.size();
  _idxToRect.push_back(boxRight);
  _parent.push_back(idx);
  _leftOrBottom.push_back(INT_MAX);
  _rightOrTop.push_back(INT_MAX);
  _leftIsNull.push_back(true);
  _rightIsNull.push_back(true);

  _leftOrBottom[idx] = idxLeft;
  _rightOrTop[idx] = idxRight;
  _leftIsNull[idx] = false;
  _rightIsNull[idx] = false;

  abkfatal(_idxToRect.size() == _parent.size(), "SlicingTree abused");
  abkfatal(_idxToRect.size() == _leftOrBottom.size(), "SlicingTree abused");
  abkfatal(_idxToRect.size() == _rightOrTop.size(), "SlicingTree abused");
  abkfatal(_idxToRect.size() == _leftIsNull.size(), "SlicingTree abused");
  abkfatal(_idxToRect.size() == _rightIsNull.size(), "SlicingTree abused");

  return;
}

void SlicingTree::eraseDecentants(slicingTreeIdxT idx) {
  abkfatal(idx < _idxToRect.size(),
           "Index out of range in slicing tree lookup");
  // TODO
}

void SlicingTree::erase(slicingTreeIdxT idx) {
  abkfatal(idx < _idxToRect.size(),
           "Index out of range in slicing tree lookup");
  // TODO
}

bool SlicingTree::isRoot(slicingTreeIdxT idx) {
  abkfatal(idx < _idxToRect.size(),
           "Index out of range in slicing tree lookup");
  return _parent[idx] == idx;
}

bool SlicingTree::isLeaf(slicingTreeIdxT idx) {
  abkfatal(idx < _idxToRect.size(),
           "Index out of range in slicing tree lookup");
  return _leftIsNull[idx] && _rightIsNull[idx];
}

}  // namespace RBPlace
