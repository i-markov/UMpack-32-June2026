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

// Created: July 25, 1997 by Igor Markov
// ... single-user development until August 13 not documented
// 970813 ilm added Constraints::getPermutation()
// 970813 mro added typedefs like rPlacement to make it compile on MSVC++
// 970815 ilm added Constraints::allSatisfied()
// 970825 ilm added Placement re-association ctor to class Constraints
// 970825 ilm changed Constraints copy-like ctors to compactify if the
//                argument was compactified
// 970828 ilm separated restriction from pull-back
// 970828 ilm reimplemented restriction with partial inverse
//                   (was with sorting and binary_search)
// 971126 ilm added Constraints::getFixedBBox()

#include <ABKCommon/abkassert.h>
#include <ABKCommon/abktypes.h>
#include <Constraints/constraint.h>
#include <Constraints/constraints.h>
#include <Placement/placement.h>
#include <Placement/subPlacement.h>

#include <algorithm>
#include <iostream>
#include <vector>

using std::cerr;
using std::endl;
using std::ostream;
using std::vector;

typedef const Placement& rPlacement;
typedef Constraint* pConstraint;
typedef FixedTypeConstraint* pFixedTypeConstraint;
typedef const Mapping& rMapping;
class NumberedPt;
typedef const NumberedPt& rNumberedPt;

Constraints::Constraints(const Constraints& from)
    : _compactified(false),
      _meat(from.getSize()),
      _fgcp(0),
      _fxgcp(0),
      _fygcp(0) {
  for (unsigned i = 0; i < from.getSize(); i++) _meat[i] = from[i].clone();
  if (from.isCompactified()) compactify();
}

Constraints& Constraints::operator=(const Constraints& c) {
  _compactified = false;
  _meat.resize(c.getSize());
  _fgcp = 0;
  _fxgcp = 0;
  _fygcp = 0;

  for (unsigned i = 0; i < c.getSize(); i++) _meat[i] = c[i].clone();
  if (c.isCompactified()) compactify();

  return *this;
}

/* (old) restriction done with sorting and binary_search
Constraints::Constraints(const Constraints& orig,const Mapping& restrictTo)
: _meat(0), _compactified(false), _fgcp(0),_fxgcp(0),_fygcp(0)
// restriction:  can result in empty constraint for non-empty original
{
  vector<unsigned> argIndices(restrictTo);
  sort(argIndices.begin(),argIndices.end());

  for(unsigned k=0;k<orig._meat.size();k++)
  {
     Constraint* tmpC=orig._meat[k]->_restriction(argIndices);
     if (tmpC)   _meat.push_back(tmpC);
  }
  if (orig.isCompactified())  compactify();
}
*/

// (new) restriction done with partial inverse
Constraints::Constraints(const Constraints& orig, const Mapping& restrictTo)
    : _compactified(false),
      _meat(0),
      _fgcp(0),
      _fxgcp(0),
      _fygcp(0)
// restriction:  can result in empty constraint for non-empty original
{
  std::vector<unsigned> partialInverse;
  restrictTo.getPartialInverse(partialInverse);

  for (unsigned k = 0; k < orig._meat.size(); k++) {
    Constraint* tmpC = orig._meat[k]->_restriction(partialInverse);
    if (tmpC) _meat.push_back(tmpC);
  }
  if (orig.isCompactified()) compactify();
}

Constraints::Constraints(const Constraints& orig,
                         const Mapping& pullBackMapping,
                         const Placement& pullBackPl)
    : _compactified(false),
      _meat(0),
      _fgcp(0),
      _fxgcp(0),
      _fygcp(0)
// Pull-back:  can result in empty constraint for non-empty original
{
  abkassert(pullBackPl.getSize() <= orig.getPlacement().getSize(),
            "Can only pullBack constraints to a smaller placement\n");
  abkfatal(pullBackPl.getSize() == pullBackMapping.getSourceSize() &&
               orig.getPlacement().getSize() == pullBackMapping.getDestSize(),
           "Mapping/Placement size mismatch during constraints pullBack\n");

  std::vector<unsigned> partialInverse;
  pullBackMapping.getPartialInverse(partialInverse);

  for (unsigned k = 0; k < orig._meat.size(); k++) {
    Constraint* tmpC = orig._meat[k]->_pullBack(partialInverse, pullBackPl);
    if (tmpC) _meat.push_back(tmpC);
  }
  if (orig.isCompactified()) compactify();
}

Constraints::Constraints(const Constraints& orig, const Placement& pl,
                         Filter flt)
    : _compactified(false),
      _meat(0),
      _fgcp(0),
      _fxgcp(0),
      _fygcp(0)
// Placement re-association (with filtering)
// can end up with an empty object for non-empty argument
{
  if (orig.isEmpty()) {
    _compactified = true;
    return;
  }
  abkfatal(orig.getPlacement().getSize() <= pl.getSize(),
           "Can not re-associate constraints to Placement of a smaller size\n");
  unsigned i;
  switch (flt) {
    case ALL:
      _meat = std::vector<Constraint*>(orig.getSize());
      for (i = 0; i < orig.getSize(); i++) _meat[i] = orig[i].clone(pl);
      break;
    case FIXEDTYPE:
      for (i = 0; i < orig.getSize(); i++)
        if (orig[i].isFixedType()) _meat.push_back(orig[i].clone(pl));
      break;
    case FIXEDONLY:
      for (i = 0; i < orig.getSize(); i++)
        if (orig[i].isFixed()) _meat.push_back(orig[i].clone(pl));
      break;
    default:
      abkfatal(0, " Unknown Constraint::Filter ");
  }
  if (orig.isCompactified()) compactify();
}

Constraints::~Constraints() {
  for (unsigned i = 0; i < getSize(); i++) delete _meat[i];
}

void Constraints::enforceAll() const {
  for (unsigned i = 0; i < _meat.size(); i++) _meat[i]->enforce();
}

bool Constraints::allSatisfied() const {
  for (unsigned k = 0; k < getSize(); k++)
    if (!_meat[k]->isSatisfied()) return false;
  return true;
}

ostream& operator<<(ostream& out, const Constraints& ctrs) {
  if (ctrs.getSize() == 0)
    out << " No spatial constraints " << endl;
  else if (ctrs.getSize() > 1)
    out << " === Total of " << ctrs.getSize()
        << " spatial constraints : " << endl;
  for (unsigned i = 0; i < ctrs.getSize(); i++)
    out << (const char*)ctrs[i].getType() << " of size " << ctrs[i].getSize()
        << endl;
  return out;
}

void Constraints::sortByType() {
  CompareByTypeSize byTypeSize;
  std::sort(_meat.begin(), _meat.end(), byTypeSize);  // maybe stable_sort ?
}

bool Constraints::inBBox(const BBox& bbx) const {
  for (unsigned i = 0; i < _meat.size(); i++)
    if (!_meat[i]->inBBox(bbx)) return false;
  return true;
}

bool Constraints::inBBoxIgnoringSoft(const BBox& bbx) const {
  for (unsigned i = 0; i < _meat.size(); i++) {
    const Constraint& c = *_meat[i];
    if (c.isSoft()) continue;
    if (!c.inBBox(bbx)) return false;
  }
  return true;
}

void Constraints::reorderSubPlWith(const Permutation& perm) {
  if (getSize())
    abkfatal(perm.getSize() == _meat[0]->_spl.getDestSize(),
             "Can't reorder items in a subplacement: wrong permutation size\n");
  for (unsigned i = 0; i < _meat.size(); i++) _meat[i]->reorder(perm);
}

class NumberedPt {
  friend void Constraints::compactify();

  Point pt;
  unsigned num;

  NumberedPt(Point p, unsigned n) : pt(p), num(n) {}

 public:
  NumberedPt() : pt(0, 0), num(0) {}
  bool operator<(const NumberedPt& npt) const {
    if (num < npt.num) return true;
    if (num > npt.num) return false;
    if (pt < npt.pt) return true;
    return false;
  }
};

void Constraints::compactify() {
  if (getSize() == 0) {
    _fgcp = NULL;
    _fxgcp = NULL;
    _fygcp = NULL;
    _compactified = true;
    return;
  }

  sortByType();
  std::vector<Constraint*>::iterator iter, beginFixedType = _meat.begin();
  std::vector<NumberedPt> fixed, fixedX, fixedY;
  unsigned i;

  // fish out FixedType  constraint and copy their contents
  // into the three vectors defined above

  for (iter = _meat.end() - 1; iter >= _meat.begin(); iter--) {
    Constraint::Type type((*iter)->getType());

    std::vector<NumberedPt>* buf = NULL;
    if (type == Constraint::FixedGroup)
      buf = &fixed;
    else if (type == Constraint::FixedXGroup)
      buf = &fixedX;
    else if (type == Constraint::FixedYGroup)
      buf = &fixedY;
    else {
      beginFixedType = iter + 1;
      break;
    }

    // now we know it's a FixedTypeConstraint, hence the casting
    pConstraint tmpc(*iter);  // Constraint*
    pFixedTypeConstraint constr(static_cast<FixedTypeConstraint*>(tmpc));
    // ptr
    rMapping mp(constr->getMapping());      // const ref
    rPlacement pl(constr->getFixedLocs());  // const ref

    for (i = 0; i < constr->getSize(); i++)
      buf->push_back(NumberedPt(pl[i], mp[i]));
  }

  // Sort, remove duplicate constraints and catch conflicts

  if (fixed.size() > 1) {
    std::sort(fixed.begin(), fixed.end());
    NumberedPt lastNuPt(fixed[0]);

    for (i = 1; i < fixed.size();
         i++)  // fixed.size() can change in the loop !!!
    {
      rNumberedPt curNuPt(fixed[i]);
      if (lastNuPt.num == curNuPt.num) {
        if (mDist(lastNuPt.pt, curNuPt.pt) > 1e-4) {
          cerr << " Conflicting FixedGroup constraints for cell " << curNuPt.num
               << " : " << endl
               << lastNuPt.pt << "        and " << curNuPt.pt << endl;
          abkfatal(0, "");
        }
        fixed.erase(fixed.begin() + i, fixed.begin() + i + 1);
        i--;
      } else
        lastNuPt = curNuPt;
    }
  }

  if (fixedX.size() > 1) {
    std::sort(fixedX.begin(), fixedX.end());
    NumberedPt lastNuPt(fixedX[0]);
    for (i = 1; i < fixedX.size();
         i++)  // fixedX.size() can change in the loop !!!
    {
      rNumberedPt curNuPt(fixedX[i]);
      if (lastNuPt.num == curNuPt.num) {
        if (fabs(lastNuPt.pt.x - curNuPt.pt.x) > 1e-4) {
          cerr << " Conflicting FixedXGroup constraints for cell "
               << curNuPt.num << " : " << lastNuPt.pt.x << "  and   "
               << curNuPt.pt.x << endl;
          abkfatal(0, "");
        }
        fixedX.erase(fixedX.begin() + i, fixedX.begin() + i + 1);
        i--;
      } else
        lastNuPt = curNuPt;
    }
  }

  if (fixedY.size() > 1) {
    std::sort(fixedY.begin(), fixedY.end());
    NumberedPt lastNuPt(fixedY[0]);
    for (i = 1; i < fixedY.size();
         i++)  // fixedY.size() can change in the loop !!!
    {
      rNumberedPt curNuPt(fixedY[i]);
      if (lastNuPt.num == curNuPt.num) {
        if (fabs(lastNuPt.pt.y - curNuPt.pt.y) > 1e-4) {
          cerr << " Conflicting FixedYGroup constraints for cell "
               << curNuPt.num << " : " << lastNuPt.pt.y << "   and    "
               << curNuPt.pt.y << endl;
          abkfatal(0, "");
        }
        fixedY.erase(fixedY.begin() + i, fixedY.begin() + i + 1);
        i--;
      } else
        lastNuPt = curNuPt;
    }
  }
  rPlacement realPl(getPlacement());  // ref
  unsigned rPlS = realPl.getSize();

  std::vector<Constraint*>::iterator citer;
  for (citer = beginFixedType; citer < _meat.end(); citer++) delete *citer;

  _meat.erase(beginFixedType, _meat.end());

  if (fixedY.size() > 0) {
    Mapping fixedYMap(fixedY.size(), rPlS);
    Placement fixedYPl(fixedY.size());
    for (i = 0; i < fixedY.size(); i++) {
      fixedYMap[i] = fixedY[i].num;
      fixedYPl[i] = fixedY[i].pt;
    }
    SubPlacement fixedYSub(realPl, fixedYMap);
    FixedYGroupConstraint fgyc(fixedYSub, fixedYPl);
    add(fgyc);
  }

  if (fixedX.size() > 0) {
    Mapping fixedXMap(fixedX.size(), rPlS);
    Placement fixedXPl(fixedX.size());
    for (i = 0; i < fixedX.size(); i++) {
      fixedXMap[i] = fixedX[i].num;
      fixedXPl[i] = fixedX[i].pt;
    }
    SubPlacement fixedXSub(realPl, fixedXMap);
    FixedXGroupConstraint fgxc(fixedXSub, fixedXPl);
    add(fgxc);
  }

  if (fixed.size() > 0) {
    Mapping fixedMap(fixed.size(), rPlS);
    Placement fixedPl(fixed.size());
    for (i = 0; i < fixed.size(); i++) {
      fixedMap[i] = fixed[i].num;
      fixedPl[i] = fixed[i].pt;
    }
    //  cerr << " XY " << realPl.getSize() << " " << fixedMap.getDestSize() <<
    //  endl;
    SubPlacement fixedSub(realPl, fixedMap);
    FixedGroupConstraint fgc(fixedSub, fixedPl);
    add(fgc);
  }

  pConstraint tmpc;
  unsigned idx = _meat.size();
  if (fixed.size() > 0) {
    tmpc = _meat[--idx];
    _fgcp = static_cast<FixedGroupConstraint*>(tmpc);
  } else
    _fgcp = NULL;

  if (fixedX.size() > 0) {
    tmpc = _meat[--idx];
    _fxgcp = static_cast<FixedXGroupConstraint*>(tmpc);
  } else
    _fxgcp = NULL;

  if (fixedY.size() > 0) {
    tmpc = _meat[--idx];
    _fygcp = static_cast<FixedYGroupConstraint*>(tmpc);
  } else
    _fygcp = NULL;

  _compactified = true;
}

Permutation& Constraints::getPermutation(Permutation& newIndices) const
// produces a permutation which sorts fixed items first
// then fixedX, then fixedY --- newIndices[i] will give the new
// index of the cell i. To reorder a vector use
// template reorderVector() from Placement/mapping.h
{
  abkfatal(isCompactified(),
           " Uncompactified constraint can't produce a permutation \n");

  unsigned plSize = getPlacement().getSize();
  Permutation perm(plSize);  // later will do newIndices=perm;

  unsigned numFixed, numFixedX, numFixedY;
  unsigned i;

  std::vector<unsigned> partInv4Fixed(plSize, unsigned(-1));
  std::vector<unsigned> partInv4FixedX(plSize, unsigned(-1));
  std::vector<unsigned> partInv4FixedY(plSize, unsigned(-1));

  if (hasFGC()) {
    getFGC().getMapping().getPartialInverse(partInv4Fixed);
    numFixed = getFGC().getSize();
  } else
    numFixed = 0;

  if (hasFXGC()) {
    getFXGC().getMapping().getPartialInverse(partInv4FixedX);
    numFixedX = getFXGC().getSize();
  } else
    numFixedX = 0;

  if (hasFYGC()) {
    getFYGC().getMapping().getPartialInverse(partInv4FixedY);
    numFixedY = getFYGC().getSize();
  } else
    numFixedY = 0;

  unsigned fixedCo = 0, fixedXCo = numFixed, fixedYCo = numFixed + numFixedX;
  unsigned movableCo = numFixed + numFixedX + numFixedY;

  for (i = 0; i < plSize; i++) {
    abkfatal(
        partInv4Fixed[i] == unsigned(-1) || partInv4FixedX[i] == unsigned(-1),
        "FixedGroup and FixedXGroup constraints are mutually exclusive for a "
        "cell\n");

    abkfatal(
        partInv4Fixed[i] == unsigned(-1) || partInv4FixedY[i] == unsigned(-1),
        "FixedGroup and FixedYGroup constraints are mutually exclusive for a "
        "cell\n");

    abkfatal(
        partInv4FixedX[i] == unsigned(-1) || partInv4FixedY[i] == unsigned(-1),
        "FixedXGroup and FixedYGroup constraints are mutually exclusive for a "
        "cell\n");

    if (partInv4Fixed[i] != unsigned(-1))
      perm[i] = fixedCo++;
    else if (partInv4FixedX[i] != unsigned(-1))
      perm[i] = fixedXCo++;
    else if (partInv4FixedY[i] != unsigned(-1))
      perm[i] = fixedYCo++;
    else
      perm[i] = movableCo++;
  }
  abkfatal(fixedCo == numFixed && fixedXCo == numFixed + numFixedX &&
               fixedYCo == numFixed + numFixedX + numFixedY &&
               movableCo == plSize,
           "Internal error in Constraints::getPermutation()");

  newIndices = perm;
  return newIndices;
}

bit_vector& Constraints::getFGCbitVec(bit_vector& isFixed) const {
  abkassert(_compactified,
            "Uncompactified constraints can't produce FGCbitVec\n");
  // by Xiaojian to fix unplaced I/O problem
  // abkfatal(!isEmpty(),"Empty constraints can't produce FGCbitVec\n");
  if (isEmpty()) {
    abkwarn(0, "Empty constraints --- I/O pads have no locations!\n");
    return isFixed;
  }

  if (!hasFGC()) {
    isFixed = bit_vector(getPlacement().getSize(), false);
    return isFixed;
  }
  // "else"
  std::vector<unsigned> tmpVec;
  getFGC().getMapping().getPartialInverse(tmpVec);
  isFixed = bit_vector(getPlacement().getSize());
  for (unsigned k = 0; k < tmpVec.size(); k++)
    isFixed[k] = (tmpVec[k] != unsigned(-1));
  return isFixed;
}

bit_vector& Constraints::getFXGCbitVec(bit_vector& isXFixed) const {
  abkassert(_compactified,
            "Uncompactified constraints can't produce FXGCbitVec\n");
  abkfatal(!isEmpty(), "Empty constraints can't produce FXGCbitVec\n");
  if (!hasFXGC()) {
    isXFixed = bit_vector(getPlacement().getSize(), false);
    return isXFixed;
  }
  // "else"
  std::vector<unsigned> tmpVec;
  getFXGC().getMapping().getPartialInverse(tmpVec);
  isXFixed = bit_vector(getPlacement().getSize());
  for (unsigned k = 0; k < tmpVec.size(); k++)
    isXFixed[k] = (tmpVec[k] != unsigned(-1));
  return isXFixed;
}

bit_vector& Constraints::getFYGCbitVec(bit_vector& isYFixed) const {
  abkassert(_compactified,
            "Uncompactified constraints can't produce FYGCbitVec\n");
  abkfatal(!isEmpty(), "Empty constraints can't produce FYGCbitVec\n");
  if (!hasFYGC()) {
    isYFixed = bit_vector(getPlacement().getSize(), false);
    return isYFixed;
  }
  // "else"
  std::vector<unsigned> tmpVec;
  getFYGC().getMapping().getPartialInverse(tmpVec);
  isYFixed = bit_vector(getPlacement().getSize());
  for (unsigned k = 0; k < tmpVec.size(); k++)
    isYFixed[k] = (tmpVec[k] != unsigned(-1));
  return isYFixed;
}

BBox Constraints::getFixedBBox() const {
  abkfatal(_compactified, "Constraints uncompactified");
  BBox bbox;

  if (hasFGC()) bbox = getFGC().getBBox();

  if (hasFXGC()) {
    Interval xinterv = getFXGC().getFixedXInterval();
    if (bbox.xMin > xinterv.low) bbox.xMin = xinterv.low;
    if (bbox.xMax < xinterv.high) bbox.xMax = xinterv.high;
  }

  if (hasFXGC()) {
    Interval yinterv = getFXGC().getFixedXInterval();
    if (bbox.yMin > yinterv.low) bbox.yMin = yinterv.low;
    if (bbox.yMax < yinterv.high) bbox.yMax = yinterv.high;
  }

  return bbox;
}
