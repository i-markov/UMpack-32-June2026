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

//! author ="Aug 16, 1997 by Igor Markov"
//! CHANGES =" fxiTypeConstr.h 970825 ilm added cloning with placement
//! re-association to all individual Constraint classes 970923 ilm added
//! simplified constructors for FixedTypeConstraints 971126 ilm added
//! FixedGroupConstraint::getBBox() FixedXGroupConstraint::getXInterval()
//! FixedYGroupConstraint::getYInterval() "

// CHANGES
// 970825 ilm added cloning with placement re-association to all
//               individual Constraint classes
// 970923 ilm added simplified constructors for FixedTypeConstraints
// 971126 ilm added  FixedGroupConstraint::getBBox()
//                   FixedXGroupConstraint::getXInterval()
//                   FixedYGroupConstraint::getYInterval()

#ifndef _FIXEDTYPECONSTRS_H_
#define _FIXEDTYPECONSTRS_H_

#include <ABKCommon/abkassert.h>
#include <Constraints/constraint.h>
#include <Geoms/interval.h>
#include <Placement/subPlacement.h>

#include <algorithm>
#include <cmath>

// descriptions follow class declarations
class FixedGroupConstraint;
// all cells must have x_coord and y_coord fixed at given points
class FixedXGroupConstraint;
// all cells must have x_coord fixed at given values
class FixedYGroupConstraint;
// all cells must have y_coord fixed at given values

/* -------------------------------------------------------------------- */

//: Constraint: all cells must have coord fixed at given values
class FixedGroupConstraint : public FixedTypeConstraint {
 protected:
  Constraint* _virtPullBack(const std::vector<unsigned>& hitIndices,
                            const SubPlacement& newSubPl) const;

 public:
  /* don't use both isSatisfied() and getPenalty() ! */
  FixedGroupConstraint(const SubPlacement& spl, const Placement& fixedLocs)
      :  // fixedLocs for constrained cells
         // spl specifies the group of constrained cells
         // (contains a reference to the placement that we're checking)
        FixedTypeConstraint(spl, fixedLocs, FixedGroup) {}

  FixedGroupConstraint(const Placement& referencePl, const Mapping& embed,
                       const Placement& fixedLocs)
      : FixedTypeConstraint(SubPlacement(referencePl, embed), fixedLocs,
                            FixedGroup) {}

  FixedGroupConstraint(const FixedGroupConstraint& fgc)
      : FixedTypeConstraint(fgc._spl, fgc._pl, FixedGroup) {}

  FixedGroupConstraint& operator=(const FixedGroupConstraint& fgc) {
    _pl = fgc._pl;
    _spl = fgc._spl;
    _type = fgc._type;
    return *this;
  }

  virtual ~FixedGroupConstraint(){};

  Constraint* clone() const { return new FixedGroupConstraint(*this); }
  Constraint* clone(const Placement& pl) const {
    SubPlacement newSubPl(pl, _spl.getSubset());
    return new FixedGroupConstraint(newSubPl, _pl);
  }
  // cloning with placement re-association

  bool isSatisfied() const {
    for (unsigned i = 0; i < _spl.getSize(); i++)
      if (_spl[i] != _pl[i]) return false;
    return true;
  }

  double getPenalty() const {
    double max = -1;
    for (unsigned i = _spl.getSize() - 1; i != unsigned(-1); i--)
      if (mDist(_spl[i], _pl[i]) > max) max = mDist(_spl[i], _pl[i]);
    return max;
  }

  void enforce() { _spl.pushForward(_pl); }

  bool inBBox(const BBox& bb) const {
    for (unsigned i = 0; i < _spl.getSize(); i++)
      if (!bb.contains(_pl[i])) return false;
    return true;
  }

  BBox getBBox() const { return _pl.getBBox(); }  // of fixed locations
};

//: Constraints: all cells must have x_coord fixed at given values
class FixedXGroupConstraint : public FixedTypeConstraint {
 protected:
  Constraint* _virtPullBack(const std::vector<unsigned>& hitIndices,
                            const SubPlacement& newSubPl) const;

 public:
  /* don't use both isSatisfied() and getPenalty() ! */
  FixedXGroupConstraint(const SubPlacement& spl, const Placement& fixedLocs)
      :  // fixedLocs for constrained cells (we ignore the y locations)
         // spl specifies the group of constrained cells
         //(contains a reference to the placement that we're checking)
        FixedTypeConstraint(spl, fixedLocs, FixedXGroup) {}

  FixedXGroupConstraint(const Placement& referencePl, const Mapping& embed,
                        const Placement& fixedLocs)
      : FixedTypeConstraint(SubPlacement(referencePl, embed), fixedLocs,
                            FixedXGroup) {}

  FixedXGroupConstraint(const FixedXGroupConstraint& fgc)
      : FixedTypeConstraint(fgc._spl, fgc._pl, FixedXGroup) {}

  FixedXGroupConstraint& operator=(const FixedXGroupConstraint& fgc) {
    _pl = fgc._pl;
    _spl = fgc._spl;
    _type = fgc._type;
    return *this;
  }

  virtual ~FixedXGroupConstraint() {}

  Constraint* clone() const { return new FixedXGroupConstraint(*this); }
  Constraint* clone(const Placement& pl) const {
    SubPlacement newSubPl(pl, _spl.getSubset());
    return new FixedXGroupConstraint(newSubPl, _pl);
  }
  // cloning with placement re-association

  bool isSatisfied() const {
    for (unsigned i = 0; i < _spl.getSize(); i++)
      if (_spl[i].x != _pl[i].x) return false;
    return true;
  }

  double getPenalty() const {
    double max = -1;
    for (unsigned i = 0; i < _spl.getSize(); i++)
      max = std::max(max, fabs(_spl[i].x - _pl[i].x));
    return max;
  }

  void enforce() { _spl.pushForwardX(_pl); }

  bool inBBox(const BBox& bb) const {
    for (unsigned i = 0; i < _spl.getSize(); i++)
      if (bb.xMin > _pl[i].x || bb.xMax < _pl[i].x) return false;
    return true;
  }

  Interval getFixedXInterval() const {
    Interval interv(_pl[0].x, _pl[0].x);
    for (unsigned i = 1; i < _spl.getSize(); i++)
      if (interv.low > _pl[i].x)
        interv.low = _pl[i].x;
      else if (interv.high < _pl[i].x)
        interv.high = _pl[i].x;
    return interv;
  }
};

//: Constraint: all cells must have y_coord fixed at given values
class FixedYGroupConstraint : public FixedTypeConstraint {
 protected:
  Constraint* _virtPullBack(const std::vector<unsigned>& hitIndices,
                            const SubPlacement& newSubPl) const;

 public:
  /* don't use both isSatisfied() and getPenalty() !*/
  FixedYGroupConstraint(const SubPlacement& spl, const Placement& fixedLocs)
      :  // fixedLocs for constrained cells  (we ignore the x locations)
         // spl specifies the group of constrained cells
         // (contains a reference to the placement that we're checking)
        FixedTypeConstraint(spl, fixedLocs, FixedYGroup) {}

  FixedYGroupConstraint(const Placement& referencePl, const Mapping& embed,
                        const Placement& fixedLocs)
      : FixedTypeConstraint(SubPlacement(referencePl, embed), fixedLocs,
                            FixedYGroup) {}

  FixedYGroupConstraint(const FixedYGroupConstraint& fgc)
      : FixedTypeConstraint(fgc._spl, fgc._pl, FixedYGroup) {}

  FixedYGroupConstraint& operator=(const FixedYGroupConstraint& fgc) {
    _pl = fgc._pl;
    _spl = fgc._spl;
    _type = fgc._type;
    return *this;
  }

  virtual ~FixedYGroupConstraint() {}

  Constraint* clone() const { return new FixedYGroupConstraint(*this); }
  Constraint* clone(const Placement& pl) const {
    SubPlacement newSubPl(pl, _spl.getSubset());
    return new FixedYGroupConstraint(newSubPl, _pl);
  }
  // cloning with placement re-association

  bool isSatisfied() const {
    for (unsigned i = 0; i < _spl.getSize(); i++)
      if (_spl[i].y != _pl[i].y) return false;
    return true;
  }

  double getPenalty() const {
    double max = -1;
    for (unsigned i = 0; i < _spl.getSize(); i++)
      max = std::max(max, fabs(_spl[i].y - _pl[i].y));
    return max;
  }

  void enforce() { _spl.pushForwardY(_pl); }

  bool inBBox(const BBox& bb) const {
    for (unsigned i = 0; i < _spl.getSize(); i++)
      if (bb.yMin > _pl[i].y || bb.yMax < _pl[i].y) return false;
    return true;
  }

  Interval getFixedYInterval() const {
    Interval interv(_pl[0].y, _pl[0].y);
    for (unsigned i = 1; i < _spl.getSize(); i++)
      if (interv.low > _pl[i].y)
        interv.low = _pl[i].y;
      else if (interv.high < _pl[i].y)
        interv.high = _pl[i].y;
    return interv;
  }
};
#endif
