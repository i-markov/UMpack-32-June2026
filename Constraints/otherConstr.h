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

//! author = "Aug 16, 1997 by Igor Markov"
//! CHANGES ="otherConstr.h  090825 ilm added cloning with placement
//! re-association to all individual Constraint classes"

// CHANGES
// 090825 ilm added cloning with placement re-association to all
//               individual Constraint classes

#ifndef _OTHERCONSTRS_H_
#define _OTHERCONSTRS_H_

#include <ABKCommon/abkassert.h>
#include <Constraints/constraint.h>
#include <Placement/subPlacement.h>

#include <cmath>

/* descriptions follow class declarations */
class TetheredGroupConstraint;
// all ref'ed cells must be tethered (to same tolerance dist) at their locations
class COGConstraint;
// referenced cells must have COG at the given center point
class EqualXConstraint;
// referenced cells must have equal x coordinates
class EqualYConstraint;
// referenced cells must have equal y coordinates
class StayTogetherConstraint;
// referenced cells must be located in a given rectangular region centered
// at their center of gravity

/* -------------------------------------------------------------------- */

//: All refered cells must be tethered (to same tolerance dist) at their
//: locations
class TetheredGroupConstraint : public Constraint {
 protected:
  Placement _pl;
  double _dist;

 protected:
  Constraint* _virtPullBack(const std::vector<unsigned>& hitIndices,
                            const SubPlacement& newSubPl) const;

 public: /* don't use both isSatisfied() and getPenalty() ! */
  TetheredGroupConstraint(const SubPlacement& spl, const Placement& pl,
                          double dist)
      :  // pl is a complete placement that contains the fixed x locations
         //   (we ignore the y locations)
         // spl specifies the group of constrained cells (contains a reference
         //  to the placement that we're checking)
        Constraint(spl, TetheredGroup),
        _pl(pl),
        _dist(dist) {
    abkfatal(pl.getSize() == spl.getSize(),
             "Can't initialize TetheredGroupConstraint: arg size mismatch");
  }

  TetheredGroupConstraint(const TetheredGroupConstraint& fgc)
      : Constraint(fgc._spl, TetheredGroup), _pl(fgc._pl), _dist(fgc._dist){};

  TetheredGroupConstraint& operator=(const TetheredGroupConstraint& fgc) {
    _pl = fgc._pl;
    _spl = fgc._spl;
    _dist = fgc._dist;
    _type = fgc._type;
    return *this;
  }

  virtual ~TetheredGroupConstraint() {}

  Constraint* clone() const { return new TetheredGroupConstraint(*this); }
  Constraint* clone(const Placement& pl) const {
    SubPlacement newSubPl(pl, _spl.getSubset());
    return new TetheredGroupConstraint(newSubPl, _pl, _dist);
  }
  // cloning with placement re-association

  bool isSatisfied() const {
    for (unsigned i = 0; i < _spl.getSize(); i++)
      if (mDist(_spl[i], _pl[i]) > _dist) return false;
    return true;
  }

  double getPenalty() const {
    double max = -1;
    for (unsigned i = 0; i < _spl.getSize(); i++) {
      double curPenalty = mDist(_spl[i], _pl[i]) - _dist;
      if (curPenalty < 0) curPenalty = 0;
      if (curPenalty > max) max = curPenalty;
    }
    return max;
  }

  void enforce() {
    for (unsigned i = 0; i < _spl.getSize(); i++)
      if (mDist(_spl[i], _pl[i]) > _dist) {
        double alpha = atan2(_spl[i].x - _pl[i].x, _spl[i].y - _pl[i].y);
        double s = sin(alpha), c = cos(alpha);
        double tau = 0.9999 / (fabs(s) + fabs(c));
        _spl[i].x = _pl[i].x + tau * _dist * s;
        _spl[i].y = _pl[i].y + tau * _dist * c;
      }
  }

  bool inBBox(const BBox& bb) const {
    for (unsigned i = 0; i < _spl.getSize(); i++)
      if (!bb.contains(_pl[i])) return false;
    return true;
  }
};

//: Center-of-Gravity constraint
//  - referenced cells must have COG at the given center point
class COGConstraint : public Constraint {
 protected:
  Point _center;

 protected:
  Constraint* _virtPullBack(const std::vector<unsigned>& hitIndices,
                            const SubPlacement& newSubPl) const;

 public: /* don't use both isSatisfied() and getPenalty() ! */
  COGConstraint(const SubPlacement& spl, const Point& center)
      : Constraint(spl, COG), _center(center) {}

  COGConstraint(const COGConstraint& cgc)
      : Constraint(cgc._spl, COG), _center(cgc._center) {}

  COGConstraint& operator=(const COGConstraint& cgc) {
    _spl = cgc._spl;
    _center = cgc._center;
    _type = cgc._type;
    return *this;
  }

  virtual ~COGConstraint() {}
  Constraint* clone() const { return new COGConstraint(*this); }
  Constraint* clone(const Placement& pl) const {
    SubPlacement newSubPl(pl, _spl.getSubset());
    return new COGConstraint(newSubPl, _center);
  }
  // cloning with placement re-association

  void enforce() { abkfatal(0, "Can't enforce COGConstraint"); }

  bool isSatisfied() const {
    Point tempPt(0, 0);
    for (unsigned i = 0; i < _spl.getSize(); i++) tempPt += _spl[i];
    tempPt.scaleBy(1.0 / _spl.getSize(), 1.0 / _spl.getSize());
    if (mDist(tempPt, _center) < 1e-5) return true;
    return false;
  }

  double getPenalty() const {
    Point tempPt(0, 0);
    for (unsigned i = 0; i < _spl.getSize(); i++) tempPt += _spl[i];
    tempPt.scaleBy(1.0 / _spl.getSize(), 1.0 / _spl.getSize());
    return mDist(tempPt, _center);
  }

  bool inBBox(const BBox& bb) const { return bb.contains(_center); }
};

//: Referenced cells must have equal x coordinates
class EqualXConstraint : public Constraint {
 protected:
  Constraint* _virtPullBack(const std::vector<unsigned>& hitIndices,
                            const SubPlacement& newSubPl) const;

 public: /* don't use both isSatisfied() and getPenalty() ! */
  EqualXConstraint(const SubPlacement& spl) : Constraint(spl, EqualX) {}

  EqualXConstraint(const EqualXConstraint& rrc)
      : Constraint(rrc._spl, EqualX){};

  EqualXConstraint& operator=(const EqualXConstraint& rrc) {
    _spl = rrc._spl;
    _type = rrc._type;
    return *this;
  }

  virtual ~EqualXConstraint() {}

  Constraint* clone() const { return new EqualXConstraint(*this); }
  Constraint* clone(const Placement& pl) const {
    SubPlacement newSubPl(pl, _spl.getSubset());
    return new EqualXConstraint(newSubPl);
  }
  // cloning with placement re-association

  bool isSatisfied() const {
    for (unsigned i = 1; i < _spl.getSize(); i++)
      if (_spl[i].x != _spl[0].x) return false;
    return true;
  }

  double getPenalty() const {
    double max = -1;
    for (unsigned i = 1; i < _spl.getSize(); i++) {
      double curPenalty = fabs(_spl[i].x - _spl[0].x);
      if (curPenalty > max) max = curPenalty;
    }
    return max;
  }

  void enforce() {
    double sum = 0;
    unsigned i;
    for (i = 0; i < _spl.getSize(); i++) sum += _spl[i].x;
    double mean = sum / _spl.getSize();
    for (i = 0; i < _spl.getSize(); i++) _spl[i].x = mean;
  }

  bool inBBox(const BBox&) const {
    abkfatal(0, "Not implemented");
    return true;
  }
};

//: Referenced cells must have equal y coordinates
class EqualYConstraint : public Constraint {
 protected:
  Constraint* _virtPullBack(const std::vector<unsigned>& hitIndices,
                            const SubPlacement& newSubPl) const;

 public: /* don't use both isSatisfied() and getPenalty() ! */
  EqualYConstraint(const SubPlacement& spl) : Constraint(spl, EqualY) {}

  EqualYConstraint(const EqualYConstraint& rrc)
      : Constraint(rrc._spl, EqualY){};

  EqualYConstraint& operator=(const EqualYConstraint& rrc) {
    _spl = rrc._spl;
    _type = rrc._type;
    return *this;
  }

  virtual ~EqualYConstraint() {}

  Constraint* clone() const { return new EqualYConstraint(*this); }
  Constraint* clone(const Placement& pl) const {
    SubPlacement newSubPl(pl, _spl.getSubset());
    return new EqualYConstraint(newSubPl);
  }
  // cloning with placement re-association

  bool isSatisfied() const {
    for (unsigned i = 1; i < _spl.getSize(); i++)
      if (_spl[i].y != _spl[0].y) return false;
    return true;
  }

  double getPenalty() const {
    double max = -1;
    for (unsigned i = 1; i < _spl.getSize(); i++) {
      double curPenalty = fabs(_spl[i].y - _spl[0].y);
      if (curPenalty > max) max = curPenalty;
    }
    return max;
  }

  void enforce() {
    double sum = 0;
    unsigned i;
    for (i = 0; i < _spl.getSize(); i++) sum += _spl[i].y;
    double mean = sum / _spl.getSize();
    for (i = 0; i < _spl.getSize(); i++) _spl[i].y = mean;
  }

  bool inBBox(const BBox&) const {
    abkfatal(0, "Not implemented");
    return true;
  }
};

//: Referenced cells must be located within the rectange centered at
// their center of gravity
class StayTogetherConstraint : public Constraint {
 protected:
  Rectangle _bb;

 protected:
  Constraint* _virtPullBack(const std::vector<unsigned>& hitIndices,
                            const SubPlacement& newSubPl) const;

 public: /* don't use both isSatisfied() and getPenalty() ! */
  StayTogetherConstraint(const SubPlacement& spl, const Rectangle& bb)
      : Constraint(spl, StayTogether), _bb(bb) {}

  StayTogetherConstraint(const StayTogetherConstraint& rrc)
      : Constraint(rrc._spl, StayTogether), _bb(rrc._bb){};

  StayTogetherConstraint& operator=(const StayTogetherConstraint& rrc) {
    _spl = rrc._spl;
    _bb = rrc._bb;
    _type = rrc._type;
    return *this;
  }

  virtual ~StayTogetherConstraint() {}

  Constraint* clone() const { return new StayTogetherConstraint(*this); }
  Constraint* clone(const Placement& pl) const {
    SubPlacement newSubPl(pl, _spl.getSubset());
    return new StayTogetherConstraint(newSubPl, _bb);
  }
  // cloning with placement re-association

  bool isSatisfied() const {
    Point cog(_spl.getCenterOfGravity());
    BBox bb(_bb);
    bb.centerAt(cog);
    for (unsigned i = 0; i < _spl.getSize(); i++)
      if (!bb.contains(_spl[i])) return false;
    return true;
  }

  double getPenalty() const {
    Point cog(_spl.getCenterOfGravity());
    BBox bb(_bb);
    bb.centerAt(cog);

    double max = -1;
    for (unsigned i = 0; i < _spl.getSize(); i++) {
      Point pt(_spl[i]);
      bb.coerce(pt);
      double curPenalty = mDist(_spl[i], pt);
      if (curPenalty > max) max = curPenalty;
    }
    return max;
  }

  void enforce() {
    Point cog(_spl.getCenterOfGravity());
    for (unsigned i = 0; i < _spl.getSize(); i++) _spl[i] = cog;
  }

  bool inBBox(const BBox&) const {
    abkfatal(0, "Not Implemented");
    return true;
  }
};
#endif
