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
//! CHANGES =" 090825 ilm added cloning with placement re-association to all
//! individual Constraint classes"

// Created: Aug 16, 1997 by Igor Markov
// CHANGES
// 090825 ilm added cloning with placement re-association to all
//               individual Constraint classes

#ifndef _RECTREGIONCONSTRS_H_
#define _RECTREGIONCONSTRS_H_

#include <ABKCommon/abkassert.h>
#include <Constraints/constraint.h>
#include <Placement/subPlacement.h>

#include <cfloat>
#include <cmath>

// descriptions follow class declarations
class SoftGroupConstraint;
// the dimensions of the BBox of given cells are bounded
// (see LEF/DEF manual p.6-51.)
class RectRegionConstraint;
// all ref'ed cells must be located in the rectangular region
class MultRectRegionConstraint;
// all ref'ed cells must be located in of of given rectangular regions
class UtilRectRegionConstraint;
// by sadya in SYNP. The region must satisfy the given utilization constraint

/* -------------------------------------------------------------------- */

//: The dimensions of the BBox of given cells are bounded
// (see LEF/DEF manual p.6-51.)
class SoftGroupConstraint : public Constraint {
 protected:
  double _maxX;
  double _maxY;
  double _halfPerim;

 protected:
  Constraint* _virtPullBack(const std::vector<unsigned>& hitIndices,
                            const SubPlacement& newSubPl) const;

 public: /* don't use both isSatisfied() and getPenalty() !*/
  SoftGroupConstraint(const SubPlacement& spl, double maxX, double maxY,
                      double halfPerim)
      :  // pl is a complete placement that contains the fixed (x,y) locations
         // spl specifies the group of constrained cells (contains a reference
         //  to the placement that we're checking)
         // maxX, maxY and halfPerim bound the dims of the BBox of given points
        Constraint(spl, SoftGroup),
        _maxX(maxX),
        _maxY(maxY),
        _halfPerim(halfPerim) {
    abkfatal(maxX >= 0 && maxY >= 0 && halfPerim >= 0,
             "Can't initialize SoftGroupConstraint: negative BBox sizes");
  }

  SoftGroupConstraint(const SoftGroupConstraint& fgc)
      : Constraint(fgc._spl, SoftGroup),
        _maxX(fgc._maxX),
        _maxY(fgc._maxY),
        _halfPerim(fgc._halfPerim) {}

  SoftGroupConstraint& operator=(const SoftGroupConstraint& fgc) {
    _spl = fgc._spl;
    _maxX = fgc._maxX;
    _maxX = fgc._maxY;
    _halfPerim = fgc._halfPerim;
    _type = fgc._type;
    return *this;
  }

  ~SoftGroupConstraint();

  Constraint* clone() const { return new SoftGroupConstraint(*this); }
  Constraint* clone(const Placement& pl) const {
    SubPlacement newSubPl(pl, _spl.getSubset());
    return new SoftGroupConstraint(newSubPl, _maxX, _maxY, _halfPerim);
  }
  // cloning with placement re-association

  bool isSatisfied() const {
    BBox bbx;
    for (unsigned i = _spl.getSize() - 1; i != unsigned(-1); i--)
      bbx += _spl[i];
    double dx = fabs(bbx.xMax - bbx.xMin);
    double dy = fabs(bbx.yMax - bbx.yMin);
    return dx <= _maxX && dy <= _maxY && dx + dy <= _halfPerim;
  }

  double getPenalty() const {
    BBox bbx;
    for (unsigned i = _spl.getSize() - 1; i != unsigned(-1); i--)
      bbx += _spl[i];
    double dx = fabs(bbx.xMax - bbx.xMin);
    double dy = fabs(bbx.yMax - bbx.yMin);
    double ddx = dx - _maxX, ddy = dy - _maxY, ddhp = dx + dy - _halfPerim;
    double penal = 0;
    if (ddx > 0) penal += ddx;
    if (ddy > 0) penal += ddy;
    if (ddhp > 0) penal += ddhp;
    return penal;
  }

  void enforce() {
    { abkfatal(0, "Can't enforce SoftGroupConstraint"); }
  }

  bool inBBox(const BBox& bbx) const {
    if (_maxX > bbx.getWidth()) return false;
    if (_maxY > bbx.getHeight()) return false;
    if (_halfPerim > bbx.getHalfPerim()) return false;
    return true;
  }

  virtual bool isSoft() const { return true; }
};

//: Referenced cells must be located in the rectangular region
class RectRegionConstraint : public Constraint {
 protected:
  BBox _bb;

 protected:
  Constraint* _virtPullBack(const std::vector<unsigned>& hitIndices,
                            const SubPlacement& newSubPl) const;

 public: /* don't use both isSatisfied() and getPenalty() ! */
  RectRegionConstraint(const SubPlacement& spl, const BBox& bb)
      : Constraint(spl, RectRegion), _bb(bb) {}

  RectRegionConstraint(const RectRegionConstraint& rrc)
      : Constraint(rrc._spl, RectRegion), _bb(rrc._bb){};

  RectRegionConstraint& operator=(const RectRegionConstraint& rrc) {
    _spl = rrc._spl;
    _bb = rrc._bb;
    _type = rrc._type;
    return *this;
  }

  ~RectRegionConstraint();

  Constraint* clone() const { return new RectRegionConstraint(*this); }
  Constraint* clone(const Placement& pl) const {
    SubPlacement newSubPl(pl, _spl.getSubset());
    return new RectRegionConstraint(newSubPl, _bb);
  }

  const BBox& getBBox() const { return (_bb); }

  bool isSatisfied() const {
    for (unsigned i = 0; i < _spl.getSize(); i++)
      if (!_bb.contains(_spl[i])) return false;
    return true;
  }

  double getPenalty() const {
    double max = -1;
    for (unsigned i = 0; i < _spl.getSize(); i++) {
      Point pt(_spl[i]);
      _bb.coerce(pt);
      double curPenalty = mDist(_spl[i], pt);
      if (curPenalty > max) max = curPenalty;
    }
    return max;
  }

  void enforce() {
    for (unsigned i = 0; i < _spl.getSize(); i++) _bb.coerce(_spl[i]);
  }

  bool inBBox(const BBox& bb) const { return bb.contains(_bb); }
};

//: Referenced cells must be located in one of given rectangular regions
class MultRectRegionConstraint : public Constraint {
 protected:
  std::vector<Rectangle> _regions;

 protected:
  Constraint* _virtPullBack(const std::vector<unsigned>& hitIndices,
                            const SubPlacement& newSubPl) const;

 public: /* don't use both isSatisfied() and getPenalty() ! */
  MultRectRegionConstraint(const SubPlacement& spl,
                           const std::vector<Rectangle>& regions)
      : Constraint(spl, RectRegion), _regions(regions) {}

  MultRectRegionConstraint(const MultRectRegionConstraint& rrc)
      : Constraint(rrc._spl, MultRectRegion), _regions(rrc._regions){};

  MultRectRegionConstraint& operator=(const MultRectRegionConstraint& rrc) {
    _spl = rrc._spl;
    _regions = rrc._regions;
    _type = rrc._type;
    return *this;
  }

  ~MultRectRegionConstraint();

  Constraint* clone() const { return new MultRectRegionConstraint(*this); }
  Constraint* clone(const Placement& pl) const {
    SubPlacement newSubPl(pl, _spl.getSubset());
    return new MultRectRegionConstraint(newSubPl, _regions);
  }

  bool isSatisfied() const {
    for (unsigned i = _spl.getSize() - 1; i != unsigned(-1); i--) {
      bool fits = false;
      for (unsigned k = _regions.size() - 1; k != unsigned(-1) && !fits; k--)
        if (_regions[k].contains(_spl[i])) fits = true;
      if (!fits) return false;
    }
    return true;
  }

  double getPenalty() const {
    double max = -1;
    for (unsigned i = 0; i < _spl.getSize(); i++) {
      double min = DBL_MAX;
      for (unsigned k = _regions.size() - 1; k != unsigned(-1); k--) {
        Point pt(_spl[i]);
        _regions[k].coerce(pt);
        double curPenalty = mDist(_spl[i], pt);
        if (curPenalty < min) min = curPenalty;
      }
      if (min > max) max = min;  // ILM ;-)
    }
    return max;
  }

  void enforce() {
    for (unsigned i = 0; i < _spl.getSize(); i++) _regions[0].coerce(_spl[i]);
  }

  bool inBBox(const BBox& bb) const {
    for (unsigned i = 0; i < _regions.size(); i++)
      if (!bb.contains(_regions[i])) return false;
    return true;
  }
};

//: The rectangular region must satisfy utilization constraint
class UtilRectRegionConstraint : public Constraint {
 protected:
  BBox _bb;
  double _utilization;

 protected:
  // folowing function has no meaning.
  Constraint* _virtPullBack(const std::vector<unsigned>& hitIndices,
                            const SubPlacement& newSubPl) const;

 public: /* don't use both isSatisfied() and getPenalty() ! */
  UtilRectRegionConstraint(const SubPlacement& spl, const BBox& bb,
                           double utilization)
      : Constraint(spl, UtilRectRegion), _bb(bb), _utilization(utilization){};

  UtilRectRegionConstraint(const UtilRectRegionConstraint& rrc)
      : Constraint(rrc._spl, UtilRectRegion),
        _bb(rrc._bb),
        _utilization(rrc._utilization){};

  UtilRectRegionConstraint& operator=(const UtilRectRegionConstraint& rrc) {
    _spl = rrc._spl;
    _bb = rrc._bb;
    _utilization = rrc._utilization;
    _type = rrc._type;
    return *this;
  }

  ~UtilRectRegionConstraint();

  Constraint* clone() const { return new UtilRectRegionConstraint(*this); }

  Constraint* clone(const Placement& pl) const {
    SubPlacement newSubPl(pl, _spl.getSubset());
    return new UtilRectRegionConstraint(newSubPl, _bb, _utilization);
  }

  const BBox& getBBox() const { return (_bb); }
  double getUtilization() const { return (_utilization); }

  bool isSatisfied() const {
    abkwarn(0, "Not Supported Yet\n");
    return false;
  }

  double getPenalty() const {
    abkwarn(0, "Not Supported Yet\n");
    return 0;
  }

  void enforce() { abkwarn(0, "Not Supported Yet\n"); }

  bool inBBox(const BBox& bb) const { return bb.contains(_bb); }
};

#endif
