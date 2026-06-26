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

//! author ="July 16, 1997 by Igor Markov"
//! CONTACTS =" Igor"

// CHANGES
// 970825 ilm  added Constraint::isFixedType()
// 090825 ilm added cloning with placement re-association to all
//               individual Constraint classes
// 970828 ilm separated restriction from pull-back
// 970828 ilm reimplemented restriction with partial inverse
//                   (was with sorting and binary_search)

#ifndef _CONSTRAINT_H_
#define _CONSTRAINT_H_

#include <ABKCommon/abkassert.h>
#include <Placement/placement.h>
#include <Placement/subPlacement.h>

#include <iostream>
#include <vector>

class Constraint;
class FixedTypeConstraint;  // abstract base class for three classes

/*================================ INTERFACES ==============================*/

//: Base class for individual constraints/groups of constraints
class Constraint {
  friend class Constraints;

 public:
  enum _Type {
    FixedGroup,
    FixedXGroup,
    FixedYGroup,
    EqualX,
    EqualY,
    COG,
    SoftGroup,
    RectRegion,
    MultRectRegion,
    TetheredGroup,
    StayTogether,
    UtilRectRegion
  };

  class Type
  //: Constraint types
  {
    _Type _type;
    static const char* _text;

   public:
    Type(_Type type) : _type(type) {}
    operator const char*() const;
    operator _Type() const { return _type; }
    bool isFixed() const { return _type == FixedGroup; }
    bool isFixedType() const {
      return _type == FixedGroup || _type == FixedXGroup ||
             _type == FixedYGroup;
    }
    bool isRectRgnType() const { return _type == RectRegion; }
    bool isUtilRectRgnType() const { return _type == UtilRectRegion; }

    bool operator<(Type t1) { return _type < t1._type; }
    bool operator<=(Type t1) { return _type <= t1._type; }
    bool operator>(Type t1) { return _type > t1._type; }
    bool operator>=(Type t1) { return _type >= t1._type; }
    bool operator==(Type t1) { return _type == t1._type; }
    bool operator!=(Type t1) { return _type != t1._type; }
  };

 protected:
  Type _type;
  SubPlacement _spl;

  virtual Constraint* _virtPullBack(const std::vector<unsigned>& hitIndices,
                                    const SubPlacement& newSubPl) const = 0;
  /* (old) restriction:  obsolete */
  //  Constraint* _restriction(const std::vector<unsigned>& sortedInx) const;

  /* new restriction */
  Constraint* _restriction(const std::vector<unsigned>& partialInv) const;
  Constraint* _pullBack(const std::vector<unsigned>& partialInverse,
                        const Placement& pullBackPl) const;
  // pull-back
 public:
  Constraint(const SubPlacement& spl, Type type) : _type(type), _spl(spl) {}
  // these three are models for children
  Constraint(const SubPlacement& spl, _Type type) : _type(type), _spl(spl) {}
  // these three are models for children
  Constraint(const Constraint& from) : _type(from._type), _spl(from._spl) {}
  Constraint& operator=(const Constraint& from) {
    abkfatal(0, "Can\'t assign objects of abstract class");
    _spl = from._spl;
    return *this;
  }
  virtual ~Constraint();

  bool isFixedType() const { return _type.isFixedType(); }
  bool isFixed() const { return _type.isFixed(); }

  const Mapping& getMapping() const { return _spl.getMapping(); }
  // synonyms
  const Subset& getSubset() const { return _spl.getSubset(); }
  // synonyms

  const Placement& getPlacement() const { return _spl.getPlacement(); }

  virtual Constraint* clone() const = 0;
  // in child classes will be: { return new <Type>Constraint(*this); }
  virtual Constraint* clone(const Placement&) const = 0;
  // cloning with placement re-association

  virtual void enforce() = 0;
  virtual bool isSatisfied() const = 0;
  virtual double getPenalty() const = 0;
  unsigned getSize() const { return _spl.getSize(); }
  Type getType() const { return _type; }
  virtual bool inBBox(const BBox&) const = 0;  // does not check _spl
  virtual void reorder(const Permutation& pm) { _spl.reorder(pm); }
  virtual bool isSoft() const { return false; }
};

inline bool operator==(Constraint::Type t1, Constraint::_Type t2)
// please Sun CC5.0 - ILM
{
  return static_cast<int>((Constraint::_Type)t1) == static_cast<int>(t2);
}

inline bool operator==(const Constraint::_Type& t1, const Constraint::_Type& t2)
// please "gcc 2.96 Red Hat Linux 7.1 2.96-81" - ILM
{
  return static_cast<int>((Constraint::Type)t2) ==
         static_cast<int>((Constraint::Type)t1);
}

//: Fixed Type Constraint
class FixedTypeConstraint : public Constraint {
 protected:
  Placement _pl;  // fixed locations
 public:
  FixedTypeConstraint(const SubPlacement& spl, const Placement& pl, Type type)
      :  // pl is a complete placement that contains the fixed (x,y) locations
         // spl specifies the group of constrained cells (contains a reference
         //  to the placement that we're checking)
        Constraint(spl, type),
        _pl(pl) {
    abkfatal(pl.getSize() == spl.getSize(),
             "Can't initialize a FixedTypeConstraint: arg size mismatch");
  }

  FixedTypeConstraint(const SubPlacement& spl, const Placement& pl, _Type type)
      : Constraint(spl, type), _pl(pl) {
    abkfatal(pl.getSize() == spl.getSize(),
             "Can't initialize a FixedTypeConstraint: arg size mismatch");
  }
  const Placement& getFixedLocs() const { return _pl; }
};

//============================= IMPLEMENTATIONS ============================

#endif
