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

//! author =" July 17, 1997 by Igor Markov"

// 970815 ilm added  support of O(1) tests whether a cell is fixed-constrained
// 970815 ilm added Constraints::allSatisfied()
// 970815 ilm added Constraints::enforceAll()
// 970828 ilm separated restriction from pull-back
// 970828 ilm reimplemented restriction with partial inverse
//                   (was with sorting and binary_search)
// 971121 ilm added Constraints::getFixedBBox()
// 971121 ilm added "add-constraint with re-association"

#ifndef _CONSTRAINTS_H_
#define _CONSTRAINTS_H_

#include <ABKCommon/abkassert.h>
#include <ABKCommon/abktypes.h>
#include <Constraints/constraint.h>
#include <Constraints/fxTypeConstr.h>  // if you need other constraints, include 2 more *.h
#include <Constraints/regionConstr.h>
#include <Placement/placement.h>
#include <Placement/subPlacement.h>
#include <vector>

class Constraints;

/*================================ INTERFACES ==============================*/

//: Keeps track of collections of constraints of different types
class Constraints {
  // ilm: first create an empy object, then throw in various individual
  // constraints with add(), then call compactify() which will sort the
  // constraints by type and size, then merge some of them
  //   add() should take an object of an individual constraint class
  //   (see realConstr.h) with local scope. It will be cloned into
  //   the given Constraints object and then deleted when it goes out of
  //   scope. No pointer in user code ! (but a small loss of efficiency,
  //   never mind)
  // Note: all individual constraints should refer to the same Placement
  //       object

 public:
  enum Filter { ALL, FIXEDTYPE, FIXEDONLY };  // more to come (?)

 protected:
  bool _compactified;
  std::vector<Constraint*> _meat;

  const FixedGroupConstraint* _fgcp;
  const FixedXGroupConstraint* _fxgcp;
  const FixedYGroupConstraint* _fygcp;

 public:
  Constraints() : _compactified(false), _fgcp(0), _fxgcp(0), _fygcp(0) {}
  Constraints(const Constraints& from);
  /* Constraints(const Constraint&  from); */

  Constraints(const Constraints& orig,
              const Subset& restrictTo);  // restriction

  Constraints(const Constraints& orig, const Mapping& pullBackMapping,
              const Placement& pullBackPl);
  // pull-back: the mapping embeds a small placement into the constrained one

  Constraints(const Constraints& orig, const Placement&, Filter flt = ALL);
  // re-association to a different Placement

  ~Constraints();

  void resize(unsigned size) { _meat.resize(size); }
  // in case you know the number of constraints
  // (or an >upper< bound within 130-140% of actual size)

  const Constraint& operator[](unsigned idx) const;
  Constraints& operator=(const Constraints& c);
  unsigned getSize() const { return _meat.size(); }
  bool isEmpty() const { return _meat.size() == 0; }
  bool allSatisfied() const;
  void enforceAll() const;

  /* Note: all individual constraints should refer to the same Placement
           object (assumed below) */

  const Placement& getPlacement() const {
    abkfatal(getSize() > 0, "Empty constraint ");
    return _meat[0]->getPlacement();
  }
  void add(const Constraint& c);
  void add(const Constraint& c, const Placement& pl);  // with re-association
  void sortByType();
  void compactify();

  bool inBBox(const BBox&) const;
  // checks that all indiv. constrs are inBBox()

  bool inBBoxIgnoringSoft(const BBox&) const;
  // checks that all indiv. non-soft constrs are inBBox()

  void reorderSubPlWith(const Permutation&);  // make sure your permutation
                                              // has the correct variance
                                              // ( direct/inverse)
                                              // for your application

  bool isCompactified() const { return _compactified; }

  BBox getFixedBBox() const;

  bool hasFGC() const { return _compactified && _fgcp; }
  const FixedGroupConstraint& getFGC() const {
    abkfatal(_compactified, "Constraints uncompactified");
    return *_fgcp;
  }

  bool hasFXGC() const { return _compactified && _fxgcp; }
  const FixedXGroupConstraint& getFXGC() const {
    abkfatal(_compactified, "Constraints uncompactified");
    return *_fxgcp;
  }

  bool hasFYGC() const { return _compactified && _fygcp; }
  const FixedYGroupConstraint& getFYGC() const {
    abkfatal(_compactified, "Constraints uncompactified");
    return *_fygcp;
  }

  Permutation& getPermutation(Permutation&) const;
  // produces a permutation which sorts fixed items first
  // then fixedX, then fixedY --- pm[i] will give the new
  // location of the cell i. To reorder a vector use
  // template reorderVector() from Placement/mapping.h

  bit_vector& getFGCbitVec(bit_vector& isFixed) const;
  bit_vector& getFXGCbitVec(bit_vector& isXFixed) const;
  bit_vector& getFYGCbitVec(bit_vector& isYFixed) const;
};

std::ostream& operator<<(std::ostream& out, const Constraints& ctrs);

//: Function object used to sort the constraints by their types.
class CompareByTypeSize {
 public:
  CompareByTypeSize() {}
  int operator()(Constraint* c1, Constraint* c2) {
    if (c1->getType() > c2->getType()) return true;
    if (c1->getType() < c2->getType()) return false;
    if (c1->getSize() < c2->getSize()) return true;
    return false;
  }
};

/*============================= IMPLEMENTATIONS ============================*/

/*
inline Constraints::Constraints(const Constraint& ctr)
:_meat(1,ctr.clone()),_compactified(false), _fgcp(0),_fxgcp(0),_fygcp(0)
{
  compactify();
}
*/

inline const Constraint& Constraints::operator[](unsigned idx) const {
  abkassert(idx < getSize(), "Constraint index out of range");
  return *(_meat[idx]);
}

inline void Constraints::add(const Constraint& c) {
  _compactified = false; /* Perhaps, this can be relaxed with an if
                             but be careful ! */
  _meat.push_back(c.clone());
}

inline void Constraints::add(const Constraint& c, const Placement& pl) {
  _compactified = false; /* Perhaps, this can be relaxed with an if
                             but be careful ! */
  _meat.push_back(c.clone(pl));
}

#endif
