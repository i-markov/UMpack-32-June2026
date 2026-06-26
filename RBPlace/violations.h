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

// 970930  Original: Paul Tucker
// 980212  dv Modified from pat's code of Legalizer. Improved interfaces
//             for violation code and oriented row-based placement.
// 980603  dv removed orientedRBPlacement class
// 990107  aec 100% re-do of the class.  Both interface and
//			class structure have significantly changed.

#ifndef _VIOLATIONS_H_
#define _VIOLATIONS_H_

#include <iosfwd>

#ifdef _MSC_VER
#pragma warning(disable : 4800)
#endif

class Violations {
  unsigned char _violations;

 public:
  enum ViolationType {
    OffSite = 1,   // cell is not on a valid site
    HOverlap = 2,  // cell overlaps another, horizontally
    VOverlap = 4,  // cell overlaps obstacles with another, vertically
    WOverlap = 8,  // cell overlaps obstacles with special net wire
    ROrient = 16,  // cell violates row orientation
    COrient = 32,  // cell violates own symmetry
    SubRow = 64,   // cell crosses a subrow boundary
    OffRow = 128   // cell not placed in any row
  };

  Violations(unsigned v = 0) : _violations(v) {}

  Violations(const Violations& v) : _violations(v._violations) {}

  Violations(const ViolationType viol) : _violations(viol) {}

  Violations(const char*);  // for mnemonic inline constants

  unsigned getNumViolations() const;

  // are there _any_violations
  bool hasViolations() const { return _violations != 0; }

  // is there a _particular_ violation
  bool hasViolation(const ViolationType viol) const {
    return _violations & viol;
  }

  void clearAll() { _violations = 0; }
  void addViolation(const ViolationType viol) { _violations += viol; }
  void clearViolation(const ViolationType viol) { _violations -= viol; }

  Violations& operator|=(const Violations& v) {
    _violations |= v._violations;
    return *this;
  }

  Violations& operator&=(const Violations& v) {
    _violations &= v._violations;
    return *this;
  }

  Violations operator|(const Violations& v) const {
    return Violations(_violations | v._violations);
  }

  Violations operator&(const Violations& v) const {
    return Violations(_violations & v._violations);
  }

  bool operator!=(const Violations& v) const {
    return _violations != v._violations;
  }

  bool operator==(const Violations& v) const {
    return _violations == v._violations;
  }

  friend std::ostream& operator<<(std::ostream& os, const Violations& v);
};

inline unsigned Violations::getNumViolations() const {
  unsigned count = 0;
  unsigned bitVec = _violations;

  if (bitVec) do {
      count++;
    } while (bitVec &= bitVec - 1);

  return count;
}

#endif
