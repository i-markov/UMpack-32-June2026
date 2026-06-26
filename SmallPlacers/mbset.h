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

#include <iostream>

const unsigned leading_bit = ~(~0u >> 1);
const unsigned bits_in_uns = sizeof(unsigned) * 8;

class mbitset {
 public:
  unsigned hi;
  unsigned lo;

 public:
  mbitset() : hi(0), lo(0) {}
  mbitset(unsigned hi_, unsigned lo_) : hi(hi_), lo(lo_) {}
  mbitset(const mbitset& bs) : hi(bs.hi), lo(bs.lo) {}
  void reset() { lo = hi = 0; }

  bool lastBitIsOne() const { return lo & 0x01; }
  bool isNonZero() const { return lo || hi; }

  void next() {
    lo >>= 1;
    if (hi & 0x01) lo |= leading_bit;
    hi >>= 1;
  }

  void set(unsigned idx) {
    if (idx < bits_in_uns)
      lo |= (1u << idx);
    else
      hi |= (1u << (idx - bits_in_uns));
  }

  void unset(unsigned idx) {
    if (idx < bits_in_uns)
      lo &= ~(1u << idx);
    else
      hi &= ~(1u << (idx - bits_in_uns));
  }

  bool operator[](unsigned idx) const {
    if (idx < bits_in_uns)
      return (lo >> idx) & 1u;
    else
      return (hi >> (idx - bits_in_uns)) & 1u;
  }

  mbitset& operator=(const mbitset& bs) {
    lo = bs.lo;
    hi = bs.hi;
    return *this;
  }
  mbitset& operator-=(const mbitset& t) {
    lo &= ~t.lo;
    hi &= ~t.hi;
    return *this;
  }
};

inline mbitset operator&(const mbitset& a, const mbitset& b) {
  return mbitset(a.hi & b.hi, a.lo & b.lo);
}

inline std::ostream& operator<<(std::ostream& fout, const mbitset& b) {
  mbitset tmp = b;
  fout << tmp.hi << ":" << tmp.lo << " =(left-to-right!) [";
  while (tmp.isNonZero()) {
    fout << (tmp.lastBitIsOne() ? "1" : "0");
    tmp.next();
  }

  fout << "]" << std::endl;
  return fout;
}
