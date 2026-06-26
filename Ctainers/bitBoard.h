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

//! author="Igor Markov Nov 11, 1997"
// 020811  ilm  ported to g++ 3.0

#ifndef _BITBOARD_H_
#define _BITBOARD_H_

#include <algorithm>
#include <vector>

#include <ABKCommon/abktypes.h>

#include <iostream>

//: bits can be set, and then all set bits can be erased in O(#bits)
// time bits can not be unset one by one in O(1) time each (search required)
// SAMPLE USAGE: mark visited nodes in graph traversals when
// multiple traversals are performed

class BitBoard {
  std::vector<unsigned> _bitIndicesSet;
  bit_vector _boardSpace;

 public:
  BitBoard(unsigned size) : _boardSpace(size, false) {}
  const std::vector<unsigned>& getIndicesOfSetBits() const {
    return _bitIndicesSet;
  }
  // To get the Indices of Bits which have been set to 1;
  const bit_vector& getBoardSpace() const { return _boardSpace; }
  // To get the bit sequence;
  unsigned getSize() const { return _boardSpace.size(); }
  // To get the size of the bitBoard;
  unsigned getNumBitsSet() const { return _bitIndicesSet.size(); }
  // To get the number of bits which have been set to 1;
  bool isBitSet(unsigned idx) const
  // To tell whether the given bit is set to 1;
  {
    /* abkassert(idx<getSize()," Index out of range"); */
    return _boardSpace[idx];
  }
  void setBit(unsigned idx)
  // To set the given bit to 1;
  {
    /*      abkassert(idx<getSize()," Index out of range"); */
    if (!_boardSpace[idx]) {
      _boardSpace[idx] = true;
      _bitIndicesSet.push_back(idx);
    }
  }
  void reset() {
    _bitIndicesSet.clear();
    std::fill(_boardSpace.begin(), _boardSpace.end(), false);
  }
  // Both functions reset and clear are to set all bits to 0;
  // difference is reset() takes O(size) time
  // but clear() takes O(# bits which are true) time

  void clear() {
    for (std::vector<unsigned>::iterator it = _bitIndicesSet.begin();
         it != _bitIndicesSet.end(); it++)
      _boardSpace[*it] = false;
    _bitIndicesSet.clear();
  }

  void reset(unsigned newSize) {
    _bitIndicesSet.clear();
    _boardSpace.assign(newSize, false);
  }

  friend std::ostream& operator<<(std::ostream& os, const BitBoard& bb);
};

#endif
